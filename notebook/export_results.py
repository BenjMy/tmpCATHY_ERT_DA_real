"""Package the results of the DA runs for a colleague: ONE small folder (or zip), nothing recomputed.

    python export_results.py OUT_DIR [--runs GLOB ...] [--only-ok] [--figures none|compare|all] [--zip] [--toml FILE]

Reads the paths of data_path.toml (registry, simulations, nb8 folder) and COPIES what already exists:
  OUT_DIR/overview.csv           ONE row per run: registry (status, duration, commit) + disk (lighten_run.py) + the key
                                 nb8 scores of summary.csv. The only new file; it recomputes nothing.
  OUT_DIR/README.txt             how to read it (run names, columns, caveats)
  OUT_DIR/tables/                nb8 summary.csv / parameters_summary.csv (one subfolder per case) + the runs registry
  OUT_DIR/figures/               nb8 figures: compare_*.png (default), all of them (--figures all) or none
  OUT_DIR/runs/<run>/            manifest.json, run_meta.json, <run>_issues.log and the tail of the run log (failed /
                                 partial runs only): exact command, seed, git state
  OUT_DIR/synthetic_twin_plan.md meaning of the scenarios
NOT copied: states (_states/, .npz, .nc), vtk, pickles: big, and unreadable without pyCATHY.
"""
import argparse
import csv
import fnmatch
import glob
import os
import shutil
import sys

import project_paths as PP
import run_registry as RR
import lighten_run as L

SCORE_COLS = ['sw_all_rmse', 'sw_all_rmse_ol', 'sw_all_skill_vs_ol', 'sw_all_spread_skill', 'sw_all_coverage',
              'sw_all_crps', 'sw_in_ERT_block_rmse', 'sw_outside_ERT_block_rmse', 'psi_all_rmse', 'psi_all_coverage',
              'states_kind']
REG_COLS = ['run', 'case', 'scenario', 'da_type', 'nens', 'seed', 'prior_bias', 'with_da', 'status', 'duration_min',
            'git_commit', 'git_branch']
TAIL_LINES = 80
OK_STATUS = ('finished', 'ok')

README = """Results of the DA twin experiments (exported by export_results.py)

START HERE
  overview.csv      one row per run. Same run name in tables/, runs/<run>/ and figures/.
  figures/          compare_*.png = scenarios / DA types side by side; per-run figures with --figures all
  tables/           the nb8 tables as written by the code (summary.csv, parameters_summary.csv), the runs registry

RUN NAMES   SYN_<case>_<scenario>_<da>_N<members>_s<seed>[_noDA][_unbiased]
  case NC = uncompacted soil, HC = compacted soil; scenario = what is perturbed / updated (see synthetic_twin_plan.md);
  da = enkf | pf | esmda; N = ensemble size; s = seed; _noDA = open loop (no assimilation, same prior);
  _unbiased = prior centred on the truth, otherwise biased.

OVERVIEW COLUMNS
  status            finished / failed / started (still running or killed): scores of a partial run cover its cycles only
  duration_min      wall time of the run
  disk_mb, vtk_mb   size of the run folder and of its vtk files; compact_store = states stored as _states/cycle_*.npz or .nc
  sw_all_rmse       RMSE of the ensemble mean water content vs the truth, all nodes (lower = better), mean over cycles
  sw_all_rmse_ol    the same for the open loop; sw_all_skill_vs_ol = 1 - rmse / rmse_ol (> 0: the DA beats the open loop)
  sw_all_coverage   share of the truth inside the 5-95 % ensemble range (0.9 = reliable; far below = overconfident)
  sw_all_spread_skill  ensemble spread / RMSE (about 1 = consistent)
  sw_all_crps       CRPS (lower = better); sw_in_ERT_block_rmse / sw_outside_ERT_block_rmse: inside / outside the ERT block
  psi_*             the same for the pressure head
CAVEATS
  * enkf / pf scores are those of the FORECAST before each update; esmda scores are posterior (smoother) states
    conditioned on all the data: not like for like (column states_kind).
  * Runs with a different N, seed or prior_bias are not directly comparable: compare one option at a time.
  * Empty score cells: run without nb8 scores (failed, no states, or not yet plotted).
  * runs/<run>/manifest.json has the exact command, seed and git commit; <run>_issues.log / log tail explain failed runs.
"""


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def list_runs(P, globs=None, only_ok=False):
    """Registry rows (dict) of the selected runs, plus run folders on disk the registry does not know."""
    _sims = os.path.abspath(P.simulations)       # the registry is shared by all folders: keep this folder's rows
    rows = {r['run']: dict(r) for r in RR.read(P.run_registry)
            if not r.get('folder') or os.path.dirname(os.path.abspath(r['folder'])) == _sims}
    if os.path.isdir(P.simulations):
        for d in sorted(os.listdir(P.simulations)):
            full = os.path.join(P.simulations, d)
            if (os.path.isdir(full) and d.startswith('SYN_') and not d.endswith('_meta') and d not in rows):
                rows[d] = dict(run=d, status='no registry row')
    out = []
    for name, r in rows.items():
        if globs and not any(fnmatch.fnmatch(name, g) for g in globs):
            continue
        if only_ok and r.get('status') not in OK_STATUS:
            continue
        out.append(r)
    return sorted(out, key=lambda r: r['run'])


def _find_one(patterns):
    for pat in patterns:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    return None


def _tail(path, n=TAIL_LINES):
    try:
        with open(path, errors='replace') as fh:
            return ''.join(fh.readlines()[-n:])
    except OSError:
        return None


def collect_run(P, row, out_dir):
    """Copy the small per-run files. Returns the list of what was copied."""
    run = row['run']
    dest = os.path.join(out_dir, 'runs', run)
    meta = os.path.join(P.simulations, run + '_meta')
    copied = []

    def cp(src, name=None):
        if src and os.path.isfile(src):
            os.makedirs(dest, exist_ok=True)
            shutil.copy2(src, os.path.join(dest, name or os.path.basename(src)))
            copied.append(name or os.path.basename(src))

    cp(os.path.join(meta, 'manifest.json'))
    cp(os.path.join(meta, 'run_meta.json'))
    cp(_find_one([os.path.join(d, run + '_issues.log') for d in (P.simulations, meta, P.logs, os.path.join(P.simulations, run))]))
    if row.get('status') != 'finished':                       # failed / started / unknown: the end of the log tells why
        log = row.get('log') or os.path.join(P.logs, run + '.log')
        txt = _tail(log) if os.path.isfile(log) else None
        if txt:
            os.makedirs(dest, exist_ok=True)
            with open(os.path.join(dest, run + '_log_tail.txt'), 'w') as fh:
                fh.write(txt)
            copied.append(run + '_log_tail.txt')
    return copied


def disk_columns(P, run):
    d = os.path.join(P.simulations, run)
    if not os.path.isdir(d):
        return dict(disk_mb='', vtk_mb='', compact_store='')
    cats, _ = L.scan(d)
    return dict(disk_mb=round(sum(cats.values()) / 1e6, 1), vtk_mb=round(cats.get('vtk', 0) / 1e6, 1),
                compact_store=L.has_compact_store(d))


def nb8_tables(P):
    """{case: folder} of the nb8 output folders holding a summary.csv (figures_nb8/<case>/, or figures_nb8/ itself)."""
    base = P.figures_synth
    found = {}
    if os.path.isfile(os.path.join(base, 'summary.csv')):
        found[''] = base
    if os.path.isdir(base):
        for d in sorted(os.listdir(base)):
            if os.path.isfile(os.path.join(base, d, 'summary.csv')):
                found[d] = os.path.join(base, d)
    return found


def read_summaries(tables):
    """{run: {score col: value}} from every summary.csv."""
    out = {}
    for folder in tables.values():
        with open(os.path.join(folder, 'summary.csv'), newline='') as fh:
            for r in csv.DictReader(fh):
                if r.get('run'):
                    out[r['run']] = {c: r[c] for c in SCORE_COLS if c in r}
    return out


def export(P, out_dir, globs=None, only_ok=False, figures='compare', make_zip=False):
    runs = list_runs(P, globs, only_ok)
    os.makedirs(out_dir, exist_ok=True)
    tables = nb8_tables(P)
    scores = read_summaries(tables)
    # tables + figures: copied as they are
    for case, folder in tables.items():
        for name in ('summary.csv', 'parameters_summary.csv'):
            if os.path.isfile(os.path.join(folder, name)):
                dst = os.path.join(out_dir, 'tables', case)
                os.makedirs(dst, exist_ok=True)
                shutil.copy2(os.path.join(folder, name), dst)
        if figures != 'none':
            for png in sorted(glob.glob(os.path.join(folder, '*.png'))):
                base = os.path.basename(png)
                if figures == 'all' and not base.startswith('compare_') and not any(
                        base.startswith(r['run'] + '_') for r in runs):
                    continue
                if figures == 'compare' and not base.startswith('compare_'):
                    continue
                dst = os.path.join(out_dir, 'figures', case)
                os.makedirs(dst, exist_ok=True)
                shutil.copy2(png, dst)
    if os.path.isfile(P.run_registry):
        os.makedirs(os.path.join(out_dir, 'tables'), exist_ok=True)
        shutil.copy2(P.run_registry, os.path.join(out_dir, 'tables', os.path.basename(P.run_registry)))
    plan = os.path.join(P.notebook, 'synthetic_twin_plan.md')
    if os.path.isfile(plan):
        shutil.copy2(plan, out_dir)
    # overview.csv
    score_cols = [c for c in SCORE_COLS if any(c in v for v in scores.values())]
    cols = REG_COLS + ['disk_mb', 'vtk_mb', 'compact_store'] + score_cols + ['files']
    with open(os.path.join(out_dir, 'overview.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in runs:
            row = {c: r.get(c, '') for c in REG_COLS}
            row.update(disk_columns(P, r['run']))
            row.update(scores.get(r['run'], {}))
            row['files'] = ' '.join(collect_run(P, r, out_dir))
            w.writerow(row)
    with open(os.path.join(out_dir, 'README.txt'), 'w') as fh:
        fh.write(README)
    result = out_dir
    if make_zip:
        result = shutil.make_archive(out_dir.rstrip(os.sep), 'zip', os.path.dirname(os.path.abspath(out_dir)),
                                     os.path.basename(os.path.abspath(out_dir)))
    return runs, result


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out_dir')
    ap.add_argument('--runs', nargs='*', default=None, help="glob(s) of run names, e.g. 'SYN_NC_*'")
    ap.add_argument('--only-ok', action='store_true', help='only runs with status finished')
    ap.add_argument('--figures', choices=['none', 'compare', 'all'], default='compare')
    ap.add_argument('--zip', action='store_true')
    ap.add_argument('--toml', default=None, help='data_path.toml (default: the project one)')
    a = ap.parse_args(argv)
    P = PP.get(a.toml) if a.toml else PP.get()
    runs, res = export(P, a.out_dir, a.runs, a.only_ok, a.figures, a.zip)
    print(f"{len(runs)} run(s) -> {res}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
