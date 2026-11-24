# =============================================================================
# ONE entry point of the synthetic-twin scenario study: runs (nb7) -> diagnostics (nb8)
# =============================================================================
# Grid = case x scenario x DA type (EnKF / PF) x seed x NENS, each with DA, plus ONE open loop (withoutDA) per
# case/scenario/seed/NENS (the open loop does not depend on the DA type). The scenarios are a ladder in the number of
# perturbed / updated parameters (make_scenarios.py): ic < ic_Ks_upd < ic_VGP_upd < ic_VGP_Ks_upd < full_upd.
#
#   python run_study.py --dry-run          # print the runs and the diagnostics call, simulate nothing
#   python run_study.py                    # QUICK: first 10 survey times, NC, EnKF + PF, DA + open loop, then diagnostics
#   python run_study.py --full             # all surveys (sim dir <root>simulations), NENS 100
#   python run_study.py --n-times 5 --nens 16 --scenarios ic ic_Ks_upd
#   python run_study.py --stage diag       # diagnostics only, on every finished run of the sim dir
# Run from anywhere (the children run in notebook/, ROOT = ../). Quick runs go to <root>simulations_n<N>/ so they never
# mix with full runs. Resumable: finished runs are skipped (--force reruns). A failed run does not stop the study
# (--stop-on-fail does); every failure is listed at the end and in <sim-dir>/study_report.json; exit code 1 if any.
# Diagnostics = nb8_plots_synthetic.py on ALL finished runs of each case (scores vs truth, parameters, DA vs open loop).
# =============================================================================
import argparse
import json
import os
import subprocess
import sys
import time

from twin_names import run_name

HERE = os.path.dirname(os.path.abspath(__file__))
ALL_SCENARIOS = ['ic', 'ic_VGP', 'ic_VGP_upd', 'ic_Ks_upd', 'ic_Archie_upd', 'ic_VGP_Ks_upd', 'full_upd']
LADDER = ['ic', 'ic_Ks_upd', 'ic_VGP_upd', 'ic_VGP_Ks_upd', 'full_upd']
QUICK_NENS, FULL_NENS = 32, 100


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cases', nargs='+', default=['NC'], choices=['NC', 'HC'])
    ap.add_argument('--scenarios', nargs='+', default=LADDER,
                    help=f"tags of make_scenarios.py ({', '.join(ALL_SCENARIOS)}) or 'all'; default = the ladder")
    ap.add_argument('--da', nargs='+', default=['enkf', 'pf'], choices=['enkf', 'pf', 'esmda'])
    ap.add_argument('--no-open-loop', action='store_true', help='skip the withoutDA runs')
    ap.add_argument('--seeds', nargs='+', type=int, default=[42])
    ap.add_argument('--nens', nargs='+', type=int, default=None, help=f'default {QUICK_NENS} (quick) / {FULL_NENS} (--full)')
    ap.add_argument('--n-times', type=int, default=10, help='quick test: first N survey times only (default 10)')
    ap.add_argument('--full', action='store_true', help='all survey times (= --n-times 0)')
    ap.add_argument('--prior-bias', choices=['biased', 'none'], default='biased')
    ap.add_argument('--extra', default='', help='extra nb7 arguments, e.g. "--obs-err 10 --inflation 1.2"')
    ap.add_argument('--diag-extra', default='', help='extra nb8 arguments, e.g. "--detail-seed 42 --t-from-h 0"')
    ap.add_argument('--stage', choices=['runs', 'diag', 'both'], default='both')
    ap.add_argument('--root', default='../', help='project folder (data/, simulations/); relative to notebook/')
    ap.add_argument('--sim-dir', default=None, help='default <root>simulations_n<N> (quick) / <root>simulations (full)')
    ap.add_argument('--log-dir', default=None, help='default <sim-dir>/logs')
    ap.add_argument('--nb7', default=os.path.join(HERE, 'nb7_DA_ERT_synthetic.py'))
    ap.add_argument('--nb8', default=os.path.join(HERE, 'nb8_plots_synthetic.py'))
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--force', action='store_true', help='rerun finished runs')
    ap.add_argument('--stop-on-fail', action='store_true')
    a = ap.parse_args(argv)
    assert a.n_times >= 0, f"--n-times must be >= 0, got {a.n_times}"
    assert all(s >= 0 for s in a.seeds), "seeds must be >= 0"
    assert a.nens is None or all(n >= 2 for n in a.nens), f"--nens must be >= 2, got {a.nens}"
    return a


def effective_n_times(a):
    return 0 if a.full else a.n_times


def resolve_root(root):
    """Absolute root WITH a trailing separator (nb7 / nb8 build paths as ROOT + 'data/...')."""
    r = root if os.path.isabs(root) else os.path.normpath(os.path.join(HERE, root))
    return r.rstrip(os.sep) + os.sep


def resolve_sim_dir(a, root):
    if a.sim_dir:
        return os.path.abspath(a.sim_dir if os.path.isabs(a.sim_dir) else os.path.join(HERE, a.sim_dir))
    n = effective_n_times(a)
    return os.path.join(root, 'simulations' if n == 0 else f'simulations_n{n}')


def normalize_scenarios(names):
    """Tags without the case prefix, 'all' expanded, duplicates dropped (order kept); unknown tag = error."""
    out = []
    for n in names:
        tags = ALL_SCENARIOS if n == 'all' else [n[3:] if n[:3] in ('NC_', 'HC_') else n]
        for t in tags:
            assert t in ALL_SCENARIOS, f"unknown scenario {n!r}; available: {ALL_SCENARIOS}"
            if t not in out:
                out.append(t)
    return out


def _run(a, root, sim_dir, case, sc, da, nens, seed, with_da):
    name = run_name(case, sc, da, nens, seed, with_da, a.prior_bias == 'biased')
    cmd = [sys.executable, a.nb7, '--case', case, '--scenario', sc, '--da-type', da, '--nens', str(nens),
           '--seed', str(seed), '--prior-bias', a.prior_bias, '--root', root, '--sim-dir', sim_dir]
    if not with_da:
        cmd.append('--no-da')
    if effective_n_times(a):
        cmd += ['--n-times', str(effective_n_times(a))]
    cmd += a.extra.split()
    return dict(name=name, cmd=cmd, case=case, scenario=sc, da=da, nens=nens, seed=seed,
                kind='DA' if with_da else 'open_loop')


def plan_runs(a, root, sim_dir):
    """The ordered list of runs (dicts with name, cmd, ...). ONE open loop per case/scenario/NENS/seed (nb7 needs --da-type enkf)."""
    nens_list = a.nens or [FULL_NENS if effective_n_times(a) == 0 else QUICK_NENS]
    runs = []
    for case in a.cases:
        for sc in normalize_scenarios(a.scenarios):
            for n in nens_list:
                for seed in a.seeds:
                    runs += [_run(a, root, sim_dir, case, sc, da, n, seed, True) for da in a.da]
                    if not a.no_open_loop:
                        runs.append(_run(a, root, sim_dir, case, sc, 'enkf', n, seed, False))
    return runs


def read_meta(sim_dir, name):
    try:
        with open(os.path.join(sim_dir, name + '_meta', 'run_meta.json')) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def read_status(sim_dir, name):
    return read_meta(sim_dir, name).get('status', 'none')


def log_tail(path):
    """Last non-empty line of a log (why a child failed); '' if unreadable."""
    try:
        with open(path, errors='replace') as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        return lines[-1][:300] if lines else ''
    except OSError:
        return ''


def call_logged(cmd, log_path, env=None):
    """Run cmd in notebook/ with stdout+stderr in log_path; returns (returncode, error text)."""
    try:
        with open(log_path, 'w') as log:
            return subprocess.call(cmd, cwd=HERE, stdout=log, stderr=subprocess.STDOUT, env=env), ''
    except OSError as e:
        return -1, f'{type(e).__name__}: {e}'


def execute_runs(runs, a, sim_dir, log_dir):
    results = []
    for k, r in enumerate(runs, 1):
        res = dict(run=r['name'], kind=r['kind'], scenario=r['scenario'], da=r['da'], nens=r['nens'], seed=r['seed'],
                   result='', minutes=0.0, reason='')
        prev = read_status(sim_dir, r['name'])
        if prev == 'finished' and not a.force:
            res['result'] = 'skipped'
            print(f"[{k}/{len(runs)}] {r['name']}: finished, skipped")
        elif a.dry_run:
            res['result'] = 'dry-run'
            print(f"[{k}/{len(runs)}] {r['name']} (previous status: {prev})\n    {' '.join(r['cmd'])}")
        else:
            print(f"[{k}/{len(runs)}] {r['name']} (previous status: {prev})", flush=True)
            log = os.path.join(log_dir, r['name'] + '.log')
            t0 = time.time()
            rc, err = call_logged(r['cmd'], log)
            res['minutes'] = round((time.time() - t0) / 60, 2)
            st = read_status(sim_dir, r['name'])
            if rc == 0 and st == 'finished':
                res['result'] = 'ok'
                res['n_updated_parameters'] = read_meta(sim_dir, r['name']).get('n_updated_parameters')
            else:
                res['result'] = 'FAILED'
                res['reason'] = err or f"rc={rc}, status={st}: {log_tail(log)}"
            print(f"    {res['result']} in {res['minutes']:.1f} min -> {log}" + (f"\n    {res['reason']}" if res['reason'] else ''))
        results.append(res)
        if res['result'] == 'FAILED' and a.stop_on_fail:
            print('stopping: --stop-on-fail')
            break
    return results


def execute_diagnostics(a, root, sim_dir, log_dir):
    """nb8 on every finished run of each case. A case with no finished run is reported as FAILED, not skipped silently."""
    results = []
    env = dict(os.environ, MPLBACKEND='Agg')
    for case in a.cases:
        cmd = [sys.executable, a.nb8, '--case', case, '--root', root, '--sim-dir', sim_dir] + a.diag_extra.split()
        res = dict(case=case, result='', reason='', summary_csv=None, log=os.path.join(log_dir, f'diagnostics_{case}.log'))
        if a.dry_run:
            res['result'] = 'dry-run'
            print(f"diagnostics {case}:\n    {' '.join(cmd)}")
        elif not os.path.isdir(sim_dir):
            res.update(result='FAILED', reason=f'sim dir not found: {sim_dir}')
        else:
            print(f"diagnostics {case} ...", flush=True)
            rc, err = call_logged(cmd, res['log'], env=env)
            csv = os.path.join(sim_dir, 'figures_nb8', case, 'summary.csv')
            if rc == 0:
                res.update(result='ok', summary_csv=csv if os.path.isfile(csv) else None)
            else:
                res.update(result='FAILED', reason=err or f"rc={rc}: {log_tail(res['log'])}")
            print(f"    {res['result']} -> {os.path.dirname(csv)}" + (f"\n    {res['reason']}" if res['reason'] else ''))
        results.append(res)
    return results


def summarize(run_res, diag_res):
    """(text, n_failures): counts and every failure, one line each."""
    n = {k: sum(r['result'] == k for r in run_res) for k in ('ok', 'skipped', 'FAILED', 'dry-run')}
    lines = [f"runs: {len(run_res)} planned | ok {n['ok']} | skipped (finished) {n['skipped']} | FAILED {n['FAILED']}"
             + (f" | dry-run {n['dry-run']}" if n['dry-run'] else '')]
    lines += [f"  FAILED {r['run']}: {r['reason']}" for r in run_res if r['result'] == 'FAILED']
    bad = [d for d in diag_res if d['result'] == 'FAILED']
    if diag_res:
        lines.append(f"diagnostics: {len(diag_res)} case(s) | FAILED {len(bad)}")
        lines += [f"  FAILED diagnostics {d['case']}: {d['reason']} (log {d['log']})" for d in bad]
        lines += [f"  summary: {d['summary_csv']}" for d in diag_res if d['summary_csv']]
    return '\n'.join(lines), n['FAILED'] + len(bad)


def main(argv=None):
    a = parse_args(argv)
    root = resolve_root(a.root)
    sim_dir = resolve_sim_dir(a, root)
    log_dir = os.path.abspath(a.log_dir) if a.log_dir else os.path.join(sim_dir, 'logs')
    runs = plan_runs(a, root, sim_dir)
    n_t = effective_n_times(a)
    print(f"study: {'FULL (all survey times)' if n_t == 0 else f'QUICK (first {n_t} survey times)'} | "
          f"{len(runs)} runs | sim dir {sim_dir}")
    if not a.dry_run and (a.stage != 'diag' or os.path.isdir(sim_dir)):     # diagnostics only never create the sim dir
        os.makedirs(log_dir, exist_ok=True)
    run_res, diag_res = [], []
    if a.stage in ('runs', 'both'):
        run_res = execute_runs(runs, a, sim_dir, log_dir)
    if a.stage in ('diag', 'both') and not (a.stop_on_fail and any(r['result'] == 'FAILED' for r in run_res)):
        diag_res = execute_diagnostics(a, root, sim_dir, log_dir)
    text, n_fail = summarize(run_res, diag_res)
    print('\n' + text)
    if not a.dry_run and os.path.isdir(sim_dir):
        with open(os.path.join(sim_dir, 'study_report.json'), 'w') as f:
            json.dump(dict(n_times=n_t, sim_dir=sim_dir, runs=run_res, diagnostics=diag_res, n_failures=n_fail,
                           argv=sys.argv, finished=time.strftime('%Y-%m-%d %H:%M:%S')), f, indent=1, default=str)
    return 1 if n_fail else 0


if __name__ == '__main__':
    sys.exit(main())
