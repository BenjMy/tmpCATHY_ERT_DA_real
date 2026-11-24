# =============================================================================
# Diagnose one DA run (filter inbreeding, dead members, parameter collapse, ERT data space)   todo.md 5.2
# =============================================================================
#   python diag_pf_run.py SYN_NC_full_upd_pf_N32_s1_unbiased            # run name: looked up in the simulations folder
#   python diag_pf_run.py ~/share/twin_results/runs/SYN_NC_full_upd_pf_N32_s1_unbiased   # or any folder with the small files
#   python diag_pf_run.py RUN --compare SYN_NC_full_upd_enkf_N32_s1_unbiased            # + same checks on a 2nd run, side by side
#   python diag_pf_run.py RUN --json diag.json --strict                                 # machine-readable; exit 1 if any BAD
#
# Reads only SMALL files (no pyCATHY / pyGIMLi needed), found by name in the folder, its sub-folders, <simulations>,
# <simulations>/<run>/, <simulations>/<run>_meta/ and the logs folder of data_path.toml:
#   <run>_run.json          options, per-parameter prior / final mean and sd, rejected members, status
#   <run>_issues.log        every warning of the run (level, kind, count, cycles, message)
#   <run>_rejections.json   members dropped by the run (CATHY stopped early, NaN, ...)
#   <run>_ERT_checks.json   per cycle: forward checks, data-space misfit / spread, normalised innovation
#   <run>_states/cycle_*.npz  (optional) genealogy of the particles, via pf_inbreeding.py (PF runs only)
# A missing file only skips the checks that need it.
#
# What is checked, in the order of the todo (each line: ok / WARN / BAD / info):
#   1  status           interrupted or failed run: the scores are partial
#   2  genealogy        founders left after the last cycle (needs <run>_states/), resampling count, N_eff warnings
#   3  dead members     rejected members, whether they stay NaN (never refilled), how much of the window CATHY ran,
#                       consecutive indices (= common cause)
#   4  parameters       spread contraction, distance of the final mean from the nominal in final SDs (only meaningful when
#                       nominal = truth, i.e. '_unbiased' runs), parameters at a bound, Ks / strictly positive
#                       parameters that went negative before the reflection (perturbed in linear space), reflection
#                       messages whose before = after (the bound printed may not be the one that acts)
#   5  ERT data space   per cycle: ensemble spread / misfit (< 0.5: under-dispersed in data space), share of the variance
#                       in the first mode (>= 0.8: redundant data, a PF treats them as independent), normalised
#                       innovation (obs error too small / too large), members left
#   6  compare          --compare B: parameters of A and B side by side; a bias present in both runs comes from the prior
#                       or the identifiability, not from the filter
#   7  next             the commands / settings to test next (same scenario, other DA, other N), generated from the run name
# The thresholds are the arguments of `diagnose()`; they are rules of thumb, not tests of significance.
# =============================================================================
import argparse
import glob
import json
import math
import os
import re
import sys

OK, WARN, BAD, INFO = 'ok', 'WARN', 'BAD', 'info'
SUFFIX = dict(issues='_issues.log', rejections='_rejections.json', run='_run.json', ert='_ERT_checks.json')
RUN_RE = re.compile(r'^SYN_(?P<case>[A-Za-z]+)_(?P<scenario>.+?)_(?P<da>enkf|pf|esmda)_N(?P<nens>\d+)_s(?P<seed>\d+)'
                    r'(?P<noda>_noDA)?(?P<unbiased>_unbiased)?(?:_(?P<tag>.+))?$')
LN10 = math.log(10.0)


# ------------------------------------------------------------------ locating the files
def _simulation_dirs(name):
    """Folders where the files of `name` may be (data_path.toml); [] when project_paths is not importable."""
    try:
        import project_paths as PP
        P = PP.get()
        sims = P.simulations
        dirs = [sims, os.path.join(sims, name), os.path.join(sims, name + '_meta')]
        logs = getattr(P, 'logs', None)
        if logs:
            dirs.append(logs)
        return dirs
    except Exception:
        return []


def _glob_first(patterns):
    for p in patterns:
        hits = sorted(glob.glob(p, recursive=True))
        if hits:
            return hits[0]
    return None


def locate(run, extra_dirs=()):
    """{'name', 'issues', 'rejections', 'run', 'ert', 'states'} (paths or None). `run` = run name or a folder.
    Folders given by the user (and their sub-folders) are searched freely; the parent folder and the
    <simulations> folders only for the exact file names of THIS run (never a sibling run)."""
    if os.path.isdir(run):
        base = os.path.abspath(run)
        name = os.path.basename(base.rstrip(os.sep))
        hit = _glob_first([os.path.join(base, '*' + SUFFIX['run']), os.path.join(base, '*', '*' + SUFFIX['run'])])
        if hit:
            name = os.path.basename(hit)[:-len(SUFFIX['run'])]
        wide = [base] + list(extra_dirs)
        parent = os.path.dirname(base)
        exact = [parent, os.path.join(parent, name + '_meta')] + _simulation_dirs(name)
    else:
        name = run
        wide = list(extra_dirs)
        exact = _simulation_dirs(name)
    out = {'name': name}
    for key, suf in SUFFIX.items():
        pats = [os.path.join(d, name + suf) for d in wide + exact]
        pats += [os.path.join(d, '*', name + suf) for d in wide]
        out[key] = _glob_first(pats)
    cand = [os.path.join(d, name + '_states') for d in wide + exact] + [os.path.join(d, '*', name + '_states') for d in wide]
    st = next((c for c in cand if glob.glob(os.path.join(c, 'cycle_*.npz'))), None)
    out['states'] = st
    return out


def _load_json(path):
    if not path:
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


# ------------------------------------------------------------------ the issues log
def _cycles(txt):
    """'3' -> (3, 3);  '0-9' -> (0, 9);  anything else -> (None, None)."""
    m = re.fullmatch(r'\s*(\d+)(?:-(\d+))?\s*', txt or '')
    if not m:
        return None, None
    a = int(m.group(1))
    return a, int(m.group(2)) if m.group(2) else a


def parse_issues(path):
    """[{level, kind, count, c0, c1, msg}] from <run>_issues.log (tab separated; whitespace fallback)."""
    out = []
    if not path:
        return out
    with open(path, errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split('\t', 4)
            if len(parts) < 5:
                parts = re.split(r'\s+', line.strip(), 4)
            if len(parts) < 5:
                continue
            level, kind, cnt, cyc, msg = parts
            c0, c1 = _cycles(cyc)
            m = re.fullmatch(r'x?(\d+)', cnt.strip())
            out.append(dict(level=level.strip(), kind=kind.strip(), count=int(m.group(1)) if m else 1,
                            c0=c0, c1=c1, msg=msg.strip()))
    return out


_REFLECT = re.compile(r'^(?P<par>\S+): (?P<n>\d+) analysed value\(s\) reflected to \[(?P<lo>[^,\]]+), (?P<hi>[^\]]+)\] '
                      r'\(members \[(?P<mem>[^\]]*)\]\): \[(?P<old>[^\]]*)\] -> \[(?P<new>[^\]]*)\]')


def _floats(txt):
    out = []
    for tok in txt.replace(',', ' ').split():
        try:
            out.append(float(tok))
        except ValueError:
            pass
    return out


def parse_reflections(issues):
    """One dict per 'reflected' message: parameter, members, old / new values, bounds, cycles."""
    out = []
    for it in issues:
        m = _REFLECT.match(it['msg'])
        if not m:
            continue
        mem = [int(x) for x in re.findall(r'\d+', m.group('mem'))]
        out.append(dict(par=m.group('par'), members=mem, old=_floats(m.group('old')), new=_floats(m.group('new')),
                        lo=float(m.group('lo')), hi=float(m.group('hi')), c0=it['c0'], c1=it['c1'],
                        n_events=it['count']))
    return out


# ------------------------------------------------------------------ the checks (each returns a list of (level, text))
def check_status(run, issues, rej):
    out = []
    status = (run or {}).get('library_status')
    ncyc = (run or {}).get('n_cycles')
    errs = [i for i in issues if i['level'].upper() == 'ERROR']
    if status is None and run is None:
        out.append((INFO, 'no <run>_run.json: status unknown'))
    elif str(status).lower() in ('ok', 'finished', 'success'):
        out.append((OK, f'library status {status}, {ncyc} cycles'))
    else:
        out.append((WARN, f'library status "{status}" after {ncyc} cycles: the scores are partial'))
    for e in errs:
        aborted = 'KeyboardInterrupt' in e['msg'] or 'aborted' in e['msg']
        out.append((WARN if aborted else BAD, f'ERROR at cycle {e["c0"]}: {e["msg"][:120]}'
                    + (' (run interrupted by hand?)' if 'KeyboardInterrupt' in e['msg'] else '')))
    return out


def check_weights(issues, da_type):
    out = []
    if da_type not in (None, 'pf'):
        return [(INFO, f'weights: not applicable to {da_type}')]
    conc = [i for i in issues if re.search(r'weight concentration|all weight on one particle', i['msg'], re.I)]
    neff = [i for i in issues if re.search(r'low N_eff|particle degeneracy', i['msg'], re.I)]
    for label, group in (('extreme weight concentration', conc), ('very low N_eff', neff)):
        if not group:
            out.append((OK, f'no "{label}" warning'))
            continue
        n = sum(g['count'] for g in group)
        c0 = min(g['c0'] for g in group if g['c0'] is not None)
        c1 = max(g['c1'] for g in group if g['c1'] is not None)
        out.append((BAD, f'"{label}" {n} times, cycles {c0}-{c1}: the weights are degenerate'))
    return out


def check_genealogy(states_dir, da_type, frac=0.25):
    if da_type not in (None, 'pf'):
        return [(INFO, f'genealogy: only for the particle filter (this run is {da_type})')], None
    if not states_dir:
        return [(INFO, 'genealogy: no <run>_states/cycle_*.npz here (pruned or not copied): run pf_inbreeding.py '
                       'on the workstation, or copy the folder')], None
    try:
        import pf_inbreeding as PI
        df, v = PI.genealogy(PI.load_cycles(states_dir), frac)
    except Exception as exc:                                        # noqa: BLE001 (report, do not crash a diagnosis)
        return [(WARN, f'genealogy: could not be computed ({type(exc).__name__}: {exc})')], None
    n = v['nens']
    txt = (f"{v['final_founders']} distinct founders left of N={n} (ESS {v['final_founders_ess']:.1f}), resampled at "
           f"{v['resampled_cycles']} of {v['n_analysis_cycles']} analysis cycles")
    if len(df):
        txt += (f"; mean spread of psi {df.spread_before.iloc[0]:.3g} -> {df.spread_after.iloc[-1]:.3g} m "
                f"(first before / last after the analysis)")
    return [(BAD if v['inbreeding'] else OK, txt)], v


def check_members(run, rej, issues, nens=None):
    """Dead members: rejections.json + NaN parameters that persist in the log."""
    out, info = [], {}
    rows = rej or []
    ids = sorted({r['member'] for r in rows})
    if run and run.get('rejected_members'):
        for lst in run['rejected_members'].values():
            ids = sorted(set(ids) | set(lst))
    info['rejected'] = ids
    if not ids:
        out.append((OK, 'no rejected member'))
    else:
        cyc = sorted({r['cycle'] for r in rows})
        stages = sorted({r['stage'] for r in rows})
        out.append((WARN, f'{len(ids)} member(s) rejected: {ids} (cycle(s) {cyc}, stage {", ".join(stages)})'))
        fr = []
        for r in rows:
            m = re.search(r't=([\d.]+) s.*last printing time ([\d.]+) s', r.get('reason', ''))
            if m and float(m.group(2)) > 0:
                fr.append(float(m.group(1)) / float(m.group(2)))
        if fr:
            out.append((INFO, f'CATHY stopped at {100 * min(fr):.0f}-{100 * max(fr):.0f} % of the window: it did not '
                              'converge (time step collapse?), look at the parameters of these members'))
        if len(ids) >= 3 and ids == list(range(ids[0], ids[0] + len(ids))):
            out.append((WARN, f'the rejected indices {ids[0]}-{ids[-1]} are consecutive: a common cause (same worker / batch '
                              '/ folder) is more likely than bad luck'))
    nan_members, last = set(), None
    for r in parse_reflections(issues):
        for mem, v in zip(r['members'], r['old'] if len(r['old']) == len(r['members']) else []):
            if math.isnan(v):
                nan_members.add(mem)
                last = max(last or 0, r['c1'] or 0)
    info['nan_members'] = sorted(nan_members)
    if nan_members:
        alive = (nens - len(nan_members)) if nens else None
        out.append((BAD, f'members {sorted(nan_members)} still have NaN parameters up to cycle {last}: dead slots are not '
                         f'refilled' + (f', the real ensemble is {alive} of {nens}' if alive else '')))
    elif ids:
        out.append((OK, 'no NaN parameter left in the log after the rejection'))
    return out, info


def check_parameters(run, issues, nominal_is_truth, contraction_bad=0.6, offset_bad=3.0):
    out, rows = [], []
    pars = (run or {}).get('parameters') or {}
    upd = {k: v for k, v in pars.items() if v.get('n_updates', 0) > 0}
    if not upd:
        return [(INFO, 'no updated parameter in <run>_run.json')], rows
    for k, v in upd.items():
        fsd, psd = v.get('final_sd'), v.get('prior_sd')
        nom = v.get('nominal')
        rows.append(dict(par=k, contraction=v.get('contraction'),
                         offset_final_sd=(v['final_mean'] - nom) / fsd if nom is not None and fsd else None,
                         offset_prior_sd=(v['final_mean'] - nom) / psd if nom is not None and psd else None,
                         final_mean=v.get('final_mean'), final_sd=fsd, nominal=nom,
                         at_bound=v.get('frac_at_bound_final', 0.0)))
    # contraction as the run reports it, or recomputed when absent
    for r in rows:
        if r['contraction'] is None and upd[r['par']].get('prior_sd'):
            r['contraction'] = 1 - upd[r['par']]['final_sd'] / upd[r['par']]['prior_sd']
    coll = [r for r in rows if r['contraction'] is not None and r['contraction'] >= contraction_bad]
    lvl = BAD if len(coll) > len(rows) / 2 else WARN if coll else OK
    out.append((lvl, f'{len(coll)} of {len(rows)} updated parameters lost >= {100 * contraction_bad:.0f} % of their spread: '
                     + ', '.join(f"{r['par']} {100 * r['contraction']:.0f} %" for r in coll)))
    away = [r for r in rows if r['offset_final_sd'] is not None and abs(r['offset_final_sd']) >= offset_bad]
    if nominal_is_truth:
        out.append((BAD if away else OK,
                    (f'{len(away)} parameters end >= {offset_bad:g} final SDs from the truth (nominal = truth): '
                     + ', '.join(f"{r['par']} {r['offset_final_sd']:+.1f}" for r in away)
                     + ' -> confident and wrong') if away else 'final means are within the final SDs of the truth'))
    else:
        out.append((INFO, 'offset from the nominal not judged (the run is biased on purpose, or --nominal-is-truth no): '
                          + ', '.join(f"{r['par']} {r['offset_final_sd']:+.1f}" for r in away)))
    bound = [r['par'] for r in rows if r['at_bound']]
    if bound:
        out.append((WARN, f'parameters with members at a bound at the end: {bound}'))
    refl = parse_reflections(issues)
    neg = {}
    for r in refl:
        spec = pars.get(r['par'], {})
        lo_clip = spec.get('clip_min', r['lo'])
        if lo_clip is not None and lo_clip > 0 and any(v < 0 for v in r['old']):
            neg.setdefault(r['par'], []).append(min(v for v in r['old'] if v < 0))
    if neg:
        out.append((WARN, 'negative analysed values for strictly positive parameters before the reflection ('
                          + ', '.join(f'{p} min {min(v):.3g}' for p, v in neg.items())
                          + '): they are perturbed in linear space; use log space (Ks, thetar, a_Archie...)'))
    same = {}
    for r in refl:
        a, b = [x for x in r['old'] if not math.isnan(x)], [x for x in r['new'] if not math.isnan(x)]
        if a and a == b:
            same.setdefault(r['par'], []).append(max(b))
    if same:
        out.append((WARN, 'reflection messages whose values before and after are identical (the log prints the same list): '
                          + ', '.join(f'{p} (max {max(v):.4g})' for p, v in same.items())
                          + '. Check that the bound printed is the one that acts (values piling up just below a round '
                            'number, e.g. 0.5, point to another bound)'))
    return out, rows


def ert_table(checks):
    """[{cycle, n_members, n_obs, misfit_log10, rel_spread, spread_over_misfit, first_mode, rms_innov, n_fwd_failed, n_fwd}]"""
    by = {}
    for ent in checks or []:
        c = ent.get('cycle')
        row = by.setdefault(c, dict(cycle=c))
        if ent.get('stage') == 'forward':
            mem = ent.get('members', {})
            row['n_fwd'] = len(mem)
            row['n_fwd_failed'] = sum(1 for m in mem.values() if not m.get('ok', True))
            st = (ent.get('ensemble') or {}).get('stats', {})
        elif ent.get('stage') == 'analysis':
            st = (ent.get('report') or {}).get('stats', {})
        else:
            continue
        for src, dst in (('n_members', 'n_members'), ('n_obs', 'n_obs'), ('median_abs_log10_misfit', 'misfit_log10'),
                         ('median_rel_spread', 'rel_spread'), ('first_mode_var_frac', 'first_mode'),
                         ('rms_normalised_innovation', 'rms_innov')):
            if st.get(src) is not None:
                row[dst] = st[src]
    rows = [r for _, r in sorted(by.items(), key=lambda kv: (kv[0] is None, kv[0]))]
    for r in rows:
        if r.get('misfit_log10') and r.get('rel_spread') is not None:
            r['spread_over_misfit'] = r['rel_spread'] / (LN10 * r['misfit_log10'])
    return [r for r in rows if 'misfit_log10' in r]


def check_ert(rows, nens=None, ratio_low=0.5, first_mode=0.8, innov=(0.3, 1.5)):
    if not rows:
        return [(INFO, 'no <run>_ERT_checks.json (or no ensemble statistics in it)')]
    out = []
    ratios = [r['spread_over_misfit'] for r in rows if 'spread_over_misfit' in r]
    if ratios:
        med = sorted(ratios)[len(ratios) // 2]
        out.append((WARN if med < ratio_low else OK,
                    f'data-space spread / misfit (rel. spread / (ln10 x log10 misfit)), median over cycles {med:.2f} '
                    f'[{min(ratios):.2f}-{max(ratios):.2f}]: '
                    + ('under-dispersed in data space' if med < ratio_low else
                       'about the size of the misfit, so the data are fitted with a spread that is not collapsed there')))
    fm = [r['first_mode'] for r in rows if 'first_mode' in r]
    if fm:
        share = sum(1 for x in fm if x >= first_mode) / len(fm)
        out.append((WARN if share >= 0.5 else OK,
                    f'first mode carries >= {100 * first_mode:.0f} % of the ensemble variance in {100 * share:.0f} % of the '
                    f'cycles (range {min(fm):.2f}-{max(fm):.2f})'
                    + (': the many ERT data are largely redundant; a PF multiplies them as if independent (weight '
                       'collapse) and they constrain one combination of parameters, not the depth profile' if share >= 0.5 else '')))
    ri = [r['rms_innov'] for r in rows if 'rms_innov' in r]
    if ri:
        lo, hi = innov
        bad = min(ri) < lo or max(ri) > hi
        out.append((WARN if bad else OK, f'RMS normalised innovation {min(ri):.2f}-{max(ri):.2f} (expect ~1; < 1 = obs error '
                                         f'larger than the misfit, > 1 = smaller)'))
    mem = [r['n_members'] for r in rows if 'n_members' in r]
    if mem and nens and min(mem) < nens:
        out.append((WARN, f'members used in the ERT analysis: {max(mem)} at cycle 0, {min(mem)} later (N = {nens})'))
    failed = [(r['cycle'], r['n_fwd_failed'], r['n_fwd']) for r in rows if r.get('n_fwd_failed')]
    if failed:
        tot = sum(f for _, f, _ in failed)
        out.append((INFO, f'forward ERT checks failed for {tot} member-cycles (|rhoa| outside the model range at a few '
                          'quadrupoles): mostly near-zero potentials, check `rhoa_implausible` if it grows'))
    return out


def suggest(name, nens, da_type):
    """Commands / settings to test next, from the run name (RUN_RE)."""
    m = RUN_RE.match(name or '')
    out = []
    if not m:
        return [(INFO, 'run name not recognised (SYN_<case>_<scenario>_<da>_N<n>_s<seed>[_noDA][_unbiased]): no command '
                       'generated')]
    d = m.groupdict()
    bias = 'none' if d['unbiased'] else 'biased'
    env = f"CASE={d['case']} BIAS={bias} SEED={d['seed']}"
    sc = d['scenario']
    out.append((INFO, f'3  same scenario / seed / N with the EnKF:   {env} DA=enkf NENS={d["nens"]} ./run_twin.sh pair {sc}'))
    out.append((INFO, f'   and ES-MDA (parameters only):           {env} DA=esmda NENS={d["nens"]} ./run_twin.sh pair {sc}'))
    out.append((INFO, '   then:  python diag_pf_run.py ' + name + ' --compare '
                + name.replace('_' + d['da'] + '_', '_enkf_', 1)))
    if da_type in (None, 'pf'):
        out.append((INFO, '4  PF settings to change one at a time (da_ert_synthetic.py, the pf_* arguments; each rerun '
                          'is compared with this one): pf_resample_threshold 1.0 -> 0.5 (resample only when N_eff < N/2, '
                          'check first that the weights are carried between cycles); pf_jitter_std_state 0 -> > 0; '
                          'pf_temper_neff_frac (does the tempering iterate up to full weights, or stop after one partial '
                          'step?); do not update a_/m_/n_Archie in this twin; use the PF for the parameters only'))
    out.append((INFO, '5  Ks and other strictly positive parameters: perturb in log space (pyCATHY perturbate / scenario '
                      'sampling), then no value needs to be reflected'))
    out.append((INFO, f'6  ensemble size, only after 2-4:   {env} DA={da_type or d["da"]} NENS_LIST="{nens} 64 100" '
                      f'./run_twin.sh compare-nens   (one diagnostic run with N=100 is enough to see whether N matters)'))
    return out


# ------------------------------------------------------------------ one run, all checks
def diagnose(run_arg, extra_dirs=(), nominal_is_truth='auto', frac=0.25):
    files = locate(run_arg, extra_dirs)
    name = files['name']
    run = _load_json(files['run'])
    rej = _load_json(files['rejections'])
    ert = _load_json(files['ert'])
    issues = parse_issues(files['issues'])
    opts = (run or {}).get('options', {})
    m = RUN_RE.match(name)
    da_type = opts.get('da_type') or (m.group('da') if m else None)
    nens = opts.get('n_members') or (int(m.group('nens')) if m else None)
    nom_truth = (('_unbiased' in name) if nominal_is_truth == 'auto' else nominal_is_truth == 'yes')
    rep = dict(name=name, files={k: v for k, v in files.items() if k != 'name'}, da_type=da_type, nens=nens,
               sections=[])

    def add(title, items):
        rep['sections'].append(dict(title=title, items=[dict(level=l, text=t) for l, t in items]))

    add('1  status', check_status(run, issues, rej))
    gen_items, gen = check_genealogy(files['states'], da_type, frac)
    add('2  genealogy and weights', check_weights(issues, da_type) + gen_items)
    mem_items, mem = check_members(run, rej, issues, nens)
    add('3  dead members', mem_items)
    par_items, par_rows = check_parameters(run, issues, nom_truth)
    add('4  parameters', par_items)
    rows = ert_table(ert)
    add('5  ERT data space', check_ert(rows, nens))
    rep['parameters'], rep['ert'], rep['members'], rep['genealogy'] = par_rows, rows, mem, gen
    rep['n_issues'] = len(issues)
    return rep


def compare(a, b):
    """Parameter table A vs B and a one-line reading."""
    pa = {r['par']: r for r in a['parameters']}
    pb = {r['par']: r for r in b['parameters']}
    common = [p for p in pa if p in pb]
    rows, both = [], []
    for p in common:
        ra, rb = pa[p], pb[p]
        rows.append(dict(par=p, a_contraction=ra['contraction'], b_contraction=rb['contraction'],
                         a_offset_prior_sd=ra['offset_prior_sd'], b_offset_prior_sd=rb['offset_prior_sd']))
        oa, ob = ra['offset_prior_sd'], rb['offset_prior_sd']
        if oa is not None and ob is not None and abs(oa) >= 1 and abs(ob) >= 1 and oa * ob > 0:
            both.append(p)
    off_a = [r['par'] for r in rows if r['a_offset_prior_sd'] is not None and abs(r['a_offset_prior_sd']) >= 1]
    if not off_a:
        msg = (OK, 'A has no parameter more than 1 prior SD from the nominal: nothing to attribute')
    elif len(both) >= 0.5 * len(off_a):
        msg = (WARN, f'{len(both)} of the {len(off_a)} parameters that A moves by >= 1 prior SD move the same way in B '
                     f'({", ".join(both)}): the bias is not specific to the filter (prior, or identifiability)')
    else:
        msg = (WARN, f'only {len(both)} of the {len(off_a)} parameters that A moves by >= 1 prior SD move the same way in B: '
                     'the offsets are specific to the filter of A')
    return dict(rows=rows, message=dict(level=msg[0], text=msg[1]))


# ------------------------------------------------------------------ output
def _fmt(x, spec='.2f'):
    return '' if x is None else format(x, spec)


def print_report(rep, out=sys.stdout):
    w = out.write
    w(f"\n=== {rep['name']}   (da = {rep['da_type']}, N = {rep['nens']}, {rep['n_issues']} issue lines) ===\n")
    miss = [k for k, v in rep['files'].items() if not v]
    if miss:
        w(f"files not found: {', '.join(SUFFIX.get(k, k) if k in SUFFIX else k + '/' for k in miss)}\n")
    for s in rep['sections']:
        w(f"\n[{s['title']}]\n")
        for it in s['items']:
            w(f"  {it['level']:<5} {it['text']}\n")
    if rep['parameters']:
        w('\n  parameters         contraction  offset/final SD  offset/prior SD   final mean      final SD\n')
        for r in rep['parameters']:
            w(f"  {r['par']:<18} {_fmt(r['contraction']):>10}  {_fmt(r['offset_final_sd'], '+.1f'):>15}  "
              f"{_fmt(r['offset_prior_sd'], '+.1f'):>15}  {_fmt(r['final_mean'], '.4g'):>12}  {_fmt(r['final_sd'], '.3g'):>10}\n")
    if rep['ert']:
        w('\n  cycle  n_mem   n_obs  misfit_log10  rel_spread  spread/misfit  first_mode  rms_innov\n')
        for r in rep['ert']:
            w(f"  {str(r['cycle']):>5}  {_fmt(r.get('n_members'), '.0f'):>5}  {_fmt(r.get('n_obs'), '.0f'):>6}  "
              f"{_fmt(r.get('misfit_log10'), '.4f'):>12}  {_fmt(r.get('rel_spread'), '.3f'):>10}  "
              f"{_fmt(r.get('spread_over_misfit')):>13}  {_fmt(r.get('first_mode')):>10}  {_fmt(r.get('rms_innov')):>9}\n")


def print_compare(a, b, cmp, out=sys.stdout):
    w = out.write
    w(f"\n=== compare  A = {a['name']}\n             B = {b['name']} ===\n")
    w('\n  parameters         contraction A   B     offset/prior SD  A      B\n')
    for r in cmp['rows']:
        w(f"  {r['par']:<18} {_fmt(r['a_contraction']):>12} {_fmt(r['b_contraction']):>6}   "
          f"{_fmt(r['a_offset_prior_sd'], '+.1f'):>14} {_fmt(r['b_offset_prior_sd'], '+.1f'):>6}\n")
    w(f"\n  {cmp['message']['level']:<5} {cmp['message']['text']}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description='Diagnose a DA run from its small files (inbreeding, dead members, parameters, '
                                             'ERT data space); optional comparison with a second run')
    ap.add_argument('run', help='run name (looked up with data_path.toml) or a folder holding the small files')
    ap.add_argument('--compare', metavar='RUN_B', help='second run (e.g. the EnKF at the same seed / N / scenario)')
    ap.add_argument('--dir', action='append', default=[], help='extra folder to search for the files (repeatable)')
    ap.add_argument('--nominal-is-truth', choices=('auto', 'yes', 'no'), default='auto',
                    help='judge the offset of the final mean from the nominal (auto: only for "_unbiased" runs)')
    ap.add_argument('--frac', type=float, default=0.25, help='INBREEDING if the founders left <= frac x N')
    ap.add_argument('--json', metavar='PATH', help='write everything as json')
    ap.add_argument('--strict', action='store_true', help='exit code 1 when any check is BAD')
    a = ap.parse_args(argv)

    rep = diagnose(a.run, a.dir, a.nominal_is_truth, a.frac)
    print_report(rep)
    cmp = None
    if a.compare:
        rep_b = diagnose(a.compare, a.dir, a.nominal_is_truth, a.frac)
        cmp = compare(rep, rep_b)
        print_report(rep_b)
        print_compare(rep, rep_b, cmp)
        rep['compare'] = dict(other=rep_b['name'], **cmp)
    print('\n[7  next]')
    for lvl, txt in suggest(rep['name'], rep['nens'], rep['da_type']):
        print(f'  {txt}')
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump(rep, fh, indent=1, default=lambda o: None if o is None else str(o))
    n_bad = sum(1 for s in rep['sections'] for it in s['items'] if it['level'] == BAD)
    print(f"\n{n_bad} BAD, {sum(1 for s in rep['sections'] for it in s['items'] if it['level'] == WARN)} WARN")
    return 1 if (a.strict and n_bad) else 0


if __name__ == '__main__':
    sys.exit(main())
