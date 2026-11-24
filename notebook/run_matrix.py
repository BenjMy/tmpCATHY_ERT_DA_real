# =============================================================================
# Scenario matrix of the twin experiment: calls nb7 for every combination
# =============================================================================
#   python run_matrix.py --dry-run                      # list the runs
#   python run_matrix.py --cases NC --seeds 42 43 44    # run, skipping the finished ones
#   python run_matrix.py --preset quick                 # NC, ic + VGP_upd, EnKF, 1 seed, NENS 32
# Run from the notebook/ folder. A run is "finished" when <sim-dir>/<run>_meta/run_meta.json says so
# (resume after a crash = run the same command again; --force reruns). Open-loop references
# (--no-da) are added per (case, scenario, seed) with --open-loop. Failed runs do not stop the matrix.
import argparse
import json
import os
import subprocess
import sys
import time

import project_paths as PP
import run_registry as RR
from twin_names import run_name

PRESETS = {
    'quick': dict(cases=['NC'], scenarios=['ic', 'ic_VGP_upd'], da=['enkf'], seeds=[42], nens=32),
    'main': dict(cases=['NC', 'HC'], scenarios=['ic', 'ic_VGP', 'ic_VGP_upd', 'ic_Ks_upd', 'ic_Archie_upd',
                                                 'full_upd'], da=['enkf', 'pf'], seeds=[42, 43, 44], nens=100),
    'ens': dict(cases=['NC'], scenarios=['ic_VGP_upd'], da=['enkf'], seeds=[42, 43, 44, 45, 46], nens=None),
}
ap = argparse.ArgumentParser()
ap.add_argument('--preset', choices=sorted(PRESETS))
ap.add_argument('--cases', nargs='+', default=None)
ap.add_argument('--scenarios', nargs='+', default=None)
ap.add_argument('--da', nargs='+', default=None, choices=['enkf', 'pf', 'esmda'])
ap.add_argument('--seeds', nargs='+', type=int, default=None)
ap.add_argument('--nens', nargs='+', type=int, default=None, help="several values = ensemble-size study")
ap.add_argument('--prior-bias', choices=['biased', 'none'], default='biased')
ap.add_argument('--open-loop', action='store_true', help='also run each combination without assimilation')
ap.add_argument('--extra', default='', help='extra nb7 arguments, e.g. "--obs-err 10 --inflation 1.2"')
ap.add_argument('--delta', type=float, default=None)
ap.add_argument('--root', default=None, help='project root (default: data_path.toml)')
ap.add_argument('--sim-dir', default=None, help='default: [outputs].simulations of data_path.toml')
ap.add_argument('--script', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nb7_DA_ERT_synthetic.py'))
ap.add_argument('--log-dir', default=None, help='default: [outputs].logs of data_path.toml')
ap.add_argument('--dry-run', action='store_true')
ap.add_argument('--force', action='store_true')
ap.add_argument('--stop-on-fail', action='store_true')
a = ap.parse_args()
PTH = PP.get()        # (P below is the preset)
a.root = a.root or PTH.root_str
a.log_dir = a.log_dir or PTH.logs

P = PRESETS.get(a.preset, dict(cases=['NC'], scenarios=['ic_VGP_upd'], da=['enkf'], seeds=[42], nens=100))
cases = a.cases or P['cases']
scenarios = a.scenarios or P['scenarios']
das = a.da or P['da']
seeds = a.seeds or P['seeds']
nens_list = a.nens or [P['nens'] or 100] if a.preset != 'ens' or a.nens else [32, 64, 128]
sim_dir = a.sim_dir or (PTH.simulations if PTH.is_root(a.root) else os.path.join(a.root, 'simulations'))
biased = a.prior_bias == 'biased'


def status(name):
    f = os.path.join(sim_dir, name + '_meta', 'run_meta.json')
    try:
        with open(f) as fh:
            return json.load(fh).get('status', '?')
    except (OSError, ValueError):
        return 'none'


runs = []


def _cmd(case, sc, da, n, seed, with_da):
    name = run_name(case, sc, da, n, seed, with_da, biased)
    cmd = [sys.executable, a.script, '--case', case, '--scenario', sc, '--da-type', da,
           '--nens', str(n), '--seed', str(seed), '--prior-bias', a.prior_bias,
           '--root', a.root] + (['--no-da'] if not with_da else []) \
          + (['--delta', str(a.delta)] if a.delta else []) \
          + (['--sim-dir', a.sim_dir] if a.sim_dir else []) + a.extra.split()
    return name, cmd


for case in cases:
    for sc in scenarios:
        for n in nens_list:
            for seed in seeds:
                for da in das:
                    runs.append(_cmd(case, sc, da, n, seed, True))
                if a.open_loop:                       # ONE open loop per scenario/N/seed (no update: same for every DA type)
                    runs.append(_cmd(case, sc, 'enkf', n, seed, False))
print(f"{len(runs)} runs")
os.makedirs(a.log_dir, exist_ok=True)
summary = []
for k, (name, cmd) in enumerate(runs, 1):
    st = status(name)
    if st == 'finished' and not a.force:
        print(f"[{k}/{len(runs)}] {name}: finished, skipped")
        summary.append((name, 'skipped'))
        continue
    print(f"[{k}/{len(runs)}] {name} (previous status: {st})")
    if a.dry_run:
        print("   ", ' '.join(cmd))
        continue
    t0 = time.time()
    with open(os.path.join(a.log_dir, name + '.log'), 'w') as log:
        rc = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT)
    res = 'ok' if rc == 0 and status(name) == 'finished' else f'FAILED (rc={rc}, status={status(name)})'
    if not res == 'ok':      # nb7 may have died before it could register itself (import error, bad argument ...)
        _row = next((r for r in RR.read() if r['run'] == name), None)
        if _row is None or _row['status'] in ('started', ''):
            RR.upsert(name, status='failed', finished=time.strftime('%Y-%m-%d %H:%M:%S'),
                      log=os.path.join(a.log_dir, name + '.log'))
    print(f"    {res} in {(time.time() - t0) / 60:.1f} min -> {a.log_dir}/{name}.log")
    summary.append((name, res))
    if res != 'ok' and a.stop_on_fail:
        break
if not a.dry_run:
    with open(os.path.join(a.log_dir, 'matrix_summary.json'), 'w') as f:
        json.dump(summary, f, indent=1)
    bad = [s for s in summary if s[1].startswith('FAILED')]
    print(f"done: {len(summary) - len(bad)} ok/skipped, {len(bad)} failed")
    sys.exit(1 if bad else 0)
