"""Protocol of todo.md section 7: measure a run (7.0) and accept / reject a lighter setting (7.4).

    python lighter_protocol.py measure OUT.json -- CMD ARGS...        # run CMD, record wall time, peak RSS, disk per category
    python lighter_protocol.py accept  BASE.csv [BASE.csv ...] --light LIGHT.csv [LIGHT.csv ...] [--cols ...]

measure   needs nothing but the standard library (peak RSS = ru_maxrss of the children, so it also works without
          /usr/bin/time). Pass --run-dir RUN_DIR (before the `--`) to add the `lighten_run.py` byte count per category
          of that folder after the command ended.
accept    each csv is a nb8 `summary.csv` of ONE seed. A lighter setting is accepted when, for every score column,
          the lighter value (mean over its seeds) is inside the seed-to-seed spread of the baseline:
              |mean_light - mean_base| <= max(k * std_base, tol)       (default k = 1, tol = 0 )
          (std over the baseline seeds, ddof=1; at least 2 baseline seeds are needed). Exit code 0 = accepted.
          Rows are matched on (scenario, da, nens): other nens/da rows of the file are compared with their own baseline.
"""
import argparse
import json
import os
import resource
import subprocess
import sys
import time

DEFAULT_COLS = ['sw_all_rmse', 'sw_all_spread', 'sw_all_coverage', 'sw_all_crps',
                'psi_all_rmse', 'psi_all_spread', 'psi_all_coverage', 'psi_all_crps']
KEYS = ['scenario', 'da', 'nens']


def measure(cmd, run_dir=None):
    """Run cmd (list). Returns dict(command, returncode, wall_s, peak_rss_mb, disk_bytes{category}, disk_total_bytes)."""
    t0 = time.time()
    rc = subprocess.call(cmd)
    wall = time.time() - t0
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss          # kB on Linux, bytes on macOS
    mb = rss / 1024 / 1024 if sys.platform == 'darwin' else rss / 1024
    out = dict(command=cmd, returncode=rc, wall_s=round(wall, 1), peak_rss_mb=round(mb, 1))
    if run_dir:
        import lighten_run as L
        cats, _ = L.scan(run_dir)
        out['disk_bytes'] = cats
        out['disk_total_bytes'] = sum(cats.values())
    return out


def _load(paths):
    import pandas as pd
    return pd.concat([pd.read_csv(p).assign(_file=os.path.basename(p)) for p in paths], ignore_index=True)


def accept(base_paths, light_paths, cols=None, k=1.0, tol=0.0):
    """Returns (ok, table DataFrame: key, column, base_mean, base_std, light_mean, diff, limit, ok)."""
    import numpy as np
    import pandas as pd
    if len(base_paths) < 2:
        raise ValueError("need >= 2 baseline seeds to know the seed-to-seed spread (todo 7.4: 3 seeds)")
    base, light = _load(base_paths), _load(light_paths)
    cols = [c for c in (cols or DEFAULT_COLS) if c in base.columns and c in light.columns]
    if not cols:
        raise ValueError("none of the score columns is in both summary.csv files")
    keys = [c for c in KEYS if c in base.columns and c in light.columns]
    rows = []
    for key, gb in base.groupby(keys, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        gl = light
        for kname, kval in zip(keys, key):
            gl = gl[gl[kname] == kval]
        if not len(gl):
            continue
        for c in cols:
            b = gb[c].astype(float).dropna()
            lv = gl[c].astype(float).dropna()
            if len(b) < 2 or not len(lv):
                continue
            sd = float(b.std(ddof=1))
            diff = float(lv.mean() - b.mean())
            lim = max(k * sd, tol)
            rows.append(dict(zip(keys, key), column=c, base_mean=float(b.mean()), base_std=sd,
                             light_mean=float(lv.mean()), diff=diff, limit=lim, ok=bool(abs(diff) <= lim)))
    tab = pd.DataFrame(rows)
    if tab.empty:
        raise ValueError("no (scenario, da, nens) row of the lighter run matches the baseline")
    return bool(tab['ok'].all()), tab


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    m = sub.add_parser('measure')
    m.add_argument('out')
    m.add_argument('--run-dir', default=None)
    a = sub.add_parser('accept')
    a.add_argument('base', nargs='+')
    a.add_argument('--light', nargs='+', required=True)
    a.add_argument('--cols', nargs='+', default=None)
    a.add_argument('-k', type=float, default=1.0, help='allowed shift in baseline standard deviations')
    a.add_argument('--tol', type=float, default=0.0, help='absolute floor of the allowed shift')
    cmd = []
    if '--' in argv:                                   # everything after the first -- is the command to measure
        i = argv.index('--')
        argv, cmd = argv[:i], argv[i + 1:]
    ns = ap.parse_args(argv)
    if ns.cmd == 'measure':
        if not cmd:
            ap.error("give the command after --")
        res = measure(cmd, ns.run_dir)
        os.makedirs(os.path.dirname(os.path.abspath(ns.out)), exist_ok=True)
        with open(ns.out, 'w') as fh:
            json.dump(res, fh, indent=1)
        print(json.dumps({k: v for k, v in res.items() if k != 'disk_bytes'}))
        return res['returncode']
    ok, tab = accept(ns.base, ns.light, ns.cols, ns.k, ns.tol)
    print(tab.round(5).to_string(index=False))
    print('ACCEPTED' if ok else 'REJECTED')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
