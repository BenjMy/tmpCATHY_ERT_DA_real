"""CSV registry of the runs: what has been done and what has not.

One row per run name in [outputs].run_registry of data_path.toml (default registry/runs_registry.csv),
updated in place. Each tracked run also gets <simulations>/<run>_meta/manifest.json (code versions + branches of
CATHY DA ERT / pyCATHY / CATHY, inputs hashes, environment: see run_manifest.py); `track(..., manifest=False)` skips it. Safe with parallel runs (file lock + atomic rewrite). A registry problem never stops a run:
it is reported in one line.

In a driver script:
    import run_registry as RR
    REG = RR.track(run_name, kind='synthetic', case='NC', scenario='NC_ic', da_type='enkf', nens=64, seed=1,
                   prior_bias='biased', with_da=True, log=...)          # status 'started'
    ...
    REG.finish('finished')            # or 'failed'; if the process ends without it, the row becomes 'failed'

Command line:
    python run_registry.py list [--status failed] [--case NC] [--kind real]
    python run_registry.py summary
    python run_registry.py backfill [--sim-dir DIR]              # register runs already on disk
"""
import atexit
import csv
import json
import os
import re
import subprocess
import sys
import time

COLUMNS = ['run', 'kind', 'case', 'scenario', 'da_type', 'nens', 'seed', 'prior_bias', 'with_da', 'status',
           'started', 'finished', 'duration_min', 'folder', 'log', 'git_commit', 'git_branch', 'pycathy_version',
           'manifest']
_NOW = '%Y-%m-%d %H:%M:%S'


def registry_path(path=None):
    if path:
        return os.path.abspath(path)
    import project_paths as PP
    return PP.get().run_registry


def _git_commit():
    try:
        out = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=os.path.dirname(os.path.abspath(__file__)),
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() if out.returncode == 0 else ''
    except Exception:
        return ''


class _Lock:
    def __init__(self, path):
        self.path = path + '.lock'

    def __enter__(self):
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        self.f = open(self.path, 'w')
        try:
            import fcntl
            fcntl.flock(self.f, fcntl.LOCK_EX)
        except ImportError:                      # Windows: no lock, rare concurrent writes
            pass
        return self

    def __exit__(self, *a):
        self.f.close()


def read(path=None):
    """List of dict rows ([] if the file does not exist)."""
    p = registry_path(path)
    if not os.path.isfile(p):
        return []
    with open(p, newline='') as f:
        return list(csv.DictReader(f))


def _write(rows, p):
    tmp = p + '.tmp'
    with open(tmp, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, p)


def upsert(run, path=None, **fields):
    """Create or update the row of `run` (unknown fields are ignored, None values keep the old value)."""
    p = registry_path(path)
    try:
        with _Lock(p):
            rows = read(p)
            row = next((r for r in rows if r['run'] == run), None)
            if row is None:
                row = {c: '' for c in COLUMNS}
                row['run'] = run
                rows.append(row)
            for k, v in fields.items():
                if k in COLUMNS and v is not None:
                    row[k] = str(v)
            _write(rows, p)
        return True
    except Exception as e:                       # never stop a run because of the registry
        print(f"[run_registry] could not update {p}: {type(e).__name__}: {e}", file=sys.stderr)
        return False


def manifest_dir(run, path=None):
    """Folder of the run's manifest: <simulations>/<run>_meta (next to nb7's run_meta.json). With a custom
    registry file (tests, other projects) it is <registry folder>/<run>_meta instead, so nothing is written
    in the project's simulations folder."""
    import project_paths as PP
    P = PP.get()
    if path is None or os.path.abspath(path) == os.path.abspath(P.run_registry):
        return os.path.join(P.simulations, run + '_meta')
    return os.path.join(os.path.dirname(os.path.abspath(path)), run + '_meta')


class Tracker:
    def __init__(self, run, path=None, manifest=True, **fields):
        self.run, self.path, self.done = run, path, False
        self.t0 = time.time()
        self.meta_dir = None
        fields.setdefault('status', 'started')
        fields.setdefault('started', time.strftime(_NOW))
        fields.setdefault('git_commit', _git_commit())
        fields.setdefault('finished', '')
        if manifest:
            fields.update(self._write_manifest(fields))
        upsert(run, path, **fields)
        atexit.register(self._at_exit)

    def _write_manifest(self, fields):
        """Write <run>_meta/manifest.json (code versions, inputs, environment); returns registry columns."""
        try:
            import run_manifest as RM
            self.meta_dir = manifest_dir(self.run, self.path)
            f = RM.write(self.run, self.meta_dir, run_dir=fields.get('folder'), fields=fields)
            if not f:
                return {}
            m = RM.read(self.meta_dir)
            return {'git_branch': (m.get('cathy_da_ert') or {}).get('branch') or 'detached',
                    'pycathy_version': (m.get('pycathy') or {}).get('version'), 'manifest': f}
        except Exception as e:                   # never stop a run because of the manifest
            print(f"[run_registry] manifest not written: {type(e).__name__}: {e}", file=sys.stderr)
            return {}

    def finish(self, status='finished', **extra):
        self.done = True
        finished, dur = time.strftime(_NOW), f"{(time.time() - self.t0) / 60:.1f}"
        if self.meta_dir:
            import run_manifest as RM
            RM.update(self.meta_dir, status=status, finished=finished, duration_min=dur)
        upsert(self.run, self.path, status=status, finished=finished, duration_min=dur, **extra)

    def _at_exit(self):
        if not self.done:                        # crashed / killed by an exception without finish()
            self.finish('failed')


def track(run, path=None, **fields):
    return Tracker(run, path, **fields)


# ---------------------------------------------------------------------------------------------------------
# backfill: register runs that already exist on disk
# ---------------------------------------------------------------------------------------------------------
_REAL = re.compile(r'^NC_ERT_(?P<m>PF_|ESMDA_)?(?P<w>withDA|withoutDA)_(?P<sc>.+?)(?P<mat>_mat)?$')


def backfill(sim_dir=None, path=None, logs_dir=None):
    """Register every run found in sim_dir: synthetic ones from <run>_meta/run_meta.json, real ones from their
    folder name (+ _finished.txt / _esmda_result.pkl). Existing rows are only completed, never downgraded.
    Returns the number of rows written."""
    import project_paths as PP
    P = PP.get()
    sim_dir = sim_dir or P.simulations
    logs_dir = logs_dir or P.logs
    if not os.path.isdir(sim_dir):
        print(f"[run_registry] no simulations folder: {sim_dir}")
        return 0
    have = {r['run']: r for r in read(path)}
    n = 0
    for d in sorted(os.listdir(sim_dir)):
        full = os.path.join(sim_dir, d)
        if d.endswith('_meta') and os.path.isfile(os.path.join(full, 'run_meta.json')):
            try:
                with open(os.path.join(full, 'run_meta.json')) as f:
                    m = json.load(f)
            except (OSError, ValueError):
                continue
            run = m.get('run') or d[:-5]
            if run in have and have[run]['status'] in ('finished', 'failed'):
                continue
            log = os.path.join(logs_dir, run + '.log')
            upsert(run, path, kind='synthetic', case=m.get('case'), scenario=m.get('scenario'),
                   da_type=m.get('da_type'), nens=m.get('nens'), seed=m.get('seed'), prior_bias=m.get('prior_bias'),
                   with_da=m.get('with_da'), status=m.get('status', 'unknown'), started=m.get('started'),
                   finished=m.get('finished'), folder=os.path.join(sim_dir, run),
                   log=log if os.path.isfile(log) else None)
            n += 1
        elif _REAL.match(d) and os.path.isdir(full):
            g = _REAL.match(d).groupdict()
            da = {'PF_': 'pf', 'ESMDA_': 'esmda', None: 'enkf'}[g['m']]
            ok = os.path.isfile(os.path.join(full, d + '_finished.txt')) or os.path.isfile(
                os.path.join(full, d + '_esmda_result.pkl'))
            if d in have and have[d]['status'] in ('finished', 'failed'):
                continue
            upsert(d, path, kind='real', case='NC', scenario=g['sc'], da_type=da, with_da=g['w'] == 'withDA',
                   status='finished' if ok else 'incomplete', folder=full)
            n += 1
    return n


# ---------------------------------------------------------------------------------------------------------
def _table(rows, cols):
    w = {c: max(len(c), *(len(str(r.get(c, ''))) for r in rows)) for c in cols} if rows else {c: len(c) for c in cols}
    print('  '.join(c.ljust(w[c]) for c in cols))
    for r in rows:
        print('  '.join(str(r.get(c, '')).ljust(w[c]) for c in cols))


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['list', 'summary', 'backfill'])
    ap.add_argument('--registry', default=None, help='csv (default: [outputs].run_registry of data_path.toml)')
    ap.add_argument('--status', default=None)
    ap.add_argument('--case', default=None)
    ap.add_argument('--kind', default=None, choices=['synthetic', 'real'])
    ap.add_argument('--sim-dir', default=None)
    a = ap.parse_args(argv)
    if a.cmd == 'backfill':
        n = backfill(a.sim_dir, a.registry)
        print(f"{n} run(s) registered in {registry_path(a.registry)}")
        return 0
    rows = read(a.registry)
    for k in ('status', 'case', 'kind'):
        v = getattr(a, k)
        if v:
            rows = [r for r in rows if r[k] == v]
    if a.cmd == 'list':
        _table(rows, ['run', 'status', 'started', 'duration_min', 'git_commit'])
        print(f"{len(rows)} run(s) in {registry_path(a.registry)}")
    else:
        from collections import Counter
        for key in ('status', 'case', 'da_type', 'kind'):
            print(f"{key:8s}", dict(Counter(r[key] or '-' for r in rows)))
        print(f"total    {len(rows)}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
