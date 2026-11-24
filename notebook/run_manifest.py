"""Reproducibility manifest of a run: which code, which versions, which inputs.

One `manifest.json` per run, written next to run_meta.json in <simulations>/<run>_meta/ when the run starts
(run_registry.track does it, so nb5 / nb6 / nb7 need no change) and completed when it finishes.

What it records
  * this repo (CATHY DA ERT): full commit, branch, upstream, remote url, `git describe`, dirty flag; when the
    working tree is dirty the uncommitted diff is saved as <run>_meta/uncommitted.patch (+ its sha256)
  * pyCATHY: version, where it is imported from, and -- when that is a git clone (editable install or
    [code].pycathy) -- commit, branch, remote, dirty flag (+ <run>_meta/pycathy_uncommitted.patch)
  * CATHY (the Fortran model): git info of its source folder if found, sha256 / size / mtime of the executable
    found in the run folder or on the search path
  * python, platform, every installed package version (pip freeze like), conda / venv, thread env vars
  * inputs: sha256 of the scenario files, data_path.toml, site_params, ERT schedule; name+size digest of the
    big input folders (cheap, no full read); copies of the small config files in <run>_meta/
  * the command line, the run fields given to the registry (seed, nens, da_type, ...), timestamps

Never stops a run: any failure is reported in one line on stderr and the manifest keeps what it could collect.

    import run_manifest as RM
    RM.write(run_name, meta_dir, fields={...})            # at start   (done by run_registry.track)
    RM.update(meta_dir, status='finished', ...)           # at the end (done by Tracker.finish)

Command line:
    python run_manifest.py collect                        # print the manifest of the current environment
    python run_manifest.py show RUN [--sim-dir DIR]       # one-screen summary of a run's manifest
    python run_manifest.py compare RUN_A RUN_B            # what differs (versions, commits, inputs)
"""
import hashlib
import json
import os
import platform
import re
import socket
import subprocess
import sys
import time

SCHEMA = 1
MANIFEST_NAME = 'manifest.json'
_NOW = '%Y-%m-%d %H:%M:%S'
_MAX_PATCH = 2_000_000          # do not store diffs larger than this (bytes); the sha256 is kept anyway
_MAX_HASH = 200_000_000         # files above this are fingerprinted by name+size+mtime only
_ENV_VARS = ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS',
             'PYTHONHASHSEED', 'DATA_PATH_TOML', 'CATHY_SRC', 'CATHY_EXE', 'CONDA_DEFAULT_ENV', 'VIRTUAL_ENV')
_KEY_PACKAGES = ('pyCATHY', 'numpy', 'scipy', 'pandas', 'xarray', 'netCDF4', 'pygimli', 'pyvista', 'vtk',
                 'matplotlib', 'scikit-learn', 'numba', 'tomli', 'pyarrow')


# ---------------------------------------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------------------------------------
def _run(cmd, cwd=None, timeout=15, strip=True):
    try:
        out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        if out.returncode != 0:
            return None
        return out.stdout.strip() if strip else out.stdout
    except Exception:
        return None


def _git(cwd, *args, strip=True):
    return _run(['git', *args], cwd=cwd, strip=strip)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def _sha256_text(text):
    return hashlib.sha256(text.encode('utf-8', 'replace')).hexdigest()


def _strip_credentials(url):
    """https://user:token@github.com/x/y.git -> https://github.com/x/y.git"""
    return re.sub(r'//[^/@]+@', '//', url) if url else url


def _stamp(path):
    st = os.stat(path)
    return {'size': st.st_size, 'mtime': time.strftime(_NOW, time.localtime(st.st_mtime))}


# ---------------------------------------------------------------------------------------------------------
# git
# ---------------------------------------------------------------------------------------------------------
def git_info(path, patch_file=None):
    """Git state of the repository containing `path` ({'git': False} when it is not one).
    A dirty tree writes its diff (tracked files vs HEAD) to patch_file."""
    if not path or not os.path.exists(path):
        return {'git': False}
    cwd = path if os.path.isdir(path) else os.path.dirname(path)
    top = _git(cwd, 'rev-parse', '--show-toplevel')
    if not top:
        return {'git': False, 'path': os.path.abspath(path)}
    commit = _git(top, 'rev-parse', 'HEAD')
    branch = _git(top, 'rev-parse', '--abbrev-ref', 'HEAD')
    info = {
        'git': True, 'toplevel': top, 'commit': commit,
        'commit_short': (commit or '')[:8],
        'branch': branch if branch != 'HEAD' else None,            # None = detached HEAD
        'detached': branch == 'HEAD',
        'describe': _git(top, 'describe', '--tags', '--always', '--dirty'),
        'upstream': _git(top, 'rev-parse', '--abbrev-ref', '--symbolic-full-name', '@{u}'),
        'remote_url': _strip_credentials(_git(top, 'remote', 'get-url', 'origin')),
        'commit_date': _git(top, 'log', '-1', '--format=%cI'),
        'commit_subject': _git(top, 'log', '-1', '--format=%s'),
    }
    if info['detached']:                                           # which branches contain this commit
        cont = _git(top, 'branch', '-a', '--contains', 'HEAD')
        info['branches_containing_head'] = [b.strip('* ').strip() for b in (cont or '').splitlines()][:10]
    status = _git(top, 'status', '--porcelain', '--untracked-files=no', strip=False)   # keep the leading ' M'
    changed = [ln[3:] for ln in (status or '').splitlines()]
    info['dirty'] = bool(changed)
    info['dirty_files'] = changed[:50]
    untracked = _git(top, 'ls-files', '--others', '--exclude-standard')
    info['untracked_files'] = (untracked or '').splitlines()[:50]
    if changed:
        diff = _git(top, 'diff', 'HEAD') or ''
        info['diff_sha256'] = _sha256_text(diff)
        if patch_file and diff and len(diff) <= _MAX_PATCH:
            os.makedirs(os.path.dirname(patch_file), exist_ok=True)
            with open(patch_file, 'w') as f:
                f.write(diff + '\n')
            info['diff_patch'] = os.path.basename(patch_file)
        elif diff:
            info['diff_patch'] = None                              # too large to store: only the hash
    return info


# ---------------------------------------------------------------------------------------------------------
# packages: pyCATHY, CATHY, python environment
# ---------------------------------------------------------------------------------------------------------
def _installed_packages():
    from importlib import metadata
    pk = {}
    for d in metadata.distributions():
        name = d.metadata['Name'] if d.metadata else None
        if name:
            pk[name] = d.version
    return dict(sorted(pk.items(), key=lambda kv: kv[0].lower()))


def _find_package_dir(name):
    """Folder of a top-level package WITHOUT importing it (respects the current sys.path order)."""
    import importlib.util
    try:
        spec = importlib.util.find_spec(name)
    except Exception:
        return None
    if spec is None:
        return None
    if spec.submodule_search_locations:
        return os.path.abspath(list(spec.submodule_search_locations)[0])
    return os.path.dirname(os.path.abspath(spec.origin)) if spec.origin else None


def _version_from_source(pkg_dir):
    """`__version__ = '...'` of the package, parsed from source (an editable install's metadata can be stale)."""
    for fn in ('_version.py', '__init__.py', 'version.py'):
        p = os.path.join(pkg_dir, fn)
        if os.path.isfile(p):
            try:
                with open(p, errors='replace') as f:
                    m = re.search(r'''__version__\s*=\s*['"]([^'"]+)['"]''', f.read())
                if m:
                    return m.group(1)
            except OSError:
                pass
    return None


def pycathy_info(meta_dir=None, packages=None):
    pkg_dir = _find_package_dir('pyCATHY')
    if not pkg_dir:
        return {'found': False}
    info = {'found': True, 'package_dir': pkg_dir,
            'version_metadata': (packages or {}).get('pyCATHY') or (packages or {}).get('pycathy'),
            'version_source': _version_from_source(pkg_dir)}
    info['version'] = info['version_metadata'] or info['version_source']
    patch = os.path.join(meta_dir, 'pycathy_uncommitted.patch') if meta_dir else None
    info['repo'] = git_info(pkg_dir, patch)
    return info


def _sha_exe(path):
    d = {'path': path}
    try:
        d.update(_stamp(path))
        d['sha256'] = _sha256(path) if d['size'] <= _MAX_HASH else None
    except OSError as e:
        d['error'] = str(e)
    return d


def cathy_info(run_dir=None, pycathy_dir=None, cathy_src=None, meta_dir=None):
    """The Fortran CATHY model: source folder (git state) and executable (sha256). Best effort: the source
    folder is looked up in $CATHY_SRC, [code].cathy, <pyCATHY>/CATHY|src and its parent; the executable in
    $CATHY_EXE, the run folder and the source folder."""
    cands = [cathy_src, os.environ.get('CATHY_SRC')]
    if pycathy_dir:
        base = os.path.dirname(pycathy_dir)
        cands += [os.path.join(pycathy_dir, 'CATHY'), os.path.join(pycathy_dir, 'src'),
                  os.path.join(base, 'CATHY'), os.path.join(base, 'cathy')]
    src = next((os.path.abspath(os.path.expanduser(c)) for c in cands if c and os.path.isdir(os.path.expanduser(c))),
               None)
    info = {'source_dir': src, 'source_repo': None, 'executable': None}
    if src:
        info['source_repo'] = git_info(src, os.path.join(meta_dir, 'cathy_uncommitted.patch') if meta_dir else None)
    exe = os.environ.get('CATHY_EXE')
    if not (exe and os.path.isfile(os.path.expanduser(exe))):
        exe = None
        for base in (run_dir, src):
            if base and os.path.isdir(base):
                for root, dirs, files in os.walk(base):
                    depth = root[len(base):].count(os.sep)
                    if depth > 2:
                        dirs[:] = []
                        continue
                    hit = [f for f in files if re.fullmatch(r'(?i)(cathy|cathy_main|cathy\.exe)', f)
                           and os.access(os.path.join(root, f), os.X_OK)]
                    if hit:
                        exe = os.path.join(root, hit[0])
                        break
            if exe:
                break
    if exe:
        info['executable'] = _sha_exe(os.path.abspath(os.path.expanduser(exe)))
    return info


# ---------------------------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------------------------
def _tree_digest(folder, max_files=200_000):
    """name+size digest of a folder (no content read): detects added / removed / resized files."""
    h, n, total = hashlib.sha256(), 0, 0
    for root, dirs, files in os.walk(folder):
        dirs.sort()
        for f in sorted(files):
            p = os.path.join(root, f)
            try:
                s = os.path.getsize(p)
            except OSError:
                continue
            h.update(f"{os.path.relpath(p, folder)}\0{s}\n".encode())
            n += 1
            total += s
            if n >= max_files:
                return {'n_files': n, 'bytes': total, 'digest': h.hexdigest(), 'truncated': True}
    return {'n_files': n, 'bytes': total, 'digest': h.hexdigest()}


def inputs_info(P, meta_dir=None):
    """Fingerprint the inputs listed in data_path.toml (P = project_paths.ProjectPaths or None)."""
    out = {'files': {}, 'folders': {}, 'copied': []}
    if P is None:
        return out
    files = {'data_path.toml': P.toml_file, 'scenario_real': P.scenario_real,
             'scenario_synthetic': P.scenario_synthetic, 'site_params': P.site_params,
             'ert_schedule': P.ert_schedule, 'mesh_forward': P.mesh_forward}
    small = ('data_path.toml', 'scenario_real', 'scenario_synthetic', 'site_params', 'ert_schedule')
    for k, p in files.items():
        if not p or not os.path.isfile(p):
            out['files'][k] = {'path': p, 'exists': False}
            continue
        d = {'path': p, 'exists': True, **_stamp(p)}
        d['sha256'] = _sha256(p) if d['size'] <= _MAX_HASH else None
        out['files'][k] = d
        if meta_dir and k in small and d['size'] < 1_000_000:          # keep a copy of the small configs
            try:
                import shutil
                os.makedirs(os.path.join(meta_dir, 'inputs'), exist_ok=True)
                shutil.copy2(p, os.path.join(meta_dir, 'inputs', os.path.basename(p)))
                out['copied'].append(os.path.join('inputs', os.path.basename(p)))
            except OSError:
                pass
    for k in ('data_3d', 'sensors_dir', 'sequence_dir', 'mesh_dir'):
        p = getattr(P, k, None)
        if p and os.path.isdir(p):
            out['folders'][k] = {'path': p, **_tree_digest(p)}
    return out


# ---------------------------------------------------------------------------------------------------------
# collect / write / update
# ---------------------------------------------------------------------------------------------------------
def collect(run=None, meta_dir=None, run_dir=None, fields=None, P=None, repo_dir=None):
    """Build the manifest dict. Every section is independent: one failing never loses the others."""
    here = repo_dir or os.path.dirname(os.path.abspath(__file__))
    m = {'schema': SCHEMA, 'run': run, 'created': time.strftime(_NOW), 'status': 'started'}
    errors = {}

    def section(name, fn):
        try:
            m[name] = fn()
        except Exception as e:                                          # noqa: BLE001 - never stop a run
            errors[name] = f"{type(e).__name__}: {e}"

    section('cathy_da_ert', lambda: git_info(here, os.path.join(meta_dir, 'uncommitted.patch') if meta_dir else None))
    pk = {}
    try:
        pk = _installed_packages()
    except Exception as e:                                              # noqa: BLE001
        errors['packages'] = f"{type(e).__name__}: {e}"
    m['packages'] = pk
    m['key_packages'] = {k: next((v for n, v in pk.items() if n.lower() == k.lower()), None) for k in _KEY_PACKAGES}
    section('pycathy', lambda: pycathy_info(meta_dir, pk))
    pyc_dir = (m.get('pycathy') or {}).get('package_dir')
    section('cathy', lambda: cathy_info(run_dir, pyc_dir, getattr(P, 'cathy', '') or None, meta_dir))
    m['python'] = {'version': sys.version.split()[0], 'executable': sys.executable, 'implementation':
                   platform.python_implementation()}
    m['platform'] = {'system': platform.system(), 'release': platform.release(), 'machine': platform.machine(),
                     'node': socket.gethostname(), 'cpu_count': os.cpu_count()}
    m['environment'] = {k: os.environ[k] for k in _ENV_VARS if k in os.environ}
    m['command'] = {'argv': list(sys.argv), 'cwd': os.getcwd()}
    m['fields'] = {k: v for k, v in (fields or {}).items() if v is not None}
    section('inputs', lambda: inputs_info(P, meta_dir))
    if errors:
        m['collect_errors'] = errors
    return m


def _dump(m, meta_dir):
    os.makedirs(meta_dir, exist_ok=True)
    final = os.path.join(meta_dir, MANIFEST_NAME)
    tmp = final + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(m, f, indent=2, default=str)
    os.replace(tmp, final)
    return final


def write(run, meta_dir, run_dir=None, fields=None, P=None):
    """Collect and write <meta_dir>/manifest.json. Returns its path, or None (and one stderr line) on failure."""
    try:
        if P is None:
            try:
                import project_paths as PP
                P = PP.get()
            except Exception:                                           # noqa: BLE001
                P = None
        return _dump(collect(run, meta_dir, run_dir, fields, P), meta_dir)
    except Exception as e:                                              # noqa: BLE001
        print(f"[run_manifest] could not write the manifest in {meta_dir}: {type(e).__name__}: {e}", file=sys.stderr)
        return None


def update(meta_dir, **fields):
    """Merge fields (status, finished, duration_min, ...) into an existing manifest; refresh the CATHY
    executable record (it only exists once the run folder is built)."""
    p = os.path.join(meta_dir, MANIFEST_NAME)
    try:
        with open(p) as f:
            m = json.load(f)
        run_dir = (m.get('fields') or {}).get('folder')
        if run_dir and not (m.get('cathy') or {}).get('executable'):
            pyc_dir = (m.get('pycathy') or {}).get('package_dir')
            try:
                m['cathy'] = cathy_info(run_dir, pyc_dir, (m.get('cathy') or {}).get('source_dir'), meta_dir)
            except Exception:                                           # noqa: BLE001
                pass
        m.update({k: v for k, v in fields.items() if v is not None})
        _dump(m, meta_dir)
        return True
    except Exception as e:                                              # noqa: BLE001
        print(f"[run_manifest] could not update {p}: {type(e).__name__}: {e}", file=sys.stderr)
        return False


def read(meta_dir):
    with open(os.path.join(meta_dir, MANIFEST_NAME)) as f:
        return json.load(f)


# ---------------------------------------------------------------------------------------------------------
# command line: collect / show / compare
# ---------------------------------------------------------------------------------------------------------
def _repo_line(label, r):
    if not r or not r.get('git'):
        return f"  {label:14s} (not a git repository)"
    where = r.get('branch') or f"detached HEAD ({', '.join(r.get('branches_containing_head') or []) or '?'})"
    dirty = f"  DIRTY: {len(r.get('dirty_files', []))} file(s)" if r.get('dirty') else ''
    return f"  {label:14s} {r.get('commit_short')}  {where}{dirty}\n  {'':14s} {r.get('remote_url') or ''}"


def summary(m):
    lines = [f"run      {m.get('run')}   status={m.get('status')}   created={m.get('created')}"]
    lines.append(_repo_line('CATHY DA ERT', m.get('cathy_da_ert')))
    pc = m.get('pycathy') or {}
    if pc.get('found'):
        lines.append(f"  {'pyCATHY':14s} version {pc.get('version')}  from {pc.get('package_dir')}")
        lines.append(_repo_line('', pc.get('repo')))
    else:
        lines.append(f"  {'pyCATHY':14s} NOT FOUND on sys.path")
    ct = m.get('cathy') or {}
    ex = ct.get('executable') or {}
    lines.append(f"  {'CATHY':14s} exe sha256 {(ex.get('sha256') or '?')[:16]}  {ex.get('path') or '(executable not found)'}")
    lines.append(_repo_line('', ct.get('source_repo')))
    lines.append(f"  {'python':14s} {m.get('python', {}).get('version')}  on {m.get('platform', {}).get('node')}")
    kp = {k: v for k, v in (m.get('key_packages') or {}).items() if v}
    lines.append(f"  {'packages':14s} " + ', '.join(f"{k}={v}" for k, v in kp.items()))
    lines.append(f"  {'command':14s} {' '.join(m.get('command', {}).get('argv', []))}")
    lines.append(f"  {'fields':14s} {m.get('fields')}")
    return '\n'.join(lines)


def _flat(d, prefix=''):
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flat(v, key + '.'))
        else:
            out[key] = v
    return out


_IGNORE_IN_COMPARE = ('created', 'finished', 'duration_min', 'command.cwd', 'platform.node', 'status')


def compare(a, b):
    """{key: (value_a, value_b)} of what differs between two manifests (timestamps / paths of the machine ignored)."""
    fa, fb = _flat(a), _flat(b)
    diff = {}
    for k in sorted(set(fa) | set(fb)):
        if k in _IGNORE_IN_COMPARE or k.endswith(('.mtime', '.toplevel')) or k.startswith(('run', 'fields.folder',
                                                                                          'fields.log')):
            continue
        if fa.get(k) != fb.get(k):
            diff[k] = (fa.get(k), fb.get(k))
    return diff


def _meta_dir(run, sim_dir=None):
    if sim_dir is None:
        import project_paths as PP
        sim_dir = PP.get().simulations
    return os.path.join(sim_dir, run + '_meta')


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['collect', 'show', 'compare'])
    ap.add_argument('runs', nargs='*')
    ap.add_argument('--sim-dir', default=None)
    ap.add_argument('--json', action='store_true', help='print the full json instead of the summary')
    a = ap.parse_args(argv)
    if a.cmd == 'collect':
        try:
            import project_paths as PP
            P = PP.get()
        except Exception:                                               # noqa: BLE001
            P = None
        m = collect('(current environment)', None, None, None, P)
        print(json.dumps(m, indent=2, default=str) if a.json else summary(m))
        return 0
    need = 1 if a.cmd == 'show' else 2
    if len(a.runs) != need:
        ap.error(f"'{a.cmd}' needs {need} run name(s)")
    ms = []
    for r in a.runs:
        p = _meta_dir(r, a.sim_dir)
        try:
            ms.append(read(p))
        except OSError:
            print(f"no manifest for run '{r}' in {p}", file=sys.stderr)
            return 1
    if a.cmd == 'show':
        print(json.dumps(ms[0], indent=2, default=str) if a.json else summary(ms[0]))
        return 0
    diff = compare(*ms)
    for k, (x, y) in diff.items():
        if k.startswith('packages.') and not a.json:
            continue                                                    # listed below, once
        print(f"{k}\n    {a.runs[0]}: {x}\n    {a.runs[1]}: {y}")
    pk = [k for k in diff if k.startswith('packages.')]
    if pk:
        print(f"{len(pk)} installed package version(s) differ: " + ', '.join(k[9:] for k in pk[:15])
              + (' ...' if len(pk) > 15 else ''))
    print('identical (code, versions, inputs)' if not diff else f"{len(diff)} difference(s)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
