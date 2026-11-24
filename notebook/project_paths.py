"""Single access point to every path of the project: reads data_path.toml.

    import project_paths as PP
    P = PP.get()                      # cached; P.simulations, P.cache, P.synthetic_dir, ...
    P.root_str                        # project root WITH trailing slash (legacy `ROOT + 'data/...'` code)

Lookup order of the toml: $DATA_PATH_TOML, <project>/data_path.toml (one level above this file),
<this folder>/data_path.toml. Without any toml the defaults reproduce the historical layout (root = ../),
so nothing breaks if the file is missing. Every path returned is ABSOLUTE (CATHY changes the cwd).

    python project_paths.py            # print every resolved path and whether it exists
    python project_paths.py --get simulations_str     # one value (used by run_twin.sh)
"""
import os
import sys
import tomllib

_HERE = os.path.dirname(os.path.abspath(__file__))

# key -> default (relative to the project root); "" = derived (see _derive)
DEFAULTS = {
    'inputs': {
        'data': 'data', 'data_3d': 'data/3d', 'site_params': 'data/site_params.json',
        'ert_schedule': 'data/3d/schedule2023_3d.csv', 'sequence_dir': 'data/3d',
        'mesh_dir': 'data/3d/mesh_Luca_Agripolis3d', 'mesh_forward': 'meshForward.vtk',
        'sensors_dir': 'data/sensors',
    },
    'scenarios': {'real': 'notebook/scenario_nc_ERT.toml', 'synthetic': 'notebook/scenarios_twin.toml'},
    'outputs': {
        'synthetic_dir': 'data/3d/synthetic', 'simulations': 'simulations', 'figures_real': '',
        'figures_synth': '', 'cache': 'notebook/cache', 'logs': 'notebook/logs',
        'run_registry': 'registry/runs_registry.csv',
    },
    'code': {'notebook': 'notebook', 'pycathy': '', 'cathy': ''},
}


def _expand(p):
    return os.path.expanduser(os.path.expandvars(str(p)))


def find_toml(explicit=None):
    """Path of the data_path.toml to use, or None."""
    chosen = explicit or os.environ.get('DATA_PATH_TOML')
    if chosen:                                   # an explicit choice must exist: no silent fallback
        if not os.path.isfile(_expand(chosen)):
            raise FileNotFoundError(f"data_path.toml not found: {chosen}")
        return os.path.abspath(_expand(chosen))
    for c in (os.path.join(os.path.dirname(_HERE), 'data_path.toml'), os.path.join(_HERE, 'data_path.toml')):
        if os.path.isfile(c):
            return os.path.abspath(c)
    return None


class ProjectPaths:
    """Resolved (absolute) paths. Attributes: root, data, data_3d, site_params, ert_schedule, sequence_dir,
    mesh_dir, mesh_forward, sensors_dir, scenario_real, scenario_synthetic, synthetic_dir, simulations,
    figures_real, figures_synth, cache, logs, run_registry, notebook, pycathy, cathy ('' if not set)."""

    def __init__(self, toml_file=None):
        self.toml_file = find_toml(toml_file)
        raw = {}
        if self.toml_file:
            with open(self.toml_file, 'rb') as f:
                raw = tomllib.load(f)
        base = os.path.dirname(self.toml_file) if self.toml_file else os.path.dirname(_HERE)
        root = _expand(raw.get('project', {}).get('root', '.' if self.toml_file else '.'))
        self.root = os.path.abspath(root if os.path.isabs(root) else os.path.join(base, root))

        def val(sec, key):
            v = raw.get(sec, {}).get(key, DEFAULTS[sec][key])
            unknown = set(raw.get(sec, {})) - set(DEFAULTS[sec])
            assert not unknown, f"data_path.toml [{sec}]: unknown key(s) {sorted(unknown)} (valid: {sorted(DEFAULTS[sec])})"
            return _expand(v)

        def res(p):
            return os.path.abspath(p if os.path.isabs(p) else os.path.join(self.root, p))

        for k in DEFAULTS['inputs']:
            if k != 'mesh_forward':
                setattr(self, k, res(val('inputs', k)))
        mf = val('inputs', 'mesh_forward')
        self.mesh_forward = mf if os.path.isabs(mf) else os.path.join(self.mesh_dir, mf)
        self.scenario_real = res(val('scenarios', 'real'))
        self.scenario_synthetic = res(val('scenarios', 'synthetic'))
        for k in ('synthetic_dir', 'simulations', 'cache', 'logs', 'run_registry'):
            setattr(self, k, res(val('outputs', k)))
        self.figures_real = res(val('outputs', 'figures_real') or os.path.join(self.simulations, 'figures_nb5'))
        self.figures_synth = res(val('outputs', 'figures_synth') or os.path.join(self.simulations, 'figures_nb8'))
        self.notebook = res(val('code', 'notebook'))
        py = val('code', 'pycathy')
        self.pycathy = res(py) if py else ''
        ca = val('code', 'cathy')
        self.cathy = res(ca) if ca else ''

    # --- legacy string forms: `ROOT + 'data/...'` and `path2prj + name` need a trailing slash -------------
    @property
    def root_str(self):
        return self.root + os.sep

    @property
    def simulations_str(self):
        return self.simulations + os.sep

    def is_root(self, root):
        """True when `root` (a '../'-style argument) is this project's root: then the toml decides the folders."""
        return os.path.abspath(_expand(root)) == self.root

    def synthetic_case_dir(self, case):
        return os.path.join(self.synthetic_dir, case)

    def ensure_dirs(self):
        for d in (self.simulations, self.cache, self.logs, os.path.dirname(self.run_registry)):
            os.makedirs(d, exist_ok=True)

    def setup_sys_path(self):
        """Put a pycathy clone (code.pycathy) in front of sys.path; no-op if unset."""
        if self.pycathy and self.pycathy not in sys.path:
            assert os.path.isdir(os.path.join(self.pycathy, 'pyCATHY')), f"code.pycathy: no pyCATHY/ in {self.pycathy}"
            sys.path.insert(0, self.pycathy)

    def as_dict(self):
        return {k: v for k, v in vars(self).items() if not k.startswith('_')}

    def report(self):
        print(f"data_path.toml: {self.toml_file or '(none found: built-in defaults)'}")
        for k, v in self.as_dict().items():
            if k == 'toml_file':
                continue
            flag = '' if not v else ('   [exists]' if os.path.exists(v) else '   [missing]')
            print(f"  {k:20s} {v}{flag}")


_CACHE = {}


def get(toml_file=None, reload=False):
    """Cached ProjectPaths (pass toml_file / reload=True to switch file, e.g. in tests)."""
    key = toml_file or os.environ.get('DATA_PATH_TOML') or ''
    if reload or key not in _CACHE:
        _CACHE[key] = ProjectPaths(toml_file)
    return _CACHE[key]


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--toml', default=None)
    ap.add_argument('--get', default=None, help='print one attribute (root_str, simulations_str, pycathy, ...)')
    a = ap.parse_args()
    P = get(a.toml)
    if a.get:
        print(getattr(P, a.get))
    else:
        P.report()
