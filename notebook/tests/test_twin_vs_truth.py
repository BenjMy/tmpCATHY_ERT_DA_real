"""Twin set-up vs truth (todo.md 5.3).

The truth (synthetic_truth_ert.py) and the DA model (da_ert_synthetic.py) both call
`model_config.setup_model(simu, case, ic_values=, extra_times=, soil_overrides=, delta=)`. Whatever differs between the
two CATHY projects therefore comes from these ARGUMENTS (or from code outside setup_model). These tests run
setup_model twice on a recording stand-in of the CATHY object (`FakeCathy`) with the arguments of the truth and of the
DA run, and compare what was written: mesh inputs, soil table, ic, atmbc, parm, boundary conditions.

No CATHY and no pyCATHY needed: `pyCATHY.meshtools` is replaced (for these tests only) by a fake with the same
contract (layer depths; mid-depth rule of map_param_to_layers, todo 5.1 (a)). What this does NOT check (needs the real
files / a run) is listed as skipped tests at the end, one per open item of todo.md 5.3.
"""
import importlib
import importlib.util
import sys
import types

import numpy as np
import pandas as pd
import pytest

CASES = ['NC', 'HC']
ERT_TIMES = [3 * 3600 + 1800.0, 3 * 3600 + 7200.0, 30 * 3600.0]       # a few surveys [s], off the hourly grid


# ------------------------------------------------------------------ stand-ins
def _fake_meshtools():
    mt = types.ModuleType('pyCATHY.meshtools')

    def get_layer_depths(dem_parameters):
        zr = np.asarray(dem_parameters['zratio'], float) * dem_parameters['base']
        bot = np.cumsum(zr)
        return bot - zr, bot

    def map_param_to_layers(ltop, lbot, ranges, values):
        mid = 0.5 * (np.asarray(ltop) + np.asarray(lbot))
        out = []
        for m in mid:
            g = next((i for i, (a, b) in enumerate(ranges) if a <= m < b), len(ranges) - 1)
            out.append(float(np.ravel(values[g])[0]))
        return np.array(out)

    mt.get_layer_depths, mt.map_param_to_layers = get_layer_depths, map_param_to_layers
    return mt


class FakeCathy:
    """Records what setup_model writes. Same method names / arguments as the CATHY object used by setup_model."""

    def __init__(self):
        self.rec = {}
        self.dem_parameters = None

    def update_prepo_inputs(self, DEM, delta_x, delta_y, dr, zratio, base):
        self.dem_parameters = dict(zratio=list(zratio), base=base)
        self.rec['prepo'] = dict(DEM=np.array(DEM), delta=(delta_x, delta_y, dr), zratio=list(zratio), base=base)

    def update_veg_map(self, m):
        self.rec['veg'] = np.array(m)

    def create_mesh_vtk(self, verbose=False):
        self.rec['mesh'] = True

    def update_zone(self, z):
        self.rec['zones'] = np.array(z)

    def update_soil(self, **kw):
        if 'SPP_map' in kw:
            self.rec['soil_df'] = kw['SPP_map'].copy()
            self.rec['soil_kw'] = {k: v for k, v in kw.items() if k != 'SPP_map'}

    def init_soil_SPP_map_df(self, nzones, nstr):
        from model_config import SOIL_COLS
        idx = pd.MultiIndex.from_product([range(1, nzones + 1), range(1, nstr + 1)], names=['zone', 'layer'])
        return pd.DataFrame(np.nan, index=idx, columns=SOIL_COLS)

    def map_prop_2mesh_markers(self, name, arr, to_nodes=True):
        n_lay, n = len(arr), self._nodes_per_level
        lev = np.minimum(np.arange(n_lay + 1), n_lay - 1)           # node level k takes the value of layer min(k, last)
        return np.repeat(np.asarray(arr, float)[lev], n)

    def update_ic(self, **kw):
        self.rec['ic'] = np.array(kw['pressure_head_ini'], float)
        self.rec['ic_kw'] = {k: v for k, v in kw.items() if k != 'pressure_head_ini'}

    def update_atmbc(self, **kw):
        self.rec['atmbc'] = dict(kw, netValue=[np.array(a) for a in kw['netValue']])

    def update_parm(self, **kw):
        self.rec['parm'] = kw

    def update_nansfneubc(self, **kw):
        self.rec['nansfneubc'] = kw

    def update_nansfdirbc(self, **kw):
        self.rec['nansfdirbc'] = kw

    def update_sfbc(self, **kw):
        self.rec['sfbc'] = kw


@pytest.fixture(scope='module')
def MC():
    """model_config (and make_scenarios) imported with a fake pyCATHY.meshtools; the real one is left untouched."""
    mp = pytest.MonkeyPatch()
    mods = ('model_config', 'nc_model', 'make_scenarios')
    had = {m: m in sys.modules for m in mods}
    mt = _fake_meshtools()
    if importlib.util.find_spec('pyCATHY') is None:
        pkg = types.ModuleType('pyCATHY')
        pkg.__path__ = []
        mp.setitem(sys.modules, 'pyCATHY', pkg)
    mp.setitem(sys.modules, 'pyCATHY.meshtools', mt)
    mc = importlib.import_module('model_config')
    yield mc
    mp.undo()
    for m in mods:
        if not had[m]:
            sys.modules.pop(m, None)


@pytest.fixture(scope='module')
def MS(MC):
    return importlib.import_module('make_scenarios')


@pytest.fixture
def mt_fake(monkeypatch):
    monkeypatch.setitem(sys.modules, 'pyCATHY.meshtools', _fake_meshtools())


def build(mc, case, **kw):
    """Run setup_model on a recording stand-in; returns (recorded dict, DEM, ic_nodes)."""
    simu = FakeCathy()
    simu._nodes_per_level = (mc.get_case(case, kw.get('delta'))['ncell'] + 1) ** 2
    DEM, ic_nodes = mc.setup_model(simu, case, **kw)
    return simu.rec, DEM, ic_nodes


def truth_build(mc, case, delta=None):
    """Arguments of synthetic_truth_ert.py: no ic_values, no soil_overrides."""
    return build(mc, case, extra_times=ERT_TIMES, delta=delta)


def da_build(mc, case, biased=False, ic_values=None, delta=None):
    """Arguments of da_ert_synthetic.py: ic prior mean per group (default ic_true), biased_overrides when --biased."""
    ic = [[v[0]] for v in mc.get_case(case)['ic_true']] if ic_values is None else ic_values
    return build(mc, case, ic_values=ic, extra_times=ERT_TIMES,
                 soil_overrides=mc.biased_overrides(case) if biased else None, delta=delta)


def assert_same(a, b, skip=()):
    assert set(a) == set(b)
    for k in sorted(set(a) - set(skip)):
        va, vb = a[k], b[k]
        if isinstance(va, pd.DataFrame):
            pd.testing.assert_frame_equal(va, vb, obj=k)
        elif isinstance(va, dict):
            assert set(va) == set(vb), k
            for kk in va:
                _eq(va[kk], vb[kk], f'{k}.{kk}')
        else:
            _eq(va, vb, k)


def _eq(x, y, what):
    if isinstance(x, (list, tuple)) and x and isinstance(x[0], np.ndarray):
        assert len(x) == len(y), what
        for i, (p, q) in enumerate(zip(x, y)):
            assert np.array_equal(p, q), f'{what}[{i}]'
    else:
        assert np.array_equal(np.asarray(x, dtype=object) if isinstance(x, str) else np.asarray(x),
                              np.asarray(y, dtype=object) if isinstance(y, str) else np.asarray(y)), what


# ------------------------------------------------------------------ mesh / geometry / time / boundary conditions
@pytest.mark.parametrize('case', CASES)
def test_unbiased_da_model_is_identical_to_the_truth_model(MC, mt_fake, case):
    """--prior-bias none, ic prior mean = ic_true: every file setup_model writes is the same as for the truth."""
    (t, _, ic_t), (d, _, ic_d) = truth_build(MC, case), da_build(MC, case)
    assert_same(t, d)
    assert np.array_equal(ic_t, ic_d)


@pytest.mark.parametrize('case', CASES)
def test_mesh_layers_and_nodes(MC, mt_fake, case):
    rec, DEM, ic = truth_build(MC, case)
    cfg = MC.get_case(case)
    n = cfg['ncell']
    assert DEM.shape == (n, n) and rec['prepo']['delta'] == (cfg['delta'],) * 3
    zr = np.asarray(rec['prepo']['zratio']) * rec['prepo']['base']
    assert len(zr) == MC.NSTR == 15 and np.isclose(zr.sum(), MC.MAXDEPTH) and np.allclose(zr, MC.ZR_THICK)
    assert ic.size == (n + 1) ** 2 * (MC.NSTR + 1)           # CATHY nodes: (n+1)^2 per level x (layers + 1) levels


def test_delta_override_changes_truth_and_da_together(MC, mt_fake):
    (t, _, _), (d, _, _) = truth_build(MC, 'NC', delta=0.25), da_build(MC, 'NC', delta=0.25)
    assert_same(t, d)
    assert t['prepo']['DEM'].shape == (20, 20)               # 5 m / 0.25 m
    with pytest.raises(AssertionError):
        assert_same(truth_build(MC, 'NC', delta=0.25)[0], da_build(MC, 'NC', delta=0.5)[0])


@pytest.mark.parametrize('case', CASES)
def test_atmbc_and_parm_times(MC, mt_fake, case):
    rec, _, _ = truth_build(MC, case)
    cfg = MC.get_case(case)
    times = np.asarray(rec['atmbc']['time'], float)
    assert np.allclose(times, MC.build_times(ERT_TIMES)) and np.all(np.diff(times) > 0)
    assert all(t in times for t in ERT_TIMES)                 # every survey has an exact CATHY time step
    assert list(rec['parm']['TIMPRTi']) == list(rec['atmbc']['time'])       # output times = atmbc times
    flux = rec['atmbc']['netValue']
    assert len(flux) == len(times)
    irr = [i for i, t in enumerate(times) if t == cfg['irr_start_hour'] * MC.HOUR]
    assert len(irr) == 1                                      # exactly one irrigation step
    y0, y1, x0, x1 = MC.irrigation_window(cfg)
    for i, f in enumerate(flux):
        expect = np.full_like(f, cfg['et'])
        if i == irr[0]:
            expect[y0:y1, x0:x1] = cfg['irr']
        assert np.array_equal(f, expect), f'atmbc step {i}'
    assert rec['atmbc']['IETO'] == 1 and rec['atmbc']['HSPATM'] == 0


def test_boundary_conditions_and_solver_settings_are_fixed(MC, mt_fake):
    rec, _, _ = truth_build(MC, 'NC')
    assert rec['nansfneubc'] == rec['nansfdirbc'] == rec['sfbc'] == dict(no_flow=True)
    for k, v in dict(DTMAX=1e3, DELTAT=1e-1, DTMIN=1e-3, IPRT1=4, VTKF=2, TRAFLAG=0).items():
        assert rec['parm'][k] == v, k
    assert rec['ic_kw'] == dict(INDP=1, IPOND=0)


def test_first_n_surveys_cuts_truth_and_da_at_the_same_survey(MC):
    times_h = {1: 0.5, 2: 2.0, 3: 26.0, 4: 8.0}
    assert MC.first_n_surveys(times_h, 2) == {1: 0.5, 2: 2.0}
    assert MC.first_n_surveys(times_h, None) == times_h and MC.first_n_surveys(times_h, 0) == times_h


# ------------------------------------------------------------------ materials / soil
@pytest.mark.parametrize('case', CASES)
def test_soil_table_rows_follow_the_material_layers(MC, mt_fake, case):
    rec, _, _ = truth_build(MC, case)
    df, cfg = rec['soil_df'], MC.get_case(case)
    assert not df[MC.SOIL_COLS].isna().any().any()
    for (zone, mat), vals in {(z, m): v for z, ms in cfg['soil'].items() for m, v in ms.items()}.items():
        for L in MC.MATERIAL_LAYERS_0[mat]:
            assert df.loc[(zone, L + 1), MC.SOIL_COLS].tolist() == [vals[c] for c in MC.SOIL_COLS], (zone, mat, L)
    assert MC.MATERIAL_LAYERS_0 == {'A': list(range(0, 6)), 'B': [6, 7, 8], 'C': list(range(9, 15))}


def test_hc_zone2_below_material_a_is_zone1_in_truth_and_da(MC, mt_fake):
    """zone2_fill_below=True: zone 2 materials B, C are copies of zone 1 (nb2 never sets them)."""
    for rec in (truth_build(MC, 'HC')[0], da_build(MC, 'HC')[0]):
        df = rec['soil_df']
        for L in MC.MATERIAL_LAYERS_0['B'] + MC.MATERIAL_LAYERS_0['C']:
            assert df.loc[(2, L + 1), MC.SOIL_COLS].tolist() == df.loc[(1, L + 1), MC.SOIL_COLS].tolist()
        assert df.loc[(2, 1), 'PERMX'] != df.loc[(1, 1), 'PERMX']                 # the compacted band differs in A


@pytest.mark.parametrize('case', CASES)
def test_scenario_nominals_equal_the_soil_table_written_to_cathy(MC, MS, mt_fake, case):
    """Independent of prior_bias_check: the unbiased prior nominal of every group = what the TRUTH project holds."""
    rec, _, _ = truth_build(MC, case)
    df = rec['soil_df']
    col = {name: MS.PRIOR[name]['col'] for name in ('SATCELL_VG', 'n_VG', 'thetar_VG', 'Ks') if name in MS.PRIOR}
    assert col, 'no soil parameter found in make_scenarios.PRIOR'
    for name, c in col.items():
        for g, nominal in zip(MS.groups(case)[2], MS.group_values(case, name, False)):
            zone, layers = (1, g) if case == 'NC' else (g[0], g[1])
            assert np.allclose(df.loc[[(zone, L + 1) for L in layers], c], nominal), (name, g)


@pytest.mark.parametrize('case', CASES)
def test_biased_da_model_differs_from_truth_only_by_biased_overrides(MC, mt_fake, case):
    t, _, ic_t = truth_build(MC, case)
    d, _, ic_d = da_build(MC, case, biased=True)
    assert_same(t, d, skip=('soil_df',))                      # mesh, ic, atmbc, parm, bc: untouched by the soil bias
    assert np.array_equal(ic_t, ic_d)
    diff = (t['soil_df'] != d['soil_df']).any(axis=0)
    assert set(diff[diff].index) == {'PERMX', 'PERMY', 'PERMZ', 'VGPSATCELL', 'VGNCELL', 'VGRMCCELL'}   # no POROS / ELSTOR
    f = MC.biased_overrides(case)
    for (zone, mat), o in f.items():
        for L in MC.MATERIAL_LAYERS_0[mat]:
            row = d['soil_df'].loc[(zone, L + 1)]
            for c, v in o.items():
                assert row[c] == v, (zone, mat, L, c)


def test_biased_hc_leaves_zone2_below_a_at_truth_values(MC, mt_fake):
    """CHARACTERISATION (not a verdict): biased_overrides only has the (zone, material) entries of CASES['HC']['soil'],
    i.e. zone 2 'A'. Zone 2 B and C (copied from zone 1 by zone2_fill_below) are NOT biased, while zone 1 B and C are.
    They are not DA parameters either (make_scenarios.groups('HC') has [2, A] only), so they stay at the truth.
    Decide whether the biased HC prior should move them too (see todo.md 5.3)."""
    t, d = truth_build(MC, 'HC')[0]['soil_df'], da_build(MC, 'HC', biased=True)[0]['soil_df']
    for L in MC.MATERIAL_LAYERS_0['B'] + MC.MATERIAL_LAYERS_0['C']:
        assert t.loc[(2, L + 1)].tolist() == d.loc[(2, L + 1)].tolist()              # zone 2: truth
        assert t.loc[(1, L + 1)].tolist() != d.loc[(1, L + 1)].tolist()              # zone 1: biased


# ------------------------------------------------------------------ initial condition
@pytest.mark.parametrize('case', CASES)
def test_ic_per_layer_prior_mean_equals_the_truth_ic_nodes(MC, mt_fake, case):
    """todo.md 5.1 (a): DA input/ic at cycle 0 = truth ic_nodes (per-layer ic_true)."""
    cfg = MC.get_case(case)
    ic_true = np.array([v[0] for v in cfg['ic_true']])
    t = truth_build(MC, case)[2]
    d = da_build(MC, case)[2]
    assert np.array_equal(t, d)
    per_level = t.reshape(MC.NSTR + 1, -1)[:, 0]              # one value per node level
    zr = np.asarray(MC.ZR_THICK)
    mid = np.cumsum(zr) - zr / 2
    expect = [ic_true[next((g for g, (a, b) in enumerate(cfg['ic_ranges']) if a <= m < b), 2)] for m in mid]
    assert np.allclose(per_level[:MC.NSTR], expect)


def test_old_scalar_ic_prior_is_the_only_difference_and_is_flagged(MC, mt_fake):
    """nb7 --ic-prior-mean -1.5: ic differs from the truth (known, +1.8 sigma in layer 0), nothing else does."""
    t, _, ic_t = truth_build(MC, 'NC')
    d, _, ic_d = da_build(MC, 'NC', ic_values=[[-1.5]] * 3)
    assert_same(t, d, skip=('ic',))
    assert not np.array_equal(ic_t, ic_d)
    top_t, top_d = ic_t.reshape(MC.NSTR + 1, -1)[0, 0], ic_d.reshape(MC.NSTR + 1, -1)[0, 0]
    assert (top_d - top_t) / 0.5 == pytest.approx(0.9 / 0.5)   # layer 0: -1.5 vs -2.4 = +0.9 m = +1.8 sigma


# ------------------------------------------------------------------ truth_parameters (what nb8 compares with)
@pytest.mark.parametrize('case', CASES)
def test_truth_parameters_describe_the_model_that_was_built(MC, case):
    tp, cfg = MC.truth_parameters(case), MC.get_case(case)
    assert tp['ic_true'] == cfg['ic_true'] and tp['ic_ranges'] == cfg['ic_ranges'] and tp['soil'] == cfg['soil']
    assert tp['material_layers'] == MC.MATERIAL_LAYERS_0 and tp['zr_thick'] == list(MC.ZR_THICK)
    assert tp['delta'] == cfg['delta'] and tp['band'] == cfg['band']


# ------------------------------------------------------------------ open items: need the real files / a run
NEEDS = 'needs the truth_<CASE>.pkl / a CATHY project / a run: not testable on config alone (todo.md 5.3)'


@pytest.mark.skip(reason=NEEDS)
def test_mesh_nodes_equal_truth_nodes():
    """truth['nodes'] vs the DA CATHY mesh points; same DEM / surface and bottom elevations."""


@pytest.mark.skip(reason=NEEDS)
def test_fwd_mesh_da_mesh_and_ert_frame_are_the_truth_ones():
    """truth['fwd_mesh'], truth['da_mesh'], ert_frame (signs, shift, origin) = the files and frame the DA reads."""


@pytest.mark.skip(reason=NEEDS)
def test_tdr_nodes_are_the_same_nodes_in_the_observation_operator():
    """truth['tdr_nodes'] / tdr_depths vs the DA observation operator."""


@pytest.mark.skip(reason=NEEDS)
def test_porosity_nodes_and_archie_nominals_equal_truth():
    """porosity_nodes and porosity_mode used by Archie; archie_truth vs archie_prior and rFluid when unbiased."""


@pytest.mark.skip(reason=NEEDS)
def test_survey_times_and_n_times_cut_match_the_truth():
    """build_times(ert_times) == truth atm_times (values, not only the assert on sizes); surveys[*]['t_h_irr'];
    --n-times N: first N surveys and the truth cut at the same last survey."""


@pytest.mark.skip(reason=NEEDS)
def test_observations_noise_and_error_model_match():
    """ERT data of the twin = truth data (sequence, geometry, noise type / % / seed); filter error model = noise added;
    tdr_sigma in the filter = the one used to perturb the virtual TDR."""


@pytest.mark.skip(reason=NEEDS)
def test_open_loop_noda_unbiased_unperturbed_reproduces_the_truth():
    """noDA, unbiased prior, N=1 at the nominal: states equal the truth sw / psi to solver precision."""
