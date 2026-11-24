"""Per-layer initial condition (todo.md 5.1 (a)): ic = one value per depth range, prior mean = ic_true.
Needs pyCATHY with the 'material' control for ic (pycathy_wrapper feat/ic-per-layer-material)."""
import os
import sys
import tomllib

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
pytest.importorskip('pyCATHY')
import make_scenarios as MS  # noqa: E402
import model_config as mc  # noqa: E402
from pyCATHY import meshtools as mt  # noqa: E402
from pyCATHY.DA import perturbate  # noqa: E402
from pyCATHY.DA.scenario import validate_scenario  # noqa: E402

pytestmark = pytest.mark.skipif('material' not in perturbate.LAYER_LIKE or
                                'material' not in __import__('pyCATHY.DA.scenario', fromlist=['x']).SUPPORTED_CONTROLS['ic'],
                                reason='pyCATHY without the per-group ic (feat/ic-per-layer-material)')
NL = mc.NSTR


def test_ic_groups_are_the_depth_ranges_of_the_truth():
    g = MS.ic_groups()
    assert sorted(l for ls in g for l in ls) == list(range(NL))          # every layer once
    # same expansion as setup_model (map_param_to_layers on the layer depths)
    bot = np.cumsum(mc.ZR_THICK)
    top = bot - np.asarray(mc.ZR_THICK)
    for case in ('NC', 'HC'):
        truth = mc.get_case(case)['ic_true']
        by_depth = mt.map_param_to_layers(top, bot, mc.IC_RANGES, truth)
        by_group = np.empty(NL)
        for gi, ls in enumerate(g):
            by_group[ls] = truth[gi][0]
        assert np.allclose(by_depth, by_group), case
    assert 9 in g[1] and 9 in mc.MATERIAL_LAYERS_0['C']                  # why ic needs its own groups


@pytest.mark.parametrize('case', ['NC', 'HC'])
@pytest.mark.parametrize('biased', [True, False])
def test_every_scenario_has_ic_centred_on_truth_and_is_valid(case, biased):
    allsc = tomllib.loads(MS.build(biased=biased))['scenario']
    truth = [v[0] for v in mc.get_case(case)['ic_true']]
    tags = [k for k in allsc if k.startswith(case + '_')]
    assert tags
    for k in tags:
        sc = allsc[k]
        i = sc['per_name'].index('ic')
        assert sc['pert_control_name'][i] == 'material'
        assert sc['per_mean'][i] == truth and sc['per_nom'][i] == truth       # prior mean = ic_true
        assert sc['ic_layers'] == MS.ic_groups()
        zone = np.ones((2, 2)) if case == 'NC' else np.array([[1., 2.], [1., 2.]])
        errs, _ = validate_scenario(sc, nlayers=NL, zone_raster=zone, raise_on_error=False)
        assert not errs, (k, errs)


def test_perturbate_builds_ic0_ic1_ic2_with_the_truth_as_mean():
    idx = pd.MultiIndex.from_product([[1], range(NL)], names=['zone', 'layer'])

    class Sim:
        dem_parameters = {'base': 2.0}
        soil = dict(PCWLT=[-150.0], OMGC=[1.0], PZ=[1.0])

        def read_inputs(self, key, **kw):
            return {'soil': (pd.DataFrame({'PERMX': 1e-6}, index=idx), pd.DataFrame({'ZROOT': [0.4]}, index=[1])),
                    'zone': (np.ones((2, 2)), {}), 'root_map': (np.ones((2, 2)), {})}[key]

    sc = tomllib.loads(MS.build(biased=False))['scenario']['NC_full_upd']
    ent = [e for e in perturbate.perturbate(Sim(), sc, 8, check_scenario=False) if e['type_parm'].startswith('ic')]
    assert [e['type_parm'] for e in ent] == ['ic0', 'ic1', 'ic2']
    assert [e['mean'] for e in ent] == [v[0] for v in mc.get_case('NC')['ic_true']]
    assert all(e['layer_groups'] == MS.ic_groups() for e in ent)


def test_esmda_still_finds_the_ic_parameters():
    from pyCATHY.DA.esmda_cathy import resolve_parameter_names
    d = {'ic0': {}, 'ic1': {}, 'ic2': {}, 'a_Archie': {}, 'Ks0': {}}
    assert resolve_parameter_names(d, ['St. var.', 'ic']) == ['ic0', 'ic1', 'ic2']
