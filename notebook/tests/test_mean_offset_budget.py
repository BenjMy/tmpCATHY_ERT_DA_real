"""mean_offset_budget: finite-N offset, van Genuchten sw, ic nonlinearity (no CATHY needed for the functions)."""
import numpy as np
import pytest

import mean_offset_budget as MB

NC_ZR = [0.03] * 4 + [0.04] * 3 + [0.06] * 2 + [0.08] * 2 + [0.10, 0.25, 0.40, 0.73]
NC_IC_RANGES = [[0.0, 0.20], [0.2, 0.44], [0.44, 2.0]]


def test_finite_n_matches_monte_carlo():
    rng = np.random.default_rng(1)
    n = 16
    means = rng.normal(size=(40_000, n)).mean(axis=1)
    exp_abs, b95 = MB.finite_n_offset(n)
    assert abs(np.abs(means).mean() - exp_abs) < 0.005
    assert abs(np.mean(np.abs(means) < b95) - 0.95) < 0.005


def test_finite_n_table_has_collapsed_rows():
    t = MB.finite_n_table([32], [2])
    assert list(t.n) == [32, 2] and t.expected_abs_sigma[1] > 0.5 > t.expected_abs_sigma[0]


def test_vg_limits_and_monotone():
    assert MB.vg_sw(0.0, 0.4, 1.4, 0.05, 0.4) == pytest.approx(0.4)          # saturated at psi = 0
    assert MB.vg_sw(1.0, 0.4, 1.4, 0.05, 0.4) == pytest.approx(0.4)          # positive head: still saturated
    assert MB.vg_sw(-1e6, 0.4, 1.4, 0.05, 0.4) == pytest.approx(0.05, abs=1e-3)
    psi = np.linspace(-5, -0.01, 50)
    assert np.all(np.diff(MB.vg_sw(psi, 0.4, 1.4, 0.05, 0.4)) > 0)


def test_layer_groups_follow_the_mid_depth_rule():
    """todo.md 5.1 (a): ic groups [0-5] [6-9] [10-14] for the NC layering."""
    g = MB.group_of_layers(NC_ZR, NC_IC_RANGES)
    assert [g.index(k) for k in (0, 1, 2)] == [0, 6, 10]
    assert g == [0] * 6 + [1] * 4 + [2] * 5


def test_nonlinearity_is_zero_for_zero_sigma_and_follows_curvature():
    soil = [dict(VGPSATCELL=0.44, VGNCELL=1.4, VGRMCCELL=0.04, POROS=0.415)]
    z = MB.ic_nonlinearity([-2.4], [1e-9], (-5, -0.05), [0], soil, n=2000)
    assert abs(z.offset_sw[0]) < 1e-6
    nl = MB.ic_nonlinearity([-2.4], [0.5], (-5, -0.05), [0], soil, n=100_000)
    assert nl.offset_sw[0] > 0          # sw(psi) is convex there: Jensen -> mean over the prior above sw(mean)
    assert abs(nl.offset_sigma[0]) < 0.3


def test_cli_on_nc():
    pytest.importorskip('pyCATHY')      # model_config / make_scenarios import pyCATHY at module level
    inp = MB.case_inputs('NC')
    assert inp['layer_ic_group'] == [0] * 6 + [1] * 4 + [2] * 5
    assert MB.main(['--case', 'NC', '--n', '32']) == 0
