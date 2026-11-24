"""prior_bias_check: sampler vs analytic moments, and what it finds in the generated (unbiased) NC scenarios."""
import numpy as np
import pytest
from scipy.stats import truncnorm

import prior_bias_check as PB


def test_truncated_normal_matches_analytic_mean():
    # Ks prior: sigma = nominal, lower bound = nominal / 100 -> truncation at -0.99 sigma
    x = PB.sample_prior('normal', 1.0, 1.0, 0.01, 100.0, n=400_000)
    exact = truncnorm((0.01 - 1.0) / 1.0, (100.0 - 1.0) / 1.0, loc=1.0, scale=1.0).mean()
    assert abs(x.mean() - exact) < 0.01
    assert 1.25 < exact < 1.33            # the mean is ~29 % above the nominal


def test_lognormal_centre_median_vs_mean():
    s = 0.7
    med = PB.sample_prior('lognormal', 2.0, s, n=400_000, lognormal_center='median')
    mean = PB.sample_prior('lognormal', 2.0, s, n=400_000, lognormal_center='mean')
    assert abs(np.median(med) / 2.0 - 1) < 0.01 and abs(med.mean() / (2.0 * np.exp(s ** 2 / 2)) - 1) < 0.02
    assert abs(mean.mean() / 2.0 - 1) < 0.02


def test_unbounded_normal_is_centred():
    x = PB.sample_prior('normal', 3.0, 0.5, n=200_000)
    assert abs(x.mean() - 3.0) < 0.01


@pytest.fixture(scope='module')
def nc_unbiased():
    pytest.importorskip('pyCATHY')        # model_config / make_scenarios import pyCATHY at module level
    return PB.prior_table('NC', 'full_upd', biased=False, n=100_000)


def test_unbiased_nominals_equal_truth(nc_unbiased):
    d = nc_unbiased[~nc_unbiased.param.str.startswith('ic')]
    assert np.allclose(d.nominal_over_truth, 1.0)


def test_unbiased_soil_and_archie_priors_are_centred_on_the_working_scale(nc_unbiased):
    d = nc_unbiased[~nc_unbiased.param.str.startswith('ic')]
    assert d.bias_sigma.abs().max() < 0.1
    ks = d[d.param == 'Ks']
    assert (ks.scale == 'log').all() and ks.mean_over_truth.between(1.2, 1.4).all()   # physical mean +29 %, log-mean ~0


def test_unbiased_ic_prior_is_off_centre_known_issue(nc_unbiased):
    """KNOWN ISSUE (todo.md s.5): ic is ONE scalar (-1.5) for 3 layers, the truth is [-2.4, -1.3, -1.15].
    When the IC prior is fixed (per-layer ic, or a centred scalar) flip this test."""
    ic = nc_unbiased[nc_unbiased.param.str.startswith('ic')].set_index('param')
    assert ic.flag.all()
    assert abs(ic.loc['ic (layer 0)', 'bias_sigma'] - 1.8) < 0.05
