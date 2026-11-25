"""Helpers called by plots_synthetic.py: nearest_index, parameter_evolution, parameter_table(params_ds=...)."""
import numpy as np
import pytest

import plots_synthetic_helpers as H


def test_nearest_index():
    t = np.array([0.0, 10.0, 20.0, 30.0])
    assert H.nearest_index(t, 11.0) == 1
    assert H.nearest_index(t, 100.0) == 3
    assert H.nearest_index([], 1.0) is None


def _entry(n_cycle=3, n_member=8, seed=0):
    rng = np.random.default_rng(seed)
    e = {'ini_perturbation': rng.normal(1.0, 0.3, n_member)}
    for k in range(1, n_cycle + 1):
        e[f'update_nb{k}'] = rng.normal(1.0, 0.3 / (1 + k), n_member)
    return e


def test_parameter_evolution_from_dict():
    dpp = {'Ks1': _entry()}
    ev = H.parameter_evolution('Ks1', dpp)
    assert ev.shape == (4, 8)
    assert np.allclose(ev[0], dpp['Ks1']['ini_perturbation'])


def test_parameter_evolution_missing_raises():
    with pytest.raises(KeyError):
        H.parameter_evolution('nope', {'Ks1': _entry()})


def test_parameter_evolution_from_dataset_member_first():
    xr = pytest.importorskip('xarray')
    arr = np.arange(8 * 3, dtype=float).reshape(8, 3)       # (member, iteration)
    ds = xr.Dataset({'Ks1': (('member', 'iteration'), arr)})
    ev = H.parameter_evolution('Ks1', {}, ds)
    assert ev.shape == (3, 8) and np.allclose(ev, arr.T)


def test_parameter_table_accepts_params_ds_and_keeps_old_call(monkeypatch):
    monkeypatch.setattr(H, 'truth_of_parameter', lambda *a, **k: 1.0)
    dpp = {'Ks1': _entry()}
    old = H.parameter_table(dpp, ['Ks1', 'absent'], {}, None, None)            # old signature still works
    new = H.parameter_table(dpp, ['Ks1', 'absent'], {}, None, None, params_ds=None)
    assert list(old.index) == ['Ks1'] and old.equals(new)
    assert {'z_prior', 'z_final', 'contraction', 'truth_in_90'} <= set(old.columns)
    assert H.parameter_table({}, ['Ks1'], {}, None, None).empty
