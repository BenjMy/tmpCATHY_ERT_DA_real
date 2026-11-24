"""summarize_run: the open loop is a separate run with its own cycles -> matched by time, not by position."""
import numpy as np
import pandas as pd

import plots_synthetic_helpers as H


def _sc(t, rmse):
    n = len(t)
    return pd.DataFrame(dict(rmse=rmse, bias=0.0, spread_skill=1.0, coverage=0.9, crps=0.1, ), index=pd.Index(t, dtype=float))


def test_open_loop_longer_than_da_run():                  # the crash: 9 DA cycles vs 28 open-loop cycles
    t_da = np.arange(9) * 24.0
    t_ol = np.arange(28) * 24.0
    da, ol = _sc(t_da, 0.02), _sc(t_ol, np.linspace(0.05, 0.08, 28))
    r = H.summarize_run(da, ol, 0.0, t_da)
    assert np.isclose(r['rmse_ol'], np.linspace(0.05, 0.08, 28)[:9].mean())
    assert np.isclose(r['skill_vs_ol'], 1 - 0.02 / r['rmse_ol'])


def test_open_loop_shorter_and_offset_times():
    t_da = np.arange(6) * 24.0
    ol = _sc([0.0, 24.001, 48.0], [0.1, 0.2, 0.3])        # 3 cycles, one 3.6 s off, rest missing
    r = H.summarize_run(_sc(t_da, 0.05), ol, 0.0, t_da)
    assert np.isclose(r['rmse_ol'], 0.2)


def test_t_from_h_and_no_common_time():
    t_da = np.arange(4) * 24.0
    ol = _sc(t_da, [0.1, 0.1, 0.3, 0.3])
    assert np.isclose(H.summarize_run(_sc(t_da, 0.05), ol, 48.0, t_da)['rmse_ol'], 0.3)
    far = _sc([1000.0, 1024.0], [0.1, 0.1])
    r = H.summarize_run(_sc(t_da, 0.05), far, 0.0, t_da)
    assert np.isnan(r['rmse_ol']) and np.isnan(r['skill_vs_ol'])


def test_without_open_loop_unchanged():
    t = np.arange(3) * 24.0
    r = H.summarize_run(_sc(t, [0.1, 0.2, 0.3]))
    assert 'rmse_ol' not in r and np.isclose(r['rmse'], 0.2) and np.isclose(r['rmse_final'], 0.3)
