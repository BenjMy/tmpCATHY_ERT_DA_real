"""lighten_run.py (report / convert / prune) and lighter_protocol.py (measure / accept)."""
import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

import lighten_run as L
import lighter_protocol as LP


def _run(tmp_path, compact=True, legacy=True):
    d = tmp_path / 'SYN_NC_x'
    (d / 'DA_Ensemble' / 'm1').mkdir(parents=True)
    (d / 'DA_Ensemble' / 'm1' / 'a.vtk').write_bytes(b'x' * 5000)
    (d / 'b.txt').write_bytes(b'y' * 100)
    if compact:
        (d / 'SYN_NC_x_states').mkdir()
        (d / 'SYN_NC_x_states' / 'cycle_0001.npz').write_bytes(b'z' * 50)
    if legacy:
        pd.DataFrame({'time': [0, 0, 1, 1], 'Ensemble_nb': [1, 1, 1, 1], 'OL': [0, 0, 0, 0],
                      'sw_bef_update_': np.linspace(0, 1, 4)}).to_pickle(d / 'SYN_NC_x_df_DA.pkl')
    return str(d)


def test_report_categories_and_vtk_hint(tmp_path, capsys):
    d = _run(tmp_path)
    assert L.report([d]) > 5000
    out = capsys.readouterr().out
    assert 'vtk' in out and 'compact store' in out and 'mostly vtk' in out and 'total:' in out


def test_prune_needs_compact_store(tmp_path):
    d = _run(tmp_path, compact=False)
    assert L.prune(d, apply=True) == []
    assert os.path.exists(os.path.join(d, 'SYN_NC_x_df_DA.pkl'))


def test_prune_dry_run_then_apply(tmp_path):
    d = _run(tmp_path)
    f = os.path.join(d, 'SYN_NC_x_df_DA.pkl')
    assert L.prune(d) == [f] and os.path.exists(f)
    assert L.prune(d, apply=True) == [f] and not os.path.exists(f)
    assert os.path.exists(os.path.join(d, 'DA_Ensemble', 'm1', 'a.vtk'))       # prune never touches the vtk


def test_convert_float32_idnode_keep(tmp_path):
    pytest.importorskip('pyarrow')
    d = _run(tmp_path)
    pkl = os.path.join(d, 'SYN_NC_x_df_DA.pkl')
    out, _, _ = L.convert(pkl, keep=True)
    df = pd.read_parquet(out)
    assert df['sw_bef_update_'].dtype == np.float32 and list(df['idnode']) == [0, 1, 0, 1]
    assert os.path.exists(pkl)
    L.convert(pkl, keep=False)
    assert not os.path.exists(pkl)


def test_main_cli(tmp_path, capsys):
    d = _run(tmp_path)
    L.main(['prune', d])
    assert 'would delete' in capsys.readouterr().out


def _summary(path, rmse, cov=0.9, nens=8):
    pd.DataFrame([dict(scenario='ic_Ks_upd', da='enkf', nens=nens, sw_all_rmse=rmse, sw_all_coverage=cov)]).to_csv(path, index=False)
    return str(path)


def test_accept_within_and_outside_spread(tmp_path):
    base = [_summary(tmp_path / f'b{i}.csv', r) for i, r in enumerate([0.040, 0.044, 0.042])]   # std = 0.002
    ok, tab = LP.accept(base, [_summary(tmp_path / 'l1.csv', 0.0435)])
    assert ok and set(tab['column']) == {'sw_all_rmse', 'sw_all_coverage'}
    ok, tab = LP.accept(base, [_summary(tmp_path / 'l2.csv', 0.050)])
    assert not ok and not tab.set_index('column').loc['sw_all_rmse', 'ok']


def test_accept_needs_two_baselines_and_matching_rows(tmp_path):
    b = _summary(tmp_path / 'b.csv', 0.04)
    with pytest.raises(ValueError, match='>= 2'):
        LP.accept([b], [b])
    b2 = _summary(tmp_path / 'b2.csv', 0.05)
    with pytest.raises(ValueError, match='matches'):
        LP.accept([b, b2], [_summary(tmp_path / 'l.csv', 0.04, nens=16)])


def test_accept_tolerance_floor_and_exit_code(tmp_path, capsys):
    base = [_summary(tmp_path / f'b{i}.csv', 0.04) for i in range(2)]        # std = 0: only the floor allows a shift
    light = _summary(tmp_path / 'l.csv', 0.041)
    assert LP.main(['accept', *base, '--light', light]) == 1
    assert LP.main(['accept', *base, '--light', light, '--tol', '0.01']) == 0
    assert 'ACCEPTED' in capsys.readouterr().out


def test_measure_records_time_rss_and_disk(tmp_path):
    d = _run(tmp_path)
    out = str(tmp_path / 'm' / 'm.json')
    rc = LP.main(['measure', out, '--run-dir', d, '--', sys.executable, '-c', 'x = bytearray(30_000_000)'])
    r = json.load(open(out))
    assert rc == 0 and r['peak_rss_mb'] > 30 and r['wall_s'] >= 0 and r['disk_bytes']['vtk'] == 5000
    assert r['disk_total_bytes'] == sum(r['disk_bytes'].values())
    assert LP.main(['measure', out, '--', sys.executable, '-c', 'raise SystemExit(3)']) == 3
