"""export_results.py: copies the small files, one overview row per run, never the big ones."""
import csv
import json
import os

import export_results as E
import project_paths as PP
import run_registry as RR

OK, BAD, NOREG = 'SYN_NC_full_upd_pf_N32_s1_unbiased', 'SYN_NC_full_upd_enkf_N32_s1_noDA_unbiased', 'SYN_NC_ic_enkf_N8_s1'


def _project(tmp_path):
    (tmp_path / 'data_path.toml').write_text('[project]\nroot = "."\n[code]\nnotebook = "nb"\n')
    P = PP.ProjectPaths(str(tmp_path / 'data_path.toml'))
    sim = tmp_path / 'simulations'
    for run in (OK, BAD, NOREG):
        (sim / run / 'DA_Ensemble').mkdir(parents=True)
        (sim / run / 'DA_Ensemble' / 'm.vtk').write_bytes(b'v' * 2_000_000)
        (sim / run / 'x.txt').write_bytes(b'x' * 100)
    for run in (OK, BAD):
        meta = sim / (run + '_meta')
        meta.mkdir()
        (meta / 'manifest.json').write_text(json.dumps({'command': 'nb7 ...'}))
        (meta / 'run_meta.json').write_text('{}')
        (sim / run / (run + '_states')).mkdir()
        (sim / run / (run + '_states') / 'cycle_0001.npz').write_bytes(b'n' * 50)
    (sim / (BAD + '_issues.log')).write_text('E: something broke\n')
    (sim / OK / (OK + '_states.big')).write_bytes(b'b' * 10)
    logs = tmp_path / 'notebook' / 'logs'
    logs.mkdir(parents=True)
    (logs / (BAD + '.log')).write_text('\n'.join(f'line {i}' for i in range(500)))
    RR.upsert(OK, P.run_registry, kind='synthetic', case='NC', scenario='NC_full_upd', da_type='pf', nens=32, seed=1,
              status='finished', duration_min=12.5)
    RR.upsert(BAD, P.run_registry, kind='synthetic', case='NC', da_type='enkf', nens=32, seed=1, status='failed',
              log=str(logs / (BAD + '.log')))
    fig = sim / 'figures_nb8' / 'NC'
    fig.mkdir(parents=True)
    with open(fig / 'summary.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['run', 'scenario', 'sw_all_rmse', 'sw_all_coverage', 'states_kind', 'unrelated'])
        w.writerow([OK, 'NC_full_upd', 0.031, 0.88, 'forecast_before_update', 'zzz'])
    (fig / 'parameters_summary.csv').write_text('run,parameter\n')
    (fig / 'compare_heatmap_NC.png').write_bytes(b'png')
    (fig / (OK + '_metrics.png')).write_bytes(b'png')
    (fig / ('SYN_OTHER_metrics.png')).write_bytes(b'png')
    return P


def _overview(out):
    with open(os.path.join(out, 'overview.csv'), newline='') as fh:
        return {r['run']: r for r in csv.DictReader(fh)}


def test_overview_joins_registry_disk_and_scores(tmp_path):
    P = _project(tmp_path)
    out = str(tmp_path / 'out')
    runs, _ = E.export(P, out)
    ov = _overview(out)
    assert set(ov) == {OK, BAD, NOREG}                              # a run on disk without registry row is kept
    assert ov[OK]['status'] == 'finished' and ov[OK]['duration_min'] == '12.5' and ov[OK]['nens'] == '32'
    assert ov[OK]['sw_all_rmse'] == '0.031' and ov[OK]['states_kind'] == 'forecast_before_update'
    assert 'unrelated' not in ov[OK] and ov[BAD]['sw_all_rmse'] == ''
    assert float(ov[OK]['disk_mb']) > 0 and float(ov[OK]['vtk_mb']) > 0 and ov[OK]['compact_store'] == 'True'
    assert ov[NOREG]['status'] == 'no registry row' and ov[NOREG]['compact_store'] == 'False'
    assert 'manifest.json' in ov[OK]['files']


def test_copies_small_files_only(tmp_path):
    P = _project(tmp_path)
    out = str(tmp_path / 'out')
    E.export(P, out)
    listing = [os.path.relpath(os.path.join(r, f), out) for r, _, fs in os.walk(out) for f in fs]
    assert not any(f.endswith(('.vtk', '.npz', '.nc', '.pkl', '.big')) for f in listing)
    assert os.path.isfile(os.path.join(out, 'runs', OK, 'manifest.json'))
    assert os.path.isfile(os.path.join(out, 'runs', BAD, BAD + '_issues.log'))
    tail = open(os.path.join(out, 'runs', BAD, BAD + '_log_tail.txt')).read().splitlines()
    assert len(tail) == E.TAIL_LINES and tail[-1] == 'line 499'
    assert not os.path.exists(os.path.join(out, 'runs', OK, OK + '_log_tail.txt'))      # finished run: no log tail
    assert os.path.isfile(os.path.join(out, 'tables', 'NC', 'summary.csv'))
    assert os.path.isfile(os.path.join(out, 'tables', 'runs_registry.csv')) and os.path.isfile(os.path.join(out, 'README.txt'))


def test_figures_modes_and_filters(tmp_path):
    P = _project(tmp_path)
    figs = lambda o: sorted(os.listdir(os.path.join(o, 'figures', 'NC'))) if os.path.isdir(os.path.join(o, 'figures')) else []
    E.export(P, str(tmp_path / 'a'))
    assert figs(str(tmp_path / 'a')) == ['compare_heatmap_NC.png']
    E.export(P, str(tmp_path / 'b'), figures='all')
    assert figs(str(tmp_path / 'b')) == sorted(['compare_heatmap_NC.png', OK + '_metrics.png'])    # other runs' figures left out
    E.export(P, str(tmp_path / 'c'), figures='none')
    assert figs(str(tmp_path / 'c')) == []
    runs, _ = E.export(P, str(tmp_path / 'd'), globs=['*_pf_*'], only_ok=True)
    assert [r['run'] for r in runs] == [OK]


def test_zip_and_cli(tmp_path):
    P = _project(tmp_path)
    out = str(tmp_path / 'share')
    _, res = E.export(P, out, make_zip=True)
    assert res.endswith('.zip') and os.path.getsize(res) > 0
    assert E.main([str(tmp_path / 'cli'), '--toml', str(tmp_path / 'data_path.toml'), '--figures', 'none']) == 0
    assert os.path.isfile(str(tmp_path / 'cli' / 'overview.csv'))


def test_empty_project(tmp_path):
    (tmp_path / 'data_path.toml').write_text('[project]\nroot = "."\n')
    runs, _ = E.export(PP.ProjectPaths(str(tmp_path / 'data_path.toml')), str(tmp_path / 'o'))
    assert runs == [] and os.path.isfile(str(tmp_path / 'o' / 'overview.csv'))
