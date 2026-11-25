"""CSV run registry: upsert, tracker (finished / crashed), backfill from disk, parallel-safe rewrite."""
import json
import os
import subprocess
import sys
import textwrap

import run_registry as RR

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_upsert_updates_in_place(tmp_path):
    p = str(tmp_path / 'reg' / 'r.csv')
    assert RR.upsert('A', p, kind='synthetic', case='NC', status='started', bogus='x')
    assert RR.upsert('B', p, status='started')
    assert RR.upsert('A', p, status='finished', duration_min=1.5)
    rows = RR.read(p)
    assert [r['run'] for r in rows] == ['A', 'B']
    a = rows[0]
    assert a['status'] == 'finished' and a['case'] == 'NC' and a['duration_min'] == '1.5'
    assert 'bogus' not in a and list(a) == RR.COLUMNS


def test_tracker_finish_and_failure(tmp_path):
    p = str(tmp_path / 'r.csv')
    t = RR.track('ok', p, kind='real')
    assert RR.read(p)[0]['status'] == 'started'
    t.finish('finished')
    r = RR.read(p)[0]
    assert r['status'] == 'finished' and r['finished'] and r['started'] and r['git_commit'] is not None
    t2 = RR.track('boom', p)
    t2._at_exit()                                  # what atexit does when the script died before finish()
    assert [x['status'] for x in RR.read(p)] == ['finished', 'failed']


def test_crashed_process_is_marked_failed(tmp_path):
    p = str(tmp_path / 'r.csv')
    code = textwrap.dedent(f"""
        import sys; sys.path.insert(0, {HERE!r})
        import run_registry as RR
        RR.track('crash', {p!r}, case='NC')
        raise SystemExit(0)
    """)
    subprocess.run([sys.executable, '-c', code], check=True)
    assert RR.read(p)[0]['status'] == 'failed'


def test_backfill_from_disk(tmp_path):
    sim = tmp_path / 'sims'
    (sim / 'SYN_NC_ic_enkf_N8_s1_unbiased_meta').mkdir(parents=True)
    (sim / 'SYN_NC_ic_enkf_N8_s1_unbiased_meta' / 'run_meta.json').write_text(json.dumps(dict(
        run='SYN_NC_ic_enkf_N8_s1_unbiased', case='NC', scenario='NC_ic', da_type='enkf_Evensen2009', nens=8, seed=1,
        prior_bias='none', with_da=True, status='finished', started='2026-10-24 10:00:00', finished='2026-10-24 10:25:00')))
    done = sim / 'NC_ERT_PF_withDA_nc_ERT_ic_mat'
    done.mkdir(parents=True)
    (done / 'NC_ERT_PF_withDA_nc_ERT_ic_mat_finished.txt').write_text('ok\n')
    (sim / 'NC_ERT_withoutDA_nc_ERT_ic_mat').mkdir()                  # no marker: incomplete
    p = str(tmp_path / 'r.csv')
    assert RR.backfill(str(sim), p, logs_dir=str(tmp_path)) == 3
    by = {r['run']: r for r in RR.read(p)}
    assert by['SYN_NC_ic_enkf_N8_s1_unbiased']['status'] == 'finished' and by['SYN_NC_ic_enkf_N8_s1_unbiased']['nens'] == '8'
    pf = by['NC_ERT_PF_withDA_nc_ERT_ic_mat']
    assert (pf['kind'], pf['da_type'], pf['status'], pf['scenario'], pf['with_da']) == ('real', 'pf', 'finished', 'nc_ERT_ic', 'True')
    assert by['NC_ERT_withoutDA_nc_ERT_ic_mat']['status'] == 'incomplete'
    assert RR.backfill(str(sim), p, logs_dir=str(tmp_path)) == 1      # finished rows are not touched again (only the incomplete one)


def test_cli_list_and_summary(tmp_path, capsys):
    p = str(tmp_path / 'r.csv')
    RR.upsert('A', p, status='finished', case='NC', da_type='enkf')
    RR.upsert('B', p, status='failed', case='HC', da_type='pf')
    assert RR.main(['list', '--registry', p, '--status', 'failed']) == 0
    out = capsys.readouterr().out
    assert 'B' in out and '1 run(s)' in out
    assert RR.main(['summary', '--registry', p]) == 0
    assert "'failed': 1" in capsys.readouterr().out


def test_same_name_in_two_folders_gives_two_rows(tmp_path):
    p = str(tmp_path / 'r.csv')
    base, logn = str(tmp_path / 'simulations_n10'), str(tmp_path / 'simulations_lognormal')
    RR.upsert('R', p, status='finished', git_commit='18430b5', folder=os.path.join(base, 'R'))
    RR.upsert('R', p, status='started', git_commit='ebe5786', folder=os.path.join(logn, 'R'))
    RR.upsert('R', p, status='finished', folder=os.path.join(logn, 'R'))
    rows = RR.read(p)
    assert len(rows) == 2
    assert RR.find('R', os.path.join(base, 'R'), p)['git_commit'] == '18430b5'      # not overwritten
    assert RR.find('R', os.path.join(logn, 'R'), p)['git_commit'] == 'ebe5786'
    assert RR.find('R', os.path.join(tmp_path, 'other', 'R'), p) is None


def test_legacy_row_without_folder_is_adopted(tmp_path):
    p = str(tmp_path / 'r.csv')
    RR.upsert('R', p, status='finished')
    RR.upsert('R', p, status='started', folder=str(tmp_path / 'sims' / 'R'))
    rows = RR.read(p)
    assert len(rows) == 1 and rows[0]['folder'] == str(tmp_path / 'sims' / 'R')


def test_tracker_finish_updates_its_own_folder_row(tmp_path):
    p = str(tmp_path / 'r.csv')
    RR.upsert('R', p, status='finished', folder=str(tmp_path / 'a' / 'R'))
    t = RR.track('R', p, manifest=False, folder=str(tmp_path / 'b' / 'R'))
    t.finish('failed')
    by = {os.path.basename(os.path.dirname(r['folder'])): r['status'] for r in RR.read(p)}
    assert by == {'a': 'finished', 'b': 'failed'}


def test_cli_list_sim_dir_filter(tmp_path, capsys):
    p = str(tmp_path / 'r.csv')
    RR.upsert('A', p, status='finished', folder=str(tmp_path / 's1' / 'A'))
    RR.upsert('B', p, status='finished', folder=str(tmp_path / 's2' / 'B'))
    assert RR.main(['list', '--registry', p, '--sim-dir', str(tmp_path / 's1')]) == 0
    out = capsys.readouterr().out
    assert 'A ' in out and 'B ' not in out and '1 run(s)' in out
