"""diag_pf_run: checks of a DA run from its small files (fixtures written here, no pyCATHY)."""
import json
import math

import numpy as np

import diag_pf_run as D

NAME = 'SYN_NC_full_upd_pf_N8_s1_unbiased'
NAME_B = 'SYN_NC_full_upd_enkf_N8_s1_unbiased'


def _par(nominal, prior_sd, final_mean, final_sd, n_updates=5, clip_min=0.0):
    return dict(nominal=nominal, prior_sd=prior_sd, final_mean=final_mean, final_sd=final_sd, n_updates=n_updates,
                clip_min=clip_min, frac_at_bound_final=0.0, contraction=1 - final_sd / prior_sd)


def _write(d, name=NAME, collapse=True, status='failed', states=False):
    run = dict(library_status=status, n_cycles=3, n_stored_cycles=3,
               options=dict(da_type=D.RUN_RE.match(name).group('da'), n_members=8, seed=1),
               rejected_members=dict(after_run=[3, 4, 5]),
               parameters=dict(
                   Ks0=_par(2e-6, 1e-6, 1.0e-6 if collapse else 1.9e-6, 1e-7 if collapse else 8e-7, clip_min=1.5e-8),
                   n_VG0=_par(1.4, 0.07, 1.5 if collapse else 1.41, 0.01 if collapse else 0.06),
                   ic0=_par(-2.4, 0.5, -2.4, 0.5, n_updates=0)))
    (d / (name + '_run.json')).write_text(json.dumps(run))
    (d / (name + '_rejections.json')).write_text(json.dumps([
        dict(cycle=1, stage='after_run', member=m, reason=f"df_mbeconv['time']{3400 + m} - CATHY stopped at t={3400 + m} s "
             "before the end of the window (last printing time 3600 s)") for m in (3, 4, 5)]))
    lines = [
        "ERROR\trun\tx1\t3\trun aborted: KeyboardInterrupt:",
        "WARNING\tlog\tx2\t0-2\t⚠️ CRITICAL: Extreme weight concentration - almost all weight on one particle",
        "WARNING\tlog\tx3\t0-2\t⚠️ WARNING: Very low N_eff - particle degeneracy",
        "WARNING\tlog\tx1\t1\tKs0: 3 analysed value(s) reflected to [1.5e-08, 0.00023] (members [1, 2, 3]): "
        "[-4.3e-07  1.9e-07  1.7e-07] -> [4.6e-07 1.9e-07 1.7e-07]",
        "WARNING\tlog\tx2\t2-3\tn_VG0: 3 analysed value(s) reflected to [1.05, 2.5] (members [3, 4, 5]): [nan nan nan] -> [nan nan nan]",
        "WARNING\tlog\tx1\t1\tSATCELL_VG0: 2 analysed value(s) reflected to [0.1, 2.0] (members [1, 2]): "
        "[0.49 0.47] -> [0.49 0.47]",
    ]
    (d / (name + '_issues.log')).write_text('\n'.join(lines) + '\n')

    def ens(c, n_mem, misfit, spread, fm):
        return dict(name='e', ok=True, stats=dict(n_members=n_mem, n_obs=100, median_abs_log10_misfit=misfit,
                                                  median_rel_spread=spread, first_mode_var_frac=fm))
    chk = []
    for c, (nm, mis, spr, fm, ri) in enumerate([(8, 0.1, 0.4, 0.9, 0.2), (5, 0.03, 0.08, 0.9, 0.6), (5, 0.05, 0.14, 0.95, 0.9)]):
        chk.append(dict(stage='forward', cycle=c, members={'1': dict(ok=True), '2': dict(ok=(c != 1))},
                        ensemble=ens(c, nm, mis, spr, fm)))
        chk.append(dict(stage='analysis', cycle=c, report=dict(name='a', ok=True, stats=dict(
            n_members=nm, n_obs=100, median_abs_log10_misfit=mis, median_rel_spread=spr, first_mode_var_frac=fm,
            rms_normalised_innovation=ri))))
    (d / (name + '_ERT_checks.json')).write_text(json.dumps(chk))
    if states:
        sd = d / (name + '_states')
        sd.mkdir()
        rng = np.random.default_rng(0)
        psi = rng.normal(-1.5, 0.5, (30, 8))
        for c in range(3):
            par = np.zeros(8, int)                               # everybody descends from member 0
            np.savez(sd / f'cycle_{c:05d}.npz', cycle=c, is_analysis=True, rejected=np.zeros(8, np.int8),
                     psi=psi.astype(np.float32), analysis=psi[:, par].astype(np.float32))
            psi = psi[:, par] + rng.normal(0, 0.01, (30, 8))


def _levels(rep, section):
    s = next(s for s in rep['sections'] if s['title'].startswith(section))
    return [(i['level'], i['text']) for i in s['items']]


def test_parse_issues_and_cycles(tmp_path):
    _write(tmp_path)
    iss = D.parse_issues(str(tmp_path / (NAME + '_issues.log')))
    assert iss[0]['level'] == 'ERROR' and iss[0]['c0'] == 3
    assert iss[1]['count'] == 2 and (iss[1]['c0'], iss[1]['c1']) == (0, 2)
    refl = D.parse_reflections(iss)
    assert [r['par'] for r in refl] == ['Ks0', 'n_VG0', 'SATCELL_VG0']
    assert refl[0]['old'][0] < 0 and math.isnan(refl[1]['old'][0])


def test_locate_by_folder_and_missing_files(tmp_path):
    _write(tmp_path)
    f = D.locate(str(tmp_path))
    assert f['name'] == NAME and f['run'] and f['issues'] and f['rejections'] and f['ert'] and f['states'] is None
    empty = tmp_path / 'empty'
    empty.mkdir()
    rep = D.diagnose(str(empty))                          # nothing found: no crash, only info
    assert all(i['level'] in (D.INFO, D.OK, D.WARN) for s in rep['sections'] for i in s['items'])


def test_collapsed_run_is_flagged(tmp_path):
    _write(tmp_path, collapse=True, states=True)
    rep = D.diagnose(str(tmp_path))
    assert any(l == D.WARN and 'partial' in t for l, t in _levels(rep, '1'))             # interrupted
    g = _levels(rep, '2')
    assert sum(1 for l, _ in g if l == D.BAD) >= 3                                       # weights x2 + inbreeding
    assert any('1 distinct founders' in t for _, t in g)
    m = _levels(rep, '3')
    assert any(l == D.BAD and 'NaN' in t and 'not refilled' in t for l, t in m)
    assert any('consecutive' in t for _, t in m)
    p = _levels(rep, '4')
    assert any(l == D.BAD and '2 of 2 updated parameters' in t for l, t in p)           # ic0 (n_updates 0) not counted
    assert any(l == D.BAD and 'confident and wrong' in t for l, t in p)
    assert any('negative analysed values' in t and 'Ks0' in t for _, t in p)
    assert any('identical' in t and 'SATCELL_VG0' in t for _, t in p)
    e = _levels(rep, '5')
    assert any('first mode' in t and l == D.WARN for l, t in e)
    assert any('members used in the ERT analysis' in t for _, t in e)


def test_healthy_run_has_no_bad(tmp_path):
    _write(tmp_path, collapse=False, status='finished')
    (tmp_path / (NAME + '_issues.log')).write_text('')
    (tmp_path / (NAME + '_rejections.json')).write_text('[]')
    rep = D.diagnose(str(tmp_path))
    run = json.loads((tmp_path / (NAME + '_run.json')).read_text())
    run['rejected_members'] = {}
    (tmp_path / (NAME + '_run.json')).write_text(json.dumps(run))
    rep = D.diagnose(str(tmp_path))
    assert not [i for s in rep['sections'] for i in s['items'] if i['level'] == D.BAD]


def test_nominal_is_truth_auto_and_override(tmp_path):
    _write(tmp_path, collapse=True)
    biased = tmp_path / 'b'
    biased.mkdir()
    _write(biased, name='SYN_NC_full_upd_pf_N8_s1', collapse=True)
    assert any('confident and wrong' in t for _, t in _levels(D.diagnose(str(tmp_path)), '4'))
    rep = D.diagnose(str(biased))                                                       # no '_unbiased' in the name
    assert not any('confident and wrong' in t for _, t in _levels(rep, '4'))
    rep = D.diagnose(str(biased), nominal_is_truth='yes')
    assert any('confident and wrong' in t for _, t in _levels(rep, '4'))


def test_ert_table_ratio():
    chk = [dict(stage='analysis', cycle=0, report=dict(stats=dict(median_abs_log10_misfit=0.05, median_rel_spread=0.23,
                                                                  n_members=8, n_obs=10, first_mode_var_frac=0.5,
                                                                  rms_normalised_innovation=1.0)))]
    rows = D.ert_table(chk)
    assert abs(rows[0]['spread_over_misfit'] - 0.23 / (math.log(10) * 0.05)) < 1e-12
    low = [dict(stage='analysis', cycle=0, report=dict(stats=dict(median_abs_log10_misfit=0.1, median_rel_spread=0.05)))]
    assert D.check_ert(D.ert_table(low))[0][0] == D.WARN                                 # under-dispersed in data space


def test_compare_attributes_common_bias(tmp_path):
    a = tmp_path / 'a'; b = tmp_path / 'b'
    a.mkdir(); b.mkdir()
    _write(a, collapse=True)
    _write(b, name=NAME_B, collapse=True)
    ra, rb = D.diagnose(str(a)), D.diagnose(str(b))
    cmp = D.compare(ra, rb)
    assert cmp['message']['level'] == D.WARN and 'not specific to the filter' in cmp['message']['text']
    _write(b, name=NAME_B, collapse=False)                                               # EnKF stays near the nominal
    cmp = D.compare(ra, D.diagnose(str(b)))
    assert 'specific to the filter' in cmp['message']['text'] and 'not specific' not in cmp['message']['text']


def test_suggest_builds_commands_from_the_run_name():
    s = '\n'.join(t for _, t in D.suggest('SYN_NC_full_upd_pf_N32_s1_unbiased', 32, 'pf'))
    assert 'CASE=NC BIAS=none SEED=1 DA=enkf NENS=32 ./run_twin.sh pair full_upd' in s
    assert '--compare SYN_NC_full_upd_enkf_N32_s1_unbiased' in s and 'pf_jitter_std_state' in s
    assert 'pf_jitter' not in '\n'.join(t for _, t in D.suggest('SYN_NC_full_upd_enkf_N32_s1_unbiased', 32, 'enkf'))
    assert 'not recognised' in D.suggest('whatever', 8, 'pf')[0][1]


def test_main_json_and_strict_exit_code(tmp_path, capsys):
    _write(tmp_path, collapse=True)
    out = tmp_path / 'diag.json'
    assert D.main([str(tmp_path), '--json', str(out)]) == 0
    assert D.main([str(tmp_path), '--strict']) == 1
    data = json.loads(out.read_text())
    assert data['name'] == NAME and data['nens'] == 8 and data['sections']
    assert 'BAD' in capsys.readouterr().out
