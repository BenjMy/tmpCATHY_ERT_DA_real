# Smoke tests of run_study.py (plan, tags, resume, failure tracking + registry, diagnostics call) and of the quick-test
# helpers of model_config (first_n_surveys, set_t_end). No pyCATHY needed: da_ert_synthetic.py / plots_synthetic.py are
# replaced by fake children and pyCATHY is stubbed when missing.
#   pytest notebook/tests/test_run_study.py     or     python notebook/tests/test_run_study.py
import hashlib
import json
import os
import pickle
import sys
import tempfile
import types

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NB = os.path.dirname(HERE)
sys.path.insert(0, NB)
import project_paths as PP  # noqa: E402
import run_registry as RR  # noqa: E402
import run_study as RS  # noqa: E402
from twin_names import run_name  # noqa: E402

FAKE_NB7 = '''
import json, os, sys
sys.path.insert(0, {nb!r})
from twin_names import run_name
a = sys.argv[1:]
g = lambda k: a[a.index(k) + 1]
tag = g('--tag') if '--tag' in a else ''
name = run_name(g('--case'), g('--scenario'), g('--da-type'), int(g('--nens')), int(g('--seed')), '--no-da' not in a,
                g('--prior-bias') == 'biased', tag)
sim = g('--sim-dir')
os.makedirs(os.path.join(sim, name + '_meta'), exist_ok=True)
fail = g('--scenario') == 'ic_Ks_upd' and g('--da-type') == 'pf'
meta = dict(run=name, status='failed' if fail else 'finished', n_updated_parameters=3, n_times=int(g('--n-times')) if '--n-times' in a else 0)
json.dump(meta, open(os.path.join(sim, name + '_meta', 'run_meta.json'), 'w'))
with open(os.path.join(sim, 'calls.txt'), 'a') as f:
    f.write(' '.join(a) + chr(10))
if fail:
    print('ValueError: fake failure'); sys.exit(1)
'''
FAKE_NB8 = '''
import os, sys
a = sys.argv[1:]
g = lambda k: a[a.index(k) + 1]
out = os.path.join(g('--sim-dir'), 'figures_nb8', g('--case'))
os.makedirs(out, exist_ok=True)
open(os.path.join(out, 'summary.csv'), 'w').write('run,rmse' + chr(10))
open(os.path.join(g('--sim-dir'), 'diag_calls.txt'), 'a').write(' '.join(a) + chr(10))
'''


def _fakes(tmp):
    n7, n8 = os.path.join(tmp, 'fake_nb7.py'), os.path.join(tmp, 'fake_nb8.py')
    open(n7, 'w').write(FAKE_NB7.format(nb=NB))
    open(n8, 'w').write(FAKE_NB8)
    return n7, n8


def _real_registry_hash():
    p = PP.get().run_registry
    return hashlib.sha1(open(p, 'rb').read()).hexdigest() if os.path.isfile(p) else None


def test_plan_quick_and_full():
    P = PP.get()
    a = RS.parse_args([])
    root = RS.resolve_root(a)
    assert root is None, "no --root: the children read data_path.toml"
    sim = RS.resolve_sim_dir(a, root)
    runs = RS.plan_runs(a, root, sim)
    assert len(runs) == 5 * (2 + 1), len(runs)                              # 5 scenarios x (EnKF + PF + ONE open loop)
    assert sim == P.simulations + '_n10' and RS.resolve_log_dir(a, root) == P.logs
    assert len({r['name'] for r in runs}) == len(runs), "run names must be unique"
    assert all(r['name'].endswith('_n10') for r in runs), "quick runs carry the tag n10 (separate registry rows)"
    ol = [r for r in runs if r['kind'] == 'open_loop']
    assert len(ol) == 5 and all('--no-da' in r['cmd'] and r['name'].endswith('_noDA_n10') for r in ol)
    assert all(r['cmd'][r['cmd'].index('--da-type') + 1] == 'enkf' for r in ol)
    assert all(r['cmd'][r['cmd'].index('--n-times') + 1] == '10' and r['cmd'][r['cmd'].index('--tag') + 1] == 'n10'
               for r in runs)
    assert all('--root' not in r['cmd'] for r in runs)
    assert all(r['cmd'][r['cmd'].index('--nens') + 1] == '32' for r in runs)
    assert all(os.path.basename(r['cmd'][1]) == 'da_ert_synthetic.py' for r in runs)
    assert [r['scenario'] for r in runs[::3]] == RS.LADDER

    f = RS.parse_args(['--full'])
    fruns = RS.plan_runs(f, None, RS.resolve_sim_dir(f, None))
    assert RS.resolve_sim_dir(f, None) == P.simulations
    assert all('--n-times' not in r['cmd'] and '--tag' not in r['cmd'] for r in fruns), "full runs = the plain call"
    assert fruns[0]['name'] == run_name('NC', 'ic', 'enkf', 100, 42, True, True), "same names as run_matrix.py"
    assert all(r['cmd'][r['cmd'].index('--nens') + 1] == '100' for r in fruns)

    r = RS.parse_args(['--root', '/tmp/otherproj'])
    assert RS.plan_runs(r, RS.resolve_root(r), '/x')[0]['cmd'][-2:] != [], "explicit root is passed on"
    assert '/tmp/otherproj/' in RS.plan_runs(r, RS.resolve_root(r), '/x')[0]['cmd']

    assert RS.normalize_scenarios(['NC_ic', 'ic', 'all'])[:2] == ['ic', 'ic_VGP']      # prefix stripped, deduplicated
    try:
        RS.normalize_scenarios(['nope'])
        raise RuntimeError('unknown scenario must be rejected')
    except AssertionError:
        pass


def test_pipeline_failure_tracking_registry_and_resume():
    reg0 = _real_registry_hash()
    with tempfile.TemporaryDirectory() as tmp:
        n7, n8 = _fakes(tmp)
        sim, reg, logs = os.path.join(tmp, 'sim'), os.path.join(tmp, 'reg', 'runs.csv'), os.path.join(tmp, 'logs')
        argv = ['--scenarios', 'ic', 'ic_Ks_upd', '--sim-dir', sim, '--log-dir', logs, '--registry', reg,
                '--nb7', n7, '--nb8', n8]
        assert RS.main(argv) == 1, "one run fails: exit code must be 1"
        rep = json.load(open(os.path.join(sim, 'study_report.json')))
        assert rep['n_failures'] == 1 and rep['n_times'] == 10
        failed = [r for r in rep['runs'] if r['result'] == 'FAILED']
        assert [r['run'] for r in failed] == ['SYN_NC_ic_Ks_upd_pf_N32_s42_n10'], failed
        assert 'fake failure' in failed[0]['reason'], failed[0]['reason']              # the reason comes from the log
        assert sum(r['result'] == 'ok' for r in rep['runs']) == 5, "a failed run must not stop the others"
        rows = {r['run']: r for r in RR.read(reg)}                                      # child died before registering
        assert rows['SYN_NC_ic_Ks_upd_pf_N32_s42_n10']['status'] == 'failed' and len(rows) == 1
        calls = open(os.path.join(sim, 'calls.txt')).read().splitlines()
        assert len(calls) == 6 and all('--n-times 10 --tag n10' in c for c in calls)
        assert rep['diagnostics'][0]['result'] == 'ok' and rep['diagnostics'][0]['summary_csv']
        dcalls = open(os.path.join(sim, 'diag_calls.txt')).read().splitlines()
        assert len(dcalls) == 1, "ONE diagnostics call per case"
        assert f'--sim-dir {sim}' in dcalls[0] and dcalls[0].split('--runs ')[1].split().__len__() == 6, dcalls[0]
        assert 'SYN_NC_ic_enkf_N32_s42_noDA_n10' in dcalls[0], "the open loops are scored too"

        assert RS.main(argv) == 1                                                       # resume: finished runs skipped
        rep = json.load(open(os.path.join(sim, 'study_report.json')))
        assert sum(r['result'] == 'skipped' for r in rep['runs']) == 5
        assert len(open(os.path.join(sim, 'calls.txt')).read().splitlines()) == 7       # only the failed run was retried
    assert _real_registry_hash() == reg0, "the project's real registry must not be touched by the tests"


def test_diag_only_without_runs_is_reported():
    with tempfile.TemporaryDirectory() as tmp:
        n7, n8 = _fakes(tmp)
        missing = os.path.join(tmp, 'missing')
        assert RS.main(['--stage', 'diag', '--sim-dir', missing, '--log-dir', os.path.join(tmp, 'l'), '--nb8', n8]) == 1
        assert RS.main(['--dry-run', '--sim-dir', missing, '--log-dir', os.path.join(tmp, 'l'), '--nb7', n7, '--nb8', n8]) == 0
        assert not os.path.exists(missing) and not os.path.exists(os.path.join(tmp, 'l')), "--dry-run creates nothing"


def _import_model_config():
    """model_config / nc_model import pyCATHY at module level: stub it when it is not installed."""
    added = []
    try:
        import pyCATHY.meshtools  # noqa: F401
    except ImportError:
        for n in ('pyCATHY', 'pyCATHY.meshtools'):
            sys.modules[n] = types.ModuleType(n)
            added.append(n)
        sys.modules['pyCATHY'].meshtools = sys.modules['pyCATHY.meshtools']
    import model_config as MC
    return MC, added


def test_quick_test_helpers_match_the_truth_time_vector():
    MC, added = _import_model_config()
    nc = MC._nc
    t_end0, case0 = nc.T_END, {k: c['t_end'] for k, c in MC.CASES.items()}
    try:
        times = {3: 5.0, 1: 1.0, 2: 3.0, 0: -0.4}                                       # shuffled on purpose
        assert list(MC.first_n_surveys(times, 2)) == [0, 1], "the EARLIEST surveys, whatever the dict order"
        assert MC.first_n_surveys(times, 0) == times and MC.first_n_surveys(times, None) == times
        assert MC.first_n_surveys(times, 99) == times

        truth_file = os.path.join(PP.get().synthetic_case_dir('NC'), 'truth_NC.pkl')
        if not os.path.isfile(truth_file):
            return
        with open(truth_file, 'rb') as f:
            tr = pickle.load(f)
        h = MC.HOUR
        stab = float(tr['stabilize_period'])
        all_h = {int(i): float(v['t_h_irr']) for i, v in tr['surveys'].items()}
        atm = np.asarray(tr['atm_times'], float)
        full = np.asarray(nc.build_times([t * h + stab for _, t in sorted(all_h.items())]))
        assert full.size == atm.size and np.allclose(full, atm), "baseline: full DA vector == truth vector"

        q = MC.first_n_surveys(all_h, 10)
        assert len(q) == 10 and max(q.values()) == sorted(all_h.values())[9]
        MC.set_t_end(0.0)
        MC.get_case('NC')                                                               # the per-case assert must still hold
        q_times = [t * h + stab for _, t in sorted(q.items())]
        vec = np.asarray(nc.build_times(q_times))
        assert vec.size < full.size / 5, "the quick model must be much shorter than the full one"
        assert vec[-1] == np.ceil(max(q_times) / h) * h, "the model stops at the last kept survey"
        pre = atm[atm <= vec[-1] + 1e-6]                                                # = the check added to the driver
        assert vec.size == pre.size and np.allclose(vec, pre), "quick DA vector must be a prefix of the truth vector"
        for t in q_times:
            assert np.min(np.abs(vec - t)) < 1e-3, "every kept survey has an exact model time"
    finally:
        nc.T_END = t_end0
        for k, v in case0.items():
            MC.CASES[k]['t_end'] = v
        for n in added + ['model_config', 'nc_model']:
            sys.modules.pop(n, None)


if __name__ == '__main__':
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            fn()
            print('ok', name)
    print('ALL OK')
