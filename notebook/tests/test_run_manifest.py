"""Run manifest: git state (branch, dirty, patch), pyCATHY discovery, CATHY executable, inputs, tracker hook."""
import json
import os
import subprocess
import sys
import textwrap

import pytest

import run_manifest as RM
import run_registry as RR
import project_paths as PP

GIT = ['git', '-c', 'user.name=T', '-c', 'user.email=t@example.org']


def _repo(path, branch='feat/x', remote='https://user:secret@github.com/o/r.git'):
    os.makedirs(path, exist_ok=True)
    subprocess.run(['git', 'init', '-q', '-b', branch, str(path)], check=True)
    (path / 'a.txt').write_text('one\n')
    subprocess.run([*GIT, '-C', str(path), 'add', '.'], check=True)
    subprocess.run([*GIT, '-C', str(path), 'commit', '-q', '-m', 'first'], check=True)
    if remote:
        subprocess.run(['git', '-C', str(path), 'remote', 'add', 'origin', remote], check=True)
    return path


def test_git_info_clean_dirty_and_credentials(tmp_path):
    r = _repo(tmp_path / 'r')
    g = RM.git_info(str(r))
    assert g['git'] and g['branch'] == 'feat/x' and len(g['commit']) == 40 and not g['dirty']
    assert g['remote_url'] == 'https://github.com/o/r.git'              # token stripped
    (r / 'a.txt').write_text('two\n')
    patch = str(tmp_path / 'meta' / 'uncommitted.patch')
    g = RM.git_info(str(r), patch)
    assert g['dirty'] and g['dirty_files'] == ['a.txt'] and g['diff_sha256']
    assert '+two' in open(patch).read()


def test_git_info_detached_and_not_a_repo(tmp_path):
    r = _repo(tmp_path / 'r')
    subprocess.run(['git', '-C', str(r), 'checkout', '-q', '--detach'], check=True)
    g = RM.git_info(str(r))
    assert g['detached'] and g['branch'] is None and 'feat/x' in ' '.join(g['branches_containing_head'])
    assert RM.git_info(str(tmp_path))['git'] is False
    assert RM.git_info(str(tmp_path / 'missing')) == {'git': False}


@pytest.fixture
def fake_pycathy(tmp_path, monkeypatch):
    clone = _repo(tmp_path / 'pycathy_wrapper', branch='dev-twin')
    pkg = clone / 'pyCATHY'
    pkg.mkdir()
    (pkg / '__init__.py').write_text("__version__ = '9.9.1'\n")
    subprocess.run([*GIT, '-C', str(clone), 'add', '.'], check=True)
    subprocess.run([*GIT, '-C', str(clone), 'commit', '-q', '-m', 'pkg'], check=True)
    monkeypatch.syspath_prepend(str(clone))
    sys.modules.pop('pyCATHY', None)
    yield clone
    sys.modules.pop('pyCATHY', None)


def test_pycathy_found_without_import_and_branch(fake_pycathy):
    info = RM.pycathy_info()
    assert info['found'] and info['version_source'] == '9.9.1'
    assert info['repo']['branch'] == 'dev-twin' and info['repo']['toplevel'] == str(fake_pycathy)
    assert 'pyCATHY' not in sys.modules                                 # discovery must not import the heavy package


def test_cathy_executable_sha_and_source(tmp_path):
    run = tmp_path / 'run'
    (run / 'bin').mkdir(parents=True)
    exe = run / 'bin' / 'cathy'
    exe.write_bytes(b'\x7fELFfake')
    exe.chmod(0o755)
    src = _repo(tmp_path / 'CATHY_src', branch='main', remote=None)
    c = RM.cathy_info(str(run), cathy_src=str(src))
    assert c['executable']['path'] == str(exe) and len(c['executable']['sha256']) == 64
    assert c['source_repo']['branch'] == 'main'
    assert RM.cathy_info(str(tmp_path / 'empty'))['executable'] is None


def test_inputs_fingerprint_and_copies(tmp_path):
    toml = tmp_path / 'data_path.toml'
    toml.write_text('[project]\nroot = "."\n')
    (tmp_path / 'data' / '3d').mkdir(parents=True)
    (tmp_path / 'data' / '3d' / 'f1.tldat').write_text('x' * 10)
    P = PP.ProjectPaths(str(toml))
    meta = tmp_path / 'meta'
    inp = RM.inputs_info(P, str(meta))
    assert inp['files']['data_path.toml']['sha256'] and (meta / 'inputs' / 'data_path.toml').is_file()
    assert inp['files']['site_params']['exists'] is False
    d1 = inp['folders']['data_3d']['digest']
    (tmp_path / 'data' / '3d' / 'f2.tldat').write_text('y')
    assert RM.inputs_info(P)['folders']['data_3d']['digest'] != d1       # an added file changes the digest


def test_collect_write_update_roundtrip(tmp_path, fake_pycathy):
    meta = str(tmp_path / 'run1_meta')
    f = RM.write('run1', meta, fields={'seed': 3, 'nens': 8, 'skipme': None})
    assert f == os.path.join(meta, 'manifest.json')
    m = RM.read(meta)
    assert m['schema'] == RM.SCHEMA and m['run'] == 'run1' and m['status'] == 'started'
    assert m['fields'] == {'seed': 3, 'nens': 8}
    assert m['cathy_da_ert']['git'] and m['pycathy']['version'] and m['python']['version']
    assert 'numpy' in {k.lower() for k in m['packages']}
    assert RM.update(meta, status='finished', duration_min='1.0')
    m2 = RM.read(meta)
    assert m2['status'] == 'finished' and m2['duration_min'] == '1.0' and m2['cathy_da_ert'] == m['cathy_da_ert']
    assert 'run1' in RM.summary(m2)


def test_compare_ignores_machine_noise_and_flags_versions():
    a = {'run': 'A', 'created': 't1', 'cathy_da_ert': {'commit': 'aaa', 'toplevel': '/x'},
         'pycathy': {'version': '1'}, 'platform': {'node': 'h1'}}
    b = {'run': 'B', 'created': 't2', 'cathy_da_ert': {'commit': 'bbb', 'toplevel': '/y'},
         'pycathy': {'version': '1'}, 'platform': {'node': 'h2'}}
    assert RM.compare(a, b) == {'cathy_da_ert.commit': ('aaa', 'bbb')}


def test_failure_never_raises(tmp_path, capsys):
    blocker = tmp_path / 'file'
    blocker.write_text('x')
    assert RM.write('r', str(blocker / 'sub')) is None                   # meta dir cannot be created
    assert 'could not write the manifest' in capsys.readouterr().err
    assert RM.update(str(tmp_path / 'nope')) is False


def test_tracker_writes_manifest_and_registry_columns(tmp_path):
    reg = str(tmp_path / 'reg' / 'r.csv')
    t = RR.track('runA', reg, kind='real', seed=7, folder=str(tmp_path / 'runA'))
    row = RR.read(reg)[0]
    meta = str(tmp_path / 'reg' / 'runA_meta')
    assert row['manifest'] == os.path.join(meta, 'manifest.json') and row['git_branch']
    assert RM.read(meta)['fields']['seed'] == 7
    t.finish('finished')
    assert RM.read(meta)['status'] == 'finished' and RR.read(reg)[0]['status'] == 'finished'


def test_tracker_without_manifest_and_crash_is_marked(tmp_path):
    reg = str(tmp_path / 'r.csv')
    RR.track('quiet', reg, manifest=False)
    assert not (tmp_path / 'quiet_meta').exists() and RR.read(reg)[0]['manifest'] == ''
    t = RR.track('boom', reg)
    t._at_exit()
    assert RM.read(str(tmp_path / 'boom_meta'))['status'] == 'failed'


def test_cli_show_and_compare(tmp_path, capsys):
    sims = tmp_path / 'sims'
    for r in ('A', 'B'):
        RM.write(r, str(sims / f'{r}_meta'), fields={'seed': 1 if r == 'A' else 2})
    assert RM.main(['show', 'A', '--sim-dir', str(sims)]) == 0
    assert 'CATHY DA ERT' in capsys.readouterr().out
    assert RM.main(['compare', 'A', 'B', '--sim-dir', str(sims)]) == 0
    assert 'fields.seed' in capsys.readouterr().out
    assert RM.main(['show', 'zzz', '--sim-dir', str(sims)]) == 1
