"""data_path.toml loader: defaults, relative/absolute resolution, overrides, unknown keys."""
import os
import pytest
import project_paths as PP


def _write(tmp_path, text):
    f = tmp_path / 'data_path.toml'
    f.write_text(text)
    return str(f)


def test_repo_toml_reproduces_legacy_layout():
    P = PP.ProjectPaths()
    nb = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    legacy_root = os.path.abspath(os.path.join(nb, '..'))
    assert P.root == legacy_root
    assert P.simulations == os.path.join(legacy_root, 'simulations')
    assert P.mesh_forward == os.path.join(legacy_root, 'data', '3d', 'mesh_Luca_Agripolis3d', 'meshForward.vtk')
    assert P.cache == os.path.join(legacy_root, 'notebook', 'cache')
    assert P.figures_real == os.path.join(P.simulations, 'figures_nb5')
    assert P.root_str.endswith(os.sep) and P.simulations_str.endswith(os.sep)
    assert P.is_root('../') or os.getcwd() != nb      # '../' is the root when run from notebook/


def test_relative_to_toml_and_absolute(tmp_path):
    f = _write(tmp_path, '[project]\nroot = "proj"\n[outputs]\nsimulations = "/abs/sims"\ncache = "c"\n'
                         '[inputs]\nmesh_forward = "/abs/m.vtk"\n')
    P = PP.ProjectPaths(f)
    assert P.root == str(tmp_path / 'proj')
    assert P.simulations == '/abs/sims'
    assert P.cache == str(tmp_path / 'proj' / 'c')
    assert P.mesh_forward == '/abs/m.vtk'
    assert P.synthetic_case_dir('NC') == str(tmp_path / 'proj' / 'data' / '3d' / 'synthetic' / 'NC')


def test_derived_figures_follow_simulations(tmp_path):
    P = PP.ProjectPaths(_write(tmp_path, '[outputs]\nsimulations = "S"\n'))
    assert P.figures_real == os.path.join(P.simulations, 'figures_nb5')
    assert P.figures_synth == os.path.join(P.simulations, 'figures_nb8')


def test_env_var_selects_toml(tmp_path, monkeypatch):
    f = _write(tmp_path, '[outputs]\nsimulations = "elsewhere"\n')
    monkeypatch.setenv('DATA_PATH_TOML', f)
    assert PP.get(reload=True).simulations == str(tmp_path / 'elsewhere')
    monkeypatch.delenv('DATA_PATH_TOML')
    PP.get(reload=True)


def test_unknown_key_and_missing_file(tmp_path):
    with pytest.raises(AssertionError):
        PP.ProjectPaths(_write(tmp_path, '[outputs]\nsimulatons = "typo"\n'))
    with pytest.raises(FileNotFoundError):
        PP.ProjectPaths(str(tmp_path / 'nope.toml'))


def test_pycathy_must_contain_package(tmp_path):
    P = PP.ProjectPaths(_write(tmp_path, '[code]\npycathy = "nowhere"\n'))
    with pytest.raises(AssertionError):
        P.setup_sys_path()
