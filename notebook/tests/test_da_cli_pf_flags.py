"""da_ert_synthetic.py exposes the PF settings as CLI flags whose defaults keep the former hard-coded values.
Checked on the source (the script needs pyCATHY and CATHY to run)."""
import ast
import os

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'da_ert_synthetic.py')


def _tree():
    with open(SRC) as fh:
        return ast.parse(fh.read())


def _flag_defaults():
    out = {}
    for n in ast.walk(_tree()):
        if isinstance(n, ast.Call) and getattr(n.func, 'attr', '') == 'add_argument' and n.args:
            kw = {k.arg: k.value for k in n.keywords}
            if 'default' in kw and isinstance(kw['default'], ast.Constant):
                out[n.args[0].value] = kw['default'].value
    return out


def test_defaults_unchanged():
    d = _flag_defaults()
    assert d['--pf-jitter-std-state'] == 0.0
    assert d['--pf-resample-threshold'] == 1.0


def test_flags_reach_analysis_config():
    for n in ast.walk(_tree()):
        if isinstance(n, ast.Call) and getattr(n.func, 'id', '') == 'AnalysisConfig':
            kw = {k.arg: ast.unparse(k.value) for k in n.keywords}
            assert kw['pf_jitter_std_state'] == 'args.pf_jitter_std_state'
            assert kw['pf_resample_threshold'] == 'args.pf_resample_threshold'
            return
    raise AssertionError('AnalysisConfig call not found')
