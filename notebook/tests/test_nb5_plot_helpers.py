import json
import numpy as np
import nb5_plot_helpers as H

SC = dict(per_name=['ic', 'SATCELL_VG', 'n_VG', 'thetar_VG'],
          per_nom=[-1.5, [.44, .44, .38, .4, .4], [1.4, 1.4, 1.42, 1.4, 1.4], [.04, .04, .065, .07, .07]])


def test_layer_groups_follow_the_scenario():
    assert H.layer_groups(SC) == [(0, 1), (2, 2), (3, 4)]
    assert H.layer_groups(dict(per_name=['ic'], per_nom=[-1.5])) == []


def test_split_name():
    assert H.split_parm_name('n_VG12') == ('n_VG', 12)
    assert H.split_parm_name('a_Archie') == ('a_Archie', None)


def _dpp():
    r = np.random.RandomState(0)
    ini = r.normal(1.4, 0.07, 32)
    return {'n_VG0': {'ini_perturbation': ini, 'update_nb1': ini + 0.01, 'update_nb2': ini * 0 + 1.4 + 0.5},
            'n_VG1': {'ini_perturbation': ini, 'update_nb1': ini * 0.1 + 1.26},
            'ic': {'x': 1}}


def test_evolution_and_shift_flags():
    d = _dpp()
    ev = H.parm_evolution(d['n_VG0'])
    assert ev.shape == (3, 32)
    t = H.parm_shift_table(d, ['n_VG0', 'n_VG1', 'ic', 'missing'])
    assert list(t.index) == ['n_VG0', 'n_VG1']
    assert 'SHIFT' in t.loc['n_VG0', 'flag'] and 'COLLAPSE' in t.loc['n_VG0', 'flag']   # all members equal
    assert t.loc['n_VG1', 'flag'] == 'COLLAPSE' and abs(t.loc['n_VG1', 'spread_ratio'] - 0.1) < 0.02


def test_rejections(tmp_path):
    p = tmp_path / 'r.json'
    p.write_text(json.dumps([dict(cycle=2, stage='after_run', member=1, reason='x'),
                             dict(cycle=2, stage='after_run', member=1, reason='x'),
                             dict(cycle=2, stage='after_analysis', member=3, reason='y')]))
    rej = H.read_rejections(p)
    t = H.rejections_per_cycle(rej, 10)
    assert t.loc[2, 'after_run'] == 0.1 and t.loc[2, 'after_analysis'] == 0.1
    assert H.read_rejections(tmp_path / 'none.json').empty
    assert H.rejections_per_cycle(H.read_rejections(None), 10).empty
