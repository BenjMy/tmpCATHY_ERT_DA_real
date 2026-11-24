"""pf_inbreeding: recover a KNOWN PF genealogy from cycle_*.npz chunks."""
import numpy as np

import pf_inbreeding as PI

N, NODES = 8, 40


def _write_run(tmp_path, parents_per_cycle, seed=0):
    """Simulate a PF: forecast psi (members drift apart by their own parameters), resample with the given
    parents, analysis = copy of the parent's psi. Returns the true founder count after each cycle."""
    rng = np.random.default_rng(seed)
    drift = rng.normal(0, 0.05, (NODES, N))            # member-specific 'parameters': copies diverge afterwards
    psi = rng.normal(-1.5, 0.5, (NODES, N))
    founder, truth = np.arange(N), []
    for c, par in enumerate(parents_per_cycle):
        par = np.asarray(par)
        ana = psi[:, par]
        np.savez(tmp_path / f'cycle_{c:05d}.npz', cycle=c, node_ids=np.arange(NODES), is_analysis=True,
                 rejected=np.zeros(N, np.int8), nens=N, psi=psi.astype(np.float32), sw=psi.astype(np.float32),
                 analysis=ana.astype(np.float32))
        founder = founder[par]
        truth.append(len(np.unique(founder)))
        psi = ana + drift                                # forward run of the analysed state (member j keeps its index)
    return truth


def test_recovers_known_genealogy(tmp_path):
    rng = np.random.default_rng(3)
    parents = [rng.choice(N, N, p=np.r_[0.4, 0.3, np.full(N - 2, 0.3 / (N - 2))]) for _ in range(6)]
    truth = _write_run(tmp_path, parents)
    df, v = PI.genealogy(PI.load_cycles(str(tmp_path)))
    assert list(df.n_founders) == truth
    assert df.n_founders.is_monotonic_decreasing
    assert df.n_parents.iloc[0] == len(np.unique(parents[0]))
    assert df.match_err.max() < 1e-6 and df.ties.sum() == 0


def test_no_resampling_keeps_every_founder(tmp_path):
    _write_run(tmp_path, [np.arange(N)] * 4)
    df, v = PI.genealogy(PI.load_cycles(str(tmp_path)))
    assert not df.resampled.any() and (df.n_founders == N).all() and not v['inbreeding']
    assert (df.n_unique_after == N).all()


def test_full_collapse_is_flagged_and_exit_code(tmp_path):
    _write_run(tmp_path, [np.zeros(N, int)] * 3)         # every cycle: all members copy member 0
    df, v = PI.genealogy(PI.load_cycles(str(tmp_path)))
    assert df.n_unique_after.iloc[0] == 1 and v['final_founders'] == 1 and v['inbreeding']
    assert PI.main([str(tmp_path)]) == 1


def test_non_analysis_cycles_are_skipped(tmp_path):
    _write_run(tmp_path, [np.arange(N)] * 2)
    z = dict(np.load(tmp_path / 'cycle_00001.npz')); z['is_analysis'] = False
    np.savez(tmp_path / 'cycle_00001.npz', **z)
    df, _ = PI.genealogy(PI.load_cycles(str(tmp_path)))
    assert list(df.cycle) == [0]
