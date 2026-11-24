# =============================================================================
# Filter inbreeding of a particle-filter run   (todo.md, section 5)
# =============================================================================
#   python pf_inbreeding.py simulations/SYN_NC_full_upd_pf_N32_s1_unbiased/<project>_states   # folder of cycle_*.npz
#   python pf_inbreeding.py <states_dir> --csv inbreeding.csv --plot inbreeding.png
#
# Inbreeding = the ensemble descends from a few particles: the PF resamples (copies) the heavy particles and
# drops the light ones, so after k cycles the members share ancestors. The resampling indices are not stored,
# but they can be recovered: with `use_enkf_update=False` and `pf_jitter_std_state=0` (nb7 settings) the
# analysed state of member j is a COPY of the pre-update state of its parent p(j)
#       analysis_c[:, j] == psi_c[:, p(j)]
# and member j of cycle c+1 is the forward run of analysis_c[:, j]. Chaining the parents gives, at each cycle,
# the founder (member of the first cycle) every member descends from.
#
# Reported per analysis cycle (CSV columns):
#   resampled         some member is not its own parent
#   n_unique_before   distinct pre-update states        (< N: clones that stayed identical, e.g. no perturbed parameter)
#   n_unique_after    distinct analysed states          (< N: the resampling duplicated particles)
#   n_parents         distinct parents                  (N - n_parents members were dropped)
#   max_offspring     largest family produced this cycle
#   n_founders        distinct founders left in the ensemble  (N at the start; 1 = full collapse)
#   founders_ess      1 / sum(share^2) of the founders  (effective number of founders)
#   spread_before / spread_after   mean over nodes of the std over members of psi (loss of spread by the analysis)
# Verdict: INBREEDING when the final n_founders <= N / 4 (override with --frac).
# Limits: if two pre-update states are bit-identical the parent is ambiguous (counted in `ties`; the lowest
# index is taken, `n_unique_*` do not depend on it). With `pf_jitter_std_state > 0` or `use_enkf_update=True`
# the copy is no longer exact: the nearest column is used and `match_err` shows how far it is.
# =============================================================================
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd


def load_cycles(states_dir):
    """[(cycle, psi[n_nodes, N], analysis[n_nodes, N], rejected[N], is_analysis)] from cycle_*.npz (StateStore chunks)."""
    files = sorted(glob.glob(os.path.join(states_dir, 'cycle_*.npz')),
                   key=lambda f: int(re.findall(r'cycle_(\d+)', f)[-1]))
    if not files:
        raise FileNotFoundError(f'no cycle_*.npz in {states_dir}')
    out = []
    for f in files:
        with np.load(f) as z:
            out.append((int(z['cycle']), z['psi'].astype(np.float64), z['analysis'].astype(np.float64),
                        z['rejected'], bool(z['is_analysis'])))
    return out


def n_unique_columns(a, tol=0.0):
    """Number of distinct columns of a [n_nodes, N] (exact when tol = 0)."""
    if tol == 0:
        return np.unique(np.ascontiguousarray(a.T), axis=0).shape[0]
    keep = []
    for j in range(a.shape[1]):
        if not any(np.max(np.abs(a[:, j] - a[:, k])) <= tol for k in keep):
            keep.append(j)
    return len(keep)


def parents(psi, analysis):
    """parent[j] = column of `psi` closest to analysis[:, j]; also (max distance, number of ambiguous ties)."""
    n = psi.shape[1]
    d = np.empty((n, n))
    for j in range(n):
        d[j] = np.max(np.abs(psi - analysis[:, [j]]), axis=0)          # Chebyshev distance to every candidate parent
    par = d.argmin(axis=1)
    best = d[np.arange(n), par]
    ties = int(np.sum((d <= best[:, None] + 1e-12).sum(axis=1) > 1))
    return par, float(best.max()), ties


def genealogy(cycles, frac=0.25):
    """DataFrame, one row per cycle that has an analysis, and the verdict."""
    n = cycles[0][1].shape[1]
    founder = np.arange(n)                                             # founder of each forecast member
    rows = []
    for c, psi, ana, rej, is_an in cycles:
        if not is_an:
            continue
        par, err, ties = parents(psi, ana)
        resampled = bool(np.any(par != np.arange(n)))
        founder = founder[par]                                         # founder of each ANALYSED member
        cnt = np.bincount(founder, minlength=n) / n
        rows.append(dict(
            cycle=c, resampled=resampled, n_unique_before=n_unique_columns(psi), n_unique_after=n_unique_columns(ana),
            n_parents=len(np.unique(par)), max_offspring=int(np.bincount(par, minlength=n).max()),
            n_founders=int(np.count_nonzero(cnt)), founders_ess=float(1.0 / np.sum(cnt ** 2)),
            spread_before=float(np.nanmean(np.std(psi, axis=1))), spread_after=float(np.nanmean(np.std(ana, axis=1))),
            match_err=err, ties=ties))
        # member j of the next cycle is the forward run of analysed member j: `founder` carries over unchanged
    df = pd.DataFrame(rows)
    final = int(df.n_founders.iloc[-1]) if len(df) else n
    verdict = dict(nens=n, n_analysis_cycles=len(df), resampled_cycles=int(df.resampled.sum()) if len(df) else 0,
                   final_founders=final, final_founders_ess=float(df.founders_ess.iloc[-1]) if len(df) else float(n),
                   inbreeding=bool(final <= frac * n))
    return df, verdict


def plot(df, n, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(df.cycle, df.n_founders, 'o-', label='distinct founders')
    ax[0].plot(df.cycle, df.founders_ess, 's--', label='founders ESS')
    ax[0].plot(df.cycle, df.n_unique_after, '^:', label='distinct analysed states')
    ax[0].axhline(n, color='gray', lw=.5)
    ax[0].set(xlabel='cycle', ylabel='members', title='genealogy', ylim=(0, n * 1.05)); ax[0].legend(fontsize=8)
    ax[1].plot(df.cycle, df.spread_before, 'o-', label='before analysis')
    ax[1].plot(df.cycle, df.spread_after, 's-', label='after analysis')
    ax[1].set(xlabel='cycle', ylabel='mean std of psi [m]', title='spread'); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description='Particle-filter genealogy (inbreeding) from the cycle_*.npz state chunks')
    ap.add_argument('states_dir', help='folder with cycle_*.npz (<project>_states/)')
    ap.add_argument('--frac', type=float, default=0.25, help='INBREEDING if final founders <= frac x N (default .25)')
    ap.add_argument('--csv'); ap.add_argument('--plot')
    a = ap.parse_args(argv)
    df, v = genealogy(load_cycles(a.states_dir), a.frac)
    pd.set_option('display.width', 200)
    print(df.to_string(index=False, float_format=lambda x: f'{x:.4g}'))
    print(f"\nN = {v['nens']}   analysis cycles = {v['n_analysis_cycles']}   resampled at {v['resampled_cycles']}   "
          f"founders left = {v['final_founders']} (ESS {v['final_founders_ess']:.1f})   ->  "
          f"{'INBREEDING' if v['inbreeding'] else 'ok'}")
    if df.match_err.max() > 1e-5:
        print(f"note: analysed states are not exact copies (max distance {df.match_err.max():.2g}): jitter / EnKF update on")
    if a.csv: df.to_csv(a.csv, index=False)
    if a.plot: plot(df, v['nens'], a.plot)
    return 1 if v['inbreeding'] else 0


if __name__ == '__main__':
    sys.exit(main())
