# =============================================================================
# Truth-aware helpers of plots_synthetic.py (pure numpy / pandas: no pyCATHY import, unit-tested)
# =============================================================================
# df_DA layout (see nb5_plot_helpers.py): LONG table, one row per (cycle 'time', member 'Ensemble_nb', node);
# the node id is the position inside each (time, member, OL) block. Cycle c = state at the END of window c =
# atm_times[c]. truth['sw'] / ['psi'] are [time, node] arrays at truth['sw_time'] (nb6).
import json
import os
import re

import numpy as np
import pandas as pd

HOUR = 3600.0
STATE_COL = {'sw': 'sw_bef_update_', 'psi': 'psi_bef_update'}
SOIL_COLUMN = {'SATCELL_VG': 'VGPSATCELL', 'n_VG': 'VGNCELL', 'thetar_VG': 'VGRMCCELL', 'Ks': 'PERMX'}


# -----------------------------------------------------------------------------
# runs
# -----------------------------------------------------------------------------
def collect_runs(sim_dir, prefix='SYN_'):
    """Table of the runs with a run_meta.json (written by nb7): one row per run."""
    rows = []
    for d in sorted(os.listdir(sim_dir)) if os.path.isdir(sim_dir) else []:
        f = os.path.join(sim_dir, d, 'run_meta.json')
        if not (d.startswith(prefix) and d.endswith('_meta') and os.path.isfile(f)):
            continue
        with open(f) as fh:
            m = json.load(fh)
        rows.append(dict(run=m.get('run', d[:-5]), case=m.get('case'), scenario=m.get('scenario'),
                         da={'pf': 'pf', 'esmda': 'esmda'}.get(str(m.get('da_type')), 'enkf'), nens=m.get('nens'), seed=m.get('seed'),
                         with_da=m.get('with_da'), prior_bias=m.get('prior_bias'), status=m.get('status'),
                         n_par=m.get('n_updated_parameters', 0), truth=m.get('truth')))
    return pd.DataFrame(rows)


def pair_open_loop(runs):
    """For every DA run the name of its open-loop run (same case/scenario/nens/seed/bias, with_da False), or None."""
    key = ['case', 'scenario', 'nens', 'seed', 'prior_bias']
    ol = runs[runs['with_da'] == False]  # noqa: E712
    lut = {tuple(r[k] for k in key): r['run'] for _, r in ol.iterrows()}
    return {r['run']: lut.get(tuple(r[k] for k in key)) for _, r in runs[runs['with_da'] == True].iterrows()}  # noqa: E712


# -----------------------------------------------------------------------------
# df_DA -> arrays
# -----------------------------------------------------------------------------
def ensemble_array(df, col, ol=False, nodes=None, mask_rejected=True):
    """df_DA -> (cycles, members, array[n_cycle, n_member, n_node]).

    ol selects the open-loop rows; nodes (optional) = node ids to keep; rejected members (DA rows only) are NaN.
    """
    keys = ['time', 'Ensemble_nb', 'OL']
    use = [c for c in keys + ['rejected', 'idnode', col] if c in df.columns]
    d = df[use].copy()
    d['OL'] = d['OL'].astype(bool)
    d = d[d['OL'] == bool(ol)]
    assert len(d), f"no rows with OL={ol}"
    if 'idnode' not in d.columns:                       # legacy table: the node id is the position in the block
        d['idnode'] = d.groupby(['time', 'Ensemble_nb'], sort=False).cumcount().to_numpy()
    cycles = np.sort(d['time'].unique())
    members = np.sort(d['Ensemble_nb'].unique())
    if nodes is None:                                    # nodes present at EVERY cycle (node subsets: store_full='analysis')
        cnt = d.groupby('idnode')['time'].nunique()
        node_ids = np.sort(cnt.index[cnt == cycles.size].to_numpy())
    else:
        node_ids = np.asarray(sorted(set(map(int, nodes))))
    d = d[d['idnode'].isin(node_ids)]
    it = np.searchsorted(cycles, d['time'].to_numpy())
    im = np.searchsorted(members, d['Ensemble_nb'].to_numpy())
    inn = np.searchsorted(node_ids, d['idnode'].to_numpy())
    arr = np.full((cycles.size, members.size, node_ids.size), np.nan)
    v = d[col].to_numpy(dtype=float)
    if mask_rejected and 'rejected' in d.columns and not ol:
        v = np.where(d['rejected'].to_numpy() > 0, np.nan, v)
    arr[it, im, inn] = v
    return cycles, members, node_ids, arr


def ensemble_array_ds(ds, field, nodes=None):
    """xarray Dataset of the states (StateStore.to_xarray / <project>_states.nc) -> (cycles, members, node_ids, array)
    in the layout of ``ensemble_array``. ``field``: 'sw' | 'psi' | 'analysis' (psi/sw are the values BEFORE the update).
    Rejected members are already NaN."""
    ds = ds if nodes is None else ds.sel(node=np.asarray(sorted(set(map(int, nodes)))))
    a = ds[field].transpose('cycle', 'member', 'node').values.astype(float)
    return ds['cycle'].values, ds['member'].values, ds['node'].values, a


def cycles_to_seconds(cycles, atm_times):
    atm = np.asarray(atm_times, dtype=float)
    c = np.asarray(cycles)
    return atm[np.round(c).astype(int)] if c.size and c.max() <= len(atm) - 1 else c.astype(float)


def truth_on_times(truth_field, truth_time, t_s, nodes=None, tol=1.0):
    """truth [time, node] at the times t_s [s] (exact output times; error if missing) -> [len(t_s), n_node]."""
    tt = np.asarray(truth_time, dtype=float)
    idx = []
    for t in np.atleast_1d(t_s):
        i = int(np.argmin(np.abs(tt - t)))
        assert abs(tt[i] - t) <= tol, f"truth has no output at t={t} s (nearest {tt[i]})"
        idx.append(i)
    out = np.asarray(truth_field)[idx]
    return out if nodes is None else out[:, np.asarray(nodes, dtype=int)]


# -----------------------------------------------------------------------------
# scores (ens: [n_t, n_member, n_node], y: [n_t, n_node])
# -----------------------------------------------------------------------------
def crps_ensemble(ens, y):
    """CRPS of an ensemble (member axis 1), NaN members ignored: mean|x-y| - 0.5 E|X-X'| (sorted form)."""
    x = np.sort(np.asarray(ens, float), axis=1)             # NaN last
    m = np.sum(~np.isnan(x), axis=1)                         # members per (t, node)
    t1 = np.nansum(np.abs(x - y[:, None, :]), axis=1) / m
    i = np.arange(1, x.shape[1] + 1)[None, :, None]
    t2 = np.nansum((2 * i - m[:, None, :] - 1) * x, axis=1) / m ** 2
    return t1 - t2


def ensemble_scores(ens, y, lo=0.05, hi=0.95):
    """Per time: rmse of the ensemble mean, bias, spread (mean member sd), spread/rmse, 90 % coverage, CRPS,
    all averaged over nodes. Returns a DataFrame indexed by the time position."""
    mean = np.nanmean(ens, axis=1)
    err = mean - y
    sd = np.nanstd(ens, axis=1, ddof=1)
    q_lo, q_hi = np.nanquantile(ens, lo, axis=1), np.nanquantile(ens, hi, axis=1)
    rmse = np.sqrt(np.nanmean(err ** 2, axis=1))
    spread = np.sqrt(np.nanmean(sd ** 2, axis=1))
    return pd.DataFrame(dict(rmse=rmse, bias=np.nanmean(err, axis=1), spread=spread,
                             spread_skill=spread / np.where(rmse > 0, rmse, np.nan),
                             coverage=np.nanmean((y >= q_lo) & (y <= q_hi), axis=1),
                             crps=np.nanmean(crps_ensemble(ens, y), axis=1)))


def nodes_in_box(points, lo, hi):
    """Boolean mask of the nodes inside the box [lo, hi] (3 coordinates each, CATHY frame)."""
    p = np.asarray(points, float)
    return np.all((p >= np.asarray(lo)) & (p <= np.asarray(hi)), axis=1)


def align_on_time(sc_ol, t_target, tol_h=0.01):
    """Per-time scores of the open loop at the times of t_target (hours, nearest instant within tol_h, else NaN).

    The open loop is a separate run with its own number of cycles (finished, crashed or still running), so it must be
    matched by time, never by position."""
    s = sc_ol.sort_index()
    return s.reindex(pd.Index(np.asarray(t_target, float)), method='nearest', tolerance=tol_h)


def summarize_run(sc_da, sc_ol=None, t_from_h=None, t_h=None):
    """One row of the summary table from per-time scores (mean over cycles >= t_from_h).

    sc_da / sc_ol are indexed by time [h]. The open loop is matched to the DA cycles by time (align_on_time):
    DA cycles without an open-loop output at the same instant are left out of rmse_ol (NaN-mean)."""
    sel = np.ones(len(sc_da), bool) if t_from_h is None else np.asarray(t_h) >= t_from_h
    row = {k: float(sc_da[k][sel].mean()) for k in ('rmse', 'bias', 'spread_skill', 'coverage', 'crps')}
    row['rmse_final'] = float(sc_da['rmse'].iloc[-1])
    if sc_ol is not None:
        t_da = sc_da.index.values if t_h is None else np.asarray(t_h)
        ol = align_on_time(sc_ol, t_da)['rmse'].to_numpy()[sel]
        row['rmse_ol'] = float(np.nanmean(ol)) if np.isfinite(ol).any() else np.nan
        row['skill_vs_ol'] = 1.0 - row['rmse'] / row['rmse_ol'] if row['rmse_ol'] > 0 else np.nan
    return row


# -----------------------------------------------------------------------------
# parameters vs truth
# -----------------------------------------------------------------------------
def split_parm_name(name):
    m = re.match(r'^(.*?)(\d+)$', name)
    return (m.group(1), int(m.group(2))) if m else (name, None)


def truth_of_parameter(name, scenario, soil_table, archie_truth):
    """True value of a scenario parameter: Archie scalars from archie_truth {'a','m','n'}; layered parameters = mean
    of the truth soil table over the rows of the group (control 'material': all zones; 'zone_material': that zone).
    Returns None when it cannot be resolved."""
    if name in ('a_Archie', 'm_Archie', 'n_Archie'):
        return float(archie_truth[name[0]])
    fam, g = split_parm_name(name)
    if fam not in SOIL_COLUMN or g is None or fam not in scenario['per_name']:
        return None
    ctl = scenario['pert_control_name'][scenario['per_name'].index(fam)]
    col = soil_table[SOIL_COLUMN[fam]].astype(float)
    if ctl == 'zone_material':
        z, layers = scenario['zone_material_layers'][g]
        rows = [(int(z), l) for l in layers]
    elif ctl == 'material':
        zones = soil_table.index.get_level_values(0).unique()
        rows = [(z, l) for z in zones for l in scenario['material_layers'][g]]
    else:
        return None
    return float(np.nanmean(col.reindex(rows).to_numpy()))


def parm_evolution(entry):
    """dict_parm_pert[p] -> [n_cycle + 1, n_member] (prior, then each analysed ensemble)."""
    rows = [np.ravel(np.asarray(entry['ini_perturbation'], float))]
    ks = sorted(int(k[len('update_nb'):]) for k in entry if str(k).startswith('update_nb'))
    rows += [np.ravel(np.asarray(entry[f'update_nb{k}'], float)) for k in ks]
    n = max(len(r) for r in rows)
    return np.vstack([np.pad(r, (0, n - len(r)), constant_values=np.nan) for r in rows])


def parameter_table(dpp, names, scenario, soil_table, archie_truth):
    """Per parameter: truth, prior / final mean and sd, bias in sd (z), contraction 1 - sd_post/sd_prior,
    and 'truth_in_90' (truth inside the final 5-95 % range). Only parameters with a stored ensemble."""
    rows = []
    for p in names:
        if p not in dpp or not isinstance(dpp[p], dict) or 'ini_perturbation' not in dpp[p]:
            continue
        ev = parm_evolution(dpp[p])
        t = truth_of_parameter(p, scenario, soil_table, archie_truth)
        m0, s0 = np.nanmean(ev[0]), np.nanstd(ev[0], ddof=1)
        m1, s1 = np.nanmean(ev[-1]), np.nanstd(ev[-1], ddof=1)
        q = np.nanquantile(ev[-1], [0.05, 0.95])
        rows.append(dict(parameter=p, truth=t, prior_mean=m0, prior_sd=s0, final_mean=m1, final_sd=s1,
                         z_prior=(m0 - t) / s0 if t is not None and s0 > 0 else np.nan,
                         z_final=(m1 - t) / s1 if t is not None and s1 > 0 else np.nan,
                         err_prior=(m0 - t) if t is not None else np.nan,
                         err_final=(m1 - t) if t is not None else np.nan,
                         contraction=1 - s1 / s0 if s0 > 0 else np.nan,
                         truth_in_90=bool(t is not None and q[0] <= t <= q[1])))
    return pd.DataFrame(rows).set_index('parameter') if rows else pd.DataFrame()
