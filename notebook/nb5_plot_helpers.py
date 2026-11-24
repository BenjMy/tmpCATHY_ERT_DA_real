# =============================================================================
# Helpers for nb5_DA_ERT_soil3D_nc_plots.py  (no pyCATHY import: testable on its own)
# =============================================================================
# df_DA (pyCATHY.DA.cathy_DA._DA_df) is a LONG table with one row per
# (cycle, ensemble member, node):
#   time (= atmbc cycle counter c: the row set of cycle c is the state at the END of the model
#   window, i.e. at atm_times[c]), Ensemble_nb, psi_bef_update, sw_bef_update_ (state before the
#   analysis), analysis, OL (open-loop flag), rejected (1 = member rejected at that cycle).
# The node id is not stored: it is the position inside each (time, member, OL) block.
import os
import re
import numpy as np
import pandas as pd

HOUR = 3600.0
STATE_COLS = {'sw': 'sw_bef_update_', 'psi': 'psi_bef_update'}


# -----------------------------------------------------------------------------
# df_DA -> ensemble statistics
# -----------------------------------------------------------------------------
def reduce_df_DA(df, nodes, cols=('sw_bef_update_', 'psi_bef_update')):
    """Add the node id and keep only the requested nodes/columns (df_DA can have millions of rows)."""
    keys = ['time', 'Ensemble_nb', 'OL']
    d = df[[c for c in keys + ['rejected'] + list(cols) if c in df.columns]].copy()
    d['OL'] = d['OL'].astype(bool)
    d['idnode'] = d.groupby(keys, sort=False).cumcount().to_numpy()
    d = d[d['idnode'].isin(list(nodes))]
    return d.reset_index(drop=True)


def cycle_to_seconds(t, atm_times):
    """df_DA['time'] (cycle counter c -> atm_times[c]); already-in-seconds values are returned as is."""
    t = np.asarray(t, dtype=float)
    atm = np.asarray(atm_times, dtype=float)
    if t.size and np.nanmax(t) <= len(atm) - 1 and np.allclose(t, np.round(t)):
        return atm[np.round(t).astype(int)]
    return t


def ensemble_stats(df_red, node, col, atm_times, stabilize_s, ol=False):
    """Mean/std/min/max/q05/q50/q95 over the members, per cycle, at one node.

    Index = hours since the START OF IRRIGATION. Rejected members are excluded (DA run only).
    """
    sub = df_red[(df_red['idnode'] == node) & (df_red['OL'] == ol)]
    if 'rejected' in sub.columns and not ol:
        sub = sub[sub['rejected'] == 0]
    if sub.empty:
        return None
    g = sub.groupby('time')[col]
    out = pd.DataFrame({'mean': g.mean(), 'std': g.std(), 'min': g.min(), 'max': g.max(),
                        'q05': g.quantile(0.05), 'q50': g.median(), 'q95': g.quantile(0.95),
                        'n': g.count()})
    out.index = (cycle_to_seconds(out.index.values, atm_times) - stabilize_s) / HOUR
    out.index.name = 't_h'
    return out.sort_index()


def ensemble_members(df_red, node, col, atm_times, stabilize_s, ol=False):
    """Pivot time x member (thin member lines)."""
    sub = df_red[(df_red['idnode'] == node) & (df_red['OL'] == ol)]
    if 'rejected' in sub.columns and not ol:
        sub = sub[sub['rejected'] == 0]
    if sub.empty:
        return None
    pv = sub.pivot_table(index='time', columns='Ensemble_nb', values=col)
    pv.index = (cycle_to_seconds(pv.index.values, atm_times) - stabilize_s) / HOUR
    return pv.sort_index()


def acquisition_times_h(obs_stamps_s, atm_times, stabilize_s):
    """Observations are stamped with the START of their model window (nb5, section 4): they are
    compared with the state at the window END = next atmbc time. Returns hours since irrigation."""
    atm = np.asarray(atm_times, dtype=float)
    out = []
    for s in obs_stamps_s:
        i = int(np.argmin(np.abs(atm - s)))
        out.append((atm[min(i + 1, len(atm) - 1)] - stabilize_s) / HOUR)
    return np.asarray(out)


# -----------------------------------------------------------------------------
# soil file -> porosity at the sensor depths
# -----------------------------------------------------------------------------
def layer_of_depth(depth_m, zr_thick):
    """0-based layer (from the top) that contains a depth [m]."""
    edges = np.cumsum(np.r_[0.0, np.asarray(zr_thick, dtype=float)])
    return int(np.clip(np.searchsorted(edges, depth_m, side='right') - 1, 0, len(zr_thick) - 1))


def porosity_at_depths(soil_df, depths_m, zr_thick, zone=1):
    """POROS of the soil layer containing each depth.

    soil_df: first output of simu.read_inputs('soil') (index zone x layer, column POROS).
    Layers are taken by POSITION from the top (the file is written top layer first), so it does
    not matter whether the index labels are 0- or 1-based.
    """
    df = soil_df
    if isinstance(df.index, pd.MultiIndex):
        zones = df.index.get_level_values(0)
        zone = zone if zone in set(zones) else zones[0]
        df = df.xs(zone, level=0)
    phi = df['POROS'].to_numpy(dtype=float)
    if len(phi) != len(zr_thick):
        raise ValueError(f"soil file has {len(phi)} layers, ZR_THICK has {len(zr_thick)}")
    return {d: float(phi[layer_of_depth(d, zr_thick)]) for d in depths_m}


# -----------------------------------------------------------------------------
# soil moisture (SMC) observations
# -----------------------------------------------------------------------------
_TIME_NAMES = ('time', 't', 'datetime', 'date', 'timestamp', 'date_time', 'hours', 't_h', 'time_h')


def load_smc(path, depths_cm, columns=None, time_col=None, time_unit='h_irr', irr_start=None,
             stabilize_s=0.0, units='auto', bin_h=0.5, verbose=True):
    """Read soil moisture data -> DataFrame (index = hours since START OF IRRIGATION, one column
    per depth in cm, volumetric water content in m3/m3).

    path       : .csv/.txt/.dat (separator auto-detected), .xlsx/.xls, .parquet, .pkl
    columns    : {depth_cm: 'col' | ['col1', 'col2']} (replicates are averaged). None = auto: the
                 column whose name contains the depth as a number (e.g. 'VWC_15', 'theta 15cm').
    time_col   : name of the time column (None = a column called time/date/... or the first one)
    time_unit  : numeric time column: 'h_irr' | 'min_irr' | 's_irr' (since irrigation start),
                 'h_cathy' | 's_cathy' (CATHY t = 0, i.e. stabilize_s before irrigation).
                 Date/time strings need irr_start (e.g. '2024-05-14 09:00').
    units      : 'auto' (values > 1.5 are taken as %), 'm3m3' or 'percent'
    bin_h      : average the raw series in bins of this width [h] (None = raw)
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.csv', '.txt', '.dat'):
        df = pd.read_csv(path, sep=None, engine='python')
    elif ext in ('.xlsx', '.xls'):
        df = pd.read_excel(path)
    elif ext == '.parquet':
        df = pd.read_parquet(path)
    elif ext in ('.pkl', '.pickle'):
        df = pd.read_pickle(path)
    else:
        raise ValueError(f"unsupported SMC file type: {ext}")
    df.columns = [str(c).strip() for c in df.columns]

    if time_col is None:
        cand = [c for c in df.columns if c.lower() in _TIME_NAMES]
        time_col = cand[0] if cand else df.columns[0]
    tcol = df[time_col]
    if pd.api.types.is_numeric_dtype(tcol):
        v = tcol.to_numpy(dtype=float)
        conv = {'h_irr': v, 'min_irr': v / 60.0, 's_irr': v / HOUR,
                'h_cathy': v - stabilize_s / HOUR, 's_cathy': (v - stabilize_s) / HOUR}
        if time_unit not in conv:
            raise ValueError(f"time_unit must be one of {list(conv)}")
        t_h = conv[time_unit]
    else:
        dt = pd.to_datetime(tcol)
        if irr_start is None:
            raise ValueError("the SMC time column holds dates: set irr_start (irrigation start date/time)")
        t0 = pd.Timestamp(irr_start)
        if getattr(dt.dt, 'tz', None) is not None:
            dt = dt.dt.tz_localize(None)
        t_h = (dt - t0).dt.total_seconds().to_numpy() / HOUR

    out = {}
    for d in depths_cm:
        if columns is not None and d in columns:
            cols = columns[d] if isinstance(columns[d], (list, tuple)) else [columns[d]]
        elif columns is not None:
            continue
        else:
            pat = re.compile(rf'(?<![\d.]){re.escape(str(d))}(?![\d.])')
            cols = [c for c in df.columns if c != time_col and pat.search(c)]
        cols = [c for c in cols if c in df.columns]
        if not cols:
            if verbose:
                print(f"  SMC: no column found for {d} cm (columns: {list(df.columns)})")
            continue
        out[d] = df[cols].apply(pd.to_numeric, errors='coerce').mean(axis=1).to_numpy()
        if verbose:
            print(f"  SMC {d:>3} cm <- {cols}")
    if not out:
        raise ValueError("no SMC depth column recognised: pass columns={depth_cm: 'column name'}")

    smc = pd.DataFrame(out, index=pd.Index(t_h, name='t_h')).sort_index()
    smc = smc[~smc.index.isna()]
    vmax = np.nanmax(smc.to_numpy())
    if units == 'percent' or (units == 'auto' and vmax > 1.5):
        smc = smc / 100.0
        if verbose:
            print(f"  SMC values read as % (max {vmax:.1f}) -> m3/m3")
    if bin_h:
        smc = smc.groupby(np.floor(smc.index.to_numpy() / bin_h) * bin_h + bin_h / 2).mean()
        smc.index.name = 't_h'
    return smc


def smc_to_saturation(smc, phi_by_depth_cm, verbose=True):
    """sw = theta / porosity (porosity from the soil file, per sensor depth)."""
    sw = smc.copy()
    for d in sw.columns:
        sw[d] = smc[d] / phi_by_depth_cm[d]
        nover = int((sw[d] > 1.0).sum())
        if verbose and nover:
            print(f"  SMC {d} cm: {nover} values above saturation (theta > porosity {phi_by_depth_cm[d]:.3f})")
    return sw


# -----------------------------------------------------------------------------
# skill of the ensemble against the (independent) SMC data
# -----------------------------------------------------------------------------
def smc_skill(stats_by_depth, sw_obs, phi_by_depth_cm, tol_h=0.75):
    """Bias / RMSE / r / coverage of the ensemble vs the SMC data, one row per depth.

    stats_by_depth: {depth_cm: ensemble_stats DataFrame}. Model values are taken at the model
    times within tol_h of an observation. Bias = model - obs. coverage = share of observations
    inside the 5-95 % ensemble band. RMSE is also given in m3/m3 (x porosity).
    """
    rows = []
    for d, st in stats_by_depth.items():
        if st is None or d not in sw_obs.columns:
            continue
        o = sw_obs[d].dropna()
        o = o[(o.index >= st.index.min()) & (o.index <= st.index.max())]
        if o.empty:
            continue
        mt = st.index.to_numpy()
        idx = np.abs(mt[None, :] - o.index.to_numpy()[:, None]).argmin(axis=1)
        keep = np.abs(mt[idx] - o.index.to_numpy()) <= tol_h
        if not keep.any():
            continue
        o, idx = o[keep], idx[keep]
        m = st['mean'].to_numpy()[idx]
        lo, hi = st['q05'].to_numpy()[idx], st['q95'].to_numpy()[idx]
        e = m - o.to_numpy()
        r = np.corrcoef(m, o.to_numpy())[0, 1] if len(o) > 2 and np.std(m) > 0 and np.std(o) > 0 else np.nan
        rows.append({'depth_cm': d, 'n_obs': len(o), 'bias_sw': e.mean(),
                     'rmse_sw': float(np.sqrt((e ** 2).mean())),
                     'rmse_vwc': float(np.sqrt((e ** 2).mean()) * phi_by_depth_cm[d]),
                     'r': r, 'coverage_5_95': float(((o.to_numpy() >= lo) & (o.to_numpy() <= hi)).mean())})
    return pd.DataFrame(rows).set_index('depth_cm') if rows else pd.DataFrame()


# -----------------------------------------------------------------------------
# plotting
# -----------------------------------------------------------------------------
def plot_ensemble(ax, st, color, label, members=None, band=('q05', 'q95'), ls='-'):
    """Median line + 5-95 % band (+ optional thin member lines) of one ensemble."""
    if st is None:
        return
    if members is not None:
        ax.plot(members.index, members.values, color=color, lw=0.4, alpha=0.15, zorder=1)
    ax.fill_between(st.index, st[band[0]], st[band[1]], color=color, alpha=0.25, lw=0, zorder=2)
    ax.plot(st.index, st['mean'], color=color, lw=1.6, ls=ls, label=label, zorder=3)


def plot_smc(ax, obs, sigma=None, label='SMC sensor'):
    """Soil-moisture observations on top of the ensemble (black dots, optional error bars)."""
    if obs is None or obs.dropna().empty:
        return
    o = obs.dropna()
    if sigma is not None:
        ax.errorbar(o.index, o.values, yerr=sigma, fmt='none', ecolor='k', elinewidth=0.6,
                    alpha=0.5, zorder=4)
    ax.plot(o.index, o.values, 'o', ms=3.2, mfc='k', mec='white', mew=0.4, label=label, zorder=5)


def mark_events(ax, ert_t_h, irr_t_h=0.0, label_ert=True):
    ax.axvline(irr_t_h, color='tab:blue', ls='--', lw=1.0, zorder=0,
               label='irrigation start' if label_ert else None)
    for k, t in enumerate(ert_t_h):
        ax.axvline(t, color='tab:red', ls=':', lw=1.0, zorder=0,
                   label='ERT survey' if (k == 0 and label_ert) else None)


def profile_at_time(df_red, nodes, t_h, col, atm_times, stabilize_s, ol=False, tol_h=1e-3):
    """Ensemble mean / 5-95 % over the members at the given nodes (a vertical column), at one time
    [h since irrigation start]. Returns a DataFrame indexed by node id (None if no data)."""
    th = (cycle_to_seconds(df_red['time'].to_numpy(), atm_times) - stabilize_s) / HOUR
    sel = df_red[(np.abs(th - t_h) <= tol_h) & (df_red['OL'] == ol) & df_red['idnode'].isin(list(nodes))]
    if 'rejected' in sel.columns and not ol:
        sel = sel[sel['rejected'] == 0]
    if sel.empty:
        return None
    g = sel.groupby('idnode')[col]
    return pd.DataFrame({'mean': g.mean(), 'q05': g.quantile(0.05), 'q95': g.quantile(0.95)})


# -----------------------------------------------------------------------------
# TDR / sensor CSV of nb3 (data/sensors/st1.csv | st5.csv, 3-row header, Port 1/2/3 = 15/30/50 cm)
# -----------------------------------------------------------------------------
def load_tdr_nb3(sensor_csv, site_params_json, exp_label, depths_cm=(15, 30, 50), bin_h=None, verbose=True):
    """Same parsing/time reference as nb3_compare_meas_vs_mod.py -> DataFrame (index = hours since
    START OF IRRIGATION, one column per depth [cm], VWC m3/m3). INFILTRATION_START[exp_label] is t = 0."""
    import json
    with open(site_params_json) as f:
        t0 = pd.Timestamp(json.load(f)['INFILTRATION_START'][exp_label])
    df = pd.read_csv(sensor_csv, header=[0, 1, 2])
    df.columns = ['_'.join([str(lv) for lv in col if lv is not None]).strip() for col in df.columns]
    df.columns = [c.replace('___', '').replace('__', '').strip() for c in df.columns]
    dt = pd.to_datetime(df.iloc[:, 0], dayfirst=True)
    t_h = (dt - t0).dt.total_seconds().to_numpy() / HOUR
    out = {}
    for port, d in zip(('Port 1', 'Port 2', 'Port 3'), depths_cm):
        cols = [c for c in df.columns if port in c and 'VWC' in c]
        if not cols:
            continue
        out[d] = pd.to_numeric(df[cols[0]], errors='coerce').to_numpy()
        if verbose:
            print(f"  TDR {d:>3} cm <- {cols[0]}")
    if not out:
        raise ValueError(f"no 'Port n ... VWC' columns in {sensor_csv}")
    tdr = pd.DataFrame(out, index=pd.Index(t_h, name='t_h')).sort_index()
    tdr = tdr[~tdr.index.isna()]
    if np.nanmax(tdr.to_numpy()) > 1.5:
        tdr = tdr / 100.0
    if bin_h:
        tdr = tdr.groupby(np.floor(tdr.index.to_numpy() / bin_h) * bin_h + bin_h / 2).mean()
        tdr.index.name = 't_h'
    return tdr


def plot_obs_vs_model(ax, stats, obs, phi, color='tab:green', tol_h=0.75):
    """1:1 scatter of ensemble mean vs TDR (saturation) with 5-95 % bars; returns n points."""
    if stats is None or obs is None:
        return 0
    o = obs.dropna()
    o = o[(o.index >= stats.index.min()) & (o.index <= stats.index.max())]
    mt = stats.index.to_numpy()
    idx = np.abs(mt[None, :] - o.index.to_numpy()[:, None]).argmin(axis=1)
    keep = np.abs(mt[idx] - o.index.to_numpy()) <= tol_h
    o, idx = o[keep], idx[keep]
    if o.empty:
        return 0
    m, lo, hi = stats['mean'].to_numpy()[idx], stats['q05'].to_numpy()[idx], stats['q95'].to_numpy()[idx]
    ax.errorbar(o.to_numpy(), m, yerr=[np.clip(m - lo, 0, None), np.clip(hi - m, 0, None)], fmt='o', ms=3, color=color, alpha=0.6, elinewidth=0.5)
    return len(o)


# =============================================================================
# Parameters, rejections (pure functions, tested in test_nb5_plot_helpers.py)
# =============================================================================
import json
import re as _re


def layer_groups(scenario, families=('SATCELL_VG', 'n_VG', 'thetar_VG')):
    """Group consecutive layers with identical nominal values in ALL layered families of a scenario
    entry (= one soil material). Returns [(first_layer, last_layer), ...]; no hand-typed
    material table: it follows the scenario, which nb5 checks against the soil table."""
    names = list(scenario['per_name'])
    noms = [np.asarray(scenario['per_nom'][names.index(f)], float) for f in families if f in names]
    noms = [n for n in noms if n.ndim == 1]
    if not noms:
        return []
    mat = np.vstack(noms)                           # (n_family, n_layers)
    groups, start = [], 0
    for i in range(1, mat.shape[1] + 1):
        if i == mat.shape[1] or not np.allclose(mat[:, i], mat[:, start]):
            groups.append((start, i - 1))
            start = i
    return groups


def split_parm_name(name):
    """'n_VG12' -> ('n_VG', 12); 'a_Archie' -> ('a_Archie', None)."""
    m = _re.match(r'^(.*?)(\d+)$', name)
    return (m.group(1), int(m.group(2))) if m else (name, None)


def parm_evolution(dpp_entry):
    """dict_parm_pert[p] -> array (n_cycles + 1, n_members): prior, then each analysed ensemble.
    Keys read: 'ini_perturbation', 'update_nb<k>' (k = 1, 2, ...), NaN = invalid member."""
    rows = [np.ravel(np.asarray(dpp_entry['ini_perturbation'], float))]
    ks = sorted(int(k[len('update_nb'):]) for k in dpp_entry if str(k).startswith('update_nb'))
    rows += [np.ravel(np.asarray(dpp_entry[f'update_nb{k}'], float)) for k in ks]
    n = max(len(r) for r in rows)
    return np.vstack([np.pad(r, (0, n - len(r)), constant_values=np.nan) for r in rows])


def parm_shift_table(dpp, names, max_shift_sd=3.0, min_spread_ratio=0.1):
    """Per parameter: prior mean/sd, last mean/sd, shift of the mean in PRIOR sd, spread ratio.
    flag = 'SHIFT' if |shift| > max_shift_sd, 'COLLAPSE' if spread ratio < min_spread_ratio
    (validation protocol V8). Parameters without a stored ensemble are skipped."""
    rows = []
    for p in names:
        if p not in dpp or not isinstance(dpp[p], dict) or 'ini_perturbation' not in dpp[p]:
            continue
        ev = parm_evolution(dpp[p])
        m0, s0 = np.nanmean(ev[0]), np.nanstd(ev[0], ddof=1)
        m1, s1 = np.nanmean(ev[-1]), np.nanstd(ev[-1], ddof=1)
        shift = (m1 - m0) / s0 if s0 > 0 else np.nan
        ratio = s1 / s0 if s0 > 0 else np.nan
        flag = ('SHIFT' if abs(shift) > max_shift_sd else '') + \
               ('COLLAPSE' if ratio < min_spread_ratio else '')
        rows.append(dict(parameter=p, family=split_parm_name(p)[0], layer=split_parm_name(p)[1],
                         prior_mean=m0, prior_sd=s0, last_mean=m1, last_sd=s1,
                         shift_prior_sd=shift, spread_ratio=ratio, n_updates=ev.shape[0] - 1, flag=flag))
    return pd.DataFrame(rows).set_index('parameter') if rows else pd.DataFrame()


def read_rejections(path):
    """<project>_rejections.json -> DataFrame (cycle, stage, member, reason); empty if absent."""
    cols = ['cycle', 'stage', 'member', 'reason']
    if not path or not os.path.isfile(path):
        return pd.DataFrame(columns=cols)
    with open(path) as f:
        data = json.load(f)
    return pd.DataFrame(data, columns=cols) if data else pd.DataFrame(columns=cols)


def rejections_per_cycle(rej, n_members):
    """Share of rejected members per (cycle, stage); a member counted once per cycle and stage."""
    if rej.empty:
        return pd.DataFrame(columns=['after_run', 'after_analysis'])
    t = (rej.drop_duplicates(['cycle', 'stage', 'member']).groupby(['cycle', 'stage']).size()
         .unstack('stage', fill_value=0) / float(n_members))
    for c in ('after_run', 'after_analysis'):
        if c not in t:
            t[c] = 0.0
    return t[['after_run', 'after_analysis']]
