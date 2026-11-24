# =============================================================================
# Plot results of nb5_DA_ERT_soil3D_nc.py (ERT assimilation, field data)
# =============================================================================
# Plotting only. There is no synthetic truth here (real ERT surveys), so the comparison is:
# ensemble with vs without DA, soil-moisture (SMC) sensors as INDEPENDENT validation data
# (they are not assimilated), ERT misfit (RMS), parameter evolution and Archie scatter.
# Scenario and method (--scenario, --da enkf|pf|esmda) are command-line arguments; path2prj / PARAM_NAMES as in nb5.
#
# Needs nb5_plot_helpers.py next to this file (like nc_model.py).
#
# What is new (see nb5_plots_suggestions.md for the reasoning):
#   - soil-moisture observations over the saturation plots: sw_obs = theta / porosity, with the
#     porosity read from the project's soil file at each sensor depth
#   - same figure in volumetric water content (theta = sw * porosity)
#   - one x axis for everything: hours since the START OF IRRIGATION (t = 0 = irrigation)
#   - markers for irrigation and for the ERT surveys (at their acquisition time)
#   - median + 5-95 % band + thin member lines, DA and open loop on the same panel
#   - saturation profiles at the survey times (ensemble vs SMC)
#   - ensemble spread vs time, SMC skill table (bias / RMSE / r / coverage) saved as CSV
#   - every section is independent: one failing plot no longer stops the others
# =============================================================================
import os
import re
import sys
import argparse
import pickle
import glob
import traceback
import warnings
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyvista as pv

from pyCATHY.DA.cathy_DA import DA
import pyCATHY.plotters.cathy_plots as cplt
from nc_model import (X, Y, ZR_THICK, IRR_START_HOUR, SENSOR_DEPTHS, SENSOR_PHI, build_times)
import nb5_plot_helpers as H

mpl.rcParams.update({"axes.grid": True, "grid.color": "gray", "grid.linestyle": '--',
                     "grid.alpha": 0.4, "figure.constrained_layout.use": True,
                     "axes.spines.top": False, "axes.spines.right": False})

# %% Settings (must match nb5)
ROOT = '../'
path2prj = "../simulations/"
PARAM_NAMES = ['ic']
# Scenario, method and VGP control come from the command line (run_twin.sh real-plots passes them);
# the defaults are the former hard-coded values:
#   python nb5_DA_ERT_soil3D_nc_plots.py --scenario nc_ERT_ic_VGP_pert --da pf
_ap = argparse.ArgumentParser(description='figures of one real-data run (nb5 EnKF / PF, nb6 ES-MDA)')
_ap.add_argument('--scenario', default='nc_ERT_ic_VGP_upd', help="e.g. nc_ERT_ic_VGP_pert, nc_ERT_ic_VGP_Archie_upd")
_ap.add_argument('--da', default='enkf', choices=['enkf', 'pf', 'esmda'], help='method of the run to plot')
_ap.add_argument('--vgp-control', default='material', choices=['material', 'layers'], help='as in nb5/nb6')
_args, _ = _ap.parse_known_args()          # unknown arguments (Jupyter, wrapper) are ignored
SCENARIO = _args.scenario
DA_METHOD = _args.da
VGP_CONTROL = _args.vgp_control            # as in nb5 ('layers' or 'material')
RUN_TAG = SCENARIO + ('_mat' if VGP_CONTROL == 'material' else '')

# project folders written by the scripts (the open loop is the same for EnKF and PF: one run per scenario, EnKF name)
_DA_PREFIX = {'enkf': 'NC_ERT_', 'pf': 'NC_ERT_PF_', 'esmda': 'NC_ERT_ESMDA_'}[DA_METHOD]
_OL_PREFIX = 'NC_ERT_ESMDA_' if DA_METHOD == 'esmda' else 'NC_ERT_'
PRJ = {'withDA': _DA_PREFIX + 'withDA_' + RUN_TAG,
       'withoutDA': _OL_PREFIX + 'withoutDA_' + RUN_TAG,     # skipped if the folder does not exist
       }
print(f'plotting: method = {DA_METHOD} | scenario = {SCENARIO} | folders = {PRJ}')
COLORS = {'withDA': 'tab:green', 'withoutDA': 'tab:orange', 'open loop': 'tab:gray'}
HOUR = 3600.0
STABILIZE_PERIOD = IRR_START_HOUR * HOUR          # CATHY t = 0 is this long before the irrigation

# virtual probes: (x, y) in the CATHY frame and depths [cm] below the soil surface.
# Default = the SMC sensor depths of nc_model (None x/y = centre of the domain).
PROBE_XY = (None, None)
PROBE_DEPTHS_CM = [int(round(d * 100)) for d in SENSOR_DEPTHS]
PROFILE_MAX_DEPTH_M = 1.0                          # depth axis of the profile figure
SHOW_MEMBERS = True                                # thin line per ensemble member

# --- soil moisture observations (SMC) ------------------------------------------
# SMC_FILE = None  -> looks for one file named *smc*/*vwc*/*tdr*/*moisture* under ROOT/data
SMC_FILE = None                    # e.g. '../data/smc/agripolis_nc_smc.csv'
SMC_COLUMNS = None                 # None = auto (column names that contain the depth, e.g. 'VWC_15');
                                   # or {15: 'TDR_15', 30: ['A_30', 'B_30'], 50: 'TDR_50'} (lists are averaged)
SMC_TIME_COL = None                # None = 'time'/'date'/... column, else the first one
SMC_TIME_UNIT = 'h_irr'            # numeric time: 'h_irr' | 'min_irr' | 's_irr' | 'h_cathy' | 's_cathy'
SMC_IRR_START = None               # irrigation start, e.g. '2024-05-14 09:00' (needed if the time column holds dates)
SMC_UNITS = 'auto'                 # 'auto' (> 1.5 means %) | 'm3m3' | 'percent'
SMC_BIN_H = 0.5                    # average the raw sensor series in bins of this width [h] (None = raw)
# TDR sensors of nb3 (data/sensors/st5.csv = NC, st1.csv = HC; time reference = INFILTRATION_START of
# data/site_params.json). Used first if the file exists; otherwise the generic SMC loader above.
COMP_LEVEL = 'no_comp'             # as in nb5: 'no_comp' -> st5/NC | 'heavy_comp' -> st1/HC
TDR_FILE = os.path.join(ROOT, 'data', 'sensors', 'st5.csv' if COMP_LEVEL == 'no_comp' else 'st1.csv')
TDR_SITE_PARAMS = os.path.join(ROOT, 'data', 'site_params.json')
TDR_LABEL = 'NC' if COMP_LEVEL == 'no_comp' else 'HC'
TDR_BIN_H = None                   # None = raw TDR readings (dots); a number = bin average [h]
SMC_SIGMA_VWC = None               # sensor error [m3/m3] -> error bars (e.g. 0.02); None = no bars

BACKUP_PREFER = ('.pkl', '.parquet')   # use ('.parquet', '.pkl') for a faster load of the 250 MB df_DA

fig_dir = os.path.join(path2prj, 'figures_nb5', f'{DA_METHOD}_{RUN_TAG}')     # one folder per method and scenario
os.makedirs(fig_dir, exist_ok=True)


def savefig(fig, name):
    fig.savefig(os.path.join(fig_dir, name), dpi=300, bbox_inches='tight')


def section(title):
    print(f'\n--- {title} ' + '-' * max(3, 70 - len(title)))


def safe(fn):
    """Run a plotting section; report the error and carry on with the next one."""
    try:
        return fn()
    except Exception as e:
        print(f'  [{fn.__name__}] skipped: {type(e).__name__}: {e}')
        traceback.print_exc(limit=2)


def load_backup(simu):
    """Read the DA backup files of a project (df_DA, df_obs, dict_parm_pert, ET_DA_xr, ...).

    Uses DA.load_backup, which reads each `<project>_<name>.pkl/.parquet` file.
    (DA has no load_pickle_backup; load_parquet_backup only reads .parquet and so
    misses the .pkl-only dict_parm_pert and ET_DA_xr.)
    """
    res = simu.load_backup(prefer=BACKUP_PREFER)
    if 'df_DA' not in res:
        raise RuntimeError(f'no df_DA readable for {simu.project_name}; got: {list(res)}')
    missing = [k for k in ('df_obs', 'df_performance', 'dict_parm_pert', 'ET_DA_xr') if k not in res]
    if missing:
        print(f'  note: {simu.project_name}: not found {missing}')
    return res


# %% Reload DA backups
section('reload DA results')
simus, results = {}, {}
for mode, prj in PRJ.items():
    if not os.path.isdir(os.path.join(path2prj, prj)):
        print(f'skipping {mode}: {prj} not found')
        continue
    simus[mode] = DA(dirName=path2prj, prj_name=prj, notebook=True)
    try:
        results[mode] = load_backup(simus[mode])
    except Exception as e:
        print(f'skipping {mode}: {type(e).__name__}: {e}')
        simus.pop(mode, None)
        continue


def esmda_summary():
    """ES-MDA is a smoother: its result is <project>_esmda_result.pkl (misfit per forward run, kept members,
    parameter ensembles), not the cycle-by-cycle df_DA of the filters. Same figure as nb6 section 7."""
    pkl = os.path.join(path2prj, PRJ['withDA'], PRJ['withDA'] + '_esmda_result.pkl')
    if not os.path.isfile(pkl):
        print(f'ES-MDA result not found: {pkl}')
        return
    from pyCATHY.plotters.esmda_plots import plot_esmda_diagnostics, plot_misfit
    with open(pkl, 'rb') as f:
        res = pickle.load(f)
    mh = np.asarray(res.misfit_history, dtype=float)
    print('misfit per forward run:', np.round(mh, 2), '| members kept:', int(res.valid.sum()), '/', len(res.valid))
    pd.DataFrame({'run': np.arange(len(mh)), 'misfit': mh}).to_csv(os.path.join(fig_dir, 'esmda_misfit.csv'), index=False)
    if len(res.param_names) <= 12:
        fig = plot_esmda_diagnostics(res)
    else:
        fig = plt.figure(figsize=(5, 3.2))
        plot_misfit(res, ax=fig.add_subplot(111))
    savefig(fig, 'esmda_diagnostics.png')


if DA_METHOD == 'esmda':
    safe(esmda_summary)
    if not results:
        print('no cycle-by-cycle backup (df_DA) in the ES-MDA project: only the ES-MDA summary was plotted')
        print(f'figures in {os.path.abspath(fig_dir)}')
        sys.exit(0)
assert results, 'no DA results found'
ref_mode = 'withDA' if 'withDA' in results else next(iter(results))
ref = simus[ref_mode]

# %% Time axis: df_DA['time'] is the atmbc CYCLE counter (state at the END of window c = atm_times[c])
try:
    ref.read_inputs('atmbc')
    atm_times = np.asarray(ref.atmbc['time'], dtype=float)
except Exception as e:
    warnings.warn(f'atmbc times not readable ({type(e).__name__}); using the hourly vector of nc_model: '
                  'times of the ERT surveys on the sub-hourly grid will be approximate')
    atm_times = build_times()
n_cycles = int(max(r['df_DA']['time'].max() for r in results.values()))
if n_cycles > len(atm_times) - 1:
    warnings.warn(f'df_DA has {n_cycles} cycles but only {len(atm_times)} atmbc times: time axis unreliable')

# %% Observations / assimilation times
observations = results[ref_mode]['df_obs']
assimilation_times = observations.index.get_level_values(1).unique().to_list()   # window STARTS [s]
ert_t_h = H.acquisition_times_h(assimilation_times, atm_times, STABILIZE_PERIOD)  # acquisition [h since irrigation]
print('ERT surveys: DA stamp [h]', np.round(np.asarray(assimilation_times) / HOUR, 2),
      '-> acquisition [h since irrigation]', np.round(ert_t_h, 2))

# %% Probe nodes (nearest CATHY node to the sensor depths)
section('probe nodes')
mesh = pv.read(os.path.join(ref.workdir, ref.project_name, 'vtk', ref.project_name + '.vtk'))
x0 = X / 2 if PROBE_XY[0] is None else PROBE_XY[0]
y0 = Y / 2 if PROBE_XY[1] is None else PROBE_XY[1]
ztop = mesh.points[:, 2].max()
idx, coords = ref.find_nearest_node([[x0, y0, ztop - d / 100] for d in PROBE_DEPTHS_CM])
nodes_of_interest = [int(i) for i in np.atleast_1d(idx)]
node_depth_cm = [float(100 * (ztop - c[2])) for c in np.atleast_2d(coords)]
print('probe nodes (sensor depth -> node, node depth):',
      {d: (n, round(z, 1)) for d, n, z in zip(PROBE_DEPTHS_CM, nodes_of_interest, node_depth_cm)})

# vertical column of nodes under the probe position (profile figure)
edges = np.cumsum(np.r_[0.0, ZR_THICK])
edges = edges[edges <= PROFILE_MAX_DEPTH_M + 1e-9]
col_idx, col_coords = ref.find_nearest_node([[x0, y0, ztop - d] for d in edges])
col_nodes = [int(i) for i in np.atleast_1d(col_idx)]
col_depth = np.array([float(ztop - c[2]) for c in np.atleast_2d(col_coords)])

# matplotlib 3D instead of a pyvista off-screen render: that one segfaults through OpenGL/VTK and
# kills the kernel (a segfault cannot be caught by try/except)
def probe_nodes_figure():
    pts = np.asarray(mesh.points)
    fig = plt.figure(figsize=(5, 5))
    ax3 = fig.add_subplot(projection='3d')
    sub = pts[::max(1, len(pts) // 4000)]
    ax3.scatter(sub[:, 0], sub[:, 1], sub[:, 2], s=1, c='gray', alpha=0.15)
    ax3.scatter(*pts[nodes_of_interest].T, c='red', s=40, depthshade=False, label='probe nodes')
    ax3.set_xlabel('x (m)')
    ax3.set_ylabel('y (m)')
    ax3.set_zlabel('z (m)')
    ax3.legend(fontsize=8)
    savefig(fig, 'probe_nodes.png')


safe(probe_nodes_figure)

# %% Porosity at the sensor depths, from the soil file of the project
section('porosity (soil file)')
depths_m = [d / 100 for d in PROBE_DEPTHS_CM]
try:
    soil_df, _ = ref.read_inputs('soil')
    phi_m = H.porosity_at_depths(soil_df, depths_m, ZR_THICK)
    phi_src = 'soil file'
except Exception as e:
    warnings.warn(f'soil file not readable ({type(e).__name__}: {e}); using SENSOR_PHI of nc_model')
    phi_m = dict(zip(depths_m, SENSOR_PHI))
    phi_src = 'nc_model.SENSOR_PHI'
PHI = {d: phi_m[d / 100] for d in PROBE_DEPTHS_CM}
print(f'porosity ({phi_src}):', PHI)

# %% Soil moisture observations -> saturation
section('soil moisture observations')


def find_smc_file():
    hits = []
    for pat in ('*smc*', '*vwc*', '*tdr*', '*moisture*', '*soil_moist*'):
        hits += glob.glob(os.path.join(ROOT, 'data', '**', pat), recursive=True)
    hits = sorted({h for h in hits if os.path.splitext(h)[1].lower() in
                   ('.csv', '.txt', '.dat', '.xlsx', '.xls', '.parquet', '.pkl')})
    if len(hits) == 1:
        return hits[0]
    if hits:
        print('  several candidate SMC files, set SMC_FILE:', *hits, sep='\n    ')
    return None


smc_vwc = smc_sw = None
smc_path = SMC_FILE or find_smc_file()
if not SMC_FILE and os.path.isfile(TDR_FILE) and os.path.isfile(TDR_SITE_PARAMS):
    try:
        print('TDR file:', TDR_FILE)
        smc_vwc = H.load_tdr_nb3(TDR_FILE, TDR_SITE_PARAMS, TDR_LABEL, PROBE_DEPTHS_CM, bin_h=TDR_BIN_H)
        smc_sw = H.smc_to_saturation(smc_vwc, PHI)
        smc_path = None
        print(f'  {len(smc_vwc)} samples, t = {smc_vwc.index.min():.1f} .. {smc_vwc.index.max():.1f} h '
              'since irrigation start')
    except Exception as e:
        smc_vwc = smc_sw = None
        print(f'  TDR not loaded: {type(e).__name__}: {e}')
if smc_vwc is None and smc_path and os.path.isfile(smc_path):
    try:
        print('SMC file:', smc_path)
        smc_vwc = H.load_smc(smc_path, PROBE_DEPTHS_CM, columns=SMC_COLUMNS, time_col=SMC_TIME_COL,
                             time_unit=SMC_TIME_UNIT, irr_start=SMC_IRR_START,
                             stabilize_s=STABILIZE_PERIOD, units=SMC_UNITS, bin_h=SMC_BIN_H)
        smc_sw = H.smc_to_saturation(smc_vwc, PHI)
        print(f'  {len(smc_vwc)} samples, t = {smc_vwc.index.min():.1f} .. {smc_vwc.index.max():.1f} h '
              'since irrigation start')
    except Exception as e:
        smc_vwc = smc_sw = None
        print(f'  SMC not loaded: {type(e).__name__}: {e}')
elif smc_vwc is None:
    print('no TDR/SMC file (set TDR_FILE / SMC_FILE): saturation plots are drawn without observations')

# %% Ensemble statistics at the probes
section('ensemble statistics')
ALL_NODES = sorted(set(nodes_of_interest) | set(col_nodes))
red, stats, members = {}, {}, {}
for mode, res in results.items():
    red[mode] = H.reduce_df_DA(res['df_DA'], ALL_NODES)
    has_ol = bool(red[mode]['OL'].any())
    for var, col in H.STATE_COLS.items():
        for node, d in zip(nodes_of_interest, PROBE_DEPTHS_CM):
            stats[(mode, var, d, False)] = H.ensemble_stats(red[mode], node, col, atm_times, STABILIZE_PERIOD)
            members[(mode, var, d, False)] = (H.ensemble_members(red[mode], node, col, atm_times,
                                                                 STABILIZE_PERIOD) if SHOW_MEMBERS else None)
            if has_ol:
                stats[(mode, var, d, True)] = H.ensemble_stats(red[mode], node, col, atm_times,
                                                               STABILIZE_PERIOD, ol=True)
    n_rej = int((red[mode].drop_duplicates(['time', 'Ensemble_nb'])['rejected'] > 0).sum()) \
        if 'rejected' in red[mode] else 0
    print(f'{mode}: {len(red[mode]):,} rows kept, open loop rows: {has_ol}, '
          f'rejected member-cycles: {n_rej}')


def ts_figure(var, name, obs=None, scale=None, ylabel=None):
    """One panel per depth; every run overlaid; SMC observations (obs) on top."""
    fig, axs = plt.subplots(1, len(PROBE_DEPTHS_CM), figsize=(3.6 * len(PROBE_DEPTHS_CM), 3.4),
                            sharex=True, sharey=True, squeeze=False)
    for ax, d, nd in zip(axs[0], PROBE_DEPTHS_CM, node_depth_cm):
        f = 1.0 if scale is None else scale[d]
        for mode in results:
            for ol in (False, True):
                st = stats.get((mode, var, d, ol))
                if st is None:
                    continue
                st = st.copy()
                for c in ('mean', 'q05', 'q95'):
                    st[c] = st[c] * f
                mem = members.get((mode, var, d, ol))
                if mem is not None:
                    mem = mem * f
                lab = 'open loop' if ol else mode
                H.plot_ensemble(ax, st, COLORS.get(lab, 'tab:blue'), lab, members=mem,
                                ls='--' if ol else '-')
        if obs is not None and d in obs.columns:
            sig = None if SMC_SIGMA_VWC is None else (SMC_SIGMA_VWC / PHI[d] if scale is None else SMC_SIGMA_VWC)
            H.plot_smc(ax, obs[d], sigma=sig)
        H.mark_events(ax, ert_t_h)
        sub = f'{d} cm' if abs(nd - d) < 0.5 else f'{d} cm (node at {nd:.0f} cm)'
        ax.set_title(f'{sub}, porosity {PHI[d]:.3f}' if obs is not None else sub, fontsize=9)
        ax.set_xlabel('time since irrigation start (h)')
    axs[0][0].set_ylabel(ylabel or var)
    h, l = axs[0][0].get_legend_handles_labels()
    for a in axs[0][1:]:
        hh, ll = a.get_legend_handles_labels()
        h += [x for x, y in zip(hh, ll) if y not in l]
        l += [y for y in ll if y not in l]
    fig.legend(h, l, loc='outside lower center', ncol=len(l), fontsize=8, frameon=False)
    savefig(fig, name)
    return fig


# %% Saturation (and theta) with the soil-moisture observations on top
section('state time series')
safe(lambda: ts_figure('sw', 'sw_states_SMC.png', obs=smc_sw, ylabel='water saturation (-)'))
if smc_vwc is not None:
    safe(lambda: ts_figure('sw', 'vwc_states_SMC.png', obs=smc_vwc, scale=PHI,
                           ylabel=r'volumetric water content (m$^3$ m$^{-3}$)'))
safe(lambda: ts_figure('psi', 'psi_states.png', ylabel=r'pressure head $\psi$ (m)'))


# %% Saturation profiles at the ERT survey times
def profiles_figure():
    times = np.unique(np.round(ert_t_h, 3))
    fig, axs = plt.subplots(1, len(times), figsize=(2.6 * len(times), 4), sharex=True, sharey=True,
                            squeeze=False)
    for ax, t in zip(axs[0], times):
        for mode in results:
            for ol in (False, True):
                if ol and not red[mode]['OL'].any():
                    continue
                p = H.profile_at_time(red[mode], col_nodes, t, 'sw_bef_update_', atm_times,
                                      STABILIZE_PERIOD, ol=ol)
                if p is None:
                    continue
                p = p.reindex(col_nodes)
                lab = 'open loop' if ol else mode
                c = COLORS.get(lab, 'tab:blue')
                ax.fill_betweenx(col_depth, p['q05'], p['q95'], color=c, alpha=0.25, lw=0)
                ax.plot(p['mean'], col_depth, color=c, lw=1.6, ls='--' if ol else '-', label=lab)
        if smc_sw is not None:
            for d in PROBE_DEPTHS_CM:
                o = smc_sw[d].dropna()
                if o.empty or np.abs(o.index - t).min() > 0.75:
                    continue
                ax.plot(o.iloc[np.abs(o.index - t).argmin()], d / 100, 'o', ms=5, mfc='k', mec='white',
                        label='TDR' if d == PROBE_DEPTHS_CM[0] else None, zorder=5)
        ax.set_title(f'{t:+.0f} h', fontsize=9)
        ax.set_xlabel('saturation (-)')
    axs[0][0].set_ylim(PROFILE_MAX_DEPTH_M, 0)
    axs[0][0].set_ylabel('depth (m)')
    h, l = axs[0][-1].get_legend_handles_labels()
    if not l:
        h, l = axs[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc='outside lower center', ncol=max(len(l), 1), fontsize=8, frameon=False)
    savefig(fig, 'sw_profiles_at_ERT_times.png')


safe(profiles_figure)


# %% Ensemble spread (std of sw) vs time: does the DA reduce it?
def spread_figure():
    fig, axs = plt.subplots(1, len(PROBE_DEPTHS_CM), figsize=(3.6 * len(PROBE_DEPTHS_CM), 3),
                            sharex=True, sharey=True, squeeze=False)
    for ax, d in zip(axs[0], PROBE_DEPTHS_CM):
        for mode in results:
            st = stats.get((mode, 'sw', d, False))
            if st is not None:
                ax.plot(st.index, st['std'], color=COLORS.get(mode), label=mode)
        H.mark_events(ax, ert_t_h, label_ert=False)
        ax.set_title(f'{d} cm', fontsize=9)
        ax.set_xlabel('time since irrigation start (h)')
    axs[0][0].set_ylabel('ensemble std of saturation (-)')
    axs[0][0].legend(fontsize=8)
    savefig(fig, 'sw_ensemble_spread.png')


safe(spread_figure)

# %% Ensemble mean vs TDR (1:1)
def scatter_figure():
    fig, axs = plt.subplots(1, len(PROBE_DEPTHS_CM), figsize=(3.4 * len(PROBE_DEPTHS_CM), 3.4),
                            sharex=True, sharey=True, squeeze=False)
    for ax, d in zip(axs[0], PROBE_DEPTHS_CM):
        for mode in results:
            n = H.plot_obs_vs_model(ax, stats.get((mode, 'sw', d, False)), smc_sw[d] if d in smc_sw else None,
                                    PHI[d], color=COLORS.get(mode, 'tab:blue'))
            ax.plot([], [], 'o', color=COLORS.get(mode, 'tab:blue'), label=f'{mode} (n={n})')
        lim = ax.get_xlim() + ax.get_ylim()
        ax.plot([min(lim), max(lim)], [min(lim), max(lim)], 'k--', lw=0.8)
        ax.set_title(f'{d} cm', fontsize=9)
        ax.set_xlabel('TDR saturation (-)')
        ax.legend(fontsize=7)
    axs[0][0].set_ylabel('ensemble mean saturation (-)')
    savefig(fig, 'sw_model_vs_TDR.png')


if smc_sw is not None:
    safe(scatter_figure)

# %% Skill against the (independent) soil-moisture sensors
section('SMC skill')
if smc_sw is not None:
    tables = {}
    for mode in results:
        st_d = {d: stats.get((mode, 'sw', d, False)) for d in PROBE_DEPTHS_CM}
        tables[mode] = H.smc_skill(st_d, smc_sw, PHI)
    if any(not t.empty for t in tables.values()):
        skill = pd.concat({m: t for m, t in tables.items() if not t.empty}, names=['run'])
        print(skill.round(3).to_string())
        skill.to_csv(os.path.join(fig_dir, 'smc_skill.csv'))
        print('bias = model - obs (saturation); coverage_5_95 = share of SMC samples inside the ensemble 5-95 % band')
    else:
        print('no overlap between SMC and model time ranges')
else:
    print('skipped (no SMC data)')

# %% Fit to ERT data (RMS)
section('ERT fit')


def rms_figure():
    fig, ax = plt.subplots(2, figsize=[11, 5])
    cplt.DA_RMS(results['withDA']['df_performance'], 'ERT', ax=ax)
    savefig(fig, 'DA_RMS_ERT.png')


if 'withDA' in results:
    safe(rms_figure)


# %% Observed ERT data per survey (apparent resistivity distribution)
def _read_tldat(fname):
    """Data block of a pygimli .tldat/.dat file -> DataFrame (fallback when df_obs lost its objects)."""
    with open(fname) as f:
        lines = f.read().split('\n')
    i0 = 2 + int(lines[0].split()[0])
    nd = int(lines[i0].split()[0])
    cols = lines[i0 + 1].lstrip('#').split()
    rows = [l.split() for l in lines[i0 + 2:i0 + 2 + nd] if l.strip()]
    return pd.DataFrame(rows, columns=cols).apply(pd.to_numeric, errors='coerce')


def ert_hist():
    fig, ax = plt.subplots()
    o = observations.xs('ERT', level=0)
    filt = sorted(glob.glob(os.path.join(path2prj, 'ert_filtered', '*.tldat')))
    for k, (t, th) in enumerate(zip(assimilation_times, ert_t_h)):
        d = (o.xs(t, level=0) if isinstance(o.index, pd.MultiIndex) else o.loc[[t]]).iloc[0]['data']
        if not isinstance(d, pd.DataFrame):          # parquet backup stores it as a string
            d = _read_tldat(filt[k])                 # same order as the surveys (sorted file names)
        col = 'rhoa' if 'rhoa' in d.columns else 'r'
        ax.hist(d[col].dropna(), bins=40, histtype='step', label=f'{th:+.0f} h')
    ax.set_xlabel(r'apparent resistivity ($\Omega$ m)' if col == 'rhoa' else f'{col} (data file column)')
    ax.set_ylabel('count')
    ax.legend(title='survey (h since irrigation)', fontsize=8)
    savefig(fig, 'ERT_data_hist.png')


safe(ert_hist)

# %% Parameter convergence (ic, Ks, Archie)
section('parameters')


def _describe(v, depth=0):
    """One-line summary of an object, to see how dict_parm_pert[p] is organised."""
    if isinstance(v, dict):
        return 'dict{' + ', '.join(f'{k}: {_describe(x, depth + 1)}' for k, x in list(v.items())[:12]) + \
            (', ...' if len(v) > 12 else '') + '}' if depth < 2 else f'dict[{len(v)}]'
    if hasattr(v, 'shape'):
        return f'{type(v).__name__}{tuple(v.shape)}'
    if isinstance(v, (list, tuple)):
        return f'{type(v).__name__}[{len(v)}]'
    return type(v).__name__


def _parm_axis(ax, p, entry):
    """Prior + analysed ensemble of one parameter vs survey time (own plot, no pyCATHY call)."""
    ev = H.parm_evolution(entry)                      # (n_updates + 1, n_members), row 0 = prior
    n_up = ev.shape[0] - 1
    if n_up == len(ert_t_h):
        x = np.r_[ert_t_h[0] - 1.0, ert_t_h]
    else:                                             # unexpected number of updates: plain index
        x = np.arange(ev.shape[0], dtype=float)
    ax.plot(x, ev, color='tab:green', lw=0.4, alpha=0.2)
    m, s = np.nanmean(ev, axis=1), np.nanstd(ev, axis=1, ddof=1)
    ax.fill_between(x, m - s, m + s, color='tab:green', alpha=0.25, lw=0)
    ax.plot(x, m, 'o-', color='tab:green', ms=3, lw=1.5)
    if 'nominal' in entry:
        ax.axhline(float(entry['nominal']), color='k', ls='--', lw=0.8)
    return n_up


def parm_figure(mode):
    """Parameter figures for one run (withDA, withoutDA, ...): not limited to the DA run."""
    dpp = results[mode].get('dict_parm_pert', {})
    parms = [p for p in dpp if isinstance(dpp[p], dict) and 'ini_perturbation' in dpp[p]]
    print(f'[{mode}] {len(parms)} parameters in dict_parm_pert, e.g. {parms[:3]}')
    if not parms:
        return
    ks = [k for k in dpp[parms[0]] if str(k).startswith('update_nb')]
    print(f'keys of {parms[0]}:', list(dpp[parms[0]].keys()))
    print(f'{len(ks)} update_nb* keys, {len(assimilation_times)} assimilation times')
    if not ks:
        print('  no update_nb<k> key: only the prior is plotted (initial ensemble)')

    # one figure per parameter family (ic, SATCELL_VG, thetar_VG, n_VG, ...), not one 45-row figure
    families = {}
    for p in parms:
        families.setdefault(re.sub(r'\d+$', '', p), []).append(p)

    for fam, plist in families.items():
        fig, axs = plt.subplots(len(plist), 2, figsize=(11, 2.2 * len(plist) + 0.5), squeeze=False,
                                gridspec_kw={'width_ratios': [3, 1]})
        for i, p in enumerate(plist):
            try:
                _parm_axis(axs[i, 0], p, dpp[p])
            except Exception as e:
                print(f'{p}: convergence plot failed: {type(e).__name__}: {e}')
            axs[i, 0].set_title(f'{p} [{mode}]', fontsize=9)
            axs[i, 0].set_xlabel('time since irrigation start (h)')
            axs[i, 1].hist(np.ravel(dpp[p]['ini_perturbation']), bins=15)
            axs[i, 1].set_title('initial ensemble', fontsize=9)
        savefig(fig, f'parm_convergence_{fam.rstrip("_")}_{mode}.png')

    tab = H.parm_shift_table(dpp, parms)
    if not tab.empty:
        print(tab.round(3).to_string())
        tab.to_csv(os.path.join(fig_dir, f'parm_shift_{mode}.csv'))


for _mode in results:                              # every loaded run, not only withDA
    safe(lambda m=_mode: parm_figure(m))

# %% Archie scatter (saturation vs converted resistivity), if stored
for mode, res in results.items():
    for key in ('df_Archie', 'Archie'):
        if key in res:
            fig, ax = plt.subplots(figsize=(6, 3))
            safe(lambda: cplt.DA_plot_Archie(res[key], ax=ax))
            savefig(fig, f'archie_{mode}.png')
            break

print(f'\nfigures in {os.path.abspath(fig_dir)}')
plt.show()
