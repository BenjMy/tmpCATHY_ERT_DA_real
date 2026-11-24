# =============================================================================
# Synthetic twin experiment: results, validation against the TRUE sw / psi / parameters, comparison of scenarios
# =============================================================================
# Plotting + scoring only (reads the runs of nb7 and the truth of nb6). Pure functions are in nb8_helpers.py
# (unit-tested); the per-run layout follows nb5_DA_ERT_soil3D_nc_plots.py (median + 5-95 % band + thin members,
# DA vs open loop, hours since irrigation, markers of irrigation / ERT surveys).
#
#   python nb8_plots_synthetic.py --case NC                       # every finished run of NC in ../simulations
#   python nb8_plots_synthetic.py --case HC --runs SYN_HC_ic_Ks_upd_enkf_N100_s42
#   python nb8_plots_synthetic.py --case NC --detail-seed 42 --t-from-h 0
#
# Output (<sim-dir>/figures_nb8/<case>/): per run (seed = --detail-seed)
#   <run>_states.png      truth vs ensemble (DA, open loop) at the sensor nodes, sw and theta (+ virtual TDR)
#   <run>_profiles.png    vertical profiles at the ERT survey times
#   <run>_metrics.png     RMSE / spread / coverage / CRPS in time (all nodes, inside / outside the ERT block)
#   <run>_params.png      parameters vs truth (prior -> final, per update)
# and, across runs:
#   summary.csv, parameters_summary.csv, compare_heatmap_<case>.png, compare_boxplot_<case>.png,
#   compare_params_<case>.png, compare_nens_<case>.png (if several NENS)
# Every section is independent: a failing plot is reported (ONE line) and the others go on; a final table lists what was
# plotted / skipped per run. DA types: enkf / pf (states = FORECAST before each update) and esmda (states = posterior
# forward run, conditioned on ALL the data): their state scores are not like for like, see 'states_kind' in summary.csv.
# =============================================================================
import argparse
import gc
import os
import pickle
import traceback
import warnings

import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

import nb8_helpers as H  # noqa: E402

mpl.rcParams.update({"axes.grid": True, "grid.color": "gray", "grid.linestyle": '--', "grid.alpha": 0.4,
                     "figure.constrained_layout.use": True, "axes.spines.top": False, "axes.spines.right": False})
COL = {'da': 'tab:green', 'ol': 'tab:orange', 'truth': 'k'}
HOUR = 3600.0

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--case', choices=['NC', 'HC'], default='NC')
ap.add_argument('--root', default='../')
ap.add_argument('--sim-dir', default=None)
ap.add_argument('--runs', nargs='+', default=None, help='run names (default: every finished DA run of the case)')
ap.add_argument('--detail-seed', type=int, default=None, help='seed of the runs plotted in detail (default: smallest)')
ap.add_argument('--no-detail', action='store_true')
ap.add_argument('--truth', default=None)
ap.add_argument('--out', default=None)
ap.add_argument('--t-from-h', type=float, default=0.0, help='scores averaged over t >= this [h since irrigation]')
ap.add_argument('--max-nodes', type=int, default=4000, help='node subsample for the scores (HC: 41.6k nodes)')
ap.add_argument('--include-running', action='store_true', help='also score runs not finished (partial results from the per-cycle chunks)')
ap.add_argument('--prefer', default='pkl', choices=['pkl', 'parquet'])
ap.add_argument('--debug', action='store_true', help='print tracebacks (default: one line per problem)')
args = ap.parse_args()

CASE, ROOT = args.case, args.root
SIM_DIR = args.sim_dir or os.path.join(ROOT, 'simulations')
OUT = args.out or os.path.join(SIM_DIR, 'figures_nb8', CASE)
os.makedirs(OUT, exist_ok=True)
TRUTH_FILE = args.truth or os.path.join(ROOT, 'data', '3d', 'synthetic', CASE, f'truth_{CASE}.pkl')


def section(t):
    print(f'\n--- {t} ' + '-' * max(3, 70 - len(t)))


def safe(fn, *a, **k):
    try:
        return fn(*a, **k)
    except Exception as e:
        print(f'  [{fn.__name__}] skipped: {type(e).__name__}: {e}')
        if args.debug:
            traceback.print_exc(limit=3)


STATUS = {}       # run -> {'scores': 'ok' | reason, 'plots': [...], 'skipped': [(what, reason)]}


def note(run, key, value):
    STATUS.setdefault(run, dict(scores='-', plots=[], skipped=[]))
    if key == 'plot':
        STATUS[run]['plots'].append(value)
    elif key == 'skip':
        STATUS[run]['skipped'].append(value)
    else:
        STATUS[run][key] = value


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'    saved {name}', flush=True)


# =============================================================================
# 1. Truth and runs
# =============================================================================
section('truth')
with open(TRUTH_FILE, 'rb') as f:
    truth = pickle.load(f)
assert truth['case'] == CASE
atm_times = np.asarray(truth['atm_times'], float)
STAB = float(truth['stabilize_period'])
nodes_xyz = np.asarray(truth['nodes'], float)
n_nodes = len(nodes_xyz)
ert_t_h = sorted(float(v['t_h_irr']) for v in truth['surveys'].values())
tdr_nodes = [int(n) for n in truth['tdr_nodes']]
tdr_depths_cm = [int(round(d * 100)) for d in truth['tdr_depths']]
print(f"truth {CASE}: {n_nodes} nodes, {len(atm_times)} output times, {len(ert_t_h)} surveys, "
      f"noise {truth['noise']['pct']} %, Archie truth {truth['archie_truth']}")

# nodes used for the scores: every node (NC) or a stride sample + probes + profile column (HC)
ztop = nodes_xyz[:, 2].max()
xc, yc = np.median(nodes_xyz[:, 0]), np.median(nodes_xyz[:, 1])
col_mask = (np.abs(nodes_xyz[:, 0] - nodes_xyz[tdr_nodes[0], 0]) < 1e-6) & (np.abs(nodes_xyz[:, 1] - nodes_xyz[tdr_nodes[0], 1]) < 1e-6)
col_nodes = np.flatnonzero(col_mask)
col_nodes = col_nodes[np.argsort(-nodes_xyz[col_nodes, 2])]            # top first
col_depth = ztop - nodes_xyz[col_nodes, 2]
stride = max(1, int(np.ceil(n_nodes / args.max_nodes)))
score_nodes = np.union1d(np.arange(0, n_nodes, stride), np.r_[tdr_nodes, col_nodes]).astype(int)
print(f"scores on {len(score_nodes)} of {n_nodes} nodes (stride {stride})")

# ERT footprint (electrode block, CATHY frame): inside / outside scores
try:
    import model_config as MC
    from pyCATHY.ERT.alignment import transform_points
    el = MC.load_ert_electrodes(MC.ert_paths(ROOT, CASE)['sequence'])
    el_c = transform_points(el, tuple(truth['ert_frame']['signs']), tuple(truth['ert_frame']['shift']))
    box_lo, box_hi = el_c.min(0) - [0.1, 0.1, 0.1], el_c.max(0) + [0.1, 0.1, 0.1]
except Exception as e:
    warnings.warn(f"electrode block unavailable ({type(e).__name__}: {e}); using a 1.2 m box around the plot centre")
    box_lo, box_hi = np.array([xc - 0.6, yc - 0.6, ztop - 1.2]), np.array([xc + 0.6, yc + 0.6, ztop + 0.1])
in_box = H.nodes_in_box(nodes_xyz[score_nodes], box_lo, box_hi)
REGIONS = {'all': np.ones(len(score_nodes), bool), 'in_ERT_block': in_box, 'outside_ERT_block': ~in_box}
print("ERT block (CATHY frame):", np.round(box_lo, 2), np.round(box_hi, 2), "|", int(in_box.sum()), "nodes inside")

section('runs')
runs = H.collect_runs(SIM_DIR)
assert len(runs), f"no run_meta.json under {SIM_DIR}"
runs = runs[(runs['case'] == CASE)]
if args.runs:
    runs = runs[runs['run'].isin(args.runs) | runs['run'].isin([r + '_noDA' for r in args.runs])]
_ok_status = ['finished', 'running'] if args.include_running else ['finished']
_unfinished = runs[(runs['with_da'] == True) & ~runs['status'].isin(_ok_status)]       # noqa: E712
UNFINISHED = {r['run']: str(r['status']) for _, r in _unfinished.iterrows()}          # reported, never silently dropped
runs = runs[runs['status'].isin(_ok_status)].reset_index(drop=True)
for _run, _st in UNFINISHED.items():
    print(f"  {_run}: NOT plotted, status '{_st}' in run_meta.json (not finished: no results expected)")
assert len(runs), "no finished run for this case: " + (f"{len(UNFINISHED)} run(s) not finished ({', '.join(UNFINISHED)})" if UNFINISHED else "none found")
ol_of = H.pair_open_loop(runs)
da_runs = runs[runs['with_da'] == True].reset_index(drop=True)  # noqa: E712
print(da_runs[['run', 'scenario', 'da', 'nens', 'seed', 'prior_bias']].to_string())


# =============================================================================
# 2. Scores per run
# =============================================================================
def load_df(run):
    from pyCATHY.DA.cathy_DA import DA
    s = DA(dirName=SIM_DIR + '/', prj_name=run, notebook=True)
    pref = ('.pkl', '.parquet') if args.prefer == 'pkl' else ('.parquet', '.pkl')
    res = s.load_backup(prefer=pref, df_DA=False)                 # states: xarray (NetCDF / chunks), no long table
    if 'states' not in res:                                       # old runs: the legacy long table
        res = s.load_backup(prefer=pref)
    assert 'states' in res or 'df_DA' in res, f"no states / df_DA (loaded: {', '.join(res) or 'nothing'})"
    return s, res


def _states_of(res):
    return res['states'] if 'states' in res else res['df_DA']


def scores_of(df, ol):
    """{field: (t_h, {region: DataFrame})} + the arrays needed for the figures."""
    out, arrs = {}, {}
    for field in ('sw', 'psi'):
        if hasattr(df, 'data_vars'):                    # xarray Dataset (compact states)
            assert not ol, "open-loop rows are a separate run with compact states"
            cyc, mem, nod, a = H.ensemble_array_ds(df, field, nodes=score_nodes)
        else:
            cyc, mem, nod, a = H.ensemble_array(df, H.STATE_COL[field], ol=ol, nodes=score_nodes)
        assert np.array_equal(nod, score_nodes), "node ids of df_DA differ from the truth nodes"
        t_s = H.cycles_to_seconds(cyc, atm_times)
        y = H.truth_on_times(truth[field], truth['sw_time'], t_s, nodes=score_nodes)
        out[field] = {r: H.ensemble_scores(a[:, :, m], y[:, m]).set_index((t_s - STAB) / HOUR)
                      for r, m in REGIONS.items()}
        arrs[field] = (t_s, a, y)
    return out, arrs


def states_kind_of(res):
    ds = res.get('states')
    return str(ds.attrs.get('states_kind', 'forecast_before_update')) if hasattr(ds, 'attrs') else 'forecast_before_update'


def one_line(e):
    return f"{type(e).__name__}: {' '.join(str(e).split())[:160]}"



def band(ax, t, a, color, label, members=True):
    """median + 5-95 % band (+ thin members) of a[t, member] ."""
    ax.fill_between(t, np.nanquantile(a, 0.05, axis=1), np.nanquantile(a, 0.95, axis=1), color=color, alpha=0.25, lw=0)
    ax.plot(t, np.nanmedian(a, axis=1), color=color, label=label)
    if members:
        ax.plot(t, a, color=color, lw=0.3, alpha=0.2)


def marks(ax):
    for k, t in enumerate(ert_t_h):
        ax.axvline(t, color='gray', ls=':', lw=0.6, label='ERT survey' if k == 0 else None)
    ax.axvline(0, color='b', ls='--', lw=0.8, label='irrigation start')


def fig_states(name):
    d = RUN[name]
    t_s, a, y = d['arrs']['sw']
    t_h = (t_s - STAB) / HOUR
    t_h_ol = (d['arrs_ol']['sw'][0] - STAB) / HOUR if d['arrs_ol'] else None   # the open loop keeps its own times
    pos = [int(np.flatnonzero(score_nodes == n)[0]) for n in tdr_nodes]
    phi = np.asarray(truth['porosity_nodes'], float)[tdr_nodes]
    tdr = pd.read_csv(truth['tdr_file']) if os.path.isfile(truth.get('tdr_file', '')) else None
    fig, axes = plt.subplots(2, len(tdr_nodes), figsize=(4.2 * len(tdr_nodes), 6.5), sharex=True, squeeze=False)
    for j, (p, dep) in enumerate(zip(pos, tdr_depths_cm)):
        ax = axes[0, j]
        if d['arrs_ol']:
            band(ax, t_h_ol, d['arrs_ol']['sw'][1][:, :, p], COL['ol'], 'open loop', members=False)
        band(ax, t_h, a[:, :, p], COL['da'], 'DA')
        ax.plot(t_h, y[:, p], color=COL['truth'], lw=1.8, label='truth')
        marks(ax)
        ax.set_title(f'{dep} cm'); ax.set_ylabel('sw [-]' if j == 0 else '')
        ax = axes[1, j]                                           # volumetric water content + virtual TDR
        if d['arrs_ol']:
            band(ax, t_h_ol, d['arrs_ol']['sw'][1][:, :, p] * phi[j], COL['ol'], 'open loop', members=False)
        band(ax, t_h, a[:, :, p] * phi[j], COL['da'], 'DA', members=False)
        ax.plot(t_h, y[:, p] * phi[j], color=COL['truth'], lw=1.8, label='truth')
        if tdr is not None and f'theta_{dep}' in tdr:
            ax.plot(tdr['t_h'], tdr[f'theta_{dep}'], '.', ms=2, color='tab:red', alpha=0.5, label='virtual TDR (not assimilated)')
        marks(ax)
        ax.set_xlabel('t since irrigation [h]'); ax.set_ylabel('theta [m3/m3]' if j == 0 else '')
    axes[0, 0].legend(fontsize=7)
    axes[1, 0].legend(fontsize=7)
    fig.suptitle(f"{name}: sw at the sensor nodes vs truth")
    save(fig, f'{name}_states.png')


def fig_profiles(name):
    d = RUN[name]
    t_s, a, y = d['arrs']['sw']
    t_h = (t_s - STAB) / HOUR
    cols = [int(np.flatnonzero(score_nodes == n)[0]) for n in col_nodes]
    sel = sorted({int(np.argmin(np.abs(t_h - t))) for t in np.linspace(min(ert_t_h), max(ert_t_h), 5)})
    fig, axes = plt.subplots(1, len(sel), figsize=(3.0 * len(sel), 5), sharey=True, squeeze=False)
    for ax, i in zip(axes[0], sel):
        i_ol = H.nearest_index(d['arrs_ol']['sw'][0], t_s[i]) if d['arrs_ol'] else None   # same instant, not same position
        if i_ol is not None:
            ol_a = d['arrs_ol']['sw'][1][i_ol]
            ax.fill_betweenx(col_depth, np.nanquantile(ol_a[:, cols], 0.05, axis=0),
                             np.nanquantile(ol_a[:, cols], 0.95, axis=0), color=COL['ol'], alpha=0.2)
        ax.fill_betweenx(col_depth, np.nanquantile(a[i][:, cols], 0.05, axis=0), np.nanquantile(a[i][:, cols], 0.95, axis=0),
                         color=COL['da'], alpha=0.3)
        ax.plot(np.nanmedian(a[i][:, cols], axis=0), col_depth, color=COL['da'], label='DA')
        ax.plot(y[i, cols], col_depth, color=COL['truth'], lw=1.8, label='truth')
        ax.set_title(f't = {t_h[i]:.1f} h'); ax.set_xlabel('sw [-]')
    axes[0, 0].set_ylabel('depth [m]'); axes[0, 0].invert_yaxis(); axes[0, 0].legend(fontsize=7)
    fig.suptitle(f"{name}: sw profile under the plot centre")
    save(fig, f'{name}_profiles.png')


def fig_metrics(name):
    d = RUN[name]
    fig, axes = plt.subplots(2, 4, figsize=(17, 7), squeeze=False)
    for r, field in enumerate(('sw', 'psi')):
        for c, (metric, lab) in enumerate((('rmse', 'RMSE of the mean'), ('spread_skill', 'spread / RMSE (1 = reliable)'),
                                           ('coverage', '90 % coverage (0.9 = reliable)'), ('crps', 'CRPS'))):
            ax = axes[r, c]
            for reg, ls in (('all', '-'), ('in_ERT_block', '--'), ('outside_ERT_block', ':')):
                s = d['sc'][field][reg]
                ax.plot(s.index, s[metric], ls=ls, color=COL['da'], label=f'DA, {reg}')
                if d['sc_ol']:
                    so = d['sc_ol'][field][reg]
                    ax.plot(so.index, so[metric], ls=ls, color=COL['ol'], label=f'open loop, {reg}')
            if metric == 'spread_skill':
                ax.axhline(1, color='k', lw=0.6)
            if metric == 'coverage':
                ax.axhline(0.9, color='k', lw=0.6)
            marks(ax)
            ax.set_title(f'{field}: {lab}', fontsize=9)
            ax.set_xlabel('t since irrigation [h]')
    axes[0, 0].legend(fontsize=6)
    fig.suptitle(f"{name}: validation against the true state")
    save(fig, f'{name}_metrics.png')


def fig_params(name):
    d = RUN[name]
    pt = d.get('ptab')
    if pt is None or not len(pt):
        return 'n/a: no updated parameter (state-only run)'      # expected, not a failure
    dpp, pds = d['res'].get('dict_parm_pert', {}), d['res'].get('params')
    n = len(pt)
    nc = min(4, n)
    nr = int(np.ceil(n / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(3.6 * nc, 2.8 * nr), squeeze=False)
    for ax, (p, r) in zip(axes.ravel(), pt.iterrows()):
        ev = H.parameter_evolution(p, dpp, pds)
        x = np.arange(ev.shape[0])
        ax.fill_between(x, np.nanquantile(ev, 0.05, axis=1), np.nanquantile(ev, 0.95, axis=1), color=COL['da'], alpha=0.3)
        ax.plot(x, np.nanmean(ev, axis=1), color=COL['da'])
        if r['truth'] is not None and np.isfinite(r['truth']):
            ax.axhline(r['truth'], color='k', lw=1.5)
        ax.set_title(f"{p}  z {r['z_prior']:.1f} -> {r['z_final']:.1f}", fontsize=8)
        ax.set_xlabel('iteration' if d['meta']['da'] == 'esmda' else 'update')
    for ax in axes.ravel()[n:]:
        ax.axis('off')
    fig.suptitle(f"{name}: parameters (mean, 5-95 %) vs truth (black); z = (mean - truth) / sd")
    save(fig, f'{name}_params.png')


# =============================================================================
# 2-4. PER RUN: load -> score -> summary row -> figures -> free memory (one run in memory at a time)
# =============================================================================
section('per run: load, score, plot, free')
RUN = {}          # run -> light results only (scores, parameter table); the big arrays are freed after each run
rows, prow = [], []
scen_cache = {}
seeds = sorted(da_runs['seed'].unique())
dseed = args.detail_seed if args.detail_seed is not None else seeds[0]
for _, r in da_runs.iterrows():
    name = r['run']
    try:
        sim, res = load_df(name)
        sc, arrs = scores_of(_states_of(res), ol=False)
        ol_name = ol_of.get(name)
        sc_ol = arrs_ol = None
        if ol_name:                                     # open-loop run (--no-da) = separate project
            try:
                _, res_ol = load_df(ol_name)
                sc_ol, arrs_ol = scores_of(_states_of(res_ol), ol=False)
                del res_ol
            except Exception as e:                      # a broken open loop must not lose the DA run
                note(name, 'skip', ('open loop', one_line(e)))
        elif 'df_DA' in res and (res['df_DA']['OL'].astype(bool)).any():
            sc_ol, arrs_ol = scores_of(res['df_DA'], ol=True)
        RUN[name] = dict(meta=r, sim=sim, res=res, sc=sc, arrs=arrs, sc_ol=sc_ol, arrs_ol=arrs_ol,
                         states_kind=states_kind_of(res))
        note(name, 'scores', 'ok' + ('' if sc_ol else ' (no open loop)'))
        print(f"  {name}: scores ok [{RUN[name]['states_kind']}]" + ("" if sc_ol else " (no open loop)"), flush=True)
    except Exception as e:
        note(name, 'scores', 'SKIPPED: ' + one_line(e))
        print(f"  {name}: skipped: {one_line(e)}", flush=True)
        if args.debug:
            traceback.print_exc(limit=3)
        continue
    d = RUN[name]
    m = d['meta']
    row = dict(run=name, scenario=m['scenario'], da=m['da'], nens=m['nens'], seed=m['seed'], n_par=m['n_par'],
               states_kind=d['states_kind'])
    for field in ('sw', 'psi'):
        for reg in REGIONS:
            t_h = d['sc'][field][reg].index.values
            s = H.summarize_run(d['sc'][field][reg], d['sc_ol'][field][reg] if d['sc_ol'] else None,
                                args.t_from_h, t_h)
            row.update({f'{field}_{reg}_{k}': v for k, v in s.items()})
    rows.append(row)
    # parameters vs truth
    try:
        import tomllib
        import make_scenarios as MS
        key = (m['scenario'], m['prior_bias'])
        if key not in scen_cache:
            scen_cache[key] = tomllib.loads(MS.build(biased=(m['prior_bias'] == 'biased')))['scenario'][m['scenario']]
        scen = scen_cache[key]
        dpp = d['res'].get('dict_parm_pert', {})
        upd = [p for p in scen['listUpdateParm'] if p != 'St. var.']
        pds = d['res'].get('params')                       # ES-MDA: every iterate (<project>_params.nc)
        pt = H.parameter_table(dpp, upd, scen, truth['soil_table'], truth['archie_truth'], params_ds=pds)
        d['ptab'], d['scen'] = pt, scen
        if len(pt):
            prow.append(dict(run=name, scenario=m['scenario'], da=m['da'], seed=m['seed'], nens=m['nens'],
                             n_par=len(pt), mean_abs_z_prior=pt['z_prior'].abs().mean(),
                             mean_abs_z_final=pt['z_final'].abs().mean(), mean_contraction=pt['contraction'].mean(),
                             frac_truth_in_90=pt['truth_in_90'].mean(),
                             mean_rel_err_prior=(pt['err_prior'].abs() / pt['truth'].abs()).mean(),
                             mean_rel_err_final=(pt['err_final'].abs() / pt['truth'].abs()).mean()))
            pt.to_csv(os.path.join(OUT, f'{name}_parameters.csv'))
    except Exception as e:
        note(name, 'skip', ('parameter table', one_line(e)))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, 'summary.csv'), index=False)      # written after every run: a crash keeps it
    # ---- per-run figures, now, while the arrays are in memory
    if not args.no_detail and r['seed'] == dseed:
        print('  plotting', name, flush=True)
        for fn in (fig_states, fig_profiles, fig_metrics, fig_params):
            try:
                why = fn(name)
                if why:
                    note(name, 'skip', (fn.__name__.replace('fig_', ''), why))
                else:
                    note(name, 'plot', fn.__name__.replace('fig_', ''))
            except Exception as e:
                note(name, 'skip', (fn.__name__.replace('fig_', ''), one_line(e)))
                print(f'  [{fn.__name__}] skipped: {one_line(e)}')
                if args.debug:
                    traceback.print_exc(limit=3)
    # ---- free the big objects: only scores (small) and the parameter table stay
    for k in ('arrs', 'arrs_ol', 'res', 'sim'):
        d[k] = None
    del sim, res, sc, arrs, sc_ol, arrs_ol
    gc.collect()

assert RUN, "no run could be scored (see the reasons above)"
section('summary')
summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT, 'summary.csv'), index=False)
psum = pd.DataFrame(prow)
if len(psum):
    psum.to_csv(os.path.join(OUT, 'parameters_summary.csv'), index=False)
_cols = ['run', 'sw_all_rmse', 'sw_all_rmse_ol', 'sw_all_skill_vs_ol', 'sw_in_ERT_block_rmse',
         'sw_outside_ERT_block_rmse', 'sw_all_coverage', 'sw_all_crps']
print(summary[[c for c in _cols if c in summary]].round(4).to_string())

# =============================================================================
# 5. Comparison of scenarios and DA types
# =============================================================================
def agg(col):
    return summary.groupby(['scenario', 'da'])[col].mean().unstack('da')


def fig_heatmap():
    items = [('sw_all_rmse', 'RMSE sw (all nodes)', 'viridis_r'), ('sw_in_ERT_block_rmse', 'RMSE sw inside the ERT block', 'viridis_r'),
             ('sw_outside_ERT_block_rmse', 'RMSE sw outside the block', 'viridis_r'),
             ('sw_all_skill_vs_ol', 'skill vs open loop (1 - RMSE/RMSE_ol)', 'RdYlGn'),
             ('sw_all_crps', 'CRPS sw', 'viridis_r'), ('sw_all_coverage', '90 % coverage sw (target 0.9)', 'coolwarm'),
             ('psi_all_rmse', 'RMSE psi [m]', 'viridis_r')]
    items = [i for i in items if i[0] in summary and summary[i[0]].notna().any()]
    nc = min(4, len(items))
    nr = int(np.ceil(len(items) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(4.3 * nc, 0.6 * summary['scenario'].nunique() * nr + 2.2 * nr), squeeze=False)
    for ax, (col, title, cmap) in zip(axes.ravel(), items):
        t = agg(col)
        im = ax.imshow(t.values, cmap=cmap, aspect='auto', **({'vmin': 0.5, 'vmax': 1.0} if 'coverage' in col else {}))
        ax.set_xticks(range(t.shape[1])); ax.set_xticklabels(t.columns)
        ax.set_yticks(range(t.shape[0])); ax.set_yticklabels([s[len(CASE) + 1:] for s in t.index], fontsize=8)
        for i in range(t.shape[0]):
            for j in range(t.shape[1]):
                if np.isfinite(t.values[i, j]):
                    ax.text(j, i, f'{t.values[i, j]:.3g}', ha='center', va='center', fontsize=7)
        ax.set_title(title, fontsize=9); ax.grid(False)
        fig.colorbar(im, ax=ax, shrink=0.8)
    for ax in axes.ravel()[len(items):]:
        ax.axis('off')
    mixed = summary['states_kind'].nunique() > 1
    fig.suptitle(f'{CASE}: scenarios x DA type (mean over seeds, t >= {args.t_from_h:g} h)'
                 + ('\nesmda = smoother (states conditioned on ALL the data); enkf/pf = forecast before each update:'
                    ' state scores not like for like' if mixed else ''))
    save(fig, f'compare_heatmap_{CASE}.png')


def fig_boxplot():
    if summary['seed'].nunique() < 2:
        print('  one seed only: no box plot')
        return
    groups = [(sc, da) for sc in sorted(summary['scenario'].unique()) for da in sorted(summary['da'].unique())]
    data = [summary[(summary['scenario'] == sc) & (summary['da'] == da)]['sw_all_rmse'].dropna().values for sc, da in groups]
    ol = summary.groupby('scenario')['sw_all_rmse_ol'].mean() if 'sw_all_rmse_ol' in summary else None
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(groups)), 4.5))
    ax.boxplot(data)
    ax.set_xticklabels([f"{sc[len(CASE) + 1:]}\n{da}" for sc, da in groups])
    if ol is not None:
        ax.axhline(ol.mean(), color=COL['ol'], ls='--', label='open loop (mean)')
        ax.legend()
    ax.set_ylabel('RMSE sw (all nodes)'); plt.setp(ax.get_xticklabels(), rotation=60, ha='right', fontsize=7)
    ax.set_title(f'{CASE}: spread over seeds')
    save(fig, f'compare_boxplot_{CASE}.png')


def fig_params_compare():
    if not len(psum):
        print('  no parameter summary')
        return
    t = psum.groupby(['scenario', 'da'])[['mean_abs_z_prior', 'mean_abs_z_final', 'mean_contraction', 'frac_truth_in_90',
                                          'mean_rel_err_prior', 'mean_rel_err_final']].mean().reset_index()
    t['label'] = [f"{s[len(CASE) + 1:]}\n{d}" for s, d in zip(t['scenario'], t['da'])]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
    x = np.arange(len(t))
    axes[0].bar(x - 0.2, t['mean_rel_err_prior'], 0.4, color='lightgray', label='prior')
    axes[0].bar(x + 0.2, t['mean_rel_err_final'], 0.4, color=COL['da'], label='final')
    axes[0].set_title('mean |relative error| of the parameters'); axes[0].legend()
    axes[1].bar(x - 0.2, t['mean_abs_z_prior'], 0.4, color='lightgray'); axes[1].bar(x + 0.2, t['mean_abs_z_final'], 0.4, color=COL['da'])
    axes[1].set_title('mean |z| = |mean - truth| / sd (final > prior = overconfident)')
    axes[2].bar(x, t['mean_contraction'], color=COL['da']); axes[2].set_title('contraction 1 - sd_post / sd_prior')
    axes[3].bar(x, t['frac_truth_in_90'], color=COL['da']); axes[3].axhline(0.9, color='k', lw=0.6)
    axes[3].set_title('share of parameters with the truth in the final 90 % range')
    for ax in axes:
        ax.set_xticks(x); ax.set_xticklabels(t['label'], rotation=60, ha='right', fontsize=7)
    fig.suptitle(f'{CASE}: parameters vs truth')
    save(fig, f'compare_params_{CASE}.png')


def fig_nens():
    if summary['nens'].nunique() < 2:
        print('  one ensemble size only')
        return
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for (sc, da), g in summary.groupby(['scenario', 'da']):
        t = g.groupby('nens')['sw_all_rmse'].agg(['mean', 'std'])
        ax.errorbar(t.index, t['mean'], t['std'].fillna(0), marker='o', label=f'{sc[len(CASE) + 1:]} {da}')
    ax.set_xlabel('ensemble size'); ax.set_ylabel('RMSE sw'); ax.legend(fontsize=7)
    save(fig, f'compare_nens_{CASE}.png')


section('comparison figures')
for fn in (fig_heatmap, fig_boxplot, fig_params_compare, fig_nens):
    safe(fn)


# =============================================================================
# 6. What was plotted / skipped
# =============================================================================
section('status per run')
skipped_any = False
for run, why in UNFINISHED.items():
    STATUS[run] = dict(scores=f"NOT FINISHED (status '{why}')", plots=[], skipped=[])
for run in list(da_runs['run']) + list(UNFINISHED):
    st = STATUS.get(run, dict(scores='not loaded', plots=[], skipped=[]))
    skipped_any |= not str(st['scores']).startswith('ok') or any(not str(w).startswith('n/a') for _, w in st['skipped'])
    print(f"  {run}\n      scores: {st['scores']} | plotted: {', '.join(st['plots']) or '-'}")
    for what, why in st['skipped']:
        print(f"      skipped {what}: {why}")
print(f"\nfigures and tables in {OUT}")
raise SystemExit(1 if skipped_any else 0)
