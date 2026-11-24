# =============================================================================
# Plot results of nb4_DA_soil3D_nc.py (SWC assimilation, twin experiment)
# =============================================================================
# Plotting only: reloads solution_obs.pkl (truth + synthetic obs, STEP A) and the
# DA backups written by nb4 (with / without DA). No simulation is run.
# Set PARAM_NAMES / SCENARIO / path2prj exactly as in nb4.
# =============================================================================
import os
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyCATHY.DA.cathy_DA import DA, dictObs_2pd
import pyCATHY.plotters.cathy_plots as cplt

mpl.rcParams.update({"axes.grid": True, "grid.color": "gray",
                     "grid.linestyle": '--'})

# %% Settings (must match nb4)
path2prj = "../simulations/"
sol_name = 'SimuID_no_comp_DAsolution'
PARAM_NAMES = ['ic']
SCENARIO = 'nc_' + '_'.join(PARAM_NAMES)
PRJ = {'withDA': 'NC_withDA_' + SCENARIO,
       'withoutDA': 'NC_withoutDA_' + SCENARIO}
POROS = None            # None = read from the soil file (must be uniform);
                        # else float, or one value per sensor [phi15, phi30, phi50]
DEPTH_LABELS = ['15 cm', '30 cm', '50 cm']
HOUR = 3600.0
fig_dir = os.path.join(path2prj, 'figures_nb4')
os.makedirs(fig_dir, exist_ok=True)

# %% Truth + synthetic observations
with open(os.path.join(path2prj, sol_name, 'solution_obs.pkl'), 'rb') as f:
    sol = pickle.load(f)
sensor_nodes = np.atleast_1d(np.asarray(sol['sensor_nodes'])).ravel().astype(int)
sw_times = np.asarray(sol['sw_times'], dtype=float)
theta_true = np.asarray(sol['theta_true'])
theta_obs = np.asarray(sol['theta_obs'])
ns = len(sensor_nodes)
labels = (DEPTH_LABELS + [f'sensor {i}' for i in range(ns)])[:ns]


# %% Reload DA backups
def load_backup(simu):
    for fn in ('load_pickle_backup', 'load_parquet_backup'):
        try:
            return getattr(simu, fn)()
        except Exception as e:
            print(f'  {fn} failed for {simu.project_name}: {type(e).__name__}')
    raise RuntimeError(f'no backup readable for {simu.project_name}')


simus, results = {}, {}
for mode, prj in PRJ.items():
    if not os.path.isdir(os.path.join(path2prj, prj)):
        print(f'skipping {mode}: {prj} not found')
        continue
    simus[mode] = DA(dirName=path2prj, prj_name=prj, notebook=True)
    results[mode] = load_backup(simus[mode])
assert results, 'no DA results found'
ref_mode = 'withDA' if 'withDA' in results else next(iter(results))

if POROS is None:
    soil_SPP, _ = simus[ref_mode].read_inputs('soil', MAXVEG=1)
    phi_u = soil_SPP['POROS'].unique()
    if len(phi_u) > 1:
        print('WARNING: several POROS values in soil file, using the first; '
              'set POROS explicitly (one per sensor).', phi_u)
    POROS = float(phi_u[0])
POROS = np.broadcast_to(np.asarray(POROS, dtype=float), (ns,))

# %% Observations / assimilation times
observations = results[ref_mode]['df_obs']
assimilation_times = observations.index.get_level_values(1).unique().to_list()
obs_names = observations.index.get_level_values(0).unique().to_list()

for res in results.values():
    df = res['df_DA']
    ut = df['time'].unique()
    n = min(len(ut), len(assimilation_times))
    df['assimilation_times'] = df['time'].map(dict(zip(ut[:n],
                                                       assimilation_times[:n])))


def savefig(fig, name):
    fig.savefig(os.path.join(fig_dir, name), dpi=300, bbox_inches='tight')


# %% Observations vs truth
fig, axs = plt.subplots(1, ns, figsize=(3.2 * ns, 3), sharey=True, squeeze=False)
for j, ax in enumerate(axs[0]):
    ax.plot(sw_times / HOUR, theta_true[:, j], 'r-', label='truth')
    ax.plot(sw_times / HOUR, theta_obs[:, j], 'k.', ms=3, alpha=.5, label='obs')
    ax.set_title(labels[j]); ax.set_xlabel('time (h)')
axs[0][0].set_ylabel('VWC (m³/m³)'); axs[0][0].legend()
savefig(fig, 'obs_vs_truth.png')


# %% Ensemble states at the sensor nodes vs truth
def plot_states(mode, var):
    fig, axs = plt.subplots(1, ns, figsize=(3.2 * ns, 3), sharey=True,
                            squeeze=False)
    for j, (ax, node) in enumerate(zip(axs[0], sensor_nodes)):
        cplt.DA_plot_time_dynamic(results[mode]['df_DA'], var, node, ax=ax,
                                  keytime='assimilation_times')
        if var == 'sw':
            ax.plot(sw_times, theta_true[:, j] / POROS[j], 'r.-', label='truth')
            sub = observations.xs(obs_names[j], level=0) \
                if j < len(obs_names) else None
            if sub is not None:
                ax.plot(sub.index.values, sub['data'].values / POROS[j], 'k*',
                        label='obs')
        ax.set_title(f'{labels[j]} (node {node})')
        if ax.get_legend():
            ax.get_legend().set_visible(False)
    fig.suptitle(f'{mode} - {var}')
    savefig(fig, f'{var}_{mode}.png')


for mode in results:
    plot_states(mode, 'sw')
    plot_states(mode, 'psi')


# %% Ensemble-mean error vs truth, with / without DA
def mean_sw(mode, node):
    df = results[mode]['df_DA']
    d = df[df['variable'] == 'sw'] if 'variable' in df.columns else None
    return d


fig, axs = plt.subplots(1, ns, figsize=(3.2 * ns, 3), sharey=True, squeeze=False)
for j, (ax, node) in enumerate(zip(axs[0], sensor_nodes)):
    for mode, res in results.items():
        df = res['df_DA']
        try:
            d = df[df['variable'] == 'sw'] if 'variable' in df else df
            col = d[node] if node in d.columns else None
            if col is None:
                continue
            m = d.groupby('assimilation_times')[node].mean()
            truth_i = np.interp(m.index.values, sw_times,
                                theta_true[:, j] / POROS[j])
            ax.plot(m.index.values / HOUR, m.values - truth_i, label=mode)
        except Exception as e:
            print('bias plot skipped:', type(e).__name__, e)
    ax.axhline(0, color='r'); ax.set_title(labels[j]); ax.set_xlabel('time (h)')
axs[0][0].set_ylabel('mean(sw) - truth (-)'); axs[0][0].legend()
savefig(fig, 'sw_bias.png')

# %% DA performance (RMSE)
if 'withDA' in results:
    fig, ax = plt.subplots(2, figsize=[11, 5])
    for name in obs_names:
        cplt.DA_RMS(results['withDA']['df_performance'], name, ax=ax)
    savefig(fig, 'DA_RMS.png')

# %% Parameter convergence
for mode in results:
    dpp = results[mode].get('dict_parm_pert', {})
    for parm in dpp:
        if parm.endswith('_ini') or 'ini_perturbation' not in dpp[parm]:
            continue
        fig, axs = plt.subplots(2, 1, figsize=[11, 5])
        cplt.DA_plot_parm_dynamic_scatter(
            parm=parm, dict_parm_pert=dpp,
            list_assimilation_times=assimilation_times, ax=axs[0])
        axs[0].set_title(f'{mode} - {parm}')
        axs[1].hist(dpp[parm]['ini_perturbation'], bins=15)
        axs[1].set_xlabel(f'{parm} initial ensemble')
        savefig(fig, f'parm_{parm}_{mode}.png')

plt.show()
