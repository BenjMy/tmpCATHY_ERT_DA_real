# =============================================================================
# Agripolis 3D hydro (NC) + Data Assimilation of soil moisture (twin experiment)
# =============================================================================
# Same structure as DA_ERT_part2_runDA.py, with the nb1 model (nc_model.py).
#   STEP A  (nb4a_DA_solution_soil3D_nc.py, run first)
#           truth run -> synthetic VWC at 15/30/50 cm (+ noise)
#   STEP B  this script: DA / open-loop project with a biased prior ensemble
#   STEP C  this script: assimilate the synthetic VWC with EnKF
#
# Set WITH_DA = False for the open-loop run (obs error = 1e99, as in the ERT
# tutorial), True for the assimilation run.
#
# NOTE: the SWC observation block follows the pyCATHY DA API as used in the ERT
# example. Check the `meta` keys / `data_type` against your pyCATHY version
# (see pyCATHY/DA/observations.py).
# =============================================================================
import os
import shutil
import pickle
import numpy as np

from pyCATHY import cathy_utils
from pyCATHY.DA.cathy_DA import DA, dictObs_2pd
from pyCATHY.DA.observations import read_observations, make_data_cov
from pyCATHY.DA import perturbate

from nc_model import IRR_START_HOUR, HOUR
from nc_model import setup_model, snap_to_times, X, Y, ELEVATION

# =============================================================================
# 0. Config
# =============================================================================
WITH_DA = True

# Parameters perturbed AND updated by the filter. Names follow the pyCATHY
# scenario convention ('ic', 'Ks' ...). To add van Genuchten parameters append
# their scenario names here (e.g. 'VGN'); check the names accepted by your
# pyCATHY version with cathy_utils.create_scenario_file(..., use_common_values=False)
# and inspect the generated json.
# PARAM_NAMES = ['ic', 'Ks']
PARAM_NAMES = ['ic']
SCENARIO = 'nc_' + '_'.join(PARAM_NAMES)     # scenario name / file name

# Prior (deliberately biased, so that DA has something to correct)
IC_PRIOR_MEAN, IC_PRIOR_STD = -1.5, 0.5      # [m] pressure head
KS_PRIOR_FACTOR, KS_PRIOR_STD = 2.0, 0.5     # multiplicative on nominal Ks
KS_NOMINAL = 2.0e-06                         # [m/s] ~ nb1 layer-1 Ks

NENS = 32
# Extra times [s] that must exist in the atmbc vector (e.g. ERT survey times).
# Must be identical to the list used in the truth run (STEP A).
# e.g. EXTRA_TIMES = [2 * HOUR, 5.5 * HOUR]
EXTRA_TIMES = []
OBS_EVERY_HOURS = 3             # assimilation frequency


import project_paths as PP                  # data_path.toml: every folder of the project
path2prj = PP.get().simulations_str          # [outputs].simulations
NB1_PRJ = 'SimuID_no_comp'      # project folder written by forward_model_nc.py
sol_name = 'SimuID_no_comp_DAsolution'
da_name = ('NC_withDA_' if WITH_DA else 'NC_withoutDA_') + SCENARIO

# =============================================================================
# 2. Load truth + synthetic observations (produced by STEP A)
# =============================================================================
with open(os.path.join(path2prj, sol_name, 'solution_obs.pkl'), 'rb') as f:
    sol = pickle.load(f)

sensor_nodes = sol['sensor_nodes']
sw_times = sol['sw_times']
theta_true = sol['theta_true']
theta_obs = sol['theta_obs']
SWC_NOISE_STD = sol['noise_std']
OBS_ERR = SWC_NOISE_STD if WITH_DA else 1e99

# =============================================================================
# 3. STEP B - DA project
# =============================================================================
da_dir = os.path.join(path2prj, da_name)
if os.path.exists(da_dir):
    shutil.rmtree(da_dir)

# DA project = copy of the nb1 project folder (mesh, soil, atmbc, bc already set)
nb1_dir = os.path.join(path2prj, NB1_PRJ)
assert os.path.isdir(nb1_dir), f"nb1 project not found: {nb1_dir} (run nb1 first)"
shutil.copytree(nb1_dir, da_dir)

simuDA = DA(dirName=path2prj, prj_name=da_name)
setup_model(simuDA, [[IC_PRIOR_MEAN]] * 3, extra_times=EXTRA_TIMES)

# =============================================================================
# 4. Observations (assimilate VWC at the 3 sensors)
# =============================================================================
# observation times: when irrigation starts and then every OBS_EVERY_HOURS
atm_times = np.unique(np.asarray(simuDA.read_inputs('atmbc').time, dtype=float))
assim_times = [t for t in atm_times
               if t >= IRR_START_HOUR * HOUR and (t / HOUR) % OBS_EVERY_HOURS == 0]

sensor_nodes = np.atleast_1d(np.asarray(sensor_nodes)).ravel()

# SWC sensors are read one by one (sensor loop > time loop), each with its own
# mesh node; read_observations returns the updated dict (see pyCATHY doc
# plot_2_prepare_SMC_obs4DA). Sensors are auto-named swc, swc2, swc3 ...
dict_obs = {}
for i_s, node in enumerate(sensor_nodes):
    for tt in assim_times:
        i_t = int(np.argmin(np.abs(sw_times - tt)))
        dict_obs = read_observations(
            dict_obs,
            obs_2_add=float(theta_obs[i_t, i_s]),
            tA=float(tt),
            mesh_nodes=int(node),
            data_type='swc',
            data_err=OBS_ERR,
            colname=' m³/m³ Water Content',
            datetime=None,
        )

data_measure_df = dictObs_2pd(dict_obs)
print(data_measure_df.head())
# df_atmbc = simuDA.read_inputs('atmbc')
# df_atmbc.time.unique()

data_cov, data_pert, stacked_data_cov = make_data_cov(
    simuDA, dict_obs, list_assimilated_obs='all',
    nb_assimilation_times=len(dict_obs))
simuDA.stacked_data_cov = stacked_data_cov

# =============================================================================
# 5. Parameter perturbation scenario (cathy_utils.create_scenario_file)
# =============================================================================
# 5.1 create the scenario file with suggested values for each parameter
scenario_file = os.path.join(simuDA.workdir, simuDA.project_name,
                             SCENARIO + '.json')
cathy_utils.create_scenario_file(
    scenario_name=os.path.join(simuDA.workdir, simuDA.project_name, SCENARIO),
    param_names=PARAM_NAMES,
    use_common_values=True,      # suggested values instead of <placeholders>
    filetype='json',
)

# 5.2 read it back and override the suggested values with our prior
scenario = cathy_utils.read_scenario_file(scenario_file)

# top-level key is the full path passed to create_scenario_file, not SCENARIO
key = next(k for k in scenario if not k.startswith('_'))
sc = scenario[key]

prior = {
    'ic': dict(per_nom=IC_PRIOR_MEAN, per_mean=IC_PRIOR_MEAN,
               per_sigma=IC_PRIOR_STD, per_bounds="None"),
    # 'Ks': dict(per_nom=KS_NOMINAL, per_mean=KS_PRIOR_FACTOR,
    #            per_sigma=KS_PRIOR_STD, per_bounds="None"),
}
for i, name in enumerate(sc['per_name']):
    for field, val in prior[name].items():
        sc[field][i] = val

# make sure the filter updates the state AND every perturbed parameter
sc['listUpdateParm'] = ['St. var.'] + PARAM_NAMES

# 5.3 inspect
print(cathy_utils.scenario_dict_to_df_list(scenario, key))

# 5.4 generate the ensembles

# simuDA.read_inputs('root_map')
# simuDA.update_veg_map()

list_pert = perturbate.perturbate(simuDA, sc, NENS)

var_per_dict_stacked = {}
for dp in list_pert:
    var_per_dict_stacked = perturbate.perturbate_parm(
        var_per_dict_stacked,
        parm=dp,
        type_parm=dp['type_parm'],
        mean=dp['mean'],
        sd=dp['sd'],
        sampling_type=dp['sampling_type'],
        ensemble_size=dp['ensemble_size'],
        per_type=dp['per_type'],
        savefig=os.path.join(simuDA.workdir, simuDA.project_name,
                             simuDA.project_name + '_' + str(dp['savefig'])),
    )

# parameters updated by the filter: state variables + all perturbed parameters
list_parm2update = ['St. var.'] + PARAM_NAMES

simuDA.update_parm(IPRT1=2, TRAFLAG=0, 
                   DTMIN=1e-3, DTMAX=1e3, DELTAT=1e2, VTKF=2)

#%%
# =============================================================================
# 6. Run EnKF
# =============================================================================
# simuDA.run_DA_sequential(
#     parallel=True,
#     dict_obs=dict_obs,
#     list_assimilated_obs='all',
#     list_parm2update=list_parm2update,
#     DA_type='enkf_Evensen2009_Sakov',
#     dict_parm_pert=var_per_dict_stacked,
#     # open_loop_run=not WITH_DA,
#     threshold_rejected=80,
# )
#%%
# =============================================================================
# 7. Results
# =============================================================================
# Reload with simuDA.load_pickle_backup() and plot with
# cplt.DA_plot_time_dynamic(results_DA['df_DA'], 'sw', node, ...) as in
# DA_ERT_part3_results.py, overlaying theta_true / phi from
# <sol_dir>/solution_obs.pkl.
print("DA run finished:", da_name)