# =============================================================================
# STEP A - DA "solution" run: nb1 model with the true parameters
# =============================================================================
# Runs the forward model, extracts the VWC at the 15/30/50 cm sensor nodes and
# writes noisy synthetic observations to <sol_dir>/solution_obs.pkl.
# Run this once, then run nb4_DA_soil3D_nc.py.
# =============================================================================
import os
import shutil
import pickle
import numpy as np

from pyCATHY import cathy_tools
from nc_model import (setup_model, IC_TRUE, X, Y, ELEVATION,
                      SENSOR_DEPTHS, SENSOR_PHI, HOUR)

# --- config ------------------------------------------------------------------
import project_paths as PP                  # data_path.toml: every folder of the project
path2prj = PP.get().simulations_str          # [outputs].simulations
sol_name = 'SimuID_no_comp_DAsolution'
SWC_NOISE_STD = 0.01            # [m3/m3] sensor noise added to synthetic data
SEED = 42
np.random.seed(SEED)

# --- truth run ---------------------------------------------------------------
sol_dir = os.path.join(path2prj, sol_name)
if os.path.exists(sol_dir):
    shutil.rmtree(sol_dir)

simu_solution = cathy_tools.CATHY(dirName=path2prj, prj_name=sol_name)
setup_model(simu_solution, IC_TRUE)
simu_solution.run_processor(IPRT1=4, DTMAX=1e3, DELTAT=1e2, DTMIN=1e-3, verbose=True)

# --- sensor nodes ------------------------------------------------------------
x0, y0 = X / 2, Y / 2
sensor_nodes = []
for d in SENSOR_DEPTHS:
    n, _ = simu_solution.find_nearest_node([x0, y0, ELEVATION - d])
    sensor_nodes.append(int(n[0]))
print("Sensor nodes:", sensor_nodes)

# --- outputs -> VWC (+ noise) -------------------------------------------------
df_sw, sw_times = simu_solution.read_outputs('sw')
sw_times = np.asarray(sw_times, dtype=float)

theta_true = np.column_stack([df_sw.iloc[:, n].to_numpy() * phi
                              for n, phi in zip(sensor_nodes, SENSOR_PHI)])
theta_obs = theta_true + np.random.normal(0, SWC_NOISE_STD, theta_true.shape)

# --- save --------------------------------------------------------------------
out = {
    'sensor_nodes': sensor_nodes,
    'sensor_depths': SENSOR_DEPTHS,
    'sensor_phi': SENSOR_PHI,
    'sw_times': sw_times,
    'theta_true': theta_true,
    'theta_obs': theta_obs,
    'noise_std': SWC_NOISE_STD,
}
outfile = os.path.join(sol_dir, 'solution_obs.pkl')
with open(outfile, 'wb') as f:
    pickle.dump(out, f)
print("Solution + synthetic observations saved in:", outfile)
