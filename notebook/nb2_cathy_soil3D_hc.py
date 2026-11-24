# =============================================================================
# Agripolis 3D hydro: COMPACTED (HC)
# =============================================================================
# This script shows how to use pyCATHY to run the hydrological model on a small scale soil scenario.

import sys
import os
import json
import shutil
import pickle
from IPython import embed
import pyvista as pv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pyCATHY.meshtools as mt
from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt
from pathlib import Path

# =============================================================================
# 0. Setup paths and experiment parameters
# =============================================================================
compLevel = 'heavy_comp'
path2prj = "../simulations/"
prj_name = 'SimuID_' + compLevel

prj_path = os.path.join(path2prj, prj_name)

if os.path.exists(prj_path):
    shutil.rmtree(prj_path)

os.makedirs(prj_path)

# Local CATHY example containing the required source files
# cathy_source = "/home/alberto/pycathy_wrapper/examples/SSHydro/soil3d_from_weill"

# # Copy the CATHY project structure
# for item in os.listdir(cathy_source):
#     src = os.path.join(cathy_source, item)
#     dst = os.path.join(prj_path, item)

#     if os.path.isdir(src):
#         shutil.copytree(src, dst)
#     else:
#         shutil.copy2(src, dst)

# initialize CATHY
simuInf3D = cathy_tools.CATHY(
    dirName=path2prj,
    prj_name=prj_name
)
# ============================================================
# 1. Mesh 3D
# ============================================================
# The mesh is a 3D grid that divides the soil volume into small cells:
# 2x2 m area with a compacted central band: maintains a
# lateral buffer that allows the flow to spread horizontally if necessary.

### 1.1 horizontal discretization
delta_x = 0.1  # [m] resolution x
delta_y = 0.1  # [m] resolution y
x = 5.0  # [m] extension x
y = 5.0  # [m] extension y
maxdepth = 2.0  # [m] depth
elevation = 1e-3  # [m]

ncellx = int(np.round(x / delta_x))
ncelly = int(np.round(y / delta_y))

DEM = np.ones((ncelly, ncellx)) * elevation
DEM[-1, -1] = elevation - 1e-4  # outlet

### 1.2 vertical discretization
zr_thickness = [
    0.03, 0.03, 0.03, 0.03,  # 0-12 cm
    0.04, 0.04, 0.04,        # 12-24 cm
    0.06, 0.06,              # 24-36 cm
    0.08, 0.08,              # 36-52 cm
    0.10,                    # 52-62 cm
    0.25,                    # 62-87 cm
    0.40,                    # 87-127 cm
    0.73                     # 127-200 cm
]
nstr = len(zr_thickness)

zratio = [z / maxdepth for z in zr_thickness]
assert np.isclose(sum(zratio), 1.0), "zratio does not sum to 1"

### 1.3 update input preprocessor
simuInf3D.update_prepo_inputs(
    DEM, delta_x=delta_x, delta_y=delta_y, dr=delta_x,
    zratio=zratio, base=maxdepth)

cell_size = delta_x * delta_y  # area
# simuInf3D.show_input('dem')

# =============================================================================
# 2. Vegetation (minimum homogeneous grass)
# =============================================================================
root_map = np.ones(DEM.shape)
simuInf3D.update_veg_map(root_map)

simuInf3D.create_mesh_vtk(verbose=False)


root_map_dict = {
    'PCANA': [-2.0],  # anoxic stress
    'PCREF': [-5.0],  # ottimo
    'PCWLT': [-15.0],  # wilting point
    'ZROOT': [0.05],  # roots depth [m]
    'PZ': [-1.0],
    'OMGC': [1.0],
}

simuInf3D.update_soil(FP_map=root_map_dict)
#simuInf3D.show_input('root_map')

# =============================================================================
# 3. Hydraulic soil properties (permeability, porosity, Mualem - Van Genuchten parameters)
# =============================================================================
# STEP 1: Create zonation
ny, nx = simuInf3D.DEM.shape
zones = np.ones(np.shape(simuInf3D.DEM))

# central compacted column (0.4 m wide)
zones[:, 23:27] = 2
simuInf3D.update_zone(zones)

#simuInf3D.show_input('zone')

# compacted layers control:
cum_depth = np.cumsum(zr_thickness)
n_compacted_layers = int(np.searchsorted(cum_depth, 0.30, side='right') + 1)
print(f"Layer compattati: 1..{n_compacted_layers} "
      f"(profondita' 0-{cum_depth[n_compacted_layers - 1]:.2f} m)")

# STEP 2: change soil parameters per layers
from pyCATHY.petro import predict_unsat_soil_hydr_param
granulo = [30, 55, 15]
SPP_map_predicted = predict_unsat_soil_hydr_param([granulo])
print("Parameters estimated from soil texture (pedotransfer):")
print(SPP_map_predicted)
df_SPP_map = simuInf3D.init_soil_SPP_map_df(nzones=2, nstr=nstr)

# volume di influenza della compattazione va a circa 30 cm, considerando l'abbassamento del suolo di circa 10cm, abbiamo quasi 1/3 (misurato in campagna)
phi, phi_comp = 0.41, 0.35 # 
elstor = 1e-05

cols = ['PERMX', 'PERMY', 'PERMZ',  # hydraulic conductivity Ks
        'ELSTOR',  # specific storage
        'POROS',  # porosity
        'VGPSATCELL',  # = 1/alpha [m]
        'VGNCELL',  # n
        'VGRMCCELL']  # residual moisture content (theta_r)
idx = pd.IndexSlice


# zona 1
#df_SPP_map.loc[idx[:, :], cols] = [2e-06, 2e-06, 2e-06, elstor, phi, 0.4, 1.5, 0.08] #: homogeneous

# Layer 1: driest
df_SPP_map.loc[idx[:, 1:6], cols] = [2.5e-06, 2.5e-06, 2.8e-06, elstor, phi, 0.3, 1.4, 0.07] 
# Layer 2: intermediate
df_SPP_map.loc[idx[:, 7:10], cols] = [3.3e-06, 3.3e-06, 3e-06, elstor, phi, 0.5, 1.5, 0.06] 
# Layer 3: wettest
df_SPP_map.loc[idx[:, 10:], cols] = [1.e-06, 1.e-06, 1.e-06, elstor, phi, 0.4, 1.4, 0.07] 

# zona 2: compacted
df_SPP_map.loc[idx[2, 1:6], cols] = [3e-07, 3e-07, .5e-08, elstor, phi_comp, 0.55, 1.28, 0.045]


simuInf3D.update_soil(
    SPP_map=df_SPP_map,
    FP_map=root_map_dict,
    PMIN=-1e12,  # Minimum pressure head (very dry condition)
)

### 3.3 mapping to layers
ltop, lbot = mt.get_layer_depths(simuInf3D.dem_parameters)
ic_layers_depth_range = [[0.0, 0.20],[0.2, 0.44],[0.44, 2.0]]
ic_layers_values = [[-1.4],[-1.13],[-1.13]]
                                                                         
ic_array = mt.map_param_to_layers(ltop, lbot,
                                   ic_layers_depth_range,
                                   ic_layers_values)

ic_nodes = simuInf3D.map_prop_2mesh_markers('ic', ic_array, to_nodes=True)
ic_nodes = np.array(ic_nodes, dtype=float).flatten()

simuInf3D.update_ic(INDP=1,
                     IPOND=0,
                     pressure_head_ini=ic_nodes,
                     )
    
# =============================================================================
# 5. Atmospheric forcing (ET + irrigation)
# =============================================================================
HOUR = 3600  # seconds
DAY = 24 * HOUR

t_end = 3 * DAY  # simulation length
times = np.arange(0, t_end + HOUR, HOUR).tolist()

et = -1e-7  # [m/s] constant ET

# Irrigazione:
rain_total_m = 0.04
irr = (rain_total_m / 2) / HOUR

irr_start_hour = 2 
# irrigation area (2m x 2m central)
nx = DEM.shape[1]
ny = DEM.shape[0]

center_x = nx // 2
center_y = ny // 2

half_size = int(2.0 / delta_x / 2) 

x_min = center_x - half_size
x_max = center_x + half_size
y_min = center_y - half_size
y_max = center_y + half_size

# costruzione forzante
atm_dem_steps = []

for i, t in enumerate(times):
    atm = np.ones((DEM.shape[0] + 1, DEM.shape[1] + 1)) * et
    if t == irr_start_hour * HOUR:   
        atm[y_min:y_max, x_min:x_max] = irr
    atm_dem_steps.append(atm)

### update atmospheric boundary conditions                                                                                                                                                              
simuInf3D.update_atmbc(
    HSPATM=0,           # prescribed flux
    IETO=1,
    time=times,
    netValue=atm_dem_steps
)

### save VTK every hour
vtk_times = np.arange(0, t_end + HOUR, HOUR).tolist()
simuInf3D.update_parm(TIMPRTi=vtk_times,
                       TRAFLAG=0,
                       VTKF=2,
                       )

### Boundary conditions
# Close lateral and bottom boundaries (no-flow), leaving only atmospheric fluxes active.
simuInf3D.update_nansfneubc(no_flow=True)
simuInf3D.update_nansfdirbc(no_flow=True)
simuInf3D.update_sfbc(no_flow=True)


# =============================================================================
# 6. Run
# =============================================================================
print("\n" + "=" * 70)
print("RUNNING INFILTRATION SIMULATION")
print("=" * 70)

simuInf3D.run_processor(
    IPRT1=4,
    DTMAX=1e3,
    DELTAT=1e-1,
    DTMIN=1e-3,
    verbose=True
)
print("Simulation complete!")


# =============================================================================
# 7. Pick-up info at nodes/sensor e salvataggio
# =============================================================================
figpath = os.getcwd() + '/' + prj_name + '/figures/'
os.makedirs(figpath, exist_ok=True)

# centro mesh
x0 = x / 2
y0 = y / 2

# profondita' target
z_15 = elevation - 0.15
z_30 = elevation - 0.30
z_50 = elevation - 0.50

node_15, pos_15 = simuInf3D.find_nearest_node([x0, y0, z_15])
node_30, pos_30 = simuInf3D.find_nearest_node([x0, y0, z_30])
node_50, pos_50 = simuInf3D.find_nearest_node([x0, y0, z_50])

print(f"15 cm: node {node_15[0]}, pos = {pos_15[0]}, depth ~ {elevation - pos_15[0][2]:.3f} m")
print(f"30 cm: node {node_30[0]}, pos = {pos_30[0]}, depth ~ {elevation - pos_30[0][2]:.3f} m")
print(f"50 cm: node {node_50[0]}, pos = {pos_50[0]}, depth ~ {elevation - pos_50[0][2]:.3f} m")

### 7.1 Output: saturazione to VWC
df_sw, sw_time = simuInf3D.read_outputs('sw')
t_mod = np.asarray(sw_time, dtype=float) 

theta_15_np = (df_sw.iloc[:, node_15[0]] * phi_comp).to_numpy(dtype=float)
theta_30_np = (df_sw.iloc[:, node_30[0]] * phi).to_numpy(dtype=float)
theta_50_np = (df_sw.iloc[:, node_50[0]] * phi).to_numpy(dtype=float)

print("Saturation output:", df_sw.shape)
print(f"Model time: {t_mod.min() / HOUR:.2f} -> {t_mod.max() / HOUR:.2f} h")

### 7.2 Output: pressure head
df_psi = simuInf3D.read_outputs('psi')
psi_15_np = df_psi.iloc[:, node_15[0]].to_numpy(dtype=float)
psi_30_np = df_psi.iloc[:, node_30[0]].to_numpy(dtype=float)
psi_50_np = df_psi.iloc[:, node_50[0]].to_numpy(dtype=float)

### 7.3 SAVE
model_output = {
    'compLevel': compLevel,
    't_mod': t_mod,                 # secondi dall'inizio SIMULAZIONE (t=0)
    'theta_15': theta_15_np,
    'theta_30': theta_30_np,
    'theta_50': theta_50_np,
    'psi_15': psi_15_np,
    'psi_30': psi_30_np,
    'psi_50': psi_50_np,
    'irr_start_hour': irr_start_hour,  # ore simulate a cui inizia l'infiltrazione
    'HOUR': HOUR,
    'phi': phi,
    'phi_comp': phi_comp,
}

outfile = figpath + f'model_output_{compLevel}.pkl'
with open(outfile, 'wb') as f:
    pickle.dump(model_output, f)

print(f"\nOutput model saved in: {outfile}")
print("Run compare_meas_vs_mod.py for comparison")
