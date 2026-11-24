#!/usr/bin/env python
# coding: utf-8

# # Agripolis 3D hydro
# This example shows how to use pyCATHY to run the hydrological model on a small scale compacted soil scenario illustrated below:

# <img src="../figures/scheme_NOcomp.png" width="400">

# ## 1. Inizializzazione simulazione CATHY
# definiamo unità temporali, inizializziamo il progetto e definiamo i dominio spaziale

# In[21]:


import sys
import os
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
# SETUP: Define paths and experiment parameters
# =============================================================================
# Soil composition: [% sand, % silt, % clay]
granulo = [25, 55, 20] # 25% sand, 55% silt, 20% clay = silty-loam soil
# Soil compaction level: how densely packed the soil is
compLevel = 'no_comp'
# Define where simulation files will be stored
path2prj = "../simulations/"

# =============================================================================
# STEP 1: Initialize the CATHY hydrological model
# =============================================================================
# Create a unique project name
prj_name = 'SimuID_' + (f'{compLevel}')
# Create the CATHY simulation object
simuInf3D = cathy_tools.CATHY(dirName=path2prj, prj_name=prj_name)
# simuInf3D.run_preprocessor(verbose=True)


# In[22]:


# ============================================================
# STEP 2: Create the 3D computational mesh
# ============================================================
# ---------- discretizzazione orizzontale
delta_x = 0.1        # [m] risoluzione x
delta_y = 0.1        # [m] risoluzione y
x = 2.0              # [m] estensione x
y = 2.0              # [m] estensione y
maxdepth = 2.0       # [m] profondità
elevation = 1e-3     # [m] quota

ncellx = int(np.round(x / delta_x))
ncelly = int(np.round(y / delta_y))

# DEM piatto con un outlet leggermente più basso
DEM = np.ones((ncelly, ncellx)) * elevation
DEM[-1, -1] = elevation - 1e-4

# Vertical discretization – focus on 0–60 cm
zr_thickness = [
    0.03, 0.03, 0.03, 0.03,   # 0–12 cm
    0.04, 0.04, 0.04,        # 12–24 cm
    0.06, 0.06,              # 24–36 cm
    0.08, 0.08,              # 36–52 cm
    0.10,                    # 52–62 cm
    0.25,                    # 62–87 cm
    0.40,                    # 87–127 cm
    0.73                     # 127–200 cm
]

zratio = [z / maxdepth for z in zr_thickness]
assert np.isclose(sum(zratio), 1.0), "zratio does not sum to 1"

# ---------- update input preprocessor
simuInf3D.update_prepo_inputs(
    DEM,
    delta_x=delta_x,
    delta_y=delta_y,
    dr=delta_x,
    zratio=zratio,
    base=maxdepth,
)

cell_size = delta_x * delta_y # area cella


simuInf3D.show_input('dem')
# simuInf3D.update_dem_parameters()
# simuInf3D.run_preprocessor(IPRT1=3)
simuInf3D.create_mesh_vtk(verbose=False)

# ## 2. Zonazione, vegetazione e irrigazione
# - non c'è vegetazione: imposto erba omogenea superficiale ma dopo posso rimuoverla
# - unica zona, a cui posso però variare i parametri nei vari layer in profondità

# In[23]:


# --- Vegetation map
root_map, _ = simuInf3D.read_inputs('root_map')
root_map[:, :] = 1 # tutta l’area coperta da erba (classe 1)

simuInf3D.update_veg_map(root_map)
# simuInf3D.update_cathyH(MAXVEG=1)

root_map_dict = {
    'PCANA': [-2.0],     # stress anossico
    'PCREF': [-5.0],     # ottimo
    'PCWLT': [-15.0],    # wilting
    'ZROOT': [0.05],      # profondità radici erba
    'PZ': [-1.0],
    'OMGC': [1.0],
}

simuInf3D.update_soil(FP_map = root_map_dict)
simuInf3D.show_input('root_map')


# definizione di permeabilità, porosità e parametri di Mualem–van Genuchten

# In[24]:


SPP, FP = simuInf3D.read_inputs('soil')#, NVEG = 1)


phi = 0.45

SPP['PERMZ'] = 4e-6
SPP['PERMX'] = 4e-6
SPP['PERMY'] = 4e-6
SPP['ELSTOR'] = 1e-5
SPP['POROS'] = phi
SPP['VGPSATCELL'] = 0.28 # 0.18 
SPP['VGNCELL'] = 1.4 # 1.6
SPP['VGRMCCELL'] = 0.09 # 0.07 

simuInf3D.show_input('zone')
SPP.head()

#%%
df_SPP_map = simuInf3D.init_soil_SPP_map_df(nzones=1, nstr=15)

cols = ['PERMX','PERMY','PERMZ','ELSTOR','POROS']
idx = pd.IndexSlice

# Layers 1–4: compacted
df_SPP_map.loc[idx[:,:], cols] = [4e-6, 4e-6, 4e-6, 1e-5, 0.33]


import pandas as pd

cols_vgp = ['VGPSATCELL','VGNCELL','VGRMCCELL']
cols = ['PERMX','PERMY','PERMZ','ELSTOR','POROS','VGPSATCELL','VGNCELL','VGRMCCELL']

idx = pd.IndexSlice

cols = ['PERMX','PERMY','PERMZ','ELSTOR','POROS','VGPSATCELL','VGNCELL','VGRMCCELL']
df_SPP_map.loc[idx[:, :], cols] = [4e-6, 4e-6, 2e-6, 1e-5, 0.43, 0.18, 1.6, 0.07]


# Layer 1: driest
df_SPP_map.loc[idx[:, 1:5], cols_vgp] = [0.20, 1.6, 0.07]

# Layer 2: intermediate
df_SPP_map.loc[idx[:, 6:10], cols_vgp] = [0.20, 1.7, 0.09]

# Layer 3: wettest
df_SPP_map.loc[idx[:, 10:], cols_vgp] = [0.20, 1.8, 0.1]



# # Layer 1: driest
# df_SPP_map.loc[idx[:, 1:5], cols_vgp] = [0.15, 1.3, 0.05]

# # Layer 2: intermediate
# df_SPP_map.loc[idx[:, 6:10], cols_vgp] = [0.15, 1.3, 0.05]

# # Layer 3: wettest
# df_SPP_map.loc[idx[:, 10:], cols_vgp] = [0.15, 1.3, 0.05]

# mt.get_layer_depths(simuInf3D.dem_parameters)


# df_SPP_map.loc[idx[:, :], cols] = [4e-6, 4e-6, 2e-6, 1e-5, 0.43, 0.18, 1.6, 0.07]
# # zona 2: compacted
# df_SPP_map.loc[idx[2, 1:4], cols] = [5e-7, 5e-7, 5e-7, 1e-5, 0.33, 0.15, 1.3, 0.08]



# ff
# In[25]:


simuInf3D.update_soil(
    SPP=df_SPP_map,
    FP_map=root_map_dict,
    #MAXVEG=1,
    PMIN=-1e12, # Minimum pressure head (very dry condition)
    # show=True,
)
# SPP.head()


# In[26]:

# simuInf3D.run_preprocessor(verbose=True)
# simuInf3D.create_mesh_vtk()

### --- Set initial conditions
# INDP=0 means we're specifying pressure head (not saturation)
# simuInf3D.run_preprocessor(verbose=True)

#%%

# ltop, lbot = mt.get_layer_depths(simuInf3D.dem_parameters)
# ic_layers_depth_range = [[0,0.20],[0.21,0.49],[0.50,2]]
# ic_layers_values = [[-5],[-5],[-2]]

# ic_array = mt.map_param_to_layers(ltop, lbot, 
#                                   ic_layers_depth_range,
#                                   ic_layers_values
#                                   )


# ic_nodes = simuInf3D.map_prop_2mesh_markers('ic',
#                                         ic_array, 
#                                         to_nodes=True
#                                         )

# simuInf3D.update_ic(INDP=1,
#                     IPOND=0,
#                     pressure_head_ini=ic_nodes,
#                     )

simuInf3D.update_ic(INDP=3, WTPOSITION=0.1)

# simuInf3D.read_inputs('ic')
# np.unique(ic_nodes)

# p = pv.Plotter()
# p.add_mesh(simuInf3D.mesh_pv_attributes, 
#            scalars='ic_nodes',
#            cmap="viridis", 
#            show_edges=True
#            )
# p.add_axes()          
# p.show_grid()         
# p.show()

# ## 3. Forzanti atmosferiche
# definizione delle atmosferic boundary cobnditions, con ET negativa (perdita) e irrigazione positiva (input)

# In[27]:


HOUR = 3600            # seconds
DAY  = 24 * HOUR

### --- simulation (1-hour steps, 2 days)
t_end = 3 * DAY
times = np.arange(0, t_end + HOUR, HOUR).tolist()

# Evapotranspiration [m/s] (negative = upward flux)
# et = -2e-7             
et = 0          

# Irrigation: 200 mm in 30 minutes
# → equivalent to 100 mm in 1 hour to match hourly BC grid
rain_total_m = 0.1       # 0.1 = 100 mm in 30 min
irr = (rain_total_m / 2) / HOUR   # = 0.05 / 3600 ≈ 1.39e-5 m/s

irr_start_hour = 1     # irrigation starts after 1 hour

### --- Irrigazione controllata
atm_dem_steps = []

for i, t in enumerate(times):
    atm = np.ones((DEM.shape[0] + 1, DEM.shape[1] + 1)) * et
    if t == irr_start_hour * HOUR:     # irrigazione al tempo desiderato
        atm[:, :] = irr

    atm_dem_steps.append(atm)

### --- Update atmospheric boundary conditions
simuInf3D.update_atmbc(
    HSPATM=0,           # prescribed flux
    IETO=1,             # Evapotranspiration option (0 = no ET)
    time=times,         # times when BC can change
    netValue=atm_dem_steps
)

### --- Save VTK every hour, only during the first day
vtk_times = np.arange(0, DAY + HOUR, HOUR).tolist()
simuInf3D.update_parm(TIMPRTi=vtk_times,
                      TRAFLAG=0,             # Transport flag (0 = no solute transport)
                      VTKF=2,                # VTK output frequency (for 3D visualization)
                     )

### --- Boundary conditions
# Close lateral and bottom boundaries (no-flow),
# leaving only atmospheric fluxes active
simuInf3D.update_nansfneubc(no_flow=True)
simuInf3D.update_nansfdirbc(no_flow=True)
simuInf3D.update_sfbc(no_flow=True)


# ## 4. Run del modello
# time stepg adattivo e output completi (pressione, saturazione, RWU)

# In[28]:



### --- Run the simulation
print("\n" + "="*70)
print("RUNNING INFILTRATION SIMULATION")
print("="*70)

simuInf3D.run_processor(
    IPRT1=4,            # Print interval option
    DTMAX=1e3,          # Maximum time step (seconds)
    DELTAT=1e-1,        # Initial time step (seconds)
    DTMIN=1e-3,         # Minimum time step (seconds)
    verbose=True        # Print progress messages
)

print("Simulation complete!")


# ## 5. Plot e visualizzazione
# estrazione info da punti nello spazio per leggere l'evoluzione temporale del modello (pressione, saturazione, contenuto d’acqua)

# In[29]:


#print(os.getcwd())

# Create a folder to save figures/plots
figpath = os.getcwd()+'/'+prj_name+'/figures/'
os.makedirs(figpath, exist_ok=True)
figpath


# In[30]:


# centro mesh
x0 = x / 2
y0 = y / 2

# profondità target
z_15 = elevation - 0.15  # 15 cm
z_30 = elevation - 0.30  # 30 cm
z_50 = elevation - 0.47  # 50 cm

# nodi più vicini a ciascuna profondità
node_15, pos_15 = simuInf3D.find_nearest_node([x0, y0, z_15])
node_30, pos_30 = simuInf3D.find_nearest_node([x0, y0, z_30])
node_50, pos_50 = simuInf3D.find_nearest_node([x0, y0, z_50])

print(f"15 cm: node {node_15[0]}, pos = {pos_15[0]}, depth ≈ {elevation - pos_15[0][2]:.3f} m")
print(f"30 cm: node {node_30[0]}, pos = {pos_30[0]}, depth ≈ {elevation - pos_30[0][2]:.3f} m")
print(f"50 cm: node {node_50[0]}, pos = {pos_50[0]}, depth ≈ {elevation - pos_50[0][2]:.3f} m")

pl = pv.Plotter(notebook=False
               )
cplt.show_vtk(unit="pressure",
              timeStep=3,
              path=Path(simuInf3D.workdir) / simuInf3D.project_name / 'vtk',
              style='wireframe',
              opacity=0.1,
              ax=pl,
              )
pl.add_points(pos_15[0], color = 'red', point_size=10, render_points_as_spheres=True)
pl.add_points(pos_30[0], color = 'k', point_size=10, render_points_as_spheres=True)
pl.add_points(pos_50[0], color = 'b', point_size=10, render_points_as_spheres=True)
pl.show()
# In[31]:


# Read water saturation results (fraction of pore space filled with water)
# sw = saturation, ranges from 0 (completely dry) to 1 (fully saturated)
df_sw, _ = simuInf3D.read_outputs('sw')
print(df_sw.shape)

# plot the saturation with time at a given mesh point
fig, ax = plt.subplots(figsize=(7, 3))

ax.plot(df_sw.iloc[:, node_15[0]], marker='o', markersize=3,
        label=f"depth: {pos_15[0][2]:.2f} m")

ax.plot(df_sw.iloc[:, node_30[0]], marker='o', markersize=3,
        label=f"depth: {pos_30[0][2]:.2f} m")

ax.plot(df_sw.iloc[:, node_50[0]], marker='o', markersize=3,
        label=f"depth: {pos_50[0][2]:.2f} m")

plt.legend(loc='best')
ax.set_xlabel('Time step', fontsize=12)
ax.set_ylabel('Saturation (-)', fontsize=12)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(figpath+'sat.png', dpi=200, bbox_inches='tight')


# In[32]:


# Read pressure head results (measure of soil water tension)
# Negative values = suction (dry), zero = saturated, positive = ponding
df_psi = simuInf3D.read_outputs('psi')
print(df_psi.shape)

# plot the pressure head with time at a given mesh point
fig, ax = plt.subplots(figsize=(7, 3))
ax.plot(df_psi.iloc[:, node_15[0]], marker='o', markersize=3,
        label=f"depth: {pos_15[0][2]:.2f} m")

ax.plot(df_psi.iloc[:, node_30[0]], marker='o', markersize=3,
        label=f"depth: {pos_30[0][2]:.2f} m")

ax.plot(df_psi.iloc[:, node_50[0]], marker='o', markersize=3,
        label=f"depth: {pos_50[0][2]:.2f} m")

plt.legend(loc='best')
ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel('Pressure head (m)', fontsize=12)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(figpath+'psi.png', dpi=200, bbox_inches='tight')


# In[33]:


from pyCATHY.plotters import cathy_plots as cplt

# timestep da visualizzare
time_index = 3
time_value = vtk_times[time_index]  # tempo corrispondente

pl = pv.Plotter()
cplt.show_vtk(
    unit="saturation",
    timeStep=time_index,
    path=Path(simuInf3D.workdir) / simuInf3D.project_name / 'vtk',
    ax=pl,
    clim=[0, 1]
)

pl.add_text(f"Timestep: {time_index}  |  {time_value/3600:.1f} hours", 
            position='upper_right', font_size=11, color='black', shadow=False)

pl.show_grid()
pl.show()


# In[34]:


figpath


# In[35]:


# SNAPSHOTS
output_dir = figpath+'vtk_snapshots/'
os.makedirs(output_dir, exist_ok=True)

vtk_path = Path(simuInf3D.workdir) / simuInf3D.project_name / 'vtk'

vtk_files = sorted(list(vtk_path.glob("1*.vtk")))  
n_steps = min(len(vtk_files), len(vtk_times))
print(f"Trovati {n_steps} timesteps")

# Loop
for t in range(n_steps):
    pl = pv.Plotter(off_screen=True)

    # Mostra il VTK
    cplt.show_vtk(
        unit="saturation",       
        timeStep=t,              
        path=vtk_path,
        ax=pl,
        clim=[0, 1],  
        #cmap='viridis'
    )

    pl.add_text(
        f"Timestep: {t}  |  {vtk_times[t]/3600:.1f} h",
        position='upper_right',
        font_size=11,
        color='black',
        shadow=False
    )

    pl.show_grid()

    # salva screenshot
    filename = output_dir + f"sat_t{t:03d}.png"
    pl.screenshot(str(filename))
    pl.close()  # chiude il plotter per liberare memoria

print("Tutti i PNG sono stati salvati in:", output_dir)


# In[36]:


from pathlib import Path
import imageio
images = sorted(Path(output_dir).glob("*.png"))
imageio.mimsave(figpath + "vwc.gif",
                [imageio.imread(img) for img in images],
                fps=2)


# ## 6. Confronto
# confronto dati VWC modellati con quelli dei sensori puntuali (FDR, tensiometri)

# In[37]:


print(os.getcwd())


# In[38]:


datadir = '../data/sensors/'
file = "st5.csv"

df = pd.read_csv(datadir+file, header=[0,1,2])

# --- Unisci i livelli delle colonne e rimuovi '_' multipli
df.columns = ['_'.join([str(lv) for lv in col if lv is not None]).strip() for col in df.columns]
df.columns = [col.replace('___','').replace('__','').strip() for col in df.columns]

# --- Estrai la colonna del tempo (la prima)
df['time'] = pd.to_datetime(df.iloc[:,0], dayfirst=True)

# --- Calcola ore dall'inizio
t0 = df['time'].iloc[0]
df['time_hours'] = (df['time'] - t0).dt.total_seconds() / 3600
df['time_seconds'] = (df['time'] - t0).dt.total_seconds()

# --- Trova le colonne VWC dei Port 1, 3 e 5
vwc_ports = []
for port in ['Port 1','Port 3','Port 5']:
    for col in df.columns:
        if port in col and 'VWC' in col:
            vwc_ports.append(col)

print("Colonne VWC selezionate:", vwc_ports)

# --- Crea un DataFrame solo con tempo e VWC dei sensori selezionati
df_vwc = df[['time_hours'] + ['time_seconds'] + vwc_ports]

# --- Controllo
df_vwc.head()


# In[39]:


fig, ax = plt.subplots(figsize=(8,4))

ax.plot(df_vwc['time_seconds'], df_vwc[vwc_ports[0]], marker='o', markersize=4, alpha = .6, label=f"Sensor 1")
ax.plot(df_vwc['time_seconds'], df_vwc[vwc_ports[1]], marker='o', markersize=4, alpha = .6, label=f"Sensor 2")
ax.plot(df_vwc['time_seconds'], df_vwc[vwc_ports[2]], marker='o', color='k', markersize=4, alpha = .6, label=f"Sensor 3")

ax.set_xlabel("Time (s)")
ax.set_ylabel("VWC (-)")
ax.grid(True, alpha=0.5)
ax.set_title("Measured VWC at sensor locations")
ax.legend()
plt.show()


# In[42]:


# porosity
phi = 0.38

# filter df
df_vwc_window = df_vwc[
    (df_vwc['time_seconds'] >= 0.482e6) &
    (df_vwc['time_seconds'] <= .6e6)
].copy()

df_vwc_window['time_seconds'] = (df_vwc_window['time_seconds'] - df_vwc_window['time_seconds'].iloc[0])

fig, ax = plt.subplots(figsize=(9, 3.5))

# measured VWC
ax.plot(df_vwc_window['time_seconds'], df_vwc_window[vwc_ports[0]],
        marker='o', markersize=6, alpha=.5, label="meas 1")

ax.plot(df_vwc_window['time_seconds'], df_vwc_window[vwc_ports[1]],
        marker='o', markersize=6, alpha=.5, label="meas 2")

ax.plot(df_vwc_window['time_seconds'], df_vwc_window[vwc_ports[2]],
        marker='o', color='k', markersize=6, alpha=.5, label="meas 3")


# modeled sat. to VWC
theta_15 = df_sw.iloc[:, node_15[0]] * phi
theta_30 = df_sw.iloc[:, node_30[0]] * phi
theta_50 = df_sw.iloc[:, node_50[0]] * phi

ax.plot(theta_15, marker='s', color='tab:blue',
        markersize=4, label="mod 1")

ax.plot(theta_30, marker='s', color='tab:orange',
        markersize=4, label="mod 2")

ax.plot(theta_50, marker='s', color='k',
        markersize=4, label="mod 3")


#ax.set_ylim(0.15, 0.45)
ax.set_xlabel("Time (s)")
ax.set_ylabel("VWC (-)")
ax.set_title("Measured vs Modeled at sensor locations: " + compLevel)
ax.grid(True, alpha=0.5)
ax.legend()

fig.savefig(figpath + 'measVSmod_' + compLevel + '.png', dpi=200, bbox_inches='tight')
plt.show()


# In[ ]:




