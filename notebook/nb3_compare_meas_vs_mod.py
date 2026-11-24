# =============================================================================
# Confronto misura vs modello
# =============================================================================
import os
import json
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================================================================
# 0. Config
# =============================================================================
compLevel = 'heavy_comp'  # change it between 'no_comp' and 'heavy_comp'

if compLevel == 'heavy_comp':
    exp_label = 'HC'
    sensor_file = 'st1.csv'
else:
    exp_label = 'NC'
    sensor_file = 'st5.csv'   

import project_paths as PP                  # data_path.toml: every folder of the project
path2prj = PP.get().simulations_str          # [outputs].simulations
prj_name = 'SimuID_' + compLevel
figpath = os.path.join(path2prj, prj_name, 'figures') + '/'
datadir = PP.get().sensors_dir + os.sep   # [inputs].sensors_dir

BASELINE_LOOKBACK_HOURS = 6
EVENT_WINDOW_HOURS = 24 * 3

# =============================================================================
# 1. Carica output modello e parametri di sito
# =============================================================================
with open(figpath + f'model_output_{compLevel}.pkl', 'rb') as f:
    model_output = pickle.load(f)

t_mod = model_output['t_mod']
theta_15_np = model_output['theta_15']
theta_30_np = model_output['theta_30']
theta_50_np = model_output['theta_50']
psi_15 = model_output['psi_15']
psi_30 = model_output['psi_30']
psi_50 = model_output['psi_50']
HOUR = model_output['HOUR']
irr_start_hour = model_output['irr_start_hour']

with open(PP.get().site_params) as f:   # [inputs].site_params
    site_params = json.load(f)

INFILTRATION_START = pd.Timestamp(site_params['INFILTRATION_START'][exp_label])
print(INFILTRATION_START)
# =============================================================================
# 2. Carica i dati sensori e allinea il tempo a quello del modello
# =============================================================================
# - 'Port 1': shallow (15cm)
# - 'Port 2': medium (30cm)
# - 'Port 3': deep (50cm)
df = pd.read_csv(datadir + sensor_file, header=[0, 1, 2])
df.columns = ['_'.join([str(lv) for lv in col if lv is not None]).strip() for col in df.columns]
df.columns = [c.replace('___', '').replace('__', '').strip() for c in df.columns]

df['datetime'] = pd.to_datetime(df.iloc[:, 0], dayfirst=True)  # real time between sensors

# RIFERIMENTO TEMPORALE DEL MODELLO
# CATHY considera t=0 come l'inizio della simulazione.
#   t = 0                    -> inizio simulazione
#   t = irr_start_hour*HOUR  -> inizio irrigazione
SIM_START_REAL = INFILTRATION_START - pd.Timedelta(seconds=irr_start_hour * HOUR)

print("Riferimento temporale:")
print(f"  Inizio simulazione : {SIM_START_REAL}")
print(f"  Inizio irrigazione : {INFILTRATION_START}")
print(f"  Irrigation model t : {irr_start_hour * HOUR:.0f} s")

df['time_seconds'] = (df['datetime'] - SIM_START_REAL).dt.total_seconds()
df['time_hours'] = df['time_seconds'] / HOUR

# FINESTRA TEMPORALE DEL CONFRONTO
t_win_start = irr_start_hour * HOUR - BASELINE_LOOKBACK_HOURS * HOUR
t_win_end = irr_start_hour * HOUR + EVENT_WINDOW_HOURS * HOUR

vwc_ports = []
for port in ['Port 1', 'Port 2', 'Port 3']:
    for col in df.columns:
        if port in col and 'VWC' in col:
            vwc_ports.append(col)
print("\nColonne VWC selezionate:", vwc_ports)

df_vwc = df[['datetime', 'time_hours', 'time_seconds'] + vwc_ports].copy()
df_vwc_window = df_vwc[(df_vwc['time_seconds'] >= t_win_start) &
                        (df_vwc['time_seconds'] <= t_win_end)].copy()

print(f"\nFinestra confronto: {t_win_start / HOUR:.2f} -> {t_win_end / HOUR:.2f} h")
print(f"Misure nella finestra: {len(df_vwc_window)}")

print("\nControllo allineamento temporale:")
print(f"  Modello: {t_mod.min() / HOUR:.2f} -> {t_mod.max() / HOUR:.2f} h")
print(f"  Misure : {df_vwc_window['time_seconds'].min() / HOUR:.2f} -> "
      f"{df_vwc_window['time_seconds'].max() / HOUR:.2f} h")

# =============================================================================
# 3. Plot misura vs modello
# =============================================================================
fig, ax = plt.subplots(figsize=(9, 3.5))

ax.plot(df_vwc_window['time_hours'], df_vwc_window[vwc_ports[0]],
        marker='o', markersize=6, alpha=.5, label="meas 1 (15cm)")
ax.plot(df_vwc_window['time_hours'], df_vwc_window[vwc_ports[1]],
        marker='o', markersize=6, alpha=.5, label="meas 2 (30cm)")
ax.plot(df_vwc_window['time_hours'], df_vwc_window[vwc_ports[2]],
        marker='o', color='k', markersize=6, alpha=.5, label="meas 3 (50cm)")

time_model_hours = t_mod / HOUR
ax.plot(time_model_hours, theta_15_np, marker='s', color='tab:blue', markersize=4, label="mod 1 (15cm)")
ax.plot(time_model_hours, theta_30_np, marker='s', color='tab:orange', markersize=4, label="mod 2 (30cm)")
ax.plot(time_model_hours, theta_50_np, marker='s', color='k', markersize=4, label="mod 3 (50cm)")

ax.axvline(irr_start_hour, color='k', linestyle='--', linewidth=1, alpha=0.6, label='inizio infiltrazione')
ax.set_xlim(t_win_start / HOUR, t_win_end / HOUR)
ax.set_xlabel("Time since simulation start (h)")
ax.set_ylabel(r"VWC ($\mathregular{m^3\,m^{-3}}$)")
ax.set_title("Measured vs Modelled at sensor locations: " + compLevel)
ax.grid(True, alpha=0.5)
ax.legend(loc='lower right', ncol=2, fontsize=8)
fig.tight_layout()
fig.savefig(figpath + 'measVSmod_' + compLevel + '.png', dpi=200, bbox_inches='tight')
plt.show()

# =============================================================================
# 4. Interpolazione modello sui tempi di misura + CSV di confronto
# =============================================================================
# NB: filtriamo esplicitamente le misure al dominio effettivo del modello [t_mod.min(), t_mod.max()] prima di interpolare.

mod_t_min, mod_t_max = t_mod.min(), t_mod.max()

in_domain_mask = (df_vwc_window['time_seconds'] >= mod_t_min) & \
                 (df_vwc_window['time_seconds'] <= mod_t_max)
n_out = (~in_domain_mask).sum()
if n_out > 0:
    print(f"\nATTENZIONE: {n_out} misure nella finestra scelta cadono fuori dal "
          f"dominio del modello [{mod_t_min/HOUR:.2f}, {mod_t_max/HOUR:.2f}] h "
          f"(finestra richiesta: [{t_win_start/HOUR:.2f}, {t_win_end/HOUR:.2f}] h). "
          "Questi punti vengono esclusi dal CSV di confronto per evitare valori "
          "modello estrapolati e non rappresentativi.")

df_vwc_compare = df_vwc_window[in_domain_mask].copy()
time_meas = df_vwc_compare['time_seconds'].to_numpy(dtype=float)

theta15_mod = np.interp(time_meas, t_mod, theta_15_np)
theta30_mod = np.interp(time_meas, t_mod, theta_30_np)
theta50_mod = np.interp(time_meas, t_mod, theta_50_np)

psi15_mod = np.interp(time_meas, t_mod, psi_15)
psi30_mod = np.interp(time_meas, t_mod, psi_30)
psi50_mod = np.interp(time_meas, t_mod, psi_50)

comparison_df = pd.DataFrame({
    "time_s": time_meas,
    "theta_meas_15": df_vwc_compare[vwc_ports[0]].to_numpy(),
    "theta_mod_15": theta15_mod,
    "theta_meas_30": df_vwc_compare[vwc_ports[1]].to_numpy(),
    "theta_mod_30": theta30_mod,
    "theta_meas_50": df_vwc_compare[vwc_ports[2]].to_numpy(),
    "theta_mod_50": theta50_mod,
    "psi_mod_15": psi15_mod,
    "psi_mod_30": psi30_mod,
    "psi_mod_50": psi50_mod,
})

comparison_df.to_csv(figpath + f"results_{compLevel}.csv", index=False, float_format="%.6f")
print(f"\nCSV saved in: {figpath}results_{compLevel}.csv")


# =============================================================================
# 5. Goodness-of-fit
# =============================================================================
def compute_gof_metrics(t_meas, y_meas, t_mod_, y_mod_):

    """
    Compute goodness-of-fit metrics between measured and modeled time series.

    RMSE: root mean square error, quantifying the overall magnitude of the
          prediction errors, with greater weight given to larger errors.
    MAE: mean absolute error, representing the average magnitude of the
         prediction errors, independent of their direction.
    Bias: mean signed error, indicating systematic over- or under-prediction.
    R2: coefficient of determination based on the squared Pearson correlation,
        describing the strength of the linear association between measured
        and modeled values.
    NSE: Nash-Sutcliffe Efficiency, evaluating model predictive performance
         relative to using the mean of the observations as a benchmark.
    """

    t_meas = np.asarray(t_meas, dtype=float)
    y_meas = np.asarray(y_meas, dtype=float)
    t_mod_ = np.asarray(t_mod_, dtype=float)
    y_mod_ = np.asarray(y_mod_, dtype=float)

    mask = (t_meas >= t_mod_.min()) & (t_meas <= t_mod_.max())
    t_meas_c, y_meas_c = t_meas[mask], y_meas[mask]
    if len(t_meas_c) < 3:
        return {"RMSE": np.nan, "MAE": np.nan, "Bias": np.nan,
                "R2": np.nan, "NSE": np.nan, "n": len(t_meas_c)}

    y_mod_interp = np.interp(t_meas_c, t_mod_, y_mod_)
    residuals = y_mod_interp - y_meas_c


    rmse = np.sqrt(np.mean(residuals ** 2))
    mae = np.mean(np.abs(residuals))
    bias = np.mean(residuals)
    r = np.corrcoef(y_meas_c, y_mod_interp)[0, 1]
    r2 = r ** 2
    denom = np.sum((y_meas_c - np.mean(y_meas_c)) ** 2)
    nse = 1 - np.sum(residuals ** 2) / denom if denom > 0 else np.nan

    return {"RMSE": rmse, "MAE": mae, "Bias": bias, "R2": r2, "NSE": nse, "n": len(t_meas_c)}


sensor_defs = [
    {"label": "15 cm", "y_mod": theta_15_np, "meas_col": vwc_ports[0]},
    {"label": "30 cm", "y_mod": theta_30_np, "meas_col": vwc_ports[1]},
    {"label": "50 cm", "y_mod": theta_50_np, "meas_col": vwc_ports[2]},
]

rows = []
for s in sensor_defs:
    m = compute_gof_metrics(
        df_vwc_window['time_seconds'].to_numpy(),
        df_vwc_window[s["meas_col"]].to_numpy(),
        t_mod,
        s["y_mod"],
    )
    rows.append({"sensor": s["label"], **m})

gof_table = pd.DataFrame(rows).set_index("sensor")[["n", "RMSE", "MAE", "Bias", "R2", "NSE"]]
print("\nFit modello vs sensori:")
print(gof_table.round(4))
gof_table.to_csv(figpath + f"gof_{compLevel}.csv")