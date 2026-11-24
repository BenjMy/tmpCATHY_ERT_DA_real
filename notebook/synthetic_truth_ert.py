# =============================================================================
# nb6: synthetic truth + synthetic ERT (+ virtual TDR) for the twin experiment
# =============================================================================
# Run from the notebook/ folder (like nb5):
#       python synthetic_truth_ert.py --case NC
#       python synthetic_truth_ert.py --case HC --delta 0.25
#
# 1. TRUTH: the hydrological solution of nb1 (NC) / nb2 (HC), built by model_config.setup_model with the
#    TRUTH parameters, on the same time vector as the DA (ERT survey times of schedule2023_3d.csv).
# 2. TRUTH -> RESISTIVITY with the SAME chain as the DA (pyCATHY.ERT.petro_Archie.SW_2_ERT_ERa_DA):
#       sw (nodes) -> Archie -> CATHY mesh -> CATHY_2_pg (trace_mesh onto the forward mesh, same
#       ERT<->CATHY frame as nb5) -> simulate_ERT.create_ERT_survey_pg (same quadrupoles as the real data).
#    Differences from the DA, on purpose (avoid the "inverse crime"):
#       - Archie parameters of the truth may differ from the DA prior (--archie offset)
#       - optional finer forward mesh for the truth (--fwd-mesh)
#       - nodal porosity is used in Archie. The UNPATCHED library applies the porosity of node 0 to every
#         node (petro_Archie.SW_2_ER0_DA); pycathy_twin_fixes.patch makes the DA nodal too.
#         --porosity-mode da_node0 reproduces the unpatched DA (only for runs without the patch)
# 3. NOISE: 5 % multiplicative log-normal noise on rhoa, seeded, added HERE (pyGIMLi noise is switched off).
# 4. FILES: one .tldat per survey in the format and with the file names of the real data (only the
#    quadrupoles that the DA keeps: err < 900 and the homogeneous-medium forward test of nb5), so
#    nb7 reads them with the nb5 code, ERT_DATA_ERR = 5.
# 5. TRUTH EXPORT (truth_<case>.pkl) and virtual TDR csv (independent validation, not assimilated).
#
# Needs next to it: model_config.py, nc_model.py, ert_schedule.py, ert_frame.py, synth_tools.py.
import argparse
import json
import os
import pickle
import shutil
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd
import pyvista as pv

import model_config as M
import synth_tools as S
import ert_schedule as ES
from ert_frame import choose_ert_frame
from nc_model import (HOUR, X, Y, ELEVATION, ZR_THICK, ERT_COMP_LEVELS, build_times, snap_to_times,
                      load_ert_electrodes, filter_tldat, unreliable_quadrupoles)

# =============================================================================
# 0. Config
# =============================================================================
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--case', default='NC', choices=list(M.CASES))
ap.add_argument('--root', default=None, help='project root (default: data_path.toml)')
ap.add_argument('--delta', type=float, default=None, help='horizontal cell size [m] (default: the case value)')
ap.add_argument('--noise-pct', type=float, default=5.0, help='relative noise on rhoa [%%]')
ap.add_argument('--seed-noise', type=int, default=1234)
ap.add_argument('--archie', default='prior', choices=['prior', 'offset'],
                help="truth Archie (a, m, n): 'prior' = the DA prior mean of nb5, 'offset' = shifted by ~1 sigma")
ap.add_argument('--porosity-mode', default='nodal', choices=['nodal', 'da_node0'])
ap.add_argument('--fwd-mesh', default=None, help='forward mesh of the TRUTH (.vtk); default = the DA mesh')
ap.add_argument('--fwd-tol', type=float, default=0.10, help='homogeneous-medium test tolerance (as nb5 ERT_FWD_TOL)')
ap.add_argument('--dense-step-h', type=float, default=None,
                help='synthetic schedule every N h from the first to the last real survey (default: real times)')
ap.add_argument('--out-dt', type=float, default=300.0, help='dense CATHY output interval [s] of the truth run')
ap.add_argument('--rhoa-ratio', type=float, default=3.0, help='drop quadrupoles with rhoa outside [rho_min/r, rho_max*r] (0 = keep)')
ap.add_argument('--tdr-sigma', type=float, default=0.02, help='virtual TDR noise [m3/m3]')
ap.add_argument('--keep-project', action='store_true', help='do not delete an existing truth project')
args = ap.parse_args()

import project_paths as PP
P = PP.get()
CASE, ROOT = args.case, (args.root or P.root_str)
_TOML = P.is_root(ROOT)
cfg = M.get_case(CASE, args.delta)
COMP_LEVEL = M.COMP_LEVEL[CASE]
_, PREFIX, _ = ERT_COMP_LEVELS[COMP_LEVEL]
IRR_START_HOUR = cfg['irr_start_hour']
STABILIZE_PERIOD = IRR_START_HOUR * HOUR
MAX_SNAP_MIN = 30

# Archie: nb5 constants (prior mean). 'offset' = +1 sigma a, +1 sigma m, -1 sigma n of make_scenario_nc_VGP.py
WATER_EC_uScm = 588
rFLUID = 1 / (WATER_EC_uScm * (1e-6 / 1e-2))
ARCHIE_PRIOR = dict(a=0.3, m=1.7, n=1.7)
ARCHIE_TRUTH = ARCHIE_PRIOR if args.archie == 'prior' else dict(a=0.4, m=1.85, n=1.55)

SITE_PARAMS = P.site_params if _TOML else ROOT + 'data/site_params.json'
SCHEDULE_FILE = P.ert_schedule if _TOML else ROOT + 'data/3d/schedule2023_3d.csv'
path2prj = P.simulations_str if _TOML else ROOT + 'simulations/'
prj_name = f'TRUTH_{CASE}' + (f'_d{args.delta}' if args.delta else '')

paths_real = M.ert_paths(ROOT, CASE)                       # real data (templates), sequence, DA mesh
paths_syn = M.ert_paths(ROOT, CASE, synthetic=True)        # synthetic data_dir
OUT_ROOT = os.path.dirname(paths_syn['data_dir'])          # data/3d/synthetic/<case>
TMPL_DIR = os.path.join(OUT_ROOT, 'template')
for d in (paths_syn['data_dir'], TMPL_DIR):
    os.makedirs(d, exist_ok=True)
DA_MESH, SEQUENCE = paths_real['mesh'], paths_real['sequence']
FWD_MESH = args.fwd_mesh or DA_MESH
for f in (DA_MESH, SEQUENCE, FWD_MESH, SCHEDULE_FILE, SITE_PARAMS):
    assert os.path.isfile(f), f"not found: {os.path.abspath(f)}"
if FWD_MESH == DA_MESH:
    warnings.warn("truth forward on the SAME mesh as the DA: mesh discretisation error is not tested "
                  "(use --fwd-mesh with a finer mesh for the final experiments)")
# create_ERT_survey_pg tries <stem>.msh, then .bms, then .vtk: a sibling file would silently be used
for m_ in {DA_MESH, FWD_MESH}:
    sib = [os.path.splitext(m_)[0] + e for e in ('.msh', '.bms') if os.path.isfile(os.path.splitext(m_)[0] + e)]
    if sib:
        raise RuntimeError(f"{sib} exist next to {m_}: simulate_ERT would read them instead of the .vtk "
                           "(different cell order/frame). Remove/rename them.")

# =============================================================================
# 1. Survey schedule -> CATHY times (same as nb5)
# =============================================================================
INFILTRATION_START = ES.load_infiltration_start(SITE_PARAMS, COMP_LEVEL)
acq = ES.read_schedule_2023(SCHEDULE_FILE, COMP_LEVEL)
ERT_TIMES_H = ES.times_h_since_irrigation(acq, INFILTRATION_START)
if args.dense_step_h:
    t0, t1 = min(ERT_TIMES_H.values()), max(ERT_TIMES_H.values())
    ERT_TIMES_H = {i: round(t0 + i * args.dense_step_h, 6) for i in range(int((t1 - t0) / args.dense_step_h) + 1)}
problems = ES.check_schedule(ERT_TIMES_H, IRR_START_HOUR)
assert not problems, "ERT schedule problems:\n  " + "\n  ".join(problems)
surveys = sorted(ERT_TIMES_H)
ert_times = [ERT_TIMES_H[i] * HOUR + STABILIZE_PERIOD for i in surveys]       # CATHY time [s]
atm_times = build_times(ert_times)
ert_times = snap_to_times(ert_times, atm_times, MAX_SNAP_MIN)
print(f"[{CASE}] {len(surveys)} surveys, CATHY times [h]: {np.round(np.asarray(ert_times) / HOUR, 2).tolist()}")

# =============================================================================
# 2. Truth run (model_config.setup_model with the TRUTH parameters, DA time vector)
# =============================================================================
import project_paths as _PP0                  # data_path.toml: [code].pycathy (clone of pycathy_wrapper) in sys.path
_PP0.get().setup_sys_path()
from pyCATHY import cathy_tools
from pyCATHY import meshtools as mt
from pyCATHY.ERT import simulate_ERT as SE
from pyCATHY.ERT.alignment import transform_points, check_ert_cathy_alignment

prj_path = os.path.join(path2prj, prj_name)
if os.path.exists(prj_path) and not args.keep_project:
    shutil.rmtree(prj_path)
os.makedirs(prj_path, exist_ok=True)
simu = cathy_tools.CATHY(dirName=path2prj, prj_name=prj_name)
DEM, ic_nodes = M.setup_model(simu, CASE, extra_times=ert_times, delta=args.delta)     # truth: no overrides
# CATHY prints at the first time STEP >= each requested time, so a continuous run does not give the state at the
# atmbc times (3627.9 s printed for 3600 s), while the DA windows end exactly there. Ask for dense output
# (every OUT_DT_S, i.e. every step) and interpolate linearly in time to the exact atmbc times.
OUT_DT_S = args.out_dt
dense = np.unique(np.round(np.r_[np.asarray(atm_times, float), np.arange(0.0, float(atm_times[-1]), OUT_DT_S)]).astype(int))  # update_parm casts to int
simu.update_parm(TIMPRTi=dense.tolist())
simu.run_preprocessor()
simu.run_processor(IPRT1=4, DTMAX=1e3, DELTAT=1e-1, DTMIN=1e-3, verbose=False)

df_sw, sw_time_raw = simu.read_outputs('sw')
df_psi = simu.read_outputs('psi')
sw_time_raw = np.asarray(sw_time_raw, dtype=float)
SW_raw, PSI_raw = df_sw.to_numpy(dtype=float), df_psi.to_numpy(dtype=float)
assert SW_raw.shape == PSI_raw.shape, f"sw {SW_raw.shape} and psi {PSI_raw.shape} differ"
SW, _gap = S.interp_in_time(SW_raw, sw_time_raw, atm_times)
PSI, _ = S.interp_in_time(PSI_raw, sw_time_raw, atm_times)
sw_time = np.asarray(atm_times, dtype=float)
print(f"CATHY output: {len(sw_time_raw)} times (largest step {_gap:.0f} s) -> interpolated to the {len(sw_time)} atmbc times")
print(f"truth: {SW.shape[0]} output times x {SW.shape[1]} nodes, sw in [{SW.min():.3f}, {SW.max():.3f}]")

_cathy_vtk = os.path.join(simu.workdir, simu.project_name, 'vtk', simu.project_name + '.vtk')
assert os.path.isfile(_cathy_vtk), f"CATHY mesh not found: {_cathy_vtk}"
meshCATHY = pv.read(_cathy_vtk)
assert meshCATHY.n_points == SW.shape[1], (meshCATHY.n_points, SW.shape[1])
cathy_bounds = np.asarray(meshCATHY.bounds, dtype=float)

# nodal porosity, built exactly as the DA does (pyCATHY.DA.cathy_DA._ert_ensemble_inputs)
soil_df, _ = simu.read_inputs('soil')
_, poro_list = simu.map_prop_2mesh_markers('POROS', soil_df['POROS'].to_list(), to_nodes=False)
poro_nodes = np.asarray(poro_list, dtype=float)
assert poro_nodes.size == SW.shape[1]
poro_archie = poro_nodes if args.porosity_mode == 'nodal' else np.full_like(poro_nodes, poro_nodes[0])
print(f"porosity: nodal range [{poro_nodes.min():.3f}, {poro_nodes.max():.3f}], mode {args.porosity_mode}")

# =============================================================================
# 3. ERT geometry: frame, meshes (same gate as nb5)
# =============================================================================
elecs = load_ert_electrodes(SEQUENCE)
mesh_da_raw = pv.read(DA_MESH)
ERT_SIGNS, ERT_SHIFT, _origin = choose_ert_frame(elecs, mesh_da_raw.points, cathy_bounds)
print(f"ERT frame: {_origin}; signs={ERT_SIGNS}, shift={np.round(ERT_SHIFT, 4).tolist()}")
mesh_fwd_raw = pv.read(FWD_MESH)
fwd_points_cathy = transform_points(mesh_fwd_raw.points, ERT_SIGNS, ERT_SHIFT)
elecs_cathy = transform_points(elecs, ERT_SIGNS, ERT_SHIFT)
check_ert_cathy_alignment(elecs_cathy, cathy_bounds, pg_points=fwd_points_cathy, on_fail='raise')

# =============================================================================
# 4. Templates = the quadrupoles the DA will keep (real data files, nb5 filters)
# =============================================================================
def real_file(i):
    return os.path.join(paths_real['data_dir'], f"{PREFIX}t{int(i):02d}.tldat")


def make_template(src, dst):
    """Valid rows of a real survey (err < 900) minus the unreliable quadrupoles (homogeneous test, DA mesh)."""
    n_valid, n_tot = filter_tldat(src, dst, min_abs_r=0.0)
    bad = unreliable_quadrupoles(DA_MESH, dst, tol=args.fwd_tol, verbose=False)
    tmp = dst + '.tmp'
    n_keep, _ = filter_tldat(dst, tmp, min_abs_r=0.0, drop=bad)
    os.replace(tmp, dst)
    return n_keep, n_valid, n_tot


templates, _cache = {}, {}
available = [i for i in surveys if os.path.isfile(real_file(i))]
assert available, f"no real survey file found in {paths_real['data_dir']}: cannot build the quadrupole templates"
for i in surveys:
    src = real_file(i) if i in available else real_file(available[0])
    if i not in available:
        warnings.warn(f"real file of survey {i} missing: the quadrupole mask of survey {available[0]} is used")
    if src not in _cache:
        dst = os.path.join(TMPL_DIR, os.path.basename(src))
        _cache[src] = (dst, make_template(src, dst))
        print(f"  template {os.path.basename(src)}: {_cache[src][1][0]} kept / {_cache[src][1][1]} valid / {_cache[src][1][2]} rows")
    templates[i] = _cache[src][0]

# =============================================================================
# 5. Forward (truth) + noise + files
# =============================================================================
rng_seed0 = args.seed_noise
tmpdir = tempfile.mkdtemp(prefix='nb6_')
tmp_path = os.path.join(tmpdir, '')                       # add_attribute_2mesh concatenates path + name
meta_fwd = {'forward_mesh_vtk_file': FWD_MESH, 'mesh_nodes_modif': fwd_points_cathy, 'data_format': 'pygimli'}
survey_out, rho_nodes_all = {}, {}
for i, t_s in zip(surveys, ert_times):
    it = S.time_index(sw_time, t_s)                       # the state at the ACQUISITION time (window end)
    rho_nodes = S.archie_rho(SW[it], poro_archie, rFLUID, ARCHIE_TRUTH['a'], ARCHIE_TRUTH['m'], ARCHIE_TRUTH['n'])
    mesh_c, name = mt.add_attribute_2mesh(rho_nodes, meshCATHY.copy(), 'ER_truth', overwrite=True, path=tmp_path)
    mesh_g, scalar_new = mt.CATHY_2_pg(mesh_c, meta_fwd, scalar=name, show=False, path=tmp_path)
    res0 = np.array(mesh_g.get_array(scalar_new), dtype=float)
    assert np.all(np.isfinite(res0)) and np.all(res0 > 0), f"survey {i}: invalid res0"

    tmpl = templates[i]
    het = SE.create_ERT_survey_pg(tmp_path, sequence=tmpl, mesh=FWD_MESH, res0=res0,
                                  fwdNoiseLevel=0, fwdNoiseAbs=0, verbose=False)   # noise: added below, seeded
    rhoa_true = np.asarray(het['rhoa'], dtype=float)
    n_tmpl = S.read_tldat_table(tmpl)['n']
    assert rhoa_true.size == n_tmpl, f"survey {i}: forward gave {rhoa_true.size} data for {n_tmpl} quadrupoles"
    # create_ERT_survey_pg replaces rhoa <= 0 by 1e-3 Ohm m: such quadrupoles are dropped, not written
    # also dropped: rhoa far outside the resistivity range of the model (near-zero-potential quadrupoles of the
    # forward, unstable: seen at 0.3 ohm.m for a model of 23-190 ohm.m); --rhoa-ratio 0 keeps them
    bad_fwd = ~(rhoa_true > 1e-2)
    if args.rhoa_ratio > 0:
        bad_fwd |= (rhoa_true < rho_nodes.min() / args.rhoa_ratio) | (rhoa_true > rho_nodes.max() * args.rhoa_ratio)
    if bad_fwd.any():
        warnings.warn(f"survey {i}: {int(bad_fwd.sum())} quadrupoles with rhoa <= 0.01 or outside the model "
                      f"resistivity range / {args.rhoa_ratio:g}: dropped")
        tmpl2 = os.path.join(tmpdir, f'tmpl_{i}.tldat')
        filter_tldat(tmpl, tmpl2, min_abs_r=0.0, drop=bad_fwd)
        tmpl, rhoa_true = tmpl2, rhoa_true[~bad_fwd]

    rhoa_obs, eps = S.add_lognormal_noise(rhoa_true, args.noise_pct / 100.0, rng_seed0 + 1000 * int(i))
    fname = os.path.join(paths_syn['data_dir'], f"{PREFIX}t{int(i):02d}.tldat")
    S.write_tldat_like(tmpl, fname, rhoa_obs, args.noise_pct / 100.0)
    tab = S.read_tldat_table(fname)
    survey_out[i] = dict(file=fname, t_s=t_s, t_h_irr=ERT_TIMES_H[i], rhoa_true=rhoa_true, rhoa_obs=rhoa_obs,
                         eps=eps, quadrupoles=np.array([[int(v) for v in r[:4]] for r in tab['rows']]),
                         k=S.column(tab, 'k'))
    rho_nodes_all[i] = rho_nodes
    print(f"  survey {i:2d} t={ERT_TIMES_H[i]:7.2f} h: {rhoa_true.size} data, rhoa_true "
          f"[{rhoa_true.min():.1f}, {rhoa_true.max():.1f}] Ohm m, rho nodes [{rho_nodes.min():.1f}, {rho_nodes.max():.1f}]")
shutil.rmtree(tmpdir, ignore_errors=True)

# sanity: compare with the real data range (good rows of the first real survey)
_real = S.read_tldat_table(templates[surveys[0]])
_r = S.column(_real, 'rhoa')
_s = np.concatenate([v['rhoa_true'] for v in survey_out.values()])
print(f"rhoa range: real (survey {surveys[0]}) [{_r.min():.1f}, {_r.max():.1f}] | synthetic truth "
      f"[{_s.min():.1f}, {_s.max():.1f}] Ohm m")
if _s.max() > 3 * _r.max() or _s.min() < _r.min() / 3:
    warnings.warn("synthetic rhoa is far outside the real range: check Archie (a, m, n), porosity and the "
                  "initial condition before assimilating")

# =============================================================================
# 6. Virtual TDR (independent validation, never assimilated)
# =============================================================================
x0, y0 = X / 2, Y / 2
ztop = meshCATHY.points[:, 2].max()
SENSOR_DEPTHS = [0.15, 0.30, 0.50]
idx, coords = simu.find_nearest_node([[x0, y0, ztop - d] for d in SENSOR_DEPTHS])
tdr_nodes = [int(i) for i in np.atleast_1d(idx)]
phi_sensor = M.sensor_phi(CASE, SENSOR_DEPTHS, args.delta)
rng = np.random.RandomState(args.seed_noise + 7)
tdr = pd.DataFrame({'t_h': (sw_time - STABILIZE_PERIOD) / HOUR})
for d, n_, phi in zip(SENSOR_DEPTHS, tdr_nodes, phi_sensor):
    th = SW[:, n_] * phi
    tdr[f'theta_true_{int(d * 100)}'] = th
    tdr[f'theta_{int(d * 100)}'] = th + rng.normal(0.0, args.tdr_sigma, th.size)
tdr_file = os.path.join(OUT_ROOT, f'tdr_synth_{CASE}.csv')
tdr.to_csv(tdr_file, index=False, float_format='%.6f')
print("virtual TDR:", tdr_file, "| nodes", tdr_nodes, "| porosity", phi_sensor)

# =============================================================================
# 7. Export the truth
# =============================================================================
truth = dict(
    case=CASE, prj_name=prj_name, delta=cfg['delta'], atm_times=np.asarray(atm_times), sw_time=sw_time,
    stabilize_period=STABILIZE_PERIOD, sw=SW, psi=PSI, nodes=np.asarray(meshCATHY.points),
    porosity_nodes=poro_nodes, soil_table=soil_df, ic_nodes=np.asarray(ic_nodes),
    params=M.truth_parameters(CASE, args.delta), archie_truth=ARCHIE_TRUTH, archie_prior=ARCHIE_PRIOR,
    rFluid=rFLUID, porosity_mode=args.porosity_mode, surveys=survey_out, rho_nodes=rho_nodes_all,
    ert_frame=dict(signs=ERT_SIGNS, shift=ERT_SHIFT, origin=_origin), fwd_mesh=FWD_MESH, da_mesh=DA_MESH,
    noise=dict(pct=args.noise_pct, type='lognormal multiplicative', seed=args.seed_noise),
    tdr_nodes=tdr_nodes, tdr_depths=SENSOR_DEPTHS, tdr_file=tdr_file, tdr_sigma=args.tdr_sigma,
)
truth_file = os.path.join(OUT_ROOT, f'truth_{CASE}.pkl')
with open(truth_file, 'wb') as f:
    pickle.dump(truth, f)
meta = {k: (v if isinstance(v, (str, int, float, bool, type(None), dict)) else str(v)) for k, v in dict(
    case=CASE, delta=cfg['delta'], noise_pct=args.noise_pct, seed_noise=args.seed_noise, archie_truth=ARCHIE_TRUTH,
    archie_prior=ARCHIE_PRIOR, porosity_mode=args.porosity_mode, fwd_mesh=FWD_MESH, da_mesh=DA_MESH,
    same_mesh=(FWD_MESH == DA_MESH), n_surveys=len(surveys), argv=' '.join(sys.argv),
    data_dir=paths_syn['data_dir']).items()}
with open(os.path.join(OUT_ROOT, f'truth_meta_{CASE}.json'), 'w') as f:
    json.dump(meta, f, indent=2)

print(f"\ntruth: {truth_file}\nsynthetic ERT: {paths_syn['data_dir']}  (assimilate with ERT_DATA_ERR = {args.noise_pct:g})")
print("nb7 must use: same case/delta, ert_paths(..., synthetic=True), STABILIZE_PERIOD =", STABILIZE_PERIOD, "s")
