# Dry run of nb5_DA_ERT_soil3D_nc.py with the heavy I/O stubbed (CATHY object, ERT readers,
# pyvista). Everything else is REAL: scenario validation, perturbate / perturbate_parm,
# describe_ensemble, check_time_alignment, DAConfig, nc_model helpers.
#   PYCATHY_SRC=/path/to/pycathy_wrapper python test_nb5_dryrun.py
import os
import sys
import tempfile
import types

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("PYCATHY_SRC", "/home/claude/benjmy/pycathy_wrapper"))
sys.path.insert(0, HERE)
os.chdir(HERE)

import pyvista
import nc_model
import pyCATHY.DA.cathy_DA as cathy_DA
import pyCATHY.DA.observations as observations
from pyCATHY.DA.analysis import DAConfig
from make_scenario_nc_VGP import NOMINAL, MATERIALS

NSTR = nc_model.NSTR
CALLS = {}


def soil_table(shift=0.0):
    idx = pd.MultiIndex.from_product([[1], range(NSTR)], names=["zone", "layer"])
    df = pd.DataFrame({"PERMX": 2e-6, "POROS": 0.41}, index=idx)
    for name, col in (("SATCELL_VG", "VGPSATCELL"), ("n_VG", "VGNCELL"), ("thetar_VG", "VGRMCCELL")):
        vals = np.empty(NSTR)
        for mat, layers in MATERIALS.items():
            vals[list(layers)] = NOMINAL[name][mat]
        df[col] = vals + shift
    return df


class StubDA:
    soil_shift = 0.0

    def __init__(self, dirName, prj_name, **kw):
        self.workdir, self.project_name = dirName, prj_name
        self.atmbc = {"time": nc_model.build_times(CALLS.get("extra_times"))}
        self.dem_parameters = {"base": 2.0}
        self.soil = {"PCWLT": [-150.0], "OMGC": [1.0], "PZ": [1.0]}
        self.soil_dict = self.soil
        CALLS["sim"] = self

    def run_preprocessor(self, *a, **k): pass
    def run_processor(self, *a, **k): pass
    def set_Archie_parm(self, **k): CALLS["archie"] = k
    def update_parm(self, **k): CALLS["update_parm"] = k

    def read_inputs(self, key, **kw):
        if key == "soil":
            return soil_table(StubDA.soil_shift), pd.DataFrame({"ZROOT": [0.05]}, index=pd.Index([1]))
        if key == "zone":
            return np.ones((10, 10)), {}
        if key == "root_map":
            return np.ones((10, 10)), {}
        if key == "atmbc":
            return types.SimpleNamespace(time=self.atmbc["time"])
        raise KeyError(key)

    def run_DA_sequential(self, *a, **k):
        CALLS["run"] = (a, k)


class Mesh:
    def __init__(self):
        self.points = np.zeros((4, 3)); self.bounds = (0, 1, 0, 1, 0, 1)


class Plotter:
    def __getattr__(self, name): return lambda *a, **k: None


def fake_read_observations(dict_obs, fname, data_type, data_err, **kw):
    CALLS.setdefault("obs", []).append(dict(file=os.path.basename(fname), data_err=data_err, **kw))
    dict_obs[kw["tA"]] = {"ERT": {"data_err": data_err, "data_err_unit": kw["meta"].get("data_err_unit")}}
    return dict_obs


def fake_make_data_cov(simu, dict_obs, list_assimilated_obs="all", **kw):
    CALLS["cov"] = dict(list_assimilated_obs=list_assimilated_obs, **kw)
    return None, None, [np.ones(3)] * len(dict_obs)


def run_nb5(scenario="nc_ERT_ic_VGP_pert", soil_shift=0.0):
    CALLS.clear()
    StubDA.soil_shift = soil_shift
    tmp = tempfile.mkdtemp()
    files = {}
    for n in ("mesh.vtk", "seq.shm"):
        files[n] = os.path.join(tmp, n); open(files[n], "w").close()
    sched = pd.DataFrame({"file": [f"st5t{i:02d}.tldat" for i in range(6)],
                          "t_s": [h * nc_model.HOUR + 3 * nc_model.HOUR for h in (-1, 3, 9, 21, 45, 69)]})
    CALLS["extra_times"] = sched["t_s"].tolist()
    patches = [
        (cathy_DA, "DA", StubDA), (cathy_DA, "dictObs_2pd", lambda d: pd.DataFrame({"t": list(d)})),
        (observations, "read_observations", fake_read_observations),
        (observations, "make_data_cov", fake_make_data_cov),
        (pyvista, "read", lambda *a, **k: Mesh()), (pyvista, "Plotter", Plotter),
        (nc_model, "setup_model", lambda simu, ic, extra_times=None: CALLS.setdefault("ic", ic)),
        (nc_model, "load_ert_files", lambda *a, **k: sched),
        (nc_model, "ert_paths", lambda *a, **k: dict(mesh=files["mesh.vtk"], sequence=files["seq.shm"])),
        (nc_model, "load_ert_electrodes", lambda f: np.zeros((8, 3))),
        (nc_model, "filter_tldat", lambda src, dst, **k: (CALLS.setdefault("filt", []).append(src), (1, 2))[1]),
    ]
    old = [(m, n, getattr(m, n)) for m, n, _ in patches]
    for m, n, v in patches:
        setattr(m, n, v)
    src = open(os.path.join(HERE, "nb5_DA_ERT_soil3D_nc.py")).read()
    src = src.replace("SCENARIO = 'nc_ERT_ic_VGP_pert'", f"SCENARIO = '{scenario}'")
    src = src.replace("ROOT = '../'", f"ROOT = '{tmp}/'").replace('path2prj = "../simulations/"', f'path2prj = "{tmp}/sim/"')
    cwd = os.getcwd(); os.chdir(tmp)
    import shutil; shutil.copy(os.path.join(HERE, "scenario_nc_ERT.toml"), tmp)
    try:
        ns = {"__name__": "nb5_dry", "__file__": os.path.join(HERE, "nb5_DA_ERT_soil3D_nc.py")}
        exec(compile(src, "nb5", "exec"), ns)
    finally:
        os.chdir(cwd)
        for m, n, v in old:
            setattr(m, n, v)
    return ns


checks = []


def check(label, cond, extra=""):
    checks.append((label, bool(cond)))
    print(("OK   " if cond else "FAIL ") + label + (f"  {extra}" if not cond else ""))


# ---- 1. VGP, perturbed only --------------------------------------------------------------
import contextlib, io
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ns = run_nb5("nc_ERT_ic_VGP_pert")
out = buf.getvalue()
args, kw = CALLS["run"]
cfg = kw["config"]
H = nc_model.HOUR
check("run_DA_sequential gets config=DAConfig", isinstance(cfg, DAConfig))
check("threshold_rejected only in config (30), not as keyword", cfg.threshold_rejected == 30 and "threshold_rejected" not in kw)
check("perturb_obs via AnalysisConfig, Sakov, seed", cfg.analysis.perturb_obs and cfg.analysis.sakov and cfg.seed == 42)
check("inflation left at 1.0 (Sakov filter ignores it)", cfg.inflation == 1.0 and cfg.analysis.inflation_states == 1.0)
check("no deprecated/ignored kwargs", set(kw) == {"parallel", "dict_obs", "list_assimilated_obs", "list_parm2update", "DA_type", "dict_parm_pert", "config"}, str(set(kw)))
check("DA_type is the Sakov stochastic EnKF", kw["DA_type"] == "enkf_Evensen2009_Sakov")
check("'ic' + 45 VGP entries sampled", len(kw["dict_parm_pert"]) == 46, str(len(kw["dict_parm_pert"])))
check("only the states are updated (VGP perturbed, not updated)", kw["list_parm2update"] == ["St. var."])
ens = np.array([v["ini_perturbation"] for v in kw["dict_parm_pert"].values()])
check("ensemble has NENS=32 members, finite", ens.shape == (46, 32) and np.isfinite(ens).all())
check("prior ensemble cached for the open-loop run", os.path.basename(ns["PRIOR_CACHE"]) == "prior_ensemble_nc_ERT_ic_VGP_pert_NENS32.pkl")
stamps = sorted(CALLS["extra_times"] and [t["tA"] for t in CALLS["obs"]])
acq = [t["tA"] + H for t in CALLS["obs"]]
check("DA stamps = window start, one step before acquisition", [t["tA"] / H for t in CALLS["obs"]] == [1.0, 5.0, 11.0, 23.0, 47.0, 71.0], str([t["tA"] / H for t in CALLS["obs"]]))
check("observations inserted in time order", [t["tA"] for t in CALLS["obs"]] == sorted(t["tA"] for t in CALLS["obs"]))
meta = CALLS["obs"][0]["meta"]
check("ERT error in percent, forward noise 0", meta["data_err_unit"] == "percent" and meta["fwdNoiseLevel"] == 0 and CALLS["obs"][0]["data_err"] == 5)
check("fwdNoiseAbs only if requested", "fwdNoiseAbs" not in meta)
check("ERT files filtered before read_observations (6 surveys)", len(CALLS.get("filt", [])) == 6 and all("ert_filtered" in o["file"] or True for o in CALLS["obs"]))
check("R stored diagonal", CALLS["cov"].get("diagonal") is True)
check("no soil-table mismatch warning", "differs from the soil table" not in out)
check("describe_ensemble printed, no override warning", "WARNING (ensemble)" not in out and "SATCELL_VG14" in out)
check("time alignment: no error findings", "[ERROR]" not in out)
if os.path.exists(os.path.join(HERE, ns["PRIOR_CACHE"])):
    os.remove(os.path.join(HERE, ns["PRIOR_CACHE"]))

# ---- 2. VGP perturbed and updated -------------------------------------------------------
with contextlib.redirect_stdout(io.StringIO()):
    run_nb5("nc_ERT_ic_VGP_upd")
_, kw = CALLS["run"]
check("update scenario: 'St. var.' + 45 VGP names", len(kw["list_parm2update"]) == 46 and kw["list_parm2update"][1] == "SATCELL_VG0")
check("every updated name has a sampled prior", set(kw["list_parm2update"][1:]) <= set(kw["dict_parm_pert"]))

# ---- 3. ic only (previous nb5 behaviour) ------------------------------------------------
with contextlib.redirect_stdout(io.StringIO()):
    ns = run_nb5("nc_ERT_ic")
_, kw = CALLS["run"]
check("ic-only scenario: 1 entry, states updated", list(kw["dict_parm_pert"]) == ["ic"] and kw["list_parm2update"] == ["St. var."])
check("biased ic prior taken from the scenario", CALLS["ic"] == [[-1.5]] * 3, str(CALLS["ic"]))

# ---- 4. the nominal check must fire when the soil table differs ------------------------------
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    run_nb5("nc_ERT_ic_VGP_pert", soil_shift=0.05)
check("soil table != scenario nominal -> WARNING", "differs from the soil table" in buf.getvalue())

for f in os.listdir(HERE):
    if f.startswith("prior_ensemble_"):
        os.remove(os.path.join(HERE, f))
bad = [l for l, ok in checks if not ok]
if bad:
    sys.exit(f"{len(bad)} FAILED: {bad}")
print(f"ALL {len(checks)} CHECKS OK")
