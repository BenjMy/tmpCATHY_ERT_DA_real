# Smoke test of scenario_nc_ERT.toml against the REAL pyCATHY scenario / perturbation code
# (stub simulation with the soil layout of nc_model: 1 zone x 15 layers, 0-based layers).
#   PYCATHY_SRC=/path/to/pycathy_wrapper python test_scenario_nc_ERT.py
import os
import sys
import tomllib

import numpy as np
import pandas as pd

sys.path.insert(0, os.environ.get("PYCATHY_SRC", "/home/claude/benjmy/pycathy_wrapper"))
from pyCATHY.DA import perturbate
from pyCATHY.DA.scenario import (validate_scenario, check_update_parameters,
                                 describe_ensemble)

HERE = os.path.dirname(os.path.abspath(__file__))
NSTR, NENS = 15, 32

with open(os.path.join(HERE, "scenario_nc_ERT.toml"), "rb") as f:
    SCEN = tomllib.load(f)["scenario"]

SOIL_INDEX = pd.MultiIndex.from_product([[1], range(NSTR)], names=["zone", "layer"])
ZONE = np.ones((10, 10))
VEG = np.ones((10, 10))


class Sim:                                    # what perturbate() reads from the DA object
    dem_parameters = {"base": 2.0}
    soil = {"PCWLT": [-150.0], "OMGC": [1.0], "PZ": [1.0]}
    soil_dict = soil

    def read_inputs(self, key, **kw):
        soil = pd.DataFrame({"PERMX": 1e-6, "POROS": 0.41}, index=SOIL_INDEX)
        feddes = pd.DataFrame({"ZROOT": [0.05]}, index=pd.Index([1]))
        return {"soil": (soil, feddes), "zone": (ZONE, {}), "root_map": (VEG, {})}[key]


def build(name):
    sc = SCEN[name]
    errors, warns = validate_scenario(sc, nlayers=NSTR, zone_raster=ZONE, veg_raster=VEG,
                                      raise_on_error=False)
    assert not errors, f"{name}: {errors}"
    np.random.seed(42)
    list_pert = perturbate.perturbate(Sim(), sc, NENS)
    out = {}
    for dp in list_pert:
        out = perturbate.perturbate_parm(
            out, parm=dp, type_parm=dp["type_parm"], mean=dp["mean"], sd=dp["sd"],
            sampling_type=dp["sampling_type"], ensemble_size=dp["ensemble_size"],
            per_type=dp["per_type"], savefig=False)
    check_update_parameters(list_pert, sc["listUpdateParm"])
    return sc, list_pert, out, warns


failures = []
for name, n_expected in (("nc_ERT_ic", 1), ("nc_ERT_ic_VGP_pert", 1 + 3 * NSTR),
                         ("nc_ERT_ic_VGP_upd", 1 + 3 * NSTR)):
    try:
        sc, lp, out, warns = build(name)
        assert len(lp) == n_expected, (len(lp), n_expected)
        for key, val in out.items():
            ens = np.asarray(val["ini_perturbation"], dtype=float)
            assert ens.shape == (NENS,), (key, ens.shape)
            assert np.all(np.isfinite(ens)), key
            lo, hi = val.get("clip_min"), val.get("clip_max")
            assert (lo is None or ens.min() >= lo - 1e-12) and (hi is None or ens.max() <= hi + 1e-12), key
            assert ens.std() > 0, f"{key}: no spread"
        if "VGP" in name:
            # prior centred on the nominal of each layer, n_VG > 1, thetar_VG >= 0
            for k in range(NSTR):
                assert abs(np.mean(out[f"n_VG{k}"]["ini_perturbation"]) - sc["per_nom"][2][k]) < 0.1
            assert min(np.min(out[f"n_VG{k}"]["ini_perturbation"]) for k in range(NSTR)) > 1.0
            assert min(np.min(out[f"thetar_VG{k}"]["ini_perturbation"]) for k in range(NSTR)) >= 0.0
            df = describe_ensemble(lp, sc["listUpdateParm"], soil_index=SOIL_INDEX,
                                   zone_raster=ZONE, veg_raster=VEG)
            assert not df.attrs["warnings"], df.attrs["warnings"]      # no override between entries
            assert df["target"].str.contains("ERROR").sum() == 0
            n_upd = int(df["updated"].sum())
            assert n_upd == (3 * NSTR if name.endswith("_upd") else 0), n_upd
        print(f"OK   {name}: {len(lp)} entries, scenario warnings: {warns}")
    except Exception as exc:                  # failure tracking: report all, then exit != 0
        failures.append((name, f"{type(exc).__name__}: {exc}"))
        print(f"FAIL {name}: {type(exc).__name__}: {exc}")

# the library must still reject what the scenario deliberately avoids
bad = dict(SCEN["nc_ERT_ic_VGP_pert"])
bad["per_name"] = ["ic", "VGP", "n_VG", "thetar_VG"]
errs, _ = validate_scenario(bad, nlayers=NSTR, zone_raster=ZONE, veg_raster=VEG, raise_on_error=False)
assert any("VGP" in e for e in errs), "validate_scenario no longer flags 'VGP' (library changed?)"

if failures:
    sys.exit(f"{len(failures)} failure(s): {failures}")
print("ALL OK")
