# =============================================================================
# Generator of scenario_nc_ERT.toml (perturbation scenarios of nb5, incl. VGP)
# =============================================================================
# Run from the notebook/ folder:   python make_scenario_nc_VGP.py
#
# Why a generator: the van Genuchten parameters (VGP) are perturbed with the
# 'layers' control, i.e. ONE list entry per soil layer (NSTR = 15 in nc_model)
# in per_nom / per_mean / per_sigma / per_type. Written by hand these lists are
# easy to get wrong (validate_scenario catches a wrong length, not a wrong value).
#
# Layout (pyCATHY.DA.scenario / perturbate, repository state of 2026-10-19):
#   - parameters built by perturbate(): SATCELL_VG (-> VGPSATCELL), n_VG
#     (-> VGNCELL), thetar_VG (-> VGRMCCELL). NOT usable: 'VGP' (its code is
#     commented out) and 'alpha_VG' (cannot be written to the soil files).
#   - keys <name><k>, k = 0..NSTR-1 is the 0-based layer from the top (the soil
#     table read back from disk numbers layers from 0).
#   - per_type is "None" (absolute sampling: the prior is N(per_mean, per_sigma)
#     truncated to per_bounds, all in the physical unit of the soil file). A
#     'multiplicative' type is avoided on purpose: with the 'layers' control the
#     prior bounds apply to the multiplicative FACTOR whereas the bounds used to
#     clip the analysed values apply to the PHYSICAL value.
#   - per_type must be a LIST (one entry per layer): a bare string is indexed
#     character by character by process_perturbation_block.
#
# FIXME (future fix, NOT implemented): mutual correlation between the VGP.
#   Every <name><k> is sampled independently of the others: alpha (SATCELL_VG),
#   n_VG and thetar_VG of one layer are uncorrelated in the prior, and so are
#   the layers of one material (layers 0-5, 6-8, 9-14 of nc_model.py are three
#   materials, each with one set of parameters). Physical VGP are correlated
#   (e.g. alpha vs n); neglecting it gives an unrealistic prior and lets the
#   EnKF update one parameter independently of the others. Botto, Belluco &
#   Camporese (2018, HESS 22, 4251-4266, https://hess.copernicus.org/articles/22/4251/2018/)
#   perturb alpha, n and thetar jointly: Johnson transforms of each parameter
#   to a normal variable (lognormal / log-ratio / hyperbolic = the transf_type
#   options already in pyCATHY), a covariance matrix between the transformed
#   parameters (Carsel & Parrish 1988), ensemble drawn through its Cholesky
#   factor, the transformed parameters in the EnKF state vector, then
#   back-transformed. To do in pyCATHY (perturbate.py / perturbate_parm_by_layers,
#   not in this script): (1) draw the VGP of one layer/material JOINTLY from a
#   multivariate normal (correlation matrix in the scenario, e.g. a new key
#   'per_corr'), (2) optionally share one draw across the layers of a material,
#   (3) apply transf_type before the draw (bounds in the transformed space),
#   (4) check the update keeps the correlation (parameter covariance in the
#   analysis, see enkf_analysis_inflation_with_parm_cov).
# =============================================================================
import os

NSTR = 15
# layer ranges (0-based, inclusive) of the three materials of nc_model.setup_model
MATERIALS = {"A": range(0, 6), "B": range(6, 9), "C": range(9, 15)}

# nominal soil-file values per material (copied from nc_model.setup_model;
# nb5 compares them with the soil table actually written, at run time)
NOMINAL = {
    "SATCELL_VG": {"A": 0.44, "B": 0.38, "C": 0.40},     # VGPSATCELL
    "n_VG":       {"A": 1.40, "B": 1.42, "C": 1.40},     # VGNCELL
    "thetar_VG":  {"A": 0.04, "B": 0.065, "C": 0.07},    # VGRMCCELL
}
# prior spread: relative (cv) or absolute (sd) std, and physical bounds
PRIOR = {
    "SATCELL_VG": dict(cv=0.25, bounds=(0.10, 2.00)),
    "n_VG":       dict(cv=0.05, bounds=(1.05, 2.50)),    # n > 1 required
    "thetar_VG":  dict(sd=0.02, bounds=(0.00, 0.15)),
}
VGP = list(NOMINAL)

# state / initial condition prior (global, no control)
IC = dict(nom=-1.5, mean=-1.5, sigma=0.5, bounds=(-5.0, -0.05))


def layer_values(name):
    out = [None] * NSTR
    for mat, layers in MATERIALS.items():
        for k in layers:
            out[k] = NOMINAL[name][mat]
    assert all(v is not None for v in out), f"{name}: layer without material"
    return out


def layer_sigma(name):
    p = PRIOR[name]
    return [round(v * p["cv"], 6) for v in layer_values(name)] if "cv" in p \
        else [p["sd"]] * NSTR


def toml_list(items):
    return "[ " + ", ".join(items) + " ]"


def fmt(v):
    return f'"{v}"' if isinstance(v, str) else repr(v)


def block(name, with_vgp, update_vgp, comment):
    """Text of one [scenario.<name>] table (lists aligned with per_name)."""
    names = ["ic"] + (VGP if with_vgp else [])
    lines = [f"# {comment}", f"[scenario.{name}]"]
    lines.append(f"per_name = {toml_list(fmt(n) for n in names)}")
    lines.append(f"pert_control_name = {toml_list(fmt('None' if n == 'ic' else 'layers') for n in names)}")

    def entry(n, key):
        if n == "ic":
            return {"nom": fmt(IC["nom"]), "mean": fmt(IC["mean"]), "sigma": fmt(IC["sigma"]),
                    "type": '"None"', "sampling": '"normal"', "transf": '"None"',
                    "bounds": "{ min = %s, max = %s }" % IC["bounds"]}[key]
        vals = layer_values(n)
        return {"nom": toml_list(map(repr, vals)), "mean": toml_list(map(repr, vals)),
                "sigma": toml_list(map(repr, layer_sigma(n))),
                "type": toml_list(['"None"'] * NSTR), "sampling": '"normal"', "transf": '"None"',
                "bounds": "{ min = %s, max = %s }" % PRIOR[n]["bounds"]}[key]

    for key, toml_key in (("type", "per_type"), ("nom", "per_nom"), ("mean", "per_mean"),
                          ("sigma", "per_sigma"), ("bounds", "per_bounds"),
                          ("sampling", "sampling_type"), ("transf", "transf_type")):
        lines.append(f"{toml_key} = [\n    " + ",\n    ".join(entry(n, key) for n in names) + "\n]")
    upd = ['"St. var."'] + ([f'"{n}{k}"' for n in VGP for k in range(NSTR)] if update_vgp else [])
    lines.append(f"listUpdateParm = {toml_list(upd)}")
    return "\n".join(lines) + "\n"


SCENARIOS = [
    ("nc_ERT_ic", False, False,
     "initial condition only, updated with the states (previous nb5 default)"),
    ("nc_ERT_ic_VGP_pert", True, False,
     "ic + VGP perturbed (prior spread of the model), only the states are updated"),
    ("nc_ERT_ic_VGP_upd", True, True,
     "ic + VGP perturbed AND updated by the EnKF (3 x 15 = 45 parameters: see the notes in nb5)"),
]

HEADER = """# scenario_nc_ERT.toml -- GENERATED by make_scenario_nc_VGP.py, edit the script, not this file.
# Read by nb5_DA_ERT_soil3D_nc.py:  tomllib.load(...)['scenario'][<name>]
# VGP = van Genuchten parameters per soil layer ('layers' control, 0-based layer keys):
#   SATCELL_VG -> VGPSATCELL, n_VG -> VGNCELL, thetar_VG -> VGRMCCELL.
# FIXME: the VGP are drawn independently (no mutual correlation), see make_scenario_nc_VGP.py.

"""

if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scenario_nc_ERT.toml")
    text = HEADER + "\n".join(block(*s) for s in SCENARIOS)
    with open(out, "w") as f:
        f.write(text)
    print("written:", out)
