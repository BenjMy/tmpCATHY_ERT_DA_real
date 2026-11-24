# =============================================================================
# Generator of the DA scenarios of the synthetic twin experiment (NC and HC)
# =============================================================================
# python make_scenarios.py [--out scenarios_twin.toml] [--kz-ratio] [--no-bias]
#
# Replaces make_scenario_nc_VGP.py. One TOML, [scenario.<CASE>_<name>], read by nb7.
#
# Controls (pyCATHY, patch pycathy_twin_fixes.patch):
#   NC : 'material'      -> 3 groups (A 0-5, B 6-8, C 9-14), written in the single zone
#   HC : 'zone_material' -> [zone, [layers]] groups: zone 1 A/B/C + zone 2 A (compacted band).
#        zone 2 below layer 5 is NOT controlled (nb2 leaves it at the library default,
#        model_config 'zone2_fill_below' = False); 'material' would overwrite the band with zone 1.
#
# Prior = the biased model of model_config.biased_overrides (Ks x5, SATCELL x1.3, n -0.08,
# thetar +0.01) unless --no-bias (prior centred on the truth: reference case).
# Sigmas are PHYSICAL (as nb5); bounds are physical. Ks: log transform, sigma = nominal.
#   TODO verify on a first run that the log-transformed Ks ensemble spans ~1 decade.
#
# Ks anisotropy ("c-lite", library): PERMZ = Ks x base kz/kx of each row, so HC keeps 0.017
# in the band. --kz-ratio adds scenarios with a Kz_ratio parameter: EXTENSION, NOT implemented in
# the library yet (the validator will reject the name) -- kept here so it can be switched on later.
#
# FIXME (as before): VGP drawn independently (Botto et al. 2018 joint draw not implemented).
# =============================================================================
import argparse
import os

import model_config as mc

LAYERS = mc.MATERIAL_LAYERS_0                       # {'A': range, 'B': ..., 'C': ...}
PRIOR = {                                           # sigma (cv: relative, sd: absolute), bounds
    "SATCELL_VG": dict(cv=0.25, bounds=(0.10, 2.00), col="VGPSATCELL"),
    "n_VG":       dict(cv=0.05, bounds=(1.05, 2.50), col="VGNCELL"),
    "thetar_VG":  dict(sd=0.02, bounds=(0.00, 0.15), col="VGRMCCELL"),
    "Ks":         dict(cv=1.00, bounds=None,         col="PERMX"),       # bounds from nominal
}
ARCHIE = dict(a_Archie=(0.3, 0.1, (0.05, 2.5)), m_Archie=(1.7, 0.15, (1.3, 3.0)),
              n_Archie=(1.7, 0.15, (1.3, 3.0)))     # (nominal, sigma, bounds) = nb5
IC_SIGMA, IC_BOUNDS = 0.5, (-5.0, -0.05)
# Ks prior. 'lognormal' (KS_PRIOR=lognormal): Ks = nominal * exp(N(0, KS_SIGMA_LOG^2)), sigma in natural-log units (0.5 = x1.65 at 1 sigma,
# x2.7 at 2 sigma, ~0.43 decade at 2 sigma): centred on the nominal in log space, always positive, no clipping.
# 'normal100' (default, unchanged): physical Gaussian with sd = 100 % of the nominal, clipped at nominal/100: ~16 % of the draws
# are negative and pile up at the lower bound (heavy-tailed prior, collapsed first updates in the N=32 run).
KS_PRIOR = os.environ.get("KS_PRIOR", "normal100")      # set KS_PRIOR=lognormal for the log-normal prior (untested)
KS_SIGMA_LOG = float(os.environ.get("KS_SIGMA_LOG", "0.5"))
VGP = ["SATCELL_VG", "n_VG", "thetar_VG"]


def groups(case):
    """(control, key, list of groups) ; a group is a list of layers (NC) or [zone, layers] (HC)."""
    mats = [list(LAYERS[m]) for m in "ABC"]
    if case == "NC":
        return "material", "material_layers", mats
    return "zone_material", "zone_material_layers", [[1, g] for g in mats] + [[2, mats[0]]]


def group_values(case, name, biased):
    """Nominal per group, from the (biased) soil table of model_config."""
    cfg = mc.get_case(case)
    src = mc.biased_overrides(case) if biased else {
        (z, m): v for z, ms in cfg["soil"].items() for m, v in ms.items()}
    col = PRIOR[name]["col"]
    out = []
    for g in groups(case)[2]:
        zone, layers = (1, g) if case == "NC" else (g[0], g[1])
        mat = next(m for m, ls in LAYERS.items() if layers[0] in ls)
        out.append(float(src[(zone, mat)][col]))
    return out


def fmt(v):
    if isinstance(v, str):
        return f'"{v}"'
    if isinstance(v, dict):
        return "{ min = %r, max = %r }" % (v["min"], v["max"])
    return repr(v)


def lst(items):
    return "[ " + ", ".join(fmt(i) for i in items) + " ]"


def entries(case, names, biased):
    """per-name dicts of the TOML keys."""
    control = groups(case)[0]
    ng = len(groups(case)[2])
    out = []
    for n in names:
        if n == "ic":
            nom = (-1.5,) * 1
            out.append(dict(name=n, control="None", nom=-1.5, mean=-1.5, sigma=IC_SIGMA, type="None",
                            bounds=dict(min=IC_BOUNDS[0], max=IC_BOUNDS[1]), transf="None"))
        elif n in ARCHIE:
            nom, sd, b = ARCHIE[n]
            out.append(dict(name=n, control="None", nom=nom, mean=nom, sigma=sd, type="None",
                            bounds=dict(min=b[0], max=b[1]), transf="None"))
        else:
            p = PRIOR[n]
            vals = group_values(case, n, biased)
            sig = [round(v * p["cv"], 10) for v in vals] if "cv" in p else [p["sd"]] * ng
            if n == "Ks" and KS_PRIOR == "lognormal":
                sig = [KS_SIGMA_LOG] * ng                    # natural-log units (np.random.lognormal(log(mean), sigma))
            if n == "Ks":
                b = dict(min=min(vals) / 100.0, max=max(vals) * 100.0)
                transf = ["log"] * ng
            else:
                b = dict(min=p["bounds"][0], max=p["bounds"][1])
                transf = ["None"] * ng
            out.append(dict(name=n, control=control, nom=vals, mean=list(vals), sigma=sig,
                            type=["None"] * ng, bounds=b, transf=transf))
    return out


def upd_keys(case, names):
    ng = len(groups(case)[2])
    keys = ["St. var."]
    for n in names:
        if n in ARCHIE:
            keys.append(n)
        elif n != "ic":
            keys += [f"{n}{k}" for k in range(ng)]
    return keys


def block(case, tag, names, updated, comment, biased=True):
    es = entries(case, names, biased)
    L = [f"# {comment}", f"[scenario.{case}_{tag}]", f"per_name = {lst(names)}",
         f"pert_control_name = {lst(e['control'] for e in es)}"]
    for key, tk in (("type", "per_type"), ("nom", "per_nom"), ("mean", "per_mean"), ("sigma", "per_sigma"),
                    ("bounds", "per_bounds")):
        L.append(f"{tk} = [\n    " + ",\n    ".join(fmt(e[key]) if not isinstance(e[key], list)
                                                   else lst(e[key]) for e in es) + "\n]")
    L.append("sampling_type = " + lst(["lognormal" if (n == "Ks" and KS_PRIOR == "lognormal") else "normal" for n in names]))
    L.append("transf_type = [\n    " + ",\n    ".join(
        lst(e["transf"]) if isinstance(e["transf"], list) else fmt(e["transf"]) for e in es) + "\n]")
    ckey = groups(case)[1]
    L.append(f"{ckey} = " + repr(groups(case)[2]).replace("range", ""))
    L.append(f"listUpdateParm = {lst(upd_keys(case, updated))}")
    return "\n".join(L) + "\n"


# (tag, perturbed names, updated names, comment)
MATRIX = [
    ("ic",          ["ic"],                                 [],            "ic only (states)"),
    ("ic_VGP",      ["ic"] + VGP,                           [],            "ic + VGP perturbed, states updated"),
    ("ic_VGP_upd",  ["ic"] + VGP,                           VGP,           "ic + VGP perturbed AND updated"),
    ("ic_Ks_upd",   ["ic", "Ks"],                           ["Ks"],        "ic + Ks perturbed AND updated"),
    ("ic_Archie_upd", ["ic"] + list(ARCHIE),                list(ARCHIE),  "ic + Archie (a, m, n) updated"),
    ("ic_VGP_Ks_upd", ["ic"] + VGP + ["Ks"],                VGP + ["Ks"],  "ic + VGP + Ks updated"),
    ("full_upd",    ["ic"] + VGP + ["Ks"] + list(ARCHIE),   VGP + ["Ks"] + list(ARCHIE),
     "ic + VGP + Ks + Archie, all updated"),
]
EXT = [("ic_Ks_Kzr_upd", ["ic", "Ks", "Kz_ratio"], ["Ks", "Kz_ratio"],
        "EXTENSION (library not ready): Ks + Kz_ratio updated; HC: prior ratio 0.1 vs true 0.017")]

HEADER = """# GENERATED by make_scenarios.py -- edit the script, not this file.
# [scenario.<CASE>_<tag>], CASE in NC / HC; read by nb7 (tomllib.load(...)['scenario'][name]).
# NC: control 'material' (material_layers); HC: 'zone_material' (zone_material_layers).
# Prior nominals: biased model (model_config.biased_overrides) unless generated with --no-bias.

"""


def build(biased=True, kz_ratio=False):
    text = HEADER
    for case in ("NC", "HC"):
        for tag, names, upd, com in MATRIX:
            text += block(case, tag, names, upd, com, biased) + "\n"
        if kz_ratio and case == "HC":
            tag, names, upd, com = EXT[0]
            text += f"# {com}\n# (not written: Kz_ratio unknown to the library; see synthetic_twin_plan.md, section 12)\n\n"
    return text


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "scenarios_twin.toml"))
    ap.add_argument("--no-bias", action="store_true")
    ap.add_argument("--kz-ratio", action="store_true")
    a = ap.parse_args()
    with open(a.out, "w") as f:
        f.write(build(not a.no_bias, a.kz_ratio))
    print("written:", a.out)
