# =============================================================================
# Why does the ensemble mean differ from the truth, other than a biased prior?   (todo.md 5.1, last item)
# =============================================================================
#   python mean_offset_budget.py                       # NC, N = 16 32 64, collapsed ensembles N_eff = 12 2
#   python mean_offset_budget.py --case HC --n 32 64 --n-eff 12 4 2
#   python mean_offset_budget.py --csv ic_nonlinearity.csv
#
# `prior_bias_check.py` answers "is the prior centred on the truth?" (an infinite ensemble). Even with a
# centred prior the ensemble mean of a FINITE ensemble is not the truth. This script puts a number on the
# two reasons that do not depend on running CATHY, so they can be subtracted before a residual is blamed
# on the filter:
#
#   1. finite N        mean of N independent draws is off by sigma/sqrt(N). Expected |offset| =
#                      sqrt(2/pi)/sqrt(N) sigma, 95 % of the runs lie within 1.96/sqrt(N) sigma. With the
#                      resampling collapse of 5.2 the number that counts is the number of independent
#                      founders / N_eff, not N (N_eff = 2 -> 0.56 sigma expected, 1.39 sigma at 95 %).
#                      Corollary: the `--tol 0.2` of prior_bias_check is BELOW the noise floor of one
#                      N = 16 run (0.20 sigma expected, 0.49 at 95 %): a single run cannot confirm a
#                      0.2 sigma bias, only the offset of the sampled prior (200 000 draws) can.
#   2. nonlinearity    mean of sw over the prior is not sw at the mean parameters. For the initial
#                      condition: sw(psi) van Genuchten, psi ~ truncated normal(ic_true, sigma). Reported
#                      per CATHY layer as the offset (mean over the prior - sw at ic_true) in sw units and
#                      in prior standard deviations of sw. This is present for N -> infinity and for an
#                      exactly centred ic prior.
#
# Not covered here (needs a run): the PF collapse itself (pf_inbreeding.py / diag_pf_run.py, 5.2) and the
# nonlinearity of the ERT forward model (sigma of the resistivity data space, see diag_pf_run.py).
#
# Assumption (stated, not checked against CATHY): VGPSATCELL is the air-entry pressure head psi_s [m]
# of  Se = [1 + (|psi| / psi_s)^n]^(-(1 - 1/n)),   sw = thetar + (POROS - thetar) Se.
# The CLI needs model_config / make_scenarios (so pyCATHY importable); the functions do not.
# =============================================================================
import argparse
import sys

import numpy as np
import pandas as pd

from prior_bias_check import sample_prior


# ------------------------------------------------------------------ 1. finite N
def finite_n_offset(n):
    """Offset of the mean of n independent draws from the truth, in prior sigmas: (expected |offset|, 95 % bound)."""
    n = np.asarray(n, dtype=float)
    return np.sqrt(2.0 / np.pi) / np.sqrt(n), 1.96 / np.sqrt(n)


def finite_n_table(ns, n_effs=()):
    rows = [dict(kind='N (independent members)', n=int(n), expected_abs_sigma=finite_n_offset(n)[0],
                 bound95_sigma=finite_n_offset(n)[1]) for n in ns]
    rows += [dict(kind='N_eff / founders (collapsed)', n=int(n), expected_abs_sigma=finite_n_offset(n)[0],
                  bound95_sigma=finite_n_offset(n)[1]) for n in n_effs]
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ 2. nonlinearity (initial condition)
def vg_sw(psi, psat, n, thetar, poros):
    """Water content of a van Genuchten soil (psi <= 0 [m]; psat = air-entry pressure head [m] = VGPSATCELL)."""
    h = np.maximum(-np.asarray(psi, dtype=float), 0.0) / psat
    se = (1.0 + h ** n) ** (-(1.0 - 1.0 / n))
    return thetar + (poros - thetar) * se


def layer_depths(zr_thick):
    """(top, bottom, mid) of every CATHY layer from the layer thicknesses [m]."""
    bot = np.cumsum(zr_thick)
    top = bot - np.asarray(zr_thick, dtype=float)
    return top, bot, 0.5 * (top + bot)


def group_of_layers(zr_thick, ic_ranges):
    """ic group of every layer by the mid-depth rule (same as the library's map_param_to_layers)."""
    mid = layer_depths(zr_thick)[2]
    return [next(g for g, (a, b) in enumerate(ic_ranges) if a <= m < b or (g == len(ic_ranges) - 1 and m >= a))
            for m in mid]


def ic_nonlinearity(ic_true, ic_sigma, ic_bounds, layer_ic_group, layer_soil, n=200_000, seed=0):
    """Per layer: offset of the prior-mean water content from the water content at the true ic.

    ic_true / ic_sigma: one value per ic group; layer_ic_group: ic group of each layer; layer_soil: dict
    (VGPSATCELL, VGNCELL, VGRMCCELL, POROS) of each layer."""
    rng = np.random.default_rng(seed)
    lo, hi = ic_bounds
    draws = {g: sample_prior('normal', ic_true[g], ic_sigma[g], lo, hi, n, rng) for g in range(len(ic_true))}
    rows = []
    for L, (g, s) in enumerate(zip(layer_ic_group, layer_soil)):
        args = (s['VGPSATCELL'], s['VGNCELL'], s['VGRMCCELL'], s['POROS'])
        sw_true = float(vg_sw(ic_true[g], *args))
        sw = vg_sw(draws[g], *args)
        rows.append(dict(layer=L, ic_group=g, ic_true=ic_true[g], sw_true=sw_true, sw_prior_mean=float(sw.mean()),
                         sw_prior_sd=float(sw.std()), offset_sw=float(sw.mean()) - sw_true,
                         offset_sigma=(float(sw.mean()) - sw_true) / float(sw.std())))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ CLI
def case_inputs(case):
    """Everything ic_nonlinearity needs for a case, from model_config / make_scenarios (needs pyCATHY)."""
    import tomllib
    import make_scenarios as MS
    import model_config as mc
    cfg = mc.get_case(case)
    sc = tomllib.loads(MS.build(biased=False))['scenario'][f'{case}_full_upd']
    i = sc['per_name'].index('ic')
    n_ic = len(cfg['ic_true'])
    sig = sc['per_sigma'][i]
    sig = list(sig) if isinstance(sig, (list, tuple)) else [sig] * n_ic
    bnd = sc['per_bounds'][i]
    zr = list(mc.ZR_THICK)
    grp = group_of_layers(zr, cfg['ic_ranges'])
    mat_of = {L: m for m, ls in mc.MATERIAL_LAYERS_0.items() for L in ls}
    soil = [{**{k: cfg['soil'][1][mat_of[L]][k] for k in ('VGPSATCELL', 'VGNCELL', 'VGRMCCELL', 'POROS')}}
            for L in range(len(zr))]
    return dict(ic_true=[v[0] for v in cfg['ic_true']], ic_sigma=sig, ic_bounds=(bnd['min'], bnd['max']),
                layer_ic_group=grp, layer_soil=soil, mat=[mat_of[L] for L in range(len(zr))], mid=layer_depths(zr)[2])


def main(argv=None):
    ap = argparse.ArgumentParser(description='Finite-N and nonlinearity offsets of the ensemble mean (not prior bias)')
    ap.add_argument('--case', default='NC', choices=['NC', 'HC'])
    ap.add_argument('--n', type=int, nargs='+', default=[16, 32, 64], help='ensemble sizes')
    ap.add_argument('--n-eff', type=float, nargs='*', default=[12, 2],
                    help='collapsed ensembles: N_eff or number of founders (see pf_inbreeding.py / diag_pf_run.py)')
    ap.add_argument('--csv', default=None, help='write the nonlinearity table here')
    a = ap.parse_args(argv)
    pd.set_option('display.width', 200)

    ft = finite_n_table(a.n, a.n_eff)
    print('1. finite N: offset of the ensemble mean of a CENTRED prior, in prior sigmas\n')
    print(ft.to_string(index=False, formatters=dict(expected_abs_sigma='{:.2f}'.format, bound95_sigma='{:.2f}'.format)))
    print("\n   -> prior_bias_check.py --tol 0.2 is below the 95 % noise of any single run with N < 96.\n")

    inp = case_inputs(a.case)
    nl = ic_nonlinearity(inp['ic_true'], inp['ic_sigma'], inp['ic_bounds'], inp['layer_ic_group'], inp['layer_soil'])
    nl.insert(1, 'mid_depth_m', inp['mid'])
    nl.insert(2, 'material', inp['mat'])
    f3 = '{:.3f}'.format
    print(f"2. nonlinearity of sw(psi) over the ic prior ({a.case}, ic ~ truncated normal(ic_true, {inp['ic_sigma']}), "
          f"N -> inf)\n")
    print(nl.to_string(index=False, formatters=dict(mid_depth_m=f3, ic_true=f3, sw_true=f3, sw_prior_mean=f3,
                                                     sw_prior_sd=f3, offset_sw=lambda v: f'{v:+.4f}',
                                                     offset_sigma=lambda v: f'{v:+.2f}')))
    worst = nl.loc[nl.offset_sigma.abs().idxmax()]
    print(f"\n   largest: layer {int(worst.layer)} ({worst.material}), mean sw {worst.offset_sw:+.4f} "
          f"({worst.offset_sigma:+.2f} sd of sw) from sw at ic_true, with an exactly centred ic prior.")
    if a.csv:
        nl.to_csv(a.csv, index=False)
        print('written', a.csv)
    return 0


if __name__ == '__main__':
    sys.exit(main())
