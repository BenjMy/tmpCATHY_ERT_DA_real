# =============================================================================
# Is the "unbiased" prior really centred on the truth?   (todo.md, section 5)
# =============================================================================
#   python prior_bias_check.py                         # NC, every scenario, unbiased prior (--prior-bias none)
#   python prior_bias_check.py --case HC --scenario full_upd
#   python prior_bias_check.py --biased                # the biased prior, for comparison
#   python prior_bias_check.py --csv prior_bias.csv --tol 0.2
#
# `--prior-bias none` only means "the prior NOMINALS are the truth values". What the filter starts from is
# the SAMPLED ensemble, and that can be off-centre even when the nominal is exact:
#   * bounded normal priors are sampled as TRUNCATED normals (pyCATHY.DA.perturbate.sampling_dist_trunc):
#     truncating a tail moves the mean (Ks: sigma = nominal, lower bound = nominal/100 -> mean x ~1.29);
#   * lognormal priors are centred on the MEDIAN by default (mean = nominal * exp(sigma^2/2));
#   * the initial condition `ic` is ONE scalar for the 3 layers (make_scenarios.entries: -1.5), whereas the
#     truth is layered (model_config ic_true) -> not centred on the truth, whatever --prior-bias says.
# This script samples each prior like the library does and reports, per parameter and group, the offset of the
# ensemble mean from the truth in units of the prior sigma (on the scale the filter works on: log for Ks).
# It does not run CATHY: it needs model_config / make_scenarios (so pyCATHY must be importable, as for nb7).
# =============================================================================
import argparse
import json
import os
import sys
import tomllib

import numpy as np
import pandas as pd
from scipy.stats import truncnorm

SOIL_NAMES = ('SATCELL_VG', 'n_VG', 'thetar_VG', 'Ks')
ARCHIE_KEYS = {'a_Archie': 'a', 'm_Archie': 'm', 'n_Archie': 'n'}


def sample_prior(sampling, mean, sigma, lo=None, hi=None, n=200_000, rng=None, lognormal_center='median'):
    """Draw the prior like pyCATHY.DA.perturbate: truncated normal when bounded, lognormal centred on
    the median (default) or the mean, uniform. `sigma` is in natural-log units for the lognormal."""
    rng = np.random.default_rng(0) if rng is None else rng
    if sampling == 'lognormal':
        mu = np.log(mean) - (0.5 * sigma ** 2 if lognormal_center == 'mean' else 0.0)
        return rng.lognormal(mu, sigma, n)
    if sampling == 'uniform':
        h = np.sqrt(3.0) * sigma
        return rng.uniform(mean - h, mean + h, n)
    if sigma == 0:
        return np.full(n, float(np.clip(mean, -np.inf if lo is None else lo, np.inf if hi is None else hi)))
    a = -np.inf if lo is None else (lo - mean) / sigma
    b = np.inf if hi is None else (hi - mean) / sigma
    if a == -np.inf and b == np.inf:
        return rng.normal(mean, sigma, n)
    return truncnorm(a, b, loc=mean, scale=sigma).rvs(n, random_state=rng)


def _as_list(v, ng):
    return list(v) if isinstance(v, (list, tuple)) else [v] * ng


def _truth_archie(MS, case):
    """Archie truth = prior nominal unless nb6 was run with an offset: read it from truth_meta_<case>.json."""
    out = {k: MS.ARCHIE[k][0] for k in MS.ARCHIE}
    try:
        import project_paths as PP
        f = os.path.join(PP.get().synthetic_dir, case, f'truth_meta_{case}.json')
        with open(f) as fh:
            tr = json.load(fh).get('archie_truth', {})
        out.update({k: tr[s] for k, s in ARCHIE_KEYS.items() if s in tr})
    except Exception:
        pass
    return out


def prior_table(case='NC', scenario='full_upd', biased=False, n=200_000, tol=0.2, seed=0):
    """DataFrame: one row per parameter x group (and per ic layer) with the offset of the sampled prior
    mean from the truth. `bias_sigma` is on the filter's working scale; `flag` = |bias_sigma| > tol."""
    import make_scenarios as MS
    import model_config as mc
    sc = tomllib.loads(MS.build(biased=biased))['scenario'][f'{case}_{scenario}']
    rng = np.random.default_rng(seed)
    cfg = mc.get_case(case)
    arch = _truth_archie(MS, case)
    ng = len(MS.groups(case)[2])
    rows = []
    for i, name in enumerate(sc['per_name']):
        samp = sc['sampling_type'][i]
        bnd = sc['per_bounds'][i]
        lo, hi = (bnd.get('min'), bnd.get('max')) if isinstance(bnd, dict) else (None, None)
        if name == 'ic':
            for layer, tv in enumerate(cfg['ic_true']):
                x = sample_prior(samp, sc['per_mean'][i], sc['per_sigma'][i], lo, hi, n, rng)
                rows.append(_row(case, scenario, f'ic (layer {layer})', None, tv[0], sc['per_mean'][i],
                                 sc['per_sigma'][i], x, 'None', tol, samp))
            continue
        means = _as_list(sc['per_mean'][i], ng if name in SOIL_NAMES else 1)
        sigmas = _as_list(sc['per_sigma'][i], len(means))
        transf = _as_list(sc['transf_type'][i], len(means))
        truth = MS.group_values(case, name, False) if name in SOIL_NAMES else [arch[name]]
        for g, (m, s, t, tv) in enumerate(zip(means, sigmas, transf, truth)):
            x = sample_prior(samp, m, s, lo, hi, n, rng)
            rows.append(_row(case, scenario, name, g if name in SOIL_NAMES else None, tv, m, s, x, t, tol, samp))
    return pd.DataFrame(rows)


def _row(case, scenario, name, group, truth, nominal, sigma, x, transf, tol, sampling):
    if transf == 'log':                                    # the filter works on ln(x): offset of the log-mean
        lx = np.log(x)
        bias_sigma = (lx.mean() - np.log(truth)) / lx.std()
        scale = 'log'
    else:
        bias_sigma = (x.mean() - truth) / sigma
        scale = 'physical'
    geo = float(np.exp(np.log(x).mean())) if np.all(x > 0) else np.nan
    return dict(case=case, scenario=scenario, param=name, group=group, truth=truth, nominal=nominal, sigma=sigma,
                sampling=sampling, scale=scale, nominal_over_truth=nominal / truth if truth else np.nan,
                mean_over_truth=x.mean() / truth if truth else np.nan,
                median_over_truth=np.median(x) / truth if truth else np.nan,
                geomean_over_truth=geo / truth if truth else np.nan,
                bias_sigma=bias_sigma, flag=bool(abs(bias_sigma) > tol))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Is the 'unbiased' prior centred on the truth? (offset of the sampled prior mean, in sigma)")
    ap.add_argument('--case', default='NC', choices=['NC', 'HC'])
    ap.add_argument('--scenario', default=None, help="scenario tag (default: all of the case)")
    ap.add_argument('--biased', action='store_true', help="the biased prior (default: --prior-bias none)")
    ap.add_argument('--tol', type=float, default=0.2, help="flag when |mean - truth| > tol x sigma (default 0.2)")
    ap.add_argument('--n', type=int, default=200_000)
    ap.add_argument('--csv', default=None)
    a = ap.parse_args(argv)
    import make_scenarios as MS
    tags = [a.scenario] if a.scenario else [t for t, *_ in MS.MATRIX]
    df = pd.concat([prior_table(a.case, t, a.biased, a.n, a.tol) for t in tags], ignore_index=True)
    # the parameters are the same in every scenario that perturbs them: print each (param, group) once
    show = df.drop_duplicates(['param', 'group']) if not a.scenario else df
    pd.set_option('display.width', 200)
    g3, g4 = (lambda v: f'{v:.3f}'), (lambda v: f'{v:.4g}')
    print(f"prior: {'BIASED' if a.biased else 'unbiased (--prior-bias none)'}   case {a.case}   tol {a.tol} sigma\n")
    print(show[['param', 'group', 'truth', 'nominal', 'scale', 'mean_over_truth', 'median_over_truth',
                'geomean_over_truth', 'bias_sigma', 'flag']].to_string(
        index=False, na_rep='-', formatters=dict(
            group=lambda v: '-' if pd.isna(v) else str(int(v)), truth=g4, nominal=g4, mean_over_truth=g3,
            median_over_truth=g3, geomean_over_truth=g3, bias_sigma=g3)))
    bad = show[show.flag]
    print(f"\n{len(bad)} / {len(show)} prior entries are off-centre by more than {a.tol} sigma:",
          ', '.join(f"{r.param}" + ('' if pd.isna(r.group) else f"[{int(r.group)}]") for r in bad.itertuples()) or 'none')
    if a.csv:
        df.to_csv(a.csv, index=False)
        print('written', a.csv)
    return 1 if len(bad) else 0


if __name__ == '__main__':
    sys.exit(main())
