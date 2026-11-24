"""Run naming shared by nb7 (writes) and run_matrix / nb8 (read): ONE place."""


def run_name(case, scenario, da_type='enkf', nens=100, seed=42, with_da=True, biased=True, tag=''):
    """SYN_<case>_<scenario tag>_<enkf|pf>_N<nens>_s<seed>[_noDA][_unbiased][_<tag>]."""
    sc = scenario[len(case) + 1:] if scenario.startswith(case + '_') else scenario
    return (f"SYN_{case}_{sc}_{da_type}_N{nens}_s{seed}" + ('' if with_da else '_noDA')
            + ('' if biased else '_unbiased') + (f'_{tag}' if tag else ''))
