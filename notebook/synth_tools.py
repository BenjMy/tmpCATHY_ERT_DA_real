# =============================================================================
# Pure helpers of the synthetic twin experiment (numpy/pandas only: unit-testable)
# =============================================================================
import numpy as np


# -----------------------------------------------------------------------------
# pygimli .tldat/.dat text format (same layout as nc_model.filter_tldat)
# -----------------------------------------------------------------------------
def read_tldat_table(path):
    """-> dict(head=[lines before the count line], n=int, header=str, cols=[names], rows=[list of tokens],
    tail=[lines after the data block]). Electrode block and column header are kept verbatim."""
    with open(path) as f:
        lines = f.read().split('\n')
    nel = int(lines[0].split()[0])
    i0 = 2 + nel
    ndat = int(lines[i0].split()[0])
    header = lines[i0 + 1]
    cols = header.lstrip('#').split()
    rows = [l.split() for l in lines[i0 + 2:i0 + 2 + ndat] if l.strip()]
    assert len(rows) == ndat, f"{path}: {len(rows)} rows read, {ndat} announced"
    return dict(head=lines[:i0], n=ndat, header=header, cols=cols, rows=rows,
                tail=lines[i0 + 2 + ndat:])


def column(tab, name):
    """Numeric column of a table returned by read_tldat_table."""
    j = tab['cols'].index(name)
    return np.array([float(r[j]) for r in tab['rows']])


def write_tldat_like(template, dst, rhoa, err_frac, r=None):
    """Write a .tldat identical to `template` (electrodes, quadrupoles, k, valid) but with the data
    columns replaced: rhoa, r = rhoa / k (unless r is given), err = err_frac (fraction, e.g. 0.05).

    Same column set and number format as the real files, so read_observations / filter_tldat /
    unreliable_quadrupoles need no change. Returns the table written."""
    tab = read_tldat_table(template)
    rhoa = np.asarray(rhoa, dtype=float)
    assert rhoa.size == tab['n'], f"rhoa has {rhoa.size} values, template has {tab['n']} rows"
    assert np.all(np.isfinite(rhoa)) and np.all(rhoa > 0), "rhoa must be finite and positive"
    cols = tab['cols']
    assert {'rhoa', 'r', 'err', 'k'} <= set(cols), f"unexpected columns {cols}"
    ik, ir, ia, ie = (cols.index(c) for c in ('k', 'r', 'rhoa', 'err'))
    k = column(tab, 'k')
    r = rhoa / k if r is None else np.asarray(r, dtype=float)
    out_rows = []
    for i, v in enumerate(tab['rows']):
        v = list(v)
        v[ia] = f"{rhoa[i]:.14e}"
        v[ir] = f"{r[i]:.14e}"
        v[ie] = f"{err_frac:.14e}"
        out_rows.append('\t'.join(v))
    with open(dst, 'w') as f:
        f.write('\n'.join(tab['head'] + [str(tab['n']), tab['header']] + out_rows + tab['tail']))
    return tab


# -----------------------------------------------------------------------------
# noise
# -----------------------------------------------------------------------------
def add_lognormal_noise(rhoa, rel_sigma, seed):
    """rhoa_obs = rhoa * exp(rel_sigma * eps), eps ~ N(0, 1), reproducible (own RandomState).

    Multiplicative log-normal noise: positive by construction (the DA assimilates log10(rhoa)), relative
    std = rel_sigma (0.05 = 5 %), and sigma_log10 = rel_sigma / ln(10) ~ 0.0217. Returns (rhoa_obs, eps)."""
    rhoa = np.asarray(rhoa, dtype=float)
    eps = np.random.RandomState(int(seed)).standard_normal(rhoa.shape)
    return rhoa * np.exp(rel_sigma * eps), eps


def archie_rho(sw, porosity, rfluid, a, m, n):
    """rho = rFluid * a * phi^-m * sw^-n, element-wise on nodal arrays (physically correct porosity map)."""
    sw = np.asarray(sw, dtype=float)
    phi = np.asarray(porosity, dtype=float)
    assert sw.shape == phi.shape, f"sw {sw.shape} and porosity {phi.shape} differ"
    assert np.all(sw > 0) and np.all(phi > 0), "sw and porosity must be > 0 for Archie"
    return rfluid * a * phi ** (-m) * sw ** (-n)


def time_index(times, t, tol=1e-3):
    """Index of the output time equal to t [s] (raise if absent)."""
    times = np.asarray(times, dtype=float)
    i = int(np.argmin(np.abs(times - t)))
    assert abs(times[i] - t) <= tol, f"time {t} s not in the output times (closest {times[i]} s)"
    return i


def interp_in_time(arr, t_src, t_tgt):
    """Linear interpolation in time of arr[time, node] from the CATHY output times t_src to t_tgt [s].

    CATHY prints at the first time STEP at or after each requested time (a continuous run does not stop
    at the atmbc times: 3627.9 s printed for a request at 3600 s), whereas the DA windows end exactly at
    the atmbc times. Returns (arr_tgt, max gap [s] between the two output times around each target)."""
    t_src = np.asarray(t_src, dtype=float)
    t_tgt = np.asarray(t_tgt, dtype=float)
    assert np.all(np.diff(t_src) > 0), "output times must be strictly increasing"
    assert t_tgt.min() >= t_src[0] - 1e-6 and t_tgt.max() <= t_src[-1] + 1e-6, "target times outside the run"
    j = np.clip(np.searchsorted(t_src, t_tgt, side='right') - 1, 0, len(t_src) - 2)
    w = ((t_tgt - t_src[j]) / (t_src[j + 1] - t_src[j]))[:, None]
    arr = np.asarray(arr, dtype=float)
    return (1 - w) * arr[j] + w * arr[j + 1], float(np.max(t_src[j + 1] - t_src[j]))
