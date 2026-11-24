# =============================================================================
# Per-case model definition (NC = uncompacted, HC = compacted): nb1 / nb2 as ONE code path
# =============================================================================
# Generalises nc_model.setup_model. nc_model.py stays untouched (nb4/nb5 keep working); every
# helper that does not depend on the case (build_times, snap_to_times, obs_window_start,
# filter_tldat, unreliable_quadrupoles, ...) is re-exported from it.
#
# Used by: nb6 (truth run, synthetic ERT), nb7 (DA), nb8 (plots).
# Truth and DA MUST call setup_model with the same case / delta / extra_times.
#
#   setup_model(simu, 'NC')                                  -> nb1 (truth)
#   setup_model(simu, 'HC', delta=0.25)                      -> nb2 on a coarser mesh
#   setup_model(simu, 'NC', ic_values=[[-1.5]]*3,
#               soil_overrides=biased_overrides('NC'))       -> biased prior for the DA
import copy
import warnings
import numpy as np
import pandas as pd

import nc_model as _nc
from nc_model import (HOUR, DAY, T_END, X, Y, MAXDEPTH, ELEVATION, ZR_THICK, NSTR, ELSTOR,  # noqa: F401
                      IC_RANGES, build_times, snap_to_times, obs_window_start, filter_tldat,
                      unreliable_quadrupoles, set_ert_boundary_markers, load_ert_electrodes,
                      load_ert_files, load_sensors_from_shm, ERT_COMP_LEVELS, scenario_to_material,
                      MATERIAL_PARAMETERS)

# 0-based layers (from the top) of each material. Disjoint: the label slices of nb1/nb2
# ([1:6], [7:10], [10:]) overlap at label 10 and the later one wins, i.e. A = 1-6, B = 7-9, C = 10-15.
MATERIAL_LAYERS_0 = {'A': list(range(0, 6)), 'B': list(range(6, 9)), 'C': list(range(9, NSTR))}
SOIL_COLS = ['PERMX', 'PERMY', 'PERMZ', 'ELSTOR', 'POROS', 'VGPSATCELL', 'VGNCELL', 'VGRMCCELL']
COMP_LEVEL = {'NC': 'no_comp', 'HC': 'heavy_comp'}


def _soil(kx, ky, kz, phi, sat, n, thr):
    return dict(PERMX=kx, PERMY=ky, PERMZ=kz, ELSTOR=ELSTOR, POROS=phi,
                VGPSATCELL=sat, VGNCELL=n, VGRMCCELL=thr)


# -----------------------------------------------------------------------------
# Cases: values copied from nb1 (NC) and nb2 (HC)
# -----------------------------------------------------------------------------
CASES = {
    'NC': dict(
        delta=0.5, rain_total=0.05, irr_start_hour=3, et=-1e-7, t_end=T_END,
        ic_ranges=IC_RANGES, ic_true=[[-2.4], [-1.3], [-1.15]],
        band=None,                                    # no zonation (1 zone)
        soil={1: {'A': _soil(2.1e-6, 2.1e-6, 2.1e-6, 0.415, 0.44, 1.40, 0.04),
                  'B': _soil(2.3e-6, 2.2e-6, 2.2e-6, 0.41, 0.38, 1.42, 0.065),
                  'C': _soil(1.5e-6, 1.5e-6, 1.5e-6, 0.41, 0.40, 1.40, 0.07)}},
        sequence='seq3d_from_st5t00tldat.shm', sequence_verified=True,
    ),
    'HC': dict(
        delta=0.1, rain_total=0.04, irr_start_hour=2, et=-1e-7, t_end=T_END,
        ic_ranges=IC_RANGES, ic_true=[[-1.4], [-1.13], [-1.13]],
        band=dict(center=X / 2, width=0.4, zone=2),   # compacted column (nb2: zones[:, 23:27] = 2)
        soil={1: {'A': _soil(2.5e-6, 2.5e-6, 2.8e-6, 0.41, 0.30, 1.40, 0.07),
                  'B': _soil(3.3e-6, 3.3e-6, 3.0e-6, 0.41, 0.50, 1.50, 0.06),
                  'C': _soil(1.0e-6, 1.0e-6, 1.0e-6, 0.41, 0.40, 1.40, 0.07)},
              # zone 2: only material A (layers 1-6 = 0-20 cm) is set in nb2
              2: {'A': _soil(3e-7, 3e-7, 0.5e-8, 0.35, 0.55, 1.28, 0.045)}},
        # nb2 never sets zone 2 below layer 6: those rows keep the init_soil_SPP_map_df defaults.
        # False = faithful to nb2 (and warns); True = copy zone 1 (materials B, C) into zone 2.
        zone2_fill_below=False,
        # Checked against st1t00.tldat (HC campaign): the 72 electrodes are IDENTICAL to NC, so the NC .shm is
        # valid as an electrode file (and for meshForward). The QUADRUPOLES differ (HC 4911 vs NC 4790, 4771
        # in common, 2203 good rows vs 2837): they always come from the HC data files themselves (templates).
        sequence='seq3d_from_st5t00tldat.shm', sequence_verified=True,
    ),
}


def get_case(case, delta=None):
    """Deep copy of CASES[case] with derived geometry. delta overrides the horizontal cell size [m]."""
    if case not in CASES:
        raise KeyError(f"unknown case {case!r}; available: {list(CASES)}")
    cfg = copy.deepcopy(CASES[case])
    cfg['name'] = case
    if delta is not None:
        cfg['delta'] = float(delta)
    d = cfg['delta']
    cfg['ncell'] = int(np.round(X / d))
    assert np.isclose(cfg['ncell'] * d, X), f"delta={d} does not divide the {X} m plot"
    assert cfg['t_end'] == _nc.T_END, "per-case t_end is not supported by build_times (global T_END)"
    cfg['irr'] = (cfg['rain_total'] / 2) / HOUR       # [m/s] for 1 h (as nb1/nb2)
    return cfg


# -----------------------------------------------------------------------------
# Pure helpers (no pyCATHY: unit-testable)
# -----------------------------------------------------------------------------
def zones_array(cfg):
    """Zone raster (ncell x ncell) of 1s, with the compacted band = band['zone']. None if no band."""
    if cfg['band'] is None:
        return None
    n, d, b = cfg['ncell'], cfg['delta'], cfg['band']
    c0 = int(np.round((b['center'] - b['width'] / 2) / d))
    c1 = int(np.round((b['center'] + b['width'] / 2) / d))
    if c1 <= c0:
        raise ValueError(f"band of {b['width']} m is thinner than one cell ({d} m)")
    if abs((c1 - c0) * d - b['width']) > 1e-9:
        warnings.warn(f"compacted band is {(c1 - c0) * d:.3f} m wide on a {d} m mesh (asked {b['width']} m)")
    z = np.ones((n, n))
    z[:, c0:c1] = b['zone']
    return z


def irrigation_window(cfg):
    """(y_min, y_max, x_min, x_max) cell slice of the 2 x 2 m central irrigated square (as nb1/nb2)."""
    n, d = cfg['ncell'], cfg['delta']
    hs = int(2.0 / d / 2)
    c = n // 2
    return c - hs, c + hs, c - hs, c + hs


def atm_steps(cfg, times):
    """List of (n+1 x n+1) net-flux rasters, one per time: ET everywhere, irrigation pulse in the square."""
    n = cfg['ncell']
    y0, y1, x0, x1 = irrigation_window(cfg)
    out = []
    for t in times:
        a = np.ones((n + 1, n + 1)) * cfg['et']
        if t == cfg['irr_start_hour'] * HOUR:
            a[y0:y1, x0:x1] = cfg['irr']
        out.append(a)
    return out


def fill_soil_df(df, cfg, overrides=None):
    """Write the soil table of a case into df (index = (zone, 1-based layer), SOIL_COLS columns).

    overrides: {'A': {col: value}} (all zones) or {(zone, 'A'): {col: value}}; absolute values.
    Returns df (modified in place)."""
    idx = pd.IndexSlice
    soil = copy.deepcopy(cfg['soil'])
    if cfg.get('zone2_fill_below') and cfg['band'] is not None:
        zb = cfg['band']['zone']
        for mat, vals in cfg['soil'][1].items():
            soil[zb].setdefault(mat, dict(vals))
    elif cfg['band'] is not None:
        zb = cfg['band']['zone']
        missing = sorted(set(MATERIAL_LAYERS_0) - set(soil[zb]))
        if missing:
            warnings.warn(f"zone {zb}: materials {missing} are not set (rows keep the init_soil_SPP_map_df "
                          "defaults), as in nb2. Set CASES['HC']['zone2_fill_below']=True to copy zone 1.")
    overrides = overrides or {}
    for zone, mats in soil.items():
        for mat, vals in mats.items():
            vals = dict(vals)
            for key in (mat, (zone, mat)):
                vals.update(overrides.get(key, {}))
            labels = [l + 1 for l in MATERIAL_LAYERS_0[mat]]
            df.loc[idx[zone, labels], SOIL_COLS] = [vals[c] for c in SOIL_COLS]
    return df


def biased_overrides(case, satcell_factor=1.3, n_shift=-0.08, thetar_shift=0.01, ks_factor=5.0,
                     n_min=1.05):
    """soil_overrides for a BIASED prior model: every (zone, material) of the truth is moved by the
    given factors (VGPSATCELL x, VGNCELL +, VGRMCCELL +, PERMX/Y/Z x). Pass to setup_model(...)."""
    cfg = get_case(case)
    out = {}
    for zone, mats in cfg['soil'].items():
        for mat, v in mats.items():
            out[(zone, mat)] = dict(
                VGPSATCELL=v['VGPSATCELL'] * satcell_factor,
                VGNCELL=max(v['VGNCELL'] + n_shift, n_min),
                VGRMCCELL=max(v['VGRMCCELL'] + thetar_shift, 0.0),
                PERMX=v['PERMX'] * ks_factor, PERMY=v['PERMY'] * ks_factor, PERMZ=v['PERMZ'] * ks_factor)
    return out


def truth_parameters(case, delta=None):
    """Everything nb8 needs to compare parameters with the truth (stored in truth_<case>.pkl)."""
    cfg = get_case(case, delta)
    return dict(case=case, delta=cfg['delta'], ic_ranges=cfg['ic_ranges'], ic_true=cfg['ic_true'],
                soil=cfg['soil'], material_layers=MATERIAL_LAYERS_0, zr_thick=list(ZR_THICK),
                band=cfg['band'], rain_total=cfg['rain_total'], irr_start_hour=cfg['irr_start_hour'],
                zone2_fill_below=cfg.get('zone2_fill_below'))


def sensor_phi(case, depths=(0.15, 0.30, 0.50), delta=None):
    """Porosity at the sensor depths at the plot centre (zone of the centre cell; rows of that zone that
    are not set fall back to zone 1, which is what nb2 does for the 30/50 cm sensors)."""
    cfg = get_case(case, delta)
    edges = np.cumsum(np.r_[0.0, ZR_THICK])
    z = zones_array(cfg)
    zone = 1 if z is None else int(z[cfg['ncell'] // 2, cfg['ncell'] // 2])
    out = []
    for d in depths:
        layer = int(np.clip(np.searchsorted(edges, d, side='right') - 1, 0, NSTR - 1))
        mat = next(m for m, ls in MATERIAL_LAYERS_0.items() if layer in ls)
        src = cfg['soil'].get(zone, {}).get(mat) or cfg['soil'][1][mat]
        out.append(src['POROS'])
    return out


def ert_paths(root='../', case='NC', synthetic=False):
    """Like nc_model.ert_paths, with the sequence file of the case; synthetic=True points data_dir to
    <root>data/3d/synthetic/<case>/procdata (written by nb6)."""
    import os
    cfg = CASES[case]
    p = _nc.ert_paths(root, COMP_LEVEL[case])
    p['sequence'] = os.path.join(root, 'data', '3d', cfg['sequence'])
    if not cfg['sequence_verified']:
        warnings.warn(f"{case}: electrode sequence {cfg['sequence']} is NOT verified for this campaign "
                      "(check CASES[...]['sequence'])")
    if synthetic:
        p['data_dir'] = os.path.join(root, 'data', '3d', 'synthetic', case, 'procdata')
    return p


def load_case_ert_files(root='../', case='NC', times_h=None, stabilize_period=None, synthetic=False):
    """Like nc_model.load_ert_files for a case; synthetic=True reads the files written by nb6.
    times_h: {survey index: hours since irrigation}. Returns a DataFrame (file, t_s, survey)."""
    import os
    cfg = CASES[case]
    stabilize_period = cfg['irr_start_hour'] * HOUR if stabilize_period is None else stabilize_period
    _, prefix, _ = ERT_COMP_LEVELS[COMP_LEVEL[case]]
    data_dir = ert_paths(root, case, synthetic)['data_dir']
    rows = []
    for i, th in sorted(times_h.items()):
        fname = os.path.join(data_dir, f"{prefix}t{int(i):02d}.tldat")
        if not os.path.isfile(fname):
            print(f"  [skip] ERT file missing: {fname}")
            continue
        rows.append(dict(file=fname, t_s=th * HOUR + stabilize_period, survey=int(i)))
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Model builder (needs pyCATHY)
# -----------------------------------------------------------------------------
def setup_model(simu, case='NC', ic_values=None, extra_times=None, soil_overrides=None, delta=None):
    """Build the model of `case` on a CATHY/DA object. Returns (DEM, ic_nodes).

    ic_values     : layered initial pressure head [[v1], [v2], [v3]] for cfg['ic_ranges'];
                    None = the truth of the case.
    extra_times   : times [s] added to the hourly atmbc/output vector (ERT surveys). Use the SAME
                    list in the truth run and in the DA run.
    soil_overrides: see fill_soil_df / biased_overrides (None = truth parameters).
    delta         : horizontal cell size override [m] (truth and DA must share it).
    """
    import pyCATHY.meshtools as mt              # lazy: the pure helpers above need no pyCATHY
    cfg = get_case(case, delta)
    n, d = cfg['ncell'], cfg['delta']
    ic_values = cfg['ic_true'] if ic_values is None else ic_values

    DEM = np.ones((n, n)) * ELEVATION
    DEM[-1, -1] = ELEVATION - 1e-4
    zratio = [z / MAXDEPTH for z in ZR_THICK]
    assert np.isclose(sum(zratio), 1.0)
    simu.update_prepo_inputs(DEM, delta_x=d, delta_y=d, dr=d, zratio=zratio, base=MAXDEPTH)

    simu.update_veg_map(np.ones(DEM.shape))
    simu.create_mesh_vtk(verbose=False)
    FP = {'PCANA': [-2.0], 'PCREF': [-5.0], 'PCWLT': [-15.0],
          'ZROOT': [0.05], 'PZ': [-1.0], 'OMGC': [1.0]}
    simu.update_soil(FP_map=FP)

    zones = zones_array(cfg)
    nzones = 1
    if zones is not None:
        simu.update_zone(zones)
        nzones = int(zones.max())
    df = simu.init_soil_SPP_map_df(nzones=nzones, nstr=NSTR)
    fill_soil_df(df, cfg, soil_overrides)
    simu.update_soil(SPP_map=df, FP_map=FP, PMIN=-1e12)

    ltop, lbot = mt.get_layer_depths(simu.dem_parameters)
    ic_array = mt.map_param_to_layers(ltop, lbot, cfg['ic_ranges'], ic_values)
    ic_nodes = np.array(simu.map_prop_2mesh_markers('ic', ic_array, to_nodes=True), dtype=float).flatten()
    nnod = (n + 1) ** 2 * (NSTR + 1)
    assert ic_nodes.size == nnod, f"ic has {ic_nodes.size} values, CATHY expects {nnod}"
    simu.update_ic(INDP=1, IPOND=0, pressure_head_ini=ic_nodes)

    times = build_times(extra_times).tolist()
    simu.update_atmbc(HSPATM=0, IETO=1, time=times, netValue=atm_steps(cfg, times))
    simu.update_parm(TIMPRTi=times, TRAFLAG=0, VTKF=2, IPRT1=4, DTMAX=1e3, DELTAT=1e-1, DTMIN=1e-3)

    simu.update_nansfneubc(no_flow=True)
    simu.update_nansfdirbc(no_flow=True)
    simu.update_sfbc(no_flow=True)
    return DEM, ic_nodes
