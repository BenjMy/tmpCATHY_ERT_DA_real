# =============================================================================
# Shared nb1 model definition (Agripolis 3D, uncompacted soil)
# =============================================================================
# Used by nb4a_DA_solution_soil3D_nc.py (truth run) and
# nb4_DA_soil3D_nc.py (DA / open-loop run) so both share exactly the same setup.
import numpy as np
import pandas as pd
import pyCATHY.meshtools as mt

HOUR = 3600
DAY = 24 * HOUR
T_END = 3 * DAY

GRANULO = [30, 55, 15]
DELTA = 0.5
X = Y = 5.0
MAXDEPTH = 2.0
ELEVATION = 1e-3
ZR_THICK = [0.03] * 4 + [0.04] * 3 + [0.06] * 2 + [0.08] * 2 + [0.10, 0.25, 0.40, 0.73]
NSTR = len(ZR_THICK)
PHI = {1: 0.415, 2: 0.41, 3: 0.41}
ELSTOR = 1e-5

IC_RANGES = [[0.0, 0.20], [0.2, 0.44], [0.44, 2.0]]
IC_TRUE = [[-2.4], [-1.3], [-1.15]]           # layered "true" initial pressure head

IRR_START_HOUR = 3
ET = -1e-7
IRR = (0.05 / 2) / HOUR


def build_times(extra_times=None):
    """Hourly atmbc/output time vector [s], optionally refined with extra times
    (e.g. ERT survey times) so each of them has an exact CATHY time step.

    If an extra time lies beyond T_END, the hourly vector is extended (to the next full
    hour) so that every survey can be assimilated; without extra times it ends at T_END."""
    t_end = T_END
    extra = None
    if extra_times is not None and len(extra_times) > 0:
        extra = np.asarray(extra_times, dtype=float).ravel()
        extra = extra[extra >= 0]
        if extra.size:
            t_end = max(T_END, float(np.ceil(extra.max() / HOUR) * HOUR))
    times = np.arange(0, t_end + HOUR, HOUR, dtype=float)
    if extra is not None and extra.size:
        times = np.unique(np.round(np.concatenate([times, extra]), 3))
    return times


def snap_to_times(obs_times, atm_times, max_snap_min=1.0):
    """Snap observation times to the closest atmbc time; raise if too far."""
    atm_times = np.asarray(atm_times, dtype=float)
    out = []
    for tt in obs_times:
        t_snap = float(atm_times[np.argmin(np.abs(atm_times - tt))])
        assert abs(t_snap - tt) <= max_snap_min * 60, (
            f"time {tt/HOUR:.2f} h is {abs(t_snap-tt)/60:.0f} min away from the "
            f"closest atmbc time ({t_snap/HOUR:.2f} h): pass it via "
            f"setup_model(..., extra_times=...)")
        out.append(t_snap)
    return out


def obs_window_start(obs_times, atm_times):
    """DA stamp (start of the model window) of observations acquired at obs_times.

    run_DA_sequential runs the model over [t, t_next] and compares an observation
    STAMPED t with the state at t_next (end of the window). Data acquired at time
    T must therefore be stamped with the atmbc time just BEFORE T. obs_times must
    already be on the atmbc grid (see snap_to_times); returns exact floats (the
    library compares the stamps with ==).
    """
    atm = np.asarray(atm_times, dtype=float)
    out = []
    for tt in obs_times:
        i = int(np.flatnonzero(atm == float(tt))[0]) if np.any(atm == float(tt)) else -1
        assert i >= 1, (f"time {tt/HOUR:.2f} h is not an atmbc time after the first one: "
                        f"no model window ends there (snap it first, and t=0 cannot be assimilated)")
        out.append(float(atm[i - 1]))
    return out


def setup_model(simu, ic_values, extra_times=None):
    """Build the nb1 model on a CATHY/DA object. Returns DEM and ic_nodes.

    extra_times: optional list of times [s] (e.g. ERT survey times) added to
    the hourly atmbc/output vector. Use the SAME list in the truth run and in
    the DA run so both share the same time grid.
    """
    ncell = int(np.round(X / DELTA))
    DEM = np.ones((ncell, ncell)) * ELEVATION
    DEM[-1, -1] = ELEVATION - 1e-4

    zratio = [z / MAXDEPTH for z in ZR_THICK]
    assert np.isclose(sum(zratio), 1.0)

    simu.update_prepo_inputs(DEM, delta_x=DELTA, delta_y=DELTA, dr=DELTA,
                             zratio=zratio, base=MAXDEPTH)

    # vegetation
    root_map = np.ones(DEM.shape)
    simu.update_veg_map(root_map)

    
    simu.create_mesh_vtk(verbose=False)

    FP = {'PCANA': [-2.0], 'PCREF': [-5.0], 'PCWLT': [-15.0],
          'ZROOT': [0.05], 'PZ': [-1.0], 'OMGC': [1.0]}
    simu.update_soil(FP_map=FP)

    # soil hydraulic properties
    df = simu.init_soil_SPP_map_df(nzones=1, nstr=NSTR)
    cols = ['PERMX', 'PERMY', 'PERMZ', 'ELSTOR', 'POROS',
            'VGPSATCELL', 'VGNCELL', 'VGRMCCELL']
    idx = pd.IndexSlice
    df.loc[idx[:, 1:6], cols] = [2.1e-06, 2.1e-06, 2.1e-06, ELSTOR, PHI[1], 0.44, 1.4, 0.04]
    df.loc[idx[:, 7:10], cols] = [2.3e-06, 2.2e-06, 2.2e-06, ELSTOR, PHI[2], 0.38, 1.42, 0.065]
    df.loc[idx[:, 10:], cols] = [1.5e-06, 1.5e-06, 1.5e-06, ELSTOR, PHI[3], 0.4, 1.4, 0.07]
    simu.update_soil(SPP_map=df, FP_map=FP, PMIN=-1e12)

    # initial condition
    ltop, lbot = mt.get_layer_depths(simu.dem_parameters)
    ic_array = mt.map_param_to_layers(ltop, lbot, IC_RANGES, ic_values)
    ic_nodes = simu.map_prop_2mesh_markers('ic', ic_array, to_nodes=True)
    ic_nodes = np.array(ic_nodes, dtype=float).flatten()
    simu.update_ic(INDP=1, IPOND=0, pressure_head_ini=ic_nodes)

    # atmospheric forcing (ET everywhere + irrigation on central 2x2 m)
    times = build_times(extra_times).tolist()
    ny, nx = DEM.shape
    cx, cy, hs = nx // 2, ny // 2, int(2.0 / DELTA / 2)
    atm_steps = []
    for t in times:
        atm = np.ones((ny + 1, nx + 1)) * ET
        if t == IRR_START_HOUR * HOUR:
            atm[cy - hs:cy + hs, cx - hs:cx + hs] = IRR
        atm_steps.append(atm)
    simu.update_atmbc(HSPATM=0, IETO=1, time=times, netValue=atm_steps)

    simu.update_parm(TIMPRTi=times, TRAFLAG=0, VTKF=2,
                     IPRT1=4, DTMAX=1e3, DELTAT=1e-1, DTMIN=1e-3)

    # closed lateral/bottom boundaries
    simu.update_nansfneubc(no_flow=True)
    simu.update_nansfdirbc(no_flow=True)
    simu.update_sfbc(no_flow=True)
    return DEM, ic_nodes


# sensor positions (centre of the plot) and the porosity used for sat -> VWC
SENSOR_DEPTHS = [0.15, 0.30, 0.50]
SENSOR_PHI = [PHI[1], PHI[2], PHI[3]]


# =============================================================================
# ERT data helpers (paths, electrodes, acquisition schedule)
# =============================================================================
# Self-contained replacement for the ERT bits of infiltration_util.py
# (set_acronyms, pathName_3dERT, compose_3d_ERT_filenames). Paths are relative to the notebook/ folder (root='../').
import os

# comp_level -> (schedule column, file prefix, folder suffix)
ERT_COMP_LEVELS = {
    'heavy_comp': ('1-HC', 'st1', 'HC'),
    'light_comp': ('3-LC', 'st3', 'LC'),
    'no_comp':    ('5-NC', 'st5', 'NC'),
}


def ert_paths(root='../', comp_level='no_comp'):
    """Return dict with data_dir, sequence, mesh (3D campaign)."""
    _, prefix, suffix = ERT_COMP_LEVELS[comp_level]
    import project_paths as _PP
    _P = _PP.get()
    if _P.is_root(root):            # data_path.toml decides the folders when root is the project root
        return dict(
            data_dir=os.path.join(_P.data_3d, prefix + suffix, 'procdata'),
            sequence=os.path.join(_P.sequence_dir, 'seq3d_from_st5t00tldat.shm'),
            mesh=_P.mesh_forward,
        )
    d3 = os.path.join(root, 'data', '3d')
    return dict(
        data_dir=os.path.join(d3, prefix + suffix, 'procdata'),
        sequence=os.path.join(d3, 'seq3d_from_st5t00tldat.shm'),
        # mesh=os.path.join(d3, 'mesh_Luca_Agripolis3d', 'meshSmall.vtk'),
        mesh=os.path.join(d3, 'mesh_Luca_Agripolis3d', 'meshForward.vtk'),
    )


def load_sensors_from_shm(sequence_file):
    sensors = []

    with open(sequence_file, "r") as f:
        n_sensors = int(f.readline().strip())
        # Skip the sensor-coordinate header, e.g. "# x y z"
        next(f)
        for _ in range(n_sensors):
            line = f.readline().strip()
            sensors.append([float(v) for v in line.split()])
    return sensors



def load_ert_electrodes(sequence_file):
    """Electrode positions [n, 3] of a pygimli sequence/data file."""
    import pygimli as pg                      # lazy: only needed for ERT runs
    
    sensors = load_sensors_from_shm(sequence_file)

    return np.array(sensors)


def load_ert_files(root='../', comp_level='no_comp', times_h=None,
                   stabilize_period=IRR_START_HOUR * HOUR):
    """ERT files + CATHY times, no schedule file needed.

    times_h: {survey index: hours since the START OF IRRIGATION} (negative for
    a baseline before irrigation). Index i -> <prefix>t<ii>.tldat, e.g. 0 ->
    st5t00.tldat. CATHY time = times_h*HOUR + stabilize_period (t=0 in CATHY is
    stabilize_period before irrigation). Missing files are skipped.
    Returns a DataFrame with columns file, t_s.
    """
    _, prefix, _ = ERT_COMP_LEVELS[comp_level]
    data_dir = ert_paths(root, comp_level)['data_dir']
    rows = []
    for i, th in sorted(times_h.items()):
        fname = os.path.join(data_dir, f"{prefix}t{int(i):02d}.tldat")
        if not os.path.isfile(fname):
            print(f"  [skip] ERT file missing: {fname}")
            continue
        rows.append(dict(file=fname, t_s=th * HOUR + stabilize_period))
    return pd.DataFrame(rows)


def filter_tldat(src, dst, err_max=900.0, min_abs_r=0.0, drop=None):
    """Copy a pygimli .tldat/.dat file keeping only usable quadrupoles.

    Drops rows flagged bad by the processing (err >= err_max, e.g. the 999.99 placeholder)
    and, if min_abs_r > 0, rows with |r| < min_abs_r [ohm] (near-zero potential, unstable
    forward response). Electrodes are kept; the data count line is rewritten.
    drop: optional bool mask (True = drop) over the rows that survive the filters.
    Returns (n_kept, n_total).
    """
    with open(src) as f:
        lines = f.read().split('\n')
    nel = int(lines[0].split()[0])
    i0 = 2 + nel                           # count line follows "n", "# x y z" and the electrodes
    ndat = int(lines[i0].split()[0])
    header = lines[i0 + 1]
    cols = header.lstrip('#').split()
    assert {'err', 'r'} <= set(cols), f"unexpected columns in {src}: {cols}"
    ie, ir = cols.index('err'), cols.index('r')
    rows = [l for l in lines[i0 + 2:i0 + 2 + ndat] if l.strip()]
    assert len(rows) == ndat, f"{src}: {len(rows)} rows read, {ndat} announced"
    keep = []
    for l in rows:
        v = l.split()
        if float(v[ie]) >= err_max or abs(float(v[ir])) < min_abs_r:
            continue
        keep.append(l)
    if drop is not None:                   # bool per row of the file AFTER the filters above
        drop = np.asarray(drop, dtype=bool)
        assert drop.size == len(keep), f"drop mask ({drop.size}) != kept rows ({len(keep)})"
        keep = [l for l, d in zip(keep, drop) if not d]
    assert keep, f"{src}: every quadrupole filtered out"
    out = lines[:i0] + [str(len(keep)), header] + keep + lines[i0 + 2 + ndat:]
    with open(dst, 'w') as f:
        f.write('\n'.join(out))
    return len(keep), ndat


def set_ert_boundary_markers(mesh):
    """Give a mesh read from .vtk the boundary markers pyGIMLi's ERT forward needs.

    A .vtk stores no boundary markers: every boundary comes back with marker 0, i.e. a fully
    closed (Neumann) box. That gives biased potentials (rhoa/rho != 1 in a homogeneous medium)
    and the 'Neumann domain without calibration (potential reference) point' warning.
    Convention: -1 = free surface (top face, z = zmax), -2 = mixed/far-field condition on all
    other outer faces. Meshes that already carry -2 (e.g. from .bms/.msh) are left unchanged.
    """
    if (np.asarray(mesh.boundaryMarkers()) == -2).any():
        return mesh
    mesh.createNeighborInfos()
    ztop = mesh.zmax()
    for b in mesh.boundaries():
        if b.outside():
            b.setMarker(-1 if abs(b.center().z() - ztop) < 1e-9 else -2)
    return mesh


def unreliable_quadrupoles(mesh_file, scheme_file, tol=0.1, rho=30.0, verbose=True):
    """Bool mask (True = unreliable) over the rows of scheme_file (pygimli .tldat/.shm).

    Homogeneous-medium test of the ERT FORWARD on EXACTLY `mesh_file` (no sibling
    .msh/.bms lookup), so it validates the same mesh the DA forward reads
    (the patched pyCATHY simulate_ERT reads only the .vtk).

    rhoa must equal rho for every quadrupole; where the finite-element answer is off by more
    than `tol` (small signals, boundaries too close, mesh too coarse near the electrodes) the
    quadrupole cannot be predicted reliably and must not be assimilated. NaN -> unreliable.
    With verbose=True it prints the mesh size, boundary markers, margin around the electrodes
    and the rhoa/rho distribution (median far from 1 -> global bias: boundaries/markers;
    median ~1 with wide tails -> coarse mesh / weak-signal quadrupoles).
    """
    import pygimli as pg                      # lazy: only needed for ERT runs
    from pygimli.physics import ert

    mesh = set_ert_boundary_markers(pg.load(mesh_file, verbose=False))
    scheme = pg.load(scheme_file)

    if verbose:
        print(f"[fwd test] mesh: {os.path.basename(mesh_file)}, {mesh.cellCount()} cells, dim={mesh.dim()}")
        print(f"[fwd test] boundary markers: {np.unique(mesh.boundaryMarkers()).tolist()} "
              "(expect -1 surface, -2 outer/mixed)")
        try:
            pos = np.array([[p.x(), p.y(), p.z()] for p in scheme.sensors()])
            lo = np.array([mesh.xmin(), mesh.ymin(), mesh.zmin()])
            hi = np.array([mesh.xmax(), mesh.ymax(), mesh.zmax()])
            print(f"[fwd test] margin around electrodes [m]: min side {np.round(pos.min(0) - lo, 2).tolist()}, "
                  f"max side {np.round(hi - pos.max(0), 2).tolist()}")
        except Exception as e:                # diagnostic only
            print(f"[fwd test] margin not computed ({e})")

    het = ert.simulate(mesh, res=np.full(mesh.cellCount(), float(rho)), scheme=scheme,
                       noiseLevel=0, noiseAbs=0, verbose=False)
    ratio = np.asarray(het['rhoa'], dtype=float) / rho
    bad = ~(np.abs(ratio - 1.0) <= tol)       # NaN -> unreliable

    if verbose:
        print(f"[fwd test] rhoa/rho: median {np.nanmedian(ratio):.3f}, "
              f"5-95% [{np.nanpercentile(ratio, 5):.3f}, {np.nanpercentile(ratio, 95):.3f}], "
              f"{int(bad.sum())}/{bad.size} outside +-{tol:.0%}")
    return bad


# =============================================================================
# Per-material VGP control (pyCATHY 'material' perturbation control)
# =============================================================================
# Layers (0-based, from the top) sharing the same soil material in setup_model (A, B, C).
# Single source of truth for the scenario conversion below (the label slices of setup_model
# overlap at layer 10 - plan item P13 - so they are NOT used here).
MATERIAL_LAYERS = [list(range(0, 6)), list(range(6, 9)), list(range(9, NSTR))]
MATERIAL_PARAMETERS = ('SATCELL_VG', 'n_VG', 'thetar_VG')


def _spec_bounds(spec):
    """(lo, hi) of a bounds spec given as [lo, hi] or {min, max}; None = open side."""
    if isinstance(spec, dict):
        return spec.get('min'), spec.get('max')
    return spec[0], spec[1]


def scenario_to_material(sc, groups=MATERIAL_LAYERS, verbose=True):
    """Copy of the per-layer VGP scenario ``sc`` with the VGP controlled per material.

    For each of SATCELL_VG / n_VG / thetar_VG with the 'layers' control: one entry per group
    (nominal, mean, sigma = mean over the layers of the group), pert_control_name 'material',
    per_bounds reduced to one spec (the intersection over the layers), listUpdateParm renamed
    (NAME0..NAME{G-1}). Every change is printed ([DA-PARAM] lines). Other parameters untouched.
    """
    import copy
    out = copy.deepcopy(sc)
    names = out['per_name']
    controls = out.setdefault('pert_control_name', ['None'] * len(names))
    out['material_layers'] = [list(g) for g in groups]
    changed, lines = [], []
    for i, name in enumerate(names):
        if name not in MATERIAL_PARAMETERS or controls[i] != 'layers':
            continue
        for key in ('per_nom', 'per_mean', 'per_sigma'):
            v = np.asarray(out[key][i], dtype=float)
            assert v.size == NSTR, (name, key, v.size, NSTR)
            gv = [v[g] for g in groups]
            if key == 'per_nom' and any(np.ptp(x) > 1e-12 * max(1.0, abs(x).max()) for x in gv):
                lines.append(f"WARNING {name}: nominals differ inside a material {[x.tolist() for x in gv]}: "
                             "the group mean is used")
            out[key][i] = [float(x.mean()) for x in gv]
        for key in ('per_type', 'transf_type', 'sampling_type'):
            if key in out and isinstance(out[key][i], list):
                first = [out[key][i][g[0]] for g in groups]
                assert all(len({out[key][i][l] for l in g}) == 1 for g in groups), (
                    f"{name}: {key} differs inside a material")
                out[key][i] = first
        pb = out.get('per_bounds')
        if pb is not None and i < len(pb) and isinstance(pb[i], (list, tuple)) and len(pb[i]) == NSTR:
            los = [b for b in (_spec_bounds(s)[0] for s in pb[i]) if b is not None]
            his = [b for b in (_spec_bounds(s)[1] for s in pb[i]) if b is not None]
            lo, hi = (max(los) if los else None), (min(his) if his else None)
            pb[i] = [lo, hi]
            lines.append(f"{name}: per-layer bounds reduced to one pair [{lo}, {hi}] (intersection)")
        controls[i] = 'material'
        changed.append(name)
        lines.append(f"{name}: 'layers' ({NSTR} params) -> 'material' ({len(groups)} params); "
                     f"nominals {[round(x, 6) for x in out['per_nom'][i]]}, "
                     f"sigmas {[round(x, 6) for x in out['per_sigma'][i]]}")
    upd, seen = [], set()
    for u in out['listUpdateParm']:
        root = u.rstrip('0123456789')
        if root in changed:
            if root not in seen:
                seen.add(root)
                upd += [f'{root}{k}' for k in range(len(groups))]
        else:
            upd.append(u)
    if verbose:
        print('[DA-PARAM] ===== VGP control converted to MATERIAL =====')
        print(f'[DA-PARAM] material_layers = {out["material_layers"]}')
        for l in lines:
            print('[DA-PARAM]', l)
        print(f'[DA-PARAM] listUpdateParm: {len(sc["listUpdateParm"])} -> {len(upd)} entries: {upd}')
    out['listUpdateParm'] = upd
    return out
