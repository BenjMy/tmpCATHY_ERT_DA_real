"""Generic ERT data helpers for DA (moved out of the Agripolis ``nc_model.py``).

``filter_tldat`` and ``read_electrodes`` need only numpy; ``unreliable_quadrupoles`` imports
pygimli lazily (homogeneous-medium test of the forward on the mesh the DA really uses).
"""
import logging
import os

import numpy as np

logger = logging.getLogger(__name__)


def read_electrodes(data_file):
    """Electrode positions [n, 3] from a pygimli .shm/.tldat/.dat header
    (line 0 = n, line 1 = '# x y z', then n lines). No pygimli needed."""
    with open(data_file) as f:
        n = int(f.readline().split()[0])
        f.readline()
        return np.array([[float(v) for v in f.readline().split()[:3]] for _ in range(n)])


def filter_tldat(src, dst, err_max=900.0, min_abs_r=0.0, drop=None):
    """Copy a pygimli data file keeping only usable quadrupoles; returns (n_kept, n_total).

    Drops rows with err >= err_max (the 999.99 placeholder of the processing) and, if
    min_abs_r > 0, rows with |r| < min_abs_r. ``drop``: optional bool mask (True = drop) over
    the rows that survive these filters. Electrodes are kept, the count line is rewritten.
    """
    with open(src) as f:
        lines = f.read().split('\n')
    nel = int(lines[0].split()[0])
    i0 = 2 + nel
    ndat = int(lines[i0].split()[0])
    header = lines[i0 + 1]
    cols = header.lstrip('#').split()
    if not {'err', 'r'} <= set(cols):
        raise ValueError(f"unexpected columns in {src}: {cols}")
    ie, ir = cols.index('err'), cols.index('r')
    rows = [l for l in lines[i0 + 2:i0 + 2 + ndat] if l.strip()]
    if len(rows) != ndat:
        raise ValueError(f"{src}: {len(rows)} rows read, {ndat} announced")
    keep = [l for l in rows
            if float(l.split()[ie]) < err_max and abs(float(l.split()[ir])) >= min_abs_r]
    if drop is not None:
        drop = np.asarray(drop, dtype=bool)
        if drop.size != len(keep):
            raise ValueError(f"drop mask ({drop.size}) != kept rows ({len(keep)})")
        keep = [l for l, d in zip(keep, drop) if not d]
    if not keep:
        raise ValueError(f"{src}: every quadrupole filtered out")
    out = lines[:i0] + [str(len(keep)), header] + keep + lines[i0 + 2 + ndat:]
    with open(dst, 'w') as f:
        f.write('\n'.join(out))
    return len(keep), ndat


def set_ert_boundary_markers(mesh):
    """-1 = free surface (top), -2 = mixed/far field on the other outer faces, for a mesh read
    from .vtk (all markers 0 = closed Neumann box, biased potentials). A mesh that already
    carries -2 is left unchanged. Same rule as ``simulate_ERT.create_ERT_survey_pg``."""
    if (np.asarray(mesh.boundaryMarkers()) == -2).any():
        return mesh
    mesh.createNeighborInfos()
    ztop = mesh.zmax()
    for b in mesh.boundaries():
        if b.outside():
            b.setMarker(-1 if abs(b.center().z() - ztop) < 1e-9 else -2)
    return mesh


def unreliable_quadrupoles(mesh_file, scheme_file, tol=0.1, rho=30.0, verbose=True):
    """Bool mask (True = unreliable) over the rows of ``scheme_file``: homogeneous-medium test
    of the forward on EXACTLY ``mesh_file``; rhoa/rho must be 1 within ``tol`` (NaN = unreliable).
    Median far from 1 = global bias (markers, boundaries); median ~1 with wide tails = coarse
    mesh / weak signals."""
    import pygimli as pg
    from pygimli.physics import ert

    mesh = set_ert_boundary_markers(pg.load(mesh_file, verbose=False))
    scheme = pg.load(scheme_file)
    if verbose:
        logger.info("fwd test: %s, %d cells, markers %s", os.path.basename(mesh_file),
                    mesh.cellCount(), np.unique(mesh.boundaryMarkers()).tolist())
    het = ert.simulate(mesh, res=np.full(mesh.cellCount(), float(rho)), scheme=scheme,
                       noiseLevel=0, noiseAbs=0, verbose=False)
    ratio = np.asarray(het['rhoa'], dtype=float) / rho
    bad = ~(np.abs(ratio - 1.0) <= tol)
    if verbose:
        logger.info("fwd test rhoa/rho: median %.3f, 5-95%% [%.3f, %.3f], %d/%d outside +-%.0f%%",
                    np.nanmedian(ratio), np.nanpercentile(ratio, 5), np.nanpercentile(ratio, 95),
                    int(bad.sum()), bad.size, tol * 100)
    return bad
