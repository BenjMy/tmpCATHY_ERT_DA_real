"""ERT FORWARD mesh: region 2 (array block) much finer than region 1 (background).
Markers: boundaries -1 surface (free), -2 other outer faces; cells 1 = background, 2 = array block.

Usage: python make_forward_mesh.py seq3d.shm meshForward [margin_xy margin_z h_reg1 h_reg2 quality n_h2 pad]
"""
import sys, numpy as np, pygimli as pg, pygimli.meshtools as mt

def load_electrodes(shm):
    with open(shm) as f:
        n = int(f.readline().split()[0]); next(f)
        return np.array([[float(v) for v in f.readline().split()[:3]] for _ in range(n)])

def build(el, margin_xy=0.8, margin_z=0.5, h_reg1=2.0, h_reg2=0.1, quality=1.8, n_h2=0, pad=0.15):
    """h_reg1 / h_reg2: max cell VOLUME [m3] in region 1 (coarse) / region 2 (fine).
    n_h2: optional uniform H2 passes (x8 cells each, keeps the region ratio)."""
    lo = el.min(0) - [margin_xy, margin_xy, margin_z]
    hi = el.max(0) + [margin_xy, margin_xy, 0.0]; hi[2] = 0.0
    world = mt.createWorld(start=lo.tolist(), end=hi.tolist(), worldMarker=True)  # -1 surface, -2 outer
    blo = el.min(0) - [pad, pad, pad]; bhi = el.max(0) + [pad, pad, 0.0]
    cube = mt.createCube(size=(bhi - blo).tolist(), pos=((blo + bhi) / 2).tolist(),
                         marker=2, area=h_reg2)
    plc = world + cube
    plc.addRegionMarker([lo[0] + 0.05, lo[1] + 0.05, lo[2] + 0.05], marker=1, area=h_reg1)
    for p in el:
        plc.createNode(p.tolist())
    mesh = mt.createMesh(plc, quality=quality)
    print("tetgen cells", mesh.cellCount())
    for _ in range(int(n_h2)):
        mesh = mesh.createH2()
    return mesh


if __name__ == "__main__":
    shm, out = sys.argv[1], sys.argv[2]
    a = [float(v) for v in sys.argv[3:]]
    if len(a) > 5: a[5] = int(a[5])           # n_h2 must be int
    el = load_electrodes(shm)
    mesh = build(el, *a)
    m = np.array(mesh.cellMarkers())
    print("cells", mesh.cellCount(), "| region1:", int((m == 1).sum()), "region2:", int((m == 2).sum()),
          "| bnd markers", np.unique(mesh.boundaryMarkers()).tolist())
    mesh.save(out + ".bms"); mesh.exportVTK(out + ".vtk")

    v = np.array([c.size() for c in mesh.cells()]); m = np.array(mesh.cellMarkers())
    for k in (1, 2):
        print("region", k, "volume %.3f m3" % v[m == k].sum(), "| cell vol max %.2e  median %.2e" % (v[m == k].max(), np.median(v[m == k])))
        
