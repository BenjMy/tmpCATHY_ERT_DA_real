# Homogeneous-medium test of the ERT forward used by pyCATHY (run in the notebook env).
#   python diag_ert_forward.py ../data/3d/mesh_Luca_Agripolis3d/meshSmall.vtk \
#        ../data/3d/seq3d_from_st5t00tldat.shm [../data/3d/st5NC/procdata/st5t00.tldat]
# A homogeneous model of RHO ohm.m MUST give rhoa == RHO for every quadrupole.
import os, sys
import numpy as np
import pygimli as pg
from pygimli.physics import ert
import pygimli.meshtools as mt

mesh_f, seq_f = sys.argv[1], sys.argv[2]
RHO = 30.0
pre, _ = os.path.splitext(mesh_f)
for ext in (".msh", ".bms", ".vtk"):          # same order as pyCATHY create_ERT_survey_pg
    if os.path.exists(pre + ext):
        mesh = mt.readGmsh(pre + ext) if ext == ".msh" else pg.load(pre + ext)
        print("mesh read from", pre + ext); break
scheme = pg.load(seq_f)
print(mesh, "\nmesh bounds:", mesh.xmin(), mesh.xmax(), mesh.ymin(), mesh.ymax(), mesh.zmin(), mesh.zmax())
pos = np.array(scheme.sensorPositions())
print("electrode bounds:", pos.min(0), pos.max(0), " n =", len(pos))
het = ert.simulate(mesh, res=np.full(mesh.cellCount(), RHO), scheme=scheme, noiseLevel=0, noiseAbs=0, verbose=False)
r = np.asarray(het["rhoa"]) / RHO
bad = np.flatnonzero(np.abs(r - 1) > 0.1)
print(f"rhoa/rho: median {np.median(r):.3f}, p5 {np.percentile(r,5):.3f}, p95 {np.percentile(r,95):.3f}; "
      f"{bad.size}/{r.size} quadrupoles off by >10%")
print("worst:", [(int(het['a'][i]), int(het['b'][i]), int(het['m'][i]), int(het['n'][i]), round(float(r[i]), 2)) for i in np.argsort(-abs(r - 1))[:10]])
if len(sys.argv) > 3:                          # observed k vs numerical k on this mesh
    obs = pg.load(sys.argv[3]); ko = np.asarray(obs["k"]); kf = np.asarray(het["k"])
    if ko.size == kf.size: print("k_fwd/k_obs median", np.median(kf / ko))
    else: print("size differs: obs", ko.size, "fwd", kf.size)
