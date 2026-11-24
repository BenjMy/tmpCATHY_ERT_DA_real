import gmsh
import numpy as np
import sys
from IPython import embed
import pandas as pd
import matplotlib.pyplot as plt
import duckdb as db

gmsh.initialize()

gmsh.model.add("lp")

model = gmsh.model
geo = model.geo
mesh = model.mesh

cl = 0.05

elecs_df = pd.read_csv('elecs72bh_def_mesh.csv')

edges = np.array((
    (-0.1, 1, 0.0),
    (-0.1, -0.1, 0.0),
    (1.1, -0.1, 0.0),
    (1.1, 1.1, 0.0)
))

edgesBottom = edges.copy()
edgesBottom[:, 2] -= 1.3

edges_df = pd.DataFrame(np.array(edges), columns=['x', 'y', 'z'])
edgesBottom_df = pd.DataFrame(np.array(edgesBottom), columns=['x', 'y', 'z'])

surf_df = db.sql(
    """
    select x, y, z, 1 as flag from elecs_df where buried = 0
    union all
    select x, y, z, 0 as flag from edges_df
    """
).df()

surf_xy = surf_df[['x', 'y']].to_numpy()

# Surface points
elecSurfTags = []  # surface electrode tags
for i, p in surf_df.iterrows():
    if p['flag'] == 1:
        pt = geo.addPoint(p['x'], p['y'], p['z'], 0.1, i + 1)
        elecSurfTags.append(pt)
    else:
        pt = geo.addPoint(p['x'], p['y'], p['z'], 1, i + 1)
    geo.synchronize()

# Bottom points
pointBottomTags = []
for i, p in edgesBottom_df.iterrows():
    pt = geo.addPoint(p['x'], p['y'], p['z'], 1, i + 1000)
    pointBottomTags.append(pt)
    geo.synchronize()

# Vertical edges
elecSurfCnt = len(elecSurfTags)
# assuming the 4 surface corners are consecutive to the electrode tags
# and that the matching bottoms start from 1000
vl1 = geo.addLine(elecSurfCnt + 1, 1000)
vl2 = geo.addLine(elecSurfCnt + 2, 1001)
vl3 = geo.addLine(elecSurfCnt + 3, 1002)
vl4 = geo.addLine(elecSurfCnt + 4, 1003)

# Bottom lines
bl1 = geo.addLine(1000, 1001)
bl2 = geo.addLine(1001, 1002)
bl3 = geo.addLine(1002, 1003)
bl4 = geo.addLine(1003, 1000)
# Bottom loop and surface
bottomCurveLoop = geo.addCurveLoop([bl1, bl2, bl3, bl4])
bottomSurf = geo.addPlaneSurface([bottomCurveLoop])

# Redefine surface edges
sl1 = geo.addLine(elecSurfCnt + 1, elecSurfCnt + 2)
sl2 = geo.addLine(elecSurfCnt + 2, elecSurfCnt + 3)
sl3 = geo.addLine(elecSurfCnt + 3, elecSurfCnt + 4)
sl4 = geo.addLine(elecSurfCnt + 4, elecSurfCnt + 1)
# Top loop and surface
topCurveLoop = geo.addCurveLoop([sl1, sl2, sl3, sl4])
topSurf = geo.addPlaneSurface([topCurveLoop])

# Side line loops
sideCurveLoop1 = geo.addCurveLoop([vl1, bl1, -vl2, -sl1])
sideCurveLoop2 = geo.addCurveLoop([vl2, bl2, -vl3, -sl2])
sideCurveLoop3 = geo.addCurveLoop([vl3, bl3, -vl4, -sl3])
sideCurveLoop4 = geo.addCurveLoop([vl4, bl4, -vl1, -sl4])
sideSurf1 = geo.addPlaneSurface([sideCurveLoop1])
sideSurf2 = geo.addPlaneSurface([sideCurveLoop2])
sideSurf3 = geo.addPlaneSurface([sideCurveLoop3])
sideSurf4 = geo.addPlaneSurface([sideCurveLoop4])
sideSurfs = [sideSurf1, sideSurf2, sideSurf3, sideSurf4]

# Volume
surfaceLoop = sideSurfs + [bottomSurf] + [topSurf]
surfloop1 = gmsh.model.geo.addSurfaceLoop(surfaceLoop)
v1 = gmsh.model.geo.addVolume([surfloop1])

# Buried electrodes
elecs_volume = elecs_df.loc[elecs_df['buried'] == 1, ['x', 'y', 'z']]
elecVolumeTags = []
for i, p in elecs_volume.iterrows():
    evt = geo.addPoint(p['x'], p['y'], p['z'], 0.1, 2000 + i)
    elecVolumeTags.append(evt)
geo.synchronize()

mesh.embed(0, elecVolumeTags, 3, v1)
mesh.embed(0, elecSurfTags, 2, topSurf)
geo.synchronize()

# Physical groups for pygimli
elecTags = elecVolumeTags + elecSurfTags
pgrp = gmsh.model.addPhysicalGroup(0, elecTags, 99)
gmsh.model.setPhysicalName(0, pgrp, "Electrodes")
geo.synchronize()
pgrp = gmsh.model.addPhysicalGroup(2, [topSurf], 1, "Neumann")
geo.synchronize()
pgrp = gmsh.model.addPhysicalGroup(2, [bottomSurf] + sideSurfs, 2, "MixedBC")
geo.synchronize()
pgrp = gmsh.model.addPhysicalGroup(3, [v1], 2)
gmsh.model.setPhysicalName(3, pgrp, "invVol")

# Reorganize the tags to handle the duplicates
geo.removeAllDuplicates()
geo.remove_all_duplicates()
geo.synchronize()

refine = True
if refine:
    mesh.field.add("Distance", 1)
    mesh.field.setNumbers(1, "NodesList", elecTags)
    mesh.field.add("Threshold", 2)
    mesh.field.setNumber(2, "IField", 1)
    mesh.field.setNumber(2, "LcMin", 0.1)
    mesh.field.setNumber(2, "LcMax", 1)
    mesh.field.setNumber(2, "DistMin", 0.1)
    mesh.field.setNumber(2, "DistMax", 3)
    mesh.field.setNumber(2, "StopAtDistMax", 1)
    mesh.field.setAsBackgroundMesh(2)

geo.synchronize()

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
