"""
Infiltration Utility Functions - Refactored and Documented
===========================================================

This module contains utility functions for infiltration experiments combined with
Electrical Resistivity Tomography (ERT). These functions help:

1. Manage experiment data and file naming conventions
2. Set up soil meshes and properties
3. Process sensor observations (water content, pressure)
4. Handle ERT measurements and forward modeling
5. Visualize and analyze results

Key Concepts:
-------------
- ERT: Electrical Resistivity Tomography - measures soil electrical properties
  to infer water content (wet soil conducts electricity better than dry)
- Infiltration: Water entering and moving through soil
- Mesh: 3D grid dividing soil into cells for simulation
- Sensors: Physical devices that measure soil conditions (moisture, pressure)
- Forward modeling: Predict what measurements would look like given a soil state

Author: [Original author]
Last modified: [Date]
"""

# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import re
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import os 
from pathlib import Path

# pyCATHY modules for hydrological modeling
from pyCATHY.plotters import cathy_plots as pltCT
import pyCATHY.meshtools as mt
from pyCATHY.petro import predict_unsat_soil_hydr_param
from pyCATHY.meshtools import CATHY_2_pg
from pyCATHY.DA.observations import read_observations

# PyGIMLi for geophysical modeling (ERT)
import pygimli as pg

# pyCATHY ERT modules
from pyCATHY.ERT import petro_Archie as Archie
from pyCATHY.ERT import simulate_ERT as simuERT

# PyVista for 3D mesh visualization
import pyvista as pv


# =============================================================================
# EXPERIMENT NAMING AND ORGANIZATION
# =============================================================================

def set_acronyms(comp_level):
    """
    Define naming conventions for different soil compaction levels.
    
    Soil compaction affects how water moves through soil:
    - Heavy compaction (HC): Soil is densely packed, slower infiltration
    - Light compaction (LC): Moderately packed soil
    - No compaction (NC): Loose soil, faster infiltration
    
    Parameters
    ----------
    comp_level : str
        Compaction level: 'heavy_comp', 'light_comp', or 'no_comp'
    
    Returns
    -------
    comp_level_acr : dict
        Acronyms for each compaction level (e.g., 'heavy_comp' -> '1-HC')
    comp_level_st : dict
        Station identifiers for each level (e.g., 'heavy_comp' -> ['st1', 'HC'])
    
    Example
    -------
    >>> acr, st = set_acronyms('heavy_comp')
    >>> print(acr['heavy_comp'])
    '1-HC'
    """
    comp_level_acr = {
        'heavy_comp': '1-HC',
        'light_comp': '3-LC',       
        'no_comp': '5-NC',
    }

    comp_level_st = {
        'heavy_comp': ['st1', 'HC'],
        'light_comp': ['st3', 'LC'],       
        'no_comp': ['st5', 'NC'],
    }
    
    return comp_level_acr, comp_level_st


def compose_2d_ERT_filenames(rootpath, index, comp_level_st, comp_level):
    """
    Generate file paths for 2D ERT measurement files.
    
    2D ERT creates a cross-section view of the soil (like a vertical slice).
    This function builds the correct file paths based on measurement timing.
    
    Parameters
    ----------
    rootpath : str
        Base directory containing ERT data
    index : array-like
        Measurement indices (e.g., [0, 1, 2, ...] for sequential measurements)
    comp_level_st : dict
        Station naming dictionary from set_acronyms()
    comp_level : str
        Compaction level identifier
    
    Returns
    -------
    ERT_time_names : list
        Complete file paths to each ERT measurement file
    """
    ERT_time_names = [
        os.path.join(
            '../data/2d/',
            comp_level_st[comp_level][0] + comp_level_st[comp_level][1],
            comp_level_st[comp_level][0] + '0' + str(r) + '.csv'
        ) 
        for r in index
    ]    
    return ERT_time_names


def compose_3d_ERT_filenames(rootpath, index, comp_level_st, comp_level):
    """
    Generate file paths for 3D ERT measurement files.
    
    3D ERT provides a complete volumetric view of the soil.
    Handles different naming conventions for measurements < 10 and >= 10.
    
    Parameters
    ----------
    rootpath : str
        Base directory containing ERT data
    index : array-like
        Measurement indices
    comp_level_st : dict
        Station naming dictionary
    comp_level : str
        Compaction level identifier
    
    Returns
    -------
    ERT_time_names : list
        Sorted list of complete file paths to 3D ERT measurement files
    
    Note
    ----
    Files are stored in the PyGIMLi ERTManager format for processing
    """
    # Build file paths for measurements 0-9 (with leading zero: 't01', 't02', etc.)
    ERT_time_names = [
        os.path.join(
            rootpath,
            comp_level_st[comp_level][0] + comp_level_st[comp_level][1],
            'gimli_inv', 
            comp_level_st[comp_level][0] + 't0' + str(r) + '/ERTManager/data.dat'
        ) 
        for r in index if int(r) < 10
    ]

    # Build file paths for measurements 10+ (no leading zero: 't10', 't11', etc.)
    ERT_time_namesSUp = [
        os.path.join(
            rootpath,
            comp_level_st[comp_level][0] + comp_level_st[comp_level][1],
            'gimli_inv',
            comp_level_st[comp_level][0] + 't' + str(r) + '/ERTManager/data.dat'
        ) 
        for r in index if int(r) >= 10
    ]
    
    # Combine and sort all file paths
    ERT_time_names = ERT_time_names + ERT_time_namesSUp
    ERT_time_names = sorted(ERT_time_names)
    
    return ERT_time_names

    
def old_compose_3d_ERT_filenames(rootpath, index, comp_level_st, comp_level):
    """
    Generate file paths for 3D ERT files using old naming convention.
    
    This function exists for backward compatibility with older data formats.
    Uses '.tldat' extension instead of the newer ERTManager format.
    
    Parameters
    ----------
    rootpath : str
        Base directory containing ERT data
    index : array-like
        Measurement indices
    comp_level_st : dict
        Station naming dictionary
    comp_level : str
        Compaction level identifier
    
    Returns
    -------
    ERT_time_names : list
        Sorted list of file paths using old format
    """
    # Files 0-9 with leading zero
    ERT_time_names = [
        os.path.join(
            rootpath,
            comp_level_st[comp_level][0] + comp_level_st[comp_level][1],
            'procdata',
            comp_level_st[comp_level][0] + 't0' + str(r) + '.tldat'
        ) 
        for r in index if int(r) < 10
    ]

    # Files 10+ without leading zero
    ERT_time_namesSUp = [
        os.path.join(
            rootpath,
            comp_level_st[comp_level][0] + comp_level_st[comp_level][1],
            'procdata',
            comp_level_st[comp_level][0] + 't' + str(r) + '.tldat'
        ) 
        for r in index if int(r) >= 10
    ]
    
    ERT_time_names = ERT_time_names + ERT_time_namesSUp
    ERT_time_names = sorted(ERT_time_names)
    
    return ERT_time_names


# =============================================================================
# EXPERIMENT SCHEDULE AND TIMING
# =============================================================================

def read_infiltration_log(comp_level, ERTdim='2d', **kwargs):
    """
    Read and process the infiltration experiment schedule.
    
    The log file contains crucial timing information:
    - When water was applied to the soil (infiltration start/stop)
    - When ERT measurements were taken
    - How much water was applied
    
    This allows us to synchronize simulations with actual measurements.
    
    Parameters
    ----------
    comp_level : str
        Soil compaction level ('heavy_comp', 'light_comp', or 'no_comp')
    ERTdim : str, optional
        ERT dimensionality: '2d' or '3d' (default: '2d')
    
    Returns
    -------
    infitration_log_selec : pd.DataFrame
        Filtered schedule with relevant experiment information:
        - Timing of each measurement
        - Infiltration duration
        - Water volume applied
    ERT_time_names : list
        File paths to corresponding ERT measurement files
    
    Notes
    -----
    The schedule file 'schedule2023_3d.csv' contains the master timing
    for the 2023 infiltration experiments at Agripolis.
    """
    # Get naming conventions for this compaction level
    comp_level_acr, comp_level_st = set_acronyms(comp_level)
    
    # Read the master schedule file
    infitration_log = pd.read_csv(
        'schedule2023_3d.csv',
        sep=',',
        date_format='%d/%m/%Y %H:%M',
        dayfirst=True  # European date format (day before month)
    )
    
    # Extract relevant columns for this compaction level
    infitration_log_selec = infitration_log[[
        'week',           # Week of experiment
        'day',            # Day within week
        'dayT',           # Total day number
        'timestep',       # Measurement timestep
        comp_level_acr[comp_level]  # Measurement datetime for this treatment
    ]]    
    
    # Infiltration duration for each treatment (in minutes)
    # These were measured during the actual experiments
    inf_len3d = {
        '5-NC': 30,   # No compaction: 30 minutes
        '3-LC': 65,   # Light compaction: 65 minutes
        '1-HC': 62,   # Heavy compaction: 62 minutes
    }
    
    # Exact start times for each infiltration event
    inf_dateStart = {
        '5-NC': '10/09/2023 11:10:00', 
        '3-LC': '10/09/2023 12:19:00',
        '1-HC': '10/09/2023 14:08:00',
    }
    
    # Add infiltration metadata to the schedule
    infitration_log_selec['infiltration_length (min)'] = inf_len3d[comp_level_acr[comp_level]]
    infitration_log_selec['inf_dateStart'] = pd.to_datetime(inf_dateStart[comp_level_acr[comp_level]])
    infitration_log_selec['infiltration_water (l/4m^2)'] = 300  # 300 liters over 4 m² area

    # Get corresponding ERT file names for these measurement times
    ERT_time_names = compose_3d_ERT_filenames(
        '../data/3d/',
        infitration_log_selec.index.values,
        comp_level_st,
        comp_level
    )
    
    # Convert measurement times to datetime objects for easier manipulation
    date_format = "%d/%m/%Y %H:%M"
    datetime2convert = infitration_log_selec[comp_level_acr[comp_level]]
    infitration_log_selec.loc[:, comp_level_acr[comp_level]] = pd.to_datetime(
        datetime2convert,
        format=date_format
    )
    
    return infitration_log_selec, ERT_time_names


def get_ERT_timesteps(infitration_log_selec, compactLevel, stabilizePeriod=3*24*60*60):
    """
    Calculate exact timing for ERT measurements relative to infiltration.
    
    This function converts absolute datetimes to elapsed seconds, which is
    what the simulation needs. It also identifies when infiltration started
    and stopped relative to the first measurement.
    
    Parameters
    ----------
    infitration_log_selec : pd.DataFrame
        Infiltration schedule from read_infiltration_log()
    compactLevel : str
        Compaction level identifier
    stabilizePeriod : int, optional
        Stabilization period in seconds before infiltration (default: 3 days)
        This allows soil to equilibrate before the experiment starts
    
    Returns
    -------
    irr_timeStart_refERT : float
        Infiltration start time in seconds relative to first ERT measurement
    irr_timeStop_refERT : float
        Infiltration stop time in seconds relative to first ERT measurement
    irr_time_s : float
        Total infiltration duration in seconds
    elapsed_time_ERT_s : np.array
        Elapsed time (seconds) for each ERT measurement
    ERT_datetimes : pd.Series
        Actual datetime stamps for each ERT measurement
    
    Example
    -------
    If first ERT at 10:00:00, infiltration at 10:05:00 for 30 minutes:
    - irr_timeStart_refERT = 300 seconds (5 minutes)
    - irr_timeStop_refERT = 2100 seconds (35 minutes)
    - irr_time_s = 1800 seconds (30 minutes)
    """
    # Get naming convention for this compaction level
    comp_level_acr, _ = set_acronyms(compactLevel)

    # Get the infiltration start time
    T0_infiltration = infitration_log_selec['inf_dateStart'][0]
    
    # Get all ERT measurement times
    ERT_datetimes = infitration_log_selec[comp_level_acr[compactLevel]]
    
    # Calculate elapsed time for each ERT measurement relative to infiltration start
    elapsed_time_ERT_s = [
        0,  # First measurement at time zero
        np.hstack([
            (Ti - T0_infiltration).total_seconds() 
            for Ti in infitration_log_selec[comp_level_acr[compactLevel]][1:]
        ])
    ]
    elapsed_time_ERT_s = np.hstack(elapsed_time_ERT_s)
    
    # Calculate infiltration duration in seconds
    irr_time_s = infitration_log_selec['infiltration_length (min)'].iloc[0] * 60
    print(f'Duration of infiltration: {irr_time_s} seconds ({irr_time_s/60:.1f} minutes)')

    # Calculate when infiltration starts relative to first ERT measurement
    irr_timeStart_refERT = (T0_infiltration - ERT_datetimes[0]).total_seconds()
    
    # Calculate when infiltration stops
    irr_timeStop_refERT = irr_timeStart_refERT + irr_time_s 
   
    return (
        irr_timeStart_refERT, 
        irr_timeStop_refERT, 
        irr_time_s, 
        elapsed_time_ERT_s, 
        ERT_datetimes
    )


# =============================================================================
# MESH CREATION - Defining the 3D Grid for Simulation
# =============================================================================

def create_hydromesh_3d(simuInfiltration):
    """
    Create a 3D computational mesh for hydrological simulation.
    
    The mesh is a 3D grid that divides the soil volume into small cells.
    Each cell tracks water content, pressure, and other properties over time.
    
    Think of it like dividing a block of soil into many small cubes - the
    simulation calculates what happens in each cube.
    
    Parameters
    ----------
    simuInfiltration : CATHY object
        CATHY simulation object to configure
    
    Returns
    -------
    cell_size : float
        Surface area of each grid cell (m²)
        Smaller cells = more accurate but slower simulation
    
    Mesh Configuration
    ------------------
    - Horizontal resolution: 0.1 m × 0.1 m cells
    - Domain size: 1.5 m × 1.5 m surface area
    - Depth: 1.5 m deep
    - Vertical layers: 12 layers with varying thickness
    
    Notes
    -----
    The 3D ERT electrode array covers a much larger area (14.2 m profile),
    but we simulate a smaller representative volume for computational efficiency.
    """
    # Define horizontal grid spacing (resolution)
    delta_x = 0.1  # meters (10 cm cells in x-direction)
    delta_y = 0.1  # meters (10 cm cells in y-direction)
    
    # Alternative coarser resolution (commented out for faster testing):
    # delta_x = 1  # 1 meter cells
    # delta_y = 1
    
    # Define total domain dimensions
    x = 1.5  # meters (domain width)
    y = 1.5  # meters (domain length)
    maxdepth = 1.5  # meters (domain depth)
    
    # Calculate number of cells in each horizontal direction
    ncellx = int(np.round(x / delta_x))  # Number of cells in x
    ncelly = int(np.round(y / delta_y))  # Number of cells in y
    
    # Create Digital Elevation Model (DEM) - surface topography
    # Using nearly flat surface at 1 mm elevation
    DEM = np.ones([ncelly, ncellx]) * 1e-3  # 1 mm elevation
    DEM[-1, -1] = 1e-3 - 1e-4  # Slight slope in corner for drainage
    
    # Adjust maximum depth (leave some space at bottom)
    maxdepth = maxdepth - 1e-1
    
    # Define vertical layer boundaries (increasing spacing with depth)
    # More layers near surface where changes happen quickly
    zb = np.linspace(0, maxdepth, 12)  # 12 vertical layers
    
    nstr = len(zb)  # Number of soil layers
    
    # Calculate layer thickness ratios (uniform in this case)
    zr = list((np.ones(len(zb))) / nstr)
    
    # Update CATHY preprocessor inputs with mesh parameters
    simuInfiltration.update_prepo_inputs(
        DEM,
        delta_x=delta_x,   # x resolution (cell size)
        delta_y=delta_y,   # y resolution
        dr=delta_x,        # radial distance for cylindrical grids
        zratio=zr,         # vertical layer thickness ratios
        base=max(zb),      # depth of model bottom
    )
    
    # Calculate cell surface area (for flux calculations)
    cell_size = delta_x * delta_y
    
    # Display the DEM
    simuInfiltration.show_input('dem')
    
    # Update DEM parameters in the model
    simuInfiltration.update_dem_parameters()
    
    # Run the CATHY preprocessor to generate the 3D mesh
    simuInfiltration.run_preprocessor(IPRT1=3)
    
    # Create VTK file for 3D visualization
    simuInfiltration.create_mesh_vtk(verbose=False)
    
    return cell_size


def create_hydromesh_2d(simuInfiltration):
    """
    Create a 2D computational mesh for hydrological simulation.
    
    2D meshes are used for cross-section simulations - like looking at a
    vertical slice through the soil. This is faster than 3D but provides
    less spatial information.
    
    Parameters
    ----------
    simuInfiltration : CATHY object
        CATHY simulation object to configure
    
    Returns
    -------
    cell_size : float
        Surface area of each grid cell (m²)
    
    Mesh Configuration
    ------------------
    - Matches 2D ERT profile: 14.2 m long with 72 electrodes
    - Electrode spacing: 0.2 m (20 cm)
    - Depth: 2.2 m
    - Vertical resolution: 15 layers
    
    Notes
    -----
    The 2D mesh is designed to exactly match the electrode positions in
    the field ERT measurements.
    """
    # Define grid spacing to match ERT electrode positions
    delta_x = 0.2  # 20 cm spacing (matches electrode spacing)
    delta_y = 0.2
    
    # Calculate number of cells to cover the 14.2 m ERT profile
    ncellx = int(np.round(14.2 / delta_x))  # ~71 cells
    ncelly = int(np.round(2 / delta_y))     # Width for 2D simulation
    
    # Create DEM with slight slope for drainage
    DEM = np.ones([ncelly, ncellx]) * 1e-1  # 10 cm elevation
    DEM[-1, -1] = 0.1 - 1e-3  # Corner slightly lower
    
    # Define depth and vertical layers
    maxdepth = 2.2 - 1e-1  # 2.1 m depth
    zb = np.linspace(0, maxdepth, 15)  # 15 vertical layers
    
    nstr = len(zb)
    zr = list((np.ones(len(zb))) / nstr)
    
    # Update mesh parameters
    simuInfiltration.update_prepo_inputs(
        DEM,
        delta_x=delta_x,
        delta_y=delta_y,  
        dr=delta_x,
        zratio=zr,
        base=max(zb),
    )
    
    cell_size = delta_x * delta_y
    
    # Generate the mesh
    simuInfiltration.show_input('dem')
    simuInfiltration.update_dem_parameters()
    simuInfiltration.run_preprocessor(IPRT1=3)
    simuInfiltration.create_mesh_vtk(verbose=False)
    
    return cell_size


# =============================================================================
# IRRIGATION ZONE DEFINITION
# =============================================================================

def bounds_2d_irrigation():
    """
    Define the irrigation zone for 2D simulations.
    
    Returns
    -------
    xb : list
        X-coordinate boundaries [start, end] in meters
    yb : list
        Y-coordinate boundaries [start, end] in meters
    
    Notes
    -----
    Irrigation is applied to a 2 m section in the middle of the 14.2 m profile
    (from 6.1 m to 8.1 m along the profile)
    """
    xb = [6.1, 8.1]  # 2 meter irrigation zone
    yb = [0, 2]      # Full width in y-direction
    return xb, yb


def bounds_3d_irrigation(meshCATHY_pv):
    """
    Define the irrigation zone for 3D simulations.
    
    Creates a 1 m × 1 m square irrigation zone centered on the mesh.
    This represents a smaller infiltration basin (similar to field setup).
    
    Parameters
    ----------
    meshCATHY_pv : pv.UnstructuredGrid
        PyVista mesh object containing the CATHY mesh
    
    Returns
    -------
    xb : list
        X-coordinate boundaries [start, end] in meters
    yb : list
        Y-coordinate boundaries [start, end] in meters
    
    Notes
    -----
    The irrigation zone is centered on the mesh, extending ±0.5 m from center
    """
    # Find center of mesh and create 1m x 1m irrigation zone
    xb = [
        meshCATHY_pv.points[:, 0].max() / 2 - 0.5,  # Center minus 0.5 m
        meshCATHY_pv.points[:, 0].max() / 2 + 0.5   # Center plus 0.5 m
    ]
    yb = [
        meshCATHY_pv.points[:, 1].max() / 2 - 0.5,
        meshCATHY_pv.points[:, 1].max() / 2 + 0.5
    ]
    return xb, yb


def get_nodes_irrigated(simuInfiltration, hydromesh=None, dim='3d'):
    """
    Identify which mesh nodes are within the irrigation zone.
    
    This function finds all surface nodes where water should be applied
    during the infiltration event.
    
    Parameters
    ----------
    simuInfiltration : CATHY object
        CATHY simulation with mesh information
    hydromesh : pv.UnstructuredGrid, optional
        PyVista mesh object (will be read if not provided)
    dim : str, optional
        Dimension: '2d' or '3d' (default: '3d')
    
    Returns
    -------
    surfnodes : np.array
        Coordinates of all surface nodes (shape: [n_nodes, 3])
    within_boundaries : np.array
        Boolean array indicating which surface nodes are in irrigation zone
    
    Example
    -------
    If you have a 10×10 surface grid and irrigate the center 2×2 cells,
    within_boundaries would be True for those 4 cells and False elsewhere.
    """
    # Read the 3D grid information
    grid3d = simuInfiltration.read_outputs('grid3d')
    nnod = grid3d['nnod']  # Number of surface nodes
    
    # Extract surface node coordinates (top layer of mesh)
    surfnodes = grid3d['mesh3d_nodes'][0:int(nnod)]
    
    # Read mesh if not provided
    if hydromesh is None:
        hydromesh = pv.read(
            os.path.join(
                simuInfiltration.workdir,
                simuInfiltration.project_name,
                'vtk',
                simuInfiltration.project_name + '.vtk'
            )
        )
    
    # Get irrigation zone boundaries
    if dim == '2d':
        xb, yb = bounds_2d_irrigation()
    else:
        xb, yb = bounds_3d_irrigation(hydromesh)
    
    # Get mesh resolution
    dx = simuInfiltration.dem_parameters['delta_x']
    dy = simuInfiltration.dem_parameters['delta_y']
    
    # Find nodes within irrigation boundaries
    condx = (surfnodes[:, 0] > xb[0]) & (surfnodes[:, 0] < xb[1])
    condy = (surfnodes[:, 1] > yb[0]) & (surfnodes[:, 1] < yb[1])
    within_boundaries = condx & condy
    
    return surfnodes, within_boundaries


# =============================================================================
# INFILTRATION RATE CALCULATION
# =============================================================================

def get_Qinfiltration(infitration_log_selec, irr_time_s, cell_size, compactLevel):
    """
    Calculate the infiltration rate (water flux) during irrigation.
    
    This converts the total volume of water applied into a rate (m/s) that
    can be used as a boundary condition in the simulation.
    
    Parameters
    ----------
    infitration_log_selec : pd.DataFrame
        Infiltration schedule with water volume information
    irr_time_s : float
        Duration of infiltration in seconds
    cell_size : float
        Surface area of mesh cells (m²)
    compactLevel : str
        Compaction level identifier
    
    Returns
    -------
    Q_infiltration_m : float
        Infiltration rate in m/s (flux applied at soil surface)
    
    Calculation Steps
    -----------------
    1. Get total water volume (liters) applied over 4 m² area
    2. For some treatments, only half the water is used
    3. Convert liters to cubic meters: L / 1000
    4. Divide by area to get water depth: m³ / 4 m² = m
    5. Divide by time to get rate: m / s
    
    Example
    -------
    300 liters over 4 m² for 30 minutes (1800 seconds):
    - Half water (if applicable): 150 L = 0.15 m³
    - Depth: 0.15 m³ / 4 m² = 0.0375 m
    - Rate: 0.0375 m / 1800 s = 2.08e-5 m/s
    """
    print(f'Infiltration duration: {irr_time_s} seconds ({irr_time_s/60:.1f} minutes)')
    
    # Get total water volume from log
    if compactLevel == 'no_comp':
        # No compaction: use half the water (soil infiltrates faster)
        infiltration_l = infitration_log_selec['infiltration_water (l/4m^2)'].iloc[0] / 2
        water_depth_m = (infiltration_l / 1000) / 4  # Convert L to m³, divide by area
        
    elif compactLevel == 'heavy_comp':
        # Heavy compaction: also use half the water
        infiltration_l = infitration_log_selec['infiltration_water (l/4m^2)'].iloc[0] / 2
        water_depth_m = (infiltration_l / 1000) / 4
    
    else:
        # Other compaction levels: use full amount
        infiltration_l = infitration_log_selec['infiltration_water (l/4m^2)'].iloc[0]
        water_depth_m = (infiltration_l / 1000) / 4
    
    # Calculate infiltration rate (m/s)
    Q_infiltration_m = water_depth_m / irr_time_s
    
    print(f'Water applied: {infiltration_l:.1f} L over 4 m²')
    print(f'Water depth: {water_depth_m:.4f} m')
    print(f'Infiltration rate: {Q_infiltration_m:.2e} m/s')
    
    return Q_infiltration_m


# =============================================================================
# ATMOSPHERIC BOUNDARY CONDITIONS (Water Application Schedule)
# =============================================================================

def prepare_atmbc_infiltration(simuInfiltration, ERT_time_names, 
                              infitration_log_selec, cell_size, 
                              compactLevel, atmbc_dim):
    """
    Prepare atmospheric boundary conditions for the infiltration simulation.
    
    Atmospheric boundary conditions (atmbc) define when and where water enters
    the soil surface. This function creates a time series of water fluxes that
    matches the experimental infiltration schedule.
    
    Think of it as programming a sprinkler system: when to turn on, how much
    water to apply, and when to turn off.
    
    Parameters
    ----------
    simuInfiltration : CATHY object
        Simulation object
    ERT_time_names : list
        ERT measurement file paths (defines observation times)
    infitration_log_selec : pd.DataFrame
        Infiltration schedule
    cell_size : float
        Mesh cell surface area (m²)
    compactLevel : str
        Soil compaction level
    atmbc_dim : str
        Dimension: '2d' or '3d'
    
    Returns
    -------
    atmbc_times : np.array
        Times (seconds) when boundary conditions change
        Includes: measurement times, infiltration start/stop, drainage period
    v_atmbc : np.array
        Water flux values (m/s) at each time
        Shape: (n_times, n_surface_nodes) for 2D or (n_times,) for 3D
    elapsed_time_ERT_s : np.array
        Elapsed time for ERT measurements
    irr_timeStart_refERT : float
        Infiltration start time (seconds)
    irr_timeStop_refERT : float
        Infiltration stop time (seconds)
    
    Time Points Included
    --------------------
    1. ERT measurement times (when we want to compare with observations)
    2. Infiltration start and stop (when water is applied)
    3. Additional points for smooth simulation:
       - 1000 s before/after infiltration start/stop
       - Multiple points during drainage period (up to 2 days)
    
    Example Timeline
    ----------------
    Time 0:     First ERT measurement
    Time 300:   Infiltration starts (water on)
    Time 2100:  Infiltration stops (water off)
    Time 3000:  ERT measurement during drainage
    Time 86400: Final measurement after 1 day
    """
    # Get naming conventions
    comp_level_acr, _ = set_acronyms(compactLevel)
    
    # Calculate timing information
    (irr_timeStart_refERT,   # When infiltration starts
     irr_timeStop_refERT,    # When infiltration stops
     irr_time_s,             # Duration of infiltration
     elapsed_time_ERT_s,     # Times of ERT measurements
     ERT_datetimes) = get_ERT_timesteps(infitration_log_selec, compactLevel)
    
    # Calculate infiltration rate
    Q_infiltration_m = get_Qinfiltration(
        infitration_log_selec,
        irr_time_s,
        cell_size,
        compactLevel
    )
    
    # Create comprehensive time array including:
    # - ERT measurement times
    # - Infiltration start/stop and nearby times
    # - Additional drainage period times
    atmbc_times = np.unique(np.sort(np.r_[
        elapsed_time_ERT_s,                      # ERT measurement times
        irr_timeStart_refERT - 1000,             # Before infiltration starts
        irr_timeStart_refERT,                    # Infiltration starts
        irr_timeStop_refERT - 1000,              # Near infiltration stop
        irr_timeStop_refERT,                     # Infiltration stops
        irr_timeStop_refERT + 1000,              # After infiltration stops
        irr_timeStop_refERT + 25000,             # Drainage period points
        irr_timeStop_refERT + 25000 * 2,         # (~7 hours intervals)
        irr_timeStop_refERT + 25000 * 3,
        irr_timeStop_refERT + 25000 * 4,
        irr_timeStop_refERT + 25000 * 5,
        irr_timeStop_refERT + 25000 * 6,
        irr_timeStop_refERT + 25000 * 7,
        irr_timeStop_refERT + 86400 * 2          # 2 days after
    ]))
    
    # Create water flux array
    if atmbc_dim == '2d':
        # For 2D: apply flux only to nodes within irrigation zone
        surfnodes, within_boundaries = get_nodes_irrigated(simuInfiltration)
        v_atmbc = np.zeros((len(atmbc_times), len(surfnodes)))
        
        # Apply infiltration rate during irrigation period
        for i, atmbci in enumerate(atmbc_times):
            if (atmbci >= irr_timeStart_refERT) & (atmbci < irr_timeStop_refERT):
                print(f'Adding infiltration at time {atmbci} s ({atmbci/60:.1f} min)')
                print(f'Cell size: {cell_size:.4f} m²')
                print(f'Infiltration rate: {Q_infiltration_m:.2e} m/s')
                v_atmbc[i, within_boundaries] = Q_infiltration_m
                
    else:
        # For 3D: simplified uniform application
        v_atmbc = np.zeros(len(atmbc_times))
        v_atmbc[(atmbc_times >= irr_timeStart_refERT) & 
                (atmbc_times < irr_timeStop_refERT)] = Q_infiltration_m
    
    return (
        atmbc_times, 
        v_atmbc, 
        elapsed_time_ERT_s, 
        irr_timeStart_refERT, 
        irr_timeStop_refERT
    )


# =============================================================================
# SENSOR POSITIONS AND OBSERVATIONS
# =============================================================================

def get_sensors_pos_ERT3d(maxDEM=1e-3):
    """
    Define sensor positions for 3D ERT monitoring.
    
    Returns coordinates of measurement points at three depths below the surface.
    These correspond to where physical sensors were installed in the field.
    
    Parameters
    ----------
    maxDEM : float, optional
        Surface elevation in meters (default: 0.001 m = 1 mm)
    
    Returns
    -------
    sensors_pos : dict
        Dictionary with sensor names as keys and [x, y, z] coordinates as values
        All sensors are at center of domain (x=0.5, y=0.5)
        Z-coordinates are depths below surface: -0.15 m, -0.30 m, -0.50 m
    
    Sensor Layout
    -------------
    - All sensors at horizontal center of domain (0.5 m, 0.5 m)
    - ERT_z1_15cm: 15 cm deep (shallow, responds quickly)
    - ERT_z2_30cm: 30 cm deep (intermediate depth)
    - ERT_z3_50cm: 50 cm deep (deeper, responds more slowly)
    
    Notes
    -----
    Negative z-values indicate depth below surface (CATHY convention)
    """
    sensors_pos = {
        'ERT_z1_15cm': [0.5, 0.5, maxDEM - 0.15],  # 15 cm depth
        'ERT_z2_30cm': [0.5, 0.5, maxDEM - 0.30],  # 30 cm depth
        'ERT_z3_50cm': [0.5, 0.5, maxDEM - 0.50],  # 50 cm depth
    }
    
    return sensors_pos


# =============================================================================
# SOIL PROPERTY CONFIGURATION
# =============================================================================

def prepare_soil_infiltration(DA_dict, simuInfiltration, filePath, 
                             porosity=0.5, soilHom=1):
    """
    Configure soil hydraulic properties for the infiltration simulation.
    
    Soil properties control how water moves through the soil. This function:
    1. Reads measured soil composition (sand, silt, clay percentages)
    2. Predicts hydraulic properties using pedotransfer functions
    3. Sets up vertical soil layers with appropriate properties
    
    Parameters
    ----------
    DA_dict : dict
        Data assimilation configuration dictionary with 'comp_level' key
    simuInfiltration : CATHY object
        Simulation object to configure
    filePath : str
        Path to soil property CSV files
    porosity : float, optional
        Soil porosity (fraction of volume that is pore space, default: 0.5)
    soilHom : int, optional
        Soil homogeneity flag (1 = homogeneous, 0 = layered, default: 1)
    
    Returns
    -------
    ordered_SPP_map : dict
        Soil Physical Properties map with hydraulic parameters for each layer
        Keys: PERMX, PERMY, PERMZ (hydraulic conductivity in each direction)
              ELSTOR (specific storage), POROS (porosity)
              VGNCELL, VGRMCCELL, VGPSATCELL (Van Genuchten parameters)
    FP_map : dict
        Feddes Parameters map for root water uptake
    layers : dict
        Layer depth boundaries (e.g., {1: [0, 0.4], 2: [0.4, 0.7]})
    zone3d_flag : list
        3D zone flags for each layer
    
    Hydraulic Parameters Explained
    ------------------------------
    - PERMX/Y/Z: How easily water flows in each direction (m/s)
      Higher values = water moves faster
    - POROS: Fraction of soil that is pore space (0-1)
      0.5 = 50% pores, 50% solid
    - VGNCELL: Van Genuchten n parameter (soil pore size distribution)
    - VGRMCCELL: Residual moisture content (water left in very dry soil)
    - VGPSATCELL: Air entry pressure (suction needed to start draining)
    
    Notes
    -----
    Van Genuchten parameters describe the soil water retention curve:
    how much water soil holds at different suction levels.
    """
    print('='*70)
    print('CONFIGURING SOIL PROPERTIES')
    print('='*70)
    
    # Read soil composition data (sand, silt, clay percentages)
    soil = pd.read_csv(filePath + 'soil.csv', sep=',')
    soil_selec = soil[soil['id'] == DA_dict['comp_level']]
    
    # Extract depth intervals and soil composition for each layer
    depths_soil = []
    VGP_predict_CATHY = []
    
    for i in range(len(soil_selec['depth(m)']) - 1):
        # Define depth interval for this layer
        if i == 0:
            depths_soil.append([0, soil_selec['depth(m)'].iloc[0]])
        else:
            depths_soil.append([
                soil_selec['depth(m)'].iloc[i],
                soil_selec['depth(m)'].iloc[i+1]
            ])
        
        # Get soil composition for this layer
        granulo = soil_selec.iloc[i][[
            'sand(%)', 
            'silt(%)', 
            'clay(%)', 
            'bd(g/cm^3)'  # Bulk density
        ]].to_numpy()
        
        print(f'Layer {i+1}: Sand={granulo[0]:.1f}%, '
              f'Silt={granulo[1]:.1f}%, Clay={granulo[2]:.1f}%')
        
        # Predict hydraulic parameters from soil texture
        # Uses pedotransfer functions (empirical relationships)
        VGP_predict_CATHY.append(
            predict_unsat_soil_hydr_param(data=[granulo])
        )
    
    # Order parameters correctly for CATHY
    orderedkeys_SPPmap = [
        'PERMX', 'PERMY', 'PERMZ',  # Hydraulic conductivity (m/s)
        'ELSTOR',                    # Specific storage (1/m)
        'POROS',                     # Porosity (-)
        'VGNCELL',                   # Van Genuchten n (-)
        'VGRMCCELL',                 # Residual moisture content (-)
        'VGPSATCELL'                 # Air entry pressure (m)
    ]
    
    # Organize predicted parameters
    SPP_map = {}
    for d in VGP_predict_CATHY:
        for key, value in d.items():
            if key not in SPP_map:
                SPP_map[key] = []
            SPP_map[key].append(value)
    
    ordered_SPP_map = {key: SPP_map[key] for key in orderedkeys_SPPmap}
    
    # Calculate mean values across all layers
    ordered_SPP_map_mean = {
        key: sum(values) / len(values) 
        for key, values in ordered_SPP_map.items()
    }
    
    # Create layer dictionary
    layers = {}
    for i, interval in enumerate(depths_soil, start=1):
        layers[i] = interval
    
    # Read Van Genuchten Parameters from measurements (if available)
    VGP = pd.read_csv(filePath + 'VGP_simple.csv', sep=',')
    VGP_selec = VGP[VGP['id'] == DA_dict['comp_level']]
    
    # Read current zone configuration
    rasterZone, headerZone = simuInfiltration.read_inputs('zone')
    
    # Configure soil properties
    if soilHom == 1:
        print('\nUsing HOMOGENEOUS soil properties (same for all layers)')
        
        nstr = simuInfiltration.dem_parameters['nstr']  # Number of layers
        
        # Use fixed, typical soil properties for all layers
        ordered_SPP_map = {
            "PERMX": [0.000188] * nstr,      # ~0.2 mm/s hydraulic conductivity
            "PERMY": [0.000188] * nstr,
            "PERMZ": [0.000188] * nstr,
            "ELSTOR": [1e-05] * nstr,        # Specific storage
            "POROS": [porosity] * nstr,      # Porosity (default 0.5)
            "VGNCELL": [1.46] * nstr,        # VG n parameter
            "VGRMCCELL": [0.15] * nstr,      # Residual moisture
            "VGPSATCELL": [0.03125] * nstr,  # Air entry pressure
        }
        
        print(f'Applied to {nstr} layers')
        print(f'Porosity: {porosity}')
        print(f'Hydraulic conductivity: {0.000188} m/s = {0.000188*100:.2f} cm/s')
    
    # Create uniform zone (all cells same soil type)
    zone_default = np.ones(np.shape(simuInfiltration.zone))
    zone_new = zone_default
    
    # Update zone in simulation
    simuInfiltration.update_zone(zone_new)
    simuInfiltration.update_soil()
    
    # Configure Feddes parameters (root water uptake)
    FP_map = simuInfiltration.soil_FP['FP_map']
    FP_map['ZROOT'] = [1e-3]  # Very shallow roots (1 mm) - minimal root uptake
    
    # Create 3D zone flag for each layer
    zone3d_flag = [zone_new] * nstr
    
    print('='*70)
    
    return ordered_SPP_map, FP_map, layers, zone3d_flag


# =============================================================================
# ERT CONFIGURATION AND DATA IMPORT
# =============================================================================

def pathName_3dERT(rootfilePath, compaction_level, synthObs, simuRef=None):
    """
    Get file paths and configuration for 3D ERT data.
    
    This function determines where ERT measurement files are located,
    whether they're synthetic (simulated) or real field measurements.
    
    Parameters
    ----------
    rootfilePath : str
        Root directory for data files
    compaction_level : str
        Soil compaction level identifier
    synthObs : bool
        True = use synthetic observations (from simulation)
        False = use real field measurements
    simuRef : CATHY object, optional
        Reference simulation (required if synthObs=True)
    
    Returns
    -------
    pathERT : str
        Path to ERT data files
    meshERT_pv : str
        Path to ERT mesh (VTK format)
    sequenceERT : str
        Path to ERT measurement sequence file
    data_format : str
        Data format identifier ('pygimli')
    elecs : np.array
        Electrode positions (x, y, z coordinates)
    
    ERT Measurement Sequence
    ------------------------
    The sequence file defines:
    - Which electrodes inject current (A and B electrodes)
    - Which electrodes measure voltage (M and N electrodes)
    - This is called a "four-point measurement" or "Wenner array"
    
    Example Measurement
    -------------------
    Electrode 1-2: inject current
    Electrode 5-6: measure voltage
    From this, calculate apparent resistivity (Ω⋅m)
    """
    _, comp_level_st = set_acronyms(compaction_level)
    
    if synthObs:
        # Use synthetic data from simulation
        print('Using SYNTHETIC ERT observations')
        pathERT = os.path.join(
            simuRef.workdir,
            simuRef.project_name,
            simuRef.project_name
        )
        sequenceERT = rootfilePath + '/assets/seq3d_from_st1t00tldat.shm'
        
    else:
        # Use real field measurements
        print('Using REAL field ERT observations')
        pathERT = rootfilePath + 'data/3d/' + ''.join(comp_level_st[compaction_level]) + '/procdata'
        
        # Sequence file varies by compaction level
        sequenceERT = rootfilePath + f'data/3d/seq3d_from_{comp_level_st[compaction_level][0]}t00tldat.shm'
        
        if compaction_level == 'light_comp':
            sequenceERT = rootfilePath + f'data/3d/seq3d_from_{comp_level_st[compaction_level][0]}t01tldat.shm'
        if compaction_level == 'no_comp':
            sequenceERT = rootfilePath + f'data/3d/seq3d_from_{comp_level_st[compaction_level][0]}t00tldat_valid.shm'
    
    # ERT mesh file (smaller mesh for efficiency)
    meshERT_pv = rootfilePath + '/assets/mesh_Luca_Agripolis3d/meshSmall.vtk'
    
    # Load sequence data to get electrode positions
    dataERT = pg.load(sequenceERT)
    data_format = 'pygimli'
    elecs = np.array(dataERT.sensors())  # Electrode coordinates
    
    print(f'ERT mesh: {meshERT_pv}')
    print(f'Number of electrodes: {len(elecs)}')
    
    return pathERT, meshERT_pv, sequenceERT, data_format, elecs


def pathName_2dERT(rootfilePath, compaction_level, synthObs, simuRef):
    """
    Get file paths and configuration for 2D ERT data.
    
    Similar to pathName_3dERT but for 2D cross-section measurements.
    
    Parameters
    ----------
    rootfilePath : str
        Root directory for data files
    compaction_level : str
        Soil compaction level identifier
    synthObs : bool
        True = synthetic, False = real data
    simuRef : CATHY object
        Reference simulation
    
    Returns
    -------
    pathERT : str
        Path to ERT data files
    meshERT : str
        Path to 2D ERT mesh (VTK format)
    sequenceERT : str
        Path to ERT measurement sequence
    data_format : str
        Data format identifier
    elecs : np.array
        Electrode positions
    """
    _, comp_level_st = set_acronyms(compaction_level)
    
    if synthObs:
        pathERT = os.path.join(
            simuRef.workdir,
            simuRef.project_name,
            simuRef.project_name
        )
    else:
        pathERT = rootfilePath + 'data/2d/' + ''.join(comp_level_st[compaction_level]) + '/procdata'
    
    # 2D sequence and mesh
    sequenceERT = rootfilePath + 'data/2d/TL_ERT_2d.shm'
    meshERT = rootfilePath + 'data/2d/mesh_k.vtk'
    
    # Load sequence
    dataERT = pg.load(sequenceERT)
    data_format = 'pygimli'
    elecs = np.array(dataERT.sensors())
    
    return pathERT, meshERT, sequenceERT, data_format, elecs


# =============================================================================
# MESH ALIGNMENT AND TRANSFORMATION
# =============================================================================

def adjust_mesh_positionning_fwd_mesh(simuInf3D_DA, mesh_OUT):
    """
    Align ERT mesh with CATHY hydrological mesh.
    
    The ERT mesh and CATHY mesh may have different:
    - Coordinate systems (origins, axes)
    - Spatial positioning
    
    This function transforms the ERT mesh to match CATHY's coordinate system
    so that we can properly transfer water content from CATHY to ERT.
    
    Parameters
    ----------
    simuInf3D_DA : CATHY object
        CATHY simulation object
    mesh_OUT : str or pv.UnstructuredGrid
        Path to ERT mesh file or PyVista mesh object
    
    Returns
    -------
    mesh_PG_new : pv.UnstructuredGrid
        Transformed ERT mesh aligned with CATHY
    mesh_PG_new.points : np.array
        Transformed node coordinates
    
    Transformations Applied
    -----------------------
    1. Read ERT mesh
    2. Apply transformation matrix (currently identity - no rotation)
    3. Translate in x-direction: +0.2 m
    4. Translate in y-direction: +0.3 m
    5. Keep z-coordinates unchanged
    
    Why This Matters
    ----------------
    For accurate ERT forward modeling, we need to:
    1. Extract water saturation from CATHY mesh nodes
    2. Map these values to ERT mesh nodes
    3. Calculate electrical resistivity at ERT nodes
    4. Simulate ERT measurements
    
    If meshes aren't aligned, we'll map water content to wrong locations!
    """
    # Read ERT mesh
    mesh_PG = pv.read(mesh_OUT)
    points = mesh_PG.points.copy()
    mesh_PG_new = mesh_PG.copy()
    
    # Define transformation matrix
    # Currently identity (no rotation), but could be modified for rotations
    transform_matrix = np.array([
        [1, 0, 0],  # X-axis: no change
        [0, 1, 0],  # Y-axis: no change
        [0, 0, 1]   # Z-axis: no change
    ])
    
    # Apply transformation
    transformed_points_mesh_PG = np.dot(points, transform_matrix)
    
    # No change in z (depth) coordinate
    transformed_points_mesh_PG[:, 2] = transformed_points_mesh_PG[:, 2]
    
    # Translate in x and y to align with CATHY mesh
    transformed_points_mesh_PG[:, 0] = transformed_points_mesh_PG[:, 0] + 0.2  # Shift +20 cm in x
    transformed_points_mesh_PG[:, 1] = transformed_points_mesh_PG[:, 1] + 0.3  # Shift +30 cm in y
    
    # Update mesh with transformed coordinates
    transformed_points_mesh_PG.points = transformed_points_mesh_PG 
    mesh_PG_new.points = transformed_points_mesh_PG.points
    
    # Verify alignment with CATHY mesh
    mesh2read = os.path.join(
        simuInf3D_DA.workdir, 
        simuInf3D_DA.project_name, 
        'vtk', 
        simuInf3D_DA.project_name + '.vtk'
    )
    
    print('Mesh alignment complete')
    print(f'ERT mesh translated by: dx=+0.2 m, dy=+0.3 m')
    
    return mesh_PG_new, mesh_PG_new.points


def adjust_mesh_positionning(mesh_OUT):
    """
    Alternative mesh positioning adjustment (swaps Y and Z axes).
    
    Used for 2D meshes where CATHY and PyGIMLi have different axis conventions.
    
    Parameters
    ----------
    mesh_OUT : str or pv.UnstructuredGrid
        Path to mesh file or PyVista mesh
    
    Returns
    -------
    mesh_PG_new : pv.UnstructuredGrid
        Transformed mesh with swapped axes
    transformed_points_mesh_PG.points : np.array
        Transformed coordinates
    
    Transformation
    --------------
    - X-axis: unchanged
    - Y-axis: becomes Z-axis
    - Z-axis: becomes Y-axis, set to constant value (1.0 m)
    
    This handles different axis conventions:
    - CATHY: X (horizontal), Y (horizontal), Z (vertical/depth)
    - PyGIMLi: X (horizontal), Y (vertical/depth), Z (horizontal)
    """
    mesh_PG = pv.read(mesh_OUT)
    points = mesh_PG.points.copy()
    mesh_PG_new = mesh_PG.copy()
    
    # Transformation matrix swaps Y and Z
    transform_matrix = np.array([
        [1, 0, 0],  # X stays X
        [0, 0, 1],  # Y becomes Z
        [0, 1, 0]   # Z becomes Y
    ])
    
    transformed_points_mesh_PG = np.dot(points, transform_matrix)
    
    # Set Y to constant (mesh is essentially 2D in X-Z plane)
    transformed_points_mesh_PG[:, 1] = 1
    
    # Adjust Z (depth) by adding 0.1 m
    transformed_points_mesh_PG[:, 2] = transformed_points_mesh_PG[:, 2] + 0.1
    
    # Update mesh
    transformed_points_mesh_PG.points = transformed_points_mesh_PG
    mesh_PG_new.points = transformed_points_mesh_PG.points
    
    return mesh_PG_new, transformed_points_mesh_PG.points


# =============================================================================
# ERT FORWARD MODELING - Converting Water Content to Electrical Measurements
# =============================================================================

def ERT_fwd_modelling_ref(simu, sc, time_of_interest, path_inputs_rhizo,
                          pathERT_results, sequenceERT_pg, meshERT_pg, 
                          sequenceERT):
    """
    Perform ERT forward modeling to generate synthetic measurements.
    
    This is the core function that converts simulated water content into
    what ERT measurements would look like. The process:
    
    1. Read water saturation from CATHY simulation
    2. Convert saturation to electrical resistivity using Archie's Law
    3. Map resistivity onto ERT mesh
    4. Simulate ERT measurements using PyGIMLi
    5. Add realistic measurement noise
    6. Save synthetic ERT data files
    
    Parameters
    ----------
    simu : CATHY object
        Completed CATHY simulation with water content results
    sc : dict
        Configuration dictionary with Archie's Law parameters:
        - rFluid_Archie: fluid electrical resistivity (Ω⋅m)
        - POROS: soil porosity (-)
        - a_Archie: tortuosity factor (-)
        - m_Archie: cementation exponent (-)
        - n_Archie: saturation exponent (-)
        - noise_level_ERT: measurement noise (%)
    time_of_interest : list
        Time points to generate ERT measurements
    path_inputs_rhizo : str
        Path to CATHY simulation inputs
    pathERT_results : str
        Path to save ERT results
    sequenceERT_pg : str
        ERT measurement sequence file
    meshERT_pg : pv.UnstructuredGrid
        ERT mesh (PyGIMLi format)
    sequenceERT : str
        ERT sequence file path
    
    Returns
    -------
    None (saves files to disk)
    
    Output Files
    ------------
    - Archie.png: Plot showing relationship between saturation and resistivity
    - Archie.csv: Saturation-resistivity data
    - ER_predicted_sol_t*.data: Synthetic ERT measurements at each time
    
    Archie's Law Explained
    ----------------------
    Archie's Law relates electrical resistivity (ρ) to water saturation (S):
    
    ρ = a * φ^(-m) * S^(-n) * ρ_fluid
    
    Where:
    - ρ = bulk soil resistivity (what ERT measures)
    - a = tortuosity factor (~0.5-2.0, accounts for pore geometry)
    - φ = porosity (fraction of soil that is pores)
    - m = cementation exponent (~1.3-2.5, describes pore connectivity)
    - S = water saturation (0 = dry, 1 = fully saturated)
    - n = saturation exponent (~1.5-2.5, describes water distribution)
    - ρ_fluid = electrical resistivity of pore water
    
    Physical Interpretation:
    - Wet soil (high S) → low resistivity → conducts electricity well
    - Dry soil (low S) → high resistivity → insulates electricity
    - Sandy soil (high φ, low m) → more sensitive to saturation changes
    - Clay soil (low φ, high m) → less sensitive to saturation changes
    
    Example
    -------
    If saturation changes from 0.5 (half-saturated) to 1.0 (fully saturated):
    - With n=2: resistivity decreases by factor of 4
    - ERT can detect this change and map water infiltration
    """
    print('\n' + '='*70)
    print('ERT FORWARD MODELING - Converting water content to ERT measurements')
    print('='*70)
    
    # Extract Archie's Law parameters
    rFluid_Archie = sc['rFluid_Archie']
    POROS = sc['POROS']
    
    # Read water saturation results from CATHY simulation
    df_sw, _ = simu.read_outputs(filename='sw')
    print(f'Read saturation data: {df_sw.shape[0]} time steps, {df_sw.shape[1]} nodes')
    
    # Create output directory if needed
    output_dir = os.path.join(pathERT_results, simu.project_name)
    isExist = os.path.exists(output_dir)
    if not isExist:
        os.makedirs(output_dir)
        print(f'Created output directory: {output_dir}')
    
    # =========================================================================
    # STEP 1: Convert saturation to electrical resistivity (Archie's Law)
    # =========================================================================
    
    print('\nStep 1: Applying Archie\'s Law')
    print(f'Porosity: {POROS}')
    print(f'Fluid resistivity: {rFluid_Archie[0]:.2f} Ω⋅m')
    print(f'Archie parameters: a={sc["a_Archie"][0]}, m={sc["m_Archie"][0]}, n={sc["n_Archie"][0]}')
    
    # Convert WITHOUT noise (for reference)
    ER_converted_ti_noNoise = Archie.Archie_rho_DA(
        rFluid_Archie=rFluid_Archie, 
        sat=[df_sw.to_numpy()],
        porosity=[POROS], 
        a_Archie=sc['a_Archie'], 
        m_Archie=sc['m_Archie'], 
        n_Archie=sc['n_Archie'],
        pert_sigma_Archie=[0]  # No noise
    )
    
    # Convert WITH noise (realistic measurements)
    ER_converted_ti = Archie.Archie_rho_DA(
        rFluid_Archie=rFluid_Archie, 
        sat=[df_sw.to_numpy()],
        porosity=[POROS], 
        a_Archie=sc['a_Archie'], 
        m_Archie=sc['m_Archie'], 
        n_Archie=sc['n_Archie'],
        pert_sigma_Archie=sc['pert_sigma_Archie']  # With noise
    )
    
    # Organize results into DataFrame
    Archie_simu_sol = pd.DataFrame([
        np.hstack(df_sw.to_numpy()),              # Saturation
        np.hstack(ER_converted_ti_noNoise[:]),    # Resistivity (no noise)
        np.hstack(ER_converted_ti[:])             # Resistivity (with noise)
    ]).transpose()
    
    Archie_simu_sol.columns = ['sw', 'rho', 'rho_noise']
    
    # =========================================================================
    # STEP 2: Plot and save Archie relationship
    # =========================================================================
    
    print('\nStep 2: Visualizing Archie\'s Law relationship')
    
    fig, ax = plt.subplots(figsize=(8, 6))
    Archie_simu_sol.plot.scatter(
        x='sw', y='rho_noise', 
        color='red', ax=ax, 
        label='With noise', alpha=0.5
    )
    Archie_simu_sol.plot.scatter(
        x='sw', y='rho',
        ax=ax, label='No noise', alpha=0.5
    )
    ax.set_xlabel('Water Saturation (-)', fontsize=12)
    ax.set_ylabel('Electrical Resistivity (Ω⋅m)', fontsize=12)
    ax.set_title('Archie\'s Law: Saturation vs Resistivity', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    save_path = os.path.join(output_dir, 'Archie.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f'Saved Archie plot: {save_path}')
    plt.close()
    
    # Save data
    csv_path = os.path.join(output_dir, 'Archie.csv')
    Archie_simu_sol.to_csv(csv_path)
    print(f'Saved Archie data: {csv_path}')
    
    # =========================================================================
    # STEP 3: Generate ERT measurements for each time step
    # =========================================================================
    
    print(f'\nStep 3: Generating ERT measurements for {len(ER_converted_ti)} time steps')
    
    for i in range(len(ER_converted_ti)):
        
        print(f'\nProcessing time step {i+1}/{len(ER_converted_ti)}...')
        
        filenameSolRoot = 'ER_predicted_sol_t' + str(i)
        
        # ---------------------------------------------------------------------
        # 3a. Read CATHY mesh data for this time step
        # ---------------------------------------------------------------------
        
        path_CATHY = os.path.join(simu.workdir, simu.project_name, 'vtk/')
        
        # Construct filename based on time index
        if i < 10:    
            filename = "10" + str(i) + ".vtk"
        elif i < 100:
            filename = "1" + str(i) + ".vtk"
        elif i < 200:
            newnb = [int(x) for x in str(i)]
            filename = "2" + str(newnb[1]) + str(newnb[2]) + ".vtk"
        elif i < 300:
            newnb = [int(x) for x in str(i)]
            filename = "3" + str(newnb[1]) + str(newnb[2]) + ".vtk"
        
        print(f'  Resistivity range: {np.min(ER_converted_ti[i, :]):.1f} - '
              f'{np.max(ER_converted_ti[i, :]):.1f} Ω⋅m')
        
        # ---------------------------------------------------------------------
        # 3b. Add converted resistivity as attribute to CATHY mesh
        # ---------------------------------------------------------------------
        
        (mesh_CATHY_new_attr, 
         mesh_CATHY_new_attr_name) = mt.add_attribute_2mesh(
             ER_converted_ti[i, :],
             path_CATHY + filename,
             'ER_converted_CATHYmsh',
             overwrite=True, 
             path=path_CATHY,
             time=i
         )
        
        # ---------------------------------------------------------------------
        # 3c. Copy attribute from CATHY mesh to ERT mesh
        # ---------------------------------------------------------------------
        
        ERT_meta_dict = sc
        mesh_pg_new_attr, scalar_new = CATHY_2_pg(
            mesh_CATHY_new_attr,
            ERT_meta_dict,
            scalar=mesh_CATHY_new_attr_name,
            show=False,
            time=i, 
            path=os.path.join(simu.workdir, simu.project_name, 'vtk/'),
        )
        
        # ---------------------------------------------------------------------
        # 3d. Simulate ERT measurements using PyGIMLi
        # ---------------------------------------------------------------------
        
        # Extract resistivity values at ERT mesh nodes
        res0 = mesh_pg_new_attr.get_array(scalar_new)
        
        # Perform ERT forward modeling
        # This simulates what the ERT instrument would measure
        ERT_predicted = simuERT.create_ERT_survey_pg(
            os.path.join(
                pathERT_results, 
                simu.project_name + '_pg_FOP',
                'pg_ER_predicted_sol_t' + str(i)
            ), 
            sequenceERT_pg, 
            meshERT_pg, 
            res0=res0,
            noise_level=sc['noise_level_ERT']  # Add realistic noise
        )
        
        # ---------------------------------------------------------------------
        # 3e. Save ERT measurements to file
        # ---------------------------------------------------------------------
        
        # Organize measurement data
        d = {
            'a': ERT_predicted['a'],      # Current electrode A
            'b': ERT_predicted['b'],      # Current electrode B
            'k': ERT_predicted['k'],      # Geometric factor
            'm': ERT_predicted['m'],      # Voltage electrode M
            'n': ERT_predicted['n'],      # Voltage electrode N
            'rhoa': ERT_predicted['rhoa'],  # Apparent resistivity (Ω⋅m)
            'valid': ERT_predicted['valid']  # Valid measurement flag
        }
        df_2_save = pd.DataFrame(data=d)
        
        print(f'  Apparent resistivity range: {min(df_2_save["rhoa"]):.1f} - '
              f'{max(df_2_save["rhoa"]):.1f} Ω⋅m')
        
        # Save to file in PyGIMLi format
        filedataext = os.path.join(output_dir, filenameSolRoot + '.data')
        ERT_predicted.save(filedataext)
        print(f'  Saved: {filedataext}')
    
    print('\n' + '='*70)
    print('ERT FORWARD MODELING COMPLETE')
    print('='*70)
    print(f'Generated {len(ER_converted_ti)} synthetic ERT datasets')
    print(f'Results saved in: {output_dir}')


# =============================================================================
# UTILITY FUNCTIONS FOR DATA MANAGEMENT AND LOGGING
# =============================================================================

def backup_simulog_DA(DAinputs, logfile='simuDAlog.csv'):
    """
    Log data assimilation simulation parameters to CSV file.
    
    This creates a record of each simulation run with all its settings,
    allowing you to track what parameters were used and when.
    
    Parameters
    ----------
    DAinputs : dict
        Dictionary of data assimilation parameters
    logfile : str, optional
        Path to log file (default: 'simuDAlog.csv')
    
    Returns
    -------
    results_df : pd.DataFrame
        Updated log DataFrame with all simulation records
    matching_index : int
        Index of current simulation in log (for reference)
    
    Notes
    -----
    If simulation with same parameters already exists, updates timestamp
    instead of creating duplicate entry.
    """
    # Read existing log file
    results_df = pd.read_csv(logfile, index_col=0)
    now = datetime.now()
    
    # Process assimilated observations field
    if len(list(DAinputs['AssObs'])) > 1:   
        AssObs_new = '/'.join(list(DAinputs['AssObs']))
        DAinputs['AssObs'] = AssObs_new
    else:
        AssObs_new = re.search(r"'(.*?)'", str(DAinputs['AssObs'])).group(1)
        DAinputs['AssObs'] = AssObs_new
    
    # Create new entry
    results_df_new = pd.DataFrame([DAinputs])
    
    # Check if simulation with same parameters already exists
    cols2check = list(DAinputs.keys())
    values = results_df_new[cols2check].values
    print((results_df[cols2check] == values).T)
    
    matching_index = results_df.index[(results_df[cols2check] == values).all(axis=1)].tolist()
    
    if matching_index:
        # Update existing entry
        now = datetime.now()
        results_df.loc[matching_index, 'datetime'] = now
        matching_index = matching_index[0]
    else:
        # Add new entry
        results_df_new['datetime'] = now
        results_df = pd.concat([results_df, results_df_new], ignore_index=True)
        matching_index = len(results_df) - 1
    
    # Save updated log
    results_df.to_csv(logfile)
    
    return results_df, matching_index


# =============================================================================
# FILE INFORMATION
# =============================================================================

print('\n' + '='*70)
print('INFILTRATION UTILITY MODULE LOADED')
print('='*70)
print('This module provides utilities for:')
print('  - Setting up infiltration experiments')
print('  - Creating computational meshes')
print('  - Managing ERT measurements')
print('  - Forward modeling ERT responses')
print('  - Data visualization and analysis')
print('='*70 + '\n')