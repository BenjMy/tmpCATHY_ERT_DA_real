"""Real acquisition times of ERT surveys -> CATHY times (generic, no pygimli/CATHY import).

Replaces the hand-typed placeholder table ERT_TIMES_H of nb5. Same time convention as
nb3_compare_meas_vs_mod.py:  site_params.json holds INFILTRATION_START per experiment
('HC', 'NC', ...) and CATHY t=0 is irr_start_hour BEFORE it.

    t_h   = (acquisition - INFILTRATION_START) [h]          (negative = before irrigation)
    t_cat = t_h * HOUR + irr_start_hour * HOUR [s]          (what load_ert_files expects)

Acquisition times file (csv, one row per survey, header required):
    survey,datetime
    0,2025-05-12 08:40
    1,2025-05-12 13:05
  'survey' = integer index of st5t<ii>.tldat; or give 'file' (st5t03.tldat) instead and the
  index is parsed from the name. Day-first dates by default (as the sensor files in nb3).
"""
import json
import re
import numpy as np
import pandas as pd

COMP_TO_LABEL = {'heavy_comp': 'HC', 'no_comp': 'NC'}     # same mapping as nb3 (+ extend here)


def load_infiltration_start(site_params_json, comp_level):
    label = COMP_TO_LABEL.get(comp_level)
    if label is None:
        raise KeyError(f"no experiment label for comp_level={comp_level!r}; add it to COMP_TO_LABEL")
    with open(site_params_json) as f:
        sp = json.load(f)
    return pd.Timestamp(sp['INFILTRATION_START'][label])


def survey_index_from_name(name):
    """'st5t03.tldat' -> 3 (last integer group of the stem)."""
    stem = re.sub(r'\.[^.]*$', '', str(name).replace('\\', '/').split('/')[-1])
    m = re.findall(r'(\d+)', stem)
    if not m:
        raise ValueError(f"cannot read a survey index from {name!r}")
    return int(m[-1])


def read_acquisition_times(csv_path, time_col='datetime', dayfirst=True):
    """csv -> {survey index: Timestamp}. Needs a 'survey' or a 'file' column."""
    df = pd.read_csv(csv_path)
    if 'survey' in df.columns:
        idx = df['survey'].astype(int)
    elif 'file' in df.columns:
        idx = df['file'].map(survey_index_from_name)
    else:
        raise ValueError(f"{csv_path}: need a 'survey' or 'file' column, found {list(df.columns)}")
    if time_col not in df.columns:
        raise ValueError(f"{csv_path}: no column {time_col!r}, found {list(df.columns)}")
    ts = pd.to_datetime(df[time_col], dayfirst=dayfirst)
    if idx.duplicated().any():
        raise ValueError(f"{csv_path}: duplicated survey index {idx[idx.duplicated()].tolist()}")
    return {int(i): pd.Timestamp(t) for i, t in zip(idx, ts)}


SCHEDULE_COL = {'heavy_comp': '1-HC', 'light_comp': '3-LC', 'no_comp': '5-NC'}   # as infiltration_utils.set_acronyms


def read_schedule_2023(csv_path, comp_level, fmt='%d/%m/%Y %H:%M'):
    """schedule2023_3d.csv (same file as infiltration_utils.read_infiltration_log) -> {survey index: Timestamp}.

    Survey index = row index of the csv (st5t<ii>.tldat), column = '5-NC' / '3-LC' / '1-HC'.
    Rows with an empty/unparsable date (survey not acquired) are skipped.
    """
    col = SCHEDULE_COL[comp_level]
    df = pd.read_csv(csv_path)
    if col not in df.columns:
        raise ValueError(f"{csv_path}: no column {col!r}, found {list(df.columns)}")
    ts = pd.to_datetime(df[col], format=fmt, errors='coerce')
    return {int(i): pd.Timestamp(t) for i, t in ts.items() if pd.notna(t)}


def times_h_since_irrigation(acq, infiltration_start):
    """{survey: Timestamp} -> {survey: hours since irrigation start} (drop-in for ERT_TIMES_H)."""
    t0 = pd.Timestamp(infiltration_start)
    return {int(k): round((pd.Timestamp(v) - t0).total_seconds() / 3600.0, 6) for k, v in acq.items()}


def check_schedule(times_h, irr_start_hour, sim_end_h=None):
    """Return a list of problems (empty = fine). t_cat = t_h + irr_start_hour must lie in (0, end]."""
    problems = []
    keys = sorted(times_h)
    vals = np.array([times_h[k] for k in keys], float)
    if np.any(np.diff(vals) <= 0):
        problems.append("survey times are not strictly increasing with the survey index: "
                        + ", ".join(f"{k}:{v:.2f}h" for k, v in zip(keys, vals)))
    t_cat = vals + irr_start_hour
    for k, t in zip(keys, t_cat):
        if t <= 0:
            problems.append(f"survey {k} is at CATHY t={t:.2f} h <= 0: before the simulation start "
                            f"(the first model time cannot be assimilated)")
        if sim_end_h is not None and t > sim_end_h:
            problems.append(f"survey {k} is at CATHY t={t:.2f} h > simulation end {sim_end_h:.2f} h")
    return problems
