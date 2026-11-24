"""Where does the disk space of a DA run go, and how to shrink it.

    python lighten_run.py report  [RUN_DIR ...]      # default: every run folder of [outputs].simulations
    python lighten_run.py convert FILE.pkl [--keep]   # legacy long df_DA pickle -> compact parquet (float32, zstd)
    python lighten_run.py prune   RUN_DIR [--apply]   # list (or delete) what nb8 does not need

What is heavy (see README, "Lighter DA folder"):
  * member CATHY output, mostly .vtk, in <run>/DA_Ensemble/   -> nb7 --cleanup-vtk all (default) removes it at the end
  * the legacy long table <run>_df_DA.pkl (float64, node ids repeated per row)
  * per-cycle chunks <run>_states/cycle_*.npz and <run>_states.nc are the compact store that nb8 reads first
`prune` only touches the legacy table and only when the compact store exists (nb8 loads `states` first and never
needs the long table then); it never touches the real-data runs of nb5/nb6 (they have no compact store).
Without --apply nothing is deleted. The parquet is written next to the pickle; the pickle is kept unless --delete.
"""
import argparse
import glob
import os
import sys

CATEGORIES = [('vtk', ('.vtk',)), ('pkl', ('.pkl', '.pickle')), ('npz', ('.npz',)), ('netcdf', ('.nc',)),
              ('parquet', ('.parquet',)), ('cathy-text', ('.txt', '.dat', '.csv', '.json', '.png'))]


def _cat(fn):
    ext = os.path.splitext(fn)[1].lower()
    for name, exts in CATEGORIES:
        if ext in exts:
            return name
    return 'other'


def fmt(n):
    for u in ('B', 'kB', 'MB', 'GB'):
        if n < 1000 or u == 'GB':
            return f"{n:.0f} {u}" if u == 'B' else f"{n:.1f} {u}"
        n /= 1000


def scan(run_dir):
    """{category: bytes}, [(bytes, path)] of one run folder (and its sibling <run>_meta)."""
    cats, files = {}, []
    for base in (run_dir, run_dir + '_meta'):
        for r, _, fs in os.walk(base):
            for f in fs:
                p = os.path.join(r, f)
                try:
                    n = os.path.getsize(p)
                except OSError:
                    continue
                cats[_cat(f)] = cats.get(_cat(f), 0) + n
                files.append((n, p))
    return cats, sorted(files, reverse=True)


def has_compact_store(run_dir):
    name = os.path.basename(os.path.normpath(run_dir))
    return bool(glob.glob(os.path.join(run_dir, name + '_states', 'cycle_*.npz'))) or os.path.isfile(
        os.path.join(run_dir, name + '_states.nc'))


def legacy_tables(run_dir):
    return sorted(glob.glob(os.path.join(run_dir, '*_df_DA.pkl')))


def report(run_dirs, top=5):
    grand = 0
    for d in run_dirs:
        cats, files = scan(d)
        tot = sum(cats.values())
        grand += tot
        print(f"\n{os.path.basename(os.path.normpath(d))}: {fmt(tot)}"
              f"  [{'compact store' if has_compact_store(d) else 'no compact store'}]")
        for c, n in sorted(cats.items(), key=lambda kv: -kv[1]):
            print(f"    {c:11s} {fmt(n):>10s}  {100 * n / max(tot, 1):4.0f} %")
        for n, p in files[:top]:
            print(f"      {fmt(n):>10s}  {os.path.relpath(p, d)}")
        if cats.get('vtk', 0) > 0.3 * tot: 
