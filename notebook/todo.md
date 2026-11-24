# TODO — tmpCATHY_ERT_DA_real

Branch: `feature/data-paths-registry-readme`
Repo: https://github.com/BenjMy/tmpCATHY_ERT_DA_real/tree/feature/data-paths-registry-readme

Legend: `[x]` done · `[ ]` to do · *italic* = note / answer

**Working rule:** changes are committed as BenjMy and shared as `.patch` files (`git format-patch`).

## Contents

1. [Next up](#1-next-up)
2. [Repository hygiene](#2-repository-hygiene)
3. [Reproducibility and run tracking](#3-reproducibility-and-run-tracking)
4. [Documentation (README)](#4-documentation-readme)
5. [Scientific checks on the DA runs](#5-scientific-checks-on-the-da-runs)
6. [Warnings to check](#6-warnings-to-check)
7. [Lighter DA runs (protocol)](#7-lighter-da-runs-protocol)
8. [Reorganise the files](#8-reorganise-the-files)
9. [Branches (pycathy_wrapper)](#9-branches-pycathy_wrapper)

---

## 1. Next up

- [x] Make this todo friendlier (sections, contents, legend)
- [x] `truth_NC.pkl` was tracked although `*.pkl` is in `.gitignore` → untracked (`git rm --cached`, file kept on disk)
- [x] Improve `.gitignore`: folders ending with `_meta` and `*.log` files (e.g. `SYN_NC_full_upd_pf_N32_s1_unbiased_issues.log`, see `tmp_treedir_typical.md`)
- [ ] Run the lighter version: which command? (see section 7)
- [ ] Check filter inbreeding in `SYN_NC_full_upd_pf_N32_s1_unbiased_issues` (see section 5)

## 2. Repository hygiene

- [x] `.gitignore`: `simulations*`, `*_meta/`, `*.log`, `*.pkl`
- [x] Untrack `data/3d/synthetic/NC/truth_NC.pkl`
- [ ] Still tracked although now ignored (a `.gitignore` never untracks): decide whether to `git rm --cached` them
  - `notebook/cache/*.pkl` (6 files, `prior_NC_*`)
  - `notebook/logs/*.log` (24 files, `SYN_NC_*`)
- [ ] Untrack `__pycache__/` and `.ipynb_checkpoints/` (40 tracked files) and add them to `.gitignore`
- [ ] Decide where the large tracked files go: `notebook/logs`, `data/**/duckdb.exe` (3 × 24 MB)

## 3. Reproducibility and run tracking

- [x] Manifest file for each run: version of CATHY DA ERT (with branch), CATHY and pyCATHY (with branch), and other useful reproducibility info → `run_manifest.py`
- [x] Register all the runs in a CSV file to track what has been done or not → `registry/runs_registry.csv`
- [ ] `data_path.toml` to control all the paths: inputs, outputs, plots, code lib, etc.

## 4. Documentation (README)

- [ ] Create a new branch on GitHub and add a README explaining the whole code architecture: input types, how to run, scenario construction, output data and plots
- [ ] Explain the synthetic ERT DA pipeline
- [ ] Explain how to use the `scenario.toml` file
- [ ] Add a "Lighter DA folder" section (`ligher_run.py` docstring already points to it)
- [ ] Add a table old name → new name once files are renamed (section 8)

## 5. Scientific checks on the DA runs

- [ ] Check filter inbreeding in `SYN_NC_full_upd_pf_N32_s1_unbiased_issues`
- [ ] Are the "unbiased" cases really unbiased?
  - Example: `SYN_NC_full_upd_pf_N32_s1_unbiased`
  - `..._profiles.png`: the truth is inside the ensemble, but the depth profile of (ensemble mean − truth) is not zero/consistent

## 6. Warnings to check

- [ ] `warnings.warn("truth forward on the SAME mesh as the DA: mesh discretisation error is not tested ...")`
- [ ] ERT frame: legacy (`farming_DA_ERT3d.py`); `signs=(1, 1, 1)`, `shift=[2.0, 2.0, 0.001]`
- [ ] pyGIMLi forward mesh margin (raw message):

```text
pycathy_wrapper/pyCATHY/ERT/alignment.py:116: UserWarning: [ERT] pg forward mesh extends only 0.50 m
beyond the electrodes (array span 1.10 m): boundary effects can bias the forward apparent resistivities
of weak-signal quadrupoles (see nc_model.unreliable_quadrupoles: homogeneous test, rhoa/rho must be 1)
```

## 7. Lighter DA runs (protocol)

*Investigation only, nothing implemented. Written 2026-11-24 from this repo + `synthetic_twin_plan.md` s.16 + a read of `pycathy_wrapper@main` (`pyCATHY/DA/cathy_DA.py`, `cathy_tools.py`). Most of the weight is on the pyCATHY side (item P).*

**Already there**
- Compact float32 `StateStore` (`_states/cycle_*.npz`): long `df_DA` 314 MB → ~35–65 MB
- `--store-full analysis`, `--store-stride`
- `--cleanup-vtk` at the END of a run
- `ligher_run.py report|convert|prune`

**Known gaps**
- README has no "Lighter DA folder" section although the `ligher_run.py` docstring points to it
- `_df_DA.pkl` is still written once at the end
- `cleanup_vtk` acts only at the end: peak disk = NENS × ~38 MB × …, vtk is ~31 of 38 MB per member
- Plan s.16 "not done": `VTKF=0` for the truth run (1.4 GB vtk, untested), `sw` as uint16

### 7.0 Baseline (one reference run, never changed)

- [ ] NC, scenario `ic_Ks_upd`, EnKF, N=8, seed 1 (22 surveys, ~25 min on 2 cores), plus the same on HC with N=8
- Record:
  - `du -sh` per category (`ligher_run.py report`)
  - peak RSS (`/usr/bin/time -v`)
  - wall time split: CATHY windows / ERT forward (pyGIMLi) / analysis (add timers)
  - the `manifest.json` of the run
- Metrics to protect: nb8 `summary.csv` (RMSE, spread, 90 % coverage, CRPS) and `parameters_summary.csv` (z prior → final)

### 7.1 Disk (no change of the result expected: states bit-identical or < 1e-6 relative)

- [ ] a. `VTKF=0` for members and truth: does `read_outputs` (sw, psi) need the vtk? If not, never write them
- [ ] b. Else delete each member's vtk after every cycle (not only at the end): peak disk becomes O(1 cycle)
- [ ] c. Per-member folder: split size into input files (mesh, dem, soil: identical for all members → symlink/hardlink instead of NENS copies?) vs outputs; keep restart files only
- [ ] d. Do not write `<run>_df_DA.pkl` when the compact store exists (flag, default off); check parameters pkl size
- [ ] e. dtype: psi float32, sw uint16 (22 MB vs 314 MB) or float16; error must stay << 5 % ERT noise
- [ ] f. Printing: `IPRT1` / `TIMPRTi` only at survey times (DTMAX effect on the printed step); the 1.4 GB truth vtk then shrinks
- [ ] g. End of run: archive as `_states.nc` (zstd) + meta + manifest, `prune --apply` the rest; document in README

### 7.2 RAM (HC: 41 616 nodes × NENS × cycles; ~85 GB with the old long table)

- [ ] Measure peak RSS vs NENS (8/16/32/64) with `--store-full all` and `analysis`
- [ ] Check that nb8 never rebuilds the long `df_DA` (property) for HC
- [ ] `--max-nodes` for scores
- [ ] One chunk per cycle must stay O(1) in memory

### 7.3 Compute time (each trade-off needs the 3-seed check of 7.4)

- [ ] a. Ensemble size N=16/32/64 (`compare-nens`): find the knee of RMSE / coverage vs cost; EnKF vs PF at equal cost
- [ ] b. HC mesh: DELTA 0.1 → 0.2 m for DA AND truth (nb6 must use the same delta; keep the inverse-crime flag in mind)
- [ ] c. ERT forward: coarser pyGIMLi mesh / fewer quadrupoles (rhoa mask already) / cache geometry factors; time per member per survey; forward mesh margin warning (0.5 m beyond electrodes) must not get worse
- [ ] d. Surveys: assimilate every k-th survey only (k=2,3) and measure the loss
- [ ] e. Parallelism: cores per member vs members in parallel (`parallel=True`, `REAL_JOBS`, `taskset`); CATHY not recompiled per member
- [ ] f. Restart: resume an interrupted run from the last `cycle_*.npz` (`load_backup`) instead of redoing it

### 7.4 Validation rule

- [ ] A lighter setting is accepted if, on 3 seeds, the nb8 scores stay within the seed-to-seed spread of the baseline (`run_manifest.py compare` checks that only the intended option changed)

### 7.5 Deliverables

- [ ] Table: option, disk / RAM / time saved, score change
- [ ] New defaults behind flags
- [ ] README section "Lighter DA folder"
- [ ] Updated `ligher_run.py` (renamed, see section 8)

### 7.P pyCATHY wrapper side

*Candidates found by reading the code; each one must be MEASURED first, no size below is measured.*

- [ ] **P1. Member folders**: `_create_subfolders_ensemble` copies `input/ output/ vtk/ prepro/` + the CATHY executable into `DA_Ensemble/cathy_origin`, then `copytree` of all of it NENS times. The `vtk/` and `output/` of the forward run are copied too (N+1 copies).
  - Test: copy only what a window needs; share exe / prepro / static input by symlink or hardlink (Windows: copy fallback); do not copy the forward-run vtk/output
- [ ] **P2. vtk**: members print vtk every window (`VTKF=2` in nb7) and `_cleanup_vtk_files` runs only at the END (`cathy_DA.py` ~865).
  - Test (a) `VTKF=0` for members: are sw/psi read from `output/` text files only? What needs the vtk (the per-member `saveMeshPath` ~2556/2648 "converted" mesh, ERT forward)?
  - Otherwise (b) a `cleanup_vtk='per_cycle'` mode in `DAConfig`
- [ ] **P3. Per-cycle backups are O(cycles²) I/O**: `_save_step_outputs` calls `backup_results_DA` (rewrites `dict_parm_pert`, `df_obs`, `ET_DA_xr`, `df_Archie` every cycle) AND `backup_simu()` (pickles the whole object, `<project>.pkl`) AND `ET_DA.nc`.
  - Measure each file size after cycle 1 / 10 / last
  - Candidates: write every k cycles + at the end + on crash (atexit / try-finally); `backup_simu` off by default (`DAConfig` flag); `dict_parm_pert` as compact arrays (npz / netcdf) instead of a pickle of nested dicts; float32
- [ ] **P4. `backup_fort777`** copies each member's `fort.777` at every cycle (`DAnb<N>fort.777`, NENS × cycles files) when `config.backup_fort777` is on: check the default, keep one stacked file per cycle instead
- [ ] **P5. Outputs printed**: the window run uses `TIMPRTi=[0, diff_time]` (start/end only, good).
  - Check CATHY's other outputs per member (mbeconv, hgraph, cumulative files) and which `IPRT1` level they need
  - `DAConfig` option to keep only psi/sw/ET and delete the rest after reading
- [ ] **P6. Compact store**: sw as uint16 / psi float32 (`state_store.py`); optional zstd netcdf chunks; `load_backup` / `states` must still read old runs
- [ ] **P7. RAM**: per-cycle objects kept in memory for every member (prediction matrices, `dict_parm_pert` growth)
  - Peak RSS vs NENS and nodes (HC 41 616 nodes)
  - Run members with a pool whose size is set from free RAM (cf. `REAL_MIN_FREE_GB`)
- [ ] **P8. Options in `DAConfig`** with documented defaults that keep today's behaviour (`light=False`), plus a `light=True` preset
  - Unit tests on stubs like `tests/test_compact_states.py`
  - Changes go to pycathy_wrapper on a branch rebased on main AFTER the `eta-assimilation-improvements` rebase (same files: `cathy_DA.py`, `cathy_tools.py`)

## 8. Reorganise the files

*Goal: better readability.*

- [ ] **Plotters in separate files** (a `plotting/` package or `*_plots.py` modules)
  - Split the figure functions out of `nb8_plots_synthetic.py` (`fig_states`, `fig_profiles`, `fig_metrics`, `fig_params`, `compare_*`), `nb5_DA_ERT_soil3D_nc_plots.py`, `nb4_DA_soil3D_nc_plots.py`, `plot_6_ESMDA_synthetic.py`, `nb5_plot_helpers.py` / `nb8_helpers.py`
  - The `nbX` scripts keep only the CLI (arguments, loop over runs, calling the plotters)
- [ ] **Rename files whose name does not say what they do.** Candidates:
  - `nb1..nb8` → e.g. `forward_model_nc.py`, `synthetic_truth_ert.py`, `da_ert_synthetic.py`, `plots_synthetic.py`
  - `fig0_infiltr_compatt` / `fig1_infiltr_scheme` notebooks
  - `plot_6_ESMDA_synthetic.py`
  - `ligher_run.py` (typo vs `lighten_run.py` in its docstring)
  - `nc_model.py` vs `model_config.py` (overlap)
  - `ert_data` / `ert_frame` / `ert_alignment` / `ert_schedule` (one `ert/` folder?)
- [ ] **When renaming:** `git mv` (keep history) and update every reference (`run_twin.sh`, `run_matrix.py`, `data_path.toml` `[scenarios]`/`[code]`, README, `synthetic_twin_plan.md`, tests, imports); check with `grep` and `pytest`; add the old → new table to the README
- [ ] **Hygiene while at it:** see section 2 (`__pycache__/`, `.ipynb_checkpoints/`, large tracked files)

## 9. Branches (pycathy_wrapper)

*Checked 2026-11-24 on `BenjMy/pycathy_wrapper`; the branches live there, not in this repo.*

- [x] Which is the last updated DA ERT branch (`feat/withIrr-on-main`, `feature/eta-assimilation-improvements` or `twin-experiment`)? Should I rebase?
  - **DA branch to use: `feature/eta-assimilation-improvements`** (4 commits not in main: `set_atmbc_components` / rain-ETp perturbation / mass-balance ETa operator, 2 mbeconv importer fixes, `plot_state_timeseries`). It is 7 commits behind main (merge-base 2026-10-26, before the mesh-vtk-properties work).
  - **Yes, rebase it on main**: `git rebase origin/main` is clean (no conflict, py_compile OK), then fast-forward/PR into main.
  - `feat/withIrr-on-main` is the most recently updated (2026-11-16) but it is the irrigation work (`irrigation.py`, `cathy_tools`, `version_config`), not DA; it is 0 behind main (main already merged into it). It touches the same files as eta (`cathy_tools.py`, `cathy_plots.py`) but merges with eta without conflict (checked): merge them in separate PRs.
  - `withIrrigation` is the older version of `withIrr-on-main` (68 behind): delete after `withIrr-on-main` is merged.
  - `twin-experiment` does not exist on origin of pycathy_wrapper, tmpCATHY_ERT_DA_real or mfc2_EA-EM_TDR: local-only? Push it or tell me where it is; not assessed.
- [x] Can I close `origin/feat/material-vgp-control` and `origin/mesh-vtk-properties`?
  - Yes, both. Each is an ancestor of `origin/main` (0 commits ahead; `git merge-base --is-ancestor` = YES; `mesh-vtk-properties` was merged by PR #33). Nothing is lost by deleting them:
    `git push origin --delete feat/material-vgp-control mesh-vtk-properties` (and `git fetch --prune` elsewhere).
