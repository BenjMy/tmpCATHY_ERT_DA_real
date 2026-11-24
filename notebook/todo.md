Repo: https://github.com/BenjMy/tmpCATHY_ERT_DA_real/tree/feature/data-paths-registry-readme

- [x] Add a manifest file for each run to check which version of CATHY DA ERT (including the branch) is used and which version of CATHY and pyCATHY as well (including the branch). And other usefuls things as you consider for reproducibility.
- [x] feat/withIrr-on-main or feature/eta-assimilation-improvements or twin-experiment? Help me to choice the last updated DA ERT branch.  Should I rebase?
  (checked 2026-11-24 on BenjMy/pycathy_wrapper; the branches live there, not in this repo)
  - DA branch to use: `feature/eta-assimilation-improvements` (4 commits not in main: set_atmbc_components / rain-ETp
    perturbation / mass-balance ETa operator, 2 mbeconv importer fixes, plot_state_timeseries). It is 7 commits behind main
    (merge-base 2026-10-26, before the mesh-vtk-properties work).
  - Yes, rebase it on main: `git rebase origin/main` is clean (no conflict, py_compile OK), then fast-forward/PR into main.
  - `feat/withIrr-on-main` is the most recently updated (2026-11-16) but it is the irrigation work (irrigation.py,
    cathy_tools, version_config), not DA; it is 0 behind main (main already merged into it). It touches the same files
    as eta (cathy_tools.py, cathy_plots.py) but merges with eta without conflict (checked): merge them in separate PRs.
  - `withIrrigation` is the older version of withIrr-on-main (68 behind): delete after withIrr-on-main is merged.
  - `twin-experiment` does not exist on origin of pycathy_wrapper, tmpCATHY_ERT_DA_real or mfc2_EA-EM_TDR: local-only?
    Push it or tell me where it is; not assessed.
- [x] By the way can I close origin/feat/material-vgp-control? and origin/mesh-vtk-properties? 
  - Yes, both. Each is an ancestor of origin/main (0 commits ahead; `git merge-base --is-ancestor` = YES;
    mesh-vtk-properties was merged by PR #33). Nothing is lost by deleting them:
    `git push origin --delete feat/material-vgp-control mesh-vtk-properties` (and `git fetch --prune` elsewhere).
- [ ] I want a data_path.toml where I can control all the path for inputs, outputs, plots, codes lib etc ...
- [ ] Register all the run in a csv files where I can keep track what has been done or not
- [ ] Investigate if it is possible to get ligher DA folder (.pkl files are heavy)
  Protocol for a lighter simulation (investigation only, nothing implemented; written 2026-11-24 from the code + synthetic_twin_plan.md s.16)
  Already there: compact float32 StateStore (`_states/cycle_*.npz`, 314 MB long df_DA -> ~35-65 MB), `--store-full analysis
  --store-stride`, `--cleanup-vtk` at the END of a run, `ligher_run.py report|convert|prune`. Known gaps: README has no
  "Lighter DA folder" section although ligher_run.py's docstring points to it; `_df_DA.pkl` is still written once at the end;
  cleanup_vtk acts only at the end (peak disk = NENS x ~38 MB x ..., vtk is ~31 of 38 MB per member); plan s.16 "not done":
  VTKF=0 for the truth run (1.4 GB vtk, untested), sw as uint16.
  - [ ] 0. Baseline (one reference run, never changed): NC, scenario ic_Ks_upd, EnKF, N=8, seed 1 (22 surveys, ~25 min on 2 cores)
        + the same on HC with N=8. Record: `du -sh` per category (`ligher_run.py report`), peak RSS (`/usr/bin/time -v`),
        wall time split CATHY windows / ERT forward (pyGIMLi) / analysis (add timers), and the manifest.json of the run.
        Metrics to protect: nb8 summary.csv (RMSE, spread, 90 % coverage, CRPS) and parameters_summary.csv (z prior -> final)
  - [ ] 1. Disk (no change of the result expected: check states are bit-identical or < 1e-6 relative)
        a. VTKF=0 for members and truth: does `read_outputs` (sw, psi) need the vtk? if not, never write them
        b. else delete each member's vtk after every cycle (not only at the end): peak disk becomes O(1 cycle)
        c. per-member folder: split size into input files (mesh, dem, soil: identical for all members -> symlink/hardlink
           instead of NENS copies?) vs outputs; keep restart files only
        d. do not write `<run>_df_DA.pkl` when the compact store exists (flag, default off); parameters pkl: check size
        e. dtype: psi float32, sw uint16 (22 MB vs 314 MB) or float16: error must stay << 5 % ERT noise
        f. printing: IPRT1 / TIMPRTi only at survey times (DTMAX effect on the printed step), then 1.4 GB truth vtk shrinks
        g. end of run: archive as `_states.nc` (zstd) + meta + manifest, `prune --apply` the rest; document in README
  - [ ] 2. RAM (HC: 41 616 nodes x NENS x cycles; ~85 GB with the old long table)
        measure peak RSS vs NENS (8/16/32/64) with `--store-full all` and `analysis`; check that nb8 never rebuilds the long
        df_DA (property) for HC; `--max-nodes` for scores; one chunk per cycle must stay O(1) in memory
  - [ ] 3. Compute time (accuracy trade-offs: each needs the 3-seed check of step 4)
        a. ensemble size: N=16/32/64 (`compare-nens`), find the knee of RMSE / coverage vs cost; EnKF vs PF at equal cost
        b. HC mesh: DELTA 0.1 -> 0.2 m for DA AND truth (nb6 must use the same delta; keep the inverse-crime flag in mind)
        c. ERT forward: coarser pyGIMLi mesh / fewer quadrupoles (rhoa mask already) / cache geometry factors; time per
           member per survey; forward mesh margin warning (0.5 m beyond electrodes) must not get worse
        d. surveys: assimilate every k-th survey only (k=2,3) and see the loss
        e. parallelism: cores per member vs members in parallel (`parallel=True`, REAL_JOBS, taskset), CATHY not recompiled per member
        f. restart: resume an interrupted run from the last `cycle_*.npz` (`load_backup`) instead of redoing it
  - [ ] 4. Validation rule: a lighter setting is accepted if, on 3 seeds, the nb8 scores stay within the seed-to-seed spread
        of the baseline (`run_manifest.py compare` to check that only the intended option changed)
  - [ ] 5. Deliverables: table (option, disk / RAM / time saved, score change), new defaults behind flags, README section
        "Lighter DA folder", updated `ligher_run.py` (renamed, see the reorganisation item)
- [ ] Create a new branch on the github and add a README that explains all the code architecture, from inputs type, how to run, scenario construction, outputs data and plots. 
- [ ] Explain synthetic ERT DA pipeline in the README
- [ ] Explain how to use scenario.toml file in the README
- [ ] Check this warming:   warnings.warn("truth forward on the SAME mesh as the DA: mesh discretisation error is not tested "
- [ ] ERT frame: legacy (farming_DA_ERT3d.py); signs=(1, 1, 1), shift=[2.0, 2.0, 0.001]
/home/ben/Nextcloud/BenCSIC/Codes/BenjMy/pycathy_wrapper/pyCATHY/ERT/alignment.py:116: UserWarning: [ERT] pg forward mesh extends only 0.50 m beyond the electrodes (array span 1.10 m): boundary effects can bias the forward apparent resistivities of weak-signal quadrupoles (see nc_model.unreliable_quadrupoles: homogeneous test, rhoa/rho must be 1)
  warnings.warn(f"[{label}] pg forward mesh extends only {margin:.2f} m beyond the 

- [ ] Reorganise the files for better readability
  - [ ] Plotters in separate files (a `plotting/` package or `*_plots.py` modules): split the figure functions out of
        nb8_plots_synthetic.py (fig_states, fig_profiles, fig_metrics, fig_params, compare_*), nb5_DA_ERT_soil3D_nc_plots.py,
        nb4_DA_soil3D_nc_plots.py, plot_6_ESMDA_synthetic.py, nb5_plot_helpers.py / nb8_helpers.py; the nbX scripts keep only
        the CLI (arguments, loop over runs, calling the plotters)
  - [ ] Rename files whose name does not say what they do. Candidates: nb1..nb8 (-> e.g. forward_model_nc.py,
        synthetic_truth_ert.py, da_ert_synthetic.py, plots_synthetic.py), fig0_infiltr_compatt / fig1_infiltr_scheme notebooks,
        plot_6_ESMDA_synthetic.py, ligher_run.py (typo vs lighten_run.py in its docstring), nc_model.py vs model_config.py
        (overlap), ert_data / ert_frame / ert_alignment / ert_schedule (one `ert/` folder?)
  - [ ] When renaming: `git mv` (keep history) and update every reference (run_twin.sh, run_matrix.py, data_path.toml
        [scenarios]/[code], README, synthetic_twin_plan.md, tests, imports); check with grep and `pytest`; add a table
        old name -> new name in the README
  - [ ] Hygiene while at it: untrack `__pycache__/` and `.ipynb_checkpoints/` (add to .gitignore); decide where the large
        tracked files go (notebook/logs, data/ duckdb.exe)

use .patch file that you commit as me as author (BenjMy)
