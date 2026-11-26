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
- [x] Run the lighter version: `CASE=NC DA=pf BIAS=none NENS=32 EXTRA="--store-full analysis --store-stride 10 --cleanup-vtk all" ./run_twin.sh pair full_upd` (wrap in `lighter_protocol.py measure` to record time / RSS / disk, see 7.0). *Scenarios run on the workstation, results not yet interpreted: see 7.6*
- [ ] Interpret the scenarios already run: send the files listed in 7.6
- [ ] Send a colleague the results: `python export_results.py ~/share/twin_results --zip`
- [ ] Run `pf_inbreeding.py` on the real PF run (the IC prior fix is done, see 5.1 (a); rerun needed)
- [ ] Write the series of tests that the twin (DA model) is set up like the truth: mesh and layers, materials, `ic`, `atmbc`, times, observations (see 5.3). *Config-level part done (22 tests, no CATHY); open: the items that need the truth pkl / a run, and two findings (HC zone 2 bias, solver `parm` differs)*

## 2. Repository hygiene

- [x] `.gitignore`: `simulations*`, `*_meta/`, `*.log`, `*.pkl`
- [x] Untrack `data/3d/synthetic/NC/truth_NC.pkl`
- [ ] Still tracked although now ignored (a `.gitignore` never untracks): decide whether to `git rm --cached` them
  - `notebook/cache/*.pkl` (6 files, `prior_NC_*`)
  - `notebook/logs/*.log` (24 files, `SYN_NC_*`)
- [x] Untrack `__pycache__/` and `.ipynb_checkpoints/` (40 tracked files) and add them to `.gitignore`
- [ ] Failing test on the untouched branch: `tests/test_nb5_plot_helpers.py::test_evolution_and_shift_flags` (`parm_shift_table` flag for `n_VG0`: expects `SHIFT` and `COLLAPSE`)
- [ ] Decide where the large tracked files go: `notebook/logs`, `data/**/duckdb.exe` (3 × 24 MB)

## 3. Reproducibility and run tracking

- [x] Manifest file for each run: version of CATHY DA ERT (with branch), CATHY and pyCATHY (with branch), and other useful reproducibility info → `run_manifest.py`
- [x] Register all the runs in a CSV file to track what has been done or not → `registry/runs_registry.csv`
- [x] `data_path.toml` to control all the paths: inputs, outputs, plots, code lib, etc. → `data_path.toml` + `project_paths.py` (commit `9f999c5`, 6 tests)
  - read by nb5–nb8, `run_matrix.py`, `run_twin.sh`; check with `python notebook/project_paths.py`
  - nb1–nb4a (incl. nb3 sensors / `site_params.json`) were still hardcoded (`../simulations/`, `../data/…`): now read the toml
  - [ ] still hardcoded: `infiltration_utils.py` (`rootfilePath + 'data/3d/…'`, `../data/2d/`, `../data/3d/`), `model_config.ert_paths` / `nc_model.ert_paths` default `root='../'` (the toml wins when `root` is the project root)

## 4. Documentation (README)

- [ ] Create a new branch on GitHub and add a README explaining the whole code architecture: input types, how to run, scenario construction, output data and plots
- [ ] Explain the synthetic ERT DA pipeline
- [ ] Explain how to use the `scenario.toml` file
- [x] Add a "Lighter DA folder" section (`lighten_run.py` docstring already points to it)
- [x] Add a table old name → new name once files are renamed (section 8; done for the files renamed so far, to extend with the next renames)

## 5. Scientific checks on the DA runs

*Tools and findings from reading the code (`make_scenarios.py`, `nb7`, `pyCATHY/DA/pf.py`, `perturbate.py`) and from simulating the priors. **Not yet confirmed on a real run**: `simulations/` is not in the repo.*

### 5.1 Are the "unbiased" cases really unbiased?

*`--prior-bias none` only means "prior nominals = truth values". Tool: `python notebook/prior_bias_check.py [--case NC] [--scenario full_upd] [--biased]` (offset of the SAMPLED prior mean from the truth, in prior sigmas, on the scale the filter works on; 6 tests).*

- [x] Tool written (`prior_bias_check.py`, `tests/test_prior_bias_check.py`)
- [x] **Found: the initial-condition prior is not centred on the truth, with or without `--prior-bias none`.** `make_scenarios.entries()` hardcodes `ic` = −1.5 (σ 0.5) and nb7 starts the DA model from `ic_values=[[−1.5]]*3`, but the truth is layered `[−2.4, −1.3, −1.15]`: offset **+1.8σ** (layer 0), −0.4σ, −0.7σ → the ensemble mean is wetter than the truth at the top and drier below, which shows up in the profile of (ensemble mean − truth), while the truth can still lie inside the ensemble (|offset| < 2σ). One scalar cannot match 3 layers.
- [x] Checked and fine (|offset| < 0.06σ): VGP (SATCELL, n, thetar) and Archie (a, m, n) nominals = truth.
- [x] Ks: sampled as a truncated normal (σ = nominal, lower bound nominal/100) → arithmetic mean ×1.29, median ×1.2 of the truth, but on the log scale the filter uses the offset is only −0.02σ. Not a filter bias; only matters when the arithmetic ensemble-mean Ks is plotted against the truth. (The comment in `make_scenarios.py` about draws piling up at the lower bound is outdated: the library now truncates.)
- [x] **(a) per-layer `ic` parameter (library change), prior mean = `ic_true`** — done (not yet run on a real DA run):
  - *library (`pycathy_wrapper`, branch `feat/ic-per-layer-material`, patch `0001-pycathy-ic-per-group-material.patch`)*: `ic` accepts `pert_control_name = "material"`: one entry per group (`ic0, ic1, ic2`) with its own mean / sigma / bounds; groups from the new scenario key `ic_layers` (fallback `material_layers`). `DA._update_ic_member` expands the group values to every layer and writes `input/ic` at cycle 0 only (the analysed state is kept afterwards); ES-MDA reads `ini_perturbation`, so it works too (the guard for `layers` is unchanged). 11 tests (`tests/test_ic_material.py`), related suites unchanged (254 passed).
  - *why `ic_layers` and not `material_layers`*: by depth (`IC_RANGES`, mid-depth rule of `map_param_to_layers`) the ic groups are `[0-5] [6-9] [10-14]`, the soil materials `[0-5] [6-8] [9-14]`: layer 9 (0.36-0.44 m) is material C but in the second ic range.
  - *this repo*: `make_scenarios.py` writes `ic` as `material` with `per_mean = per_nom = ic_true` (NC `[-2.4, -1.3, -1.15]`, HC `[-1.4, -1.13, -1.13]`), sigma 0.5 per layer, bounds `[-5, -0.05]` and the `ic_layers` key, in every scenario and **with or without** `--prior-bias` (the soil bias never touched ic). nb7 starts the DA model from the group means; `--ic-prior-mean M` is still there and now means "one scalar M for every group" (`-1.5` = the old prior, to compare). `prior_bias_check.py` reports the per-layer offset (now ~0σ) and has `--ic-scalar M` for the old one. 5 new tests (`tests/test_ic_per_layer.py`), `test_prior_bias_check.py` flipped.
  - *needs the library branch*: `code.pycathy` in `data_path.toml` (or the installed pyCATHY) must contain the patch, otherwise the validator refuses `ic` + `material`.
  - *consequence*: results of earlier runs are **not comparable** (the `ic` prior changed); `ic0..ic2` replace the single `ic` key in `dict_parm_pert` / `df_DA` parameters. nb4-nb6 (real data) still use their own scalar `ic` scenario.
- ~~(b) centred scalar for `none` mode (mean of `ic_true` = −1.62)~~ and ~~(c) keep it and rename the "unbiased" label~~: superseded by (a); (c) still applies to a run started with `--ic-prior-mean`.
- [ ] Rerun the unbiased PF case with the per-layer ic (and once with `--ic-prior-mean -1.5` for the comparison) and compare the `*_profiles.png` (mean − truth by depth)
- [x] Other (non-bias) reasons the mean differs from the truth, to separate from the above: finite N (σ/√N), nonlinearity (mean of sw ≠ sw at the mean parameters), PF collapse (5.2). Tool: `python notebook/mean_offset_budget.py [--case NC] [--n 16 32 64] [--n-eff 12 2]` (no CATHY run; `tests/test_mean_offset_budget.py`, 6 tests). *Results (NC, per-layer `ic` prior centred on `ic_true`):*
  - *finite N: expected |offset| of the ensemble mean = 0.20σ (N=16), 0.14σ (32), 0.10σ (64); 95 % of the runs within 0.49 / 0.35 / 0.24σ. With the PF collapse of 5.2 the relevant number is the founders / N_eff: 12 → 0.23σ (95 %: 0.57σ), **2 founders → 0.56σ (95 %: 1.39σ)**. So `prior_bias_check.py --tol 0.2` is below the noise of a single run (N < 96): only the sampled prior (200 000 draws) can show a 0.2σ bias, a run cannot.*
  - *nonlinearity of the initial state: with an exactly centred `ic` prior the mean of sw over the prior is 0.002 (layers 0-5) to 0.006 (layers 6-14) above sw at `ic_true`, i.e. ≤ 0.19 sd of sw: small. A mean − truth profile of the initial ensemble larger than that is NOT nonlinearity.*
  - *reading rule for the `*_profiles.png`: (mean − truth) within ~σ_sw/√N_eff of 0 (+ the small offset above) is sampling noise; a persistent offset beyond that points to the prior (5.1) or to the filter (5.2). Assumption: `VGPSATCELL` = air-entry head ψ_s of the van Genuchten curve (not checked against CATHY).*
  - *not covered (needs a run): the PF collapse itself (5.2) and the nonlinearity of the ERT forward model.*

### 5.2 Filter inbreeding (`SYN_NC_full_upd_pf_N32_s1_unbiased_issues`)

*Tool: `python notebook/pf_inbreeding.py <project>_states/ [--csv ..] [--plot ..]` (reads `cycle_*.npz`, no pyCATHY; 4 tests). It recovers the resampling parents because the analysed state is an exact copy of the parent's state (`use_enkf_update=False`, `pf_jitter_std_state=0`), chains them and reports distinct founders, founders ESS, distinct states and spread per cycle.*

- [x] Tool written (`pf_inbreeding.py`, `tests/test_pf_inbreeding.py`: known genealogies recovered exactly)
- [x] Why it is expected from the nb7 settings (`pf_resample_threshold=1.0`, `pf_jitter_std_state=0.0`, `pf_temper_neff_frac=0.3`, `pf_liu_west_a=0.98`):
  - tempering only guarantees N_eff ≥ 0.3·N, but the threshold 1.0·N is above any N_eff → **resampling at every analysis** (22 surveys → 22 times)
  - states get no jitter, so copies of a particle are identical at the analysis; only the parameters are rejuvenated (Liu–West, noise std √(1−a²) ≈ 0.2 × weighted std)
  - *simulation, not your run*: N=32, tempered N_eff ≈ 12, resampling every cycle → 2 founders left after 22 cycles. Static parameters (Ks, VGP, Archie) then descend from 2–3 initial draws and can only move by the Liu–West noise
- [x] `diag_pf_run.py RUN [--compare RUN_B]`: one report for the checks below from the small files (status, weights, founders, dead members, parameter collapse, ERT data space, next commands); `tests/test_diag_pf_run.py` (9 tests). *First reading of `SYN_NC_full_upd_pf_N32_s1_unbiased` (run interrupted at cycle 10): N_eff very low in every cycle; 12 of 15 parameters lost 60-86 % of their spread and 4 end >= 3 final SDs from the truth; members 15-19 rejected at cycle 1 (CATHY stopped at 95-98 % of the window) and still NaN at cycle 10 (ensemble of 27, slots never refilled); ERT data-space spread is about the size of the misfit but 80-96 % of the variance is in one mode (redundant data, weights collapse, depth profile not constrained)*
- [ ] **Run `pf_inbreeding.py` on the real N32 run** (and N64) and confirm: founders left, distinct states, spread loss; add `--plot` (`diag_pf_run.py` does it when `<run>_states/` is next to the small files)
- [ ] Find why members 15-19 died at cycle 1 (consecutive indices: same worker / folder?) and make the PF refill dead slots from live parents
- [x] `export_results.py` now copies `<run>_run.json`, `<run>_rejections.json`, `<run>_ERT_checks.json` (when they exist) to `runs/<run>/`, so `python diag_pf_run.py runs/<run>` works on a shared folder (`tests/test_export_results.py::test_copies_the_json_files_diag_pf_run_needs`). *The `_states/` folder is still not exported (big): `pf_inbreeding.py` needs the real run folder*
- [ ] Log of the reflections prints identical before / after values (SATCELL_*, a_Archie) and the values pile up just below 0.5 while the printed bound is [0.1, 2.0]: check which bound acts
- [ ] Perturb Ks (and other strictly positive parameters) in log space: negative analysed values (Ks2 -4.3e-07 at cycle 1) are reflected
  - *Found 2026-11-25: the default Ks prior `normal100` (Gaussian, sd = 100 % of nominal, `make_scenarios.py`) gives ~16 % negative draws piled at the lower bound. `KS_PRIOR=lognormal KS_SIGMA_LOG=0.5` (log-normal, always positive) already exists but is untested: run it, README "If the ensembles collapse", test 3*
- [ ] Compare with EnKF at the same seed / N / scenario (`--da-type enkf`): same bias profile? then the cause is the prior (5.1), not the PF
- [ ] Run the 5 tests of README "If the ensembles collapse" (N 32/64, inflation, log-normal Ks, PF jitter + resample threshold, unbiased prior) and compare `summary.csv` / `parameters_summary.csv`. *Seen at N=16 (4 scenarios, seed 42): state coverage <= 0.03 (PF `ic_VGP_Ks_upd` 0.31), parameter contraction 0.94-0.9998, `frac_truth_in_90` ~ 0, z final 89-14668; Ks relative error 5.2 -> 15.5 (EnKF) / 5.4 (PF) in `ic_Ks_upd`*
- [ ] Candidate remedies, each to TEST (nothing changed): `pf_jitter_std_state` > 0 (now `--pf-jitter-std-state`); larger N; resample only when N_eff < N/2 (check first that the weights are carried between cycles when not resampling: not verified); use the PF only for parameters
- [ ] pyCATHY side: store the resampling indices per cycle (tiny change in `pf.py` / `AnalysisResult.extra`) so the genealogy is exact instead of reconstructed

### 5.3 Twin set-up must match the truth run (series of tests)

*Why: the truth (`synthetic_truth_ert.py`, nb6) and the DA model (`da_ert_synthetic.py` -> `model_config.setup_model`) are built by two different scripts. Anything that differs besides the declared prior (soil bias, `ic` prior, Ks prior) is an unintended bias and would be read as a filter problem. Already checked in the code: `delta` (assert in `da_ert_synthetic.py`) and the noise % (warning). Everything else below is NOT checked yet.*

*Reference = the `truth_<CASE>.pkl` (keys `delta, atm_times, nodes, soil_table, ic_nodes, params, archie_truth, fwd_mesh, da_mesh, ert_frame, surveys, tdr_nodes, ...`). Target = the DA project built by `setup_model` (CATHY input files in `<run>/input/`) with `--prior-bias none` and the truth parameters. Idea: build the DA model WITHOUT running CATHY (`dry-run` / preprocessor only) and compare it with the truth, for NC and HC. Put them in `tests/test_twin_vs_truth.py`; the ones that need CATHY/pyCATHY are skipped when it is not installed (like `test_nb5_dryrun.py`).*

*Done (config level, **no CATHY needed**): both projects are built by the same `model_config.setup_model`, so `tests/test_twin_vs_truth.py` runs it twice on a recording stand-in of the CATHY object (`FakeCathy`; `pyCATHY.meshtools` faked for these tests only) with the arguments of `synthetic_truth_ert.py` and of `da_ert_synthetic.py`, and compares what is written (mesh inputs, soil table, `ic`, `atmbc`, `parm`, boundary conditions). 22 tests pass for NC and HC; 3 deliberate breaks of `model_config.py` (irrigation hour, ic ignored, bias on POROS) are caught. What still needs the real files / a run is in the unticked items and as 7 skipped tests with the same names.*

- [ ] **Mesh / geometry**
  - [ ] same 3D mesh: number of nodes and elements, node coordinates (`truth['nodes']` vs DA mesh points), same DEM / surface and bottom elevations
  - [x] number of layers (z-layers of the CATHY mesh) equal in truth and DA; layer thickness / depth of each layer identical (`test_mesh_layers_and_nodes`, `test_delta_override_changes_truth_and_da_together`; `ic` node count = (n+1)² × (layers+1))
  - [ ] ERT forward mesh (`fwd_mesh`) and DA mesh (`da_mesh`) are the files the truth records, and the ERT frame (`signs`, `shift`, `origin`) is the same
  - [ ] TDR / sensor nodes (`tdr_nodes`, `tdr_depths`) are the same nodes in the DA observation operator
- [ ] **Materials / soil**
  - [x] number of soil zones and materials (`soil_table`) equal; mapping layer -> material (NC `[0-5] [6-8] [9-14]`, see 5.1 (a)) identical (`test_soil_table_rows_follow_the_material_layers`, `test_hc_zone2_below_material_a_is_zone1_in_truth_and_da`)
  - [x] with `--prior-bias none`: every (zone, material) of the DA soil table = truth value (Ks, SATCELL / porosity, n, VGP, thetar, ...); with `--biased`: differs exactly by `biased_overrides()` and nothing else (`test_unbiased_da_model_is_identical_to_the_truth_model`, `test_scenario_nominals_equal_the_soil_table_written_to_cathy`, `test_biased_da_model_differs_from_truth_only_by_biased_overrides`: only PERMX/Y/Z, VGPSATCELL, VGNCELL, VGRMCCELL change, not POROS / ELSTOR). *The factor VALUES (×1.3, −0.08, +0.01, ×5) are not re-checked: they are the definition of the bias*
  - [ ] **Decide (found by `test_biased_hc_leaves_zone2_below_a_at_truth_values`)**: with `--biased` on HC, `biased_overrides()` only has the entries of `CASES['HC']['soil']`, i.e. zone 2 material A. Zone 2 B and C (copied from zone 1 by `zone2_fill_below`) stay at the TRUTH while zone 1 B and C are biased; they are not DA parameters either (`groups('HC')` = zone 1 A/B/C + zone 2 A). Is the HC biased prior meant to leave the compacted band below 20 cm unbiased? Either bias them too (add the keys in `biased_overrides`) or document it; the test pins the current behaviour
  - [ ] porosity per node (`porosity_nodes`) and the porosity used by Archie in the ERT forward model (`porosity_mode`) are the same
  - [ ] Archie truth vs prior (`archie_truth`, `archie_prior`, `rFluid`): the DA nominals = truth when unbiased
- [ ] **Initial condition**
  - [x] DA `input/ic` at cycle 0 = `truth['ic_nodes']` mapped by layer (per-layer `ic_true`, NC `[-2.4, -1.3, -1.15]`, HC `[-1.4, -1.13, -1.13]`), with the old scalar `--ic-prior-mean` flagged as a known difference (`test_ic_per_layer_prior_mean_equals_the_truth_ic_nodes`, `test_old_scalar_ic_prior_is_the_only_difference_and_is_flagged`: only `ic` differs, +0.9 m = +1.8σ in layer 0). *`truth['ic_nodes']` from the pkl itself is not read yet*
  - [ ] same unit / convention (pressure head vs saturation) and same spin-up / `stabilize_period`
- [ ] **Atmospheric boundary conditions and time**
  - [x] `atmbc` of the DA run = truth `atm_times` and values (rain / irrigation, ETp): same number of rows, same times, same units (`test_atmbc_and_parm_times`: ET everywhere, irrigation pulse only in the central 2 × 2 m at `irr_start_hour`; output times `TIMPRTi` = atmbc times). *Against the `atm_times` stored in the pkl: still open, next item*
  - [ ] time vector: `build_times(ert_times)` == `atm_times` (already used for an assert in `da_ert_synthetic.py` line ~285: extend to the values), same `sw_time`, same ERT survey times (`surveys[*]['t_h_irr']`)
  - [ ] with `--n-times N`: the first N surveys only, and the truth is cut at the same last survey
- [ ] **Other CATHY inputs**: boundary conditions (nodal / surface), root map / ET parameters, `parm` file (solver tolerances, time stepping, `DELTAT`, `TMAX`), `soil` file ordering, output times (`IPRT`), `delta` -> all identical between truth and DA input folders. *Done for everything `setup_model` writes (no-flow bc, `DTMAX/DELTAT/DTMIN`, `IPRT1`, `VTKF`, `TRAFLAG`, `ic` flags: `test_boundary_conditions_and_solver_settings_are_fixed`). Open: the files as written in the two real `input/` folders.*
  - [ ] **Found by reading the scripts (not by a test, the block is in the script): the DA run overwrites the `parm` settings after `setup_model`** (`da_ert_synthetic.py` l.536-538) and they differ from the truth (`synthetic_truth_ert.py` l.150, `setup_model`):
    | | truth | DA (EnKF / PF) | DA (ES-MDA) |
    |---|---|---|---|
    | `DELTAT` (first time step) | 0.1 s | `min(1e2, 0.5 × shortest window)` = up to 100 s | `min(that, esmda_dtmax)` |
    | `DTMAX` | 1e3 s | 1e3 s | `--esmda-dtmax` = 60 s |
    | `IPRT1` | 4 | 2 | 2 |
    | print times | dense (every `--out-dt`) + interpolated to the atmbc times | atmbc times only | atmbc times |
    *Why it may matter: each DA window restarts CATHY, so the first time step of every window is 100 s in the twin and 0.1 s in the truth after the irrigation pulse; the Richards solution can differ by more than the discretisation error you want to test (truth and DA also share the mesh, see the warning in 6). Size of the effect unknown: the open-loop sanity test below measures it (noDA, N=1 at the nominal). If it is not negligible: run the truth with the DA's `DELTAT / DTMAX`, or the DA with the truth's.*
- [ ] **Observations**
  - [ ] ERT data of the twin = the truth data (same electrode sequence, geometry, noise type / % / seed: `truth['noise']` vs `ERT_DATA_ERR`); error model used in the filter = the noise added to the truth
  - [ ] TDR sigma (`tdr_sigma`) in the filter = the one used to perturb the virtual TDR
- [ ] **Open-loop sanity test (cheap)**: run the DA model with `noDA`, unbiased prior, no perturbation (N=1 / all members at the nominal) and compare its states with the truth `sw` / `psi`: they must agree to solver precision. If not, the set-up differs somewhere in the list above, before any filter is involved
- [ ] Report the diff as a table (item, truth, DA, OK / DIFF) written to `<run>_meta/twin_vs_truth.json`, and register it in the manifest so a run with a DIFF is visible in the registry
- [ ] Same checks for the HC case (the HC runs failed at 13:02 before the pycathy patch: first `CASE=HC ./run_twin.sh truth` must pass)

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
- `lighten_run.py report|convert|prune` (*the committed file was cut after `report`, 77 lines, did not compile: completed and tested, 2026-11-24*)

**Known gaps**
- ~~README has no "Lighter DA folder" section~~ *(added)*
- `_df_DA.pkl` is still written once at the end
- `cleanup_vtk` acts only at the end: peak disk = NENS × ~38 MB × …, vtk is ~31 of 38 MB per member
- Plan s.16 "not done": `VTKF=0` for the truth run (1.4 GB vtk, untested), `sw` as uint16

### 7.0 Baseline (one reference run, never changed)

- [ ] NC, scenario `ic_Ks_upd`, EnKF, N=8, seed 1 (22 surveys, ~25 min on 2 cores), plus the same on HC with N=8
- *Tool ready: `lighter_protocol.py measure OUT.json --run-dir RUN -- CMD` (wall, peak RSS, disk per category). The runs themselves need CATHY + pyGIMLi: to do on the workstation.*
- Record:
  - `du -sh` per category (`lighten_run.py report`)
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
  - *Read of `nb8_plots_synthetic.py` only: it calls `load_backup(prefer=..., df_DA=False)` and uses the long table only when `states` is missing (old runs). Whether `load_backup` itself rebuilds it inside pyCATHY is NOT checked (pyCATHY not available); measure peak RSS on an HC run.*
- [x] `--max-nodes` for scores: already in `nb8_plots_synthetic.py` (default 4000 nodes, stride on the node axis)
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
  - *Rule implemented: `lighter_protocol.py accept BASE*.csv --light LIGHT*.csv` (shift <= k x baseline std, default k=1, exit code 0/1). Still to run on real seeds.*

### 7.5 Deliverables

- [ ] Table: option, disk / RAM / time saved, score change
- [ ] New defaults behind flags
- [x] README section "Lighter DA folder"
- [x] `lighten_run.py` completed (report / convert / prune) + `lighter_protocol.py measure|accept` with tests (`tests/test_lighter_protocol.py`); renamed `lighten_run.py` (section 8)

### 7.6 Comparison analysis of the runs already done

*What is needed to interpret the runs (small files only; no `_states/`, `.npz`, `.nc`, `.vtk`, `.pkl`). Paths relative to the project root; the simulations folder is `[outputs].simulations` of `data_path.toml`.*

**Files to provide**
- [ ] `registry/runs_registry.csv`: which runs exist, status, `duration_min`, commit, branch
- [ ] Output of `CASE=NC ./run_twin.sh status` (one line per run: ok / partial / no states / failed)
- [ ] nb8 tables: `<simulations>/figures_nb8/summary.csv` and `parameters_summary.csv` (RMSE, spread, 90 % coverage, CRPS; z prior → final). Figures (`compare_*.png`) only if wanted
- [ ] Per run `<run>_meta/manifest.json` and `run_meta.json`: exact command (incl. `EXTRA`), seed, nens, git state
- [ ] `python lighten_run.py report <run>` text of each run: disk per category
- [ ] `lighter_protocol.py measure` json of each run (wall time, peak RSS, disk), if the run was wrapped
- [ ] The `*_issues.log` of runs marked `failed` / `partial` (e.g. `SYN_NC_full_upd_enkf_N32_s1_noDA_unbiased`: 28 chunks but status failed, 1 error)
- [ ] State which run is the baseline and which one is the lighter variant (same seed, same scenario, only one option changed)
- *One command to collect them (adjust the paths): `zip send.zip registry/runs_registry.csv simulations/figures_nb8/*.csv simulations/*_meta/manifest.json simulations/*_meta/run_meta.json simulations/*_issues.log`*

**Comparison analysis: what exists / what is missing**
- Exists: `lighter_protocol.py accept` (scores of lighter vs baseline seeds: within the baseline seed-to-seed std, needs >= 2 baseline seeds), `run_manifest.py compare <run_a> <run_b>` (what differs between two runs), nb8 `compare_*` figures and tables (`compare-da`, `compare-scenarios`, `compare-bias`, `compare-nens`), `lighten_run.py report`
- [x] `export_results.py OUT_DIR [--runs GLOB] [--only-ok] [--figures none|compare|all] [--zip]`: packages the files above + `overview.csv` (registry + disk + key scores per run) + `README.txt` for a colleague; copies, recomputes nothing (README section "Sharing the results")
- [ ] Table of 7.5 (option, disk / RAM / time saved, score change): a `lighter_protocol.py table` subcommand joining the measure json, the report and the summary.csv of baseline vs lighter. *Not implemented*
- [ ] Run baseline and lighter on 3 seeds each (`SEED=1 2 3`), otherwise `accept` cannot decide
- [ ] The PF run `SYN_NC_full_upd_pf_N32_s1_unbiased` had only 9 cycles when plotted (`INCLUDE_RUNNING=1`): rerun `plots` once finished, its scores are partial until then
- [ ] Check why the open loop `..._noDA_unbiased` is `failed` before using its `skill_vs_ol`

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

*Status 2026-11-24: no run is in progress any more, so the five renames of the synthetic-twin chain are done (below). The split of the plotters and the names still to decide are open, see "Still open".*

**Done**
- [x] `ligher_run.py` → `lighten_run.py` (`git mv`; README, tests, `lighter_protocol.py`, `export_results.py` updated; `grep -r ligher` finds only this note)
- [x] The synthetic-twin chain renamed with `git mv` (history kept), every reference updated in `README.md`, `run_twin.sh`, `run_matrix.py`, `nb4_DA_soil3D_nc.py` (comment), `synthetic_twin_plan.md`, the scripts' own docstrings and `tests/test_nb8_summarize.py` (import):

| Old | New |
|---|---|
| `nb1_cathy_soil3D_nc.py` | `forward_model_nc.py` |
| `nb6_synth_truth_ERT.py` | `synthetic_truth_ert.py` |
| `nb7_DA_ERT_synthetic.py` | `da_ert_synthetic.py` |
| `nb8_plots_synthetic.py` | `plots_synthetic.py` |
| `nb8_helpers.py` | `plots_synthetic_helpers.py` |

- [x] Old → new table added to the README ("Renamed files")
- *Checked: `git grep` finds no old name left outside this todo, the historical logs and `truth_meta_NC.json`; `bash -n run_twin.sh`, `py_compile` of all scripts, `pytest` (60 tests: 53 passed, 6 skipped, 1 failed = the already known `test_evolution_and_shift_flags` of section 2; `test_nb5_dryrun.py` and `test_scenario_nc_ERT.py` cannot be collected here, no pyCATHY). **No CATHY / pyGIMLi run was done after the rename**: do a `./run_twin.sh ... --dry-run` on the workstation first.*
- *Left as it is on purpose: the old names inside `notebook/logs/*.log`, `data/3d/synthetic/NC/truth_meta_NC.json`, run `manifest.json` / `run_meta.json` and the registry (historical record of what was run); the short names `nb6` / `nb7` / `nb8` in comments; the output folder `figures_nb8/` (a path written by the code, other tools read it); the test file names (`test_nb8_summarize.py`, `test_nb5_plot_helpers.py`).*

**Still open: names to decide (not proposed without reading what each one is used for)**

- `nb2_cathy_soil3D_hc.py`, `nb3_compare_meas_vs_mod.py`, `nb4_*`, `nb4a_*`
- `nb5_DA_ERT_soil3D_nc*.py`, `nb6_ESMDA_ERT_soil3D_nc.py`, `plot_6_ESMDA_synthetic.py`, `nb5_plot_helpers.py`
- `nc_model.py` / `model_config.py`, `ert_*.py`

**Why the plotters are not split yet**
- `plots_synthetic.py` (541 lines), `nb5_DA_ERT_soil3D_nc_plots.py` (597) are scripts with module-level state (`args`, `truth`, `RUN`, `scores_of`...) that the figure functions read as globals: splitting them into a `plotting/` package means passing that state explicitly, which needs a real `plots` run to verify.
- None of the `nbX` scripts can be run without CATHY + pyGIMLi, so only `grep` + `pytest` (stubs) can check a change here.

- [ ] **Plotters in separate files** (a `plotting/` package or `*_plots.py` modules)
  - Split the figure functions out of `plots_synthetic.py` (`fig_states`, `fig_profiles`, `fig_metrics`, `fig_params`, `compare_*`), `nb5_DA_ERT_soil3D_nc_plots.py`, `nb4_DA_soil3D_nc_plots.py`, `plot_6_ESMDA_synthetic.py`, `nb5_plot_helpers.py` / `plots_synthetic_helpers.py`
  - The `nbX` scripts keep only the CLI (arguments, loop over runs, calling the plotters)
- [ ] **Rename files whose name does not say what they do.** Candidates:
  - ~~`nb1`, `nb6_synth`, `nb7`, `nb8`~~ *done: `forward_model_nc.py`, `synthetic_truth_ert.py`, `da_ert_synthetic.py`, `plots_synthetic.py`*
  - still to do: `nb2`, `nb3`, `nb4`, `nb4a`, `nb5`, `nb6_ESMDA` (real-data and soil-moisture-only twin scripts)
  - `fig0_infiltr_compatt` / `fig1_infiltr_scheme` notebooks
  - `plot_6_ESMDA_synthetic.py`
  - ~~`ligher_run.py` (typo vs `lighten_run.py` in its docstring)~~ *done: `git mv` to `lighten_run.py`*
  - `nc_model.py` vs `model_config.py` (overlap)
  - `ert_data` / `ert_frame` / `ert_alignment` / `ert_schedule` (one `ert/` folder?)
- [ ] **When renaming:** `git mv` (keep history) and update every reference (`run_twin.sh`, `run_matrix.py`, `data_path.toml` `[scenarios]`/`[code]`, README, `synthetic_twin_plan.md`, tests, imports); check with `grep` and `pytest`; add the old → new table to the README
- [x] **Hygiene while at it:** `__pycache__/` and `.ipynb_checkpoints/` untracked and ignored (see section 2). Still open: tracked `notebook/cache/*.pkl`, `notebook/logs/*.log`, `duckdb.exe`

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
