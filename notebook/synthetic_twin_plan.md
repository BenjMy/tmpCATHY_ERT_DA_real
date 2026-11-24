# Synthetic twin experiment: ERT data assimilation with pyCATHY (NC and HC)

Date: 2026-10-23. Status: step 0 (config), nb6 (truth + synthetic ERT), library patch and scenario generator written and unit-tested; **NC chain run for real (nb6 -> nb7 -> nb8, section 15); HC and the parameter scenarios not run yet**; nb7, run_matrix and nb8 written and tested on stubs only (see section 14).

**Goal.** Use the hydrological models of nb1 (NC, uncompacted) and nb2 (HC, compacted) as the *truth*. Generate synthetic ERT data (5 % noise) with the real electrode sequence. Assimilate them with the nb5 scheme. Compare scenarios and DA types against the known truth (saturation `sw`, pressure head `psi`, parameters).

---

## 0. What the supplied files tell us (facts that shape the design)

| Source | Fact | Consequence |
|---|---|---|
| `st5t00.tldat` | 72 electrodes: 4 vertical lines at the corners of a **1 × 1 m** block (x,y ∈ {0,1}, z = -1.1 … 0, 12 each) + 24 surface electrodes (5×5 grid at 0.1…0.9 m, centre removed) | ERT footprint is only 1 × 1 m, centred on the 5 × 5 m CATHY plot (legacy shift (2.0, 2.0, +1e-3)). For HC it sees only the central part of the compacted band. |
| `st5t00.tldat` | 4790 quadrupoles; columns `a b m n err i ip iperr k r rhoa u valid`; ~40 % have `err = 999.99` (flagged bad); good rows have `err` ≈ 0.054-0.073 (a **fraction**, ~5-7 %) | Synthetic files must keep this format. Keep the same quadrupoles and the same bad-row mask so `filter_tldat` and `unreliable_quadrupoles` behave as with real data. Write both `r` and `rhoa` (the plot helper reads `rhoa`, falling back to `r`). |
| `st5t00.tldat` | `k` takes both signs; `r = rhoa / k` | Compute `r` from `rhoa_obs` and `k`; never add noise to `r` independently. |
| `site_params.json` | `INFILTRATION_START` NC 2023-10-09 11:00, HC 14:00; Ks experiment NC 2e-6, HC 5e-8 m/s | Reuse the same time reference (`ert_schedule.py`). |
| `ert_schedule.py` | `COMP_TO_LABEL` has `heavy_comp` and `no_comp` only; `check_schedule` requires strictly increasing times and CATHY t > 0 | OK for NC and HC. |
| `scenario_nc_ERT.toml` | VGP prior nominals **equal the truth** (0.44/0.38/0.40, 1.4/1.42/1.4, 0.04/0.065/0.07); IC prior -1.5 ± 0.5 vs truth [-2.4, -1.3, -1.15] | The VGP prior is currently unbiased. A twin experiment needs an explicit, documented prior bias (section 4). |
| `make_scenario_nc_VGP.py` | NC only; no Ks entry; no zones | Generalise for HC (2 zones) and add Ks. |
| `simulate_ERT.create_ERT_survey_pg` | `ert.simulate(..., noiseLevel=fwdNoiseLevel)`; `noiseAbs` only passed if `fwdNoiseAbs` is given | See issue I1: the DA forward may add unseeded absolute noise. |
| `nb2` vs `nb1` | Different cell size (0.1 vs 0.5 m), rain (0.04 vs 0.05 m), irrigation start (2 vs 3 h), IC | Globals in `nc_model.py` (`DELTA`, `IRR`, `IRR_START_HOUR`, `IC_TRUE`, …) must become per-case. |

---

## 1. File layout

| File | Role |
|---|---|
| `model_config.py` | Generalises `nc_model.py`: `CASES = {'NC': {...}, 'HC': {...}}`; `setup_model(simu, case, soil=…, ic_values=…, extra_times=…)`; ERT paths per case. Keeps all `nc_model` helpers (`build_times`, `snap_to_times`, `filter_tldat`, …). |
| `make_scenarios.py` | Replaces `make_scenario_nc_VGP.py`: writes one toml per case, including Ks, zones, Archie, layered IC. |
| `nb6_synth_truth_ERT.py` | Truth run, truth export, synthetic ERT `.tldat`, virtual TDR. |
| `nb7_DA_ERT_synthetic.py` | nb5 driver, parametrised (CLI), reads the synthetic files, writes `run_meta.json`. |
| `run_matrix.py` | Loops over the scenario matrix and the seeds. |
| `nb8_helpers.py` | Truth-aware helpers on top of `nb5_plot_helpers.py`. |
| `nb8_plots_synthetic.py` | Comparison and validation figures and tables. |

---

## 2. Step 0: `model_config.py`

1. One dict per case with: grid (`DELTA`, `X`, `Y`, `MAXDEPTH`), zonation (HC: central band), soil table per zone and layer (Ks, `POROS`, VGP), IC layers, `RAIN_TOTAL`, `IRR_START_HOUR`, ET, `T_END`, ERT paths, `COMP_LEVEL`.
2. Values copied from nb1 (NC) and nb2 (HC). nb1/nb2 stay as the reference runs; add a test that `model_config` reproduces their `sw` output (max abs difference < 1e-8) before using it.
3. `setup_model(..., soil_overrides=None)`: allows a biased prior soil table. The truth uses no override.
4. Truth parameters in one object (`TRUTH[case]`) that is pickled with the truth.
5. **HC mesh size.** At 0.1 m the HC mesh has 51×51×16 = 41 616 nodes (NC: 1 936). `df_DA` (already ~250 MB for NC) would be ~20× larger. Decide before running:
   - (a) 0.1 m with only a small ensemble, or
   - (b) coarser mesh with the band aligned to the grid (0.25 m → 2-cell band = 0.5 m wide; 0.2 m cannot centre a 0.4 m band), truth and DA on the **same** mesh.
   Recommendation: (b) for the matrix, (a) for one final run.
6. Layer indexing: the label slices (`1:6`, `7:10`, `10:`) overlap at layer 10; the later assignment wins, so materials are A = 0-5, B = 6-8, C = 9-14 (0-based), matching `MATERIAL_LAYERS`. Document it, do not change it.

---

## 3. Step 1: `nb6_synth_truth_ERT.py`

For each case (`NC`, `HC`):

1. **Time grid.** `ERT_TIMES_H` from `schedule2023_3d.csv` (as nb5), `ert_times` → `extra_times`. The truth and the DA must use the same `build_times(extra_times)` vector. Optional: denser synthetic schedule (e.g. every 6 h) as a separate experiment.
2. **Truth run.** `setup_model` with truth parameters and truth IC (layered). Output at every atmbc time.
3. **Save `truth_{case}.pkl`:**
   - `sw`, `psi` at all nodes and all times (array `[time, node]`), `atm_times`, node coordinates, mesh bounds;
   - truth parameters: Ks, VGP per layer and zone, `POROS`, layered IC profile, Archie a, m, n, rFluid;
   - random seed and noise settings.
   Node ordering is the CATHY node order; `df_DA` uses the position inside each (time, member) block, so both match if the mesh is identical (assert node count and coordinates).
4. **Resistivity of the truth.**
   - Same CATHY → pygimli mapping as the DA (reuse the pyCATHY petro/Archie path, same ERT↔CATHY frame transform, same alignment gate), with **truth** Archie parameters.
   - Archie is part of the experiment design (factor `archie_truth`): (i) truth = prior mean (a, m, n = 0.3, 1.7, 1.7), isolates hydrology; (ii) truth offset by about 1σ of the prior (e.g. 0.4, 1.85, 1.55), tests Archie error.
   - Stage 2: forward on a **finer** mesh than `meshForward.vtk` to avoid the inverse crime (ask how `meshForward.vtk` was built).
5. **Forward.** `ert.simulate(mesh, res=res0, scheme=scheme_from_st5t00, noiseLevel=0, noiseAbs=0)` with the boundary markers of `set_ert_boundary_markers`. Noise is added **by us**, not by pyGIMLi (pyGIMLi's own noise is `noiseLevel + noiseAbs/|u|`, i.e. much more than 5 % on weak signals).
6. **5 % noise.** `rhoa_obs = rhoa_true · exp(0.05 · ε)`, ε ~ N(0,1), seeded `RandomState(seed_noise)`. Log-normal noise keeps `rhoa > 0` and is consistent with `USE_LOG_OBS=True` (σ_log10 ≈ 0.0217). Store `rhoa_true`, `rhoa_obs`, seed.
7. **Write `.tldat`** with the real file as template (same header, electrodes, row order): `rhoa = rhoa_obs`, `r = rhoa_obs / k`, `err = 0.05` on the rows that are good in the real file, `err = 999.99` on the others (option `MASK=same|all_valid`). File names `st5t<ii>.tldat` (and `st1t<ii>.tldat` for HC) in a `synthetic/` folder, so `load_ert_files` only needs a different `data_dir`.
8. **Virtual TDR** at 15, 30, 50 cm at the plot centre: truth θ = sw · φ plus N(0, 0.02 m³/m³), written in the format of `st5.csv`/`st1.csv` (or a simplified csv + loader). Not assimilated: independent validation, as in nb5.
9. **Sanity prints:** `rhoa` range vs the real files (real data are mostly 10-100 Ω·m; with a = 0.3, m = n = 1.7 and φ ≈ 0.41, dry cells give rhoa of several 100 Ω·m, so check that the truth is not systematically out of the real range), count of rows per survey, share of good rows.

---

## 4. Step 2: `nb7_DA_ERT_synthetic.py`

Copy of nb5 with these changes only:

1. CLI (`argparse`): `--case`, `--scenario`, `--da-type`, `--nens`, `--seed`, `--inflation`, `--obs-err`, `--log-obs`, `--with-da`, `--prior-bias`.
2. `ert_paths` → synthetic folder; `ERT_DATA_ERR = 5` (also run 10 as an error-misspecification test).
3. Project name `SYN_{case}_{scenario}_{da}_{nens}_s{seed}[_OL]`: nb5's `RUN_TAG` does not contain NENS, seed or DA type, and `shutil.rmtree(da_dir)` would delete earlier runs.
4. `run_meta.json` in the project folder with all options (scenario, DA type, NENS, seed, inflation, obs error, log-obs, Archie truth, prior bias). The plot script groups on it, not on file names.
5. Move `PRIOR_CACHE` to `cache/`; the hash already contains scenario, NENS, SEED. Open loop and DA runs share the prior (same seed).
6. **Prior bias** (explicit, logged), applied through `soil_overrides` and the toml means:
   - VGP: SATCELL_VG × 1.3, n_VG − 0.08, thetar_VG + 0.01 relative to truth (about 1.2σ of the current priors: cv 25 %, σ_n ≈ 0.07, σ_θr = 0.02). Keep the truth inside ±2σ.
   - Ks: factor 5 on a log scale (σ_log10 ≈ 0.5).
   - IC: ~~homogeneous prior -1.5 ± 0.5 (truth is layered: structural error)~~ → per-layer prior, mean = `ic_true`, ± 0.5 (todo.md 5.1 (a)); `--ic-prior-mean -1.5` gives the old homogeneous prior.
   - A zero-bias reference run per scenario (prior = truth) to separate DA error from prior bias.
7. **Scenario toml** (`make_scenarios.py`): add
   - `Ks` (PERMX/Y/Z together, log-normal, per material; physical bounds): confirm the pyCATHY parameter name and write rule in `perturbate` (the commented plot call in nb5 uses `parm='Ks'`);
   - `ic_layered` (3 layers, 'layers' control) so the IC can match the truth structure;
   - HC: zone-aware parameters (compacted zone 2: layers 0-5 with their own Ks, φ, VGP). Check in `describe_ensemble` which zone each key writes to.

---

## 5. Step 3: scenario matrix (`run_matrix.py`)

| Factor | Levels |
|---|---|
| Case | NC, HC |
| Updated parameters | ic (states only); ic + VGP-material (9); ic + VGP-layers (45); ic + Ks; ic + VGP + Ks; ic + Archie (a, m, n); ic + VGP + Ks + Archie |
| DA type | EnKF (inflation 1.0, 1.1), PF (tempered, Liu-West) |
| Ensemble | NENS 64, 100 |
| Obs | `USE_LOG_OBS` True/False; obs error 5 % (correct) vs 10 % |
| Archie truth | = prior mean / offset |
| Prior | zero bias / biased |
| Seeds | ≥ 3 (one run is anecdotal) |
| Open loop | one per (case, scenario, seed) |

Do not run the full product. Order: (1) ic-only EnKF with and without bias, NC; (2) + VGP-material; (3) + Ks; (4) PF; (5) HC with the winners; (6) seeds. Each run writes `run_meta.json`, `_ERT_checks.json`, `_rejections.json`.

---

## 6. Step 4: `nb8_plots_synthetic.py` (+ `nb8_helpers.py`)

Reuse `nb5_plot_helpers.py` (ensemble stats, probes, profiles, `parm_evolution`, `parm_shift_table`, rejections). Add:

1. **States vs truth.** At the probe nodes (15/30/50 cm), `sw` and `psi`: ensemble median + 5-95 % band + member lines, truth as a black line, virtual TDR as dots, ERT survey times as dotted vertical lines. Use `sw_bef_update_` (forecast) and the analysis separately. Truth aligned to `atm_times[c]` (state at the END of window c).
2. **Profiles at survey times** (nb5 `profiles_figure`) with the truth profile.
3. **Metrics against the truth, all nodes** (ERT footprint and whole domain separately, per layer, per cycle): bias, RMSE (`sw`; `psi` and log10(-psi)), spread/RMSE ratio, coverage of the 5-95 % band, CRPS, rank histogram. Table `metrics_state.csv`.
4. **Parameters vs truth.** For each parameter family (ic, SATCELL_VG, n_VG, thetar_VG, Ks, a/m/n): prior and analysed ensemble over cycles, true value as a horizontal line; table `(mean - truth)/prior_sd`, posterior/prior sd ratio, flags SHIFT/COLLAPSE (existing `parm_shift_table`) plus `|z| > 2` = TRUE-MISS. For `ic` compare the scalar with the depth-weighted truth profile and show the structural error.
5. **ERT fit.** χ² (should be ≈ 1) and normalised RMS per cycle from `df_performance`; share of rejected members per cycle/stage; number of quadrupoles kept per survey.
6. **Comparison panels.** Scenarios as columns, DA types as line styles, open loop in grey; one heat map (rows = scenario, columns = case/seed) of final state RMSE and of mean |z| of parameters; box plots over seeds.
7. **Virtual TDR skill** (bias, RMSE, r, coverage) as `smc_skill`, now also against the truth (checks that the TDR noise model is not the limiting factor).
8. Each section stays inside `safe(...)` so one failing figure does not stop the rest.

---

## 7. Checks to build in

1. Perfect-twin test: no noise, zero prior bias, truth in the ensemble → RMSE and χ² ≈ 0.
2. Truth forward with a homogeneous medium: `rhoa/rho = 1` ± 10 % on the kept quadrupoles (reuse `unreliable_quadrupoles`).
3. Truth lies inside every `per_bounds` entry; truth VGP within ±2σ of the biased prior.
4. Node count/coordinates identical between truth and DA mesh; same time vector.
5. Synthetic `rhoa` range compared with real `rhoa` (section 3.9).
6. Unit test of `model_config` against nb1/nb2 output.
7. Memory estimate (nodes × NENS × cycles × columns) before each HC run.

---

## 8. Issues found while reading the code

| ID | Where | Issue | Suggested fix |
|---|---|---|---|
| I1 | nb5 / `simulate_ERT.py` | `ERT_NOISE_FLOOR_V = None` → `fwdNoiseAbs` not passed → `ert.simulate` keeps its default `noiseAbs` (1e-4 V in recent pyGIMLi). With `fwdNoiseLevel = 0` this may still add unseeded absolute noise to every ensemble prediction (verify in your pyGIMLi version). `unreliable_quadrupoles` already passes `noiseAbs=0`. | Always pass `noiseAbs = 0` from the DA (new key or default 0 in `create_ERT_survey_pg`). |
| I2 | `nc_model.ert_paths` | `sequence` hard-coded to `seq3d_from_st5t00tldat.shm` for every comp level: HC would silently use the NC electrodes. | Per-case sequence file; assert electrodes equal the first lines of the data files. |
| I3 | nb5 | `da_name` ignores NENS, seed, DA type; `rmtree` deletes the previous run. | Unique run name (section 4.3). |
| I4 | nb5 | `PRIOR_CACHE` written to the current directory. | `cache/` folder. |
| I5 | nb2 | `n_compacted_layers` (`searchsorted(..., 'right') + 1`) over-counts by one and is not used; the compacted layers are actually 1:6 (0-20 cm) although the comment says ~30 cm. | Decide the true depth, set the slice from it, delete the print. |
| I6 | nb1/nb2/`nc_model` | Overlapping label slices at layer 10 (harmless, see 2.6). | Use disjoint slices or comment it. |
| I7 | nb1 vs nb2 | `read_outputs('sw')` returns `(df, time)`, `read_outputs('psi')` only a `df` in the scripts: check against the source (needed for the truth export). | Verify in pyCATHY. |
| I8 | `make_scenario_nc_VGP.py` | VGP drawn independently (existing FIXME). In the twin the truth VGP are exact material constants, so correlation matters less, but keep the FIXME. | Later: joint draw (Botto et al. 2018). |

---

## 9. Open points / files still needed

1. Header (first 3 rows) of `schedule2023_3d.csv`.
2. HC: ERT data folder layout (`st1HC/procdata/st1t<ii>.tldat`), the electrode sequence file for HC (same as NC?), and the forward mesh (same `meshForward.vtk`?).
3. How `meshForward.vtk` was generated (needed for the finer truth mesh).
4. pyCATHY source at commit `083f0ac` to verify: `read_outputs` return types, `perturbate` Ks naming and zone handling, the petro/Archie mapping function (to reuse for the truth), columns read by `read_observations`, `df_DA` node ordering, pygimli noise defaults. Repository to clone: please confirm the owner/name (I assume `BenjMy/pycathy_wrapper`).
5. Decision on the HC mesh (section 2.5) and on the prior bias magnitudes (section 4.6).

---

## 10. Order of work

1. `model_config.py` + reproduction test of nb1/nb2.
2. `nb6` for NC (truth, synthetic ERT, virtual TDR) + checks of section 7.
3. `nb7` for NC with `ic` only (zero bias → biased).
4. `nb8` core figures (states vs truth, metrics, ERT fit).
5. Parameters: VGP-material, Ks, Archie; PF.
6. HC (after answers to section 9).
7. Seeds and matrix, summary figures.


---

## 11. What has been done (changes)

All new files are in `twin/`; library changes are in `pycathy_twin_fixes.patch` (git diff against pyCATHY commit 526056a, "refactor pf"; applies cleanly to a fresh checkout).

### 11.1 Twin files

| File | Content | Verified |
|---|---|---|
| `model_config.py` | `CASES` NC/HC copied from nb1/nb2 (soil tables, IC, band, irrigation, rain, ET); `setup_model`, `biased_overrides`, `truth_parameters`, `sensor_phi`, `ert_paths(synthetic=True)`; re-exports the `nc_model` helpers | `test_model_config.py` (tables equal nb1/nb2, band, irrigation, bias, sensor porosity) |
| `synth_tools.py` | `.tldat` read/write keeping the real file layout, seeded multiplicative log-normal noise, Archie, time lookup | `test_synth_tools.py` on the real `st5t00.tldat` (2837 good rows, sigma_log10 ~ 0.0217 for 5 %, round trip, accepted by `filter_tldat`) |
| `ert_frame.py` | ERT <-> CATHY frame choice factored out of nb5 | compiled |
| `nb6_synth_truth_ERT.py` | truth run -> sw/psi/nodal porosity -> Archie -> `trace_mesh` -> pyGIMLi forward -> 5 % noise -> `.tldat` with the real quadrupoles; virtual TDR; `truth_<case>.pkl` + json | py_compile only (needs CATHY + pyGIMLi) |
| `make_scenarios.py` | 7 scenarios per case (ic, ic+VGP, +VGP upd, +Ks upd, Archie upd, VGP+Ks, full) with `material` (NC) / `zone_material` (HC); biased or unbiased prior; `Kz_ratio` stub (off) | `test_make_scenarios.py` (all 14 pass the library validator) |

### 11.2 Library changes (`pycathy_twin_fixes.patch`)

| # | File | Change | Why |
|---|---|---|---|
| L1 | `ERT/petro_Archie.py` | new `nodal_porosity()`; `SW_2_ER0_DA` passes one porosity per node (ensemble container, nodal list, scalar, with a size check) | **Bug**: `Archie_rho_DA` indexed `porosity[0]` (node 0) for every node, so the nodal porosity map was ignored. |
| L2 | `DA/cathy_DA.py` | porosity log line reads the first valid member; console line handles `zone_material` | follow-up of L1 / L3 |
| L3 | `DA/param_tables.py`, `DA/perturbate.py`, `DA/scenario.py` | new control `zone_material` (`zone_material_layers = [[zone, [layers]], ...]`), validated | `layers` / `material` write into **all zones**, so in HC the prior of zone 1 overwrote the compacted band. |
| L4 | `DA/perturbate.py`, `DA/scenario.py` | **Ks** goes through the shared builder (`process_perturbation_block`): all controls, log / log-ratio transforms, per-layer bounds; `Ks` added to `MATERIAL_PARAMETERS` | Ks was a special branch without these. |
| L5 | `DA/param_tables.py`, `DA/cathy_DA.py` | Ks "c-lite": `PERMX = PERMY = Ks`, `PERMZ = Ks x base kz/kx` of each soil row (ratio read from the member soil file before it is emptied, kept in `df.attrs["kz_ratio"]`) | HC band has kz/kx ~ 0.017; the old code forced it to 1 whenever Ks was perturbed. |
| T | `tests/` | `test_archie_porosity_nodal.py`, `test_zone_material_control.py`, `test_ks_anisotropy_clite.py` | all pass |

**Regression check** (full suite, `--continue-on-collection-errors`): baseline 10 failed / 284 passed / 15 collection errors; with the patch 10 failed / 292 passed (+8 new) / 15 errors; re-run after adding c-lite (4 more tests): identical failing set. The pre-existing failures are missing files/modules (`scenario.toml`, `pyCATHY.ERT.ert_data`, `utils_CATHY`, `ERT_checks`) and are unrelated.

### 11.3 Design decisions taken

- **Truth = continuous run, DA = per-window restarts** (as in nb5).
- **Observation time**: obs stamped with the start of the model window and compared to the state at the window end (`obs_window_start`), unchanged from nb5.
- **Noise**: only in the data (log-normal, 5 %, seeded, seed + 1000 x survey index); pyGIMLi noise is disabled (`noiseLevel = 0`, `noiseAbs = 0`). The DA uses `ERT_DATA_ERR = 5 %` on log(rhoa).
- **Quadrupoles**: those of the real files (valid rows), so the synthetic data have the real geometry; rows with rhoa <= 0.01 after the forward are dropped.
- **Prior bias** (to avoid a trivially perfect prior): Ks x5, SATCELL x1.3, n -0.08, thetar +0.01, IC: per layer, centred on the layered truth (was -1.5 +/- 0.5, todo.md 5.1 (a)). A zero-bias reference run is kept.
- **Archie truth**: prior (a, m, n) = (0.3, 1.7, 1.7), or `--archie offset` = (0.4, 1.85, 1.55) so the DA has a petrophysical error to recover; rFluid = 1/(588e-4) ~ 17 ohm.m.
- **Ks anisotropy**: option (c-lite) chosen, i.e. the truth stays equal to nb2 and the prior keeps the base kz/kx per row; a separate `Kz_ratio` parameter is an extension (section 13).

---

## 12. Limits (what is not covered or not verified)

### 12.1 Not verified

1. **nb6 and the patch were never run with real CATHY / pyGIMLi** (not installed in the working environment). Unit tests cover the pure-Python parts only. First run = debugging run; expect small path / API fixes.
2. HC electrode sequence: electrodes verified equal to NC; quadrupoles come from the HC data files (see "Decisions taken").
3. `schedule2023_3d.csv` (first rows) never seen: the ERT times are assumed to be the ones of nb5.
4. Ks sigma / bounds in `make_scenarios.py` (log transform, sigma = nominal, bounds nominal /100 ... x100): check on the first run that the ensemble spans about one decade.

### 12.2 Structural limits of the experiment

| Limit | Effect | Mitigation |
|---|---|---|
| Truth and DA forward use the **same ERT mesh** (`meshForward.vtk`) unless `--fwd-mesh` | inverse crime on discretisation, optimistic ERT fit | give a finer truth mesh (needs to know how `meshForward.vtk` was built) |
| Truth and DA share Archie **form**; with `--archie prior` also the parameters | optimistic petrophysics | use `--archie offset` for at least one scenario |
| **Multi-zone Archie unsupported** in the library (`get_Archie_ens_i` reads `len(rFluid) > 1` as "per member"; `Archie_rho_DA` loops zones without a mask) | one (a, m, n, rFluid) for the whole domain, also for the compacted band | truth uses one set too, so the twin is consistent; real-data runs would need it (est. half a day) |
| Ks isotropy | `PERMZ = Ks x base ratio`: the **ratio itself is not estimated** | extension `Kz_ratio` (section 13) |
| VGP drawn independently (FIXME Botto et al. 2018) | unrealistic prior correlations, EnKF may update parameters independently | joint draw |
| HC zone 2 below layer 5 keeps library defaults (as in nb2) | faithful to nb2 but physically odd; not controlled by `zone_material` | option `zone2_fill_below = True` copies zone 1 (decision needed) |
| Porosity at zone / layer interfaces is a cell -> node average | small smoothing of the truth porosity | none needed; document |
| HC mesh: 0.1 m = 41.6k nodes (20x NC); 0.25 m gives a 0.5 m band | memory (`df_DA`), band width mismatch with nb2 | decision (section 2.5) |
| ERT footprint 1 x 1 m block (72 electrodes) | sw/psi outside the block are constrained only through the model, not the data | validation metrics must be reported **inside** and **outside** the block |
| Real-data features not simulated | electrode contact errors, temperature, borehole effects, systematic errors | add a "model error" scenario (e.g. noise 10 %, 1 % systematic bias) |

### 12.3 Open issues from section 8 (status)

I1 (`noiseAbs` default): mitigated in nb6, **must be set to 0 in nb7**. I2 (hard-coded NC sequence): fixed in `model_config.ert_paths`, HC file still to supply. I3, I4: to do in nb7 (unique run name, `cache/`). I5, I6: documented, not changed in nb1/nb2. I7: `read_outputs` return types handled in nb6 defensively, to confirm at the first run. I8: unchanged.

---

## 13. Perspectives of improvement

### 13.1 Ks anisotropy: should the kz/kx ratio be perturbed too?

Currently (c-lite) Ks is perturbed, the ratio is fixed to the base table. If the ratio is also uncertain (it is, in compacted soil: ~0.017 measured vs 1 in an uncompacted soil):

- add a parameter `Kz_ratio` (log-normal, same groups as Ks), with `PERMX = PERMY = Ks`, `PERMZ = Ks x Kz_ratio`; the soil table must be assembled from both parameters before it is written (about half a day, library change in `param_tables` / `perturbate` / `scenario`);
- staging: main matrix with the fixed ratio; **one** HC scenario with a perturbed and deliberately biased ratio (e.g. prior 0.1 vs truth 0.017) to test recoverability;
- expectation: ERT constrains the ratio weakly (it acts on vertical flux through the band), so the posterior may stay close to the prior; the twin experiment is precisely what tells this. Keep Ks and ratio in log space and use bounds.
- `make_scenarios.py --kz-ratio` is a placeholder; it only writes a comment until the library supports it.

### 13.2 Library

1. Multi-zone Archie (mask by zone, per-zone rFluid / a / m / n; distinguish "per member" from "per zone" lists).
2. Joint VGP draw per material (multivariate normal in transformed space, optional shared draw across layers of a material).
3. Make `PERMZ` handling explicit in the scenario (`kz_ratio` key) rather than via `df.attrs`.
4. `create_ERT_survey_pg`: read `.vtk` before `.msh` (nb6 currently aborts if `.msh` / `.bms` siblings exist).
5. PR-ready branch from `pycathy_twin_fixes.patch` (cannot push to `BenjMy/pycathy_wrapper` from here; send the patch or give a fork).

### 13.3 Experiment design

1. Denser ERT schedule (`--dense-step-h`) to test information content vs cost.
2. Ensemble-size and seed matrix (e.g. NENS 32/64/128 x 5 seeds) with box plots; EnKF vs PF at equal cost.
3. Model-error scenarios: wrong IC structure, wrong irrigation timing, wrong Archie form, finer truth mesh.
4. Joint assimilation of ERT + virtual TDR (already written by nb6) vs ERT only.
5. Localisation / inflation sweeps (`inflation`, obs error 3-10 %).
6. Final HC run at 0.1 m with a small ensemble.

### 13.4 Metrics and plots (nb8)

States vs truth at probes and profiles (RMSE, bias), spread-skill ratio, 90 % coverage, CRPS; parameters vs truth with z-scores and prior-to-posterior contraction; ERT chi^2 per survey; heat maps of scenario x DA type; box plots over seeds; virtual TDR skill (independent validation, not assimilated).

---

## 14. Order of work (status)

| # | Item | Status |
|---|---|---|
| 1 | `model_config.py` + tests | done (reproduces nb1/nb2 tables; the `sw` reproduction test needs CATHY) |
| 2 | `synth_tools.py`, `ert_frame.py`, nb6 | written, nb6 not run |
| 3 | Library patch L1-L5 + tests | done, no regression |
| 4 | `make_scenarios.py` | done |
| 5 | `nb7_DA_ERT_synthetic.py` | written (py_compile + pyflakes clean, `--help` works, not run): nb5 with CLI, case config, synthetic ERT (times read from `truth_<case>.pkl`), generated scenarios, biased prior, `ERT_DATA_ERR = 5`, `fwdNoiseAbs = 0`, unique run name `SYN_<case>_<scenario>_<da>_N<nens>_s<seed>`, `cache/`, `run_meta.json` (config, checks, status, traceback), asserts that DA and truth share atm times and delta |
| 6 | `run_matrix.py`, `twin_names.py` | written: presets quick / main / ens, `--open-loop` adds the no-DA reference, resumes (skips `finished` runs), one log per run, `--dry-run` tested |
| 7 | `nb8_helpers.py`, `nb8_plots_synthetic.py` | written; helpers unit-tested (`test_nb8_helpers.py`), whole script tested end to end on fake truth + fake runs with a stub DA class (`test_nb8_e2e.py`); real df_DA / dict_parm_pert layout still to confirm at the first run. Per run: states at the sensor nodes (+ virtual TDR), profiles, metrics in time (all / inside / outside the ERT block: RMSE, spread-skill, 90 % coverage, CRPS), parameters vs truth. Across runs: summary.csv, heat maps scenario x DA, box plots over seeds, parameter z / contraction / coverage, ensemble-size study |
| 8 | First real run NC (nb6 -> nb7 -> nb8), then HC | **next, needs CATHY + pyGIMLi**: `python nb6_synth_truth_ERT.py --case NC`, then `python run_matrix.py --preset quick --dry-run`, then without `--dry-run`, then `python nb8_plots_synthetic.py --case NC` |
| 9 | `Kz_ratio` extension, multi-zone Archie | optional |

### Decisions taken (23 Oct 2026, from your answers)

1. **HC zone 2 below layer 5**: as nb2 (library defaults kept, `zone2_fill_below = False`).
2. **HC mesh**: the question was about the CATHY grid (cell size 0.1 m in nb2 vs 0.5 m in nb1; the ERT mesh is separate). Default = **0.1 m as in nb2** (41.6k nodes, heavy for `df_DA`: 20x NC). `--delta` stays available to run the matrix on a coarser grid (0.25 m gives a 0.5 m band instead of 0.4 m) if memory becomes the limit.
3. **HC electrodes / quadrupoles** (`st1t00.tldat`): same 72 electrodes as NC (so the NC `.shm` is valid for electrodes and `meshForward`); **different quadrupoles** (4911 vs 4790, 4771 in common, 2203 good rows vs 2837, good rhoa 2.1-298 ohm.m). nb6 therefore takes the quadrupoles of each case's own data files. HC ERT times come from the `1-HC` column of `schedule2023_3d.csv` (23 surveys, first 09/10/2023 12:54, i.e. **before** the HC irrigation start at 14:00; NC first survey 10:34 vs 11:00).
4. **Fine truth mesh**: `make_forward_mesh.py` builds `meshForward` (region 2 = electrode block, max cell volume `h_reg2`, optional `n_h2` H2 refinement). Truth mesh = same call with one more refinement, e.g. `python make_forward_mesh.py seq3d_from_st5t00tldat.shm meshForwardTruth 0.8 0.5 2.0 0.1 1.8 1 0.15` (8x more cells in the block), passed to nb6 with `--fwd-mesh`.
5. **Archie truth**: prior values (a, m, n) = (0.3, 1.7, 1.7) (`--archie prior`, default).

### Still open

- Confirm the HC ERT survey count and indices used by nb5 (the schedule has 23 rows for each campaign).


---

## 15. First real run (NC, 23-24 Oct 2026): what it showed

Environment of the test: pyGIMLi 1.6.1 + `tetgen` + `pyarrow` (pip), `gfortran`, `liblapack-dev`, `libblas-dev` (apt); CATHY sources fetched and compiled by pyCATHY. Data tree rebuilt from the files you sent (`st5t00` / `st1t00`, schedule, `site_params.json`); **the forward mesh was a structured hex substitute** (0.1 m grid around the array, 10 816 cells; `make_forward_mesh.py` failed here because pyGIMLi's tetgen wrapper could not read its `.poly` file), so the real `meshForward.vtk` must be used for your runs.

Run: `nb6 --case NC` (about 3 min), `nb7 --case NC --scenario ic --nens 8` (22 surveys, about 25 min), `nb8 --case NC`. All three completed.

### 15.1 Bugs found by the real run and fixed

| # | Where | Problem | Fix |
|---|---|---|---|
| R1 | nb6 | **A continuous CATHY run does not print at the atmbc times** (3627.9 s printed for a request at 3600 s: first step after the request), so truth and DA states were not at the same times; the old assertion failed. | Dense output (`--out-dt`, default 300 s, integer seconds because `update_parm` casts to int) and **linear interpolation in time to the exact atmbc times** (`synth_tools.interp_in_time`; largest CATHY step 1000 s). `truth['sw_time']` = atmbc times. |
| R2 | nb6 | Near-zero-potential quadrupoles gave rhoa far outside the model range (0.3 ohm.m for a model of 23-190 ohm.m, up to 12 rows per survey). | Dropped when outside [rho_min / r, rho_max x r], `--rhoa-ratio` (default 3, 0 = keep). 1-6 rows per survey are dropped, 17 of 23 surveys keep all 2837. |
| R3 | nb7 | `cache/` relative path broke because CATHY changes the working directory. | absolute path. |
| R4 | env | `.shm` first line must be the electrode count alone. | (test data only) |

### 15.2 Results of the smoke run (not a scientific result)

`ic` scenario, biased prior (Ks x5, SATCELL x1.3, ...), 8 members: RMSE sw 0.062 over all nodes, 0.104 inside the ERT block, 90 % coverage 0.2 %. The DA trajectory drains faster than the truth (Ks x5) and ERT cannot correct it because only the initial condition is perturbed: **a structural bias with an under-dispersed ensemble**. That is the expected behaviour and the reason for the `ic_VGP_upd`, `ic_Ks_upd` and `full_upd` scenarios.

### 15.3 Still to do after this run

1. Run `ic_VGP_upd`, `ic_Ks_upd`, `full_upd` and the open loop (`--no-da`): check that the backup of `dict_parm_pert` is written and read by nb8 (the `ic`-only run has none, so `nb8` could not exercise the parameter figures on real output).
2. HC run (nb6 HC at 0.1 m, then nb7): zone / `zone_material` / Ks c-lite paths are tested on unit tests only.
3. Re-run with the real `meshForward.vtk` (and a finer truth mesh for the final experiments).
4. Expected runtime: nb7 ~ 25 min for NENS = 8 on NC with 2 cores; NENS = 100 needs a multi-core machine (`parallel=True`) and HC is about 20x heavier in memory.

---

## 16. Memory cost of the DA (24 Oct 2026)

Measured on the real NC run (8 members, 362 cycles, 1936 nodes): `df_DA` = 5.6M rows x 7 float64 = **314 MB** (5 of the 7 columns repeat per node); rebuilt by `pd.concat` at every member and every cycle; NENS = 100 -> ~4 GB, HC at 0.1 m -> ~85 GB. Also `backup_results_DA` wrote the whole table again (parquet) every cycle.

Changes (library, `pycathy_twin_fixes`):

| # | Change |
|---|---|
| M1 | `DA/state_store.py` (compact store, float32, one chunk per cycle) is **wired into the DA class**: `DAConfig(compact_states=True, store_full='all'|'analysis', store_nodes=[...])`; `_DA_df` appends to the store (no concat); `df_DA` is a property rebuilding the legacy long table on demand (compact dtypes + `idnode`); `<project>_states/cycle_*.npz` written per cycle (O(1)); `<project>_df_DA.pkl` written ONCE at the end; `load_backup` rebuilds `df_DA` from the chunks if a run was interrupted. |
| M2 | **Bug**: `backup_results_DA` raised `KeyError: 0` on `df_obs['data_cov'][0]` (MultiIndex) and the `except` skipped ALL attributes, so `dict_parm_pert`, `df_obs`, `df_Archie` were never saved (nb8 could not draw the parameters). Now one try per attribute, `.iloc[0]`, `df_DA` not rewritten with the store. |
| M3 | nb7: `--store-full analysis --store-stride 10` keeps every node only at the survey cycles, probes + profile + a sample elsewhere (HC). nb8 reads `idnode` and uses the nodes present at every cycle. |

Measured alternatives for the same content: long float32 table, parquet zstd 65 MB; arrays [cycle, member, node] float32 compressed 35 MB; sw as uint16 22 MB (vs 314 MB). Tests: `tests/test_compact_states.py` (7); full suite: same failing tests as before (10), 304 passed.

Not done: delete member vtk files during the run (`DAConfig.cleanup_vtk` only acts at the end), `VTKF = 0` for the truth run (1.4 GB of vtk, untested whether `read_outputs` needs them), writing `sw` as uint16.
