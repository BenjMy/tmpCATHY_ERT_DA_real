# tmpCATHY_ERT_DA_real

Hydrological simulation of the Agripolis infiltration test (UNIPD, Oct 2023) with
[pyCATHY](https://github.com/BenjMy/pycathy_wrapper) and **data assimilation (DA) of 3D ERT data**, in two flavours:

| Flavour | Observations | Truth known? | Scripts |
|---|---|---|---|
| **Real data** | measured ERT (`st5t<ii>.tldat`, `st1t<ii>.tldat`) | no | `nb5_DA_ERT_soil3D_nc.py` (EnKF / PF), `nb6_ESMDA_ERT_soil3D_nc.py` (ES-MDA), `nb5_DA_ERT_soil3D_nc_plots.py` |
| **Synthetic twin** | ERT simulated from a known truth (5 % noise) | yes | `nb6_synth_truth_ERT.py` -> `nb7_DA_ERT_synthetic.py` -> `nb8_plots_synthetic.py` |

Two soil cases are modelled: **NC** (uncompacted, `nb1`) and **HC** (heavily compacted band, `nb2`).
Everything is driven from `notebook/`; `notebook/run_twin.sh` is the single entry point.

> The scripts need pyCATHY (with the twin patches), CATHY (Fortran, compiled by pyCATHY) and pyGIMLi.
> Parts of this README (paths, run registry) describe code that was tested without CATHY/pyGIMLi
> (unit tests and stubs only); see [Status and limits](#status-and-limits).

---

## 1. Repository layout

```
data_path.toml                 <- ALL paths (inputs, outputs, plots, cache, logs, code lib)
data/
  site_params.json             infiltration start per experiment, Ks, van Genuchten reference values
  sensors/                     TDR / soil moisture series (st1.csv = HC, st3.csv = LC, st5.csv = NC)
  3d/
    schedule2023_3d.csv        acquisition time of every ERT survey (columns 1-HC, 3-LC, 5-NC)
    seq3d_from_st5t00tldat.shm electrode positions + measurement sequence (72 electrodes)
    mesh_Luca_Agripolis3d/     meshForward.vtk = pyGIMLi forward mesh used by the DA
    st5NC/procdata/            real NC ERT: st5t00.tldat, st5t01.tldat, ...
    st1HC/procdata/            real HC ERT: st1t00.tldat, ...
    synthetic/<CASE>/          written by nb6: truth_<CASE>.pkl, procdata/*.tldat, tdr_synth_<CASE>.csv
  2d/                          2D time-lapse inversions (pre-processing, not used by the DA)
notebook/
  project_paths.py             reads data_path.toml (every script imports it)
  run_registry.py, run_manifest.py  run registry (csv) and per-run reproducibility manifest
  nb1_*, nb2_*, nb3_*          forward models (NC, HC) and measurement-vs-model comparison
  nb4*, nb4a*                  twin experiment on soil moisture (TDR) only
  nb5_*, nb6_ESMDA_*           real-data ERT DA
  nb6_synth_*, nb7_*, nb8_*    synthetic ERT twin: truth -> DA -> scores/figures
  model_config.py, nc_model.py case configuration (soil tables, grid, IC, ERT paths) and helpers
  make_scenarios.py            generates the synthetic scenarios (NC + HC)
  make_scenario_nc_VGP.py      generates scenario_nc_ERT.toml (real-data scenarios)
  scenario_nc_ERT.toml         real-data scenarios (generated, do not edit by hand)
  run_matrix.py, run_twin.sh   batch drivers
  ert_*.py, synth_tools.py     ERT schedule / frame alignment / .tldat I/O / noise
  tests/                       unit tests (pytest)
  synthetic_twin_plan.md       design notes, decisions and known limits of the twin experiment
simulations/                   (git-ignored) one folder per run
registry/runs_registry.csv     one row per run: what has been done or not
```

---

## 2. Paths: `data_path.toml`

Every folder the code reads or writes is defined in **`data_path.toml`** (repository root) and resolved by
`notebook/project_paths.py`. Relative paths are relative to `[project].root`, which is relative to the toml
itself; absolute paths, `~` and `$ENV` work everywhere. An empty string means "derive from the parent folder".

| Section | Keys |
|---|---|
| `[project]` | `root` |
| `[inputs]` | `data`, `data_3d`, `site_params`, `ert_schedule`, `sequence_dir`, `mesh_dir`, `mesh_forward`, `sensors_dir` |
| `[scenarios]` | `real` (scenario_nc_ERT.toml), `synthetic` (scenarios_twin.toml) |
| `[outputs]` | `synthetic_dir`, `simulations`, `figures_real`, `figures_synth`, `cache`, `logs`, `run_registry` |
| `[code]` | `notebook`, `pycathy` (a clone of `pycathy_wrapper`; empty = installed package) |

```bash
python notebook/project_paths.py                    # print every resolved path and whether it exists
DATA_PATH_TOML=/other/data_path.toml ./notebook/run_twin.sh status   # use another file (all scripts)
```

Precedence: command-line flags (`--root`, `--sim-dir`) and the shell variables of `run_twin.sh`
(`ROOT`, `SIM_DIR`, `SYN_DIR`, `PYCATHY`) override the toml; the toml overrides the built-in defaults, which
reproduce the historical layout (scripts run from `notebook/`, root = `../`). A typo in a key name is an error.

---

## 3. Inputs

| Input | Format | Where | Used by |
|---|---|---|---|
| ERT surveys | pyGIMLi `.tldat` (a b m n err i ip iperr k r rhoa u valid; `err = 999.99` = flagged bad) | `data/3d/st5NC/procdata`, `st1HC/procdata` | nb5, nb6 (real), template of nb6 (synthetic) |
| Electrodes / sequence | `.shm` (first line = number of electrodes) | `[inputs].sequence_dir` | all ERT scripts |
| Forward mesh | pyGIMLi `.vtk` (region 2 = electrode block) | `[inputs].mesh_dir/mesh_forward` | nb5-nb7 (nb6 truth: optional finer `--fwd-mesh`, see `make_forward_mesh.py`) |
| Survey times | `schedule2023_3d.csv` + `INFILTRATION_START` in `site_params.json` | `[inputs]` | `ert_schedule.py` (hours since irrigation) |
| TDR / soil moisture | csv | `[inputs].sensors_dir` | nb3, plots (independent validation, never assimilated) |
| Case definition | python dict (`CASES['NC' / 'HC']`: grid, soil table per zone/layer, IC, rain, irrigation, ET) | `model_config.py` | all models |
| Scenario | toml (section 6) | `[scenarios]` | nb5, nb6 (real), nb7 (synthetic) |

---

## 4. Pipelines

### 4.1 Forward models and checks
- `nb1_cathy_soil3D_nc.py` (NC) and `nb2_cathy_soil3D_hc.py` (HC): CATHY infiltration models, reference runs.
- `nb3_compare_meas_vs_mod.py`: measured TDR vs model.
- `nb4a_DA_solution_soil3D_nc.py` -> `nb4_DA_soil3D_nc.py` (+ `_plots`): small twin experiment assimilating synthetic soil moisture only.

### 4.2 Real-data ERT DA (`nb5`, `nb6_ESMDA`)
Biased prior ensemble (from a scenario) -> CATHY forecast per ERT window -> saturation to resistivity (Archie)
-> pyGIMLi forward on `meshForward.vtk` re-positioned on the CATHY mesh -> update with the measured
`log10(rhoa)`. `nb5` is sequential (EnKF or particle filter, `REAL_USE_PF=1`); `nb6_ESMDA` is the smoother
(ES-MDA, parameters only, experimental). `WITH_DA = False` (`REAL_WITH_DA=0`) gives the **open loop**
(data error 1e99, nothing is assimilated).

### 4.3 Synthetic ERT DA pipeline (twin experiment)

```
make_scenarios.py ──► scenarios (NC_*, HC_*)                              ┐
model_config.py   ──► truth parameters / biased prior soil table          │
                                                                           ▼
nb6_synth_truth_ERT.py   truth run ─► sw, psi ─► Archie ─► pyGIMLi forward ─► 5 % noise ─► .tldat
        └─ writes data/3d/synthetic/<CASE>/{truth_<CASE>.pkl, procdata/*.tldat, tdr_synth_<CASE>.csv}
                                                                           │
nb7_DA_ERT_synthetic.py  prior ensemble ─► CATHY window ─► ERT forward ─► EnKF | PF | ES-MDA update
        └─ writes <simulations>/<run>/ and <simulations>/<run>_meta/run_meta.json
                                                                           │
nb8_plots_synthetic.py   compares every run with the TRUE sw, psi and parameters
        └─ writes <figures_synth>/<CASE>/ (figures, summary.csv, parameters_summary.csv)
```

**Step 1 - `nb6_synth_truth_ERT.py --case NC|HC`** (about 3 min for NC):
1. Truth = the model of `nb1`/`nb2` with the **true** parameters, on the same time vector as the DA
   (survey times of the case). CATHY does not print at the requested times, so the output is dense
   (`--out-dt`, default 300 s) and linearly interpolated to the exact atmbc times.
2. `sw` -> resistivity with the same Archie chain as the DA (truth Archie parameters:
   `--archie prior` (default, a = 0.3, m = n = 1.7) or `--archie offset`), nodal porosity.
3. pyGIMLi forward with the **real quadrupoles** of each case, no pyGIMLi noise; optional finer truth mesh
   (`--fwd-mesh`) to avoid an inverse crime. Rows with rhoa outside `[rho_min/r, rho_max*r]`
   (`--rhoa-ratio`, default 3) are dropped.
4. Noise added by the script: `rhoa_obs = rhoa_true * exp(0.05 * N(0,1))` (`--noise-pct`, seeded).
5. Files with the format and names of the real data (`st5t<ii>.tldat`), the virtual TDR series (15, 30,
   50 cm, not assimilated) and `truth_<CASE>.pkl` (sw, psi, parameters, times, noise settings).

**Step 2 - `nb7_DA_ERT_synthetic.py`** (one run = one case/scenario/DA type/N/seed):
- reads the synthetic `.tldat`; checks that DA and truth share the atmbc times and cell size;
- builds the prior ensemble from the scenario (`--prior-bias biased` applies Ks x5, SATCELL_VG x1.3,
  n_VG -0.08, thetar_VG +0.01; `none` centres the prior on the truth); the prior is cached so that the DA
  and the open loop of one seed share the same ensemble;
- `--da-type enkf | pf` (sequential, one CATHY window per survey, observation stamped at the window start,
  compared with the state at its end) or `esmda`; `--no-da` is the open loop;
- `ERT_DATA_ERR = 5 %` on `log10(rhoa)` (`--obs-err`, `--log-obs`), forward noise always 0;
- run name: `SYN_<CASE>_<scenario>_<enkf|pf|esmda>_N<nens>_s<seed>[_noDA][_unbiased]`;
- states are stored compactly (`compact_states`, `--store-full all|analysis`, `--cleanup-vtk`).

**Step 3 - `nb8_plots_synthetic.py --case NC`** (see section 7 for the outputs).

More detail, decisions and known limits (inverse crime, single Archie set, Ks anisotropy, memory of HC at
0.1 m): `notebook/synthetic_twin_plan.md`.

---

## 5. Running

```bash
cd notebook
./run_twin.sh help                       # all commands and variables (full documentation: header of the file)
python project_paths.py                  # check the resolved folders first
```

**Synthetic twin**
```bash
CASE=NC ./run_twin.sh truth                          # nb6: truth + synthetic ERT
CASE=NC ./run_twin.sh da ic_VGP_upd                  # nb7: one EnKF run (N=64, unbiased prior by default)
CASE=NC DA=pf BIAS=biased ./run_twin.sh da ic_Ks_upd # particle filter, biased prior
CASE=NC ./run_twin.sh pair ic_Ks_upd                 # DA + its open loop
CASE=NC ./run_twin.sh status                         # one line per run (ok / partial / failed ...)
CASE=NC ./run_twin.sh plots                          # nb8
./run_twin.sh everything --dry-run                   # list the whole matrix, simulate nothing
CASES="NC" SCENARIOS="ic ic_Ks_upd" DAS="enkf" SEEDS="1" NENS=32 ./run_twin.sh everything
python run_matrix.py --preset quick --dry-run        # presets: quick | main | ens
```
`run_matrix.py` skips runs whose `<run>_meta/run_meta.json` says `finished` (resume = same command; `--force` reruns).

**Real data**
```bash
./run_twin.sh real-everything nc_ERT_ic nc_ERT_ic_VGP_upd   # PF, EnKF, ES-MDA per scenario, then figures
REAL_DAS="enkf" ./run_twin.sh real-seq nc_ERT_ic_VGP_pert   # one scenario, EnKF only
NENS_REAL=16 NITER_REAL=2 ./run_twin.sh real-esmda nc_ERT_ic   # small ES-MDA test
```
Tests: `cd notebook && python -m pytest tests -q` (some need pyCATHY, see `PYCATHY_SRC` in the tests).

---

## 6. Scenarios and the scenario `.toml` files

A **scenario** states which quantities are perturbed (prior ensemble) and which are updated by the DA.
Two files, same schema:

| File | Generated by | Used by | Names |
|---|---|---|---|
| `notebook/scenario_nc_ERT.toml` | `make_scenario_nc_VGP.py` (one soil layer = one entry, 15 layers) | nb5, nb6 (`REAL_SCENARIO=<name>`) | `nc_ERT_ic`, `nc_ERT_ic_VGP_pert`, `nc_ERT_ic_VGP_upd`, `nc_ERT_ic_VGP_Archie_pert`, `nc_ERT_ic_VGP_pert_Archie_upd`, `nc_ERT_ic_VGP_Archie_upd` |
| synthetic (built in memory by `make_scenarios.build()`; written to `scenarios_twin.toml` by `python make_scenarios.py`) | `make_scenarios.py` | nb7 (`--scenario <tag>` or `--scenario-file`), nb8 | `<CASE>_ic`, `_ic_VGP`, `_ic_VGP_upd`, `_ic_Ks_upd`, `_ic_Archie_upd`, `_ic_VGP_Ks_upd`, `_full_upd` |

The files are **generated: edit the generator (`MATRIX` in `make_scenarios.py`, or `make_scenario_nc_VGP.py`),
not the toml.** Every table is `[scenario.<name>]`; all lists have one entry per name in `per_name`; an
entry is itself a list when the control has several groups.

| Key | Meaning |
|---|---|
| `per_name` | perturbed quantities: `ic` (initial pressure head, m), `SATCELL_VG`, `n_VG`, `thetar_VG` (van Genuchten), `Ks`, `a_Archie`, `m_Archie`, `n_Archie` |
| `pert_control_name` | how a quantity maps onto the soil: `None` (one scalar), `layers` (one value per soil layer), `material` (NC: groups A = layers 0-5, B = 6-8, C = 9-14), `zone_material` (HC: `[zone, layers]` groups, zone 2 = compacted band) |
| `material_layers` / `zone_material_layers` | the groups of the control above |
| `ic_layers` | groups of the initial condition (`ic` uses `pert_control_name = "material"`, one value per depth range of the truth, prior mean = `ic_true`): `[0-5] [6-9] [10-14]`, not the soil materials (layer 9 differs); needs the pyCATHY `ic` + `material` support |
| `per_type` | sampling type; `"None"` = absolute (prior = N(`per_mean`, `per_sigma`) truncated to `per_bounds`, physical units) |
| `per_nom`, `per_mean` | nominal / prior mean (equal in the generated files) |
| `per_sigma` | prior standard deviation (physical units; for a `log` transform of `Ks` see `KS_PRIOR` in `make_scenarios.py`) |
| `per_bounds` | `{ min, max }` physical bounds, also used to clip the analysed values |
| `sampling_type` | `"normal"` (or `"lognormal"` with `KS_PRIOR=lognormal`) |
| `transf_type` | transform in the update space: `"None"`, `"log"`, `"log-ratio"` (Archie a / m, n are set by nb5-nb7) |
| `listUpdateParm` | what the filter **updates**: `"St. var."` (states) plus the parameter keys: `<name><k>` for `layers` / `material` groups (e.g. `Ks0, Ks1, Ks2`, `SATCELL_VG0 ... SATCELL_VG14`) or the bare name for scalars (`a_Archie`) |

Reading the table: a quantity that is in `per_name` but not in `listUpdateParm` is **perturbed only** (prior
spread, open ensemble); one that is in both is **estimated**. `["St. var."]` alone = states only.

Example (generated `NC_ic_Ks_upd`: initial condition and Ks perturbed, states and Ks updated):
```toml
[scenario.NC_ic_Ks_upd]
per_name          = [ "ic", "Ks" ]
pert_control_name = [ "material", "material" ]
per_type          = [ [ "None", "None", "None" ], [ "None", "None", "None" ] ]
per_mean          = [ [ -2.4, -1.3, -1.15 ], [ 1.05e-05, 1.15e-05, 7.5e-06 ] ]
per_sigma         = [ [ 0.5, 0.5, 0.5 ],  [ 1.05e-05, 1.15e-05, 7.5e-06 ] ]
per_bounds        = [ { min = -5.0, max = -0.05 }, { min = 7.5e-08, max = 0.00115 } ]
sampling_type     = [ "normal", "normal" ]
transf_type       = [ [ "None", "None", "None" ], [ "log", "log", "log" ] ]
material_layers   = [ [0,1,2,3,4,5], [6,7,8], [9,10,11,12,13,14] ]
ic_layers         = [ [0,1,2,3,4,5], [6,7,8,9], [10,11,12,13,14] ]
listUpdateParm    = [ "St. var.", "Ks0", "Ks1", "Ks2" ]
```
(nominals shown are those of the biased prior; `--prior-bias none` / `make_scenarios.py --no-bias` centres them on the truth.)

**Adding a scenario**: add a row to `MATRIX` in `make_scenarios.py` (tag, perturbed names, updated names,
comment), then `python make_scenarios.py` to inspect it; it is available to `nb7 --scenario <tag>`,
`run_matrix.py --scenarios <tag>` and `nb8` without further change. For the real data add it in
`make_scenario_nc_VGP.py` and regenerate `scenario_nc_ERT.toml`. Use `--scenario-file my.toml` (nb7) or
`[scenarios].real` in `data_path.toml` (nb5/nb6) to point to another file.
Known limit: the van Genuchten parameters are drawn independently (no mutual correlation).

---

## 7. Outputs and plots

| Output | Location | Written by |
|---|---|---|
| Truth, synthetic ERT, virtual TDR | `<synthetic_dir>/<CASE>/` | nb6 (synthetic) |
| Run folder (CATHY, ensemble, states) | `<simulations>/<run>/` (`<run>_states/cycle_*.npz`, `<run>_states.nc`, `<run>_df_DA.pkl` once at the end, `<run>_dict_parm_pert.pkl`, `<run>_ERT_checks.json`, `<run>_rejections.json`) | nb7 / nb5 |
| Run metadata | `<simulations>/<run>_meta/run_meta.json` (all options, checks, status, traceback) | nb7 |
| Run manifest | `<simulations>/<run>_meta/manifest.json` (+ `uncommitted.patch`, `pycathy_uncommitted.patch`, `inputs/`) : code versions, branches, inputs, environment (section 8) | nb5-nb7 (run_registry) |
| Filtered ERT files | `<simulations>/ert_filtered_<run>/` | nb7 |
| Prior ensemble cache | `<cache>/prior_*.pkl` | nb5-nb7 |
| Real-data run folders | `<simulations>/NC_ERT_[PF_|ESMDA_]{withDA|withoutDA}_<scenario>[_mat]/` with `_finished.txt` (ES-MDA: `_esmda_result.pkl`) | nb5 / nb6 |
| Batch logs | `<logs>/<run>.log`, `matrix_summary.json` | run_matrix.py |
| Synthetic figures and tables | `<figures_synth>/<CASE>/` | nb8 |
| Real-data figures | `<figures_real>/<method>_<scenario>[_mat]/` | nb5 plots |
| Run registry | `<run_registry>` | nb5-nb7, run_matrix.py (section 8) |

nb8 per run: `<run>_states.png` (truth vs ensemble at the 15/30/50 cm nodes, DA vs open loop, virtual TDR),
`<run>_profiles.png` (profiles at survey times), `<run>_metrics.png` (RMSE, spread, 90 % coverage, CRPS in
time; all nodes, inside / outside the ERT block), `<run>_params.png` (parameters vs truth, prior -> final).
Across runs: `summary.csv`, `parameters_summary.csv`, `compare_heatmap_<CASE>.png`,
`compare_boxplot_<CASE>.png`, `compare_params_<CASE>.png`, `compare_nens_<CASE>.png`. Every figure is
independent: a failing one is reported in one line and the others are still made. Note that for ES-MDA the
states are the posterior forward run (conditioned on all data), so their scores are not like for like with
the EnKF/PF forecasts (`states_kind` in `summary.csv`).

---

## 8. Run registry

Every run is recorded in the CSV `[outputs].run_registry` (default `registry/runs_registry.csv`) so you can
see what has been done and what has not. One row per run name, updated in place
(`started` -> `finished` | `failed`; a run that dies without reporting is closed as `failed` when the process
exits):

`run, kind (synthetic|real), case, scenario, da_type, nens, seed, prior_bias, with_da, status, started, finished, duration_min, folder, log, git_commit, git_branch, pycathy_version, manifest`

```bash
python notebook/run_registry.py list [--status failed] [--case NC]   # table of the runs
python notebook/run_registry.py summary                              # counts per status / case / DA type
python notebook/run_registry.py backfill                             # register past runs found in <simulations>
```

### Run manifest (reproducibility)

When a run starts, `run_registry.track()` also writes `<simulations>/<run>_meta/manifest.json` (completed with
status / duration when the run ends; nb5, nb6 and nb7 need no change; `track(..., manifest=False)` disables it).
It answers "exactly what produced this run?":

| Section | Content |
|---|---|
| `cathy_da_ert` | this repo: full commit, **branch**, upstream, remote url (credentials stripped), `git describe`, dirty flag + list of modified files. A dirty tree is saved as `uncommitted.patch` (and its sha256) so the run can be rebuilt |
| `pycathy` | version, folder it is imported from (without importing it), and if it is a git clone (editable install or `[code].pycathy`): commit, **branch**, remote, dirty flag (`pycathy_uncommitted.patch`), e.g. your local twin fixes |
| `cathy` | the Fortran model: git commit/branch of its source folder (`[code].cathy`, `$CATHY_SRC`, or next to pyCATHY) and sha256 / size / date of the executable found in the run folder (`$CATHY_EXE` to force it) |
| `packages`, `key_packages` | every installed package version (pip freeze like), and numpy / scipy / pandas / xarray / pygimli ... on top |
| `python`, `platform`, `environment` | interpreter, OS, host, CPU count, conda / venv, thread variables (`OMP_NUM_THREADS`, ...) |
| `inputs` | sha256 of `data_path.toml`, the scenario `.toml` files, `site_params.json`, the ERT schedule, the forward mesh; name+size digest of `data/3d`, `sensors`, `mesh` folders (no full read); copies of the small config files in `<run>_meta/inputs/` |
| `command`, `fields` | command line + working directory; seed, nens, da_type, scenario, ... as given to the registry |

```bash
python notebook/run_manifest.py show    <run>            # one-screen summary of a run
python notebook/run_manifest.py compare <run_a> <run_b>  # what differs: commits, branches, versions, inputs, seed
python notebook/run_manifest.py collect                  # the manifest of the current environment (before a run)
```

To rebuild a run: check out `cathy_da_ert.commit`, `git apply uncommitted.patch` if it was dirty, do the same in
the pyCATHY clone, then run again with the recorded `command` and `fields.seed`. A manifest problem never stops a
run (one line on stderr). `git_branch` and `pycathy_version` are also in the registry csv.

---

## Lighter DA folder

Where the disk goes, and what is already available (no change of the DA result expected):

| What | Size order | How to reduce |
|---|---|---|
| member CATHY output, mostly `.vtk`, in `<run>/DA_Ensemble/` | ~31 of ~38 MB per member | `nb7 --cleanup-vtk all` (default; acts at the END of the run, so the peak is still NENS x ~38 MB) |
| legacy long table `<run>_df_DA.pkl` (float64, node ids repeated) | ~314 MB | compact store (default, `compact_states=True`): `<run>_states/cycle_*.npz`, float32, ~35-65 MB |
| compact store | ~35-65 MB | `--store-full analysis` + `--store-stride N` (probe nodes only) |

```bash
python notebook/ligher_run.py report [RUN_DIR ...]      # bytes per category (vtk / pkl / npz / nc / ...), biggest files
python notebook/ligher_run.py prune RUN_DIR [--apply]   # list (or delete) the legacy long table, only if the compact store exists
python notebook/ligher_run.py convert FILE_df_DA.pkl [--delete]   # legacy pickle -> parquet (float32, zstd)
```

`prune` never touches the vtk and never touches runs without a compact store (nb8 then needs the long table).
The real-data runs of nb5/nb6 have no compact store.

**Measuring and accepting a lighter setting** (protocol of `notebook/todo.md`, section 7):

```bash
# 7.0 baseline: wall time, peak RSS and disk per category of one run, as json
python notebook/lighter_protocol.py measure baseline_s1.json --run-dir <RUN_DIR> -- ./notebook/run_twin.sh ...

# 7.4 a lighter setting is accepted if, for each score of nb8 summary.csv, its mean over the seeds is
#     within the seed-to-seed spread of the baseline (>= 2 baseline seeds, 3 recommended); exit code 0 = accepted
python notebook/lighter_protocol.py accept base_s1.csv base_s2.csv base_s3.csv --light light_s1.csv light_s2.csv light_s3.csv
```

`measure` records the peak RSS of the child processes, so run it on the command that starts the whole DA.
`run_manifest.py compare <run_a> <run_b>` checks that only the intended option changed between the two runs.
Nothing was measured on a real DA run yet: the sizes above are those of `notebook/synthetic_twin_plan.md` section 16.

**Sharing the results with a colleague**

```bash
python notebook/export_results.py ~/share/twin_results --zip                  # all runs, compare_*.png figures
python notebook/export_results.py ~/share/nc_pf --runs 'SYN_NC_*_pf_*' --only-ok --figures all
```

It copies what already exists (nb8 `summary.csv` / `parameters_summary.csv`, figures, the registry, per-run
`manifest.json`, `run_meta.json`, issues log and, for failed runs, the end of the run log) and writes one
`overview.csv` (a row per run: status, duration, disk, key scores) and a `README.txt` that explains the columns.
Nothing is recomputed; states, vtk and pickles are never copied.

---

## Status and limits
- This branch was checked with unit tests (pytest, stubs for pyCATHY/pyGIMLi) and `py_compile`/`bash -n`;
  **no CATHY/pyGIMLi run was done with these changes**. Run `python notebook/project_paths.py` and a small
  `./run_twin.sh everything --dry-run` first.
- The NC synthetic chain was run end to end once (see `notebook/synthetic_twin_plan.md`, section 15);
  HC and most parameter scenarios were not.
- `data/` and `notebook/logs/` contain large tracked files (3 x `duckdb.exe`, run logs); see the storage
  notes in `notebook/synthetic_twin_plan.md` section 16 for the lighter DA output options.
