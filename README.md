# tmpCATHY_ERT_DA_real

Hydrological simulation of the Agripolis infiltration test (UNIPD, Oct 2023) with
[pyCATHY](https://github.com/BenjMy/pycathy_wrapper) and **data assimilation (DA) of 3D ERT data**, in two flavours:

| Flavour | Observations | Truth known? | Scripts |
|---|---|---|---|
| **Real data** | measured ERT (`st5t<ii>.tldat`, `st1t<ii>.tldat`) | no | `nb5_DA_ERT_soil3D_nc.py` (EnKF / PF), `nb6_ESMDA_ERT_soil3D_nc.py` (ES-MDA), `nb5_DA_ERT_soil3D_nc_plots.py` |
| **Synthetic twin** | ERT simulated from a known truth (5 % noise) | yes | `synthetic_truth_ert.py` -> `da_ert_synthetic.py` -> `plots_synthetic.py` |

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
  forward_model_nc.py, nb2_*, nb3_*  forward models (NC, HC) and measurement-vs-model comparison
  nb4*, nb4a*                  twin experiment on soil moisture (TDR) only
  nb5_*, nb6_ESMDA_*           real-data ERT DA
  synthetic_truth_ert.py, da_ert_synthetic.py, plots_synthetic.py (+ plots_synthetic_helpers.py)
                               synthetic ERT twin: truth -> DA -> scores/figures
  model_config.py, nc_model.py case configuration (soil tables, grid, IC, ERT paths) and helpers
  make_scenarios.py            generates the synthetic scenarios (NC + HC)
  make_scenario_nc_VGP.py      generates scenario_nc_ERT.toml (real-data scenarios)
  scenario_nc_ERT.toml         real-data scenarios (generated, do not edit by hand)
  run_matrix.py, run_twin.sh   batch drivers
  run_study.py                 quick scenario study (runs + diagnostics in one command)
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
- `forward_model_nc.py` (NC) and `nb2_cathy_soil3D_hc.py` (HC): CATHY infiltration models, reference runs.
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
synthetic_truth_ert.py   truth run ─► sw, psi ─► Archie ─► pyGIMLi forward ─► 5 % noise ─► .tldat
        └─ writes data/3d/synthetic/<CASE>/{truth_<CASE>.pkl, procdata/*.tldat, tdr_synth_<CASE>.csv}
                                                                           │
da_ert_synthetic.py  prior ensemble ─► CATHY window ─► ERT forward ─► EnKF | PF | ES-MDA update
        └─ writes <simulations>/<run>/ and <simulations>/<run>_meta/run_meta.json
                                                                           │
plots_synthetic.py   compares every run with the TRUE sw, psi and parameters
        └─ writes <figures_synth>/<CASE>/ (figures, summary.csv, parameters_summary.csv)
```

**Step 1 - `synthetic_truth_ert.py --case NC|HC`** (about 3 min for NC):
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

**Step 2 - `da_ert_synthetic.py`** (one run = one case/scenario/DA type/N/seed):
- reads the synthetic `.tldat`; checks that DA and truth share the atmbc times and cell size;
- builds the prior ensemble from the scenario (`--prior-bias biased` applies Ks x5, SATCELL_VG x1.3,
  n_VG -0.08, thetar_VG +0.01; `none` centres the prior on the truth); the prior is cached so that the DA
  and the open loop of one seed share the same ensemble;
- `--da-type enkf | pf` (sequential, one CATHY window per survey, observation stamped at the window start,
  compared with the state at its end) or `esmda`; `--no-da` is the open loop;
- `ERT_DATA_ERR = 5 %` on `log10(rhoa)` (`--obs-err`, `--log-obs`), forward noise always 0;
- run name: `SYN_<CASE>_<scenario>_<enkf|pf|esmda>_N<nens>_s<seed>[_noDA][_unbiased]`;
- states are stored compactly (`compact_states`, `--store-full all|analysis`, `--cleanup-vtk`).

**Step 3 - `plots_synthetic.py --case NC`** (see section 7 for the outputs).

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

**Quick scenario study + diagnostics in one command** (`run_study.py`, run from `notebook/`)

```bash
# 1. Preview: lists the 15 runs and the diagnostics call, without running simulations
python run_study.py --dry-run

# 2. Quick test: first 10 survey times, N=32, EnKF + PF, with DA + open loop, then diagnostics
python run_study.py

# 3. Same as above, but only 2 scenarios and a smaller ensemble (faster)
python run_study.py --scenarios ic ic_Ks_upd --nens 16

# 4. Use a different number of survey times
python run_study.py --n-times 5

# 5. Full run on all 23 surveys (N=100, same run names and folder as run_matrix.py)
python run_study.py --full

# 6. Run diagnostics only, on simulations that have already finished
python run_study.py --stage diag

# 7. Rerun everything, ignoring previously finished runs
python run_study.py --force
```

Quick runs (`--n-times N` of `da_ert_synthetic.py`: the first N surveys, the model stops at the last one) are tagged `n<N>` in the run name and go to `<simulations>_n<N>/`, so they stay apart from full runs and from their rows in `registry/runs_registry.csv`. Failures are listed at the end and in `<sim-dir>/study_report.json`.

<details>
<summary>More <code>run_study.py</code> examples</summary>

Every example below can be previewed first by adding `--dry-run`. Number of runs =
cases × scenarios × `--nens` values × seeds × (number of `--da` types + 1 open loop).
The open loop is run once per case / scenario / N / seed, always with `--da-type enkf`,
whatever `--da` is (`--no-open-loop` removes it).

```bash
# Scenarios (ladder: ic < ic_Ks_upd < ic_VGP_upd < ic_VGP_Ks_upd < full_upd)
python run_study.py --scenarios ic_Ks_upd ic_VGP_Ks_upd        # a chosen subset
python run_study.py --scenarios all                            # the 7 scenarios (adds ic_VGP, ic_Archie_upd): 21 runs
python run_study.py --scenarios NC_ic_Ks_upd                   # a NC_ / HC_ prefix is accepted and ignored

# Filter and open loop
python run_study.py --da pf --scenarios ic_VGP_Ks_upd          # PF + its open loop (2 runs)
python run_study.py --da enkf --no-open-loop                   # EnKF with DA only (5 runs)

# Ensemble size and seeds (the grid is multiplied)
python run_study.py --nens 16 --scenarios ic_Ks_upd ic_VGP_Ks_upd --seeds 1 2 3   # 18 runs: 3 seeds -> box plots
python run_study.py --nens 8 16 32 --scenarios ic              # 9 runs: N comparison (compare_nens_NC.png)

# Prior: unbiased prior mean (run names get the suffix _unbiased)
python run_study.py --prior-bias none --scenarios ic

# Case: HC is much heavier than NC (41.6k vs 1936 nodes): start small
python run_study.py --cases HC --scenarios ic --nens 8 --dry-run
python run_study.py --cases NC HC --scenarios ic ic_Ks_upd --nens 16   # diagnostics run once per case

# Survey times
python run_study.py --n-times 5 --nens 16                      # first 5 surveys (folder <simulations>_n5/)
python run_study.py --full --scenarios ic --nens 100           # all surveys (folder <simulations>/, no _n tag)

# Extra arguments passed to the children (quote them)
python run_study.py --scenarios ic_VGP_Ks_upd --extra "--inflation 1.2"           # to da_ert_synthetic.py
python run_study.py --stage diag --diag-extra "--detail-seed 42 --t-from-h 0 --debug"   # to plots_synthetic.py

# Where things are written (default: data_path.toml)
python run_study.py --sim-dir ~/tmp/study_test --log-dir ~/tmp/study_test/logs
python run_study.py --root /path/to/project/

# Run only part of the pipeline
python run_study.py --stage runs --scenarios ic                # simulations only, no diagnostics
python run_study.py --stage diag --scenarios ic                # diagnostics only (see the note below)

# Failures
python run_study.py --stop-on-fail                             # stop at the first failed run (diagnostics skipped)
python run_study.py --force --scenarios full_upd --da pf       # rerun these runs even if finished
```

Typical study: the 4 scenarios below at N=16 (8 EnKF/PF DA runs + 4 open loops).
Rerun the same command to resume it, or add `--stage diag` to redo only the figures.

```bash
python run_study.py --scenarios ic ic_Ks_upd ic_VGP_Ks_upd full_upd --nens 16
```

**Note on `--stage diag`:** the runs that are scored are the ones the arguments describe, so give
the same `--scenarios`, `--nens`, `--seeds`, `--da`, `--cases`, `--n-times` / `--full` and `--prior-bias`
as for the simulations. Otherwise the planned run names do not match the finished ones.

</details>

<details>
<summary><code>run_study.py</code> options</summary>

| Option | Default | Description |
|---|---|---|
| `--cases NC HC` | `NC` | Case(s) to run. |
| `--scenarios` | the 5-step ladder | `ic ic_Ks_upd ic_VGP_upd ic_VGP_Ks_upd full_upd`, plus `ic_VGP`, `ic_Archie_upd`, or `all` (7). |
| `--da` | `enkf pf` | Filter(s): `enkf`, `pf`, `esmda`. |
| `--no-open-loop` | off | Skip the without-DA runs. |
| `--seeds 42 43 44` | `42` | One run per seed. |
| `--nens 16 32` | 32 (quick) / 100 (`--full`) | Ensemble size(s), each >= 2. |
| `--n-times N` | `10` | Quick test: first N survey times (`0` = all). |
| `--full` | off | All survey times (same as `--n-times 0`). |
| `--prior-bias biased\|none` | `biased` | Biased or unbiased prior mean. |
| `--stage runs\|diag\|both` | `both` | Simulations, diagnostics, or both. |
| `--dry-run` | off | Print the runs and the diagnostics call, simulate nothing. |
| `--force` | off | Rerun finished runs. |
| `--stop-on-fail` | off | Stop after the first failed run. |
| `--extra "..."` | empty | Extra arguments for `da_ert_synthetic.py`, e.g. `"--obs-err 10 --inflation 1.2"`. |
| `--diag-extra "..."` | empty | Extra arguments for `plots_synthetic.py`, e.g. `"--detail-seed 42 --t-from-h 0"`. |
| `--root`, `--sim-dir`, `--log-dir`, `--registry` | `data_path.toml` | Override the project root, simulations folder, logs folder, `runs_registry.csv`. |
| `--nb7`, `--nb8` | `da_ert_synthetic.py`, `plots_synthetic.py` | Alternative simulation / diagnostics scripts. |

</details>

Resuming and error handling:

- **Automatic resume:** Re-running the same command resumes the study. Runs that have already finished successfully are skipped.
- **`skipped` in the report** means the run was already `finished` (its results exist), not that it was not run.
- **Force execution:** Use `--force` to rerun everything, including previously completed runs.
- **Exit code:** The script returns exit code `1` if any simulation or diagnostic step fails.
- **Diagnostics `rc=1` does not always mean no figures:** `plots_synthetic.py` exits with `1` when any run was skipped or not finished (for example a failed run), or when a figure or table was skipped, even though it wrote the others. Read the `status per run` block at the end of `notebook/logs/diagnostics_<case>.log`.
- **Failed runs:** the reason is in `study_report.json` and in `<run>_rejections.json` (rejected ensemble members) next to the run. Fix the cause, then rerun the same command: only the failed runs are redone.

**If the ensembles collapse** (state coverage near 0, parameter `mean_contraction` near 1 and `frac_truth_in_90` near 0 in `parameters_summary.csv`), test one change at a time against the baseline study, N=16 here to stay cheap:
```bash
# 1. Larger ensemble (cost grows about linearly with N)
python run_study.py --scenarios ic_Ks_upd ic_VGP_Ks_upd --nens 32 64

# 2. Stronger state inflation (default 1.1; whether the parameters are inflated too is not checked here)
python run_study.py --scenarios ic_Ks_upd ic_VGP_Ks_upd --nens 16 --da enkf --extra "--inflation 1.3" --sim-dir ../simulations_infl13

# 3. Log-normal Ks prior instead of the default 'normal100' (Gaussian, sd = 100 % of the nominal: ~16 % negative draws
#    piled at the lower bound). Positive by construction; KS_SIGMA_LOG is in natural-log units (0.5 = x1.65 at 1 sigma). Untested option.
KS_PRIOR=lognormal KS_SIGMA_LOG=0.5 python run_study.py --scenarios ic_Ks_upd ic_VGP_Ks_upd --nens 16 --sim-dir ../simulations_lognormal

# 4. PF: jitter on the resampled states and resampling only when N_eff < N/2 (defaults 0 and 1.0 = as before).
#    With a threshold < 1 the weights must be carried between cycles: not verified, see todo.md 5.2.
python run_study.py --da pf --scenarios ic_VGP_Ks_upd --nens 16 --extra "--pf-jitter-std-state 0.01 --pf-resample-threshold 0.5" --sim-dir ../simulations_pfjitter

# 5. Separate the effect of the biased prior (Ks x5, SATCELL x1.3, ...) from the filter
python run_study.py --prior-bias none --scenarios ic_Ks_upd ic_VGP_Ks_upd --nens 16 --sim-dir ../simulations_unbiased
```
Run names do not encode `--extra` or `KS_PRIOR`, so give every variant its own `--sim-dir` (otherwise its runs are skipped as already finished); the environment variable also has to be set for `--stage diag`. Each variant reruns its open loop. Compare the variants side by side:
```python
import pandas as pd
dirs = ['simulations_n10', 'simulations_lognormal']          # sibling folders of notebook/
cols = ['sw_all_rmse', 'sw_all_coverage', 'sw_all_spread_skill', 'sw_all_crps']
print(pd.concat({d: pd.read_csv(f'../{d}/figures_nb8/NC/summary.csv').set_index('run')[cols] for d in dirs}, axis=1).round(3))
```
Good signs: coverage towards 0.9, spread_skill towards 1, `mean_contraction` below 1 with `frac_truth_in_90` above 0.

<details>
<summary>Cookbook: diagnostics, plots, redoing runs, HC, which code ran</summary>

All commands are run from `notebook/`. Add `--dry-run` to any `run_study.py` call to see what it would do.

**Where are my runs?** Quick studies write to `<simulations>_n<N>/`, every variant to the `--sim-dir` you gave.
The registry (`registry/runs_registry.csv`) has **one row per run name and run folder**: the same name run
in two folders (baseline and `--sim-dir ../simulations_lognormal`, say) gives two rows, and the `sim_dir` column of
`run_registry.py list` tells them apart (`--sim-dir DIR` keeps one folder). Rows written before this change have no
folder and are completed by the next run of that name. The per-folder truth is still
`<sim-dir>/<run>_meta/run_meta.json` (this is what `run_study.py` reads to skip a run).

```bash
ls -d ../simulations*                                    # candidate folders
python run_registry.py list --status finished --case NC  # what the registry says
python run_registry.py list --status failed              # what to fix
python run_registry.py summary                           # counts per status / case / DA type
ls ../simulations_n10 | grep _meta | head                # runs present in one folder
```

**Diagnostics and plots through `run_study.py --stage diag`** (no simulation is started). Give the same
selectors as for the simulations (`--nens`, `--scenarios`, `--cases`, ...) and the `--sim-dir` of the folder:

```bash
# every NC run of the study in the default quick folder
python run_study.py --stage diag --nens 16 --scenarios ic ic_Ks_upd ic_VGP_Ks_upd full_upd

# one DA type at a time (the open loop is always included: it is shared by the DA types)
python run_study.py --stage diag --da enkf --nens 16 --scenarios ic ic_Ks_upd ic_VGP_Ks_upd full_upd
python run_study.py --stage diag --da pf   --nens 16 --scenarios ic ic_Ks_upd ic_VGP_Ks_upd full_upd

# one scenario, one folder (variant runs: set the same environment variable as for the simulations)
KS_PRIOR=lognormal KS_SIGMA_LOG=0.5 python run_study.py --stage diag --scenarios ic_Ks_upd --nens 16 \
    --sim-dir ../simulations_lognormal

# options of plots_synthetic.py through --diag-extra
python run_study.py --stage diag --scenarios ic_Ks_upd --nens 16 --diag-extra "--detail-seed 42 --t-from-h 24"
python run_study.py --stage diag --scenarios ic_Ks_upd --nens 16 --diag-extra "--no-detail"   # tables + comparison figures only
python run_study.py --stage diag --scenarios ic_Ks_upd --nens 16 --diag-extra "--debug"       # tracebacks instead of one line
```

Scenarios you have not run are reported by the diagnostics, not silently skipped: list only the ones that exist.
Results: `<sim-dir>/figures_nb8/<CASE>/` (`summary.csv`, `parameters_summary.csv`, `compare_*.png`, per-run figures).

**Diagnostics and plots with `plots_synthetic.py` directly** (any selection of runs, partial runs included):

```bash
python plots_synthetic.py --case NC --sim-dir ../simulations_n10               # every finished DA run of the case
python plots_synthetic.py --case NC --sim-dir ../simulations_n10 --include-running   # also runs not finished (partial cycles)
python plots_synthetic.py --case NC --sim-dir ../simulations_n10 \
    --runs SYN_NC_ic_Ks_upd_enkf_N16_s42_n10 SYN_NC_ic_Ks_upd_enkf_N16_s42_noDA_n10   # one DA run + its open loop

# every finished EnKF run of N=16 (the open loops are named ..._enkf_..._noDA, so they are included)
RUNS=$(python run_registry.py list --sim-dir ../simulations_n10 | awk '$2=="finished" && $1 ~ /^SYN_NC_/ && $1 ~ /_enkf_/ && $1 ~ /_N16_/ {print $1}')
python plots_synthetic.py --case NC --sim-dir ../simulations_n10 --runs $RUNS

# every finished PF run of N=16, plus the open loops they are compared with
RUNS=$(python run_registry.py list --sim-dir ../simulations_n10 | awk '$2=="finished" && $1 ~ /^SYN_NC_/ && $1 ~ /_N16_/ && ($1 ~ /_pf_/ || $1 ~ /_noDA/) {print $1}')
python plots_synthetic.py --case NC --sim-dir ../simulations_n10 --runs $RUNS --out ../figures_pf_only
```

The registry lists the runs of all folders: add `--sim-dir ../simulations_n10` to `run_registry.py list` to keep the ones of the folder you plot.
Look at the table of a study:

```bash
column -s, -t < ../simulations_n10/figures_nb8/NC/summary.csv | less -S
tail -n 40 logs/diagnostics_NC.log          # the "status per run" block explains an rc=1
```

**Redoing runs.** A run whose status is `failed` or `started` is rerun by the same command, no `--force` needed;
`--force` is only for runs that are `finished`. The selectors cannot pick one run, but they can isolate it:

```bash
python run_study.py --scenarios ic_Ks_upd --nens 16 --force --dry-run               # preview: what a forced rerun would do
python run_study.py --scenarios ic_Ks_upd --nens 16 --force                         # EnKF + PF + open loop again
python run_study.py --scenarios full_upd --nens 16 --da pf --no-open-loop --force   # the PF run only
python run_study.py --scenarios full_upd --nens 16 --da enkf --no-open-loop --force # the EnKF run only
cp -r ../simulations_lognormal ../simulations_lognormal_backup                      # before a forced rerun you may want to compare with
```

To find out why a run failed:

```bash
python run_registry.py list --status failed
tail -n 30 logs/SYN_NC_full_upd_pf_N16_s42_n10.log
python run_manifest.py show SYN_NC_full_upd_pf_N16_s42_n10      # code versions, inputs, command, traceback
python diag_pf_run.py SYN_NC_full_upd_pf_N16_s42_n10            # inbreeding, dead members, parameter collapse
```

**HC case.** HC has 2 soil zones (compacted band): the soil table has `nzone x nstr` rows (zone first, then
layer) and the 3D property mapping of pyCATHY must give one marker per (zone, layer). With a pyCATHY that only knows
layer markers, `./run_twin.sh truth` stops with `Mismatch between property map (len=30) and unique mesh cell markers
(len=15)`. In `model_config.py` zone 2 defines only material A (layers 1-6); the rows of B and C are NaN unless
`CASES['HC']['zone2_fill_below'] = True` copies zone 1 into them (NaN porosity gives `nodal range [nan, nan]` in the
truth log and `sw and porosity must be > 0 for Archie`). Always do the truth first and check it:

```bash
python project_paths.py                              # folders and [code].pycathy
CASE=HC ./run_twin.sh truth                          # look for: porosity: nodal range [finite, finite]
ls ../data/3d/synthetic/HC/                          # truth_HC.pkl, procdata/, tdr_synth_HC.csv
python run_study.py --cases HC --scenarios ic --nens 8 --dry-run
python run_study.py --cases HC --scenarios ic --nens 8                 # small HC test (41.6k nodes)
python run_study.py --cases HC --scenarios ic ic_Ks_upd --nens 16 --force   # redo HC runs that failed earlier
python run_study.py --stage diag --cases HC --scenarios ic ic_Ks_upd --nens 16
```

HC runs that fail within seconds (registry `duration_min` near 0) did not get to simulate: read
`logs/<run>.log`, and rule out a missing `truth_HC.pkl` first. A truth made with another `zone2_fill_below` than the
DA is not comparable: change the flag, rerun the truth, then the DA.

**Which code ran?** The pyCATHY that is imported is the clone in `[code].pycathy` of `data_path.toml` (the branch
is whatever is checked out there; `data_path.toml` cannot select a branch).

```bash
git -C ~/path/to/pycathy_wrapper branch --show-current
git -C ~/path/to/pycathy_wrapper log --oneline -1
git -C ~/path/to/pycathy_wrapper status --short                       # local uncommitted changes
python -c "import project_paths as P; P.get().setup_sys_path(); import pyCATHY; print(pyCATHY.__file__)"
git branch --show-current && git log --oneline -1                      # this repo
python run_manifest.py show    SYN_NC_ic_Ks_upd_enkf_N16_s42_n10       # branch + commit of both repos for one run
python run_manifest.py compare SYN_NC_ic_enkf_N16_s42_n10 SYN_NC_ic_Ks_upd_enkf_N16_s42_n10   # what differs between two runs
```

Merging a pyCATHy feature branch: `git fetch origin && git log --oneline --left-right origin/<a>...origin/<b>` shows what is
on each side (no `<` lines: `<b>` is ahead, so `git merge --ff-only origin/<b>` on `<a>` is a plain fast-forward).
Check `run_manifest.py show <run>` afterwards: runs made before the merge keep their old commit.

</details>

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
| Study report and logs | `<simulations>_n10/study_report.json` (`_n<N>` with `--n-times N`, none with `--full`), `<logs>/<run>.log`, `<logs>/diagnostics_<case>.log` | run_study.py |
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
see what has been done and what has not. One row per run name and run folder, updated in place
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
python notebook/lighten_run.py report [RUN_DIR ...]      # bytes per category (vtk / pkl / npz / nc / ...), biggest files
python notebook/lighten_run.py prune RUN_DIR [--apply]   # list (or delete) the legacy long table, only if the compact store exists
python notebook/lighten_run.py convert FILE_df_DA.pkl [--delete]   # legacy pickle -> parquet (float32, zstd)
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

**Diagnosing a run (filter inbreeding, dead members, parameter collapse)**

```bash
python notebook/diag_pf_run.py SYN_NC_full_upd_pf_N32_s1_unbiased                         # run name, or a folder with the small files
python notebook/diag_pf_run.py SYN_NC_full_upd_pf_N32_s1_unbiased --compare SYN_NC_full_upd_enkf_N32_s1_unbiased
python notebook/diag_pf_run.py <run> --json diag.json --strict                            # exit code 1 if any check is BAD
```

It reads `<run>_run.json`, `<run>_issues.log`, `<run>_rejections.json`, `<run>_ERT_checks.json` and, when present,
`<run>_states/` (genealogy through `pf_inbreeding.py`); no pyCATHY needed, a missing file only skips its checks.
Each line is `ok` / `WARN` / `BAD` / `info`: interrupted run, degenerate weights, founders left, rejected members that stay NaN
(never refilled), spread lost by each parameter and distance of the final mean from the nominal (judged only for `_unbiased`
runs), negative values of strictly positive parameters before the reflection, ERT data-space spread / misfit, redundancy of the
data (first mode), normalised innovation. `--compare` puts two runs side by side to tell a bias of the prior / identifiability
(present in both) from one of the filter. It ends with the commands to run next (EnKF / ES-MDA at the same seed and N, the PF
settings to change one at a time, N). The thresholds are rules of thumb. Tests: `tests/test_diag_pf_run.py`.

**Renamed files**

| Old name | New name |
|---|---|
| `notebook/ligher_run.py` | `notebook/lighten_run.py` |
| `notebook/nb1_cathy_soil3D_nc.py` | `notebook/forward_model_nc.py` |
| `notebook/nb6_synth_truth_ERT.py` | `notebook/synthetic_truth_ert.py` |
| `notebook/nb7_DA_ERT_synthetic.py` | `notebook/da_ert_synthetic.py` |
| `notebook/nb8_plots_synthetic.py` | `notebook/plots_synthetic.py` |
| `notebook/nb8_helpers.py` | `notebook/plots_synthetic_helpers.py` |

Done with `git mv`, so `git log --follow <new name>` shows the full history. Runs made before the rename (their
`manifest.json`, `run_meta.json`, logs, `truth_meta_<CASE>.json`, the registry) record the command with the OLD names:
that is a historical record and was left as it is. The short names `nb6` / `nb7` / `nb8` in comments, and the output
folder `figures_nb8/`, are unchanged.

---

## Status and limits
- This branch was checked with unit tests (pytest, stubs for pyCATHY/pyGIMLi) and `py_compile`/`bash -n`;
  **no CATHY/pyGIMLi run was done with these changes**. Run `python notebook/project_paths.py` and a small
  `./run_twin.sh everything --dry-run` first.
- The NC synthetic chain was run end to end once (see `notebook/synthetic_twin_plan.md`, section 15);
  HC and most parameter scenarios were not.
- `data/` and `notebook/logs/` contain large tracked files (3 x `duckdb.exe`, run logs); see the storage
  notes in `notebook/synthetic_twin_plan.md` section 16 for the lighter DA output options.
