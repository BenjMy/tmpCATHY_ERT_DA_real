# Running the Synthetic DA Study

Run the study orchestration script from the `notebook` directory.

## Basic usage

```bash
cd notebook
git pull    # Get run_study.py (commit e749412)

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

## More examples

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

## Available options

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

## Output locations

| Output | Location |
|---|---|
| Quick runs | `<simulations>_n10/` (run names end in `_n10`; `_n<N>` with `--n-times N`) |
| Full runs | `<simulations>/` (no tag) |
| Study report | `<simulations>_n10/study_report.json` |
| Execution logs | `notebook/logs/<run>.log`, `notebook/logs/diagnostics_<case>.log` |
| Figures and summary | `<simulations>_n10/figures_nb8/<case>/` (`NC` or `HC`) |
| Summary table | `<simulations>_n10/figures_nb8/<case>/summary.csv` |
| Parameter table | `<simulations>_n10/figures_nb8/<case>/parameters_summary.csv` |

## Resuming and error handling

- **Automatic resume:** Re-running the same command resumes the study. Runs that have already finished successfully are skipped.
- **`skipped` in the report** means the run was already `finished` (its results exist), not that it was not run.
- **Force execution:** Use `--force` to rerun everything, including previously completed runs.
- **Exit code:** The script returns exit code `1` if any simulation or diagnostic step fails.
- **Diagnostics `rc=1` does not always mean no figures:** `plots_synthetic.py` exits with `1` when any run was skipped or not finished (for example a failed run), or when a figure or table was skipped, even though it wrote the others. Read the `status per run` block at the end of `notebook/logs/diagnostics_<case>.log`.
- **Failed runs:** the reason is in `study_report.json` and in `<run>_rejections.json` (rejected ensemble members) next to the run. Fix the cause, then rerun the same command: only the failed runs are redone.
