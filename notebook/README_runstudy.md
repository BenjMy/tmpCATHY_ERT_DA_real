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

## Available options

| Option | Description |
|---|---|
| `--scenarios` | Select scenarios: `ic ic_Ks_upd ic_VGP_upd ic_VGP_Ks_upd full_upd`, or `all`. |
| `--da enkf` | Run only the EnKF filter. |
| `--da pf` | Run only the PF filter. |
| `--no-open-loop` | Skip the without-DA runs. |
| `--seeds 42 43 44` | Add multiple random seeds. |
| `--cases NC HC` | Include the HC case. |
| `--stop-on-fail` | Stop execution after the first failure. |
| `--extra "--obs-err 10 --inflation 1.2"` | Pass additional arguments to `da_ert_synthetic.py`. |
| `--diag-extra "--detail-seed 42 --t-from-h 0"` | Pass additional arguments to `plots_synthetic.py`. |

## Output locations

| Output | Location |
|---|---|
| Quick runs | `<simulations>_n10/` (run names end in `_n10`) |
| Study report | `<simulations>_n10/study_report.json` |
| Execution logs | `notebook/logs/<run>.log` |
| Figures and summary | `<simulations>_n10/figures_nb8/NC/` |
| Summary table | `<simulations>_n10/figures_nb8/NC/summary.csv` |

## Resuming and error handling

- **Automatic resume:** Re-running the same command resumes the study. Runs that have already finished successfully are skipped.
- **Force execution:** Use `--force` to rerun everything, including previously completed runs.
- **Exit code:** The script returns exit code `1` if any simulation or diagnostic step fails.
