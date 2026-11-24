#!/usr/bin/env bash
# =============================================================================================
# run_twin.sh - ONE entry point of the two DA experiments on ERT data
#
#   SYNTHETIC TWIN  nb6 (synthetic truth + noisy ERT) -> nb7 (DA: EnKF | PF | ES-MDA) -> nb8 (figures / tables)
#   REAL DATA       nb5 (EnKF | PF, sequential)  +  nb6_ESMDA (ES-MDA smoother)  -> nb5_..._plots (figures)
#
# PUT THIS FILE NEXT TO THE .py FILES (nb5_..., nb6_..., nb7_..., nb8_..., run_matrix.py, make_scenarios.py,
# nc_model.py, ert_schedule.py, nb5_plot_helpers.py, scenario_nc_ERT.toml, ...).
# It runs from that folder (the scripts expect ROOT=../ = the project folder with data/ and simulations/).
#
# ---- QUICK START: THE TWO CALLS THAT RUN EVERYTHING ---------------------------------------
#   SYNTHETIC   ./run_twin.sh everything         # NC + HC; 6 scenarios x (EnKF, PF, ES-MDA) + open loop; then figures
#               ./run_twin.sh everything --dry-run          # print the list of runs first (nothing is simulated)
#   REAL DATA   ./run_twin.sh real-everything    # real ERT, scenario by scenario: PF, EnKF, then ES-MDA (N=6); then figures
#   both        WITH_REAL=1 ./run_twin.sh everything        # synthetic matrix first, then the real-data runs
#   Both are RESUMABLE: runs that already finished are skipped when you launch the command again.
#   Long jobs: nohup ./run_twin.sh real-everything > real.log 2>&1 &     then   tail -f real.log
# ---- SYNTHETIC TWIN: commands and examples ----------------------------------------------
#   (variables before the command; defaults: CASE=NC NENS=64 SEED=1 DA=enkf BIAS=none)
#   restrict 'everything':  CASES="NC" SCENARIOS="ic ic_Ks_upd" DAS="enkf esmda" SEEDS="1" NENS=32 ./run_twin.sh everything
#   chmod +x run_twin.sh
#   ./run_twin.sh syn-pf                                 # 'everything' with the particle filter only (PF + open loops; CASES, SCENARIOS, NENS apply)
#   ./run_twin.sh plan                                   # what 'everything' would run, and the meaning of each scenario
#   CASE=NC ./run_twin.sh truth                          # synthetic truth + ERT (5 % noise) of the uncompacted case
#   CASE=HC ./run_twin.sh truth                          # same for the compacted case
#   CASE=NC ./run_twin.sh da ic_Ks_upd                   # one EnKF run (N=64, prior centred on the truth)
#   CASE=NC DA=pf ./run_twin.sh da ic_VGP_upd            # particle filter
#   CASE=NC DA=esmda ./run_twin.sh da ic_Ks_upd          # smoother ES-MDA (parameters only, 4 iterations)
#   CASE=NC BIAS=biased ./run_twin.sh da full_upd        # biased prior
#   EXTRA="--obs-err 10 --inflation 1.2" ./run_twin.sh da ic_VGP_upd
#   CASE=HC NENS=32 ./run_twin.sh pair ic_Ks_upd         # DA + open loop of one scenario
#   CASE=NC ./run_twin.sh compare-da                     # EnKF vs PF vs ES-MDA (same scenario/seed/N)
#   CASE=NC ./run_twin.sh compare-scenarios              # all 6 scenarios, EnKF, DA + open loop
#   CASE=NC ./run_twin.sh compare-bias | compare-nens    # biased vs unbiased prior | ensemble size (NENS_LIST="16 32 64")
#   CASE=NC ./run_twin.sh status                         # one line per run: ok / partial / no states / no results / failed
#   CASE=NC ./run_twin.sh plots                          # status table, then nb8 figures + tables (REAL=0: without the real-data ones)
# ---- REAL DATA: commands and examples ---------------------------------------------------
#   Scripts (defaults = same folder as this file; override with the variable):
#     REAL_NB5    nb5_DA_ERT_soil3D_nc.py       sequential DA, EnKF or PF (REAL_USE_PF=1)
#     REAL_NB6    nb6_ESMDA_ERT_soil3D_nc.py    ES-MDA smoother (parameters only, 'ic' estimated as a parameter)
#     REAL_PLOTS  nb5_DA_ERT_soil3D_nc_plots.py figures (--scenario, --da enkf|pf|esmda); no truth: no score figures
#   Commands ([sc...] = scenario names; none given = every name of REAL_SCENARIOS):
#     ./run_twin.sh real-everything [sc...]    for each scenario: PF (noDA, withDA) -> EnKF (noDA, withDA) -> ES-MDA
#                                              (withDA, withoutDA; ES-MDA uses a small ensemble, NENS_ESMDA_REAL,
#                                              default 6); then real-plots (PF / EnKF are skipped if absent from
#                                              REAL_DAS; the noDA runs are skipped with REAL_OPEN_LOOP=0)
#                                              (per scenario, in the order of REAL_DAS; EnKF runs withoutDA first)
#     ./run_twin.sh real-pf [sc...]            PF only: per scenario withoutDA + withDA (own PF open loop), parallel
#                                              (REAL_JOBS), then PF figures; no EnKF, no ES-MDA
#     ./run_twin.sh real-seq [sc...]           nb5: every method of REAL_DAS ("enkf pf") x scenario, + open loop
#     ./run_twin.sh real-esmda [sc...]         nb6: ES-MDA x scenario, + open loop
#     ./run_twin.sh real-plots [sc...]         figures of every FINISHED run (method x scenario); missing runs are skipped
#     ./run_twin.sh real                       ONE run of nb5 with the settings written in the script (no loop)
#   Examples:
#     ./run_twin.sh real-everything nc_ERT_ic nc_ERT_ic_VGP_upd      # two scenarios only
#     REAL_DAS="pf" ./run_twin.sh real-seq                           # PF only (no open loop: PF compares with the EnKF one)
#     REAL_DAS="enkf" ./run_twin.sh real-seq nc_ERT_ic_VGP_pert      # EnKF, one scenario
#     REAL_USE_PF=1 ./run_twin.sh real                               # a single PF run with the nb5 settings
#     NENS_REAL=16 NITER_REAL=2 ./run_twin.sh real-esmda nc_ERT_ic   # small test (16 members, 2 ES-MDA iterations)
#     REAL_DRY_RUN=1 ./run_twin.sh real-esmda nc_ERT_ic              # nb6 prints the plan, simulates nothing
#     REAL_PLOT_DAS="pf esmda" ./run_twin.sh real-plots              # figures of some methods only
#   Variables of the real case (all optional):
#     REAL_SCENARIOS  scenario names of scenario_nc_ERT.toml (default: the 6 names below)
#     REAL_DAS        methods of real-seq                       default "enkf pf"
#     REAL_PLOT_DAS   methods of real-plots                     default "enkf pf esmda"
#     REAL_OPEN_LOOP  1 = one open loop per scenario            default 1  (0 = skip: halves the cost)
#     REAL_JOBS       scenarios run at the same time (real-everything)   default 2  (1 = sequential)
#     REAL_RESERVE_CORES  cores left free; the rest is split between the jobs, each pinned with taskset   default 2
#     REAL_MIN_FREE_GB    a new job starts only if this much RAM is available                              default 4
#     REAL_STAGGER    seconds between two job launches                                                      default 60
#                     (parallel logs: ../simulations/real_logs/<scenario>.log)
#     REAL_FORCE      1 = redo runs that are already finished   default 0
#     REAL_DRY_RUN    1 = dry run (nb6 only)                    default 0
#     NENS_REAL SEED_REAL NITER_REAL   members | seed | ES-MDA iterations (default: the values written in nb5 / nb6)
#     REAL_RUN_SUFFIX "_mat" if VGP_CONTROL = 'material' in nb5/nb6 (default), "" for 'layers'
#     WITH_REAL=1     'everything' also runs real-seq and real-esmda (the real figures are always attempted at the
#                     end of 'everything', for the real runs that exist)
#   Scenarios of the real case (scenario_nc_ERT.toml; NO Ks scenario here, unlike the synthetic one):
#     nc_ERT_ic                       initial condition only
#     nc_ERT_ic_VGP_pert              + van Genuchten perturbed, NOT updated (prior spread only)
#     nc_ERT_ic_VGP_upd               + van Genuchten updated
#     nc_ERT_ic_VGP_Archie_pert       + Archie perturbed, not updated
#     nc_ERT_ic_VGP_pert_Archie_upd   VGP perturbed only; Archie (a, m, n) updated
#     nc_ERT_ic_VGP_Archie_upd        VGP + Archie updated
#   What each run leaves in ../simulations/ (folder = <prefix>_<withDA|withoutDA>_<scenario>_mat):
#     EnKF  NC_ERT_<..>             PF  NC_ERT_PF_<..>             ES-MDA  NC_ERT_ESMDA_<..>
#     'finished' marker: <folder>_finished.txt (EnKF, PF) | <folder>_esmda_result.pkl (ES-MDA)
#     The open loop (R = 1e99, nothing is assimilated) is shared: EnKF and PF use the EnKF one; ES-MDA has its own.
#     Figures: ../simulations/figures_nb5/<method>_<scenario>_mat/      Logs: <folder>/console_log.txt
#   Notes:
#     - a run restarts from scratch: its folder is deleted and rebuilt (nb5/nb6 behaviour); finished runs are
#       skipped, a crashed or interrupted one is redone.
#     - the real runs are not parallel to each other: they follow one after the other (each one is already
#       parallel inside, parallel=True). Cost = scenarios x methods (+ open loops) x one full run.
#     - ES-MDA is EXPERIMENTAL (see nb6): test a small case first (NENS_REAL=16 NITER_REAL=2).
#     - the figures need the run folders: ES-MDA has no cycle-by-cycle backup (df_DA) unless the run wrote one;
#       then only its misfit / diagnostics figure is made.
#     - matplotlib runs with the Agg backend in the loops (no window blocks the batch): MPLBACKEND=... to change.
#   Other real-data helpers:
#     REAL_PLOTS=/path/script.py ./run_twin.sh real-plots          # another plotting script
#     REAL_NB5=/path/script.py ./run_twin.sh real                  # another nb5 script
# =============================================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"          # folder of this file = folder of the .py files
NB_DIR="${NB_DIR:-$HERE}"                                      # working directory of the runs (scripts use ROOT=../)
ROOT="${ROOT:-../}"
PY="${PY:-python}"

# pyCATHY: already importable, or a clone given by PYCATHY, or a clone found next to this folder
if [ -z "${PYCATHY:-}" ] && ! $PY -c "import pyCATHY" 2>/dev/null; then
  for c in "$HERE/../pycathy_wrapper" "$HERE/../../pycathy_wrapper" "$HOME/benjmy/pycathy_wrapper" "$HOME/pycathy_wrapper"; do
    [ -d "$c/pyCATHY" ] && { PYCATHY="$c"; break; }
  done
fi
[ -n "${PYCATHY:-}" ] && export PYTHONPATH="$HERE:$PYCATHY${PYTHONPATH:+:$PYTHONPATH}" || export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"

CASE="${CASE:-NC}"                                             # NC (uncompacted, nb1) | HC (compacted, nb2)
NENS="${NENS:-64}"
SEED="${SEED:-1}"
DA="${DA:-enkf}"                                               # enkf | pf | esmda (smoother)
BIAS="${BIAS:-none}"                                           # none (prior centred on the truth) | biased
EXTRA="${EXTRA:-}"                                             # extra nb7 arguments
REAL_NB5="${REAL_NB5:-$HERE/nb5_DA_ERT_soil3D_nc.py}"
REAL_PLOTS="${REAL_PLOTS:-$HERE/nb5_DA_ERT_soil3D_nc_plots.py}"                                   # real-data plotting script (no truth: no score figures)
REAL_NB6="${REAL_NB6:-$HERE/nb6_ESMDA_ERT_soil3D_nc.py}"        # real-data ES-MDA script (smoother)
REAL_SCENARIOS="${REAL_SCENARIOS:-nc_ERT_ic nc_ERT_ic_VGP_pert nc_ERT_ic_VGP_upd nc_ERT_ic_VGP_Archie_pert nc_ERT_ic_VGP_pert_Archie_upd nc_ERT_ic_VGP_Archie_upd}"  # names in scenario_nc_ERT.toml
REAL_DAS="${REAL_DAS:-enkf pf}"                                 # sequential filters of real-seq (nb5): enkf | pf
REAL_PLOT_DAS="${REAL_PLOT_DAS:-enkf pf esmda}"                 # methods plotted by real-plots (folders that do not exist are skipped)
REAL_OPEN_LOOP="${REAL_OPEN_LOOP:-1}"                          # 1 = also one open loop per scenario (like the synthetic matrix)
REAL_RUN_SUFFIX="${REAL_RUN_SUFFIX:-_mat}"                     # "_mat" when VGP_CONTROL = material in nb6, "" for layers
WITH_REAL="${WITH_REAL:-0}"                                    # 1 = 'everything' also runs the real-data case
REAL_JOBS="${REAL_JOBS:-2}"                                    # scenarios run at the same time by real-everything (1 = sequential)
REAL_RESERVE_CORES="${REAL_RESERVE_CORES:-2}"                  # cores left free for the machine; the others are split between the jobs
REAL_MIN_FREE_GB="${REAL_MIN_FREE_GB:-4}"                      # a new job starts only if at least this much RAM is available
REAL_STAGGER="${REAL_STAGGER:-60}"                             # seconds between two job launches (lets RAM ramp up before the next check)

# variables of 'everything'
CASES="${CASES:-NC HC}"
DAS="${DAS:-enkf pf esmda}"
SEEDS="${SEEDS:-$SEED}"
ALL_SCENARIOS="ic ic_VGP ic_VGP_upd ic_Ks_upd ic_Archie_upd full_upd"
SCENARIOS="${SCENARIOS:-$ALL_SCENARIOS}"

usage() {
  sed -n '/^# ---- QUICK START/,/^# ---- SYNTHETIC/p' "${BASH_SOURCE[0]}" | sed '$d' | sed 's/^# \{0,1\}//'
  echo "  Full documentation (every command, variable and scenario): the header of this file."
cat <<USG

Commands, synthetic twin:
  everything | syn-pf | plan | scenarios | truth | da <sc> | openloop <sc> | pair <sc> | matrix <preset>
  compare-scenarios | compare-da | compare-bias | compare-nens | status | plots
Commands, real data:
  real-everything [sc...] | real-pf [sc...] | real-seq [sc...] | real-esmda [sc...] | real-plots [sc...] | real
  help
Scenarios, synthetic (names of make_scenarios.py):
  ic             initial condition only (perturbed, estimated)
  ic_VGP         + van Genuchten perturbed, NOT updated (prior spread only)
  ic_VGP_upd     + van Genuchten perturbed AND updated
  ic_Ks_upd      + Ks updated
  ic_Archie_upd  + Archie (a, m, n) updated
  full_upd       ic + VGP + Ks + Archie, all updated
Scenarios, real data (scenario_nc_ERT.toml; override with REAL_SCENARIOS):
  nc_ERT_ic | nc_ERT_ic_VGP_pert | nc_ERT_ic_VGP_upd | nc_ERT_ic_VGP_Archie_pert
  nc_ERT_ic_VGP_pert_Archie_upd | nc_ERT_ic_VGP_Archie_upd
DA types: enkf, pf (sequential, states + parameters) | esmda (smoother, parameters only; 'ic' estimated as a parameter)
USG
}

cd "$NB_DIR"
cmd="${1:-help}"; shift || true
case "$cmd" in help|-h|--help) usage; exit 0;; esac
if ! $PY -c "import pyCATHY" 2>/dev/null; then
  echo "pyCATHY not importable: set PYCATHY=/path/to/pycathy_wrapper"; exit 1
fi

nb7() {  # nb7 <scenario> [more args]
  local sc="$1"; shift
  $PY "$HERE/nb7_DA_ERT_synthetic.py" --case "$CASE" --scenario "$sc" --da-type "$DA" --nens "$NENS" \
      --seed "$SEED" --prior-bias "$BIAS" --root "$ROOT" $EXTRA "$@"
}
run_status() {  # table of the runs of the case: what each one left on disk (pyCATHY.DA.run_status)
  $PY -m pyCATHY.DA.run_status "${ROOT}simulations" --glob "${STATUS_GLOB:-SYN_${CASE}_*}" "$@" \
    || echo "(pyCATHY.DA.run_status not available: update pyCATHY)"
}
real_plots() {  # real_plots [scenario ...]: figures of every finished real run (method x scenario); missing folders are skipped
  [ -f "$REAL_PLOTS" ] || { echo "set REAL_PLOTS=/path/to/the real-data plotting script"; return 1; }
  local scs=("$@") sc da pre failed=0
  [ ${#scs[@]} -gt 0 ] || read -r -a scs <<< "$REAL_SCENARIOS"
  for sc in "${scs[@]}"; do
    for da in $REAL_PLOT_DAS; do
      case "$da" in esmda) pre="NC_ERT_ESMDA_";; pf) pre="NC_ERT_PF_";; *) pre="NC_ERT_";; esac
      if [ ! -d "${ROOT}simulations/${pre}withDA_${sc}${REAL_RUN_SUFFIX}" ]; then
        echo "[skip] plots $da $sc: no run folder"; continue
      fi
      echo "================ real plots: $da  $sc ================"
      MPLBACKEND="${MPLBACKEND:-Agg}" $PY "$REAL_PLOTS" --scenario "$sc" --da "$da" \
        --vgp-control "$([ "$REAL_RUN_SUFFIX" = "_mat" ] && echo material || echo layers)" \
        || { echo "FAILED plots: $da $sc"; failed=1; }
    done
  done
  return $failed
}
real_loop() {  # real_loop <kind: seq|pf|enkf|esmda> [scenario ...]: one run per scenario (+ open loop) on the real ERT data
  local kind="$1"; shift
  local script scs=("$@") dalist sc da w worder tag prj marker pf failed=0 envv
  case "$kind" in
    esmda)   script="$REAL_NB6"; dalist="esmda" ;;
    pf|enkf) script="$REAL_NB5"; dalist="$kind" ;;
    *)       script="$REAL_NB5"; dalist="$REAL_DAS" ;;
  esac
  [ -f "$script" ] || { echo "script not found: $script (REAL_NB5 / REAL_NB6)"; return 1; }
  [ ${#scs[@]} -gt 0 ] || read -r -a scs <<< "$REAL_SCENARIOS"
  for sc in "${scs[@]}"; do
    for da in $dalist; do
      worder="1 0"; [ "$da" = "enkf" ] && worder="0 1"     # EnKF: open loop (withoutDA) first, then withDA
      for w in $worder; do
        # one open loop per scenario (it does not depend on the DA method): with the first DA of the list only
        if [ "$w" = "0" ]; then
          [ "$REAL_OPEN_LOOP" = "1" ] || continue
          [ "$da" = "$(echo $dalist | cut -d' ' -f1)" ] || continue
          [ "$da" = "pf" ] && continue
        fi
        if [ "$w" = "1" ]; then tag=withDA; else tag=withoutDA; fi
        pf=0
        case "$da" in
          esmda) prj="NC_ERT_ESMDA_${tag}_${sc}${REAL_RUN_SUFFIX}"; marker="_esmda_result.pkl" ;;
          pf)    prj="NC_ERT_PF_${tag}_${sc}${REAL_RUN_SUFFIX}";    marker="_finished.txt"; pf=1 ;;
          *)     prj="NC_ERT_${tag}_${sc}${REAL_RUN_SUFFIX}";       marker="_finished.txt" ;;
        esac
        if [ -f "${ROOT}simulations/$prj/${prj}${marker}" ] && [ "${REAL_FORCE:-0}" != "1" ] && [ "${REAL_DRY_RUN:-0}" != "1" ]; then
          echo "[skip] $da $tag $sc: finished (REAL_FORCE=1 to redo)"; continue
        fi
        echo "================ real $da: $tag  $sc ================"
        envv=(MPLBACKEND="${MPLBACKEND:-Agg}" REAL_SCENARIO="$sc" REAL_WITH_DA="$w" REAL_USE_PF="$pf")
        [ -n "${NENS_REAL:-}" ]  && envv+=(REAL_NENS="$NENS_REAL")
        [ -n "${SEED_REAL:-}" ]  && envv+=(REAL_SEED="$SEED_REAL")
        [ -n "${NITER_REAL:-}" ] && envv+=(REAL_NITER="$NITER_REAL")
        [ "${REAL_DRY_RUN:-0}" = "1" ] && envv+=(REAL_DRY_RUN=1)
        if env "${envv[@]}" $PY "$script"; then :
        elif [ "${REAL_DRY_RUN:-0}" = "1" ]; then :        # nb6 ends a dry run with SystemExit
        else echo "FAILED: $da $tag $sc"; failed=1; fi
      done
    done
  done
  return $failed
}
real_run() {  # real_run <da: enkf|pf|esmda> <w: 1 = withDA | 0 = withoutDA> <scenario>: ONE real run (skipped if finished)
  local da="$1" w="$2" sc="$3" script tag prj marker pf=0 envv
  if [ "$da" = "esmda" ]; then script="$REAL_NB6"; else script="$REAL_NB5"; fi
  [ -f "$script" ] || { echo "script not found: $script (REAL_NB5 / REAL_NB6)"; return 1; }
  if [ "$w" = "1" ]; then tag=withDA; else tag=withoutDA; fi
  case "$da" in
    esmda) prj="NC_ERT_ESMDA_${tag}_${sc}${REAL_RUN_SUFFIX}"; marker="_esmda_result.pkl" ;;
    pf)    prj="NC_ERT_PF_${tag}_${sc}${REAL_RUN_SUFFIX}";    marker="_finished.txt"; pf=1 ;;
    *)     prj="NC_ERT_${tag}_${sc}${REAL_RUN_SUFFIX}";       marker="_finished.txt" ;;
  esac
  if [ -f "${ROOT}simulations/$prj/${prj}${marker}" ] && [ "${REAL_FORCE:-0}" != "1" ] && [ "${REAL_DRY_RUN:-0}" != "1" ]; then
    echo "[skip] $da $tag $sc: finished (REAL_FORCE=1 to redo)"; return 0
  fi
  echo "================ real $da: $tag  $sc ================"
  envv=(MPLBACKEND="${MPLBACKEND:-Agg}" REAL_SCENARIO="$sc" REAL_WITH_DA="$w" REAL_USE_PF="$pf")
  [ -n "${NENS_REAL:-}" ]  && envv+=(REAL_NENS="$NENS_REAL")
  [ -n "${SEED_REAL:-}" ]  && envv+=(REAL_SEED="$SEED_REAL")
  [ -n "${NITER_REAL:-}" ] && envv+=(REAL_NITER="$NITER_REAL")
  [ "${REAL_DRY_RUN:-0}" = "1" ] && envv+=(REAL_DRY_RUN=1)
  if env "${envv[@]}" $PY "$script"; then return 0
  elif [ "${REAL_DRY_RUN:-0}" = "1" ]; then return 0        # nb6 ends a dry run with SystemExit
  else echo "FAILED: $da $tag $sc"; return 1; fi
}
real_phases() {  # real_phases "<da:w> <da:w> ..." [scenario ...]: for each scenario, the phases in the given order
  local spec="$1"; shift
  local scs=("$@") ph da w sc failed=0
  [ ${#scs[@]} -gt 0 ] || read -r -a scs <<< "$REAL_SCENARIOS"
  for sc in "${scs[@]}"; do
    echo "################ scenario $sc ################"
    for ph in $spec; do
      da="${ph%%:*}"; w="${ph##*:}"
      if [ "$w" = "0" ] && [ "$REAL_OPEN_LOOP" != "1" ]; then echo "[skip] $ph $sc (REAL_OPEN_LOOP=0)"; continue; fi
      real_run "$da" "$w" "$sc" || failed=1
    done
  done
  return $failed
}
real_mem_ok() {  # enough available RAM (>= REAL_MIN_FREE_GB) to start another job?
  local kb; kb=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo 2>/dev/null || true)
  [ -z "$kb" ] || [ "$kb" -ge $(( REAL_MIN_FREE_GB * 1048576 )) ]
}
real_parallel() {  # real_parallel "<da:w> ..." [scenario ...]: scenarios in parallel (REAL_JOBS at a time), phases in order inside each
  local spec="$1"; shift
  local scs=("$@") sc k slot usable per cpus logdir failed=0 ncpu
  local -a slot_sc=()
  [ ${#scs[@]} -gt 0 ] || read -r -a scs <<< "$REAL_SCENARIOS"
  if [ "$REAL_JOBS" -le 1 ] || [ ${#scs[@]} -le 1 ]; then real_phases "$spec" "${scs[@]}"; return $?; fi
  ncpu=$(nproc); usable=$(( ncpu - REAL_RESERVE_CORES )); [ "$usable" -ge 1 ] || usable=1
  [ "$REAL_JOBS" -le "$usable" ] || { echo "REAL_JOBS=$REAL_JOBS > $usable usable cores: using $usable"; REAL_JOBS=$usable; }
  per=$(( usable / REAL_JOBS ))
  logdir="${ROOT}simulations/real_logs"; mkdir -p "$logdir"
  echo "parallel: $REAL_JOBS jobs x $per cores (cores 0-$(( REAL_JOBS * per - 1 )) of $ncpu; $(( ncpu - REAL_JOBS * per )) left free), logs: $logdir/<scenario>.log"
  for sc in "${scs[@]}"; do
    while :; do                                                  # wait for a free slot and enough RAM
      slot=-1
      for (( k=0; k<REAL_JOBS; k++ )); do
        if [ -z "${slot_sc[$k]:-}" ] || [ -f "$logdir/${slot_sc[$k]}.rc" ]; then slot=$k; break; fi
      done
      if [ "$slot" -ge 0 ] && real_mem_ok; then break; fi
      sleep 10
    done
    slot_sc[$slot]="$sc"
    cpus="$(( slot * per ))-$(( slot * per + per - 1 ))"
    rm -f "$logdir/$sc.rc"
    echo "[start] $sc on cores $cpus ($(date +%H:%M:%S))"
    (
      taskset -cp "$cpus" "$BASHPID" >/dev/null 2>&1 || echo "(taskset unavailable: no CPU pinning)"
      renice -n 10 -p "$BASHPID" >/dev/null 2>&1 || true         # lower priority: the machine stays responsive
      export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$per}" MKL_NUM_THREADS="${MKL_NUM_THREADS:-$per}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-$per}"
      rc=0; real_phases "$spec" "$sc" > "$logdir/$sc.log" 2>&1 || rc=$?
      echo "$rc" > "$logdir/$sc.rc"
      echo "[finished] $sc rc=$rc ($(date +%H:%M:%S))"          # console line as soon as this scenario ends
    ) &
    sleep "$REAL_STAGGER"
  done
  wait
  for sc in "${scs[@]}"; do
    if [ "$(cat "$logdir/$sc.rc" 2>/dev/null || echo 1)" = "0" ]; then echo "[done] $sc"; else echo "[FAILED] $sc (see $logdir/$sc.log)"; failed=1; fi
  done
  return $failed
}
real_has() { case " $REAL_DAS " in *" $1 "*) return 0;; esac; return 1; }   # is <enkf|pf> in REAL_DAS ?
truth_file() { echo "${ROOT}data/3d/synthetic/$1/truth_$1.pkl"; }
matrix() { $PY "$HERE/run_matrix.py" --root "$ROOT" --prior-bias "$BIAS" ${EXTRA:+--extra "$EXTRA"} "$@"; }

case "$cmd" in
  plan)
    echo "everything = for CASES=[$CASES]:  scenarios + truth (if missing) ; then every run of"
    echo "  SCENARIOS=[$SCENARIOS] x DAS=[$DAS] x SEEDS=[$SEEDS], N=$NENS, prior=$BIAS, + one open loop per scenario ; then plots"
    for c in $CASES; do matrix --cases "$c" --scenarios $SCENARIOS --da $DAS --seeds $SEEDS --nens "$NENS" --open-loop --dry-run | grep -E "^\[|runs" ; done
    usage | tail -12 ;;
  everything)
    for c in $CASES; do
      echo "================ case $c ================"
      if [ ! -f "$(truth_file "$c")" ] && [[ " $* " != *" --dry-run "* ]]; then
        CASE="$c" "$0" truth
      fi
      matrix --cases "$c" --scenarios $SCENARIOS --da $DAS --seeds $SEEDS --nens "$NENS" --open-loop "$@" || echo "some runs of $c failed (see logs/)"
    done
    if [[ " $* " != *" --dry-run "* ]]; then
      if [ "$WITH_REAL" = "1" ]; then
        "$0" real-seq || echo "real-data EnKF/PF run(s) failed"
        "$0" real-esmda || echo "real-data ES-MDA run(s) failed"
      fi
      for c in $CASES; do REAL=0 CASE="$c" "$0" plots || true; done      # synthetic figures (REAL=0: real-data ones are made once, below)
      if [ -n "$REAL_PLOTS" ]; then "$0" real-plots || true; fi
    fi ;;
  syn-pf)            # synthetic twin, PF only: 'everything' restricted to DAS=pf (truth if missing, runs + open loops, figures)
    DAS="pf" "$0" everything "$@" ;;
  scenarios)
    flag=""; [ "$BIAS" = "none" ] && flag="--no-bias"
    $PY "$HERE/make_scenarios.py" $flag ${KZ_RATIO:+--kz-ratio} "$@" ;;
  truth)
    $PY "$HERE/nb6_synth_truth_ERT.py" --case "$CASE" --root "$ROOT" --noise-pct "${NOISE:-5}" "$@" ;;
  da)        [ $# -ge 1 ] || { usage; exit 1; }; sc="$1"; shift; nb7 "$sc" "$@" ;;
  openloop)  [ $# -ge 1 ] || { usage; exit 1; }; sc="$1"; shift; DA=enkf nb7 "$sc" --no-da "$@" ;;
  pair)      [ $# -ge 1 ] || { usage; exit 1; }; sc="$1"; shift; nb7 "$sc" "$@"; DA=enkf nb7 "$sc" --no-da "$@" ;;
  matrix)    p="${1:-quick}"; shift || true; matrix --preset "$p" --open-loop "$@" ;;
  compare-scenarios)
    matrix --cases "$CASE" --scenarios ${SCENARIOS:-$ALL_SCENARIOS} --da "$DA" --seeds "$SEED" --nens "$NENS" --open-loop "$@" ;;
  compare-da)
    matrix --cases "$CASE" --scenarios "${SCEN:-ic_Ks_upd}" --da enkf pf esmda --seeds "$SEED" --nens "$NENS" --open-loop "$@" ;;
  compare-bias)
    for b in none biased; do
      BIAS="$b" matrix --cases "$CASE" --scenarios "${SCEN:-ic_Ks_upd}" --da "$DA" --seeds "$SEED" --nens "$NENS" --open-loop "$@"
    done ;;
  compare-nens)
    matrix --cases "$CASE" --scenarios "${SCEN:-ic_Ks_upd}" --da "$DA" --seeds "$SEED" --nens ${NENS_LIST:-16 32 64} "$@" ;;
  status) run_status "$@" ;;
  plots)
    run_status || true                                          # which runs can be plotted, and why not (one line each)
    a=(--case "$CASE" --root "$ROOT")
    [ -n "${RUNS:-}" ] && a+=(--runs $RUNS)
    [ -n "${INCLUDE_RUNNING:-}" ] && a+=(--include-running)
    rc=0; $PY "$HERE/nb8_plots_synthetic.py" "${a[@]}" "$@" || rc=$?
    if [ -n "$REAL_PLOTS" ] && [ "${REAL:-1}" != "0" ]; then real_plots || true; fi
    exit $rc ;;
  real-plots) real_plots "$@" ;;
  real-everything)   # real data only, per scenario: PF, then EnKF, last ES-MDA (small ensemble); then the figures
    rc=0
    phases=()
    real_has pf   && phases+=(pf:0 pf:1)
    real_has enkf && phases+=(enkf:0 enkf:1)
    [ ${#phases[@]} -gt 0 ] && real_parallel "${phases[*]}" "$@" || true
    # ES-MDA last, with a smaller ensemble (default 6; override with NENS_ESMDA_REAL)
    NENS_REAL="${NENS_ESMDA_REAL:-6}" real_parallel "esmda:1 esmda:0" "$@" || { echo "some real runs failed"; rc=1; }
    if [ -n "$REAL_PLOTS" ]; then real_plots "$@" || rc=1; fi
    exit $rc ;;
  real-pf)           # real data, PF only, per scenario: withoutDA (open loop) then withDA; scenarios REAL_JOBS at a time; no ES-MDA
    rc=0
    real_parallel "pf:0 pf:1" "$@" || { echo "some real PF runs failed"; rc=1; }
    if [ -n "$REAL_PLOTS" ]; then REAL_PLOT_DAS="pf" real_plots "$@" || rc=1; fi
    exit $rc ;;
  real-esmda) real_loop esmda "$@" ;;
  real-seq)   real_loop seq "$@" ;;
  real)
    [ -f "$REAL_NB5" ] || { echo "set REAL_NB5=/path/to/the real-data notebook script (nb5 EnKF/PF or nb6_ESMDA)"; exit 1; }
    $PY "$REAL_NB5" "$@" ;;
  *) echo "unknown command: $cmd"; usage; exit 1 ;;
esac