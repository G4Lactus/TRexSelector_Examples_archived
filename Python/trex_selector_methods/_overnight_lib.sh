#!/usr/bin/env bash
# ==============================================================================
# _overnight_lib.sh  —  shared runner for the per-suite overnight scripts.
# ==============================================================================
#
# Sourced by each suite's run_overnight.sh. Runs a list of demos end-to-end to
# (re)generate their simulation_results/ data, then renders any correlation
# heatmaps present in the suite. Designed for unattended overnight use:
#
#   * continue-on-error — one failing demo never aborts the batch;
#   * per-demo timestamped logs + a run summary under <suite>/_overnight_logs/;
#   * SKIP_HEAVY=1 skips the demos listed in HEAVY_DEMOS (RAM/time hogs).
#
# Not sourced by anything type-checked; pure bash.
# ==============================================================================

_find_venv() {
  local d="$1"
  while [[ "$d" != "/" ]]; do
    if [[ -x "$d/.venv/bin/python" ]]; then echo "$d/.venv/bin/python"; return 0; fi
    d="$(dirname "$d")"
  done
  return 1
}

# run_suite <suite_dir> <demo_name> [<demo_name> ...]
# Each <demo_name> is a subdir containing <demo_name>.py.
run_suite() {
  local suite_dir="$1"; shift
  local suite_name; suite_name="$(basename "$suite_dir")"

  local venv; venv="$(_find_venv "$suite_dir")" || {
    echo "error: no .venv/bin/python found above $suite_dir" >&2
    echo "       create it and install the demo deps, then re-run." >&2
    return 1
  }

  local logdir="$suite_dir/_overnight_logs"; mkdir -p "$logdir"
  local stamp; stamp="$(date +%Y%m%d_%H%M%S)"
  local summary="$logdir/summary_${stamp}.txt"
  local skip_heavy="${SKIP_HEAVY:-0}"

  {
    echo "=============================================================="
    echo " overnight run: $suite_name"
    echo " started : $(date)"
    echo " venv    : $venv"
    echo " SKIP_HEAVY=$skip_heavy   HEAVY_DEMOS='${HEAVY_DEMOS:-}'"
    echo "=============================================================="
    echo
  } | tee "$summary"

  local name dir log t0 rc ok=0 fail=0 skip=0
  for name in "$@"; do
    dir="$suite_dir/$name"
    if [[ ! -f "$dir/$name.py" ]]; then
      echo "[MISS ] $name  (no $name.py)" | tee -a "$summary"; continue
    fi
    if [[ "$skip_heavy" == "1" && " ${HEAVY_DEMOS:-} " == *" $name "* ]]; then
      echo "[SKIP ] $name  (HEAVY, SKIP_HEAVY=1)" | tee -a "$summary"
      skip=$((skip + 1)); continue
    fi
    log="$logdir/${name}_${stamp}.log"
    echo "[START] $name  $(date +%H:%M:%S)" | tee -a "$summary"
    t0=$SECONDS
    ( cd "$dir" && "$venv" "$dir/$name.py" ) >"$log" 2>&1; rc=$?
    if [[ $rc -eq 0 ]]; then
      echo "[ OK  ] $name  $((SECONDS - t0))s  -> $(basename "$log")" | tee -a "$summary"
      ok=$((ok + 1))
    else
      echo "[FAIL ] $name  $((SECONDS - t0))s  rc=$rc  -> $(basename "$log")" | tee -a "$summary"
      tail -n 4 "$log" 2>/dev/null | sed 's/^/         | /' | tee -a "$summary"
      fail=$((fail + 1))
    fi
  done

  echo | tee -a "$summary"
  echo "--- design-correlation heatmaps (one per DGP) ---" | tee -a "$summary"
  local corr_sh="$suite_dir/generate_corr_plots.sh"
  if [[ -f "$corr_sh" ]]; then
    if bash "$corr_sh" >>"$summary" 2>&1; then
      echo "[corr ] $(basename "$suite_dir") DGP heatmaps rendered" | tee -a "$summary"
    else
      echo "[corr!] $(basename "$suite_dir") (failed — see summary)" | tee -a "$summary"
    fi
  else
    echo "(no generate_corr_plots.sh in this suite yet)" | tee -a "$summary"
  fi

  {
    echo
    echo "=== $suite_name done: $ok ok, $fail failed, $skip skipped ==="
    echo " finished: $(date)"
    echo " summary : $summary"
  } | tee -a "$summary"

  return 0
}
