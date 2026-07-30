#!/usr/bin/env bash
#
# Regenerate this demo's scalability dashboard via the suite plotting module
#   ../trex_plt_utils.R  (ggplot2, 'scalability' mode)
# Figures are written to simulation_results/plots/.
#
# Demo 08's CSV uses the wider n,p,is_barrier,solver,num_mc,snr,metric,value
# schema rather than the suite's tidy solver,metric,snr,value, so the plotter
# handles it in a dedicated 'scalability' mode: a 2x2 dashboard of runtime,
# memory, the out-of-core money plot and FDR control vs problem size.
#
# Usage:
#   ./generate_plots.sh                          # this demo's default CSV
#   ./generate_plots.sh <other.csv>              # explicit CSV (e.g. smoke)
#   ./generate_plots.sh --formats png            # extra args pass through to
#   ./generate_plots.sh --tfdr 0.05              # trex_plt_utils.R
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

rscript="${RSCRIPT:-$(command -v Rscript || true)}"
if [[ -z "$rscript" ]]; then
  echo "error: no Rscript found on PATH (set RSCRIPT=/path/to/Rscript)" >&2
  exit 1
fi

data="$here/simulation_results/data"

# First non-flag argument overrides this demo's default CSV; fall back to the
# smoke-run CSV when the full benchmark has not been recorded yet.
csv="$data/demo_trex_08_scalability_results.csv"
if [[ $# -gt 0 && "$1" != -* ]]; then
  csv="$1"
  shift
elif [[ ! -f "$csv" && -f "$data/demo_trex_08_scalability_smoke.csv" ]]; then
  echo "note: full benchmark CSV not found - plotting the smoke-run CSV."
  csv="$data/demo_trex_08_scalability_smoke.csv"
fi

stem="$(basename "$csv" .csv)_dashboard"

exec "$rscript" "$here/../trex_plt_utils.R" scalability "$csv" \
  --title 'T-Rex scalability on the memory-mapped pipeline (TLARS / TOMP / TAFS)' \
  --stem "$stem" \
  --outdir "$here/simulation_results/plots" \
  "$@"
