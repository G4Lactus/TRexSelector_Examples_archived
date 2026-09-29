#!/usr/bin/env bash
#
# Overnight: (re)generate the T-Rex SPCA suite's result data.
#
# Command:
#   ./run_overnight.sh              # run everything (default)
#   SKIP_HEAVY=1 ./run_overnight.sh # no-op here (no HEAVY demos in this suite)
#
# The demo takes ~2-4 h at num_MC=200 over the 9-point SNR grid with 6 workers
# (mclapply, one fork per trial). Figures are regenerated separately via the
# demo's generate_plots.sh.
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$here/../../_overnight_lib.sh"

export HEAVY_DEMOS=""

run_suite "$here" \
  demo_trex_spca_01_mc_sim
