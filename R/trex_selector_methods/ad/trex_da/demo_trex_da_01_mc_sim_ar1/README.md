# Demo 01: T-Rex+DA+AR1 on AR(1) Toeplitz Data

## Purpose

The foundational dependency-aware scenario: recover a sparse support under AR(1)
Toeplitz column correlation using the `AR1` DA method. The demo runs Monte Carlo
sweeps over SNR and rho, plus a 2D gap-vs-rho study that probes exactly when the
correction window starts suppressing true signals.

## Data-generating process

`dgp_ar1_snr` builds a Toeplitz design where `Sigma[j,k] = rho^|j-k|` via the
recursion `X[i,j] = rho * X[i,j-1] + sqrt(1 - rho^2) * eta[i,j]`, then calibrates
noise to a target SNR. Monte Carlo config: `n=300, p=1000, s=10, amplitude=3.0,
rho=0.7, SNR=2.0, tFDR=0.2, K=20, num_MC=200, seed=2026, max_gap=20`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "AR1"),
control = trex_control(...))` with the default solver set from
`da_sim_common.R`: `TREX-DA+AR1: TLARS` / `TREX-DA+AR1: TAFS` /
`TREX-DA+AR1: TOMP` (greedy solvers with exchangeable tie-breaking
`exch_tie_alpha=0.25` and stagnation control). Parts 1-2 also run the base
(no-DA) `TREX: TLARS/TAFS/TOMP` comparison rows.

## What it runs

- Part 1: SNR sweep `{0.1, 0.2, 0.5, 1.0, 2.0, 5.0}`, CappedSpread(max_gap=20) and Random support, DA + base T-Rex rows.
- Part 2: rho sweep `seq(0.0, 0.9, by = 0.1)` at SNR=2.0, CappedSpread and Random support, DA + base T-Rex rows.
- Part 3: 2D gap x rho sweep; `gap_grid = {100, 50, 20, 15, 10, 5, 1}`, `rho_grid = {0.0..0.9}`, with the kappa-boundary `kappa = ceiling(log(0.02)/log(rho))`.

Result stems (`simulation_results/data/`): `da_trex_mc_da_ar1_snr_{capped,random}`,
`da_trex_mc_da_ar1_rho_{capped,random}`, `da_trex_mc_da_ar1_gap_rho`
(each as `.txt` + `.csv`).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_01_mc_sim_ar1/demo_trex_da_01_mc_sim_ar1.R      # full run
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_01_mc_sim_ar1/demo_trex_da_01_mc_sim_ar1.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` as TXT and CSV files; the
script sources the consolidated `../trex_da_dgps.R`, the shared
`../da_sim_common.R`, and `../../support_generators.R`.

## Interpretation

- The AR1 correction substantially reduces FDR relative to the classical (no-DA) T-Rex on the same correlated data, at a real TPR cost as SNR grows since some correlated neighbors of true positives are actively excluded.
- The 2D gap x rho study is the key diagnostic: once active variables sit closer together than the kappa-boundary implied by rho, they fall inside each other's correction windows and TPP collapses. This is a structural property, not a bug.
- Comparing CappedSpread vs Random support shows whether deterministic even spacing behaves differently from random placement at the same nominal spacing.
