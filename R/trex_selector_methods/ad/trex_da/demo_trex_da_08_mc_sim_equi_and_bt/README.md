# Demo 08: T-Rex+DA+EQUI on Equicorrelated Data

## Purpose

Cover the equicorrelation (EQUI), or compound-symmetry, dependency structure, where every predictor is tied to the same
shared latent factor, using the `EQUI` dependency-aware (DA) method. Unlike first-order autoregressive (AR(1))
correlation, equicorrelation has no notion of column distance, so the correction acts uniformly across all variables.
The single 2D sweep locates the DA suppression cliff in rho and shows its signal-to-noise ratio (SNR) independence.

## Data-generating process

`dgp_equi_snr` builds a pure equicorrelated design via a single latent factor,
`X[i,j] = sqrt(rho) * f[i] + sqrt(1 - rho) * eta[i,j]`. Monte Carlo config:
`n=300, p=1000, s=10, amplitude=3.0, tFDR=0.2, K=20, num_MC=200, seed=2026`,
Random support (redrawn per trial). Per-SNR-column seed blocks:
`col_seed = 2026 + si * 1000000`.

## DA method

`TRexDASelector$new(..., da_control = trex_da_control(da_method = "EQUI"), control = trex_control(...))` with the
default solver set from `da_sim_common.R`: `TREX-DA+EQUI: TLARS` / `TREX-DA+EQUI: TAFS` / `TREX-DA+EQUI: TOMP`, plus
the base (no-DA) `TREX: TLARS/TAFS/TOMP` comparison rows.

## What it runs

- Part 1: 2D SNR x rho sweep — for each SNR column in `{0.5, 1, 2, 5}`, a full rho sweep over
  `{0, 0.025, 0.05, 0.1, 0.2, 0.4, 0.7, 0.9}` (dense near 0 to locate the DA suppression cliff), 200 Monte Carlo trials
  per point.

Earlier revisions of this demo also carried binary-tree (BT) hierarchical-block SNR and linkage sweeps. They were
dropped as redundant — the BT method is covered on better-suited block designs in Demos 02-05. (The folder keeps the
cpp name `equi_and_bt` for parity.)

Result stems (`simulation_results/data/`): one `.txt`/`.csv` pair per SNR column,
`da_trex_mc_da_equi_rho_snr{0p5,1,2,5}`.

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_08_mc_sim_equi_and_bt/demo_trex_da_08_mc_sim_equi_and_bt.R      # full run
Rscript R/trex_selector_methods/ad/trex_da/demo_trex_da_08_mc_sim_equi_and_bt/demo_trex_da_08_mc_sim_equi_and_bt.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` as `.txt` and `.csv` files. The
script sources the consolidated `../trex_da_dgps.R`, the shared
`../da_sim_common.R`, and `../../../support_generators.R`.

## Interpretation

- Equicorrelation ties every predictor to one shared factor, so DA-EQUI suppresses all candidates as soon as rho is
  noticeably above 0, independent of signal strength, while base Terminating Random Experiments (T-Rex) loses false
  discovery rate (FDR) control at the same tiny rho.
- The rho sweep isolates how the strength of the shared factor drives selection difficulty. At rho=0 the design is
  effectively independent, and the correction has little to do — the rho=0 column is the working anchor.
- Reading the four SNR columns side by side shows the cliff location and its SNR independence in a single figure.
