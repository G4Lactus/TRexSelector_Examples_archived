# Demo 06: T-Rex+GVS on AR(1) Blocks

## Purpose

Evaluate the Terminating Random Experiments (T-Rex) selector with Grouped Variable Selection (T-Rex+GVS) when
within-block correlation follows a **first-order autoregressive (AR(1)) time-series structure** instead of
equicorrelation. This uses the same four-block active/trap geometry as Demo 03, but replaces the
shared-latent-factor model with an AR(1)/Toeplitz covariance within each block, testing robustness to realistic,
decaying correlation patterns.

## Data-generating process

`dgp_ar1_blocks_snr` — four blocks of sizes `{20, 50, 80, 65}`: three active (`beta = 3`, `s = 150`) plus one
inactive AR(1) trap, shuffled into `p = 500` with independent and identically distributed (i.i.d.) `N(0,1)` gaps.
Within-block correlation is AR(1)/Toeplitz: `Cor(X_j, X_k) = rho^|j-k|`, so correlation decays with column distance
and a block no longer behaves like one tight cluster. The AR(1) parameter `rho` is passed directly to the
data-generating process (DGP). Dimensions: `n = 200`, `p = 500`.

## Methods compared

- **TENET** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **TENET_AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **TIENET_AUG** — `GVS_type = "IEN"` (`en_solver` not applicable, left at the library default)

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "Single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — 2-D signal-to-noise ratio (SNR) x rho grid (5 SNR x 6 rho): `SNR in {0.2, 0.5, 1.0, 2.0, 5.0}`,
  `rho in {0.30, 0.50, 0.70, 0.85, 0.90, 0.99}` (rho enters the DGP directly).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_06_mc_sim_ar1_blocks/demo_trex_gvs_06_mc_sim_ar1_blocks.R      # all parts
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_06_mc_sim_ar1_blocks/demo_trex_gvs_06_mc_sim_ar1_blocks.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` (`gvs_ar1_blocks_2d.txt` + `.csv`). The script sources the shared
`../trex_gvs_dgps.R`, `../trex_gvs_sim_common.R`, and `../../../simulation_utils.R`.

## Interpretation

- The false discovery proportion (FDP) should stay low (often well below `tFDR = 0.1`) across the SNR axis for all
  three methods, but the true positive proportion (TPP) should recover much more slowly with SNR than in the
  equicorrelated Hastie/Mixed-Blocks demos — AR(1) decay means only nearby columns are strongly correlated, weakening
  group-level evidence aggregation.
- TIENET_AUG is expected to sit between its near-zero TPP on the equicorrelated Gaussian demos and the methods based
  on the elastic net (EN), showing real (if smaller) power gains as SNR increases.
- The rho axis is the most informative view, because it shows how quickly recovery degrades as `rho` decreases and
  the AR(1) decay steepens.
