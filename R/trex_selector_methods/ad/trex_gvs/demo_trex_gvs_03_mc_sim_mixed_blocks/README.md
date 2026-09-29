# Demo 03: T-Rex+GVS on Mixed Blocks

## Purpose

Evaluate the Terminating Random Experiments (T-Rex) selector with Grouped Variable Selection (T-Rex+GVS) on a more
realistic grouped design. Four contiguous equicorrelated blocks of unequal size are placed in **random order** with
**random white-noise gaps** between and around them, including one **inactive equicorrelated trap block**. A step up
in realism from Demo 01's tidy contiguous Hastie blocks.

## Data-generating process

`dgp_mixed_blocks_snr` — four equicorrelated blocks of sizes `{20, 50, 80, 65}`: three active (`beta = 3`,
`s = 150`) plus one inactive trap block that shares the active blocks' correlation structure. Blocks are shuffled into
a random order and separated by random-width white-noise gaps within `p = 500`. Each block column is
`X_j = Z_g + sd_x * xi`, `rho = 1/(1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`. Default `sd_x = sqrt(0.01)` gives
`rho ~ 0.99`.

## Methods compared

- **TENET** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **TENET_AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **TIENET_AUG** — `GVS_type = "IEN"` (`en_solver` not applicable, left at the library default)

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "Single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — 2-D signal-to-noise ratio (SNR) x rho grid (5 SNR x 7 rho): `SNR in {0.2, 0.5, 1.0, 2.0, 5.0}`,
  `rho in {0.30, 0.50, 0.70, 0.85, 0.90, 0.95, 0.99}` (`sd_x` derived from rho, and the rho grid is the union of the
  two former presets' grids).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_03_mc_sim_mixed_blocks/demo_trex_gvs_03_mc_sim_mixed_blocks.R      # all parts
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_03_mc_sim_mixed_blocks/demo_trex_gvs_03_mc_sim_mixed_blocks.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` (`gvs_mixed_blocks_2d.txt` + `.csv`). The script sources the
shared `../trex_gvs_dgps.R`, `../trex_gvs_sim_common.R`, and `../../../simulation_utils.R`.

## Interpretation

- The false discovery proportion (FDP) should stay controlled near `tFDR = 0.1` for TENET / TENET_AUG across the
  grid, since the inactive trap block should be excluded once its group is identified as null.
- The true positive proportion (TPP) should climb with SNR, while the random gaps and block order mainly test that the
  grouping logic is invariant to block placement rather than materially changing detection power.
- Comparing the high-rho columns (`rho ~ 0.99`) against the mid-rho columns (`rho = 0.70`–`0.85`) shows how
  performance degrades as within-block correlation decreases, all else equal.
- As in Demos 01–02, expect TIENET_AUG to be substantially more conservative than the methods based on the elastic
  net (EN).
