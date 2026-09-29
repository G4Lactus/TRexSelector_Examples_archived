# Demo 05: T-Rex+GVS on Heavy-Tailed (Student-t(3)) Blocks

## Purpose

Stress-test T-Rex+GVS on non-Gaussian data. This demo reuses the same four-block active/trap geometry as Demo 03, but replaces every Gaussian component of the data-generating process — latent factors, within-block perturbations, and response noise — with scaled Student-t(3) variables, testing robustness to heavy-tailed departures from the Gaussian assumption.

## Data-generating process

`dgp_t3_equi_blocks_snr` — four equicorrelated blocks of sizes `{20, 50, 80, 65}`: three active (`beta = 3`, `s = 150`) plus one inactive heavy-tailed trap, shuffled into `p = 500` with white-noise gaps. Columns follow `X_j = sqrt(rho) * Z_g + sqrt(1 - rho) * E`, where the factors `Z`, the perturbations `E`, and the response noise are drawn from `t(df) / sqrt(df)` with `df = 3` (rescaled to unit variance). Within-block equicorrelation `rho` is passed directly to the DGP. Dimensions: `n = 200`, `p = 500`.

## Methods compared

- **TENET** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **TENET_AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **TIENET_AUG** — `GVS_type = "IEN"` (`en_solver` not applicable; left at the library default)

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "Single"`, `df = 3`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — 2-D SNR x rho grid (5 SNR x 6 rho): `SNR in {0.2, 0.5, 1.0, 2.0, 5.0}`, `rho in {0.30, 0.50, 0.70, 0.85, 0.90, 0.99}` (rho enters the DGP directly).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_05_mc_sim_t3_blocks/demo_trex_gvs_05_mc_sim_t3_blocks.R      # all parts
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_05_mc_sim_t3_blocks/demo_trex_gvs_05_mc_sim_t3_blocks.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` (`gvs_t3_blocks_2d.txt` + `.csv`); the script sources the shared `../trex_gvs_dgps.R`, `../trex_gvs_sim_common.R`, and `../../simulation_utils.R`.

## Interpretation

- TENET and TENET_AUG should remain reasonably FDP-controlled with TPP rising toward the mid-0.9s as SNR grows, though slightly below the Gaussian equicorrelated case (Demo 03) at matched SNR — heavy tails inject occasional large outliers that blur group boundaries.
- Unlike the Gaussian equicorrelated demos, TIENET_AUG is expected to be much less conservative here: FDP tightening below target as SNR grows while TPP climbs into the 0.7–0.8 range — a more competitive pattern.
- Key question over the SNR axis: does heavy-tailed noise mainly cost power (lower TPP), or also compromise FDP control (drifting above 0.1)?
