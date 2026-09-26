# Demo 07: T-Rex+GVS on Heterogeneous ARMA Blocks

## Purpose

Stress-test T-Rex+GVS when different blocks follow **different time-series structures** rather than a single uniform correlation model. This is the most heterogeneous of the four-block scenarios: each active block has its own ARMA specification, testing whether the selector can still recover grouped signal when "the group correlation model" is not one fixed thing.

## Data-generating process

`dgp_arma_blocks_snr` — four blocks of sizes `{20, 50, 80, 65}`, each a distinct time-series process, shuffled into `p = 500` with white-noise gaps:

- **Block 1** (20 vars, active): AR(1), coefficient `ar_coef`.
- **Block 2** (50 vars, active): MA(3), `ma_coefs = c(0.5, 0.3, 0.1)` (fixed).
- **Block 3** (80 vars, active): ARMA(2,1) (`s = 150` over the three active blocks).
- **Block 4** (65 vars, inactive trap): AR(1), coefficient `ar_coef`.

The swept correlation knob is `ar_coef` (shared by the AR/ARMA blocks) — **not** `rho`, since the MA(3) block has no AR component. Dimensions: `n = 200`, `p = 500`.

## Methods compared

- **TENET** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **TENET_AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **TIENET_AUG** — `GVS_type = "IEN"` (`en_solver` not applicable; left at the library default)

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "Single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — 2-D SNR x `ar_coef` grid (5 SNR x 6 ar_coef): `SNR in {0.2, 0.5, 1.0, 2.0, 5.0}`, `ar_coef in {0.30, 0.50, 0.70, 0.80, 0.90, 0.95}` (labeled `ac=` in the result tables).

## Running

```bash
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_07_mc_sim_arma_blocks/demo_trex_gvs_07_mc_sim_arma_blocks.R      # all parts
Rscript R/trex_selector_methods/trex_gvs/demo_trex_gvs_07_mc_sim_arma_blocks/demo_trex_gvs_07_mc_sim_arma_blocks.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` (`gvs_arma_blocks_2d.txt` + `.csv`); the script sources the shared `../trex_gvs_dgps.R`, `../trex_gvs_sim_common.R`, and `../../simulation_utils.R`.

## Interpretation

- Because the three active blocks have genuinely different correlation structures (AR(1), MA(3), ARMA(2,1)), recovery is expected to be uneven across blocks — aggregate TPP is a blend of these different per-block difficulties.
- FDP should stay reasonably controlled near or below `tFDR = 0.1`; the trap block (AR(1), same `ar_coef` as the active AR(1) block) is structurally similar to an active block and is the main test of false-group exclusion.
- Compare against Demo 06 (uniform AR(1) blocks) to see whether cross-block heterogeneity costs additional power beyond what AR(1) decay alone costs.
