# Demo 01: T-Rex+GVS on Hastie Equicorrelated Blocks

## Purpose

Evaluate grouped-signal recovery under strong equicorrelation, using the classical three-group design of Zou & Hastie (2005). This is the canonical starting point for the T-Rex+GVS suite: every variable inside an active group is itself active, so the question is purely **group-level recovery**, not within-group sparse-support identification.

## Data-generating process

`dgp_hastie_snr` — three equicorrelated latent-factor groups of 50 variables each (columns 1–50, 51–100, 101–150), all active with `beta = 3`, plus 350 i.i.d. background columns. Each column is `X_j = Z_g + sd_x * xi`, giving within-group correlation `rho = 1/(1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`, `s = 150`. Default `sd_x = sqrt(0.01)` gives `rho ~ 0.99`; SNR is calibrated through `sigma_y`.

## Methods compared

- **TENET** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **TENET_AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **TIENET_AUG** — `GVS_type = "IEN"` (`en_solver` not applicable; left at the library default)

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "Single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — 2-D SNR x rho grid (5 SNR x 6 rho): `SNR in {0.2, 0.5, 1.0, 2.0, 5.0}`, `rho in {0.30, 0.50, 0.70, 0.90, 0.95, 0.99}` (`sd_x` derived from rho).
- **Part 2** — `lambda2_method` comparison, `CV_1SE_SVD` vs `CV_1SE_CCD` (TENET and TIENET_AUG), fixed `SNR = 2.0`, `rho = 0.7`.
- **Part 3** — `hc_linkage` comparison, Single / Complete / Average / WPGMA (TENET and TIENET_AUG), same operating point, `lambda2 = CV_1SE_CCD`.

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_01_mc_sim_hastie_en_blocks/demo_trex_gvs_01_mc_sim_hastie_en_blocks.R      # all parts
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_01_mc_sim_hastie_en_blocks/demo_trex_gvs_01_mc_sim_hastie_en_blocks.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` (`gvs_Hastie_2d`, `gvs_Hastie_lambda2_method`, `gvs_Hastie_hc_linkage`, each as `.txt` + `.csv`); the `run_part_*` flags at the top toggle individual parts; the script sources the shared `../trex_gvs_dgps.R`, `../trex_gvs_sim_common.R`, and `../../simulation_utils.R`.

## Interpretation

- FDP should stay near the `tFDR = 0.1` target across the SNR axis for TENET and TENET_AUG, with TPP rising from weak recovery at low SNR toward near-complete recovery once `SNR >= 1`.
- TIENET_AUG is expected to behave much more conservatively on this all-active design (very low FDP, slowly rising TPP); this pattern recurs across the other equicorrelated-block demos.
- The rho axis should be fairly stable across correlation levels, since HAC-based grouping aggregates evidence across strongly correlated columns.
- Parts 2–3 check that the qualitative conclusions are not sensitive to the `lambda2` selection rule or the HAC linkage choice.
