# Demo 02: T-Rex+GVS on Scattered Blocks

## Purpose

Test whether spatially **scattering** the members of a correlated group across the design matrix — rather than keeping them contiguous — affects T-Rex+GVS group recovery. This isolates the effect of column layout from correlation strength: the group structure is identical to Demo 01, only the column positions are permuted.

## Data-generating process

`dgp_scattered_grouped_snr` — `n_groups = 3` equicorrelated latent-factor groups of `group_size = 50` variables each (all active, `beta = 3`, `s = 150`), whose members are assigned by a random column permutation instead of contiguous blocks. Background columns are i.i.d. `N(0,1)`. Each grouped column is `X_j = Z_g + sd_x * xi`, with within-group correlation `rho = 1/(1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`. Default `sd_x = sqrt(0.01)` gives `rho ~ 0.99`.

## Methods compared

- **TENET** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **TENET_AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **TIENET_AUG** — `GVS_type = "IEN"` (`en_solver` not applicable; left at the library default)

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "Single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — 2-D SNR x rho grid (5 SNR x 6 rho): `SNR in {0.2, 0.5, 1.0, 2.0, 5.0}`, `rho in {0.30, 0.50, 0.70, 0.90, 0.95, 0.99}` (`sd_x` derived from rho).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_02_mc_sim_scattered_blocks/demo_trex_gvs_02_mc_sim_scattered_blocks.R      # all parts
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_02_mc_sim_scattered_blocks/demo_trex_gvs_02_mc_sim_scattered_blocks.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` (`gvs_scattered_grouped_2d.txt` + `.csv`); the script sources the shared `../trex_gvs_dgps.R`, `../trex_gvs_sim_common.R`, and `../../simulation_utils.R`.

## Interpretation

- Because HAC-based grouping clusters columns by **correlation**, not by column adjacency, scattering group members is expected to leave FDP/TPP patterns very close to Demo 01's contiguous-block results.
- Best read side-by-side with Demo 01: similar results confirm that group discovery depends on correlation structure rather than spatial layout; a noticeable divergence would point to a layout-sensitive bug or clustering limitation.
- As in Demo 01, expect TIENET_AUG to remain conservative (very low FDP and TPP) on this all-active-groups design.
