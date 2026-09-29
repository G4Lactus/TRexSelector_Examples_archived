# Demo 04: T-Rex+GVS on Negative-Correlation Traps

## Purpose

Test whether T-Rex+GVS can recover an **active group with negative within-group correlation** (sign-flipped halves) while **excluding two spatially separated inactive trap groups** that share the same sign-flipped pattern. This stresses the grouping logic's ability to distinguish signal from confounding correlation geometry, not just correlation strength.

## Data-generating process

`dgp_neg_traps_snr` — each correlated group loads half its members `+Z` and half `-Z`, inducing negative within-group correlation:

- **Group 1** (cols 1–100): ACTIVE, `beta = +3 / -3`, `s = 100`.
- **Trap 1** (cols 101–200): inactive, same sign-flipped structure.
- **Noise** (cols 201–300): white noise.
- **Trap 2** (cols 301–360): second, smaller inactive sign-flipped trap.
- **Noise** (cols 361–500): white noise.

Within-group correlation follows `rho = 1/(1 + sd_x^2)`. Dimensions: `n = 200`, `p = 500`. Default `sd_x = sqrt(0.01)` (`rho ~ 0.99`).

## Methods compared

- **TENET** — `GVS_type = "EN"`, `en_solver = "TENET"`
- **TENET_AUG** — `GVS_type = "EN"`, `en_solver = "TENET_AUG"`
- **TIENET_AUG** — `GVS_type = "IEN"` (`en_solver` not applicable; left at the library default)

Shared controls: `K = 20`, `tFDR = 0.1`, `corr_max = 0.98`, `hc_dist = "Single"`, `num_MC = 200`, `seed = 2026`.

## Parts

- **Part 1** — 2-D SNR x rho grid (5 SNR x 6 rho): `SNR in {0.2, 0.5, 1.0, 2.0, 5.0}`, `rho in {0.30, 0.50, 0.70, 0.90, 0.95, 0.99}` (`sd_x` derived from rho).

## Running

```bash
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_04_mc_sim_neg_traps/demo_trex_gvs_04_mc_sim_neg_traps.R      # all parts
Rscript R/trex_selector_methods/ad/trex_gvs/demo_trex_gvs_04_mc_sim_neg_traps/demo_trex_gvs_04_mc_sim_neg_traps.R 8    # 8 worker cores
```

Results are written to `simulation_results/data/` (`gvs_Neg-Traps_2d.txt` + `.csv`); the script sources the shared `../trex_gvs_dgps.R`, `../trex_gvs_sim_common.R`, and `../../simulation_utils.R`.

## Interpretation

- TENET and TENET_AUG should track each other closely with FDP near the `tFDR = 0.1` target and TPP reaching near-complete recovery at moderate SNR — the sign-flipped active group is still identified as one group via the (unsigned) correlation-based HAC step, and both traps should be excluded.
- TIENET_AUG is expected to stay conservative here (low FDP, slowly rising TPP), similar to the uniformly-positive-correlation scenarios; the sign flips alone do not change IEN behavior much.
- The two spatially separated traps are a stronger test of false-group exclusion than a single trap — watch whether FDP creeps up with SNR.
