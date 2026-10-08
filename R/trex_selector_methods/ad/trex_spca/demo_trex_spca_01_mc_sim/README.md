# Demo 01: T-Rex SPCA Monte Carlo Simulation (R)

## Purpose

Compare **T-Rex SPCA**, the Terminating Random Experiments (T-Rex) selector for sparse principal component analysis
(SPCA), in 4 solver/mode combinations against **ordinary PCA** and **oracle-thresholded PCA**, the two principal
component analysis (PCA) baselines, on a synthetic sparse 3-factor model. The demo evaluates the loading-support
recovery of the first principal component (PC1) by the false discovery rate (FDR) and the true positive rate (TPR), and
the cumulative proportion of explained variance (PEV). This is the R counterpart of the C++ demo
`demo_trex_spca_01_mc_sim.cpp` and runs the identical simulation configuration.

---

## Data Generation Parameters (`generate_sparse_factor_model()`)

```math
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E}.
```

- $\boldsymbol{Z}$ ($n \times M$): factor columns $\mathcal{N}(0, \sigma_m^2)$ with $\sigma_m \in \lbrace5, 3, 1\rbrace$
  for $m=1,2,3$.
- $\boldsymbol{V}$ ($p \times M$): each column has exactly `p1` nonzero entries (value $0.9$), drawn without replacement
  from a shared candidate pool of size `overlap_pool_size` — factor supports can partially overlap.
- $\boldsymbol{E}$: independent and identically distributed (i.i.d.) Gaussian, scaled to hit a target signal-to-noise
  ratio (SNR) in dB.
- $n=50$, $p=100$, $p_1=5$, $M=3$, `overlap_pool_size=30`.

All methods share one **center-only** $\boldsymbol{X}$ (see `center_columns()` in the suite-level
`../trex_spca_sim_utils.R`: mean subtraction, no column scaling). This puts OrdPCA, OraclePCA, and TRexSPCA on a common
covariance-PCA footing, matching the legacy Comprehensive R Archive Network (CRAN) R reference. The column scales must
**not** be normalized — the factor amplitude signal lives in the column variances, and z-scoring destroys it (see the
History note below).

---

## Control Parameters

```text
tFDR = 0.10
lambda_2 = -1         # < 0 triggers k-fold CV for ridge-penalty selection (lambda2_method = "CV_1SE_CCD")
num_MC = 200          # per SNR point of the grid {-10, -7, -5, -3, 0, 3, 5, 7, 10} dB
base_seed = 42        # seeds the per-trial DATA DGP deterministically:
                      #   trial mc uses base_seed_offset + mc*1000, where
                      #   base_seed_offset = base_seed + round(snr)
```

The selector itself runs with `seed = -1L` (hardware entropy for the dummies) — required for a valid Monte Carlo (MC)
FDR estimate.

---

## Methods Compared

| Method | Sparsity mechanism | `mode` | `en_solver` |
|---|---|---|---|
| **OrdPCA** | None (all $p$ loadings retained) | — | — |
| **OraclePCA** | $\text{Top-}p_1$ by magnitude (true support size known) | — | — |
| **TRexSPCA-EN-Act** | T-Rex+GVS selection | `"ActiveSet"` | `"TENET"` |
| **TRexSPCA-ENaug-Act** | T-Rex+GVS selection | `"ActiveSet"` | `"TENET_AUG"` |
| **TRexSPCA-EN-Thr** | T-Rex+GVS selection | `"Thresholded"` | `"TENET"` |
| **TRexSPCA-ENaug-Thr** | T-Rex+GVS selection | `"Thresholded"` | `"TENET_AUG"` |

The PCA baselines use base-R `svd()`. Each T-Rex variant is invoked as

```r
sel <- TRexSPCASelector$new(
  X + 0, M = M, tFDR = tFDR,
  control = trex_spca_control(
    mode = "ActiveSet",          # or "Thresholded"
    gvs_type = "EN",
    en_solver = "TENET",         # or "TENET_AUG"
    lambda_2 = -1,
    lambda2_method = "CV_1SE_CCD"
  ),
  seed = -1L
)
```

with result fields `$V` (sparse loadings), `$Z` (scores), and `$active_sets` (**1-based**).

> **Binding note**: `trex_spca_control()` does not expose the per-PC selector's internal `scaling_mode`, but this is
> moot — the C++ demo and the Python port also use the library default (unit-L2 column scaling), which matches the
> legacy `lm_dummy.R` centering + L2 normalization.

---

## Metrics

- **FDR / TPR** — evaluated strictly on **PC1's** loading support only. PCs 2–3 are excluded because ordinary PCA's
  orthogonality constraint mixes their supports across the true factors, making a per-component ground-truth
  comparison ambiguous beyond PC1.
- **PEV** — cumulative proportion of explained variance across all $M$ components, via QR decomposition of the
  estimated score matrix.

---

## Output Files

Written to `simulation_results/data/`:

- `demo_trex_spca_01_mc_sim.txt` — human-readable table
  (methods × FDR/TPR/PEV under three normalizations × SNR columns).
- `demo_trex_spca_01_mc_sim.csv` — tidy long format: `method, metric, snr_db, value`
  (metrics `FDR`, `TPR`, `PEV`, `PEVsig`, `PEVsigmix`).
- `demo_trex_spca_01_mc_sim_union_fdr.csv` — Section 2, tidy long format:
  `pc, metric, snr_db, value` (metrics `FDRunion`, `k`).

Figures (via `./generate_plots.sh`, written to `simulation_results/plots/`):

- `demo_trex_spca_01_mc_sim.{png,pdf}` — TPR/FDR vs SNR.
- `demo_trex_spca_01_mc_sim_pev.{png,pdf}` — cumulative PEV under both normalizations
  (total variance vs Definition 1, incl. the OrdPCA Sig+Mix reference curve).
- `demo_trex_spca_01_mc_sim_union_fdr.{png,pdf}` — union-support FDR and support size per PC.

---

## Running the Demo

```bash
Rscript R/trex_selector_methods/ad/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim.R       # default 6 worker cores
Rscript R/trex_selector_methods/ad/trex_spca/demo_trex_spca_01_mc_sim/demo_trex_spca_01_mc_sim.R 8     # 8 worker cores
TREX_SPCA_PARTS=2 Rscript .../demo_trex_spca_01_mc_sim.R      # only Section 2 (union-FDR)
TREX_SPCA_NUM_MC=10 Rscript .../demo_trex_spca_01_mc_sim.R    # pilot-run MC override
```

Dependencies: `TRexSelectorNeo` only (figures additionally need `ggplot2` + `patchwork`). The Monte Carlo trials run
on spawned R worker processes (a PSOCK cluster via `mc_map()` in [../../../mc_runtime.R](../../../mc_runtime.R)), each
with single-threaded solvers.

The demo is **resumable**. Every finished cell (one method at one SNR point, or one SNR point of Section 2) is appended
to `simulation_results/cache/demo_trex_spca_01_mc_sim.cells.tsv`, and a restarted run reads those cells back instead
of recomputing them. The cache is removed once all results are written, and a change of the configuration or of the
run tag `TREX_RUN_TAG` sets an old cache aside.

---

## Simulation Results

No R results are committed yet. They follow with the first full run of the R suite. The
[C++ README](../../../../../cpp/trex_selector_methods/ad/trex_spca/demo_trex_spca_01_mc_sim/README.md) shows the
reference results over the same grid, and the R run is expected to agree with them within Monte Carlo noise. Exact
numbers cannot match, because R's generator draws other data from the same seeds, and the dummies use hardware entropy
(selector `seed = -1L`), so even an R rerun reproduces its table only within Monte Carlo noise ($\approx \pm 0.01$ at
200 trials).

> **History note (z-scoring episode, fixed 2026-07-08)**: an earlier version of this pipeline z-score-standardized
> $\boldsymbol{X}$, which converts the covariance PCA into a correlation PCA and destroys the factor amplitude signal.
> The measured T-Rex FDR inflated to $\approx 0.53$ at $-10$ dB. The C++ validation program
> `validation_trex_spca_06_handrolled_comparison` isolated the cause, and the pipeline is back to center-only.
