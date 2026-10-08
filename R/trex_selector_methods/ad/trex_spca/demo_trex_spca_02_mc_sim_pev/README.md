# Demo 02: Cumulative Explained Variance of T-Rex SPCA (R)

## Purpose

While demo 01 asks whether the false discovery rate (FDR) of the loading support is controlled, this demo asks two
other questions. **What proportion of variance does the selected support actually explain, and where does a
non-sparse method's apparent advantage come from?**

It runs all four sweeps of **Fig. 3** of [[1]](#references) for the Terminating Random Experiments (T-Rex) selector for
sparse principal component analysis (SPCA), with the cumulative explained variance (EV) under the paper's Definition 1,
and reports the conventional total-variance normalization alongside it. Under one denominator ordinary principal
component analysis (PCA) looks best, under the other it looks worst, and both statements come from the very same
explained variance. A fifth part maps the realized FDR and true positive rate (TPR) of the first principal component
(PC1) over the full grid of target FDR and signal-to-noise ratio (SNR).

This is the R counterpart of the C++ demo `demo_trex_spca_02_mc_sim_pev.cpp` and the Python demo of the same name, with
the identical configuration. The [C++ README](../../../../../cpp/trex_selector_methods/ad/trex_spca/demo_trex_spca_02_mc_sim_pev/README.md)
holds the full discussion of the metric and the reference results.

---

## Data Generation Parameters (`generate_sparse_factor_model()`)

The sparse $M$-factor model of [demo 01](../demo_trex_spca_01_mc_sim/README.md), with $p_1 = 10$ as in Fig. 3:

```math
\boldsymbol{X} = \boldsymbol{Z}\boldsymbol{V}^\top + \boldsymbol{E}
```

- $n = 50$, $p = 100$, $M = 3$ true factors, $\sigma_m \in \lbrace5, 3, 1\rbrace$.
- $p_1 = 10$ active loadings of value $0.9$ per factor, drawn without replacement from a **shared pool** of
   `overlap_pool_size` $= 30$ candidates, so factor supports partially coincide.
- All methods see the same **center-only** $\boldsymbol{X}$ (covariance-PCA footing).

---

## Control Parameters

```text
tFDR = 0.10           # Target FDR (swept in Part 4)
lambda_2 = -1         # < 0 selects the ridge penalty by k-fold CV (lambda2_method = "CV_1SE_CCD")
MC = 200              # Monte Carlo repetitions per grid point
base_seed = 42        # Per-grid-point seed band, derived from the swept value
```

Only the **data** is seeded deterministically, while each trial's dummies come from hardware entropy (selector
`seed = -1L`), which is what makes the Monte Carlo estimate valid. The seeds keep the values of the C++ and Python demos,
while R's generator draws other numbers from them, so the three languages agree within Monte Carlo noise only.

---

## Methods Compared

Identical to [demo 01](../demo_trex_spca_01_mc_sim/README.md#methods-compared): **OrdPCA**, **OraclePCA**, and the four
T-Rex SPCA variants crossing `TENET`/`TENET_AUG` with `ActiveSet`/`Thresholded`. The proportion of explained variance
(PEV) is reported as `PEV` (total-variance share), `PEVsig` (Definition 1, the EV over the Signal + Mixed EV of the
data) and `PEVsigmix` (the method's own Signal + Mixed part over the same denominator, for OrdPCA the paper's
"Ordinary PCA (Sig + Mix)" reference curve).

---

## The Sweeps

Each of Parts 1 to 4 varies **one** axis and holds the rest at their defaults (#PCs $= 3$, SNR $= 0$ dB, $p_1 = 10$,
$\mathrm{tFDR} = 0.10$):

| Part | Axis | Grid |
| --- | --- | --- |
| 1 | Number of **extracted** PCs | $\lbrace1, 3, 5, 10, 20, 30, 40, 49\rbrace$ |
| 2 | SNR (dB) | $\lbrace-10, -7, -5, -3, 0, 3, 5, 7, 10\rbrace$ |
| 3 | Number of true active loadings $p_1$ | $\lbrace1, 5, 10, 15, 20, 25, 30\rbrace$ |
| 4 | Target FDR | $\lbrace0.01, 0.05, 0.10, 0.15, \ldots, 0.50\rbrace$ |
| 5 | Realized FDR / TPR **heatmaps** over the full target-FDR $\times$ SNR grid | $11 \times 9$ points, PC1 support, two selector responses |

- **Part 1 varies only how many components are extracted.** The data keeps its $M = 3$ true factors, so every
   component beyond the third can contribute nothing but Null EV. The grid stops at 49, because the centered
   $\boldsymbol{X}$ has rank at most $n - 1 = 49$.
- **Only the T-Rex variants respond to Part 4.** The two PCA baselines never see `tFDR` and are flat references.
- **Part 5** runs the per-PC selector for grouped variable selection (GVS) with the elastic net, T-Rex+GVS(EN), twice
   on the same data, once on the **plug-in response** $y = \boldsymbol{X} v_1$ (what T-Rex SPCA regresses on) and once
   on the **oracle response** $y = z_1$, the true factor scores. The pairing isolates the mechanism behind the FDR
   overshoot at loose targets.

---

## Output Files

Written to `simulation_results/data/`:

- `demo_trex_spca_02_mc_sim_pev_pcs.txt` / `.csv`
- `demo_trex_spca_02_mc_sim_pev_snr.txt` / `.csv`
- `demo_trex_spca_02_mc_sim_pev_loadings.txt` / `.csv`
- `demo_trex_spca_02_mc_sim_pev_tfdr.txt` / `.csv`
- `demo_trex_spca_02_mc_sim_pev_fdr_heatmap.csv` — Part 5 (`response,metric,tfdr,snr_db,value`)

Each comma-separated values (CSV) file records FDR, TPR, `PEV`, `PEVsig` and `PEVsigmix`, while the `snr`, `loadings`
and `tfdr` tables show only the cumulative-PEV rows.

Figures in Portable Network Graphics (PNG) and Portable Document Format (PDF) go to `simulation_results/plots/`,
produced by `./generate_plots.sh`, which combines the four sweeps into one panel figure per normalization:

- `demo_trex_spca_02_mc_sim_pev_signal.png` — cumulative PEV under Definition 1.
- `demo_trex_spca_02_mc_sim_pev_total.png` — cumulative PEV as a share of the total variance.
- `demo_trex_spca_02_mc_sim_pev_fdr_heatmap.png` — Part 5, realized FDR / TPR heatmaps for both responses.

---

## Running the Demo

```bash
Rscript R/trex_selector_methods/ad/trex_spca/demo_trex_spca_02_mc_sim_pev/demo_trex_spca_02_mc_sim_pev.R      # 6 workers
Rscript R/trex_selector_methods/ad/trex_spca/demo_trex_spca_02_mc_sim_pev/demo_trex_spca_02_mc_sim_pev.R 8    # 8 workers
TREX_SPCA_PARTS=5 Rscript .../demo_trex_spca_02_mc_sim_pev.R      # only the Part 5 heatmaps
TREX_SPCA_NUM_MC=10 Rscript .../demo_trex_spca_02_mc_sim_pev.R    # pilot-run MC override
./generate_plots.sh                                               # figures from the saved CSVs
```

The Monte Carlo trials run on spawned R worker processes (a PSOCK cluster via `mc_map()` in
[../../../mc_runtime.R](../../../mc_runtime.R)), each with single-threaded solvers. The demo is **resumable**. Every
finished cell (one method at one grid point, or one heatmap cell) is appended to
`simulation_results/cache/demo_trex_spca_02_mc_sim_pev.cells.tsv`, and a restarted run reads those cells back instead
of recomputing them. The cache is removed once all results are written, and a change of the configuration or of the
run tag `TREX_RUN_TAG` sets an old cache aside. Part 1 dominates the runtime, because extracting 49 components means 49
T-Rex selections per trial.

---

## Simulation Results

No R results are committed yet. They follow with the first full run of the R suite. The
[C++ README](../../../../../cpp/trex_selector_methods/ad/trex_spca/demo_trex_spca_02_mc_sim_pev/README.md) shows the
reference results, and the R run is expected to agree with them within Monte Carlo noise.

---

## References

1. Machkour, J., Breloy, A., Muma, M., Palomar, D. P., & Pascal, F., "Sparse PCA with False Discovery Rate Controlled
   Variable Selection.", IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2024,
   pp. 9716–9720, IEEE.
   [DOI-Link](https://doi.org/10.1109/ICASSP48485.2024.10448237)
