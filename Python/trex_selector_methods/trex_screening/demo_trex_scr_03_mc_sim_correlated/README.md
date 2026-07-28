# Demo 03: Screen-TRex under Correlated Designs (Python)

## Purpose

The demo measures what **predictor correlation does to Screen-TRex**, and whether the **dependency-aware (DA)**
 variants repair the damage.
 Three correlation structures are put in front of the screener — AR(1), equicorrelated, and
 block-equicorrelated — and in each case the plain Ordinary and Bootstrap-CI rules are compared against
 their DA counterparts, which pre-group correlated variables before the votes are counted.
 Demo 01 established the i.i.d. baseline; everything reported here is the cost of leaving it.
 Screening returns a *candidate set*, and FDR/TPR are evaluated on the individual selected variables
 (see [What is actually measured](../README.md#what-is-actually-measured-in-these-demos)).

This is the Python port of the C++ demo
[demo_trex_scr_03_mc_sim_correlated](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_03_mc_sim_correlated/README.md),
using `trex_selector_neo`'s `TRexScreeningSelector`.

---

## Data Generation Parameters (`dgp_ar1`, `dgp_equi`, `dgp_block_equi`)

We consider the linear model:

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon},
\qquad \boldsymbol{\epsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma_{\varepsilon}^2 \boldsymbol{I}_n)
$$

- $\boldsymbol{y} \in \mathbb{R}^n$ is the response vector.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ is the design matrix.
- $\boldsymbol{\beta} \in \mathbb{R}^p$ is the coefficient vector, with $s$ nonzero entries.
- $\boldsymbol{\epsilon}$ is the noise vector, i.i.d. standard normal.
- $\sigma_{\varepsilon}^2$ is the noise variance, calibrated to achieve a target linear signal-to-noise ratio (SNR).
- $n = 300$, $p = 1000$, $s = p_1 = 10$ (high-dimensional, $p > n$).

Three correlation structures are used. **AR(1)**, a decaying-correlation chain:

$$
x_j = \rho\, x_{j-1} + \sqrt{1 - \rho^2}\, w_j
\qquad \Longrightarrow \qquad
\mathrm{corr}(x_j, x_{j'}) = \rho^{|j - j'|}
$$

**Equicorrelated**, one shared latent factor $z$ driving every column, so *all* pairs correlate equally:

$$
x_j = \sqrt{\rho}\, z + \sqrt{1 - \rho}\, w_j
\qquad \Longrightarrow \qquad
\mathrm{corr}(x_j, x_{j'}) = \rho \quad \forall\, j \neq j'
$$

**Block-equicorrelated**, the same construction but with one latent factor $z_k$ per block: columns
 correlate at $\rho$ *within* a block and are independent across blocks. Five contiguous blocks of 200
 columns each are used.

- All columns are **standardised to zero mean and unit variance** after generation, so $\rho$ is the
   realized column correlation rather than a nominal parameter.
- $\boldsymbol{\beta}$ has $p_1 = 10$ **evenly spaced** active entries, all of equal magnitude
   ($\beta_j = 1$ for the AR(1) parts, $\beta_j = 3$ for the equicorrelated and block-equicorrelated
   parts, matching the C++ source); the design and the noise are drawn afresh in every Monte Carlo trial.

![Design-matrix correlation — AR(1), rho = 0.5](../dgp_corr_plots/scr_ar1.png)

![Design-matrix correlation — equicorrelated, rho = 0.4](../dgp_corr_plots/scr_equi.png)

![Design-matrix correlation — block-equicorrelated, 5 blocks, rho = 0.9](../dgp_corr_plots/scr_block_equi.png)

Empirical Pearson correlation of exemplary designs from the three DGPs ($n=300$, $p=1000$, seed 2026), in
**natural variable order — no hierarchical clustering**. The AR(1) panel is cropped to the first 150
columns so the narrow decaying band around the diagonal stays visible at print size; the equicorrelated
panel is uniformly correlated everywhere (one shared factor); the block panel shows the five independent
$200 \times 200$ blocks. One figure per DGP lives in [`../dgp_corr_plots/`](../dgp_corr_plots/);
regenerate the whole suite with [`../generate_corr_plots.sh`](../generate_corr_plots.sh) (shared plotter
[`../../corr_plt_utils.py`](../../corr_plt_utils.py)).

---

## Control Parameters

```text
K = 20                       # Random experiments per T-loop iteration
rho_thr_DA = 0.02            # Correlation threshold for dependency-aware grouping
n_blocks = 5                 # Blocks assumed by the DA-BLOCK-EQUI variant (Part 5)
R_boot = 1000                # Bootstrap replicates (Bootstrap-CI rule only)
ci_grid_step = 0.001         # Bootstrap-CI threshold grid granularity
solver = TLARS               # T-Rex solver backend
MC = 200                     # Monte Carlo repetitions per grid point
```

Note that Screen-TRex has **no target-FDR parameter**: unlike the classical T-Rex selector, screening
 thresholds the voting statistic instead of calibrating to a user-specified level.

---

## Methods Compared

Every part compares four methods — the two thresholding rules of Demo 01, plain and dependency-aware
 [[1]](#references):

- **Screen-TRex Ordinary: TLARS** — majority vote $\{ j : \Phi_j > 0.5 \}$, no correlation handling
   (`ScreenTRexMethod.TREX`).
- **Screen-TRex Bootstrap-CI: TLARS** — bootstrap confidence band around the estimated FDR curve
   (`R_boot = 1000` replicates), no correlation handling.
- **Screen-TRex Ordinary+DA / Bootstrap-CI+DA: TLARS** — the same two rules on top of a
   **dependency-aware** method (`ScreenTRexMethod.TREX_DA_AR1`, `.TREX_DA_EQUI`, or
   `.TREX_DA_BLOCK_EQUI`), which groups variables whose correlation exceeds `rho_thr_DA = 0.02` *before*
   the votes are aggregated, so a cluster of near-duplicate predictors cannot inflate the selection count
   on its own.

The DA variant is always matched to the design that generated the data — this demo therefore reports the
 **best case for DA**, not a misspecification study. All methods also report an **estimated FDR**
 alongside the realized FDR/TPR; under correlation the gap between the two becomes the main story.

---

## The Five Parts

| Part | Design | Fixed | Sweep | DA variant |
| :--- | :----- | :---- | :---- | :--------- |
| 1 | AR(1) | $\rho = 0.5$ | SNR $\in \{0.01, 0.1, 0.2, 0.5, 0.6, 1, 2, 5\}$ | DA-AR1 |
| 2 | Equicorrelated | $\rho = 0.4$ | same SNR grid | DA-EQUI |
| 3 | AR(1) | SNR $= 1$ | $\rho \in \{0.1, \dots, 0.9\}$ | DA-AR1 |
| 4 | Equicorrelated | SNR $= 1$ | $\rho \in \{0.1, \dots, 0.9\}$ | DA-EQUI |
| 5 | Block-equicorrelated (5 blocks) | SNR $= 1$ | $\rho \in \{0.1, \dots, 0.9\}$ | DA-BLOCK-EQUI |

Parts 1–2 ask *how much signal is needed at a fixed correlation level*; Parts 3–5 ask *how much
 correlation the screener survives at a fixed signal level*. 200 MC trials per grid point throughout,
 with per-cell trial seeds following the C++ source (`base_seed = 42 + idx * 1000`).

---

## Output Files

Written to `simulation_results/data/`, each as `.txt` (human-readable table) and `.csv` (tidy long):

- `scr_ar1_snr_n300_p1000_rho0.50` — Part 1: AR(1), SNR sweep.
- `scr_equi_snr_n300_p1000_rho0.40` — Part 2: equicorrelated, SNR sweep.
- `scr_ar1_rho_n300_p1000_snr1.00` — Part 3: AR(1), $\rho$ sweep.
- `scr_equi_rho_n300_p1000_snr1.00` — Part 4: equicorrelated, $\rho$ sweep.
- `scr_block_equi_rho_n300_p1000_blocks5_snr1.00` — Part 5: block-equicorrelated, $\rho$ sweep.

Each file holds FDR, TPR, and estimated FDR per method and grid point. The file names and CSV schemas
match the C++ demo, and the suite plotter [`../trex_scr_plt_utils.py`](../trex_scr_plt_utils.py) is the
same module the C++ suite uses. Figures (PNG + PDF) go to `simulation_results/plots/`, produced by
`./generate_plots.sh`.

---

## Running the Demo

```bash
python demo_trex_scr_03_mc_sim_correlated.py   # full run (~26 min with 6 workers)
./generate_plots.sh                            # render the figures below from the saved CSVs
```

Monte Carlo trials are parallelized with `ProcessPoolExecutor`; the worker count defaults to
`min(6, usable CPUs)` and is overridden with `SCR_NUM_WORKERS` (or the suite-wide `TREX_NUM_WORKERS`).
This is the longest of the four non-biobank demos: five parts $\times$ four methods $\times$ 8–9 grid
points $\times$ 200 trials.

---

## Simulation Results

The numbers and figures below are this Python run's own results (200 MC trials per grid point). They agree
with the [C++ reference](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_03_mc_sim_correlated/README.md#simulation-results)
up to Monte-Carlo noise — the two ports draw from independent RNG streams (NumPy here, the C++ engine
there).

### Part 1 — AR(1) Design, SNR Sweep ($\rho = 0.5$)

- **AR(1) correlation is survivable.** At $\rho = 0.5$ the Ordinary rule behaves much as in Demo 01: FDR
   falls $0.43 \to 0.11$ across the sweep while TPR climbs $0.002 \to 1.00$, recovering all ten active
   variables at $\mathrm{SNR} = 5$. Local, decaying dependence does not break the screener.
- **DA lowers the error rate at a modest power cost.** Ordinary+DA-AR1 roughly *halves* the plain
   Ordinary FDR for $\mathrm{SNR} \ge 1$ ($0.053$ vs. $0.107$ at $\mathrm{SNR} = 1$; $0.042$ vs. $0.093$
   at $\mathrm{SNR} = 2$) while giving up only a few points of TPR ($0.756$ vs. $0.824$).
- **Bootstrap-CI+DA-AR1 is extremely conservative.** Its FDR never exceeds $0.021$ once
   $\mathrm{SNR} \ge 0.5$, but its TPR tops out at $0.430$ at $\mathrm{SNR} = 5$ — it buys near-zero
   false discoveries by leaving more than half the true signal on the table.
- **The internal estimate is still well-behaved here.** For $\mathrm{SNR} \ge 0.5$ it brackets the
   realized FDR closely and often exceeds it (Ordinary+DA-AR1: $0.131$ estimated vs. $0.053$ realized at
   $\mathrm{SNR} = 1$). This is the last part where that holds.

TPR (left) and FDR (right) vs. SNR, one line per method; on the FDR panel the solid line is
the realized FDR and the dashed line the procedure's own estimated FDR.

![Screen-TRex on an AR(1) design: TPR/FDR vs SNR](simulation_results/plots/scr_ar1_snr_n300_p1000_rho0.50.png)

### Part 2 — Equicorrelated Design, SNR Sweep ($\rho = 0.4$)

- **Equicorrelation breaks screening outright.** Even at a moderate $\rho = 0.4$, the Ordinary rule's FDR
   sits between $0.80$ and $0.96$ for every $\mathrm{SNR} \ge 0.1$ — four out of five selections are
   false — and it is still $0.802$ at $\mathrm{SNR} = 5$, where TPR has reached $0.994$. More signal
   recovers the true variables but does *not* clean up the candidate set.
- **DA-EQUI does not rescue the Ordinary rule.** Ordinary+DA-EQUI tracks plain Ordinary almost exactly
   ($0.846$ vs. $0.855$ at $\mathrm{SNR} = 1$). Matching the DA variant to the true design structure is
   not enough when every column shares one latent factor.
- **Bootstrap-CI+DA-EQUI helps partially, and pays for it.** It is the one method that bends the curve
   ($0.680$ vs. $0.845$ at $\mathrm{SNR} = 1$), but $0.60$–$0.68$ FDR is still unusable, and its TPR
   trails the others throughout ($0.819$ vs. $0.994$ at $\mathrm{SNR} = 5$).
- **A non-monotone artifact at the bottom of the sweep.** The realized Ordinary FDR is *lower* at
   $\mathrm{SNR} = 0.01$ ($0.567$) than at $\mathrm{SNR} = 0.1$ ($0.959$). This is not a sweet spot: at
   $\mathrm{SNR} = 0.01$ the rule selects almost nothing (TPR $= 0.002$), so there are barely any
   discoveries left to be false. Low-SNR FDR numbers must be read alongside the TPR column.
- **The FDR estimate has already come apart.** For $\mathrm{SNR} \ge 0.1$ it reports $0.020$–$0.113$
   while the realized FDR is $0.80$–$0.96$ — off by an order of magnitude, in the optimistic direction.

TPR (left) and FDR (right) vs. SNR, one line per method; solid = realized FDR, dashed =
estimated FDR.

![Screen-TRex on an equicorrelated design: TPR/FDR vs SNR](simulation_results/plots/scr_equi_snr_n300_p1000_rho0.40.png)

### Part 3 — AR(1) Design, $\rho$ Sweep (SNR $= 1$)

- **Degradation is gradual and only bites at high $\rho$.** The Ordinary rule holds FDR near
   $0.067$–$0.090$ up to $\rho = 0.4$, then deteriorates to $0.124$, $0.148$, $0.218$, $0.330$, and
   finally $0.511$ at $\rho = 0.9$. TPR is remarkably stable ($\approx 0.82$) until $\rho = 0.6$ and drops
   to $0.637$ only at $\rho = 0.9$.
- **DA-AR1 buys roughly one step of the $\rho$ grid.** Ordinary+DA-AR1 lowers FDR at every $\rho$
   ($0.417$ vs. $0.511$ at $\rho = 0.9$; $0.234$ vs. $0.330$ at $\rho = 0.8$) at a cost of a few to
   twelve points of TPR. It is a real improvement, but it delays the failure rather than preventing it.
- **Bootstrap-CI+DA-AR1 controls the FDR by declining to select.** Its FDR stays within $0.007$–$0.010$
   across the *entire* $\rho$ range — nominally the best result in this demo — but its TPR falls from
   $0.177$ at $\rho = 0.1$ to $0.093$ at $\rho = 0.6$ and $0.022$ at $\rho = 0.9$. At high correlation it
   returns essentially an empty set. This is a bias/variance trade, not a free fix, and a near-zero FDR
   is only meaningful when read next to its TPR.
- **The estimate degrades in step with the method.** At $\rho = 0.9$ Ordinary reports $0.078$ against a
   realized $0.511$ — the estimate is least trustworthy exactly where the user most needs a warning.

TPR (left) and FDR (right) vs. $\rho$, one line per method; solid = realized FDR, dashed = estimated FDR.

![Screen-TRex on an AR(1) design: TPR/FDR vs rho](simulation_results/plots/scr_ar1_rho_n300_p1000_snr1.00.png)

### Part 4 — Equicorrelated Design, $\rho$ Sweep (SNR $= 1$)

- **Failure from the very first grid point.** Ordinary FDR starts at $0.661$ for $\rho = 0.1$ and rises
   monotonically to $0.965$ at $\rho = 0.9$, while TPR collapses $0.904 \to 0.050$. Unlike AR(1), there
   is no usable region of the sweep.
- **No DA variant rescues it.** Ordinary+DA-EQUI ends at $0.964$ and Bootstrap-CI at $0.965$ —
   indistinguishable from plain Ordinary. Bootstrap-CI+DA-EQUI is again the only one that differs at low
   $\rho$ ($0.007$ FDR at $\rho = 0.1$), and again only because it selects almost nothing there
   (TPR $= 0.010$); by $\rho = 0.5$ it has converged to the others at $0.761$.
- **This is the demo's central caution: the FDR estimate is wildly optimistic exactly where the method
   fails.** Across the whole sweep the estimate reports $0.027$–$0.077$ while the realized FDR is
   $0.66$–$0.97$ — a gap of a factor of twenty or more, and it is *smallest* at the mildest correlation.
   Under equicorrelation the procedure's self-assessment cannot be trusted at all. Contrast this with
   Demo 01, where the same estimate sat *above* the realized FDR and erred on the safe side.

TPR (left) and FDR (right) vs. $\rho$, one line per method; solid = realized FDR, dashed = estimated FDR.

![Screen-TRex on an equicorrelated design: TPR/FDR vs rho](simulation_results/plots/scr_equi_rho_n300_p1000_snr1.00.png)

### Part 5 — Block-Equicorrelated Design, $\rho$ Sweep (SNR $= 1$, 5 blocks)

- **Block structure sits between the two extremes — but closer to the bad one.** Ordinary FDR grows
   $0.193 \to 0.915$ over the sweep. It is far better than the fully equicorrelated case at low $\rho$
   ($0.193$ vs. $0.661$ at $\rho = 0.1$), because independence *across* blocks survives, but by
   $\rho = 0.9$ the two designs are equally hopeless.
- **DA-BLOCK-EQUI helps meaningfully at low correlation only.** Ordinary+DA-BLOCK-EQUI cuts FDR from
   $0.193$ to $0.046$ at $\rho = 0.1$ and from $0.408$ to $0.182$ at $\rho = 0.2$ — the clearest DA win
   anywhere in this demo — but the advantage has eroded by $\rho = 0.5$ ($0.653$ vs. $0.718$) and is gone
   at $\rho = 0.9$ ($0.911$ vs. $0.915$). The TPR cost is real throughout ($0.594$ vs. $0.851$ at
   $\rho = 0.1$).
- **Bootstrap-CI+DA-BLOCK-EQUI is degenerate here.** Its FDR never exceeds $0.118$, but its TPR never
   exceeds $0.047$ and is exactly $0.000$ at $\rho = 0.2$. The method returns empty or near-empty sets
   across the sweep; its FDR column carries no usable information.
- **The estimate is again too optimistic, though less catastrophically.** Ordinary reports
   $0.037$–$0.097$ against a realized $0.19$–$0.91$. As in Part 4, it gives no signal that anything has
   gone wrong.

TPR (left) and FDR (right) vs. $\rho$, one line per method; solid = realized FDR, dashed = estimated FDR.

![Screen-TRex on a block-equicorrelated design: TPR/FDR vs rho](simulation_results/plots/scr_block_equi_rho_n300_p1000_blocks5_snr1.00.png)

---

## References

1. Machkour, J., Muma, M., & Palomar, D. P., "False Discovery Rate Control for Fast Screening of
   Large-Scale Genomics Biobanks.", IEEE Statistical Signal Processing Workshop (SSP), 2023,
    pp. 666–670, IEEE.
    [DOI-Link](https://doi.org/10.1109/SSP53291.2023.10207957)

---

**Last updated**: 2026-07-27
