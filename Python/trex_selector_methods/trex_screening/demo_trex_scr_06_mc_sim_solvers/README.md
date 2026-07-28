# Demo 06: Solver Backends for Screen-TRex — TLARS, TAFS, TOMP (Python)

## Purpose

The demo asks whether the **choice of T-Rex solver backend** changes screening performance, or whether it
 is merely a speed knob.
 The same i.i.d. Gaussian design as [Demo 01](../demo_trex_scr_01_mc_sim_screen_trex/README.md) is used, so
 the numbers here are directly comparable to the baseline established there.
 Three backends are swept over SNR — **TLARS** (terminating LARS), **TAFS** with $\rho = 0.3$ (terminating
 adaptive forward selection), and **TOMP** (terminating orthogonal matching pursuit) — first under the
 Ordinary thresholding rule alone, then crossed with the Bootstrap-CI rule.
 Screening returns a *candidate set*, and FDR/TPR are evaluated on the individual selected variables
 (see [What is actually measured](../README.md#what-is-actually-measured-in-these-demos)).

This is the Python port of the C++ demo
[demo_trex_scr_06_mc_sim_solvers](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_mc_sim_solvers/README.md),
using `trex_selector_neo`'s `TRexScreeningSelector`.

---

## Data Generation Parameters (`dgp_iid_snr`)

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
- $n = 300$, $p = 1000$, $s = 10$ (high-dimensional, $p > n$).

The design matrix has **no correlation structure**:

$$
X_{ij} \sim \mathcal{N}(0,1) \quad \text{i.i.d.}
$$

- The active support is drawn uniformly at random *per Monte Carlo trial*, so the results are not tied to
   one support pattern.
- All active coefficients are $\beta_j = 1$; the support-selection and coefficient RNGs are offset from the
   trial seed so they stay independent of the design and noise draws.
- This is the **same DGP as Demo 01**, so the two demos are directly comparable — any difference in the
   curves below is attributable to the solver backend, not to the data.

![Design-matrix correlation — i.i.d. Gaussian design (no correlation structure)](../dgp_corr_plots/scr_iid.png)

Empirical Pearson correlation of an exemplary design drawn from this DGP (`dgp_iid_snr`, $n=300$,
$p=1000$, seed 2026), in **natural variable order — no hierarchical clustering**. One figure per DGP lives
in [`../dgp_corr_plots/`](../dgp_corr_plots/); regenerate the whole suite with
[`../generate_corr_plots.sh`](../generate_corr_plots.sh) (shared plotter
[`../../corr_plt_utils.py`](../../corr_plt_utils.py)).

---

## Control Parameters

```text
K = 20                       # Random experiments per T-loop iteration
R_boot = 1000                # Bootstrap replicates (Bootstrap-CI rule only)
ci_grid_step = 0.001         # Bootstrap-CI threshold grid granularity
solver = TLARS / TAFS / TOMP # T-Rex solver backends compared
rho_afs = 0.3                # AFS penalty parameter (TAFS only)
tloop_stagnation_stop = True # Stagnation stop enabled for every backend here
MC = 200                     # Monte Carlo repetitions per grid point
```

Note that Screen-TRex has **no target-FDR parameter**: unlike the classical T-Rex selector, screening
 thresholds the voting statistic instead of calibrating to a user-specified level.

---

## Methods Compared

Three T-Rex solver backends [[1]](#references), all driving `ScreenTRexMethod.TREX`:

- **TLARS** — terminating LARS, the default backend; follows the least-angle regression path until the
   dummy-based stopping criterion fires.
- **TAFS-0.3** — terminating adaptive forward selection with penalty $\rho = 0.3$; a greedy forward path
   with an adaptive penalty on the already-selected directions.
- **TOMP** — terminating orthogonal matching pursuit; the greediest of the three, taking the predictor most
   correlated with the current residual at each step.

Each backend is combined with the Screen-TRex thresholding rules:

- **Ordinary** — selects $\{ j : \Phi_j > 0.5 \}$, a simple majority vote of the random experiments.
- **Bootstrap-CI** — builds a bootstrap confidence band around the estimated FDR curve (`R_boot = 1000`
   replicates) and picks its threshold from that band.

---

## The Two Parts

- **Part 1 — solvers under the Ordinary rule.** Three series (TLARS, TAFS-0.3, TOMP), Ordinary screening
   only, over $\mathrm{SNR} \in \{0.01, 0.1, 0.2, 0.5, 0.6, 1, 2, 5\}$ with 200 MC trials per grid point.
- **Part 2 — solvers crossed with both rules.** The same three backends under both Ordinary and
   Bootstrap-CI thresholding, giving six series over the same SNR grid and MC budget.

Every trial draws a fresh design, support, and noise realization; per-cell trial seeds follow the C++
source (`base_seed = 24 + snr_idx * 1000`).

---

## Output Files

Written to `simulation_results/data/`:

- `scr_solvers_snr_n300_p1000_s10.txt` / `.csv` — Part 1: FDR, TPR, and estimated FDR per solver and SNR
   level.
- `scr_solver_method_snr_n300_p1000_s10.txt` / `.csv` — Part 2: the same metrics for all six
   solver × rule combinations.

The file names and CSV schemas match the C++ demo, and the suite plotter
[`../trex_scr_plt_utils.py`](../trex_scr_plt_utils.py) is the same module the C++ suite uses.

Figures (PNG + PDF) go to `simulation_results/plots/`, produced by `./generate_plots.sh`.

---

## Running the Demo

```bash
python demo_trex_scr_06_mc_sim_solvers.py   # full run (~12 min with 6 workers)
./generate_plots.sh                         # render the figures below from the saved CSVs
```

Monte Carlo trials are parallelized with `ProcessPoolExecutor`; the worker count defaults to
`min(6, usable CPUs)` and is overridden with `SCR_NUM_WORKERS` (or the suite-wide `TREX_NUM_WORKERS`).

---

## Simulation Results

The numbers and figures below are this Python run's own results (200 MC trials per grid point). They agree
with the [C++ reference](../../../../cpp/trex_selector_methods/trex_screening/demo_trex_scr_06_mc_sim_solvers/README.md#simulation-results)
up to Monte-Carlo noise — the two ports draw from independent RNG streams (NumPy here, the C++ engine
there).

### Part 1 — Three Solvers, Ordinary Screening

- **The solver is not a pure speed knob.** At mid SNR the greedier backends are markedly more powerful than
   TLARS at comparable FDR: at $\mathrm{SNR} = 1$ the TPR is $0.938$ (TAFS-0.3) and $0.953$ (TOMP) against
   $0.829$ (TLARS), while the realized FDR is $0.058$, $0.062$, and $0.066$ respectively. The gap is already
   visible at $\mathrm{SNR} = 0.5$–$0.6$ ($0.477$/$0.475$ vs. $0.431$, and $0.618$/$0.627$ vs. $0.541$).
- **All three converge at high SNR.** By $\mathrm{SNR} = 2$ TOMP reaches TPR $1.000$, TAFS-0.3 $0.999$ and
   TLARS $0.994$; at $\mathrm{SNR} = 5$ all three recover every active variable at FDR $0.057$–$0.059$. The
   choice of backend therefore only matters where the signal is neither absent nor obvious.
- **At very low SNR the backends are indistinguishable — and all uncontrolled.** For $\mathrm{SNR} \le 0.2$
   the realized FDR runs $0.214$–$0.470$ for all three with TPR $\le 0.126$, reproducing the Demo 01 finding
   that screening is not FDR-controlling when there is no recoverable signal.
- **The internal FDR estimate behaves the same way for every backend.** It peaks around $0.45$–$0.55$ at
   $\mathrm{SNR} = 0.1$–$0.2$ and then settles at $0.094$ for all three at $\mathrm{SNR} = 5$, sitting
   *above* the realized FDR once signal is present.

TPR (left) and FDR (right) vs. SNR, one line per solver backend; on the FDR panel the solid line
is the realized FDR and the dashed line the procedure's own estimated FDR.

![Solver comparison under Ordinary screening: TPR/FDR vs SNR](simulation_results/plots/scr_solvers_snr_n300_p1000_s10.png)

### Part 2 — Three Solvers × Two Thresholding Rules

- **Under Bootstrap-CI the backend matters far more than under the Ordinary rule.** At $\mathrm{SNR} = 1$
   the bootstrap TPR spans $0.358$ (TLARS), $0.568$ (TAFS-0.3), and $0.899$ (TOMP) — more than half the
   achievable power — against a spread of roughly $0.12$ for the same three backends under the Ordinary
   rule ($0.838$ / $0.937$ / $0.957$).
- **TLARS is by far the most conservative backend.** Its bootstrap FDR stays at or below $0.058$ from
   $\mathrm{SNR} = 0.5$ onwards, but its TPR reaches only $0.895$ at $\mathrm{SNR} = 5$. TOMP under
   Bootstrap-CI, in contrast, behaves almost like Ordinary screening: at $\mathrm{SNR} = 1$ its TPR is
   $0.899$ vs. $0.957$ ordinary, at FDR $0.059$ vs. $0.065$.
- **One mild non-monotonicity at high SNR, reported as observed.** TAFS-0.3 under Bootstrap-CI peaks at
   TPR $0.705$ ($\mathrm{SNR} = 2$) and then *drops* to $0.630$ at $\mathrm{SNR} = 5$. This is a
   non-monotonicity in the bootstrap threshold at high SNR; no mechanism for it has been verified here, so
   it is flagged rather than explained away. The C++ reference reports a weaker version of the same dip for
   TOMP as well; in this run TOMP under Bootstrap-CI rises monotonically instead
   ($0.899 \to 0.940 \to 0.965$), which puts the TOMP dip within Monte Carlo noise rather than establishing
   it as a property of the backend.
- **The low-SNR regime is unaffected by either choice.** For $\mathrm{SNR} \le 0.2$ all six series have
   TPR $\le 0.123$ and FDR between $0.121$ and $0.470$; neither the backend nor the thresholding rule
   rescues screening when there is nothing to find.
- **Practical takeaway.** TLARS is the safe conservative default; TOMP and TAFS-0.3 buy substantial power at
   mid SNR at a comparable error rate. Under Bootstrap-CI the backend choice moves the operating point more
   than the thresholding rule does — pick the solver first, the rule second.

TPR (left) and FDR (right) vs. SNR, one line per solver × rule combination; on the FDR panel the
solid line is the realized FDR and the dashed line the procedure's own estimated FDR.

![Solver × thresholding rule: TPR/FDR vs SNR](simulation_results/plots/scr_solver_method_snr_n300_p1000_s10.png)

---

## References

1. Machkour, J., Muma, M., & Palomar, D. P., "False Discovery Rate Control for Fast Screening of
   Large-Scale Genomics Biobanks.", IEEE Statistical Signal Processing Workshop (SSP), 2023,
    pp. 666–670, IEEE.
    [DOI-Link](https://doi.org/10.1109/SSP53291.2023.10207957)

---

**Last updated**: 2026-07-27
