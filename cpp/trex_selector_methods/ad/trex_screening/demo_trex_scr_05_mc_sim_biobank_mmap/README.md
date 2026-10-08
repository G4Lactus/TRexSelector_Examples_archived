# Demo 05: Biobank Screen-TRex (Algorithm 1) — Memory-Mapped

## Purpose

The demo repeats the **Biobank Screen-TRex workflow** ("Algorithm 1") of demo 04 with
 **memory-mapped dummy matrices** (`use_memory_mapping = true`), the storage mode intended for
 biobank-scale problems where holding the dummy matrices in random-access memory (RAM) is infeasible.
 The algorithm is unchanged. Per phenotype it tries **Screen-TRex Ordinary** first. If the estimated false discovery
 rate (FDR) falls outside the acceptable window $\lbrack 0.05, 0.15\rbrack$ it retries with
 **Screen-TRex Bootstrap-CI**, and if that also lands outside the window it falls back to the
 **classical Terminating Random Experiments (T-Rex) selector** at a target FDR of $0.10$.
 Demo 04 is the in-memory reference for every number reported here — the point of this demo is to show
 that moving the dummy matrices to disk changes the memory footprint and not the statistics.
 Screening returns a *candidate set*, and the FDR and the true positive rate (TPR) are evaluated on the individual
 selected variables
 (see [What is actually measured](../README.md#what-is-actually-measured-in-these-demos)).

---

## Data Generation Parameters (`make_iid_dgp`)

We consider the linear model:

```math
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon},
\qquad \boldsymbol{\epsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma_{\varepsilon}^2 \boldsymbol{I}_n)
```

- $\boldsymbol{y} \in \mathbb{R}^n$ is the response vector.
- $\boldsymbol{X} \in \mathbb{R}^{n \times p}$ is the design matrix.
- $\boldsymbol{\beta} \in \mathbb{R}^p$ is the coefficient vector, with $s$ nonzero entries.
- $\boldsymbol{\epsilon}$ is the noise vector, independent and identically distributed (i.i.d.) standard normal.
- $\sigma_{\varepsilon}^2$ is the noise variance, calibrated to achieve a target linear signal-to-noise ratio (SNR).
- $n = 300$, $p = 1{,}000$, with $s = 10$ in Part 1 and $s = 5$ in Part 2 (high-dimensional, $p > n$).

The design matrix has **no correlation structure**:

```math
X_{ij} \sim \mathcal{N}(0,1) \quad \text{i.i.d.}
```

- The active support is drawn uniformly at random *per Monte Carlo (MC) trial*, so the results are not tied to
   one support pattern.
- All active coefficients are $\beta_j = 1$, and the support-selection and coefficient random number generators
   (RNGs) are offset from the trial seed so they stay independent of the design and noise draws.
- In Part 2 all $q = 5$ phenotypes of a trial share **one design matrix** $\boldsymbol{X}$ (a fixed
   `X_seed` per trial), while each phenotype receives its **own random support** and its own noise draw —
   the biobank setting, where a single genotype matrix is screened against many traits.

---

## Control Parameters

```text
K = 20                       # Random experiments per T-loop iteration
lower_bound_FDR = 0.05       # Lower edge of the acceptable estimated-FDR window
upper_bound_FDR = 0.15       # Upper edge of the acceptable estimated-FDR window
target_FDR_trex = 0.10       # Target FDR of the classical T-Rex fallback
R_boot = 1000                # Bootstrap replicates (Bootstrap-CI rule only)
ci_grid_step = 0.001         # Bootstrap-CI threshold grid granularity
solver = TLARS               # T-Rex solver backend
use_memory_mapping = true    # Dummy matrices stored on disk (Boost memory-mapped files)
MC = 200                     # Monte Carlo repetitions per grid point
```

The only difference from demo 04 is `use_memory_mapping = true`. The window $\lbrack 0.05, 0.15\rbrack$ is the
 routing criterion. A screening result is *accepted* only if the procedure's own estimated FDR lands inside it.

---

## Methods Compared

The three routing outcomes of Algorithm 1 [[1]](#references), all using `ScreenTRexMethod::TREX`:

- **Screen-TRex Ordinary** — the majority-vote rule ($\Phi_j > 0.5$), tried first for every phenotype.
- **Screen-TRex Bootstrap-CI** — the bootstrap confidence-band rule (`R_boot = 1000` replicates), tried
   only when the Ordinary estimate falls outside the window.
- **T-Rex fallback: TLARS** — the classical T-Rex selector at `target_FDR_trex = 0.10`, used when neither
   screening rule produces an acceptable estimate.

All reported FDR and TPR values are **conditional on routing**. They average only over the phenotypes that
 were actually handled by the respective method. The additional **Usage %** row is unconditional and states
 what fraction of all screenings each method took. A metric row for a method with $0\%$ usage is therefore
 an average over an empty set, not a performance statement.

---

## The Two Parts

- **Part 1 — single phenotype.** One response per trial, $s = 10$ active variables, SNR sweep over
   $\lbrace0.1, 0.5, 1, 2, 5\rbrace$, 200 MC trials per grid point.
- **Part 2 — $q = 5$ phenotypes sharing one $\boldsymbol{X}$.** $s = 5$ active variables per phenotype, the
   same SNR grid and 200 MC trials, hence $5 \times 200 = 1{,}000$ phenotype screenings per SNR level.

---

## Output Files

Written to `simulation_results/data/`:

- `scr_biobank_mmap_snr_n300_p1000_s10.txt` / `.csv` — Part 1: Usage %, FDR, TPR, and estimated FDR
   per method and SNR level.
- `scr_biobank_mmap_multi_n300_p1000_q5_s5.txt` / `.csv` — Part 2: the same metrics for the
   multi-phenotype study.

Figures in Portable Network Graphics (PNG) and Portable Document Format (PDF) go to `simulation_results/plots/`,
produced by `./generate_plots.sh`.

---

## Running the Demo

```bash
./build/release/bin/trex_selector_methods/ad/trex_screening/demo_trex_scr_05_mc_sim_biobank_mmap/demo_trex_scr_05_mc_sim_biobank_mmap
./generate_plots.sh   # render the figures below from the saved CSVs
```

---

## Simulation Results

### Part 1 — Single Phenotype ($s = 10$)

- **Routing behaves exactly as in the in-memory run.** At $\mathrm{SNR} = 0.1$ all phenotypes ($100.0\%$)
   go to the classical T-Rex fallback, both screening estimates being far above the window ($0.460$
   Ordinary, $0.382$ Bootstrap-CI). At $\mathrm{SNR} = 1$ the split is $58.0\%$ Ordinary / $35.0\%$
   Bootstrap-CI / $7.0\%$ fallback, and at $\mathrm{SNR} = 5$ Bootstrap-CI handles $89.0\%$ of all
   screenings with the fallback never invoked.
- **Conditional error control on the screening routes is good once signal is present.** For
   $\mathrm{SNR} \ge 1$ the realized FDR is $0.052\text{–}0.070$ for Ordinary and $0.057 \to 0.050$ for
   Bootstrap-CI, with TPR reaching $1.000$ for both at $\mathrm{SNR} = 5$. As in demo 04,
   $\mathrm{SNR} = 0.5$ is the weak point, where the accepted screenings realize an FDR of $0.208$ (Ordinary) and
   $0.177$ (Bootstrap-CI).
- **The fallback's TPR of $0.000$ at $\mathrm{SNR} = 2$ and $5$ is not a failure.** Its usage is $0.0\%$ at
   those points — nothing is routed to it, so the conditional average is over an empty set and prints as
   zero. Where it is used it is the conservative option. At $\mathrm{SNR} = 1$ its realized FDR is $0.012$
   against a target of $0.10$, with a TPR of $0.486$.
- **The estimates mostly err on the safe side.** Ordinary's estimate settles at $0.095\text{–}0.118$ against a
   realized $0.052\text{–}0.070$ for $\mathrm{SNR} \ge 1$, and Bootstrap-CI's at $0.043\text{–}0.085$ against a
   realized $0.050\text{–}0.057$, below the realized value only at $\mathrm{SNR} = 1$.

TPR (conditional on routing), FDR (solid = realized, dashed = the procedure's own estimate, with the
$\lbrack 0.05, 0.15\rbrack$ routing window shaded), and per-method Usage % vs. SNR.

![Biobank Screen-TRex (memory-mapped), single phenotype: TPR/FDR/Usage vs SNR](simulation_results/plots/scr_biobank_mmap_snr_n300_p1000_s10.png)

### Part 2 — Five Phenotypes Sharing One Design ($q = 5$, $s = 5$)

- **The fallback dominates, and more signal does not change that.** The classical T-Rex selector handles
   $84.7\%\text{–}100.0\%$ of all screenings ($100.0\%$ at $\mathrm{SNR} = 0.1$, $86.6\%$ at $\mathrm{SNR} = 5$),
   with Ordinary at $1.7\%\text{–}3.3\%$ and Bootstrap-CI at $7.6\%\text{–}13.2\%$ for $\mathrm{SNR} \ge 0.5$.
- **The rarely-used screening routes are poorly calibrated.** For $\mathrm{SNR} \ge 0.5$ the realized FDR of
   the accepted screening results is $0.306\text{–}0.343$ (Ordinary) and $0.305\text{–}0.315$ (Bootstrap-CI), roughly
   twice the upper edge of the routing window, against estimates of only $0.179\text{–}0.194$ and
   $0.144\text{–}0.160$. With just $s = 5$ active variables per phenotype a single false selection moves the false
   discovery proportion (FDP) a long way, and the estimate does not track it — even though TPR is $1.000$ for both
   rules from $\mathrm{SNR} = 1$ onward.
- **The fallback remains reliable on the bulk of the traffic.** Its realized FDR is $0.012\text{–}0.019$ against
   a target of $0.10$, with TPR rising $0.013 \to 1.000$ across the sweep.
- **Memory mapping changes the footprint, not the statistics.** The Part 2 Ordinary route here reads
   Usage $0.0/3.3/3.1/1.7/2.7\%$, FDR $0.000/0.343/0.306/0.325/0.319$ and estimate
   $0.524/0.194/0.179/0.180/0.181$, against demo 04's $0.0/3.3/2.9/3.4/2.2\%$,
   $0.000/0.338/0.307/0.306/0.316$ and $0.531/0.192/0.180/0.180/0.181$ — agreement within Monte Carlo
   noise. Exact equality is not expected, because the selector is constructed with `seed = -1`, which the library
   resolves via `std::random_device`, so the two runs use different random experiments.

TPR (conditional on routing), FDR (solid = realized, dashed = the procedure's own estimate, with the
$\lbrack 0.05, 0.15\rbrack$ routing window shaded), and per-method Usage % vs. SNR.

![Biobank Screen-TRex (memory-mapped), five phenotypes: TPR/FDR/Usage vs SNR](simulation_results/plots/scr_biobank_mmap_multi_n300_p1000_q5_s5.png)

---

## References

1. Machkour, J., Muma, M., & Palomar, D. P., "False Discovery Rate Control for Fast Screening of
   Large-Scale Genomics Biobanks.", IEEE Statistical Signal Processing Workshop (SSP), 2023,
    pp. 666–670, IEEE.
    [DOI-Link](https://doi.org/10.1109/SSP53291.2023.10207957)
