# Demo 01: Single-Run T-Rex Selector (Python)

## Purpose

Demonstrate the classical T-Rex selector as a single-run example — correct API
usage and the resulting console/file output — in two contrasting scenarios:

1. **High-dimensional** ($p > n$): moderate sample, many predictors,
2. **Low-dimensional** ($n > p$): large sample, moderate predictors.

This is not a Monte Carlo study; the results come from a single random seed and
should not be read as population averages (see Demos 02–07 for that). Python port
of
[cpp/trex_selector_methods/trex/demo_trex_01_single_run](../../../../cpp/trex_selector_methods/trex/demo_trex_01_single_run/),
using `trex_selector_neo`'s `TRexSelector`.

---

## Data Generation Parameters

Both scenarios use the shared Gaussian generator `dgp_gauss_snr`:

$$
\boldsymbol{y} = \boldsymbol{X}\boldsymbol{\beta} + \boldsymbol{\epsilon},
\qquad \boldsymbol{\epsilon} \sim \mathcal{N}(\boldsymbol{0}, \sigma^2 \boldsymbol{I}_n),
$$

with an i.i.d. standard-normal design $\boldsymbol{X}$ (independent columns), so
the design carries **no correlation structure** — the empirical column-correlation
matrix is essentially the identity:

![Design-matrix correlation of the shared i.i.d. Gaussian DGP (natural order, no clustering)](../dgp_corr_plots/trex_gauss.png)

- Fixed true support (0-based): $\mathcal{S}^\ast = \{4, 27, 42, 43, 128, 149\}$ ($s = 6$)
- Fixed coefficients $\beta_j = 1$ (`rnd_coef=False`; as in C++, the `rnd_coef`
  path selects a fixed hardcoded vector $\{-0.4, -0.25, -0.8, 1.1, 2.5, -1.2\}$
  that mimics random coefficients — it is not a random draw)
- SNR $= 1.0$, tFDR $= 0.1$
- **Scenario A — High-dimensional** ($p > n$): $n = 150$, $p = 300$ (run first, as in C++)
- **Scenario B — Low-dimensional** ($n > p$): $n = 5000$, $p = 1000$
- Seeds mirror C++ structurally (data seed `58`, selector seed `42`); the RNG
  streams differ across languages, so the drawn data — and hence the exact
  selection — are not bit-identical to the C++ output.

---

## Control Parameters

```text
solver_name           = TLARS      # Base solver (least-angle regression)
K                     = 20         # Random experiments per T-loop iteration
max_dummy_multiplier  = 10         # Max L = 10 x p
use_max_T_stop        = True       # Cap T <= ceil(n/2)
dummy_distribution    = Normal     # Dummy predictors drawn from N(0,1)
lloop_strategy        = STANDARD   # Fresh i.i.d. dummies each L-loop iteration
tloop_stagnation_stop = False      # No stagnation early-exit
tFDR                  = 0.1        # Target FDR control level
```

---

## Output Files

Dual console + file output (mirrors the C++ `print_dual` pattern). One `.txt`
per scenario in `simulation_results/data/`, stem `d01_trex_basic_n{n}_p{p}.txt`:

- `d01_trex_basic_n150_p300.txt` (high-dimensional)
- `d01_trex_basic_n5000_p1000.txt` (low-dimensional)

Each file records the data/calibration header, then the C++-style selection
block: selected indices, FDP/TPP, the sparsely-printed `phi_prime` and `phi_mat`
(nonzero entries with indices), `fdp_hat_mat`, `r_mat`, and the voting grid.

---

## Recorded Results

From the committed recordings (single seed — not an average):

| Scenario | Calibration | Selected | TP | FP | FDP | TPP |
| --- | --- | --- | --- | --- | --- | --- |
| A — high-dim ($n{=}150$, $p{=}300$) | $T_{\text{stop}} = 41$, $L = 10$ | 6 | 5 | 1 | 0.1667 | 0.8333 |
| B — low-dim ($n{=}5000$, $p{=}1000$) | $T_{\text{stop}} = 2$, $L = 2$ | 6 | 6 | 0 | 0.0000 | 1.0000 |

- The **low-dimensional** case is easy: all six actives recovered with no false
  positive after only $T_{\text{stop}} = 2$ T-loop steps.
- The **high-dimensional** case recovers five of the six actives and admits a
  single false positive (index 64), for FDP $= 1/6 \approx 0.167$ on this seed.
  Being a one-shot run, FDP can exceed the tFDR $= 0.1$ target on an individual
  draw; the Monte Carlo demos show the *average* FDR stays controlled.
- Numbers differ from the C++ recording (which recovers three of six on its seed)
  because NumPy's PRNG produces a different draw — the qualitative behaviour is
  the same.

---

## Running

```bash
python Python/trex_selector_methods/trex/demo_trex_01_single_run/demo_trex_01_single_run.py
```

Single-seed demo — no `multiprocessing`.

---

**Last updated**: 2026-07-25
