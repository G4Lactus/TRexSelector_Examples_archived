# Demo 08: Scalability Benchmark

## Purpose

Characterize **Terminating Random Experiments (T-Rex) selector runtime and memory usage** as $n$ and $p$ scale to
large sizes, and demonstrate that the selector goes **high-dimensional and
out-of-core** — into design volumes that exceed physical random-access memory
(RAM). Everything runs on the **fully memory-mapped pipeline**. Every Monte Carlo
(MC) trial generates its $\mathbf{X}/\mathbf{y}$ into memory-mapped backing files
(`SyntheticDataMapped` over `memmap::MemoryMappedMatrix`, $\mathbf{X}$ never
materialized in RAM), and `use_memory_mapping = true` additionally memory-maps the
internal dummy matrices $\mathbf{D}$ and serializes solver warm-start state to
disk. The benchmark reports:

- **Runtime scaling** — wall-clock per phase (data generation, `select()`,
  $\text{per-}T\text{-iteration}$) as a function of $n$ and $p$
- **Memory scaling** — resident set size (RSS) after selection, process peak RSS,
  physical footprint, and the $\mathbf{X}/\mathbf{y}$ backing sizes on disk
- **Solver comparison** — TLARS vs. TOMP vs. TAFS at scale
- **Selection quality at scale** — false discovery rate (FDR), true positive rate (TPR), average $L$, average $T$ as
  a sanity check

---

## The Grid

A 2-D sweep over $n$ and $p$, `num_MC` sequential trials per $(n, p, \text{solver})$, default 10 and settable
through the command-line interface (CLI):

- $n \in \lbrace1{,}000, 1{,}500, 2{,}000, 5{,}000\rbrace$
- $p \in \lbrace1{,}000, 5{,}000, 10{,}000, 20{,}000, 50{,}000, 100{,}000, 150{,}000, 200{,}000\rbrace$

$\Rightarrow$ 32 grid points × 3 solvers. The logical design volume is
$[\mathbf{X}\,|\,\mathbf{D}] = n \times (1 + L) \times p \times 8$ bytes with
$L = 10$:

| $n$ \\ $p$ | 50k | 100k | 150k | 200k |
| --- | --- | --- | --- | --- |
| 2,000 | 8.2 GiB | 16.4 GiB | 24.6 GiB | 32.8 GiB |
| 5,000 | 20.5 GiB | 41.0 GiB | 61.5 GiB | **82.0 GiB** |

The **top-right corner** ($n = 5{,}000$, $p = 200{,}000$) is the **out-of-core
barrier point**. Its $[\mathbf{X}\,|\,\mathbf{D}] \approx 82$ GiB exceeds a 64 GiB
machine's physical RAM, so it is feasible only because the pipeline is fully
memory-mapped. Barrier points (design $>$ physical RAM) are flagged with `*` in
the output table.

The signal-to-noise ratio (SNR) is fixed at $1.0$ — it is deliberately **not** a scaling axis here (sweeping
it only burns time). That sweep belongs to Demos 02–07.

A tiny smoke grid ($n \in \lbrace150, 200\rbrace$, $p \in \lbrace500, 1{,}000\rbrace$) is available via
the `smoke` argument for end-to-end pipeline verification.

### Memory behaviour

The whole pipeline is memory-mapped, so $\mathbf{X}/\mathbf{y}$ live on disk,
the dummy matrices $\mathbf{D}$ are memory-mapped (mmap), and solver warm-start state is
serialized between T-loop iterations. The **resident** memory therefore stays
bounded even at the largest points. Measured max RSS is a few GiB while the
$\mathbf{X}$ backing at the corner is ~8 GiB, most of it reclaimable page cache of
the touched mmap regions.

The **physical footprint** is a different story, and is why the `Foot_MiB` column
exists. The per-experiment coefficient-path matrices
($(p + L \cdot p) \times \text{path-steps}$, allocated per experiment per T-iteration
and never touched again) accumulate as cold anonymous heap. At the largest scales this reaches tens
of GiB — macOS's memory compressor swaps it out, so the process does **not** run out of memory (OOM),
but on systems without generous swap/compression this allocation pattern is the
scalability limit. It is the main upstream optimization opportunity this benchmark
surfaces (free or stream the beta paths after each T-iteration's $\Phi$
computation).

For disk sizing, the $\mathbf{X}$ backing file is $n \cdot p \cdot 8$ bytes per trial
and the mmap dummy matrix $\mathbf{D}$ grows up to `max_dummy_multiplier` × that,
so the corner needs hundreds of GB of free disk for the dummy backing.

---

## Data Generation Parameters

The data-generating process (DGP) mirrors **Demo 03** (variable support), with Normal predictors and Normal noise.

- **True support cardinality**: fixed at $s = 10$, while the **indices are redrawn
  per trial** via `draw_support` (uniform without replacement,
  `mt19937(seed + 500000)`)
- **True coefficients**: fixed $\beta_j = 1$ (`rnd_coef = false`)
- **SNR**: fixed at $1.0$
- Every trial writes its own mmap-backed $\mathbf{X}/\mathbf{y}$
  (`SyntheticDataMapped`), fed to the selector through `Eigen::Map` views and
  removed at end of scope by a Resource Acquisition Is Initialization (RAII)
  `MmapFileGuard` (declared before the data
  object so the mmap handles close before the files are deleted — exception-safe).

Trials run **sequentially** (no Open Multi-Processing (OpenMP) parallel MC loop) so
wall-clock and RSS measurements are undistorted, while `omp_set_num_threads(6)` still applies to the K
random experiments inside `select()` (inner OpenMP) — the parallelism this
benchmark relies on.

---

## Control Parameters

Every control option is set **explicitly and identically** in the C++ and Python
demos (nothing left to a struct/helper default — that is what made the two
diverge):

```text
K                        = 20        # Random experiments per T-loop iteration
max_dummy_multiplier     = 10        # Max dummies L = 10p (drives the D volume)
use_max_T_stop           = true      # Cap T <= ceil(n/2)
opt_threshold            = 0.75      # Voting/optimization threshold (C++ 0.75, Python default 0.02 — pinned)
dummy_distribution       = Normal    # Dummy predictors drawn from N(0,1)
lloop_strategy           = STANDARD  # Fresh i.i.d. dummy matrix per L-loop iteration
tloop_stagnation_stop    = true      # Early exit when R_mat stagnates
tloop_max_stagnant_steps = 5         # Stagnation window
parallel_rnd_experiments = false     # Sequential experiments (Eigen-parallel inside); one D working set at a time
use_memory_mapping       = true      # All points: X/y mmap + D mmap + solver serialization
tFDR                     = 0.1       # Target FDR control level
```

`parallel_rnd_experiments = false` is essential here, because the alternative (parallel
across the K experiments) holds up to *K* dummy-augmented working sets in RAM at
once, which defeats the out-of-core goal. With it off, the K experiments run
sequentially and each one's linear algebra uses the inner OpenMP threads
(`omp_set_num_threads(6)`).

---

## Solvers Compared

Three base solvers spanning the least-angle regression (LARS) path and greedy families (Demo 05 configuration):

- **TLARS** — equiangular LARS path that terminates on its own
- **TOMP** — greedy orthogonal matching pursuit
- **TAFS** — greedy adaptive forward selection, `rho_afs = 0.3` (TLARS/TOMP take `0`)

---

## Measurements

Per $(n, p)$ × solver, averaged/aggregated over the cached trials:

| Column | Meaning |
| --- | --- |
| `DesignGiB` | Logical design volume $[\mathbf{X}\,\|\,\mathbf{D}]$ (GiB) |
| `SelMed_s` / `SelMin_s` / `SelMax_s` | Median / min / max `select()` wall-clock (s) |
| `DGP_s` | Data-generation wall-clock (s) |
| `s/Titer` | `select()` time per T-loop iteration (s) |
| `PeakRSS_MiB` | Process peak RSS (running maximum) |
| `RSS_MiB` | Post-`select()` resident set size |
| `Foot_MiB` | Physical footprint (anonymous heap incl. compressed/swapped pages) |
| `X+y_MiB` | $\mathbf{X}/\mathbf{y}$ mmap backing-file sizes on disk |
| `FDR` / `TPR` / `Avg_L` / `Avg_T` | Selection-quality sanity metrics |

RSS is read via `mach task_info(MACH_TASK_BASIC_INFO)` (macOS) / `/proc/self/statm`
(Linux), the physical footprint via `task_info(TASK_VM_INFO).phys_footprint`
(macOS) / `VmRSS + VmSwap` (Linux), and peak RSS via `getrusage`. The
`SelMed`/`SelMin`/`SelMax` triple (rather than a single mean) is deliberate,
because per-trial `select()` times span a wide range at scale.

---

## Output Files

Written to `simulation_results/data/` with stem `demo_trex_08_scalability_results`,
**re-derived from the trial cache and re-saved after each completed grid point**
(a crash or OOM kill on a later, larger point preserves everything finished so far):

- **`demo_trex_08_scalability_results_trials.csv`** — the **raw per-trial cache**
  (source of truth): one row per $(n, p, \text{solver}, \text{mc})$ with its seed
  and every measured quantity, appended as each trial finishes.
- **`demo_trex_08_scalability_results.txt`** — per-solver aligned table, one row
  per $(n, p)$ grid point across the 14 metric columns, with barrier points marked `*`.
- **`demo_trex_08_scalability_results.csv`** — tidy long summary, header column
  order `n,p,is_barrier,solver,num_mc,snr,metric,value` (`num_mc` is the actual
  cached trial count per cell).

---

## Expandability

The benchmark is designed to grow incrementally rather than be re-run from
scratch. Every trial is keyed by $(n, p, \text{solver}, \text{mc})$ in the raw
cache, and a run computes **only the trials it is missing**:

- **num_MC top-up** — run at `num_MC = 2` now, then `num_MC = 10` later computes
  just the 8 new trials per cell and re-aggregates, and the first 2 are reused.
- **Grid growth** — append new $n$ / $p$ values to `N_VALUES` / `P_VALUES` later,
  so existing cells are reused and only the new ones are computed.
- **Trial-level resume** — a crash or OOM mid-grid-point keeps every trial that
  already finished (each is flushed to the cache immediately).

Per-trial seeds are a stable function of $(n, p, \text{mc})$ —
`24 + n·100003 + p·13 + mc` (mod $2^{31}{-}1$), independent of grid position and
**identical to the Python demo's `_trial_seed`** — so reuse and reproducibility
survive reordering or extending the grid. (The two languages seed identically, but the
*draws* still differ because `std::mt19937` and NumPy's pseudo-random number generator
(PRNG) are different generators.)

---

## Results Visualization

Demo 08's comma-separated values (CSV) file uses the wider `n,p,is_barrier,solver,num_mc,snr,metric,value` schema
rather than the suite's tidy `solver,metric,snr,value`, so the suite plotter
handles it in a dedicated `scalability` mode
([`../trex_plt_utils.py`](../trex_plt_utils.py) `scalability <csv> …`). The wide
CSV is read by `read_scalability()` and rendered by `plot_scalability_dashboard()`
— a 2×2 dashboard (runtime vs $p$, memory vs $p$, the out-of-core money plot with
the RAM line and the barrier point ringed, and FDR control at scale). Render it
with `./generate_plots.sh`.

---

## Status of the Recorded Data

The full 32-point × 3-solver grid is a **long-running benchmark**. A **first pass
at `num_MC = 2`** is run to get early numbers over the whole grid, and `num_MC` is then
topped up to 10 later as a pure cache extension (only the 8 missing trials per cell
are computed). The committed `demo_trex_08_scalability_smoke.{txt,csv}` is the
end-to-end smoke recording that exercises the pipeline and table/CSV format.

---

## Building and Running

Release build is recommended (the large points are heavy):

```bash
cd TRexSelector_Examples/cpp
cmake --build build/release --target demo_trex_08_mc_sim_scalability

bin=build/release/bin/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability/demo_trex_08_mc_sim_scalability
"$bin"           # full grid, num_MC = 10 (long-running)
"$bin" 2         # full grid, num_MC = 2  (first pass; top up to 10 later)
"$bin" smoke     # tiny end-to-end smoke check (num_MC = 2)
"$bin" smoke 3   # smoke grid, num_MC = 3
```

The optional integer argument sets `num_MC`, and re-running with a larger value tops
up the trial cache rather than recomputing. Grid points run smallest-first and the
summaries are re-derived from the cache after each point. The benchmark sweeps its
own memory-mapped orphans (`demo08_scalability_*`, `trex_memmap_*`) at startup and
on exit / SIGINT / SIGTERM (SIGKILL cannot be trapped, so the next startup sweep
reclaims it).
