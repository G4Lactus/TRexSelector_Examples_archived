# Demo 08: Scalability Benchmark (Python)

## Purpose

Characterize Terminating Random Experiments (T-Rex) selector **runtime and memory
usage** as $n$ and $p$ scale to large sizes, and demonstrate that the selector
goes **high-dimensional and out-of-core** — into design volumes that exceed
physical random-access memory (RAM). Everything runs on the **fully
memory-mapped pipeline**. Every Monte Carlo trial generates $\mathbf{X}$ directly
into a memory-mapped (mmap) backing file (chunked column generation over
`MemoryMappedMatrix`, so $\mathbf{X}$ is never materialized in RAM), and
`use_memory_mapping = True` additionally memory-maps the internal dummy matrices
$\mathbf{D}$ and serializes solver warm-start state to disk. Python port of
[cpp/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability](../../../../../cpp/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability/),
using `trex_selector_neo`.

The benchmark reports, per $(n, p)$ × solver grid point:

- **Runtime scaling** — wall-clock per phase (data generation, `select()`,
  $\text{per-}T\text{-iteration}$, total)
- **Memory scaling** — post-`select()` resident set size (RSS), physical
  footprint, process peak RSS, and the $\mathbf{X}/\mathbf{y}$ backing size on disk
- **Selection quality at scale** — false discovery rate (FDR), true positive rate
  (TPR), average $L$, average $T$ as a sanity check

---

## Data Generation Parameters

The data-generating process (DGP) mirrors **Demo 03** (variable support) — an
independent and identically distributed (i.i.d.) standard-normal design with
independent columns, so the column-correlation is essentially the identity:

![Design-matrix correlation of the shared i.i.d. Gaussian DGP (shown at n=300 and p=1,000 in natural order, no clustering)](../dgp_corr_plots/trex_gauss.png)

- **True support**: cardinality $s = 10$, and the **indices are redrawn per
  trial** via `make_support_random`
- **True coefficients**: fixed $\beta_j = 1$ (`rnd_coef = False`)
- **Signal-to-noise ratio (SNR)**: fixed at $1.0$ — SNR is deliberately **not** a
  scaling axis here (that sweep belongs to Demos 02–07), so this demo scales $n$
  and $p$ only
- **Seeding**: `base_seed = 24 + 1000 * (grid index)`, unique per grid point and
  shared across solvers so every solver sees identical data

Trials run **sequentially** (no `ProcessPoolExecutor`) so wall-clock and RSS
readings are undistorted, while the K random experiments inside `select()` still
run in parallel via the selector's inner Open Multi-Processing (OpenMP),
**pinned to 6 threads**
(`set_num_threads(6)`, `_NUM_THREADS`) to match the C++ demo's
`omp_set_num_threads(6)`. Left at the library default (all logical cores) the
$\text{large-}p$ points run 6 vs. 12 threads' worth of concurrent dummy matrices near
the RAM limit, which distorts the timing — so both languages use the same thread
budget for a fair comparison.

---

## The Grid

A 2-D sweep over $n$ and $p$, `num_MC = 10` sequential trials per $(n, p, \text{solver})$:

- $n \in \lbrace1{,}000, 1{,}500, 2{,}000, 5{,}000\rbrace$
- $p \in \lbrace1{,}000, 5{,}000, 10{,}000, 20{,}000, 50{,}000, 100{,}000, 150{,}000, 200{,}000\rbrace$

$\Rightarrow$ 32 grid points × 3 solvers. The logical design volume is
$[\mathbf{X}\,|\,\mathbf{D}] = n \times (1 + L) \times p \times 8$ bytes with
$L = 10$, which is what the memory-mapped pipeline offloads to disk:

| $n$ \\ $p$ | 50k | 100k | 150k | 200k |
| --- | --- | --- | --- | --- |
| 2,000 | 8.2 GiB | 16.4 GiB | 24.6 GiB | 32.8 GiB |
| 5,000 | 20.5 GiB | 41.0 GiB | 61.5 GiB | **82.0 GiB** |

The **top-right corner** ($n = 5{,}000$, $p = 200{,}000$) is the **out-of-core
barrier point**. Its $[\mathbf{X}\,|\,\mathbf{D}] \approx 82$ GiB exceeds the
machine's 64 GiB of physical RAM, so it is feasible only because the pipeline is
fully memory-mapped. Barrier points (design $>$ physical RAM) are flagged with
`*` in the output table.

---

## Control Parameters

Every control option is set **explicitly and identically** to the C++ demo
(nothing left to the Python helper's defaults — which differ from the C++ struct
and are what made the two demos diverge):

```text
K                        = 20        # Random experiments per T-loop iteration
max_dummy_multiplier     = 10        # Max dummies L = 10p (drives the D volume)
use_max_T_stop           = True      # Cap T <= ceil(n/2)
opt_threshold            = 0.75      # Voting/optimization threshold (Python default is 0.02 — pinned to C++ 0.75)
dummy_distribution       = Normal    # Dummy predictors drawn from N(0,1)
lloop_strategy           = STANDARD  # Fresh i.i.d. dummy matrix per L-loop iteration
tloop_stagnation_stop    = True      # Early exit when R_mat stagnates
tloop_max_stagnant_steps = 5         # Stagnation window
parallel_rnd_experiments = False     # Sequential experiments (Python default is True — pinned to False)
use_memory_mapping       = True      # All points: X mmap + D mmap + solver serialization
tFDR                     = 0.1       # Target FDR control level
```

Two of these — `opt_threshold` and `parallel_rnd_experiments` — have **different
defaults in the Python bindings than in the C++ struct** (0.02 vs 0.75, and True
vs false). Leaving them implicit is exactly what made the first Python pass
incomparable to C++, so all options are pinned here (matching demos 02/03/07).
`parallel_rnd_experiments = False` also keeps only one experiment working set in
RAM at a time — essential for the out-of-core corner.

---

## Solvers Compared

Three base solvers spanning the least-angle regression (LARS) path and greedy
families (Demo 05 configuration):

- **TLARS** — equiangular LARS path that terminates on its own
- **TOMP** — greedy orthogonal matching pursuit
- **TAFS** — greedy adaptive forward selection, `rho_afs = 0.3` (TLARS/TOMP take `0`)

---

## Measurements

Per $(n, p)$ × solver, averaged/aggregated over the `num_MC` sequential trials:

| Column | Meaning |
| --- | --- |
| `DesignGiB` | Logical design volume $[\mathbf{X}\,\|\,\mathbf{D}]$ (GiB) |
| `SelMed_s` / `SelMin_s` / `SelMax_s` | Median / min / max `select()` wall-clock (s) |
| `DGP_s` | Data-generation wall-clock (s) |
| `s/Titer` | `select()` time per T-loop iteration (s) |
| `PeakRSS_MiB` | Process peak RSS (running maximum) |
| `RSS_MiB` | Post-`select()` resident set size |
| `Foot_MiB` | Physical footprint (anonymous heap incl. compressed/swapped pages) |
| `X+y_MiB` | $\mathbf{X}$ backing-file size + in-RAM $\mathbf{y}$ |
| `FDR` / `TPR` / `Avg_L` / `Avg_T` | Selection-quality sanity metrics |

RSS is read via `ps -o rss=` (macOS) / `/proc/self/statm` (Linux), the physical
footprint via `task_info(TASK_VM_INFO).phys_footprint` through `ctypes` (macOS) /
`VmRSS + VmSwap` (Linux), and peak RSS via `resource.getrusage`.

The `SelMed`/`SelMin`/`SelMax` triple (rather than a single mean) is deliberate.
Per-trial `select()` times span a wide range at the larger grid points, so the
median is the robust summary and the min–max band shows the spread.

---

## Output Files

Written to `simulation_results/data/` with stem `demo_trex_08_scalability_results`,
**re-derived from the trial cache and re-saved after each completed grid point**
(a crash or out-of-memory (OOM) kill on a later, larger point preserves
everything finished so far):

- **`demo_trex_08_scalability_results_trials.csv`** — the **raw per-trial cache**
  (source of truth): one row per $(n, p, \text{solver}, \text{mc})$ with its seed
  and every measured quantity. Appended to as each trial finishes.
- **`demo_trex_08_scalability_results.txt`** — per-solver aligned table, one row
  per $(n, p)$ grid point across the 14 metric columns, with barrier points
  marked `*`.
- **`demo_trex_08_scalability_results.csv`** — tidy long summary, header column
  order `n,p,is_barrier,solver,num_mc,snr,metric,value` (`num_mc` is the actual
  cached trial count per cell). Both summaries are derived from the cache.

---

## Expandability

The benchmark is designed to grow incrementally rather than be re-run from
scratch. Every trial is keyed by $(n, p, \text{solver}, \text{mc})$ in the raw
cache, and a run computes **only the trials it is missing**:

- **num_MC top-up** — run at `num_MC = 2` now, then `num_MC = 10` later computes
  just the 8 new trials per cell and re-aggregates, reusing the first 2.
- **Grid growth** — append new $n$ / $p$ values to `_N_VALUES` / `_P_VALUES`
  later. Existing cells are reused, and only the new ones are computed.
- **Trial-level resume** — a crash or OOM mid-grid-point keeps every trial that
  already finished (each is flushed to the cache immediately).

Per-trial seeds are a stable function of $(n, p, \text{mc})$ —
`24 + n·100003 + p·13 + mc` (mod $2^{31}{-}1$), independent of grid position and
identical to the C++ demo — so reuse and reproducibility survive reordering or
extending the grid. (The two languages seed identically, but the *draws* still
differ because NumPy's pseudo-random number generator (PRNG) and C++'s
`std::mt19937` are different generators.)

---

## Results Visualization

Demo 08's comma-separated values (CSV) file uses the wider
`n,p,is_barrier,solver,…` schema rather than the suite's tidy
`solver,metric,snr,value`, so the C++ suite plotter handles it in a dedicated
`scalability` mode
([`../../../../../cpp/trex_selector_methods/ad/trex_core/trex_plt_utils.py`](../../../../../cpp/trex_selector_methods/ad/trex_core/trex_plt_utils.py)
`scalability <csv> …`). The wide CSV is read by `read_scalability()` and rendered
by `plot_scalability_dashboard()` as a 2×2 dashboard:

1. **Runtime vs $p$** (colour = $n$ and linestyle = solver, with a median + min–max band).
2. **Memory vs $p$** — peak RSS per $n$ against the $\mathbf{X}/\mathbf{y}$ backing
   size, with the physical-RAM line far above the resident set.
3. **Out-of-core money plot** — resident peak RSS vs design volume (log–log), the
   64 GiB RAM line drawn in, and the out-of-core barrier point ringed, which sits
   at $\approx 82$ GiB of design while staying resident.
4. **FDR control at scale** — per solver, against the target FDR (tFDR) line.

Since the Python demo emits the identical CSV schema, the same plotter renders the
Python results, but a per-demo `generate_plots.sh` is not committed on the Python
side.

---

## Status of the Recorded Data

The full 32-point × 3-solver grid is a **long-running benchmark** (the larger
$(n, p)$ points and the out-of-core corner take minutes each, sequentially). A
**first pass at `num_MC = 2`** is run to get early numbers over the whole grid,
and `num_MC` is then topped up to 10 later as a pure cache extension (only the 8
missing trials per cell are computed). The committed
`demo_trex_08_scalability_smoke.{txt,csv}` is the end-to-end **smoke recording**
(tiny $n \in \lbrace150, 200\rbrace$, $p \in \lbrace500, 1{,}000\rbrace$) that
exercises the full memory-mapped pipeline and the table/CSV format — a pipeline
check, not a scientific result.

---

## Running

```bash
demo=Python/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability/demo_trex_08_mc_sim_scalability.py

python "$demo"            # full grid, num_MC = 10 (long-running)
python "$demo" 2         # full grid, num_MC = 2  (first pass; top up to 10 later)
python "$demo" smoke     # tiny end-to-end smoke check (num_MC = 2)
python "$demo" smoke 3   # smoke grid, num_MC = 3
```

The optional integer argument sets `num_MC`, and re-running with a larger value
tops up the trial cache rather than recomputing. Grid points run smallest-first and the
summaries are re-derived from the cache after each point, so an OOM/crash on a
later point preserves everything finished so far. The benchmark installs cleanup
handlers (`atexit` + `SIGINT`/`SIGTERM`) that sweep any orphaned memory-mapped
backing files on exit.
