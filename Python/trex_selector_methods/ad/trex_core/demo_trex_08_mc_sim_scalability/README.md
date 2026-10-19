# Demo 08: Scalability Benchmark (Python)

## Purpose

Characterize Terminating Random Experiments (T-Rex) selector **runtime and memory
usage** as $n$ and $p$ scale to large sizes, and demonstrate that the selector
goes **high-dimensional and out-of-core**, into data volumes that exceed
physical random-access memory (RAM). Everything runs on the **fully
memory-mapped pipeline**. Every Monte Carlo (MC) trial generates $\mathbf{X}$ directly
into a memory-mapped (mmap) backing file (chunked column generation over
`MemoryMappedMatrix`, so $\mathbf{X}$ is never materialized in RAM), and
`use_memory_mapping = True` additionally memory-maps the internal dummy matrices
$\mathbf{D}$ and serializes solver warm-start state to disk. The selector receives a zero-copy view of the
mapped file, so the binding never copies $\mathbf{X}$. Python port of
[cpp/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability](../../../../../cpp/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability/README.md),
using `trex_selector_neo`. The C++ README explains the memory behaviour and the expected cost in detail.

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
- **Seeding**: a stable function of $(n, p, \text{mc})$ shared across solvers (see
  [Expandability](#expandability)), so every solver sees identical data

Trials run **sequentially** (no `ProcessPoolExecutor`) so wall-clock and RSS
readings are undistorted. The $K$ random experiments inside `select()` also run one after another, and the
selector's Open Multi-Processing (OpenMP) threads parallelize the column loops inside each experiment. The
thread count is **pinned to `TREX_NUM_WORKERS` (default 8)** (`set_num_threads(_NUM_THREADS)`) to match the C++
and R demos. Left at the library default (all logical cores) the threads would compete with the page cache the
out-of-core points live on.

---

## The Grid

Two blocks with separate output files and trial caches, `num_MC = 10` sequential trials per
$(n, p, \text{solver})$:

- **Base block** with `max_dummy_multiplier = 10`, the full cross of
  $n \in \lbrace300, 500, 1{,}000, 1{,}500, 2{,}000, 5{,}000\rbrace$ and
  $p \in \lbrace1{,}000, 5{,}000, 10{,}000, 20{,}000, 50{,}000, 100{,}000, 150{,}000, 200{,}000\rbrace$
- **Large-p block** with `max_dummy_multiplier = 5`, $n \in \lbrace300, 500\rbrace$ and
  $p \in \lbrace1{,}000{,}000, 5{,}000{,}000\rbrace$

$\Rightarrow$ 48 + 4 grid points × 3 solvers. A point runs out-of-core when the data volume one T-iteration touches
exceeds RAM. Each T-iteration runs all $K = 20$ random experiments, and each reads $\mathbf{X}$ and its own
memory-mapped dummy file of $L \cdot p$ columns:

```math
V = n \cdot p \cdot 8 \cdot (1 + K L) \ \text{bytes}.
```

At $L = 1$ this gives the following volumes in GiB (**bold** = above the 64 GiB RAM of the reference machine,
the columns $p \le 20{,}000$ stay below 16 GiB and are omitted):

| $n$ \\ $p$ | 50,000 | 100,000 | 150,000 | 200,000 | 1,000,000 | 5,000,000 |
| --- | --- | --- | --- | --- | --- | --- |
| 300 | 2.3 | 4.7 | 7.0 | 9.4 | 46.9 | **234.7** |
| 500 | 3.9 | 7.8 | 11.7 | 15.6 | **78.2** | **391.2** |
| 1,000 | 7.8 | 15.6 | 23.5 | 31.3 | — | — |
| 1,500 | 11.7 | 23.5 | 35.2 | 46.9 | — | — |
| 2,000 | 15.6 | 31.3 | 46.9 | 62.6 | — | — |
| 5,000 | 39.1 | **78.2** | **117.3** | **156.5** | — | — |

The output tables report the volume at the realized mean $L$ (`VolGiB`) and mark points above physical RAM with
`*`. The worst-case dummy disk at the cap reaches about 2 terabytes (TB) at $n = 500$, $p = 5{,}000{,}000$. Before
each point the demo checks that this worst case plus a 50 GiB reserve fits the free space of the temp directory
and skips the point with a message otherwise, because the library does not check free disk and a write into a
memory-mapped file on a full disk ends the process with `SIGBUS`.

---

## Control Parameters

Every control option is set **explicitly and identically** to the C++ demo
(nothing left to the Python helper's defaults — which differ from the C++ struct
and are what made the two demos diverge):

```text
K                        = 20        # Random experiments per T-loop iteration
max_dummy_multiplier     = 10 / 5    # Cap on L (base block / large-p block); drives the D volume
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
RAM at a time, which is essential out-of-core.

---

## Solvers Compared

Three base solvers spanning the least-angle regression (LARS) path and greedy
families (Demo 05 configuration):

- **TLARS** — terminating LARS, an equiangular path that stops on its own
- **TOMP** — terminating orthogonal matching pursuit (OMP), a greedy solver
- **TAFS** — terminating adaptive forward selection (AFS), `rho_afs = 0.3` (TLARS/TOMP take `0`)

---

## Measurements

Per $(n, p)$ × solver, averaged/aggregated over the `num_MC` sequential trials:

| Column | Meaning |
| --- | --- |
| `VolGiB` | Data volume per T-iteration $V$ at the mean realized $L$ (GiB) |
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
median is the robust summary and the min–max band shows the spread. `PeakRSS_MiB` is the running maximum of the
whole process, so the two blocks are best run in separate invocations.

---

## Output Files

Written to `simulation_results/data/`, **re-derived from the trial cache and re-saved after each completed grid
point** (a crash or out-of-memory (OOM) kill on a later, larger point preserves everything finished so far). The
base block uses the stem `demo_trex_08_scalability_results`, the large-p block
`demo_trex_08_scalability_large_p`.

- **`<stem>_trials.csv`** — the **raw per-trial cache**
  (source of truth): one row per $(n, p, \text{solver}, \text{mc})$ with its seed
  and every measured quantity. Appended to as each trial finishes.
- **`<stem>.txt`** — per-solver aligned table, one row
  per $(n, p)$ grid point across the 14 metric columns, with out-of-core points
  marked `*`.
- **`<stem>.csv`** — tidy long summary, header column
  order `n,p,is_barrier,solver,num_mc,snr,metric,value` (`num_mc` is the actual
  cached trial count per cell). Both summaries are derived from the cache.

The smoke grid writes `demo_trex_08_scalability_smoke.{txt,csv}` and keeps its trial cache in the temp directory,
where the startup cleanup removes it, so every smoke run computes fresh.

---

## Expandability

The benchmark is designed to grow incrementally rather than be re-run from
scratch. Every trial is keyed by $(n, p, \text{solver}, \text{mc})$ in the raw
cache, and a run computes **only the trials it is missing**:

- **num_MC top-up** — run at `num_MC = 2` now, then `num_MC = 10` later computes
  just the 8 new trials per cell and re-aggregates, reusing the first 2.
- **Grid growth** — append new $n$ / $p$ values to a block's value lists
  later. Existing cells are reused, and only the new ones are computed.
- **Trial-level resume** — a crash or OOM mid-grid-point keeps every trial that
  already finished (each is flushed to the cache immediately).

Per-trial seeds are a stable function of $(n, p, \text{mc})$ —
`24 + n·100003 + p·13 + mc` (mod $2^{31}{-}1$), independent of grid position and
identical to the C++ and R demos — so reuse and reproducibility survive reordering or
extending the grid. (The languages seed identically, but the *draws* still
differ because NumPy's pseudo-random number generator (PRNG) and C++'s
`std::mt19937` are different generators.)

---

## Results Visualization

Demo 08's comma-separated values (CSV) files use the wider
`n,p,is_barrier,solver,…` schema rather than the suite's tidy
`solver,metric,snr,value`, so the suite plotter handles them in a dedicated
`scalability` mode ([`../trex_plt_utils.py`](../trex_plt_utils.py), identical to the C++ copy,
`scalability <csv> [<csv> …]`). The CSVs of both blocks are read by `read_scalability()`, concatenated, and
rendered by `plot_scalability_dashboard()` as a 2×2 dashboard:

1. **Runtime vs $p$** (colour = $n$ and linestyle = solver, with a median + min–max band).
2. **Memory vs $p$** — peak RSS per $n$ against the $\mathbf{X}/\mathbf{y}$ backing
   size, with the physical-RAM line.
3. **Out-of-core money plot** — peak RSS against the volume per T-iteration (log–log), the
   64 GiB RAM line drawn in, and the out-of-core points ringed.
4. **FDR control at scale** — per solver, against the target FDR (tFDR) line.

`./generate_plots.sh` renders it into `simulation_results/plots/`.

---

## Status of the Recorded Data

The committed `demo_trex_08_scalability_results.{txt,csv}` is a partial **first pass at `num_MC = 2`** from July
on an older library commit ($n \in \lbrace1{,}000, 1{,}500, 2{,}000\rbrace$ only, with the old `DesignGiB`
column). A full run on the current library pin replaces it. The committed
`demo_trex_08_scalability_smoke.{txt,csv}` is the end-to-end **smoke recording** on the current code
(tiny $n \in \lbrace150, 200\rbrace$, $p \in \lbrace500, 1{,}000\rbrace$, caps 10 and 5) that
exercises the full memory-mapped pipeline and the table/CSV format — a pipeline
check, not a scientific result.

---

## Running

```bash
demo=Python/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability/demo_trex_08_mc_sim_scalability.py

python "$demo" base      # base block, num_MC = 10 (long-running)
python "$demo" large     # large-p block, num_MC = 10 (needs up to 2 TB free disk)
python "$demo" large 2   # large-p block, num_MC = 2 (first pass; top up later)
python "$demo"           # both blocks in one process, base first
python "$demo" smoke     # tiny end-to-end smoke check (num_MC = 2)
python "$demo" smoke 3   # smoke grid, num_MC = 3
```

The optional integer argument sets `num_MC`, and re-running with a larger value
tops up the trial cache rather than recomputing. Grid points run smallest-first and the
summaries are re-derived from the cache after each point, so an OOM/crash on a
later point preserves everything finished so far. The benchmark installs cleanup
handlers (`atexit` + `SIGINT`/`SIGTERM`) that sweep orphaned memory-mapped files
(`demo08_scalability_*`, `trex_memmap_*`, `trex_warmstart_*`) on exit. `SIGKILL` and `SIGBUS` cannot be handled
from Python, so the next startup sweep reclaims what they leave behind.
