# Demo 08: Scalability Benchmark

## Purpose

Characterize **Terminating Random Experiments (T-Rex) selector runtime and memory usage** as $n$ and $p$ scale to
large sizes, and demonstrate that the selector goes **high-dimensional and
out-of-core**, into data volumes that exceed physical random-access memory
(RAM). Everything runs on the **fully memory-mapped pipeline**. Every Monte Carlo
(MC) trial generates its $\mathbf{X}/\mathbf{y}$ into memory-mapped (mmap) backing files
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

Two blocks with separate output files and trial caches, `num_MC` sequential trials per $(n, p, \text{solver})$,
default 10 and settable through the command-line interface (CLI):

- **Base block** with `max_dummy_multiplier = 10`, the full cross of
  - $n \in \lbrace300, 500, 1{,}000, 1{,}500, 2{,}000, 5{,}000\rbrace$
  - $p \in \lbrace1{,}000, 5{,}000, 10{,}000, 20{,}000, 50{,}000, 100{,}000, 150{,}000, 200{,}000\rbrace$
- **Large-p block** with `max_dummy_multiplier = 5`, few observations and millions of predictors
  - $n \in \lbrace300, 500\rbrace$
  - $p \in \lbrace1{,}000{,}000, 5{,}000{,}000\rbrace$

$\Rightarrow$ 48 + 4 grid points × 3 solvers. The quantity that decides whether a point runs in RAM or out-of-core
is the **data volume one T-iteration touches**. Each T-iteration runs all $K = 20$ random experiments, and each
experiment reads $\mathbf{X}$ and its own memory-mapped dummy file of $L \cdot p$ columns, so

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

A larger $L$ multiplies the dummy part, so a point with $L = 2$ already needs roughly twice the volume. The output
tables report the volume at the realized mean $L$ (`VolGiB`) and mark points above physical RAM with `*`.

The cells marked — are not part of the grid. At $p = 5{,}000{,}000$ every $n \ge 1{,}000$ needs at least
8 terabytes (TB) of dummy disk at the cap, and at $p = 1{,}000{,}000$ the volume reaches 156 GiB or more.

The **disk** the memory-mapped dummy files may occupy in the worst case (all $K$ files filled up to the cap,
plus $\mathbf{X}$ and the solver checkpoints) is, in GiB,

| $n$ \\ $p$ | 50,000 | 100,000 | 150,000 | 200,000 | 1,000,000 | 5,000,000 |
| --- | --- | --- | --- | --- | --- | --- |
| 300 | 23 | 46 | 68 | 91 | 229 | 1,147 |
| 500 | 38 | 76 | 113 | 151 | 380 | **1,899** |
| 1,000 | 75 | 150 | 226 | 301 | — | — |
| 1,500 | 113 | 225 | 338 | 451 | — | — |
| 2,000 | 150 | 300 | 450 | 600 | — | — |
| 5,000 | 375 | 749 | 1,124 | 1,499 | — | — |

The cap of 5 in the large-p block keeps $n = 500$, $p = 5{,}000{,}000$ at about 2 TB. With the base cap of 10 the
same point would need 4 TB, more than the free disk of the reference machine.

The signal-to-noise ratio (SNR) is fixed at $1.0$ — it is deliberately **not** a scaling axis here (sweeping
it only burns time). That sweep belongs to Demos 02–07.

A tiny smoke grid ($n \in \lbrace150, 200\rbrace$, $p \in \lbrace500, 1{,}000\rbrace$, caps 10 and 5) is
available via the `smoke` argument for end-to-end pipeline verification.

### Memory behaviour

The whole pipeline is memory-mapped, so $\mathbf{X}/\mathbf{y}$ live on disk, the dummy matrices $\mathbf{D}$ are
memory-mapped, and solver warm-start state is serialized between T-loop iterations. Neither binding nor library
copies $\mathbf{X}$. The selector centers and normalizes it in place on the mapped file and restores it at the
end of `select()`.

The **anonymous heap** stays small. The coefficient path of every experiment is stored sparsely (one index and
value per active coefficient and step) and is reduced to the voting statistics right after the experiment, so
nothing of size $(p + L p) \times \text{steps}$ is ever allocated. What remains is about 16–32 bytes per column of
$[\mathbf{X}\,|\,\mathbf{D}]$ for the running experiment (correlations, inactive set, compensation terms) and a
few $p \times T$ accumulators of the voting statistics. The `Foot_MiB` column tracks this heap.

The **resident** memory is mostly page cache of the mapped files. While the volume $V$ fits into the page cache,
every T-iteration reads the dummy files from RAM. Once it does not, the dummy files are read back from the
solid-state drive (SSD) on every T-iteration and `select()` slows down by a factor of 5–10. In the July first
pass this began between about 35 and 47 GiB of volume at $L = 1$, below the 64 GiB of RAM, because the page cache
shares RAM with the rest of the system and trials with $L = 2$ double the volume.

The library does not check free disk before it creates the dummy files, and a write into a memory-mapped file on
a full disk ends the process with `SIGBUS` instead of an error. The demo therefore checks every point before its
trials. A point whose worst-case disk need plus a 50 GiB reserve does not fit the free space of the temp
directory is skipped with a message.

### Expected cost of the large-p block

An estimate extrapolated from the July first pass (base grid only, older library), not a measurement. At
$L = 1$ one `select()` should take about 2–6 min at $n = 300$, $p = 1{,}000{,}000$ and about 30–60 min at
$n = 500$, $p = 5{,}000{,}000$, and a point with $L = 2$ costs about three times as much. The decisive unknown is
$L$. With $n \le 500$ the true signals at SNR $= 1$ stand only slightly above the largest of millions of null
correlations, and the smoke grid ($n \le 200$) already climbs to $L = 4$–10. A calibration run at
$p = 1{,}000{,}000$ measures $L$ and $T$ before the $p = 5{,}000{,}000$ points are started.

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

Trials run **sequentially** so wall-clock and RSS measurements are undistorted. The $K$ random experiments inside
`select()` also run one after another, and the `TREX_NUM_WORKERS` Open Multi-Processing (OpenMP) threads
(default 8) parallelize the column loops inside each experiment (correlations, dummy generation, normalization).

---

## Control Parameters

Every control option is set **explicitly and identically** in the C++, Python and R
demos (nothing left to a struct/helper default — that is what made the ports
diverge):

```text
K                        = 20        # Random experiments per T-loop iteration
max_dummy_multiplier     = 10 / 5    # Cap on L (base block / large-p block); drives the D volume
use_max_T_stop           = true      # Cap T <= ceil(n/2)
opt_threshold            = 0.75      # Voting/optimization threshold (C++ 0.75, Python default 0.02 — pinned)
dummy_distribution       = Normal    # Dummy predictors drawn from N(0,1)
lloop_strategy           = STANDARD  # Fresh i.i.d. dummy matrix per L-loop iteration
tloop_stagnation_stop    = true      # Early exit when R_mat stagnates
tloop_max_stagnant_steps = 5         # Stagnation window
parallel_rnd_experiments = false     # Sequential experiments; one D working set at a time
use_memory_mapping       = true      # All points: X/y mmap + D mmap + solver serialization
tFDR                     = 0.1       # Target FDR control level
```

`parallel_rnd_experiments = false` is essential here, because the alternative (parallel
across the K experiments) holds up to *K* dummy-augmented working sets in RAM at
once, which defeats the out-of-core goal. `STANDARD` writes one memory-mapped dummy file per experiment and
regenerates all of them in every L-loop iteration, while the T-loop only reads them.

---

## Solvers Compared

Three base solvers spanning the least-angle regression (LARS) path and greedy families (Demo 05 configuration):

- **TLARS** — terminating LARS, an equiangular path that stops on its own
- **TOMP** — terminating orthogonal matching pursuit (OMP), a greedy solver
- **TAFS** — terminating adaptive forward selection (AFS), `rho_afs = 0.3` (TLARS/TOMP take `0`)

---

## Measurements

Per $(n, p)$ × solver, averaged/aggregated over the cached trials:

| Column | Meaning |
| --- | --- |
| `VolGiB` | Data volume per T-iteration $V$ at the mean realized $L$ (GiB) |
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

`PeakRSS_MiB` is the running maximum of the whole process, so every point after the largest one so far repeats
its value. Running the two blocks in separate invocations keeps the base block's peak out of the large-p table.

---

## Output Files

Written to `simulation_results/data/`, **re-derived from the trial cache and re-saved after each completed grid
point** (a crash or out-of-memory kill on a later, larger point preserves everything finished so far). The base
block uses the stem `demo_trex_08_scalability_results`, the large-p block `demo_trex_08_scalability_large_p`.

- **`<stem>_trials.csv`** — the **raw per-trial cache** (source of truth): one row per
  $(n, p, \text{solver}, \text{mc})$ with its seed and every measured quantity, appended as each trial finishes.
- **`<stem>.txt`** — per-solver aligned table, one row per $(n, p)$ grid point across the 14 metric columns,
  with out-of-core points marked `*`.
- **`<stem>.csv`** — tidy long summary, header column order `n,p,is_barrier,solver,num_mc,snr,metric,value`
  (`num_mc` is the actual cached trial count per cell, `is_barrier` = 1 where `VolGiB` exceeds RAM).

The smoke grid writes `demo_trex_08_scalability_smoke.{txt,csv}` and keeps its trial cache in the temp directory,
where the startup cleanup removes it, so every smoke run computes fresh.

---

## Expandability

The benchmark is designed to grow incrementally rather than be re-run from
scratch. Every trial is keyed by $(n, p, \text{solver}, \text{mc})$ in the raw
cache, and a run computes **only the trials it is missing**:

- **num_MC top-up** — run at `num_MC = 2` now, then `num_MC = 10` later computes
  just the 8 new trials per cell and re-aggregates, and the first 2 are reused.
- **Grid growth** — append new $n$ / $p$ values to a block's value lists later,
  so existing cells are reused and only the new ones are computed.
- **Trial-level resume** — a crash or out-of-memory kill mid-grid-point keeps every trial that
  already finished (each is flushed to the cache immediately).

Per-trial seeds are a stable function of $(n, p, \text{mc})$ —
`24 + n·100003 + p·13 + mc` (mod $2^{31}{-}1$), independent of grid position and
**identical to the Python and R demos' `_trial_seed`** — so reuse and reproducibility
survive reordering or extending the grid. (The languages seed identically, but the
*draws* still differ because `std::mt19937`, NumPy's and R's pseudo-random number generators
(PRNGs) are different generators.)

---

## Results Visualization

Demo 08's comma-separated values (CSV) files use the wider `n,p,is_barrier,solver,num_mc,snr,metric,value` schema
rather than the suite's tidy `solver,metric,snr,value`, so the suite plotter
handles them in a dedicated `scalability` mode
([`../trex_plt_utils.py`](../trex_plt_utils.py) `scalability <csv> [<csv> …]`). The CSVs of both blocks are read
by `read_scalability()`, concatenated, and rendered by `plot_scalability_dashboard()`. The result is a 2×2
dashboard (runtime vs $p$, memory vs $p$, peak RSS against the volume with the RAM line and the out-of-core
points ringed, and FDR control at scale), in which the $n = 300$ and $n = 500$ lines run from
$p = 1{,}000$ to $5{,}000{,}000$. Render it with `./generate_plots.sh`.

---

## Status of the Recorded Data

The committed `demo_trex_08_scalability_results.{txt,csv}` and the dashboard are a **first pass at
`num_MC = 2`** from July on an older library commit, recorded before the $n = 300$ and $n = 500$ rows and the
large-p block existed. They still carry the old `DesignGiB` column (logical $[\mathbf{X}\,|\,\mathbf{D}]$ at the
cap), which the plotter reads as a fallback. A full run on the current library pin replaces them. The committed
`demo_trex_08_scalability_smoke.{txt,csv}` is the end-to-end smoke recording on the current code.

---

## Building and Running

Release build is recommended (the large points are heavy):

```bash
cd TRexSelector_Examples/cpp
cmake --build build/release --target demo_trex_08_mc_sim_scalability

bin=build/release/bin/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability/demo_trex_08_mc_sim_scalability
"$bin" base      # base block, num_MC = 10 (long-running)
"$bin" large     # large-p block, num_MC = 10 (needs up to 2 TB free disk)
"$bin" large 2   # large-p block, num_MC = 2 (first pass; top up later)
"$bin"           # both blocks in one process, base first
"$bin" smoke     # tiny end-to-end smoke check (num_MC = 2)
"$bin" smoke 3   # smoke grid, num_MC = 3

TREX_SOLVERS=TLARS "$bin" calib   # calibration: large-p block at p = 1,000,000, TLARS, num_MC = 1
```

The optional integer argument sets `num_MC`, and re-running with a larger value tops
up the trial cache rather than recomputing. `calib` runs only the $p = 1{,}000{,}000$ points of the large-p block
(default `num_MC = 1`), and its trials stay in the large-p cache as regular trials. `TREX_SOLVERS` (comma
separated) limits which solvers compute new trials, while the cached trials of the other solvers still enter the
tables. Grid points run smallest-first and the
summaries are re-derived from the cache after each point. The benchmark sweeps its
own memory-mapped orphans (`demo08_scalability_*`, `trex_memmap_*`, `trex_warmstart_*`) at startup and
on exit / `SIGINT` / `SIGTERM` / `SIGBUS` (`SIGKILL` cannot be trapped, so the next startup sweep
reclaims it). The memory-mapped files go to the temp directory (`TMPDIR`), which must sit on the volume with the
free disk shown in the console header.
