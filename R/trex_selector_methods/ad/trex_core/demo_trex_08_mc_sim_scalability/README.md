# Demo 08: Scalability Benchmark

## Purpose

Characterize Terminating Random Experiments (T-Rex) selector **runtime and memory usage** as $n$ and $p$ scale to
large sizes, and demonstrate that the selector goes **high-dimensional and
out-of-core**, into data volumes that exceed physical random-access memory (RAM). Everything runs on
the **fully memory-mapped pipeline**. Every Monte Carlo (MC) trial generates $\mathbf{X}$
directly into a memory-mapped (mmap) backing file (`dgp_chunked_gauss_snr` with
`chunk_cols = 1000`, so $\mathbf{X}$ is never materialized in RAM), and
`use_memory_mapping = TRUE` additionally memory-maps the internal dummy matrices
$\mathbf{D}$ and serializes solver warm-start state to disk. The selector receives the mapped file through its own
memory-mapped entry point, so the binding never copies $\mathbf{X}$. R port of the
[C++ reference study](../../../../../cpp/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability/README.md),
whose README explains the memory behaviour and the expected cost in detail.

## Data-generating process

The data-generating process (DGP) mirrors **Demo 03** (variable support). The design is independent and
identically distributed (i.i.d.) standard-normal, the $s = 10$ support indices are **redrawn per trial** via
`make_support_random`, the coefficients are fixed at $\beta_j = 1$, the signal-to-noise ratio (SNR) is fixed at
$1.0$ (deliberately not a scaling axis here), and $\mathrm{tFDR} = 0.1$. Trials run **sequentially** (no worker
pool) so wall-clock and resident set size (RSS) readings are undistorted. The selector runs with
`seed = -1` (hardware entropy) as in all MC demos.

Per-trial seeds are a stable function of $(n, p, \text{mc})$ —
`(24 + n*100003 + p*13 + mc) mod (2^31 - 1)` — independent of grid position and
identical to the C++/Python demos, so cache reuse and reproducibility survive
reordering or extending the grid.

The design family is the suite's shared one — independent columns, so the
empirical column-correlation matrix is essentially the identity at every grid
point. While $n$ and $p$ scale, the correlation structure does not:

![Design-matrix correlation of the shared i.i.d. Gaussian DGP (shown at n=300, p=1,000, natural order, no clustering)](../dgp_corr_plots/trex_gauss.png)

The $K$ random experiments inside `select()` run one after another, and the selector's Open Multi-Processing
(OpenMP) threads parallelize the column loops inside each experiment. The thread count is pinned to
`TREX_NUM_WORKERS` (default 8, `0` = all cores) through `trex_set_num_threads()`, as in the C++/Python demos.

## The grid

Two blocks with separate output files and trial caches, `num_MC = 10` sequential trials per
$(n, p, \text{solver})$, three solvers (**TLARS**, **TOMP**, **TAFS** with `rho_afs = 0.3`, and `lambda2 = 0`
everywhere):

- **Base block** with `max_dummy_multiplier = 10`, the full cross of
  $n \in \lbrace300, 500, 1{,}000, 1{,}500, 2{,}000, 5{,}000\rbrace$ and
  $p \in \lbrace1{,}000, 5{,}000, 10{,}000, 20{,}000, 50{,}000, 100{,}000, 150{,}000, 200{,}000\rbrace$
- **Large-p block** with `max_dummy_multiplier = 5`, $n \in \lbrace300, 500\rbrace$ and
  $p \in \lbrace1{,}000{,}000, 5{,}000{,}000\rbrace$

A point runs out-of-core when the data volume one T-iteration touches exceeds RAM. Each T-iteration runs all
$K = 20$ random experiments, and each reads $\mathbf{X}$ and its own memory-mapped dummy file of $L \cdot p$
columns:

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
each point the demo checks (via `df -Pk`) that this worst case plus a 50 GiB reserve fits the free space of the
temp directory and skips the point with a message otherwise, because the library does not check free disk and a
write into a memory-mapped file on a full disk ends the process with `SIGBUS`.

## Control

Every option explicit, identical to the C++ reference: `K = 20`,
`max_dummy_multiplier = 10` (base block) or `5` (large-p block), `use_max_T_stop = TRUE`, `opt_threshold = 0.75`,
`dummy_distribution = Normal`, `lloop_strategy = "STANDARD"`,
`tloop_stagnation_stop = TRUE`, `tloop_max_stagnant_steps = 5`,
`use_openmp = FALSE` (sequential random experiments — one experiment working set
in RAM at a time, essential out-of-core),
`use_memory_mapping = TRUE`.

## Measurements

Per $(n, p)$ x solver, aggregated over the cached trials: `VolGiB` (data volume per T-iteration at the mean
realized $L$), `SelMed_s`/`SelMin_s`/`SelMax_s` (median/min/max `select()` wall-clock), `DGP_s`,
`s/Titer`, `PeakRSS_MiB`, `RSS_MiB`, `Foot_MiB`, `X+y_MiB`, and
`FDR`/`TPR`/`Avg_L`/`Avg_T` as selection-quality sanity checks.

Memory metrics in R (documented divergences from the Python/C++ instrumentation):

- `RSS_MiB` — `ps -o rss= -p <pid>` (macOS) / `/proc/self/statm` (Linux).
- `PeakRSS_MiB` — Linux: exact `VmHWM`. macOS: running maximum of the RSS
  **samples** taken at the instrumentation points — a lower bound of the true
  kernel high-water mark (base R exposes no `getrusage()`).
- `Foot_MiB` — Linux: `VmRSS + VmSwap`. macOS: **not available (NA)**. The `phys_footprint`
  counter needs mach `task_info()`, which is unreachable from base R, so NA is recorded
  rather than a fake number.

## Output files and expandability

Written to `simulation_results/data/`, re-derived from the trial cache and re-saved after each completed grid
point. The base block uses the stem `demo_trex_08_scalability_results`, the large-p block
`demo_trex_08_scalability_large_p`, and smoke mode `demo_trex_08_scalability_smoke`.

- `<stem>_trials.csv` — the **raw per-trial cache** (source of truth): one row per
  $(n, p, \text{solver}, \text{mc})$ with columns
  `n,p,solver,mc,trial_seed,dgp_s,select_s,rss_mib,fp_mib,peak_rss_mib,data_mib,fdp,tpp,L,T`.
  Appended to as each trial finishes. Smoke mode keeps it in the temp directory, where the startup cleanup
  removes it, so every smoke run computes fresh.
- `<stem>.txt` — per-solver aligned table, one row per $(n, p)$ grid point, with
  out-of-core points marked `*`.
- `<stem>.csv` — tidy long summary, header `n,p,is_barrier,solver,num_mc,snr,metric,value`
  (`num_mc` is the actual cached trial count per cell).

A run computes **only the trials it is missing** from the cache, so `num_MC` top-up,
grid growth, and trial-level resume all reuse finished trials. The script sweeps
orphaned `demo08_scalability_*` backing files, `trex_memmap_*` dummy dirs and `trex_warmstart_*` checkpoint dirs
from the temp dir at startup and on exit. `SIGKILL` and `SIGBUS` cannot be handled from R, so the next startup
sweep reclaims what they leave behind.

## Running

```bash
demo=R/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability/demo_trex_08_mc_sim_scalability.R

Rscript "$demo" base       # base block, num_MC = 10 (long-running)
Rscript "$demo" large      # large-p block, num_MC = 10 (needs up to 2 TB free disk)
Rscript "$demo" large 2    # large-p block, num_MC = 2 (first pass; top up later)
Rscript "$demo"            # both blocks in one process, base first
Rscript "$demo" smoke      # tiny end-to-end smoke check (num_MC = 2)
Rscript "$demo" smoke 3    # smoke grid, num_MC = 3
```

The block, `smoke` and integer arguments are order-independent, and the smoke grid is
$(n, p) \in \lbrace(150, 500), (200, 1{,}000)\rbrace$ with dummy caps 10 and 5. The script sources
`../../../support_generators.R` and `../trex_sim_utils.R`.

Unlike the other demos in this suite, this one is **run manually** (it is not part
of the overnight runner), and so far only the smoke recording is committed. Run the script
first, then `./generate_plots.sh` to render the scalability dashboard of runtime, memory and false discovery rate
(FDR) versus problem size from the wide summary comma-separated values (CSV) files of both blocks via the shared
[`../trex_plt_utils.R`](../trex_plt_utils.R) into `simulation_results/plots/`.
