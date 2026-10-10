# ==============================================================================
# demo_trex_08_mc_sim_scalability.R
# ==============================================================================
#
# Scalability benchmark for the T-Rex Selector: runtime and memory usage as the
# problem scales, across three base solvers (TLARS, TOMP, TAFS).
#
# Mirrors:
# cpp/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability.cpp
# and
# Python/trex_selector_methods/ad/trex_core/demo_trex_08_mc_sim_scalability.py
#
# Two aims:
#   (a) high-dimensionality  — sweep p up to 200,000;
#   (b) volumes beyond RAM   — the largest p columns push the design [X|D] past
#                              physical RAM, feasible only because the pipeline
#                              is fully out-of-core (memory-mapped).
#
# Design of the study:
#   * 2-D grid  n in {1000,1500,2000,5000} x p in {1k,5k,10k,20k,50k,100k,150k,200k}
#   * SNR = 1 fixed (SNR is not a scaling axis; sweeping it only burns time),
#     demo-03 DGP (variable support s=10 redrawn per trial, unit coefficients,
#     Normal predictors/noise).
#   * num_MC = 10 SEQUENTIAL trials per (n, p, solver). Trials are NOT run in
#     parallel: this is a timing/memory benchmark, and concurrent trials would
#     (i) corrupt per-trial wall-clock, (ii) confound RSS, and (iii) multiply
#     the on-disk dummy backing.
#
# Fully memory-mapped pipeline in every point: X is generated column-chunk by
# column-chunk into an mmap backing file via dgp_chunked_gauss_snr() (never
# materialized in RAM), and use_memory_mapping = TRUE memory-maps the dummy
# matrices D and serializes solver warm-start state between T-loop iterations.
#
# Per grid point we record: timing (DGP, select() median/min/max, per-T-iter),
# memory (post-select RSS, physical footprint, process peak RSS, X+y bytes,
# theoretical design volume), and FDR / TPR / Avg L / Avg T.
#
# Memory metrics in R (see the profiling section):
#   * rss_mib      — `ps -o rss= -p <pid>` (current resident set size).
#   * peak_rss_mib — Linux: exact VmHWM from /proc/self/status. Elsewhere
#                    (macOS): running maximum of the ps-based RSS SAMPLES taken
#                    at the instrumentation points — a lower bound of the true
#                    kernel high-water mark (base R exposes no getrusage()).
#   * fp_mib       — Linux: VmRSS + VmSwap. macOS: NA (the phys_footprint
#                    counter needs mach task_info(), unreachable from base R;
#                    NA rather than a fake number).
#
# Self-cleanup: the process sweeps its own memory-mapped orphans
# (demo08_scalability_* backing files and trex_memmap_* dummy dirs in the temp
# dir) at startup and on exit, so a cancelled run does not strand tens of GB.
# (A SIGKILL cannot be trapped; the next startup sweep reclaims it.)
#
# Expandable by design. Every trial is appended to a raw per-trial cache
# (<stem>_trials.csv, one row per (n, p, solver, mc) keyed by its seed). A run
# computes only the trials it is missing and derives the summary table/CSV from
# the cache, so:
#   * num_MC top-up  — run at num_MC=2 now, later num_MC=10 computes just the 8
#     new trials per cell and re-aggregates (no recompute of the first 2);
#   * grid growth    — append new n/p values later; existing cells are reused;
#   * trial-level resume — a crash/OOM mid-point keeps every finished trial.
# Per-trial seeds are a stable function of (n, p, mc) (grid-position
# independent), so reuse and reproducibility survive reordering or extending
# the grid.
#
# Results go to simulation_results/data/ as an aligned table (.txt) and a tidy
# summary CSV, both re-derived from the cache and re-saved after each (n, p)
# point.
#
# Run:
#   Rscript demo_trex_08_mc_sim_scalability.R           # full benchmark, num_MC=10
#   Rscript demo_trex_08_mc_sim_scalability.R 2         # full benchmark, num_MC=2
#   Rscript demo_trex_08_mc_sim_scalability.R smoke     # tiny end-to-end check
#   Rscript demo_trex_08_mc_sim_scalability.R smoke 3   # smoke grid, num_MC=3
#
# ==============================================================================

library(TRexSelectorNeo)

# ==============================================================================
# Setup
# ==============================================================================
this_dir_ <- tryCatch(
  dirname(normalizePath(sys.frame(1L)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("--file=", "", file_arg[1]))) else "."
  }
)
source(file.path(this_dir_, "..", "..", "..", "support_generators.R"))
source(file.path(this_dir_, "..", "trex_sim_utils.R"))

# ==============================================================================
# Configuration
# ==============================================================================

OUT_DIR <- file.path(this_dir_, "simulation_results", "data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

S_TRUE     <- 10L      # cardinality of the true support (redrawn per trial)
TFDR       <- 0.1
SNR        <- 1.0      # fixed; SNR is not a scaling axis
CHUNK_COLS <- 1000L    # column-chunk size for the out-of-core X generation
DUMMY_MULT <- 10L      # max_dummy_multiplier; drives the D volume (10 x p)

# NOTE (thread pinning): the C++/Python demos pin the selector's inner OpenMP
# pool to TREX_NUM_WORKERS threads (default 8). The R
# package exports no equivalent, so the library's own defaults apply here —
# a documented divergence, not a config difference in trex_control().

# 2-D scaling grid. The largest p columns push the design [X|D] past physical
# RAM — the out-of-core regime that only the memory-mapped pipeline can reach.
N_VALUES <- c(1000L, 1500L, 2000L, 5000L)
P_VALUES <- c(1000L, 5000L, 10000L, 20000L, 50000L, 100000L, 150000L, 200000L)
NUM_MC   <- 10L

# Tiny end-to-end smoke grid (pipeline verification only).
SMOKE_GRID <- list(c(150L, 500L), c(200L, 1000L))
SMOKE_MC   <- 2L

# Base solvers to compare (rho_afs only used by TAFS, which requires rho in
# (0, 1]; mirrors the C++ DemoSolverInfo list — TOMP default, TAFS rho_afs=0.3).
SOLVERS <- list(
  list(solver = "TLARS", name = "TLARS", lambda2 = 0.0, rho_afs = 0.0),
  list(solver = "TOMP",  name = "TOMP",  lambda2 = 0.0, rho_afs = 0.0),
  list(solver = "TAFS",  name = "TAFS",  lambda2 = 0.0, rho_afs = 0.3)
)

# Temp dir used for the X backing files and searched for mmap orphans.
# Matches Python tempfile.gettempdir() / C++ std::filesystem::temp_directory_path().
TMP_DIR <- {
  td <- Sys.getenv("TMPDIR", unset = "")
  if (!nzchar(td)) td <- "/tmp"
  normalizePath(td, mustWork = FALSE)
}


.design_gib <- function(n, p) {
  # Theoretical design volume [X|D] in GiB (X plus 10x-p dummies, float64).
  n * (1 + DUMMY_MULT) * as.numeric(p) * 8.0 / 1024^3
}


.physical_ram_gib <- function() {
  # Physical RAM in GiB (NA if unavailable) — used to flag out-of-core points.
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    out <- tryCatch(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE),
                    error = function(e) NA_character_)
    val <- suppressWarnings(as.numeric(out[[1]]))
    if (is.finite(val)) return(val / 1024^3)
  } else if (file.exists("/proc/meminfo")) {
    line <- grep("^MemTotal:", readLines("/proc/meminfo"), value = TRUE)
    val  <- suppressWarnings(as.numeric(sub("^MemTotal:\\s+(\\d+).*$", "\\1",
                                            line[[1]])))
    if (is.finite(val)) return(val / 1024^2)   # kB -> GiB
  }
  NA_real_
}

RAM_GIB <- .physical_ram_gib()


.make_scalability_ctrl <- function(solver_entry) {
  # T-Rex control (base config of demos 02/03/07) for one grid point.
  # EVERY control option is set explicitly and identically to the C++ demo's
  # make_scalability_trex_ctrl() — no reliance on trex_control() defaults.
  TRexSelectorNeo::trex_control(
    solver                   = solver_entry$solver,
    lambda2                  = solver_entry$lambda2,
    rho_afs                  = solver_entry$rho_afs,
    K                        = 20L,
    max_dummy_multiplier     = DUMMY_MULT,
    use_max_T_stop           = TRUE,
    opt_threshold            = 0.75,
    dummy_distribution       = trex_dummy_distribution("Normal"),
    lloop_strategy           = "STANDARD",
    tloop_stagnation_stop    = TRUE,
    tloop_max_stagnant_steps = 5L,
    use_openmp               = FALSE,          # sequential random experiments
    use_memory_mapping       = TRUE            # X mmap + D mmap + solver serialization
  )
}


.benchmark_grid <- function(num_MC = NUM_MC) {
  # List of (n, p, num_MC, is_barrier) points over the full n x p grid.
  # is_barrier flags points whose logical design [X|D] exceeds this machine's
  # physical RAM (feasible only out-of-core).
  grid <- list()
  for (n in N_VALUES) {
    for (p in P_VALUES) {
      grid[[length(grid) + 1L]] <- list(
        n = n, p = p, num_MC = num_MC,
        is_barrier = isTRUE(.design_gib(n, p) > RAM_GIB)
      )
    }
  }
  grid
}


.smoke_grid <- function(num_MC = SMOKE_MC) {
  lapply(SMOKE_GRID, function(np) {
    list(n = np[[1]], p = np[[2]], num_MC = num_MC, is_barrier = FALSE)
  })
}


# ==============================================================================
# Expandable per-trial cache (top-up num_MC / grid; trial-level resume)
# ==============================================================================

# Deterministic, grid-position-independent, 32-bit-safe seed. Distinct (n, p)
# already yield different data (the X draw has a different shape), so we only
# need determinism in (n, p, mc), not cross-cell seed uniqueness. Matches the
# C++/Python _trial_seed so the languages seed identically.
SEED_MOD <- 2147483647   # 2^31 - 1

.trial_seed <- function(n, p, mc) {
  as.integer((24 + as.numeric(n) * 100003 + as.numeric(p) * 13 + mc) %% SEED_MOD)
}

# Raw per-trial store: one row per (n, p, solver, mc). The summary .txt/.csv
# are derived from it, so a later run tops up missing trials without
# recomputation. mc is 0-based (mirrors the Python cache format).
RAW_FIELDS <- c("n", "p", "solver", "mc", "trial_seed", "dgp_s", "select_s",
                "rss_mib", "fp_mib", "peak_rss_mib", "data_mib",
                "fdp", "tpp", "L", "T")

.raw_cache_path <- function(file_stem) {
  file.path(OUT_DIR, paste0(file_stem, "_trials.csv"))
}

.cell_key <- function(n, p, solver_name) sprintf("%d|%d|%s", n, p, solver_name)

.load_trial_cache <- function(raw_path) {
  # Load the raw per-trial store into cache[[key]][[as.character(mc)]] = row.
  cache <- list()
  if (!file.exists(raw_path)) return(cache)
  df <- utils::read.csv(raw_path, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(df))) {
    key <- .cell_key(df$n[[i]], df$p[[i]], df$solver[[i]])
    if (is.null(cache[[key]])) cache[[key]] <- list()
    cache[[key]][[as.character(df$mc[[i]])]] <- as.list(df[i, RAW_FIELDS])
  }
  cache
}

.append_trial_row <- function(raw_path, row) {
  # Append one trial to the raw store (write header if the file is new).
  df     <- as.data.frame(row, stringsAsFactors = FALSE)[, RAW_FIELDS]
  is_new <- !file.exists(raw_path)
  utils::write.table(df, raw_path, sep = ",", row.names = FALSE,
                     col.names = is_new, append = !is_new, quote = FALSE)
}

.aggregate_cell <- function(cell_cache, n, p) {
  # Aggregate all cached trials of one (n, p, solver) cell into a summary list.
  # Aggregates over every cached trial (monotone: more trials -> a better
  # estimate), so re-running at a lower num_MC never discards work.
  mcs  <- sort(as.integer(names(cell_cache)))
  rows <- lapply(as.character(mcs), function(k) cell_cache[[k]])
  num  <- function(field) vapply(rows, function(r) as.numeric(r[[field]]), 0.0)

  fdp   <- num("fdp")
  tpp   <- num("tpp")
  L     <- num("L")
  T_    <- num("T")
  dgp   <- num("dgp_s")
  sel   <- num("select_s")
  titer <- sel / pmax(T_, 1.0)
  total <- dgp + sel

  list(
    design_gib   = .design_gib(n, p),
    sel_med_s    = stats::median(sel),
    sel_min_s    = min(sel),
    sel_max_s    = max(sel),
    avg_dgp_s    = mean(dgp),
    avg_titer_s  = mean(titer),
    avg_total_s  = mean(total),
    avg_select_s = mean(sel),
    peak_rss_mib = max(num("peak_rss_mib")),
    max_rss_mib  = max(num("rss_mib")),
    max_fp_mib   = if (all(is.na(num("fp_mib")))) NA_real_ else max(num("fp_mib"), na.rm = TRUE),
    data_mib     = as.numeric(rows[[length(rows)]][["data_mib"]]),
    avg_fdr      = mean(fdp),
    avg_tpr      = mean(tpp),
    avg_L        = mean(L),
    avg_T        = mean(T_),
    num_mc       = length(rows)
  )
}


# ==============================================================================
# Self-cleanup of memory-mapped orphans
# ==============================================================================

ORPHAN_GLOBS <- c("demo08_scalability_*", "trex_memmap_*")

sweep_mmap_orphans <- function(reason = "") {
  # Remove this demo's stray mmap backing in the temp dir (best-effort).
  # Targets `demo08_scalability_*` X backing files and `trex_memmap_*` dummy
  # dirs left behind by a killed/cancelled run. Assumes SOLO execution (this is
  # a sequential benchmark). Returns bytes freed, invisibly.
  freed <- 0
  for (pattern in ORPHAN_GLOBS) {
    for (path in Sys.glob(file.path(TMP_DIR, pattern))) {
      info <- file.info(path)
      if (isTRUE(info$isdir)) {
        files <- list.files(path, recursive = TRUE, full.names = TRUE)
        freed <- freed + sum(file.size(files), na.rm = TRUE)
        unlink(path, recursive = TRUE, force = TRUE)
      } else {
        if (is.finite(info$size)) freed <- freed + info$size
        unlink(path)
      }
    }
  }
  if (freed > 0) {
    cat(sprintf("[cleanup] swept %.1f GiB of mmap orphans%s\n",
                freed / 1024^3,
                if (nzchar(reason)) sprintf(" (%s)", reason) else ""))
  }
  invisible(freed)
}

# Sweep on normal interpreter exit (best-effort; SIGKILL cannot be trapped —
# the next startup sweep reclaims whatever is left).
.cleanup_env_ <- new.env(parent = emptyenv())
invisible(reg.finalizer(.cleanup_env_, function(e) sweep_mmap_orphans("exit"),
                        onexit = TRUE))


# ==============================================================================
# Memory profiling hooks (ps / /proc based — see file header for caveats)
# ==============================================================================

.rss_state_ <- new.env(parent = emptyenv())
.rss_state_$peak_mib <- 0.0

current_rss_mib <- function() {
  # Current resident set size (RSS) of this process in MiB (0.0 if unavailable).
  # Every call also feeds the sampled peak (see sampled_peak_rss_mib()).
  mib <- 0.0
  if (file.exists("/proc/self/statm")) {
    parts <- strsplit(readLines("/proc/self/statm", n = 1L), " ")[[1]]
    pages <- suppressWarnings(as.numeric(parts[[2]]))
    if (is.finite(pages)) mib <- pages * 4096 / 1024^2
  } else {
    out <- tryCatch(
      system2("ps", c("-o", "rss=", "-p", Sys.getpid()), stdout = TRUE),
      error = function(e) NA_character_
    )
    kib <- suppressWarnings(as.numeric(trimws(out[[1]])))
    if (is.finite(kib)) mib <- kib / 1024.0   # ps reports KiB
  }
  .rss_state_$peak_mib <- max(.rss_state_$peak_mib, mib)
  mib
}

peak_rss_mib <- function() {
  # Peak RSS of this process in MiB.
  # Linux: exact VmHWM (kernel high-water mark) from /proc/self/status.
  # Elsewhere (macOS): running max of the ps-based samples taken at the
  # instrumentation points — a LOWER BOUND of the true peak (documented
  # divergence; base R exposes no getrusage()).
  if (file.exists("/proc/self/status")) {
    line <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
    if (length(line) > 0L) {
      kib <- suppressWarnings(as.numeric(sub("^VmHWM:\\s+(\\d+).*$", "\\1",
                                             line[[1]])))
      if (is.finite(kib)) return(kib / 1024.0)
    }
  }
  current_rss_mib()          # refresh the sample before reporting
  .rss_state_$peak_mib
}

current_footprint_mib <- function() {
  # Physical footprint of this process in MiB.
  # Linux approximation: VmRSS + VmSwap. macOS: NA — the phys_footprint
  # counter (Activity Monitor's "Memory") needs mach task_info(), which base R
  # cannot reach; NA is recorded rather than a fake number.
  if (file.exists("/proc/self/status")) {
    lines <- readLines("/proc/self/status")
    hits  <- grep("^(VmRSS|VmSwap):", lines, value = TRUE)
    kib   <- suppressWarnings(sum(as.numeric(sub("^\\S+\\s+(\\d+).*$", "\\1",
                                                 hits))))
    if (is.finite(kib)) return(kib / 1024.0)
  }
  NA_real_
}


# ==============================================================================
# Single-trial runner (sequential; timing + memory instrumented)
# ==============================================================================

run_scalability_trial <- function(n, p, ctrl, tFDR, snr, seed, support, coefs,
                                  mmap_filestem) {
  # One MC trial on the fully memory-mapped pipeline. The X backing file is
  # removed via on.exit() (exception-safe; mirrors the C++ RAII MmapFileGuard).
  # Returns a list of per-trial measurements.
  m <- list()
  x_path <- paste0(mmap_filestem, "_X.dat")
  on.exit(unlink(x_path), add = TRUE)

  current_rss_mib()   # baseline sample for the sampled-peak tracker

  t0 <- proc.time()[["elapsed"]]
  # Chunked out-of-core DGP: X goes straight to the mmap backing file,
  # working memory is O(n * CHUNK_COLS).
  dat <- dgp_chunked_gauss_snr(n, p, support, x_path,
                               chunk_cols = CHUNK_COLS, coefs = coefs,
                               snr = snr, seed = seed)
  m$dgp_s <- proc.time()[["elapsed"]] - t0

  m$data_mib <- (file.size(x_path) + length(dat$y) * 8) / 1024^2
  current_rss_mib()   # post-DGP sample

  t1 <- proc.time()[["elapsed"]]
  sel <- TRexSelector$new(
    X       = dat$X,
    y       = dat$y,
    tFDR    = tFDR,
    seed    = -1L,     # hardware entropy per trial, as in all MC demos
    verbose = FALSE,
    control = ctrl
  )
  sel$select()
  m$select_s <- proc.time()[["elapsed"]] - t1

  m$rss_mib      <- current_rss_mib()
  m$fp_mib       <- current_footprint_mib()
  m$peak_rss_mib <- peak_rss_mib()   # process-lifetime max up to this trial

  m$fdp <- TRexSelectorNeo::trex_fdp(sel$selected_indices, support)
  m$tpp <- TRexSelectorNeo::trex_tpp(sel$selected_indices, support)
  m$L   <- as.numeric(sel$L)
  m$T   <- as.numeric(sel$T_stop)

  # Drop the selector and mmap view before the backing file is removed.
  rm(sel, dat)
  m
}


# ==============================================================================
# Grid-point loop (n, p, solver): sequential MC trials
# ==============================================================================

run_scalability_grid_point <- function(n, p, num_MC, is_barrier, solver_entry,
                                       tFDR, snr, rnd_coef,
                                       cardinality_true_support,
                                       cell_cache, raw_path) {
  # Ensure >= num_MC trials exist for this (n, p, solver) cell, then aggregate.
  # Only the missing trial indices [0, num_MC) are computed; each is appended
  # to the raw store immediately (trial-level durability).
  ctrl <- .make_scalability_ctrl(solver_entry)
  tag  <- if (is_barrier) "BARRIER" else sprintf("n%d_p%d", n, p)
  mmap_filestem <- file.path(TMP_DIR, paste0("demo08_scalability_", tag))

  design <- .design_gib(n, p)
  over   <- if (is_barrier) "  (> RAM — out-of-core)" else ""
  have   <- length(cell_cache)
  todo   <- setdiff(0L:(num_MC - 1L), as.integer(names(cell_cache)))
  reuse_note <- if (have > 0L) sprintf("  (reusing %d cached)", have) else ""
  cat(sprintf(
    "  [%s] n=%d, p=%d  (design [X|D] = %.1f GiB%s) — %d new / %d MC trials (mmap)%s ...\n",
    solver_entry$name, n, p, design, over, length(todo), num_MC, reuse_note
  ))

  for (mc in todo) {
    trial_seed <- .trial_seed(n, p, mc)

    # Variable support (demo 03 pattern; make_support_random adds the
    # +500000 offset internally, mirroring C++ draw_support).
    support <- make_support_random(cardinality_true_support, p,
                                   seed = trial_seed)
    coefs <- if (rnd_coef) {
      set.seed(trial_seed + 600000L)     # isolated coefficient seed stream
      stats::rnorm(cardinality_true_support)
    } else {
      rep(1.0, cardinality_true_support)
    }

    m <- run_scalability_trial(n, p, ctrl, tFDR, snr, trial_seed,
                               support, coefs, mmap_filestem)

    row <- list(n = n, p = p, solver = solver_entry$name, mc = mc,
                trial_seed   = trial_seed,
                dgp_s        = round(m$dgp_s, 6L),
                select_s     = round(m$select_s, 6L),
                rss_mib      = round(m$rss_mib, 3L),
                fp_mib       = round(m$fp_mib, 3L),
                peak_rss_mib = round(m$peak_rss_mib, 3L),
                data_mib     = round(m$data_mib, 3L),
                fdp = round(m$fdp, 6L), tpp = round(m$tpp, 6L),
                L = round(m$L, 3L), T = round(m$T, 3L))
    .append_trial_row(raw_path, row)          # durable before we continue
    cell_cache[[as.character(mc)]] <- row

    cat(sprintf(
      "    trial %d/%d: dgp=%.2fs  select=%.2fs  FDP=%.3f  TPP=%.3f  L=%.0f  T=%.0f  RSS=%.1f MiB  FP=%.1f MiB\n",
      mc + 1L, num_MC, m$dgp_s, m$select_s, m$fdp, m$tpp, m$L, m$T,
      m$rss_mib, m$fp_mib
    ))
  }

  agg <- .aggregate_cell(cell_cache, n, p)
  cat(sprintf(
    "  [%s] n=%d, p=%d done (%d trials). select median=%.2fs [%.2f, %.2f]  peakRSS=%.0f MiB  FDR=%.3f  TPR=%.3f\n\n",
    solver_entry$name, n, p, agg$num_mc,
    agg$sel_med_s, agg$sel_min_s, agg$sel_max_s,
    agg$peak_rss_mib, agg$avg_fdr, agg$avg_tpr
  ))
  list(agg = agg, cell_cache = cell_cache)
}


# ==============================================================================
# Output: aligned table (console + .txt) and tidy CSV
# ==============================================================================

# Per-grid-point metric columns (table order and CSV metric set):
# label, decimals, aggregate key.
METRICS <- list(
  list("DesignGiB",   1L, "design_gib"),
  list("SelMed_s",    3L, "sel_med_s"),
  list("SelMin_s",    3L, "sel_min_s"),
  list("SelMax_s",    3L, "sel_max_s"),
  list("DGP_s",       3L, "avg_dgp_s"),
  list("s/Titer",     4L, "avg_titer_s"),
  list("PeakRSS_MiB", 1L, "peak_rss_mib"),
  list("RSS_MiB",     1L, "max_rss_mib"),
  list("Foot_MiB",    1L, "max_fp_mib"),
  list("X+y_MiB",     1L, "data_mib"),
  list("FDR",         4L, "avg_fdr"),
  list("TPR",         4L, "avg_tpr"),
  list("Avg_L",       2L, "avg_L"),
  list("Avg_T",       2L, "avg_T")
)

save_and_print_scalability_results <- function(file_stem, tFDR, snr, grid,  # nolint: object_length_linter. C++ name
                                               solvers, results,
                                               num_MC_target) {
  # Print the scalability results as an aligned per-solver table (rows = (n, p)
  # grid points, columns = metrics) to console + .txt, and write a tidy
  # long-format CSV. `results[[solver_name]][[key]]` holds the agg lists.
  # Mirrors the C++/Python function of the same name.
  txt_path <- file.path(OUT_DIR, paste0(file_stem, ".txt"))
  con      <- file(txt_path, open = "wt")
  on.exit(close(con), add = TRUE)
  out      <- .make_dual_out(con)

  fmt_val <- function(v, prec, width) {
    # (R's sprintf allows at most one '*' per conversion — build the format.)
    if (is.na(v)) sprintf("%*s", width, "NA") else sprintf(sprintf("%%%d.%df", width, prec), v)
  }

  out(sprintf("\n%s\n", strrep("=", 70L)))
  out(sprintf("=== T-Rex Scalability Results (tFDR = %.2f, SNR = %.1f, num_MC = %d) ===\n",
              tFDR, snr, num_MC_target))
  out(sprintf("%s\n\n", strrep("=", 70L)))
  out("Fully memory-mapped pipeline (X mmap + D mmap + solver serialization).\n")
  out(sprintf("Physical RAM: %.0f GiB. Grid: n x p; the barrier row (p*) is the out-of-core point.\n\n",
              RAM_GIB))

  n_w   <- 8L
  p_w   <- 9L
  met_w <- 13L
  for (solver_entry in solvers) {
    grid_res <- results[[solver_entry$name]]
    out(paste0(solver_entry$name, "\n"))
    hdr <- sprintf("%-*s%-*s", n_w, "n", p_w, "p")
    for (mspec in METRICS) hdr <- paste0(hdr, sprintf("%*s", met_w, mspec[[1L]]))
    out(paste0(hdr, "\n", strrep("-", nchar(hdr)), "\n"))
    for (gp in grid) {
      agg <- grid_res[[.cell_key(gp$n, gp$p, solver_entry$name)]]
      if (is.null(agg)) next
      lbl_p <- paste0(gp$p, if (gp$is_barrier) "*" else "")
      row   <- sprintf("%-*d%-*s", n_w, gp$n, p_w, lbl_p)
      for (mspec in METRICS) {
        row <- paste0(row, fmt_val(agg[[mspec[[3L]]]], mspec[[2L]], met_w))
      }
      out(paste0(row, "\n"))
    }
    out("\n")
  }
  out("(* = out-of-core barrier point, design [X|D] > physical RAM)\n\n")

  close(con)
  on.exit(NULL)

  # Tidy long-format CSV: n, p, is_barrier, solver, num_mc, snr, metric, value
  # (num_mc is the actual cached trial count per cell, not the run target).
  rows <- list()
  for (solver_entry in solvers) {
    grid_res <- results[[solver_entry$name]]
    for (gp in grid) {
      agg <- grid_res[[.cell_key(gp$n, gp$p, solver_entry$name)]]
      if (is.null(agg)) next
      for (mspec in METRICS) {
        rows[[length(rows) + 1L]] <- data.frame(
          n = gp$n, p = gp$p, is_barrier = as.integer(gp$is_barrier),
          solver = solver_entry$name, num_mc = agg$num_mc,
          snr = sprintf("%.6f", snr), metric = mspec[[1L]],
          value = round(agg[[mspec[[3L]]]], 6L),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  csv_path <- file.path(OUT_DIR, paste0(file_stem, ".csv"))
  utils::write.csv(do.call(rbind, rows), csv_path,
                   row.names = FALSE, quote = FALSE)
  cat(sprintf("[Info] CSV results saved to:          %s\n", csv_path))
  cat(sprintf("[Info] Table results saved to:        %s\n\n", txt_path))
}


# ==============================================================================
# Scalability benchmark driver
# ==============================================================================

demo_trex_scalability <- function(grid, file_stem, num_MC_target,
                                  rnd_coef = FALSE) {
  # Clear any orphans stranded by a previous cancelled run before we start.
  sweep_mmap_orphans("startup")

  raw_path <- .raw_cache_path(file_stem)
  cache    <- .load_trial_cache(raw_path)   # cache[[key]][[mc]] = row
  cached_trials <- sum(vapply(cache, length, 0L))

  cat("\n", strrep("=", 70L), "\n", sep = "")
  cat("Demo: T-Rex Selector Scalability Benchmark  |  TLARS / TOMP / TAFS\n")
  cat(strrep("=", 70L), "\n\n", sep = "")

  max_design <- .design_gib(N_VALUES[[length(N_VALUES)]],
                            P_VALUES[[length(P_VALUES)]])
  # Out-of-core count over the FULL n x p grid (header info; mirrors Python).
  n_over <- sum(vapply(
    N_VALUES,
    function(n) {
      sum(vapply(P_VALUES,
                 function(p) isTRUE(.design_gib(n, p) > RAM_GIB),
                 FALSE))
    },
    0L
  ))
  cat("Fully memory-mapped pipeline (X mmap + D mmap + solver serialization)\n")
  cat(sprintf("Physical RAM: %.0f GiB\n", RAM_GIB))
  cat(sprintf("Grid: n in {%s} x p in {%s}  (%d points)\n",
              paste(N_VALUES, collapse = ", "),
              paste(P_VALUES, collapse = ", "),
              length(N_VALUES) * length(P_VALUES)))
  cat(sprintf("Largest design [X|D] = %.1f GiB (n=%d, p=%d); %d point(s) exceed RAM (out-of-core)\n",
              max_design, N_VALUES[[length(N_VALUES)]],
              P_VALUES[[length(P_VALUES)]], n_over))
  cat(sprintf("SNR = %.1f (fixed),  s = %d,  tFDR = %.2f,  num_MC = %d sequential per (n, p, solver)\n",
              SNR, S_TRUE, TFDR, num_MC_target))
  cat(sprintf("Trial cache: %s  (%d trial(s) already present — only missing trials are computed)\n\n",
              raw_path, cached_trials))

  # results[[solver_name]][[cell_key]] = agg
  results <- setNames(vector("list", length(SOLVERS)),
                      vapply(SOLVERS, `[[`, "", "name"))
  for (nm in names(results)) results[[nm]] <- list()

  # (n, p)-outer loop; every solver finishes the current point before the
  # next (larger) one starts, and results are re-saved after each point so a
  # crash on a bigger point preserves the finished ones.
  for (gi in seq_along(grid)) {
    gp <- grid[[gi]]
    cat(strrep("=", 50L), "\n", sep = "")
    cat(sprintf("Grid point: n=%d, p=%d%s\n", gp$n, gp$p,
                if (gp$is_barrier) "   [BARRIER — out-of-core]" else ""))
    cat(sprintf("  design [X|D] = %.1f GiB  (X backing %.2f GiB + 10x-p dummies)\n",
                .design_gib(gp$n, gp$p),
                gp$n * as.numeric(gp$p) * 8.0 / 1024^3))
    cat(strrep("=", 50L), "\n\n", sep = "")

    for (solver_entry in SOLVERS) {
      key        <- .cell_key(gp$n, gp$p, solver_entry$name)
      cell_cache <- if (is.null(cache[[key]])) list() else cache[[key]]
      point <- run_scalability_grid_point(
        gp$n, gp$p, gp$num_MC, gp$is_barrier, solver_entry, TFDR, SNR,
        rnd_coef, S_TRUE, cell_cache, raw_path
      )
      cache[[key]] <- point$cell_cache
      results[[solver_entry$name]][[key]] <- point$agg
    }

    # Incremental save over completed points (derived from the trial cache).
    save_and_print_scalability_results(file_stem, TFDR, SNR,
                                       grid[seq_len(gi)], SOLVERS, results,
                                       num_MC_target)
  }
}


# ==============================================================================
# Main
# ==============================================================================

# CLI:  [smoke] [num_MC]  (order-independent; num_MC is any integer token)
cli_args_ <- commandArgs(trailingOnly = TRUE)
smoke_    <- "smoke" %in% cli_args_
mc_tokens_ <- suppressWarnings(
  as.integer(cli_args_[grepl("^\\+?[0-9]+$", cli_args_)])
)
num_MC_cli_ <- if (length(mc_tokens_) > 0L) mc_tokens_[[1L]] else NA_integer_

if (smoke_) {
  mc_ <- if (!is.na(num_MC_cli_)) num_MC_cli_ else SMOKE_MC
  demo_trex_scalability(.smoke_grid(mc_),
                        "demo_trex_08_scalability_smoke", mc_)
} else {
  mc_ <- if (!is.na(num_MC_cli_)) num_MC_cli_ else NUM_MC
  demo_trex_scalability(.benchmark_grid(mc_),
                        "demo_trex_08_scalability_results", mc_)
}

sweep_mmap_orphans("shutdown")
cat("\ndemo_trex_08_mc_sim_scalability.R complete.\n")
