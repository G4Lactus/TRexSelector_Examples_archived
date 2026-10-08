# ==============================================================================
# mc_runtime.R
# ==============================================================================
#
# Monte Carlo runtime shared by the R T-Rex demos:
#
#   mc_map()           — runs a function over tasks on a persistent PSOCK
#                        cluster (spawned R processes, no fork), results in
#                        task order. Every worker loads the package and runs
#                        its solvers single-threaded (6-core rule: the
#                        parallelism lives in the workers).
#   mc_stop_cluster()  — stops the persistent cluster (also done at exit).
#   cell_cache_open()  — resumable grid cells: every finished cell (method x
#   cell_cached()        grid point) is appended to a cache file, and a
#   cell_cache_done()    restarted demo reads it back instead of recomputing.
#
# The cluster is created on the first mc_map() call and receives the globals
# of the demo script at that moment, so the demo sources its helpers and
# defines its functions before the first Monte Carlo loop runs.
#
# The cache lives in <demo>/simulation_results/cache/ (gitignored). Its first
# line holds a fingerprint of the demo configuration and the run tag
# TREX_RUN_TAG (set by the overnight runners to the library pin). A cache with
# another fingerprint is set aside as <file>.stale. A demo that finishes calls
# cell_cache_done(), which removes the cache, so a later run starts fresh.
# ==============================================================================

.mc_rt <- new.env(parent = emptyenv())


# ==============================================================================
# Persistent PSOCK cluster
# ==============================================================================

#' Persistent cluster with `workers` processes, created on first use.
.mc_cluster <- function(workers) {
  if (!is.null(.mc_rt$cl) && identical(.mc_rt$workers, workers)) return(.mc_rt$cl)
  mc_stop_cluster()
  cl <- parallel::makePSOCKcluster(workers)
  globals <- setdiff(ls(globalenv(), all.names = TRUE), ".mc_rt")
  parallel::clusterExport(cl, globals, envir = globalenv())
  parallel::clusterEvalQ(cl, {
    suppressMessages(library(TRexSelectorNeo))
    TRexSelectorNeo::trex_set_num_threads(1L)
  })
  .mc_rt$cl <- cl
  .mc_rt$workers <- workers
  if (is.null(.mc_rt$finalizer)) {
    reg.finalizer(.mc_rt, function(e) mc_stop_cluster(), onexit = TRUE)
    .mc_rt$finalizer <- TRUE
  }
  cl
}

#' Stop the persistent cluster (no-op when none is running).
mc_stop_cluster <- function() {
  if (!is.null(.mc_rt$cl)) {
    try(parallel::stopCluster(.mc_rt$cl), silent = TRUE)
    .mc_rt$cl <- NULL
    .mc_rt$workers <- NULL
  }
  invisible(NULL)
}

#' Apply `fun` to every element of `tasks` on `workers` processes.
#'
#' One worker (or a single task) runs in this process. Errors in a worker stop
#' the demo with the worker's message. Load balancing hands out one task at a
#' time, because trial run times vary with the selected support.
#'
#' @param tasks   List or vector of task arguments.
#' @param fun     Function of one task.
#' @param workers Number of worker processes.
#' @return List of results in task order.
mc_map <- function(tasks, fun, workers = 1L) {
  TRexSelectorNeo::trex_set_num_threads(1L)
  if (workers <= 1L || length(tasks) <= 1L) return(lapply(tasks, fun))
  parallel::parLapplyLB(.mc_cluster(as.integer(workers)), tasks, fun)
}


# ==============================================================================
# Resumable grid cells
# ==============================================================================

#' Open (or start) the cell cache of one demo.
#'
#' @param dir         Cache folder, usually <demo>/simulation_results/cache.
#' @param stem        Demo stem; the file is <stem>.cells.tsv.
#' @param config      Demo configuration (any R object); its deparsed text and
#'                    TREX_RUN_TAG form the fingerprint.
#' @return Cache handle (environment) for cell_cached() / cell_cache_done().
cell_cache_open <- function(dir, stem, config) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, paste0(stem, ".cells.tsv"))
  header <- paste0("# ", paste(deparse(config), collapse = " "),
                   " | tag=", Sys.getenv("TREX_RUN_TAG"))
  cells <- list()
  if (file.exists(path)) {
    lines <- readLines(path, warn = FALSE)
    if (length(lines) > 0L && identical(lines[1L], header)) {
      for (ln in lines[-1L]) {
        parts <- strsplit(ln, "\t", fixed = TRUE)[[1L]]
        n <- length(parts)
        # A line cut off by a hard stop lacks the end marker and is recomputed.
        if (n < 3L || parts[n] != ".") next
        kv <- strsplit(parts[2L:(n - 1L)], "=", fixed = TRUE)
        if (any(lengths(kv) != 2L)) next
        vals <- suppressWarnings(as.numeric(vapply(kv, `[`, "", 2L)))
        cells[[parts[1L]]] <- setNames(vals, vapply(kv, `[`, "", 1L))
      }
      cat(sprintf("[cache] %d finished cell(s) restored from %s\n\n",
                  length(cells), path))
    } else {
      file.rename(path, paste0(path, ".stale"))
      cat(sprintf("[cache] configuration changed, old cache kept as %s.stale\n\n",
                  basename(path)))
    }
  }
  if (!file.exists(path)) writeLines(header, path)
  cache <- new.env(parent = emptyenv())
  cache$path <- path
  cache$cells <- cells
  cache
}

#' Result of one grid cell: from the cache when present, else computed.
#'
#' @param cache   Handle from cell_cache_open(), or NULL to disable caching.
#' @param key     Unique cell label without tabs or newlines.
#' @param compute Function without arguments returning a named numeric vector
#'                or a named list of numeric scalars.
#' @return Named list of the cell's numbers.
cell_cached <- function(cache, key, compute) {
  if (is.null(cache)) return(as.list(unlist(compute())))
  hit <- cache$cells[[key]]
  if (!is.null(hit)) {
    cat(sprintf("  [cache] %s\n", key))
    return(as.list(hit))
  }
  vals <- unlist(compute())
  line <- paste(c(key, sprintf("%s=%.17g", names(vals), vals), "."), collapse = "\t")
  cat(line, "\n", sep = "", file = cache$path, append = TRUE)
  cache$cells[[key]] <- vals
  as.list(vals)
}

#' Remove the cache after the demo has written all its results.
cell_cache_done <- function(cache) {
  if (!is.null(cache) && file.exists(cache$path)) unlink(cache$path)
  invisible(NULL)
}
