# ==============================================================================
# ts_perf_solvers.R
# ==============================================================================
#
# The T-solvers compared by the performance demos (demo_ts_perf_*), each with
# the settings of its own solver demo. Mirrors cpp/tsolvers/ts_perf_solvers.hpp.
#
# Every solver runs as a reader (normalize = FALSE, intercept = FALSE): the
# caller centers and scales X and D and centers y once in front, and all
# solvers share the same R matrices. The data are marked precentered, so the
# row-bounded active sets stop at the rank n - 1 of the centered design
# (T-ENET with lambda2 > 0 has no such cap). T-AFS appears once per rho.
#
# ==============================================================================

#' Name, family and factory of one compared solver; make(X, D, y) builds it as a reader.
.solver_spec <- function(cls, name, family, ...) {
  extra <- list(...)
  make <- function(X, D, y) {
    options <- do.call(tsolver_options, c(list(normalize = FALSE, intercept = FALSE, precentered = TRUE), extra))
    cls$new(X, D, y, options = options)
  }
  list(name = name, family = family, make = make)
}

#' The 12 demo solvers with their demo settings, T-AFS for rho in {0.1, 0.3, 0.5, 1}.
demo_solvers <- function() {
  lars <- "LARS-based"
  omp <- "OMP-based"
  afs <- "AFS-based"
  list(
    .solver_spec(TLARS_Solver, "T-LARS", lars),
    .solver_spec(TLASSO_Solver, "T-LASSO", lars),
    .solver_spec(TSTEPWISE_Solver, "T-Stepwise", lars),
    .solver_spec(TENET_Solver, "T-ENET (lambda2 0.5)", lars, lambda2 = 0.5),
    .solver_spec(TSTAGEWISE_Solver, "T-Stagewise", lars),
    .solver_spec(TOMP_Solver, "T-OMP", omp),
    .solver_spec(TGP_Solver, "T-GP", omp),
    .solver_spec(TACGP_Solver, "T-ACGP", omp),
    .solver_spec(TMP_Solver, "T-MP", omp),
    .solver_spec(TNCGMP_Solver, "T-NCGMP", omp, ncgmp_variant = "line_search"),
    .solver_spec(TOOLS_Solver, "T-OLS", omp),
    .solver_spec(TAFS_Solver, "T-AFS (rho 0.1)", afs, rho_afs = 0.1),
    .solver_spec(TAFS_Solver, "T-AFS (rho 0.3)", afs, rho_afs = 0.3),
    .solver_spec(TAFS_Solver, "T-AFS (rho 0.5)", afs, rho_afs = 0.5),
    .solver_spec(TAFS_Solver, "T-AFS (rho 1)", afs, rho_afs = 1.0)
  )
}
