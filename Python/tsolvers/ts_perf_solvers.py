# ==============================================================================
# ts_perf_solvers.py
# ==============================================================================
#
# The T-solvers compared by the performance demos (demo_ts_perf_*), each with
# the settings of its own solver demo. Mirrors cpp/tsolvers/ts_perf_solvers.hpp.
#
# Every solver runs as a reader (normalize=False, intercept=False): the caller
# centers and scales X and D and centers y once in front, and all solvers share
# the Fortran-ordered arrays without copies. The data are marked precentered,
# so the row-bounded active sets stop at the rank n - 1 of the centered design
# (T-ENET with lambda2 > 0 has no such cap). T-AFS appears once per rho.
#
# ==============================================================================

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from trex_selector_neo.tsolvers.ad import SolverOptions
from trex_selector_neo.tsolvers.ad.afs_based import TAFS_Solver
from trex_selector_neo.tsolvers.ad.lars_based import (
    TENET_Solver,
    TLARS_Solver,
    TLASSO_Solver,
    TSTAGEWISE_Solver,
    TSTEPWISE_Solver,
)
from trex_selector_neo.tsolvers.ad.omp_based import (
    NCGMPVariant,
    TACGP_Solver,
    TGP_Solver,
    TMP_Solver,
    TNCGMP_Solver,
    TOMP_Solver,
    TOOLS_Solver,
)


@dataclass(frozen=True)
class SolverSpec:
    """Name, family and factory of one compared solver."""

    name: str                               # label in tables and figures, e.g. "T-LARS"
    family: str                             # "LARS-based", "OMP-based" or "AFS-based"
    make: Callable[[Any, Any, Any], Any]    # builds the solver as a reader on (X, D, y)


def _spec(cls, name, family, **extra):
    """Spec for solver class `cls` with the reader options plus `extra`."""
    def make(X, D, y):
        return cls(X, D, y, SolverOptions(normalize=False, intercept=False, precentered=True, **extra))
    return SolverSpec(name, family, make)


def demo_solvers():
    """The 12 demo solvers with their demo settings, T-AFS for rho in {0.1, 0.3, 0.5, 1}."""
    lars, omp, afs = "LARS-based", "OMP-based", "AFS-based"
    return [
        _spec(TLARS_Solver, "T-LARS", lars),
        _spec(TLASSO_Solver, "T-LASSO", lars),
        _spec(TSTEPWISE_Solver, "T-Stepwise", lars),
        _spec(TENET_Solver, "T-ENET (lambda2 0.5)", lars, lambda2=0.5),
        _spec(TSTAGEWISE_Solver, "T-Stagewise", lars),
        _spec(TOMP_Solver, "T-OMP", omp),
        _spec(TGP_Solver, "T-GP", omp),
        _spec(TACGP_Solver, "T-ACGP", omp),
        _spec(TMP_Solver, "T-MP", omp),
        _spec(TNCGMP_Solver, "T-NCGMP", omp, ncgmp_variant=NCGMPVariant.LineSearch),
        _spec(TOOLS_Solver, "T-OLS", omp),
        _spec(TAFS_Solver, "T-AFS (rho 0.1)", afs, rho_afs=0.1),
        _spec(TAFS_Solver, "T-AFS (rho 0.3)", afs, rho_afs=0.3),
        _spec(TAFS_Solver, "T-AFS (rho 0.5)", afs, rho_afs=0.5),
        _spec(TAFS_Solver, "T-AFS (rho 1)", afs, rho_afs=1.0),
    ]
