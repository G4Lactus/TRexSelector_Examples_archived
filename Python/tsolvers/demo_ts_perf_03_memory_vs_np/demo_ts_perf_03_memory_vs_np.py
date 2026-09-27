# ==============================================================================
# demo_ts_perf_03_memory_vs_np.py
# ==============================================================================
#
# Solver memory vs. the ratio n/p, for all demo T-solvers.
#
# Mirrors cpp/tsolvers/demo_ts_perf_03_memory_vs_np/demo_ts_perf_03_memory_vs_np.cpp:
# p = 1,000, L = p, n in {50, 100, 200, 500, 1,000, 2,000}, 10 true
# coefficients equal to 3 on a fixed random support, noise sigma about 1,
# T = 10. Per run and n the parent draws and scales the data once; every solver
# then runs in a forked child as a precentered reader. The child first reads
# every page of X, D and y and solves a small warm-up problem, so its baseline
# counts the data, the solver code and the interpreter, and reports its peak
# resident set minus that baseline as the solver memory. In Python this also
# includes what the binding allocates around the C++ solver.
#
# Usage: demo_ts_perf_03_memory_vs_np.py [num_MC] [parallel children]   (default 500, 6)
#
# ==============================================================================

import os
import struct
import sys
import time

import numpy as np

from trex_selector_neo.utils import set_num_threads

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_THIS_DIR)
if _PARENT_DIR not in sys.path:
    sys.path.insert(0, _PARENT_DIR)

from ts_demo_utils import display_path, output_dir, print_section_header
from ts_perf_solvers import demo_solvers
from ts_perf_utils import (
    draw_support,
    parse_positive,
    peak_rss_mib,
    reader_data,
)

# Setup
p = 1_000
num_dummies = p
num_active = 10
beta_active = 3.0
snr = 90.0                    # signal variance 10 * 3^2 over noise variance 1
T_stop = 10
support_seed = 12345
data_seed_offset = 1000       # run mc uses seed data_seed_offset + mc
n_grid = [50, 100, 200, 500, 1_000, 2_000]
REPORT = struct.Struct("ddi")  # solver memory (MiB), steps, ok

SUPPORT = draw_support(p, num_active, support_seed)
COEFS = [beta_active] * num_active


def warm_up(spec):
    """Solve a small problem, so the baseline holds the solver's code pages and a warm allocator."""
    X, D, y = reader_data(20, 30, 30, [0, 1], [1.0, 1.0], 10.0, seed=1)
    solver = spec.make(X, D, y)
    solver.setTieSeed(1)
    solver.executeStep(2, True)


def launch(spec, s, X, D, y, seed):
    """Run one solver in a forked child and return (pid, read end of its pipe, solver index)."""
    r, w = os.pipe()
    pid = os.fork()
    if pid == 0:
        os.close(r)
        memory, steps, ok = 0.0, 0.0, 0
        _ = float(X.sum() + D.sum() + y.sum())  # touch every data page
        try:
            warm_up(spec)
        except (RuntimeError, ValueError):
            pass  # the measured run below reports a failing solver
        baseline = peak_rss_mib()
        try:
            solver = spec.make(X, D, y)
            solver.setTieSeed(seed)
            solver.executeStep(T_stop, True)
            steps = float(solver.getNumSteps())
            memory = peak_rss_mib() - baseline
            ok = 1
        except (RuntimeError, ValueError):
            ok = 0
        os.write(w, REPORT.pack(memory, steps, ok))
        os._exit(0)
    os.close(w)
    return pid, r, s


def collect(running):
    """Wait for one of the children and read its report: (solver index, memory, steps, ok)."""
    pid, _ = os.wait()
    k = next(i for i, c in enumerate(running) if c[0] == pid)
    _, fd, s = running.pop(k)
    data = os.read(fd, REPORT.size)
    os.close(fd)
    if len(data) != REPORT.size:
        return s, 0.0, 0.0, 0
    return (s, *REPORT.unpack(data))


def main():
    num_MC = parse_positive(sys.argv, 1, 500)
    parallel = parse_positive(sys.argv, 2, 6)
    set_num_threads(1)  # single-threaded parent and solvers, so fork() stays safe with OpenMP
    solvers = demo_solvers()
    S, K = len(solvers), len(n_grid)

    print()
    print_section_header("T-Solver Performance: Memory vs. n/p")
    print(f"p = {p}, L = p, n in {{50, 100, 200, 500, 1000, 2000}}, T = {T_stop}, true coefficients: "
          f"{num_active} x {beta_active:g}, SNR = {snr:g}")
    print(f"{S} solvers, {num_MC} runs, {parallel} children at a time (one forked child per solver)\n")

    memory = [[[] for _ in range(K)] for _ in range(S)]
    steps = [[[] for _ in range(K)] for _ in range(S)]
    failures = 0
    t_start = time.perf_counter()
    for mc in range(num_MC):
        for k, n in enumerate(n_grid):
            X, D, y = reader_data(n, p, num_dummies, SUPPORT, COEFS, snr, seed=data_seed_offset + mc)
            running = []

            def store(result, k=k):
                nonlocal failures
                s, mem, st, ok = result
                if ok:
                    memory[s][k].append(mem)
                    steps[s][k].append(st)
                else:
                    failures += 1

            for s, spec in enumerate(solvers):
                if len(running) >= parallel:
                    store(collect(running))
                running.append(launch(spec, s, X, D, y, data_seed_offset + mc))
            while running:
                store(collect(running))
        if (mc + 1) % 50 == 0 or mc + 1 == num_MC:
            print(f"  runs done: {mc + 1} of {num_MC}")
    wall = time.perf_counter() - t_start

    # Data size per n: X, D and y as doubles
    data_mib = [8.0 * n * (p + num_dummies + 1) / (1024.0 * 1024.0) for n in n_grid]

    def header(title):
        return [title, f"{'solver':<22}" + "".join(f"{'n=' + str(n):>9}" for n in n_grid), "-" * (22 + 9 * K)]

    out = [""] + header("Data X, D and y (MiB)")
    out.append(f"{'all solvers':<22}" + "".join(f"{v:>9.2f}" for v in data_mib))
    out.append("")

    def table(title, value, precision):
        out.extend(header(title))
        for s, spec in enumerate(solvers):
            out.append(f"{spec.name:<22}" + "".join(f"{value(s, k):>9.{precision}f}" for k in range(K)))
        out.append("")

    def largest(s, k):
        return max(memory[s][k]) if memory[s][k] else 0.0

    table("Median solver memory (MiB)", lambda s, k: float(np.median(memory[s][k])), 2)
    table("Largest solver memory over the runs (MiB)", largest, 2)
    table("Median solver memory relative to the data (%)",
          lambda s, k: 100.0 * float(np.median(memory[s][k])) / data_mib[k], 1)
    table("Median steps to T = 10", lambda s, k: float(np.median(steps[s][k])), 1)
    out.append(f"Failed solver runs: {failures}")
    out.append(f"Wall time: {wall:.1f} s")
    text = "\n".join(out) + "\n"
    print(text, end="")

    stem = f"memory_vs_np_mc{num_MC}"
    data_dir = output_dir(_THIS_DIR)
    csv_path = os.path.join(data_dir, stem + ".csv")
    with open(csv_path, "w") as f:
        f.write("solver,family,n,n_over_p,data_mib,solver_mib,solver_max_mib,steps\n")
        for s, spec in enumerate(solvers):
            for k, n in enumerate(n_grid):
                row = [n / p, data_mib[k], float(np.median(memory[s][k])), largest(s, k), float(np.median(steps[s][k]))]
                f.write(f"{spec.name},{spec.family},{n}," + ",".join(f"{v:.6g}" for v in row) + "\n")
    txt_path = os.path.join(data_dir, stem + ".txt")
    with open(txt_path, "w") as f:
        f.write(text)

    print(f"\nResults written to:\n  {display_path(csv_path)}\n  {display_path(txt_path)}\n")


if __name__ == "__main__":
    main()
