#!/usr/bin/env python3
"""
ts_plt_utils.py — shared plotting utilities for the T-solver demo suite.

Suite-level companion of ts_perf_utils.hpp. Each demo's generate_plots.sh
supplies its CSV files and titles, so no per-demo configuration is needed.

Plot patterns hosted here:
    plot_terminated_paths()  — coefficient paths vs. number of included dummies,
                               one panel per T_stop of a warm-started solver
                               (the figure of the R tlars vignette).

CLI mode (typically driven from a demo's generate_plots.sh):
    python ts_plt_utils.py paths <stem>_path.csv --solver T-LARS --title "..."

The path mode reads <stem>_path.csv (T_stop,step,index,is_dummy,beta) and the
sibling <stem>_actions.csv (T_stop,step,index,is_dummy,added) written by
ts_perf::PathRecorder. Figures default to the simulation_results/plots/
sibling of the data folder; override with --outdir.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless / batch-safe: no display required

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

# Colours and line styles of the tlars vignette path plot
COLOR_ORIGINAL = "black"
COLOR_DUMMY = "#d62728"
COLOR_STEP = "0.65"
MAX_TOP_LABELS = 40  # beyond this, the top axis labels get small and rotated


# ==============================================================================
# Shared helpers
# ==============================================================================

def default_plots_dir(csv_path: Path) -> Path:
    """simulation_results/plots/ next to simulation_results/data/, else the CSV's folder."""
    parent = csv_path.resolve().parent
    if parent.name == "data" and parent.parent.name == "simulation_results":
        return parent.parent / "plots"
    return parent


def save_figure(fig: plt.Figure, outdir: Path, stem: str, formats: list[str]) -> None:
    """Write the figure to outdir/<stem>.<ext> for each requested format."""
    outdir.mkdir(parents=True, exist_ok=True)
    for extension in formats:
        output_path = outdir / f"{stem}.{extension}"
        fig.savefig(
            output_path,
            dpi=300 if extension == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )
        print(f"Wrote: {output_path}")


# ==============================================================================
# Terminated solution paths
# ==============================================================================

def read_paths(path_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read <stem>_path.csv and its sibling <stem>_actions.csv."""
    actions_csv = path_csv.with_name(path_csv.name.replace("_path.csv", "_actions.csv"))
    if actions_csv == path_csv or not actions_csv.exists():
        raise FileNotFoundError(f"no actions file next to {path_csv}")
    return pd.read_csv(path_csv), pd.read_csv(actions_csv)


def _draw_path_panel(ax: plt.Axes, path: pd.DataFrame, actions: pd.DataFrame, T_stop: int) -> None:
    """One panel: coefficient paths up to the T_stop-th included dummy."""
    num_steps = int(actions["step"].max())
    steps = range(num_steps + 1)

    # Every variable that is ever nonzero, as a dense line over all steps
    wide = path.pivot_table(index="step", columns="index", values="beta", aggfunc="first")
    wide = wide.reindex(steps).fillna(0.0)
    dummy_of = path.groupby("index")["is_dummy"].first()

    dummy_steps = actions.loc[(actions["is_dummy"] == 1) & (actions["added"] == 1), "step"].tolist()

    for step in range(1, num_steps + 1):
        ax.axvline(step, color=COLOR_STEP, linestyle=(0, (5, 2, 1, 2)), linewidth=0.6, zorder=0)
    for step in dummy_steps:
        ax.axvline(step, color=COLOR_DUMMY, linewidth=1.3, zorder=1)
    ax.axhline(0.0, color=COLOR_ORIGINAL, linewidth=0.8, zorder=1)

    for index in wide.columns:
        is_dummy = bool(dummy_of[index])
        ax.plot(
            wide.index,
            wide[index],
            color=COLOR_DUMMY if is_dummy else COLOR_ORIGINAL,
            linestyle="--" if is_dummy else "-",
            linewidth=1.2,
            zorder=3 if is_dummy else 2,
        )

    # Margin on both sides, so the curves visibly start from 0 at step 0
    margin = max(0.5, 0.02 * num_steps)
    ax.set_xlim(-margin, num_steps + margin)
    ax.set_xticks(dummy_steps)
    ax.set_xticklabels([str(k) for k in range(1, len(dummy_steps) + 1)])
    ax.set_xlabel("Number of included dummies")
    ax.set_ylabel("Coefficients")
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    # Top axis: the variable (index from the CSV) or dummy (D) entering at each step
    added = actions[actions["added"] == 1]
    top = ax.secondary_xaxis("top")
    top.set_xticks(added["step"].tolist())
    labels = ["D" if d else str(i) for i, d in zip(added["index"], added["is_dummy"])]
    if len(labels) > MAX_TOP_LABELS:
        top.set_xticklabels(labels, fontsize=6, rotation=90)
    else:
        top.set_xticklabels(labels, fontsize=8)
    top.set_xlabel("Index of the entering variable (0-based, D = dummy)", fontsize=9)

    num_dummies = len(dummy_steps)
    dummy_word = "dummy" if num_dummies == 1 else "dummies"
    ax.set_title(f"$T = {T_stop}$: {num_steps} steps, {num_dummies} {dummy_word} included",
                 fontsize=11, loc="left")


def plot_terminated_paths(path: pd.DataFrame, actions: pd.DataFrame, solver: str,
                          title: str | None = None) -> plt.Figure:
    """Coefficient paths vs. number of included dummies, one panel per T_stop."""
    T_values = sorted(actions["T_stop"].unique())
    fig, axes = plt.subplots(len(T_values), 1, figsize=(10, 4.2 * len(T_values)), squeeze=False)

    for ax, T_stop in zip(axes[:, 0], T_values):
        _draw_path_panel(ax, path[path["T_stop"] == T_stop], actions[actions["T_stop"] == T_stop],
                         int(T_stop))

    handles = [
        Line2D([], [], color=COLOR_ORIGINAL, linestyle="-", label="Original variables"),
        Line2D([], [], color=COLOR_DUMMY, linestyle="--", label="Dummies"),
    ]
    axes[0, 0].legend(handles=handles, loc="upper left", frameon=True, facecolor="white",
                      edgecolor="0.7")
    fig.suptitle(title or f"{solver}: terminated solution paths", fontsize=13)
    fig.tight_layout()
    return fig


# ==============================================================================
# CLI
# ==============================================================================

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)

    paths = sub.add_parser("paths", help="terminated solution paths of one solver")
    paths.add_argument("path_csv", type=Path, help="<stem>_path.csv written by PathRecorder")
    paths.add_argument("--solver", required=True, help="solver name for the default title")
    paths.add_argument("--title", default=None, help="figure title")
    paths.add_argument("--stem", default=None, help="output stem (default: CSV stem without _path)")
    paths.add_argument("--outdir", type=Path, default=None, help="output folder")
    paths.add_argument("--formats", nargs="+", default=["png", "pdf"], help="output formats")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    if args.mode == "paths":
        path, actions = read_paths(args.path_csv)
        fig = plot_terminated_paths(path, actions, args.solver, args.title)
        stem = args.stem or args.path_csv.stem.removesuffix("_path")
        save_figure(fig, args.outdir or default_plots_dir(args.path_csv), stem, args.formats)
        plt.close(fig)


if __name__ == "__main__":
    main()
