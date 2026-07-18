from __future__ import annotations

import argparse
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import polars as pl
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from nrf1_proteomics._adjust_text import adjust_text_labels
from nrf1_proteomics.analysis import DEFAULT_RAW_DATA, PROJECT_ROOT, analyze_raw_data

DEFAULT_FIGURE_DIR = PROJECT_ROOT / "figures"
FOLD_CHANGE_THRESHOLD = 1.0
PVALUE_THRESHOLD = 0.05

BACKGROUND = "not significant"
PVALUE_ONLY = "p < 0.05"
FOLD_CHANGE_ONLY = "abs(log2 FC) > 1"
PVALUE_AND_FOLD_CHANGE = "p < 0.05 and abs(log2 FC) > 1"
CATEGORY_ORDER = (
    BACKGROUND,
    PVALUE_ONLY,
    FOLD_CHANGE_ONLY,
    PVALUE_AND_FOLD_CHANGE,
)
CATEGORY_STYLES: Mapping[str, Mapping[str, object]] = {
    BACKGROUND: {"size": 24, "alpha": 0.58},
    PVALUE_ONLY: {"size": 30, "alpha": 0.82},
    FOLD_CHANGE_ONLY: {"size": 30, "alpha": 0.82},
    PVALUE_AND_FOLD_CHANGE: {"size": 44, "alpha": 0.9},
}


@dataclass(frozen=True)
class VolcanoPlotSpec:
    key: str
    filename: str
    x_column: str
    pvalue_column: str
    title: str
    xlabel: str


VOLCANO_PLOTS = (
    VolcanoPlotSpec(
        key="ha_chol",
        filename="results-hachol.png",
        x_column="log2_FC_HAchol_HAchow",
        pvalue_column="pvalue",
        title="Nrf1 AP-MS: HA cholesterol vs HA chow",
        xlabel="log2(FC HA chol/HA chow)",
    ),
    VolcanoPlotSpec(
        key="ha_bort",
        filename="results-habort.png",
        x_column="log2_FC_HAbort_HAchow",
        pvalue_column="pvaluebort",
        title="Nrf1 AP-MS: HA bortezomib vs HA chow",
        xlabel="log2(FC HA bort/HA chow)",
    ),
    VolcanoPlotSpec(
        key="ha_chow_lacz",
        filename="results-hachow-lacz.png",
        x_column="log2_FC_HAchow_lacZ",
        pvalue_column="pvaluechow",
        title="Nrf1 AP-MS: HA chow vs lacZ",
        xlabel="log2(FC HA chow/lacZ)",
    ),
)


def prepare_volcano_data(
    analysis: pl.DataFrame,
    spec: VolcanoPlotSpec,
    fold_change_threshold: float = FOLD_CHANGE_THRESHOLD,
    pvalue_threshold: float = PVALUE_THRESHOLD,
) -> pl.DataFrame:
    """Add plotting columns for one volcano plot from an analyzed DataFrame."""

    fold_change_significant = pl.col("log2_fold_change").abs() > fold_change_threshold
    pvalue_significant = pl.col("pvalue") < pvalue_threshold

    return (
        analysis.select(
            "Gene",
            pl.col(spec.x_column).alias("log2_fold_change"),
            pl.col(spec.pvalue_column).alias("pvalue"),
        )
        .with_columns(
            negative_log10_pvalue=pl.col("pvalue").log(base=10).neg(),
            fold_change_significant=fold_change_significant,
            pvalue_significant=pvalue_significant,
        )
        .with_columns(
            volcano_category=pl.when(
                pl.col("fold_change_significant") & pl.col("pvalue_significant")
            )
            .then(pl.lit(PVALUE_AND_FOLD_CHANGE))
            .when(pl.col("pvalue_significant"))
            .then(pl.lit(PVALUE_ONLY))
            .when(pl.col("fold_change_significant"))
            .then(pl.lit(FOLD_CHANGE_ONLY))
            .otherwise(pl.lit(BACKGROUND)),
            point_color=pl.when(
                pl.col("fold_change_significant") & pl.col("pvalue_significant")
            )
            .then(pl.lit("green"))
            .when(pl.col("pvalue_significant"))
            .then(pl.lit("red"))
            .when(pl.col("fold_change_significant"))
            .then(pl.lit("orange"))
            .otherwise(pl.lit("black")),
        )
        .select(
            "Gene",
            "log2_fold_change",
            "pvalue",
            "negative_log10_pvalue",
            "volcano_category",
            "point_color",
        )
    )


def plot_volcano(
    analysis: pl.DataFrame,
    spec: VolcanoPlotSpec,
    label_genes: Sequence[str] | None = None,
    fold_change_threshold: float = FOLD_CHANGE_THRESHOLD,
    pvalue_threshold: float = PVALUE_THRESHOLD,
) -> Figure:
    """Create a volcano plot from an analyzed Polars DataFrame."""

    plot_data = prepare_volcano_data(
        analysis=analysis,
        spec=spec,
        fold_change_threshold=fold_change_threshold,
        pvalue_threshold=pvalue_threshold,
    )
    fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=150, constrained_layout=True)

    for category in CATEGORY_ORDER:
        category_data = plot_data.filter(pl.col("volcano_category") == category)
        if category_data.is_empty():
            continue
        style = CATEGORY_STYLES[category]
        _ = ax.scatter(
            _float_column(category_data, "log2_fold_change"),
            _float_column(category_data, "negative_log10_pvalue"),
            s=cast(float, style["size"]),
            color=_string_column(category_data, "point_color"),
            alpha=cast(float, style["alpha"]),
            edgecolors="none",
            label=category,
        )

    _style_volcano_axes(
        ax=ax,
        spec=spec,
        fold_change_threshold=fold_change_threshold,
        pvalue_threshold=pvalue_threshold,
    )
    _annotate_significant_genes(ax, plot_data, label_genes)
    return fig


def write_volcano_plots(
    analysis: pl.DataFrame | None = None,
    output_dir: str | Path = DEFAULT_FIGURE_DIR,
    raw_path: str | Path = DEFAULT_RAW_DATA,
    specs: Sequence[VolcanoPlotSpec] = VOLCANO_PLOTS,
) -> dict[str, Path]:
    """Write all standard volcano plots and return their output paths."""

    data = analyze_raw_data(raw_path) if analysis is None else analysis
    figure_dir = Path(output_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)

    output_paths: dict[str, Path] = {}
    for spec in specs:
        fig = plot_volcano(data, spec)
        output_path = figure_dir / spec.filename
        fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        output_paths[spec.key] = output_path

    return output_paths


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    _ = parser.add_argument("--raw", type=Path, default=DEFAULT_RAW_DATA)
    _ = parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    args = parser.parse_args(argv)
    raw_path = cast(Path, args.raw)
    output_dir = cast(Path, args.output_dir)

    _ = write_volcano_plots(raw_path=raw_path, output_dir=output_dir)
    return 0


def _style_volcano_axes(
    ax: Axes,
    spec: VolcanoPlotSpec,
    fold_change_threshold: float,
    pvalue_threshold: float,
) -> None:
    _ = ax.axvline(-fold_change_threshold, color="#6f7378", linewidth=1, linestyle="--")
    _ = ax.axvline(fold_change_threshold, color="#6f7378", linewidth=1, linestyle="--")
    _ = ax.axhline(
        -math.log10(pvalue_threshold),
        color="#6f7378",
        linewidth=1,
        linestyle=":",
    )
    _ = ax.set_title(spec.title, fontweight="bold")
    _ = ax.set_xlabel(spec.xlabel)
    _ = ax.set_ylabel("-log10(pvalue)")
    ax.grid(True, color="#d7dce1", linewidth=0.6, alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _ = ax.legend(frameon=False, loc="best", fontsize=8)


def _annotate_significant_genes(
    ax: Axes,
    plot_data: pl.DataFrame,
    label_genes: Sequence[str] | None = None,
) -> None:
    label_data = plot_data.filter(
        (pl.col("volcano_category") == PVALUE_AND_FOLD_CHANGE)
        & pl.col("Gene").is_not_null()
        & (pl.col("Gene") != "")
    )
    if label_genes is not None:
        label_data = label_data.filter(pl.col("Gene").is_in(label_genes))
    if label_data.is_empty():
        return

    label_x = _float_column(label_data, "log2_fold_change")
    label_y = _float_column(label_data, "negative_log10_pvalue")
    font_size = max(6.0, 10.0 - label_data.height / 10)
    labels = [
        ax.text(
            x,
            y,
            gene,
            fontsize=font_size,
            fontweight="bold",
            color="#111827",
            zorder=4,
        )
        for x, y, gene in zip(label_x, label_y, _string_column(label_data, "Gene"))
    ]
    adjust_text_labels(
        labels,
        x=_float_column(plot_data, "log2_fold_change"),
        y=_float_column(plot_data, "negative_log10_pvalue"),
        objects=[legend] if (legend := ax.get_legend()) is not None else None,
        ax=ax,
        iter_lim=500,
        min_arrow_len=4,
        arrowprops={"arrowstyle": "-", "color": "#111827", "linewidth": 0.8},
    )


def _float_column(data: pl.DataFrame, column: str) -> list[float]:
    values = cast(list[object], data.get_column(column).to_list())
    return [cast(float, value) for value in values]


def _string_column(data: pl.DataFrame, column: str) -> list[str]:
    values = cast(list[object], data.get_column(column).to_list())
    return [cast(str, value) for value in values]
