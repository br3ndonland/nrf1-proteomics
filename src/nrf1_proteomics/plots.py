# pyright: reportUnknownMemberType=false

from __future__ import annotations

import argparse
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias, cast

import numpy as np
import polars as pl
import scipy.stats
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.text import Annotation
from numpy.typing import NDArray

from nrf1_proteomics.analysis import (
    CONTRASTS,
    DEFAULT_RAW_DATA,
    FDR_THRESHOLD,
    GROUP_LABELS,
    GROUP_ORDER,
    PROJECT_ROOT,
    TREAT_LFC,
    AnalysisArtifacts,
    analyze_raw_data,
    calculate_model_residuals,
    read_raw_data,
)
from nrf1_proteomics.moderated import FloatArray

ObjectArray: TypeAlias = NDArray[np.object_]

DEFAULT_FIGURE_DIR = PROJECT_ROOT / "figures"
PVALUE_THRESHOLD = 0.05
VOLCANO_LABEL_LOG2_FOLD_CHANGE = 1.0
BACKGROUND = "not discovered (q >= 0.05)"
INCREASED = "increased (q < 0.05)"
DECREASED = "decreased (q < 0.05)"
CATEGORY_ORDER = (BACKGROUND, DECREASED, INCREASED)
CATEGORY_COLORS: Mapping[str, str] = {
    BACKGROUND: "#d1d5db",
    INCREASED: "#b91c1c",
    DECREASED: "#1d4ed8",
}
GROUP_COLORS: Mapping[str, str] = {
    "lacz": "#6b7280",
    "ha_chow": "#2563eb",
    "ha_chol": "#d97706",
    "ha_bort": "#b91c1c",
}
RUN_MARKERS: Mapping[str, str] = {"run_1": "o", "run_2": "s"}


@dataclass(frozen=True)
class VolcanoPlotSpec:
    key: str
    filename: str
    title: str
    xlabel: str


VOLCANO_PLOTS = (
    VolcanoPlotSpec(
        key="ha_chol_vs_ha_chow",
        filename="results-cholesterol",
        title="Nrf1 proteomics - cholesterol vs. chow",
        xlabel="Adjusted log2 fold change (cholesterol/chow)",
    ),
    VolcanoPlotSpec(
        key="ha_bort_vs_ha_chow",
        filename="results-bortezomib",
        title="Nrf1 proteomics - bortezomib vs. chow",
        xlabel="Adjusted log2 fold change (bortezomib/chow)",
    ),
    VolcanoPlotSpec(
        key="ha_chow_vs_lacz",
        filename="results-chow-vs-lacz",
        title="Nrf1 proteomics - chow vs. lacZ",
        xlabel="Adjusted log2 fold change (chow/lacZ)",
    ),
)


def prepare_volcano_data(
    analysis: pl.DataFrame,
    spec: VolcanoPlotSpec,
) -> pl.DataFrame:
    """Add volcano plotting fields without removing diagnostic warnings."""

    data = analysis.filter(pl.col("contrast") == spec.key)
    return (
        data.with_columns(
            negative_log10_p_value=pl.col("p_value").clip(1e-300, 1.0).log10().neg(),
            volcano_category=pl.when(
                pl.col("discovery") & (pl.col("adjusted_log2_fold_change") > 0)
            )
            .then(pl.lit(INCREASED))
            .when(pl.col("discovery"))
            .then(pl.lit(DECREASED))
            .otherwise(pl.lit(BACKGROUND)),
        )
        .with_columns(
            point_color=pl.col("volcano_category").replace_strict(CATEGORY_COLORS)
        )
        .sort("raw_row")
    )


def plot_volcano(
    analysis: pl.DataFrame,
    spec: VolcanoPlotSpec,
) -> Figure:
    """Create an FDR-coded volcano plot containing every primary protein."""

    plot_data = prepare_volcano_data(analysis, spec)
    fig, ax = plt.subplots(figsize=(7.6, 5.5), constrained_layout=True)
    for category in CATEGORY_ORDER:
        category_data = plot_data.filter(pl.col("volcano_category") == category)
        if category_data.is_empty():
            continue
        _ = ax.scatter(
            _float_column(category_data, "adjusted_log2_fold_change"),
            _float_column(category_data, "negative_log10_p_value"),
            s=30 if category == BACKGROUND else 42,
            color=CATEGORY_COLORS[category],
            alpha=0.65 if category == BACKGROUND else 0.88,
            edgecolors="none",
            label=category,
            zorder=2,
        )

    warnings = plot_data.filter(pl.col("diagnostic_warning"))
    if not warnings.is_empty():
        _ = ax.scatter(
            _float_column(warnings, "adjusted_log2_fold_change"),
            _float_column(warnings, "negative_log10_p_value"),
            s=48,
            facecolors="none",
            edgecolors="#6b7280",
            linewidths=0.75,
            alpha=0.9,
            label="diagnostic flag (outline)",
            zorder=3,
        )

    _style_volcano_axes(ax, spec)
    _annotate_significant_genes(ax, plot_data)
    return fig


def plot_qc_overview(
    artifacts: AnalysisArtifacts,
    raw_data: pl.DataFrame,
) -> Figure:
    """Plot channel distributions, sample similarity, PCA, and pooled-mix QC."""

    eligible_ids = artifacts.normalized_data.get_column("protein_id").to_list()
    eligible_raw = raw_data.filter(pl.col("protein_id").is_in(eligible_ids)).sort(
        "raw_row"
    )
    metadata = artifacts.sample_metadata
    sample_columns = _string_column(metadata, "sample")
    mix_columns = list(dict.fromkeys(_string_column(metadata, "mix_column")))
    raw_log2 = cast(
        FloatArray,
        np.log2(eligible_raw.select(sample_columns).to_numpy()),
    )
    normalized = cast(
        FloatArray,
        artifacts.normalized_data.select(sample_columns).to_numpy(),
    )
    fig, raw_axes = cast(
        tuple[Figure, object],
        plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True),
    )
    (
        raw_distribution_ax,
        normalized_distribution_ax,
        correlation_ax,
        pca_ax,
        mix_ax,
        channel_shift_ax,
    ) = _axes_tuple(raw_axes, 6)

    _distribution_panel(raw_distribution_ax, raw_log2, "Raw log2 reporter S/N")
    _distribution_panel(
        normalized_distribution_ax,
        normalized,
        "Median-centered log2 abundance",
    )

    correlation = np.corrcoef(normalized.T)
    image = correlation_ax.imshow(correlation, vmin=0, vmax=1, cmap="viridis")
    _ = correlation_ax.set_title("Sample correlation")
    _ = correlation_ax.set_xticks([])
    _ = correlation_ax.set_yticks([])
    _ = fig.colorbar(image, ax=correlation_ax, fraction=0.046, pad=0.04)

    _pca_panel(pca_ax, normalized, metadata)

    mix = cast(
        FloatArray,
        np.log2(eligible_raw.select(mix_columns).to_numpy()),
    )
    _ = mix_ax.scatter(
        mix[:, 0],
        mix[:, 1],
        s=18,
        color="#4b5563",
        alpha=0.65,
    )
    limits = (
        float(np.min(mix)) - 0.2,
        float(np.max(mix)) + 0.2,
    )
    _ = mix_ax.plot(
        limits,
        limits,
        color="#9ca3af",
        linestyle="--",
        linewidth=1,
    )
    _ = mix_ax.set_xlim(limits)
    _ = mix_ax.set_ylim(limits)
    correlation_mix = cast(
        float,
        scipy.stats.pearsonr(mix[:, 0], mix[:, 1]).statistic,
    )
    _ = mix_ax.set_title(f"Pooled mix cross-run QC (r = {correlation_mix:.2f})")
    _ = mix_ax.set_xlabel("Run 1 mix log2 S/N")
    _ = mix_ax.set_ylabel("Run 2 mix log2 S/N")

    channel_shifts = np.median(raw_log2, axis=0) - np.median(raw_log2)
    _ = channel_shift_ax.bar(
        np.arange(len(sample_columns)),
        channel_shifts,
        color=[GROUP_COLORS[group] for group in _string_column(metadata, "treatment")],
    )
    _ = channel_shift_ax.axhline(0, color="#4b5563", linewidth=1)
    _ = channel_shift_ax.set_title("Removed channel median shifts")
    _ = channel_shift_ax.set_xlabel("Biological reporter channel")
    _ = channel_shift_ax.set_ylabel("Log2 shift")
    _ = channel_shift_ax.set_xticks([])

    _ = fig.suptitle(
        "Nrf1 proteomics quality-control overview",
        fontweight="bold",
    )
    return fig


def plot_model_diagnostics(artifacts: AnalysisArtifacts) -> Figure:
    """Plot variance moderation, residual, influence, and p-value diagnostics."""

    diagnostics = artifacts.diagnostics
    _, standardized_residuals, _ = calculate_model_residuals(
        artifacts.normalized_data, artifacts.sample_metadata
    )
    fig, raw_axes = cast(
        tuple[Figure, object],
        plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True),
    )
    (
        variance_ax,
        qq_ax,
        cooks_ax,
        cholesterol_ax,
        bortezomib_ax,
        chow_ax,
    ) = _axes_tuple(raw_axes, 6)

    mean_abundance = _float_array(diagnostics, "mean_log2_abundance")
    raw_variance = _float_array(diagnostics, "raw_residual_variance")
    prior_variance = _float_array(diagnostics, "prior_variance")
    _ = variance_ax.scatter(
        mean_abundance,
        np.log10(raw_variance),
        s=17,
        color="#9ca3af",
        alpha=0.6,
        label="raw",
    )
    order = np.argsort(mean_abundance)
    _ = variance_ax.plot(
        mean_abundance[order],
        np.log10(prior_variance[order]),
        color="#b91c1c",
        linewidth=1.4,
        label="fitted prior",
    )
    _ = variance_ax.set_title("Abundance-variance trend")
    _ = variance_ax.set_xlabel("Mean normalized log2 abundance")
    _ = variance_ax.set_ylabel("log10 residual variance")
    _ = variance_ax.legend(frameon=False)

    ordered_residuals = np.sort(standardized_residuals.ravel())
    probabilities = (np.arange(ordered_residuals.size) + 0.5) / ordered_residuals.size
    expected = scipy.stats.norm.ppf(probabilities)
    _ = qq_ax.scatter(
        expected,
        ordered_residuals,
        s=8,
        color="#4b5563",
        alpha=0.45,
    )
    limits = [
        float(min(np.min(expected), np.min(ordered_residuals))),
        float(max(np.max(expected), np.max(ordered_residuals))),
    ]
    _ = qq_ax.plot(
        limits,
        limits,
        color="#b91c1c",
        linestyle="--",
        linewidth=1,
    )
    _ = qq_ax.set_title("Standardized residual Q-Q plot")
    _ = qq_ax.set_xlabel("Expected normal quantile")
    _ = qq_ax.set_ylabel("Observed residual quantile")

    max_cooks = _float_column(diagnostics, "max_cooks_distance")
    threshold = _float_column(diagnostics, "cooks_distance_threshold")[0]
    _ = cooks_ax.hist(max_cooks, bins=35, color="#6b7280", alpha=0.8)
    _ = cooks_ax.axvline(
        threshold,
        color="#b91c1c",
        linestyle="--",
        linewidth=1.3,
    )
    _ = cooks_ax.set_title("Maximum Cook distance per protein")
    _ = cooks_ax.set_xlabel("Maximum Cook distance")
    _ = cooks_ax.set_ylabel("Proteins")

    p_value_axes = (cholesterol_ax, bortezomib_ax, chow_ax)
    for ax, contrast in zip(p_value_axes, CONTRASTS, strict=True):
        contrast_data = artifacts.contrast_results.filter(
            pl.col("contrast") == contrast.key
        )
        p_values = _float_array(contrast_data, "p_value")
        _ = ax.hist(
            p_values,
            bins=[index / 20 for index in range(21)],
            color="#2563eb",
            alpha=0.75,
        )
        _ = ax.set_title(contrast.label)
        _ = ax.set_xlabel("Moderated raw p-value")
        _ = ax.set_ylabel("Proteins")
    _ = fig.suptitle("Blocked model diagnostics", fontweight="bold")
    return fig


def plot_c1q_observations(artifacts: AnalysisArtifacts) -> Figure:
    """Plot every normalized C1q mouse observation with run and group summaries."""

    metadata = artifacts.sample_metadata
    c1q_genes = cast(
        list[str],
        artifacts.normalized_data.filter(pl.col("gene").str.contains(r"(?i)^c1q[abc]$"))
        .get_column("gene")
        .unique()
        .sort()
        .to_list(),
    )
    if len(c1q_genes) != 3:
        raise ValueError("Expected one normalized row for each C1q subunit")
    fig, raw_axes = cast(
        tuple[Figure, object],
        plt.subplots(
            1,
            3,
            figsize=(14, 4.8),
            sharey=True,
            constrained_layout=True,
        ),
    )
    axes = _axes_tuple(raw_axes, 3)
    sample_columns = _string_column(metadata, "sample")
    group_positions = np.arange(len(GROUP_ORDER), dtype=float)
    offsets = np.linspace(-0.12, 0.12, 6)
    treatments = _string_column(metadata, "treatment")
    runs = _string_column(metadata, "run")
    for ax, gene in zip(axes, c1q_genes, strict=True):
        protein = artifacts.normalized_data.filter(pl.col("gene") == gene)
        if protein.height != 1:
            continue
        values = cast(
            FloatArray,
            np.asarray(
                [_float_column(protein, sample)[0] for sample in sample_columns],
                dtype=np.float64,
            ),
        )
        for group_index, group in enumerate(GROUP_ORDER):
            sample_indices = [
                index
                for index, treatment in enumerate(treatments)
                if treatment == group
            ]
            group_values = np.asarray(
                [_array_scalar_1d(values, index) for index in sample_indices],
                dtype=np.float64,
            )
            run_values = [runs[index] for index in sample_indices]
            for point_index, (value, run) in enumerate(
                zip(group_values, run_values, strict=True)
            ):
                _ = ax.scatter(
                    _array_scalar_1d(group_positions, group_index)
                    + _array_scalar_1d(offsets, point_index),
                    float(value),
                    marker=RUN_MARKERS[run],
                    s=42,
                    color=GROUP_COLORS[group],
                    edgecolors="white",
                    linewidths=0.6,
                    zorder=3,
                )
            mean = float(np.mean(group_values))
            critical = cast(
                float,
                scipy.stats.t.ppf(0.975, len(group_values) - 1),
            )
            half_interval = critical * cast(
                float,
                scipy.stats.sem(group_values),
            )
            _ = ax.errorbar(
                _array_scalar_1d(group_positions, group_index),
                mean,
                yerr=half_interval,
                color="#111827",
                marker="_",
                markersize=15,
                capsize=4,
                linewidth=1.4,
                zorder=4,
            )
        _ = ax.set_title(gene)
        _ = ax.set_xticks(
            group_positions,
            [GROUP_LABELS[group] for group in GROUP_ORDER],
        )
        ax.tick_params(axis="x", rotation=30)
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.7)
    _ = axes[0].set_ylabel("Median-centered log2 reporter S/N")
    legend = [
        Line2D(
            [0], [0], marker=marker, color="none", markerfacecolor="#6b7280", label=run
        )
        for run, marker in RUN_MARKERS.items()
    ]
    _ = axes[-1].legend(handles=legend, frameon=False, title="TMT run")
    _ = fig.suptitle(
        "C1q normalized mouse-level observations",
        fontweight="bold",
    )
    return fig


def plot_effect_sizes(
    artifacts: AnalysisArtifacts,
    proteins_per_contrast: int = 8,
) -> Figure:
    """Plot confidence intervals for proteins with the smallest q-values."""

    fig, raw_axes = cast(
        tuple[Figure, object],
        plt.subplots(1, 3, figsize=(16, 7), constrained_layout=True),
    )
    axes = _axes_tuple(raw_axes, 3)
    for ax, contrast in zip(axes, CONTRASTS, strict=True):
        data = (
            artifacts.contrast_results.filter(pl.col("contrast") == contrast.key)
            .sort("q_value", "p_value")
            .head(proteins_per_contrast)
            .sort("adjusted_log2_fold_change")
        )
        y = np.arange(data.height)
        effects = _float_array(data, "adjusted_log2_fold_change")
        lower = _float_array(data, "ci_lower_95")
        upper = _float_array(data, "ci_upper_95")
        colors = [
            CATEGORY_COLORS[INCREASED if value > 0 else DECREASED]
            for value in _float_column(data, "adjusted_log2_fold_change")
        ]
        _ = ax.errorbar(
            effects,
            y,
            xerr=np.vstack([effects - lower, upper - effects]),
            fmt="none",
            ecolor="#9ca3af",
            elinewidth=1.5,
            capsize=3,
            zorder=1,
        )
        _ = ax.scatter(effects, y, color=colors, s=44, zorder=2)
        labels: list[str] = []
        for row in data.iter_rows(named=True):
            result = cast(Mapping[str, object], row)
            name = str(result["gene"] or result["protein_id"])
            treat_status = "TREAT" if bool(result["treat_discovery"]) else "no TREAT"
            q_value = float(cast(float, result["q_value"]))
            labels.append(f"{name}  q={q_value:.3g}; {treat_status}")
        _ = ax.set_yticks(y, labels, fontsize=8)
        _ = ax.axvline(0, color="#4b5563", linewidth=1)
        _ = ax.axvline(
            -TREAT_LFC,
            color="#9ca3af",
            linestyle="--",
            linewidth=1,
        )
        _ = ax.axvline(
            TREAT_LFC,
            color="#9ca3af",
            linestyle="--",
            linewidth=1,
        )
        discoveries = data.filter(pl.col("discovery")).height
        _ = ax.set_title(
            f"{contrast.label}\n{discoveries} of top {data.height} have q < 0.05"
        )
        _ = ax.set_xlabel("Adjusted log2 fold change (95% CI)")
        ax.grid(axis="x", color="#e5e7eb", linewidth=0.7)
    _ = fig.suptitle(
        "Leading effects ranked by Benjamini-Hochberg q-value",
        fontweight="bold",
    )
    return fig


def plot_sensitivities(artifacts: AnalysisArtifacts) -> Figure:
    """Compare discovery counts and primary-effect concordance by scenario."""

    data = artifacts.sensitivities
    scenarios = cast(
        list[str],
        data.get_column("scenario").unique().sort().to_list(),
    )
    ordered_scenarios = ["primary"] + [
        scenario for scenario in scenarios if scenario != "primary"
    ]
    fig, raw_axes = cast(
        tuple[Figure, object],
        plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True),
    )
    discovery_ax, effect_ax = _axes_tuple(raw_axes, 2)

    width = 0.23
    x = np.arange(len(ordered_scenarios))
    for contrast_index, contrast in enumerate(CONTRASTS):
        counts = [
            data.filter(
                (pl.col("scenario") == scenario)
                & (pl.col("contrast") == contrast.key)
                & pl.col("discovery")
            ).height
            for scenario in ordered_scenarios
        ]
        _ = discovery_ax.bar(
            x + (contrast_index - 1) * width,
            counts,
            width,
            label=contrast.label,
        )
    _ = discovery_ax.set_xticks(
        x,
        [_scenario_label(name) for name in ordered_scenarios],
    )
    discovery_ax.tick_params(axis="x", rotation=30)
    _ = discovery_ax.set_ylabel("Benjamini-Hochberg discoveries")
    _ = discovery_ax.set_title("Discovery counts")
    _ = discovery_ax.legend(frameon=False, fontsize=8)

    alternatives = ordered_scenarios[1:]
    effect_change = np.zeros((len(alternatives), len(CONTRASTS)))
    for row_index, scenario in enumerate(alternatives):
        for column_index, contrast in enumerate(CONTRASTS):
            primary = data.filter(
                (pl.col("scenario") == "primary") & (pl.col("contrast") == contrast.key)
            ).select("protein_id", pl.col("adjusted_log2_fold_change").alias("primary"))
            alternative = data.filter(
                (pl.col("scenario") == scenario) & (pl.col("contrast") == contrast.key)
            ).select(
                "protein_id", pl.col("adjusted_log2_fold_change").alias("alternative")
            )
            paired = primary.join(alternative, on="protein_id", how="inner")
            effect_change[row_index, column_index] = cast(
                float,
                paired.select(
                    (pl.col("alternative") - pl.col("primary")).abs().median()
                ).item(),
            )
    image = effect_ax.imshow(effect_change, vmin=0, cmap="Blues", aspect="auto")
    _ = effect_ax.set_xticks(
        np.arange(len(CONTRASTS)),
        [contrast.label for contrast in CONTRASTS],
        rotation=25,
    )
    _ = effect_ax.set_yticks(
        np.arange(len(alternatives)),
        [_scenario_label(name) for name in alternatives],
    )
    _ = effect_ax.set_title("Median absolute effect change from primary")
    for row_index in range(len(alternatives)):
        for column_index in range(len(CONTRASTS)):
            _ = effect_ax.text(
                column_index,
                row_index,
                f"{_array_scalar_2d(effect_change, row_index, column_index):.3f}",
                ha="center",
                va="center",
                color="#111827",
            )
    _ = fig.colorbar(image, ax=effect_ax, fraction=0.046, pad=0.04)
    _ = fig.suptitle("Prespecified sensitivity analyses", fontweight="bold")
    return fig


def write_figures(
    artifacts: AnalysisArtifacts | None = None,
    output_dir: str | Path = DEFAULT_FIGURE_DIR,
    raw_path: str | Path = DEFAULT_RAW_DATA,
    formats: Sequence[str] = ("png", "svg"),
) -> dict[str, Path]:
    """Write all analysis figures in the requested publication formats."""

    invalid_formats = set(formats) - {"png", "svg"}
    if invalid_formats:
        raise ValueError(f"Unsupported figure format: {sorted(invalid_formats)}")
    analysis = analyze_raw_data(raw_path) if artifacts is None else artifacts
    raw_data = read_raw_data(raw_path)
    figures: dict[str, Figure] = {
        "qc-overview": plot_qc_overview(analysis, raw_data),
        "model-diagnostics": plot_model_diagnostics(analysis),
        "results-c1q": plot_c1q_observations(analysis),
        "results-effect-sizes": plot_effect_sizes(analysis),
        "results-sensitivity": plot_sensitivities(analysis),
    }
    for spec in VOLCANO_PLOTS:
        figures[spec.filename] = plot_volcano(analysis.contrast_results, spec)

    figure_dir = Path(output_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for name, figure in figures.items():
        for figure_format in formats:
            path = figure_dir / f"{name}.{figure_format}"
            figure.savefig(
                path,
                dpi=300 if figure_format == "png" else None,
                bbox_inches="tight",
                facecolor="white",
            )
            paths[f"{name}.{figure_format}"] = path
        plt.close(figure)
    return paths


def write_volcano_plots(
    analysis: pl.DataFrame | None = None,
    output_dir: str | Path = DEFAULT_FIGURE_DIR,
    raw_path: str | Path = DEFAULT_RAW_DATA,
    specs: Sequence[VolcanoPlotSpec] = VOLCANO_PLOTS,
    formats: Sequence[str] = ("png", "svg"),
) -> dict[str, Path]:
    """Write the three standard volcano plots without regenerating other figures."""

    contrast_results = (
        analyze_raw_data(raw_path).contrast_results if analysis is None else analysis
    )
    figure_dir = Path(output_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for spec in specs:
        figure = plot_volcano(contrast_results, spec)
        for figure_format in formats:
            path = figure_dir / f"{spec.filename}.{figure_format}"
            figure.savefig(
                path,
                dpi=300 if figure_format == "png" else None,
                bbox_inches="tight",
                facecolor="white",
            )
            paths[f"{spec.key}.{figure_format}"] = path
        plt.close(figure)
    return paths


def textxy(
    ax: Axes,
    x: Sequence[float],
    y: Sequence[float],
    labels: Sequence[str],
    *,
    origin: tuple[float, float] = (0.0, 0.0),
    offset_points: float = 2.0,
    font_size: float = 8.0,
) -> list[Annotation]:
    """Label points with an offset away from the supplied origin."""

    if not (len(x) == len(y) == len(labels)):
        raise ValueError("x, y, and labels must have the same length")
    origin_x, origin_y = origin
    annotations: list[Annotation] = []
    for x_value, y_value, label in zip(x, y, labels, strict=True):
        x_direction = 1.0 if x_value >= origin_x else -1.0
        y_direction = 1.0 if y_value >= origin_y else -1.0
        annotation = ax.annotate(
            label,
            xy=(x_value, y_value),
            xytext=(x_direction * offset_points, y_direction * offset_points),
            textcoords="offset points",
            horizontalalignment="left" if x_direction > 0 else "right",
            verticalalignment="bottom" if y_direction > 0 else "top",
            fontsize=font_size,
            color="#111827",
            zorder=4,
        )
        annotation.set_in_layout(False)
        annotations.append(annotation)
    return annotations


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    _ = parser.add_argument("--raw", type=Path, default=DEFAULT_RAW_DATA)
    _ = parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    _ = parser.add_argument(
        "--format", dest="formats", action="append", choices=("png", "svg")
    )
    args = parser.parse_args(argv)
    formats = cast(list[str] | None, args.formats) or ["png"]
    _ = write_figures(
        raw_path=cast(Path, args.raw),
        output_dir=cast(Path, args.output_dir),
        formats=formats,
    )
    return 0


def _distribution_panel(ax: Axes, values: FloatArray, title: str) -> None:
    channel_count = cast(int, values.shape[1])
    _ = ax.boxplot(
        [values[:, index] for index in range(channel_count)],
        showfliers=False,
        patch_artist=True,
        boxprops={"facecolor": "#dbeafe", "edgecolor": "#6b7280"},
        medianprops={"color": "#b91c1c"},
    )
    _ = ax.set_title(title)
    _ = ax.set_xlabel("Biological reporter channel")
    _ = ax.set_ylabel("Log2 S/N")
    _ = ax.set_xticks([])


def _pca_panel(ax: Axes, normalized: FloatArray, metadata: pl.DataFrame) -> None:
    sample_matrix = normalized.T
    sample_matrix = sample_matrix - np.mean(sample_matrix, axis=0)
    u, singular_values, _ = np.linalg.svd(sample_matrix, full_matrices=False)
    scores = u[:, :2] * singular_values[:2]
    explained = singular_values**2 / np.sum(singular_values**2)
    treatments = _string_column(metadata, "treatment")
    runs = _string_column(metadata, "run")
    for treatment in GROUP_ORDER:
        for run, marker in RUN_MARKERS.items():
            sample_indices = [
                index
                for index, (sample_treatment, sample_run) in enumerate(
                    zip(treatments, runs, strict=True)
                )
                if sample_treatment == treatment and sample_run == run
            ]
            _ = ax.scatter(
                [_array_scalar_2d(scores, index, 0) for index in sample_indices],
                [_array_scalar_2d(scores, index, 1) for index in sample_indices],
                color=GROUP_COLORS[treatment],
                marker=marker,
                s=62,
                edgecolors="white",
                linewidths=0.7,
            )
    _ = ax.set_title("PCA of normalized samples")
    _ = ax.set_xlabel(f"PC1 ({_array_scalar_1d(explained, 0) * 100:.1f}%)")
    _ = ax.set_ylabel(f"PC2 ({_array_scalar_1d(explained, 1) * 100:.1f}%)")
    treatment_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=GROUP_COLORS[group],
            label=GROUP_LABELS[group],
        )
        for group in GROUP_ORDER
    ]
    run_handles = [
        Line2D(
            [0],
            [0],
            marker=marker,
            color="none",
            markerfacecolor="#6b7280",
            label=run,
        )
        for run, marker in RUN_MARKERS.items()
    ]
    _ = ax.legend(
        handles=[*treatment_handles, *run_handles],
        frameon=False,
        fontsize=7,
    )


def _style_volcano_axes(ax: Axes, spec: VolcanoPlotSpec) -> None:
    _ = ax.axvline(-TREAT_LFC, color="#6b7280", linewidth=1, linestyle="--")
    _ = ax.axvline(TREAT_LFC, color="#6b7280", linewidth=1, linestyle="--")
    _ = ax.axhline(
        -math.log10(PVALUE_THRESHOLD),
        color="#9ca3af",
        linewidth=1,
        linestyle=":",
    )
    _ = ax.set_title(spec.title, fontweight="bold")
    _ = ax.set_xlabel(spec.xlabel)
    _ = ax.set_ylabel("-log10(raw moderated p-value)")
    ax.grid(True, color="#e5e7eb", linewidth=0.6, alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    _ = ax.legend(frameon=False, loc="best", fontsize=8)


def _annotate_significant_genes(
    ax: Axes,
    plot_data: pl.DataFrame,
) -> None:
    label_data = plot_data.filter(
        (pl.col("q_value") < FDR_THRESHOLD)
        & (pl.col("adjusted_log2_fold_change").abs() > VOLCANO_LABEL_LOG2_FOLD_CHANGE)
        & pl.col("gene").is_not_null()
        & (pl.col("gene") != "")
    )
    if label_data.is_empty():
        return
    font_size = max(6.0, 9.0 - label_data.height / 10)
    _ = textxy(
        ax,
        _float_column(label_data, "adjusted_log2_fold_change"),
        _float_column(label_data, "negative_log10_p_value"),
        _string_column(label_data, "gene"),
        offset_points=4,
        font_size=font_size,
    )


def _float_column(data: pl.DataFrame, column: str) -> list[float]:
    values = cast(list[object], data.get_column(column).to_list())
    return [float(cast(float, value)) for value in values]


def _float_array(data: pl.DataFrame, column: str) -> FloatArray:
    return np.asarray(_float_column(data, column), dtype=np.float64)


def _string_column(data: pl.DataFrame, column: str) -> list[str]:
    values = cast(list[object], data.get_column(column).to_list())
    return [str(value) for value in values]


def _scenario_label(name: str) -> str:
    return name.replace("_", " ")


def _axes_tuple(value: object, expected: int) -> tuple[Axes, ...]:
    if not isinstance(value, np.ndarray):
        raise TypeError("Expected Matplotlib to return an axes array")
    axes = cast(ObjectArray, value).ravel()
    count = axes.size
    if count != expected:
        raise ValueError(f"Expected {expected} axes, found {count}")
    return tuple(cast(Axes, axes[index]) for index in range(count))


def _array_scalar_1d(values: FloatArray, index: int) -> float:
    return float(cast(np.float64, values[index]))


def _array_scalar_2d(values: FloatArray, row: int, column: int) -> float:
    return float(cast(np.float64, values[row, column]))


if __name__ == "__main__":
    raise SystemExit(main())
