from functools import lru_cache
from pathlib import Path

import polars as pl
from matplotlib import pyplot as plt

import nrf1_proteomics.plots as plots
from nrf1_proteomics.analysis import AnalysisArtifacts, analyze_raw_data

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DATA = PROJECT_ROOT / "data" / "nrf1-proteomics-raw.csv"


@lru_cache(maxsize=1)
def _artifacts() -> AnalysisArtifacts:
    return analyze_raw_data(RAW_DATA)


def test_prepare_volcano_data_assigns_bh_direction_colors() -> None:
    spec = plots.VOLCANO_PLOTS[0]
    analysis = pl.DataFrame(
        {
            "raw_row": [1, 2, 3],
            "gene": ["background", "up", "down"],
            "contrast": [spec.key] * 3,
            "adjusted_log2_fold_change": [0.5, 1.5, -1.5],
            "p_value": [0.01, 0.01, 0.01],
            "discovery": [False, True, True],
            "diagnostic_warning": [False, False, True],
        }
    )

    plot_data = plots.prepare_volcano_data(analysis, spec)

    assert plot_data.select("volcano_category", "point_color").rows() == [
        (plots.BACKGROUND, plots.CATEGORY_COLORS[plots.BACKGROUND]),
        (plots.INCREASED, plots.CATEGORY_COLORS[plots.INCREASED]),
        (plots.DECREASED, plots.CATEGORY_COLORS[plots.DECREASED]),
    ]


def test_prepare_volcano_data_keeps_all_diagnostic_warnings() -> None:
    artifacts = _artifacts()
    spec = plots.VOLCANO_PLOTS[0]

    plot_data = plots.prepare_volcano_data(artifacts.contrast_results, spec)

    assert plot_data.height == 252
    assert plot_data.filter(pl.col("diagnostic_warning")).height > 0
    assert plot_data.filter(pl.col("volcano_category") != plots.BACKGROUND).is_empty()


def test_plot_volcano_labels_bh_discoveries_with_twofold_effects() -> None:
    spec = plots.VOLCANO_PLOTS[0]
    analysis = pl.DataFrame(
        {
            "raw_row": [1, 2, 3, 4],
            "gene": ["label up", "label down", "small effect", "not discovered"],
            "contrast": [spec.key] * 4,
            "adjusted_log2_fold_change": [1.5, -1.5, 0.9, 1.5],
            "p_value": [0.001, 0.001, 0.001, 0.2],
            "q_value": [0.01, 0.01, 0.01, 0.2],
            "discovery": [True, True, True, False],
            "diagnostic_warning": [False, False, False, False],
        }
    )
    figure = plots.plot_volcano(
        analysis,
        spec,
    )

    assert {text.get_text() for text in figure.axes[0].texts} == {
        "label up",
        "label down",
    }
    plt.close(figure)


def test_plot_volcano_does_not_label_nominal_c1q_points() -> None:
    figure = plots.plot_volcano(
        _artifacts().contrast_results,
        plots.VOLCANO_PLOTS[0],
    )

    assert not figure.axes[0].texts
    plt.close(figure)


def test_effect_size_title_spells_out_benjamini_hochberg() -> None:
    figure = plots.plot_effect_sizes(_artifacts())

    assert {text.get_text() for text in figure.texts} == {
        "Leading effects ranked by Benjamini-Hochberg q-value"
    }
    plt.close(figure)


def test_write_volcano_plots_creates_png_and_svg_files(tmp_path: Path) -> None:
    artifacts = _artifacts()

    output_paths = plots.write_volcano_plots(
        analysis=artifacts.contrast_results,
        output_dir=tmp_path,
    )

    assert len(output_paths) == 6
    assert {path.suffix for path in output_paths.values()} == {".png", ".svg"}
    assert all(
        path.exists() and path.stat().st_size > 1000 for path in output_paths.values()
    )


def test_write_figures_creates_complete_figure_set(tmp_path: Path) -> None:
    output_paths = plots.write_figures(
        artifacts=_artifacts(),
        raw_path=RAW_DATA,
        output_dir=tmp_path,
        formats=("png",),
    )

    assert set(output_paths) == {
        "model-diagnostics.png",
        "qc-overview.png",
        "results-bortezomib.png",
        "results-c1q.png",
        "results-cholesterol.png",
        "results-chow-vs-lacz.png",
        "results-effect-sizes.png",
        "results-sensitivity.png",
    }
    assert all(
        path.exists() and path.stat().st_size > 10_000 for path in output_paths.values()
    )


def test_textxy_places_labels_away_from_origin() -> None:
    figure, ax = plt.subplots()

    annotations = plots.textxy(
        ax,
        x=[1.0, 1.0, -1.0, -1.0],
        y=[1.0, -1.0, 1.0, -1.0],
        labels=["upper right", "lower right", "upper left", "lower left"],
        offset_points=7.5,
        font_size=11.0,
    )

    assert [annotation.get_position() for annotation in annotations] == [
        (7.5, 7.5),
        (7.5, -7.5),
        (-7.5, 7.5),
        (-7.5, -7.5),
    ]
    assert [annotation.get_horizontalalignment() for annotation in annotations] == [
        "left",
        "left",
        "right",
        "right",
    ]
    assert [annotation.get_verticalalignment() for annotation in annotations] == [
        "bottom",
        "top",
        "bottom",
        "top",
    ]
    assert all(annotation.get_fontsize() == 11.0 for annotation in annotations)
    plt.close(figure)
