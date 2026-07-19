from pathlib import Path

import polars as pl
from matplotlib import pyplot as plt

import nrf1_proteomics.analysis
import nrf1_proteomics.plots

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DATA = PROJECT_ROOT / "data" / "nrf1-proteomics-raw.csv"


def _passing_assumptions(
    spec: nrf1_proteomics.plots.VolcanoPlotSpec, row_count: int
) -> dict[str, list[float]]:
    return {column: [0.5] * row_count for column in spec.assumption_pvalue_columns}


def test_prepare_volcano_data_assigns_original_point_colors() -> None:
    spec = nrf1_proteomics.plots.VOLCANO_PLOTS[0]
    analysis = pl.DataFrame(
        {
            "Gene": ["background", "pvalue", "fold_change", "both"],
            spec.x_column: [0.5, 0.5, 1.5, 1.5],
            spec.pvalue_column: [0.5, 0.01, 0.5, 0.01],
            **_passing_assumptions(spec, 4),
        }
    )

    plot_data = nrf1_proteomics.plots.prepare_volcano_data(analysis, spec)

    assert plot_data.select("volcano_category", "point_color").rows() == [
        (nrf1_proteomics.plots.BACKGROUND, "black"),
        (nrf1_proteomics.plots.PVALUE_ONLY, "red"),
        (nrf1_proteomics.plots.FOLD_CHANGE_ONLY, "orange"),
        (nrf1_proteomics.plots.PVALUE_AND_FOLD_CHANGE, "green"),
    ]


def test_prepare_volcano_data_excludes_invalid_significant_results() -> None:
    spec = nrf1_proteomics.plots.VOLCANO_PLOTS[0]
    assumptions = _passing_assumptions(spec, 5)
    for row_index, column in enumerate(spec.assumption_pvalue_columns):
        assumptions[column][row_index] = 0.01
    assumptions[spec.assumption_pvalue_columns[0]][3] = 0.01

    analysis = pl.DataFrame(
        {
            "Gene": [
                "levene_failure",
                "first_shapiro_failure",
                "second_shapiro_failure",
                "nonsignificant_failure",
                "valid_significant",
            ],
            spec.x_column: [0.5, 0.5, 1.5, 0.5, 1.5],
            spec.pvalue_column: [0.01, 0.01, 0.01, 0.5, 0.01],
            **assumptions,
        }
    )

    plot_data = nrf1_proteomics.plots.prepare_volcano_data(analysis, spec)

    assert plot_data.select("Gene", "volcano_category").rows() == [
        ("nonsignificant_failure", nrf1_proteomics.plots.BACKGROUND),
        ("valid_significant", nrf1_proteomics.plots.PVALUE_AND_FOLD_CHANGE),
    ]


def test_prepare_volcano_data_highlights_c1q_hits() -> None:
    analysis = nrf1_proteomics.analysis.analyze_raw_data(RAW_DATA)
    plot_data = nrf1_proteomics.plots.prepare_volcano_data(
        analysis, nrf1_proteomics.plots.VOLCANO_PLOTS[0]
    )

    c1q_hits = plot_data.filter(pl.col("Gene").is_in(["C1qa", "C1qb", "C1qc"]))

    assert c1q_hits.height == 3
    assert c1q_hits.select(
        pl.col("volcano_category").unique()
    ).to_series().to_list() == [nrf1_proteomics.plots.PVALUE_AND_FOLD_CHANGE]


def test_plot_volcano_labels_points_meeting_both_thresholds() -> None:
    spec = nrf1_proteomics.plots.VOLCANO_PLOTS[0]
    analysis = pl.DataFrame(
        {
            "Gene": ["background", "pvalue", "fold_change", "both_a", "both_b"],
            spec.x_column: [0.5, 0.5, 1.5, 1.5, -1.5],
            spec.pvalue_column: [0.5, 0.01, 0.5, 0.01, 0.01],
            **_passing_assumptions(spec, 5),
        }
    )

    fig = nrf1_proteomics.plots.plot_volcano(analysis, spec)

    assert {text.get_text() for text in fig.axes[0].texts} == {"both_a", "both_b"}
    plt.close(fig)


def test_write_volcano_plots_creates_png_files(tmp_path: Path) -> None:
    analysis = nrf1_proteomics.analysis.analyze_raw_data(RAW_DATA)

    output_paths = nrf1_proteomics.plots.write_volcano_plots(
        analysis=analysis, output_dir=tmp_path
    )

    assert set(output_paths) == {"ha_chol", "ha_bort", "ha_chow_lacz"}
    for path in output_paths.values():
        assert path.exists()
        assert path.suffix == ".png"
        assert path.stat().st_size > 1000


def test_textxy_places_labels_away_from_origin() -> None:
    fig, ax = plt.subplots()

    annotations = nrf1_proteomics.plots.textxy(
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
    assert all(annotation.get_fontweight() == "normal" for annotation in annotations)
    plt.close(fig)
