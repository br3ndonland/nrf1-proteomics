from pathlib import Path

import polars as pl
from matplotlib import pyplot as plt

from nrf1_proteomics.analysis import analyze_raw_data
from nrf1_proteomics.plots import (
    BACKGROUND,
    FOLD_CHANGE_ONLY,
    PVALUE_AND_FOLD_CHANGE,
    PVALUE_ONLY,
    VOLCANO_PLOTS,
    plot_volcano,
    prepare_volcano_data,
    write_volcano_plots,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DATA = PROJECT_ROOT / "data" / "nrf1-proteomics-raw.csv"


def test_prepare_volcano_data_assigns_original_point_colors() -> None:
    spec = VOLCANO_PLOTS[0]
    analysis = pl.DataFrame(
        {
            "Gene": ["background", "pvalue", "fold_change", "both"],
            spec.x_column: [0.5, 0.5, 1.5, 1.5],
            spec.pvalue_column: [0.5, 0.01, 0.5, 0.01],
        }
    )

    plot_data = prepare_volcano_data(analysis, spec)

    assert plot_data.select("volcano_category", "point_color").rows() == [
        (BACKGROUND, "black"),
        (PVALUE_ONLY, "red"),
        (FOLD_CHANGE_ONLY, "orange"),
        (PVALUE_AND_FOLD_CHANGE, "green"),
    ]


def test_prepare_volcano_data_highlights_c1q_hits() -> None:
    analysis = analyze_raw_data(RAW_DATA)
    plot_data = prepare_volcano_data(analysis, VOLCANO_PLOTS[0])

    c1q_hits = plot_data.filter(pl.col("Gene").is_in(["C1qa", "C1qb", "C1qc"]))

    assert c1q_hits.height == 3
    assert c1q_hits.select(
        pl.col("volcano_category").unique()
    ).to_series().to_list() == [PVALUE_AND_FOLD_CHANGE]


def test_plot_volcano_labels_points_meeting_both_thresholds() -> None:
    spec = VOLCANO_PLOTS[0]
    analysis = pl.DataFrame(
        {
            "Gene": ["background", "pvalue", "fold_change", "both_a", "both_b"],
            spec.x_column: [0.5, 0.5, 1.5, 1.5, -1.5],
            spec.pvalue_column: [0.5, 0.01, 0.5, 0.01, 0.01],
        }
    )

    fig = plot_volcano(analysis, spec)

    assert {text.get_text() for text in fig.axes[0].texts} == {"both_a", "both_b"}
    plt.close(fig)


def test_write_volcano_plots_creates_png_files(tmp_path: Path) -> None:
    analysis = analyze_raw_data(RAW_DATA)

    output_paths = write_volcano_plots(analysis=analysis, output_dir=tmp_path)

    assert set(output_paths) == {"ha_chol", "ha_bort", "ha_chow_lacz"}
    for path in output_paths.values():
        assert path.exists()
        assert path.suffix == ".png"
        assert path.stat().st_size > 1000
