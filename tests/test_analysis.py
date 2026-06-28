import math
from collections.abc import Mapping
from pathlib import Path
from typing import cast

import polars as pl

from nrf1_proteomics.analysis import analyze_raw_data

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DATA = PROJECT_ROOT / "data" / "nrf1-proteomics-raw.csv"


def test_analysis_returns_historical_shape_and_columns() -> None:
    result = analyze_raw_data(RAW_DATA)

    assert result.shape == (252, 7)
    assert result.columns == [
        "Gene",
        "log2_FC_HAchol_HAchow",
        "log2_FC_HAchow_lacZ",
        "pvalue",
        "pvaluechow",
        "log2_FC_HAbort_HAchow",
        "pvaluebort",
    ]
    assert result.filter(pl.col("Gene") == "").height == 3


def test_analysis_matches_known_historical_hits() -> None:
    result = analyze_raw_data(RAW_DATA)

    assert_gene_values(
        result,
        "C1qc",
        log2_FC_HAchol_HAchow=1.257259532,
        log2_FC_HAchow_lacZ=0.089171005,
        pvalue=0.014190301,
        pvaluechow=0.860112069,
        log2_FC_HAbort_HAchow=0.094604099,
        pvaluebort=0.840348994,
    )
    assert_gene_values(
        result,
        "Sfn",
        log2_FC_HAchol_HAchow=0.413544935,
        log2_FC_HAchow_lacZ=-1.042830185,
        pvalue=0.403472758,
        pvaluechow=0.209067802,
        log2_FC_HAbort_HAchow=2.599948172,
        pvaluebort=0.013774412,
    )


def test_analysis_finds_expected_cholesterol_hits() -> None:
    result = analyze_raw_data(RAW_DATA)
    result = result.with_columns(
        delta_cholesterol=pl.col("log2_FC_HAchol_HAchow")
        - pl.col("log2_FC_HAchow_lacZ")
    )

    assert result.filter(pl.col("delta_cholesterol") > 0.58).height == 33


def assert_gene_values(
    result: pl.DataFrame,
    gene: str,
    **expected_values: float,
) -> None:
    rows = result.filter(pl.col("Gene") == gene)
    assert rows.height == 1
    row = cast(Mapping[str, object], rows.row(0, named=True))

    for column, expected_value in expected_values.items():
        actual_value = cast(float, row[column])
        assert math.isclose(actual_value, expected_value, abs_tol=1e-9)
