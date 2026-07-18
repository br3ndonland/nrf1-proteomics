from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import cast

import polars as pl
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RAW_DATA = PROJECT_ROOT / "data" / "nrf1-proteomics-raw.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "nrf1-proteomics-analyzed.csv"

RAW_COLUMNS = (
    "protein_id",
    "gene",
    "description",
    "peptides_run_1",
    "peptides_run_2",
    "mix_run_1",
    "lacz_1",
    "lacz_2",
    "lacz_3",
    "ha_chow_1",
    "ha_chow_2",
    "ha_chol_1",
    "ha_chol_2",
    "ha_bort_1",
    "ha_bort_2",
    "mix_run_2",
    "lacz_4",
    "lacz_5",
    "lacz_6",
    "ha_chow_3",
    "ha_chow_4",
    "ha_chol_3",
    "ha_chol_4",
    "ha_bort_3",
    "ha_bort_4",
)
TEXT_COLUMNS = ("protein_id", "gene", "description")
PEPTIDE_COLUMNS = ("peptides_run_1", "peptides_run_2")
SIGNAL_COLUMNS = RAW_COLUMNS[5:]

SAMPLE_MIX_COLUMNS = {
    "lacz_1": "mix_run_1",
    "lacz_2": "mix_run_1",
    "lacz_3": "mix_run_1",
    "ha_chow_1": "mix_run_1",
    "ha_chow_2": "mix_run_1",
    "ha_chol_1": "mix_run_1",
    "ha_chol_2": "mix_run_1",
    "ha_bort_1": "mix_run_1",
    "ha_bort_2": "mix_run_1",
    "lacz_4": "mix_run_2",
    "lacz_5": "mix_run_2",
    "lacz_6": "mix_run_2",
    "ha_chow_3": "mix_run_2",
    "ha_chow_4": "mix_run_2",
    "ha_chol_3": "mix_run_2",
    "ha_chol_4": "mix_run_2",
    "ha_bort_3": "mix_run_2",
    "ha_bort_4": "mix_run_2",
}

GROUP_RATIO_COLUMNS = {
    "lacz": (
        "lacz_1_ratio",
        "lacz_2_ratio",
        "lacz_3_ratio",
        "lacz_4_ratio",
        "lacz_5_ratio",
        "lacz_6_ratio",
    ),
    "ha_chow": (
        "ha_chow_1_ratio",
        "ha_chow_2_ratio",
        "ha_chow_3_ratio",
        "ha_chow_4_ratio",
    ),
    "ha_chol": (
        "ha_chol_1_ratio",
        "ha_chol_2_ratio",
        "ha_chol_3_ratio",
        "ha_chol_4_ratio",
    ),
    "ha_bort": (
        "ha_bort_1_ratio",
        "ha_bort_2_ratio",
        "ha_bort_3_ratio",
        "ha_bort_4_ratio",
    ),
}

OUTPUT_COLUMNS = (
    "Gene",
    "log2_FC_HAchol_HAchow",
    "log2_FC_HAchow_lacZ",
    "pvalue",
    "pvaluechow",
    "log2_FC_HAbort_HAchow",
    "pvaluebort",
)


def read_raw_data(path: str | Path = DEFAULT_RAW_DATA) -> pl.DataFrame:
    """Read the TCMP mass spectrometry export into typed Polars columns."""

    return (
        pl.read_csv(
            path,
            has_header=False,
            skip_rows=4,
            new_columns=list(RAW_COLUMNS),
            infer_schema_length=0,
        )
        .select(
            pl.col(TEXT_COLUMNS).fill_null(""),
            pl.col(PEPTIDE_COLUMNS).cast(pl.Int64),
            pl.col(SIGNAL_COLUMNS).cast(pl.Float64),
        )
        .with_row_index("_row")
    )


def filter_confident_proteins(raw_data: pl.DataFrame) -> pl.DataFrame:
    """Keep proteins with at least two quantified peptides in both runs."""

    return raw_data.filter(
        pl.col("peptides_run_1") >= 2,
        pl.col("peptides_run_2") >= 2,
    )


def add_normalized_ratios(protein_data: pl.DataFrame) -> pl.DataFrame:
    """Normalize each reporter-channel signal to the mix channel for its run."""

    return protein_data.with_columns(
        (pl.col(sample_column) / pl.col(mix_column)).alias(f"{sample_column}_ratio")
        for sample_column, mix_column in SAMPLE_MIX_COLUMNS.items()
    )


def analyze_raw_data(path: str | Path = DEFAULT_RAW_DATA) -> pl.DataFrame:
    """Return the analyzed Nrf1 proteomics result table as a Polars DataFrame."""

    ratios = add_normalized_ratios(filter_confident_proteins(read_raw_data(path)))
    fold_changes = _calculate_fold_changes(ratios)
    pvalues = _calculate_pvalues(ratios)

    return (
        fold_changes.join(pvalues, on="protein_id", how="inner")
        .sort("_row")
        .select(OUTPUT_COLUMNS)
    )


def write_analysis(
    output_path: str | Path = DEFAULT_OUTPUT,
    raw_path: str | Path = DEFAULT_RAW_DATA,
) -> pl.DataFrame:
    """Write the analyzed result table to CSV and return it."""

    result = analyze_raw_data(raw_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    result.write_csv(output_path)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    _ = parser.add_argument("--raw", type=Path, default=DEFAULT_RAW_DATA)
    _ = parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    raw_path = cast(Path, args.raw)
    output_path = cast(Path, args.output)

    _ = write_analysis(output_path=output_path, raw_path=raw_path)
    return 0


def _calculate_fold_changes(ratios: pl.DataFrame) -> pl.DataFrame:
    return ratios.with_columns(
        lacz_mean=_mean_expr(GROUP_RATIO_COLUMNS["lacz"]),
        ha_chow_mean=_mean_expr(GROUP_RATIO_COLUMNS["ha_chow"]),
        ha_chol_mean=_mean_expr(GROUP_RATIO_COLUMNS["ha_chol"]),
        ha_bort_mean=_mean_expr(GROUP_RATIO_COLUMNS["ha_bort"]),
    ).select(
        "_row",
        "protein_id",
        pl.col("gene").alias("Gene"),
        (pl.col("ha_chol_mean") / pl.col("ha_chow_mean"))
        .log(base=2)
        .alias("log2_FC_HAchol_HAchow"),
        (pl.col("ha_chow_mean") / pl.col("lacz_mean"))
        .log(base=2)
        .alias("log2_FC_HAchow_lacZ"),
        (pl.col("ha_bort_mean") / pl.col("ha_chow_mean"))
        .log(base=2)
        .alias("log2_FC_HAbort_HAchow"),
    )


def _calculate_pvalues(ratios: pl.DataFrame) -> pl.DataFrame:
    rows: list[dict[str, str | float]] = []
    ratio_columns = (
        "protein_id",
        *GROUP_RATIO_COLUMNS["lacz"],
        *GROUP_RATIO_COLUMNS["ha_chow"],
        *GROUP_RATIO_COLUMNS["ha_chol"],
        *GROUP_RATIO_COLUMNS["ha_bort"],
    )

    for raw_row in ratios.select(ratio_columns).iter_rows(named=True):
        row = cast(Mapping[str, object], raw_row)
        rows.append(
            {
                "protein_id": str(row["protein_id"]),
                "pvalue": _t_test_pvalue(
                    _row_values(row, GROUP_RATIO_COLUMNS["ha_chol"]),
                    _row_values(row, GROUP_RATIO_COLUMNS["ha_chow"]),
                ),
                "pvaluechow": _t_test_pvalue(
                    _row_values(row, GROUP_RATIO_COLUMNS["ha_chow"]),
                    _row_values(row, GROUP_RATIO_COLUMNS["lacz"]),
                ),
                "pvaluebort": _t_test_pvalue(
                    _row_values(row, GROUP_RATIO_COLUMNS["ha_bort"]),
                    _row_values(row, GROUP_RATIO_COLUMNS["ha_chow"]),
                ),
            }
        )

    return pl.DataFrame(rows)


def _mean_expr(columns: Sequence[str]) -> pl.Expr:
    if not columns:
        msg = "At least one column is required."
        raise ValueError(msg)

    total = pl.col(columns[0])
    for column in columns[1:]:
        total = total + pl.col(column)
    return total / len(columns)


def _row_values(row: Mapping[str, object], columns: Sequence[str]) -> list[float]:
    return [cast(float, row[column]) for column in columns]


def _t_test_pvalue(sample: Sequence[float], reference: Sequence[float]) -> float:
    return cast(
        float,
        stats.ttest_ind(
            sample,
            reference,
            equal_var=True,
            alternative="two-sided",
        )[1],
    )
