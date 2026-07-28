from __future__ import annotations

import argparse
import csv
import json
import math
import re
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import numpy as np
import polars as pl
import scipy.stats

from nrf1_proteomics.moderated import (
    FloatArray,
    ModeratedModelFit,
    fit_moderated_model,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
_PROJECT_RAW_DATA = PROJECT_ROOT / "data" / "nrf1-proteomics-raw.csv"
_PACKAGED_RAW_DATA = (
    Path(__file__).resolve().parent / "data" / "nrf1-proteomics-raw.csv"
)
DEFAULT_RAW_DATA = (
    _PROJECT_RAW_DATA if _PROJECT_RAW_DATA.is_file() else _PACKAGED_RAW_DATA
)
DEFAULT_OUTPUT_DIR = Path(tempfile.gettempdir()) / "nrf1-proteomics" / "derived"

TEXT_COLUMNS = (
    "protein_id",
    "gene",
    "description",
)
PEPTIDE_COLUMNS = (
    "peptides_run_1",
    "peptides_run_2",
)
RAW_COLUMN_COUNT = 25

GROUP_LABELS: dict[str, str] = {
    "lacz": "lacZ",
    "ha_chow": "chow",
    "ha_chol": "cholesterol",
    "ha_bort": "bortezomib",
}
GROUP_ORDER = tuple(GROUP_LABELS)
TREATMENT_BY_EXPORT_LABEL = {
    "lacz": "lacz",
    "hacont": "ha_chow",
    "hachol": "ha_chol",
    "habor": "ha_bort",
}
FDR_THRESHOLD = 0.05
TREAT_LFC = math.log2(1.5)


@dataclass(frozen=True)
class ExportHeader:
    sample_labels: tuple[str, ...]
    reporter_tags: tuple[str, ...]
    run_ids: tuple[str, ...]


@dataclass(frozen=True)
class ContrastSpec:
    numerator: str
    denominator: str

    @property
    def key(self) -> str:
        return f"{self.numerator}_vs_{self.denominator}"

    @property
    def label(self) -> str:
        return f"{GROUP_LABELS[self.numerator]} - {GROUP_LABELS[self.denominator]}"

    @property
    def treatment_weights(self) -> tuple[float, ...]:
        return tuple(
            float(group == self.numerator) - float(group == self.denominator)
            for group in GROUP_ORDER
        )


CONTRASTS = (
    ContrastSpec("ha_chol", "ha_chow"),
    ContrastSpec("ha_bort", "ha_chow"),
    ContrastSpec("ha_chow", "lacz"),
)


@dataclass(frozen=True)
class AnalysisArtifacts:
    sample_metadata: pl.DataFrame
    normalized_data: pl.DataFrame
    contrast_results: pl.DataFrame
    omnibus_results: pl.DataFrame
    diagnostics: pl.DataFrame
    filter_audit: pl.DataFrame
    sensitivities: pl.DataFrame


@dataclass(frozen=True)
class ModelResult:
    contrast_results: pl.DataFrame
    omnibus_results: pl.DataFrame
    diagnostics: pl.DataFrame


def read_raw_data(path: str | Path = DEFAULT_RAW_DATA) -> pl.DataFrame:
    """Read and validate the 25-column TCMP mass spectrometry export."""

    raw_path = Path(path)
    header = _read_export_header(raw_path)
    signal_columns = _canonical_signal_columns(header)
    raw_columns = (*TEXT_COLUMNS, *PEPTIDE_COLUMNS, *signal_columns)
    data = (
        pl.read_csv(
            raw_path,
            has_header=False,
            skip_rows=4,
            new_columns=list(raw_columns),
            infer_schema_length=0,
        )
        .select(
            pl.col(TEXT_COLUMNS).fill_null(""),
            pl.col(PEPTIDE_COLUMNS).cast(pl.Int64, strict=True),
            pl.col(signal_columns).cast(pl.Float64, strict=True),
        )
        .with_row_index("raw_row", offset=1)
    )
    if data.width != RAW_COLUMN_COUNT + 1:
        message = f"Expected {RAW_COLUMN_COUNT} raw columns, found {data.width - 1}."
        raise ValueError(message)
    if data.get_column("protein_id").n_unique() != data.height:
        raise ValueError("protein_id must be unique for protein-wise inference.")
    if data.select(pl.any_horizontal(pl.col(signal_columns).is_null()).any()).item():
        raise ValueError("Signal columns must contain finite numeric values.")
    return data


def build_sample_metadata(
    path: str | Path = DEFAULT_RAW_DATA,
) -> pl.DataFrame:
    """Read sample labels, reporter tags, and TMT runs from the export header."""

    header = _read_export_header(Path(path))
    signal_columns = _canonical_signal_columns(header)
    channels = tuple(
        zip(
            signal_columns,
            header.sample_labels,
            header.reporter_tags,
            header.run_ids,
            strict=True,
        )
    )
    run_order = tuple(dict.fromkeys(header.run_ids))
    run_numbers = {run_id: index + 1 for index, run_id in enumerate(run_order)}
    mix_columns = {
        run_id: column
        for column, sample_label, _, run_id in channels
        if sample_label.casefold() == "mix"
    }
    if set(mix_columns) != set(run_order):
        raise ValueError("Expected one pooled mix channel in each TMT run.")

    rows: list[dict[str, object]] = []
    for sample, sample_label, reporter_tag, run_id in channels:
        if sample.startswith("mix_"):
            continue
        treatment, replicate_text = sample.rsplit("_", maxsplit=1)
        if treatment not in GROUP_LABELS or not replicate_text.isdigit():
            raise ValueError(f"Invalid canonical biological channel: {sample}")
        rows.append(
            {
                "sample": sample,
                "sample_label": sample_label,
                "treatment": treatment,
                "treatment_label": GROUP_LABELS[treatment],
                "run": f"run_{run_numbers[run_id]}",
                "tmt_run_id": run_id,
                "replicate": int(replicate_text),
                "reporter_tag": reporter_tag,
                "mix_column": mix_columns[run_id],
                "biological_replicate": True,
            }
        )
    return pl.DataFrame(rows)


def build_filter_audit(raw_data: pl.DataFrame) -> pl.DataFrame:
    """Record primary and one-peptide sensitivity eligibility for every protein."""

    audit_rows: list[dict[str, object]] = []
    biological_columns = _biological_columns(raw_data)
    columns = (*TEXT_COLUMNS, *PEPTIDE_COLUMNS, *biological_columns, "raw_row")
    for raw_row in raw_data.select(columns).iter_rows(named=True):
        row = cast(Mapping[str, object], raw_row)
        positive = all(
            float(cast(float, row[column])) > 0 for column in biological_columns
        )
        peptides_1 = int(cast(int, row["peptides_run_1"]))
        peptides_2 = int(cast(int, row["peptides_run_2"]))
        primary_reasons = _exclusion_reasons(
            peptides_1, peptides_2, positive, minimum=2
        )
        sensitivity_reasons = _exclusion_reasons(
            peptides_1, peptides_2, positive, minimum=1
        )
        audit_rows.append(
            {
                "raw_row": int(cast(int, row["raw_row"])),
                "protein_id": str(row["protein_id"]),
                "gene": str(row["gene"]),
                "description": str(row["description"]),
                "peptides_run_1": peptides_1,
                "peptides_run_2": peptides_2,
                "all_biological_signals_positive": positive,
                "included_primary": not primary_reasons,
                "primary_exclusion_reason": ";".join(primary_reasons) or "included",
                "included_one_peptide_sensitivity": not sensitivity_reasons,
                "sensitivity_exclusion_reason": ";".join(sensitivity_reasons)
                or "included",
            }
        )
    return pl.DataFrame(audit_rows).sort("raw_row")


def build_design_matrix(
    sample_metadata: pl.DataFrame | None = None,
) -> pl.DataFrame:
    """Build the full-rank four-group cell-means design with a fixed run block."""

    metadata = build_sample_metadata() if sample_metadata is None else sample_metadata
    run_levels = cast(
        list[str],
        metadata.get_column("run").unique(maintain_order=True).to_list(),
    )
    return metadata.select(
        *(
            (pl.col("treatment") == group).cast(pl.Float64).alias(group)
            for group in GROUP_ORDER
        ),
        *((pl.col("run") == run).cast(pl.Float64).alias(run) for run in run_levels[1:]),
    )


def normalize_expression(
    raw_data: pl.DataFrame,
    filter_audit: pl.DataFrame,
    sample_metadata: pl.DataFrame | None = None,
    *,
    minimum_peptides: int = 2,
    median_center: bool = True,
    mix_reference: bool = False,
) -> tuple[pl.DataFrame, FloatArray]:
    """Log2 transform eligible mouse channels and optionally normalize them."""

    metadata = build_sample_metadata() if sample_metadata is None else sample_metadata
    sample_columns = cast(
        list[str],
        metadata.get_column("sample").to_list(),
    )
    mix_columns = cast(
        list[str],
        metadata.get_column("mix_column").to_list(),
    )
    eligibility_column = (
        "included_primary"
        if minimum_peptides == 2
        else "included_one_peptide_sensitivity"
    )
    eligible_ids = (
        filter_audit.filter(pl.col(eligibility_column))
        .get_column("protein_id")
        .to_list()
    )
    proteins = raw_data.filter(pl.col("protein_id").is_in(eligible_ids)).sort("raw_row")
    expression = cast(
        FloatArray,
        np.asarray(
            proteins.select(sample_columns).to_numpy(),
            dtype=np.float64,
        ),
    )
    if not np.all(np.isfinite(expression)) or np.any(expression <= 0):
        raise ValueError("Modeled reporter signals must be positive and finite.")
    expression = np.log2(expression)
    if mix_reference:
        mix_values = np.asarray(
            np.column_stack(
                [
                    proteins.get_column(mix_column).to_numpy()
                    for mix_column in mix_columns
                ]
            ),
            dtype=np.float64,
        )
        if not np.all(np.isfinite(mix_values)) or np.any(mix_values <= 0):
            raise ValueError(
                "Pooled mix signals must be positive for ratio sensitivity."
            )
        expression = expression - np.log2(mix_values)
    if median_center:
        channel_medians = np.median(expression, axis=0)
        grand_median = float(cast(np.float64, np.median(expression)))
        expression = expression - channel_medians + grand_median

    normalized = proteins.select("raw_row", *TEXT_COLUMNS).with_columns(
        pl.Series(sample, expression[:, index])
        for index, sample in enumerate(sample_columns)
    )
    return normalized, expression


def fit_protein_models(
    normalized_data: pl.DataFrame,
    expression: FloatArray,
    sample_metadata: pl.DataFrame,
    *,
    trend: bool = True,
    robust: bool = True,
) -> ModelResult:
    """Fit the blocked protein-wise model and its planned moderated contrasts."""

    design = cast(
        FloatArray,
        np.asarray(
            build_design_matrix(sample_metadata).to_numpy(),
            dtype=np.float64,
        ),
    )
    design_columns = cast(int, design.shape[1])
    if np.linalg.matrix_rank(design) != design_columns:
        raise ValueError("The treatment and run design matrix must be full rank.")
    if expression.shape != (normalized_data.height, sample_metadata.height):
        raise ValueError(
            "Expression dimensions do not match protein and sample metadata."
        )

    identifiers = normalized_data.select("raw_row", *TEXT_COLUMNS)
    contrast_matrix = np.asarray(
        [
            (
                *contrast.treatment_weights,
                *([0.0] * (design_columns - len(GROUP_ORDER))),
            )
            for contrast in CONTRASTS
        ],
        dtype=np.float64,
    ).T
    moderated_fit = fit_moderated_model(
        expression,
        design,
        contrast_matrix,
        trend=trend,
        robust=robust,
        treat_lfc=TREAT_LFC,
    )

    contrast_results = _contrast_table(identifiers, moderated_fit)
    omnibus_results = _omnibus_table(identifiers, moderated_fit)
    diagnostics = _diagnostic_table(
        identifiers,
        expression,
        design,
        sample_metadata,
        moderated_fit,
    )
    return ModelResult(contrast_results, omnibus_results, diagnostics)


def analyze_raw_data(path: str | Path = DEFAULT_RAW_DATA) -> AnalysisArtifacts:
    """Run the primary analysis and prespecified sensitivity analyses."""

    raw_path = Path(path)
    raw_data = read_raw_data(raw_path)
    sample_metadata = build_sample_metadata(raw_path)
    filter_audit = build_filter_audit(raw_data)
    normalized_data, expression = normalize_expression(
        raw_data,
        filter_audit,
        sample_metadata,
    )
    primary = fit_protein_models(normalized_data, expression, sample_metadata)
    contrast_results = primary.contrast_results.join(
        primary.diagnostics.select("protein_id", "diagnostic_warning"),
        on="protein_id",
        how="left",
        validate="m:1",
    )

    sensitivity_frames = [
        _sensitivity_frame("primary", contrast_results),
        _run_sensitivity(
            "no_median_centering",
            raw_data,
            filter_audit,
            sample_metadata,
            median_center=False,
        ),
        _run_sensitivity(
            "mix_reference_ratios",
            raw_data,
            filter_audit,
            sample_metadata,
            mix_reference=True,
        ),
        _run_sensitivity(
            "standard_ebayes",
            raw_data,
            filter_audit,
            sample_metadata,
            trend=False,
            robust=False,
        ),
        _run_sensitivity(
            "one_peptide_positive",
            raw_data,
            filter_audit,
            sample_metadata,
            minimum_peptides=1,
        ),
    ]
    sensitivities = pl.concat(sensitivity_frames, how="vertical").sort(
        "scenario", "contrast", "raw_row"
    )
    return AnalysisArtifacts(
        sample_metadata=sample_metadata,
        normalized_data=normalized_data,
        contrast_results=contrast_results,
        omnibus_results=primary.omnibus_results,
        diagnostics=primary.diagnostics,
        filter_audit=filter_audit,
        sensitivities=sensitivities,
    )


def write_analysis(
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
    raw_path: str | Path = DEFAULT_RAW_DATA,
) -> AnalysisArtifacts:
    """Write all public analysis artifacts and return them."""

    artifacts = analyze_raw_data(raw_path)
    derived_dir = Path(output_dir)
    derived_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "sample_metadata": artifacts.sample_metadata,
        "normalized": artifacts.normalized_data,
        "contrasts": artifacts.contrast_results,
        "omnibus": artifacts.omnibus_results,
        "diagnostics": artifacts.diagnostics,
        "filter_audit": artifacts.filter_audit,
        "sensitivities": artifacts.sensitivities,
    }
    for name, frame in outputs.items():
        frame.write_csv(derived_dir / f"nrf1-proteomics-{name}.csv")
    summary_path = derived_dir / "nrf1-proteomics-summary.json"
    _ = summary_path.write_text(
        json.dumps(analysis_summary(artifacts), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return artifacts


def analysis_summary(artifacts: AnalysisArtifacts) -> dict[str, object]:
    """Return a JSON-compatible summary of analysis settings and results."""

    discovery_counts = {
        contrast.key: artifacts.contrast_results.filter(
            (pl.col("contrast") == contrast.key) & pl.col("discovery")
        ).height
        for contrast in CONTRASTS
    }
    treat_counts = {
        contrast.key: artifacts.contrast_results.filter(
            (pl.col("contrast") == contrast.key) & pl.col("treat_discovery")
        ).height
        for contrast in CONTRASTS
    }
    sensitivity_counts: dict[str, dict[str, int]] = {}
    scenarios = cast(
        list[str],
        artifacts.sensitivities.get_column("scenario").unique().sort().to_list(),
    )
    for scenario in scenarios:
        scenario_data = artifacts.sensitivities.filter(pl.col("scenario") == scenario)
        sensitivity_counts[scenario] = {
            contrast.key: scenario_data.filter(
                (pl.col("contrast") == contrast.key) & pl.col("discovery")
            ).height
            for contrast in CONTRASTS
        }
    design = build_design_matrix(artifacts.sample_metadata)
    design_values = cast(
        FloatArray,
        np.asarray(design.to_numpy(), dtype=np.float64),
    )
    return {
        "raw_proteins": artifacts.filter_audit.height,
        "primary_proteins": artifacts.normalized_data.height,
        "one_peptide_sensitivity_proteins": artifacts.filter_audit.filter(
            pl.col("included_one_peptide_sensitivity")
        ).height,
        "biological_samples": artifacts.sample_metadata.height,
        "design_rank": int(cast(int, np.linalg.matrix_rank(design_values))),
        "residual_degrees_of_freedom": artifacts.sample_metadata.height - design.width,
        "model": "log2_abundance ~ 0 + treatment + TMT_run",
        "empirical_bayes": {"trend": True, "robust": True},
        "fdr_method": "Benjamini-Hochberg within each contrast",
        "fdr_threshold": FDR_THRESHOLD,
        "treat_log2_threshold": TREAT_LFC,
        "primary_discoveries": discovery_counts,
        "treat_discoveries": treat_counts,
        "omnibus_discoveries": artifacts.omnibus_results.filter(
            pl.col("discovery")
        ).height,
        "diagnostic_warnings": artifacts.diagnostics.filter(
            pl.col("diagnostic_warning")
        ).height,
        "sensitivity_discoveries": sensitivity_counts,
    }


def calculate_model_residuals(
    normalized_data: pl.DataFrame,
    sample_metadata: pl.DataFrame,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return residuals, standardized residuals, and leverage for plotting."""

    sample_columns = cast(
        list[str],
        sample_metadata.get_column("sample").to_list(),
    )
    expression = cast(
        FloatArray,
        np.asarray(
            normalized_data.select(sample_columns).to_numpy(),
            dtype=np.float64,
        ),
    )
    design = cast(
        FloatArray,
        np.asarray(
            build_design_matrix(sample_metadata).to_numpy(),
            dtype=np.float64,
        ),
    )
    coefficients = np.linalg.lstsq(design, expression.T, rcond=None)[0]
    residuals = expression.T - design @ coefficients
    hat = design @ np.linalg.inv(design.T @ design) @ design.T
    leverage = np.diag(hat)
    degrees_of_freedom = cast(int, design.shape[0]) - cast(int, design.shape[1])
    mean_squared_error = np.sum(residuals**2, axis=0) / degrees_of_freedom
    scale = cast(
        FloatArray,
        np.sqrt(np.outer(1.0 - leverage, mean_squared_error)),
    )
    standardized = cast(
        FloatArray,
        np.divide(
            residuals,
            scale,
            out=np.zeros_like(residuals),
            where=scale > 0,
        ),
    )
    return residuals, standardized, leverage


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    _ = parser.add_argument("--raw", type=Path, default=DEFAULT_RAW_DATA)
    _ = parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(argv)
    _ = write_analysis(
        output_dir=cast(Path, args.output_dir),
        raw_path=cast(Path, args.raw),
    )
    return 0


def _read_export_header(path: Path) -> ExportHeader:
    with path.open(encoding="utf-8-sig", newline="") as file:
        reader = csv.reader(file)
        try:
            rows = [next(reader) for _ in range(4)]
        except StopIteration as error:
            raise ValueError("Expected four TCMP export header rows.") from error
    if any(len(row) != RAW_COLUMN_COUNT for row in rows):
        raise ValueError(f"Expected a {RAW_COLUMN_COUNT}-column TCMP export.")

    sample_labels = tuple(value.strip() for value in rows[2][5:])
    reporter_tags = tuple(value.strip() for value in rows[3][5:])
    if not all(sample_labels):
        raise ValueError("Reporter-channel sample labels must not be empty.")
    if not all(reporter_tags):
        raise ValueError("TMT reporter tags must not be empty.")

    run_ids: list[str] = []
    current_run = ""
    for value in rows[1][5:]:
        current_run = value.strip() or current_run
        if not current_run:
            raise ValueError("Each reporter channel must belong to a TMT run.")
        run_ids.append(current_run)
    return ExportHeader(sample_labels, reporter_tags, tuple(run_ids))


def _canonical_signal_columns(header: ExportHeader) -> tuple[str, ...]:
    run_order = tuple(dict.fromkeys(header.run_ids))
    run_numbers = {run_id: index + 1 for index, run_id in enumerate(run_order)}
    columns: list[str] = []
    for sample_label, run_id in zip(
        header.sample_labels,
        header.run_ids,
        strict=True,
    ):
        if sample_label.casefold() == "mix":
            columns.append(f"mix_run_{run_numbers[run_id]}")
            continue
        match = re.fullmatch(r"(.+?)[\s_-]*(\d+)", sample_label)
        if match is None:
            raise ValueError(f"Cannot parse biological sample label: {sample_label}")
        export_treatment = re.sub(r"[^a-z]", "", match.group(1).casefold())
        treatment = TREATMENT_BY_EXPORT_LABEL.get(export_treatment)
        if treatment is None:
            raise ValueError(f"Unknown treatment label: {sample_label}")
        columns.append(f"{treatment}_{int(match.group(2))}")
    if len(columns) != len(set(columns)):
        raise ValueError("Reporter-channel sample identifiers must be unique.")
    return tuple(columns)


def _biological_columns(raw_data: pl.DataFrame) -> tuple[str, ...]:
    fixed_columns = {"raw_row", *TEXT_COLUMNS, *PEPTIDE_COLUMNS}
    return tuple(
        column
        for column in raw_data.columns
        if column not in fixed_columns and not column.startswith("mix_")
    )


def _exclusion_reasons(
    peptides_1: int,
    peptides_2: int,
    positive: bool,
    *,
    minimum: int,
) -> list[str]:
    reasons: list[str] = []
    if peptides_1 < minimum:
        reasons.append(f"peptides_run_1_lt_{minimum}")
    if peptides_2 < minimum:
        reasons.append(f"peptides_run_2_lt_{minimum}")
    if not positive:
        reasons.append("nonpositive_biological_signal")
    return reasons


def _array_scalar_1d(values: FloatArray, index: int) -> float:
    return float(cast(np.float64, values[index]))


def _array_scalar_2d(values: FloatArray, row: int, column: int) -> float:
    return float(cast(np.float64, values[row, column]))


def _contrast_table(
    identifiers: pl.DataFrame,
    moderated_fit: ModeratedModelFit,
) -> pl.DataFrame:
    coefficients = moderated_fit.coefficients
    stdev_unscaled = moderated_fit.stdev_unscaled
    posterior_variance = moderated_fit.posterior_variance
    standard_errors = stdev_unscaled * np.sqrt(posterior_variance)[:, np.newaxis]
    t_statistic = moderated_fit.moderated_t
    p_values = moderated_fit.p_value
    total_df = moderated_fit.total_degrees_of_freedom
    treat_t = moderated_fit.treat_t
    treat_p = moderated_fit.treat_p_value
    annotations = identifiers.to_dicts()
    rows: list[dict[str, object]] = []
    for column, contrast in enumerate(CONTRASTS):
        q_values = adjust_bh(p_values[:, column])
        treat_q_values = adjust_bh(treat_p[:, column])
        critical_t = np.asarray(
            scipy.stats.t.ppf(0.975, total_df),
            dtype=np.float64,
        )
        for row_index, annotation in enumerate(annotations):
            log2_fold_change = _array_scalar_2d(coefficients, row_index, column)
            standard_error = _array_scalar_2d(standard_errors, row_index, column)
            discovery = _array_scalar_1d(q_values, row_index) < FDR_THRESHOLD
            rows.append(
                {
                    **annotation,
                    "contrast": contrast.key,
                    "contrast_label": contrast.label,
                    "numerator": contrast.numerator,
                    "denominator": contrast.denominator,
                    "adjusted_log2_fold_change": log2_fold_change,
                    "linear_fold_change": float(2**log2_fold_change),
                    "moderated_standard_error": standard_error,
                    "ci_lower_95": log2_fold_change
                    - _array_scalar_1d(critical_t, row_index) * standard_error,
                    "ci_upper_95": log2_fold_change
                    + _array_scalar_1d(critical_t, row_index) * standard_error,
                    "moderated_t": _array_scalar_2d(t_statistic, row_index, column),
                    "total_degrees_of_freedom": _array_scalar_1d(total_df, row_index),
                    "p_value": _array_scalar_2d(p_values, row_index, column),
                    "q_value": _array_scalar_1d(q_values, row_index),
                    "discovery": discovery,
                    "direction": (
                        "increased"
                        if discovery and log2_fold_change > 0
                        else "decreased"
                        if discovery
                        else "not_discovered"
                    ),
                    "treat_log2_threshold": TREAT_LFC,
                    "treat_t": _array_scalar_2d(treat_t, row_index, column),
                    "treat_p_value": _array_scalar_2d(treat_p, row_index, column),
                    "treat_q_value": _array_scalar_1d(treat_q_values, row_index),
                    "treat_discovery": (
                        _array_scalar_1d(treat_q_values, row_index) < FDR_THRESHOLD
                    ),
                }
            )
    return pl.DataFrame(rows).sort("contrast", "raw_row")


def _omnibus_table(
    identifiers: pl.DataFrame,
    moderated_fit: ModeratedModelFit,
) -> pl.DataFrame:
    f_statistic = moderated_fit.moderated_f
    p_values = moderated_fit.f_p_value
    q_values = adjust_bh(p_values)
    total_df = moderated_fit.total_degrees_of_freedom
    return identifiers.with_columns(
        pl.Series("moderated_f", f_statistic),
        pl.lit(len(CONTRASTS)).alias("numerator_degrees_of_freedom"),
        pl.Series("total_degrees_of_freedom", total_df),
        pl.Series("p_value", p_values),
        pl.Series("q_value", q_values),
        pl.Series("discovery", q_values < FDR_THRESHOLD),
    ).sort("raw_row")


def _diagnostic_table(
    identifiers: pl.DataFrame,
    expression: FloatArray,
    design: FloatArray,
    sample_metadata: pl.DataFrame,
    moderated_fit: ModeratedModelFit,
) -> pl.DataFrame:
    coefficients = np.linalg.lstsq(design, expression.T, rcond=None)[0]
    residuals = expression.T - design @ coefficients
    hat = design @ np.linalg.inv(design.T @ design) @ design.T
    leverage = np.diag(hat)
    design_columns = cast(int, design.shape[1])
    residual_df = cast(int, design.shape[0]) - design_columns
    mean_squared_error = np.sum(residuals**2, axis=0) / residual_df
    cooks_threshold = 4.0 / residual_df
    groups = np.asarray(
        sample_metadata.get_column("treatment").to_numpy(),
        dtype=np.str_,
    )
    rows: list[dict[str, object]] = []
    posterior_variance = moderated_fit.posterior_variance
    prior_variance = moderated_fit.prior_variance
    prior_df = moderated_fit.prior_degrees_of_freedom
    annotations = identifiers.to_dicts()
    for protein_index, annotation in enumerate(annotations):
        protein_residuals = residuals[:, protein_index]
        shapiro = scipy.stats.shapiro(protein_residuals)
        brown_forsythe = scipy.stats.levene(
            *(protein_residuals[groups == group] for group in GROUP_ORDER),
            center="median",
        )
        shapiro_statistic = cast(float, shapiro.statistic)
        shapiro_p_value = cast(float, shapiro.pvalue)
        brown_forsythe_statistic = cast(float, brown_forsythe.statistic)
        brown_forsythe_p_value = cast(float, brown_forsythe.pvalue)
        mse = _array_scalar_1d(mean_squared_error, protein_index)
        if mse > 0:
            cooks = (
                protein_residuals**2
                / (design_columns * mse)
                * leverage
                / (1.0 - leverage) ** 2
            )
            max_cooks = float(np.max(cooks))
        else:
            max_cooks = 0.0
        shapiro_warning = bool(shapiro_p_value < 0.05)
        brown_forsythe_warning = bool(brown_forsythe_p_value < 0.05)
        cooks_warning = bool(max_cooks > cooks_threshold)
        rows.append(
            {
                **annotation,
                "mean_log2_abundance": float(np.mean(expression[protein_index, :])),
                "raw_residual_variance": mse,
                "posterior_residual_variance": _array_scalar_1d(
                    posterior_variance, protein_index
                ),
                "prior_variance": _array_scalar_1d(prior_variance, protein_index),
                "prior_degrees_of_freedom": _array_scalar_1d(prior_df, protein_index),
                "residual_shapiro_statistic": shapiro_statistic,
                "residual_shapiro_p_value": shapiro_p_value,
                "brown_forsythe_statistic": brown_forsythe_statistic,
                "brown_forsythe_p_value": brown_forsythe_p_value,
                "max_cooks_distance": max_cooks,
                "cooks_distance_threshold": cooks_threshold,
                "shapiro_warning": shapiro_warning,
                "brown_forsythe_warning": brown_forsythe_warning,
                "cooks_distance_warning": cooks_warning,
                "diagnostic_warning": (
                    shapiro_warning or brown_forsythe_warning or cooks_warning
                ),
            }
        )
    return pl.DataFrame(rows).sort("raw_row")


def adjust_bh(p_values: Sequence[float] | FloatArray) -> FloatArray:
    """Adjust one prespecified family of p-values with Benjamini-Hochberg."""

    values = np.asarray(p_values, dtype=np.float64)
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        raise ValueError("BH adjustment requires one finite p-value vector.")
    return np.asarray(
        scipy.stats.false_discovery_control(values, method="bh"),
        dtype=np.float64,
    )


def _run_sensitivity(
    name: str,
    raw_data: pl.DataFrame,
    filter_audit: pl.DataFrame,
    sample_metadata: pl.DataFrame,
    *,
    minimum_peptides: int = 2,
    median_center: bool = True,
    mix_reference: bool = False,
    trend: bool = True,
    robust: bool = True,
) -> pl.DataFrame:
    normalized, expression = normalize_expression(
        raw_data,
        filter_audit,
        sample_metadata,
        minimum_peptides=minimum_peptides,
        median_center=median_center,
        mix_reference=mix_reference,
    )
    model = fit_protein_models(
        normalized,
        expression,
        sample_metadata,
        trend=trend,
        robust=robust,
    )
    return _sensitivity_frame(name, model.contrast_results)


def _sensitivity_frame(name: str, contrast_results: pl.DataFrame) -> pl.DataFrame:
    return contrast_results.select(
        pl.lit(name).alias("scenario"),
        "raw_row",
        *TEXT_COLUMNS,
        "contrast",
        "adjusted_log2_fold_change",
        "p_value",
        "q_value",
        "discovery",
    )


if __name__ == "__main__":
    raise SystemExit(main())
