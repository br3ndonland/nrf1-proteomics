# pyright: reportAny=false, reportExplicitAny=false, reportMissingTypeStubs=false, reportUnknownArgumentType=false, reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnusedCallResult=false
import json
import math
import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np
import polars as pl

from nrf1_proteomics.analysis import (
    CONTRASTS,
    DEFAULT_OUTPUT_DIR,
    FDR_THRESHOLD,
    AnalysisArtifacts,
    adjust_bh,
    analyze_raw_data,
    build_design_matrix,
    build_filter_audit,
    build_sample_metadata,
    fit_protein_models,
    normalize_expression,
    read_raw_data,
    write_analysis,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DATA = PROJECT_ROOT / "data" / "nrf1-proteomics-raw.csv"


@lru_cache(maxsize=1)
def _artifacts() -> AnalysisArtifacts:
    return analyze_raw_data(RAW_DATA)


def test_default_output_directory_is_temporary() -> None:
    assert DEFAULT_OUTPUT_DIR == (
        Path(tempfile.gettempdir()) / "nrf1-proteomics" / "derived"
    )
    assert not DEFAULT_OUTPUT_DIR.is_relative_to(PROJECT_ROOT)


def test_raw_schema_and_sample_mapping() -> None:
    raw = read_raw_data(RAW_DATA)
    metadata = build_sample_metadata(RAW_DATA)

    assert raw.shape == (633, 26)
    assert raw.get_column("protein_id").n_unique() == 633
    assert metadata.shape == (18, 10)
    assert metadata.group_by("treatment").len().sort("treatment").rows() == [
        ("ha_bort", 4),
        ("ha_chol", 4),
        ("ha_chow", 4),
        ("lacz", 6),
    ]
    assert metadata.group_by("treatment", "run").len().height == 8
    assert set(metadata.get_column("sample")) < set(raw.columns)
    assert metadata.get_column("sample_label")[0] == "Lac Z 1"
    assert metadata.get_column("reporter_tag")[0] == "129N"
    assert metadata.get_column("tmt_run_id").n_unique() == 2


def test_sample_labels_are_read_from_export_header(tmp_path: Path) -> None:
    header_path = tmp_path / "renamed-sample.csv"
    header = "\n".join(RAW_DATA.read_text(encoding="utf-8-sig").splitlines()[:4])
    _ = header_path.write_text(
        header.replace("Lac Z 1", "LAC-Z 1") + "\n",
        encoding="utf-8-sig",
    )

    metadata = build_sample_metadata(header_path)

    assert metadata.get_column("sample_label")[0] == "LAC-Z 1"
    assert metadata.get_column("sample")[0] == "lacz_1"


def test_filter_audit_records_primary_and_sensitivity_sets() -> None:
    raw = read_raw_data(RAW_DATA)
    audit = build_filter_audit(raw)

    assert audit.height == 633
    assert audit.filter(pl.col("included_primary")).height == 252
    assert audit.filter(pl.col("included_one_peptide_sensitivity")).height == 386
    assert audit.filter(
        pl.col("included_primary") & ~pl.col("all_biological_signals_positive")
    ).is_empty()
    assert audit.filter(pl.col("primary_exclusion_reason") == "").is_empty()


def test_normalization_equalizes_channel_medians() -> None:
    raw = read_raw_data(RAW_DATA)
    audit = build_filter_audit(raw)
    normalized, expression = normalize_expression(raw, audit)

    assert normalized.shape == (252, 22)
    assert np.all(np.isfinite(expression))
    assert np.allclose(np.median(expression, axis=0), np.median(expression), atol=1e-12)


def test_primary_model_matches_regression_baseline() -> None:
    artifacts = _artifacts()
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

    assert artifacts.contrast_results.shape == (756, 25)
    assert artifacts.omnibus_results.shape == (252, 10)
    assert artifacts.diagnostics.shape == (252, 19)
    assert discovery_counts == {
        "ha_chol_vs_ha_chow": 0,
        "ha_bort_vs_ha_chow": 47,
        "ha_chow_vs_lacz": 0,
    }
    assert treat_counts == {contrast.key: 0 for contrast in CONTRASTS}
    assert artifacts.omnibus_results.filter(pl.col("discovery")).height == 46


def test_sensitivity_discovery_counts_match_regression_baseline() -> None:
    counts = (
        _artifacts()
        .sensitivities.group_by("scenario", "contrast")
        .agg(pl.col("discovery").sum().alias("discoveries"))
    )
    nonzero = {
        (str(scenario), str(contrast)): int(discoveries)
        for scenario, contrast, discoveries in counts.filter(
            pl.col("discoveries") > 0
        ).iter_rows()
    }

    assert counts.height == 15
    assert nonzero == {
        ("mix_reference_ratios", "ha_bort_vs_ha_chow"): 49,
        ("mix_reference_ratios", "ha_chol_vs_ha_chow"): 1,
        ("no_median_centering", "ha_bort_vs_ha_chow"): 39,
        ("one_peptide_positive", "ha_bort_vs_ha_chow"): 40,
        ("primary", "ha_bort_vs_ha_chow"): 47,
        ("standard_ebayes", "ha_bort_vs_ha_chow"): 46,
    }


def test_c1q_effects_are_nominal_not_bh_discoveries() -> None:
    results = _artifacts().contrast_results.filter(
        (pl.col("contrast") == "ha_chol_vs_ha_chow")
        & pl.col("gene").is_in(["C1qa", "C1qb", "C1qc"])
    )
    effects = dict(results.select("gene", "adjusted_log2_fold_change").iter_rows())

    assert results.height == 3
    assert math.isclose(effects["C1qa"], 1.1216268973, abs_tol=1e-9)
    assert math.isclose(effects["C1qb"], 0.9902542863, abs_tol=1e-9)
    assert math.isclose(effects["C1qc"], 1.2703860284, abs_tol=1e-9)
    assert results.filter(pl.col("p_value") < 0.05).height == 3
    assert results.filter(pl.col("q_value") < FDR_THRESHOLD).is_empty()


def test_bh_adjustment_uses_one_vector_at_a_time() -> None:
    actual = adjust_bh([0.01, 0.04, 0.03, 0.002])

    assert np.allclose(actual, [0.02, 0.04, 0.04, 0.008])


def test_design_and_contrast_directions_with_known_effects() -> None:
    metadata = build_sample_metadata()
    design = build_design_matrix(metadata)
    rng = np.random.default_rng(20260719)
    protein_count = 30
    expression = np.zeros((protein_count, metadata.height))
    treatment = metadata.get_column("treatment").to_numpy()
    run = metadata.get_column("run").to_numpy()
    sample_columns = metadata.get_column("sample").to_list()
    for protein_index in range(protein_count):
        expression[protein_index] = (
            10
            + protein_index / 20
            + 0.7 * (run == "run_2")
            + rng.normal(0, 0.08 + protein_index / 1000, metadata.height)
        )
    expression[0] += (
        0.5 * (treatment == "ha_chow")
        + 2.5 * (treatment == "ha_chol")
        - 0.5 * (treatment == "ha_bort")
    )
    normalized = pl.DataFrame(
        {
            "raw_row": np.arange(1, protein_count + 1),
            "protein_id": [f"protein_{index}" for index in range(protein_count)],
            "gene": [f"gene_{index}" for index in range(protein_count)],
            "description": [""] * protein_count,
            **{
                sample: expression[:, sample_index]
                for sample_index, sample in enumerate(sample_columns)
            },
        }
    )
    model = fit_protein_models(normalized, expression, metadata)
    known = model.contrast_results.filter(pl.col("protein_id") == "protein_0")
    effects = dict(known.select("contrast", "adjusted_log2_fold_change").iter_rows())

    assert design.shape == (18, 5)
    assert np.linalg.matrix_rank(design.to_numpy()) == 5
    assert 1.7 < effects["ha_chol_vs_ha_chow"] < 2.3
    assert -1.3 < effects["ha_bort_vs_ha_chow"] < -0.7
    assert 0.3 < effects["ha_chow_vs_lacz"] < 0.7


def test_diagnostics_never_change_model_rows_or_discoveries() -> None:
    artifacts = _artifacts()
    protein_ids = set(artifacts.normalized_data.get_column("protein_id"))

    assert artifacts.diagnostics.filter(pl.col("diagnostic_warning")).height > 0
    for contrast in CONTRASTS:
        results = artifacts.contrast_results.filter(pl.col("contrast") == contrast.key)
        assert results.height == 252
        assert set(results.get_column("protein_id")) == protein_ids
        assert results.get_column("diagnostic_warning").null_count() == 0


def test_write_analysis_creates_all_public_outputs(tmp_path: Path) -> None:
    artifacts = write_analysis(output_dir=tmp_path, raw_path=RAW_DATA)
    expected = {
        "nrf1-proteomics-sample_metadata.csv",
        "nrf1-proteomics-normalized.csv",
        "nrf1-proteomics-contrasts.csv",
        "nrf1-proteomics-omnibus.csv",
        "nrf1-proteomics-diagnostics.csv",
        "nrf1-proteomics-filter_audit.csv",
        "nrf1-proteomics-sensitivities.csv",
        "nrf1-proteomics-summary.json",
    }
    summary = json.loads(
        (tmp_path / "nrf1-proteomics-summary.json").read_text(encoding="utf-8")
    )

    assert {path.name for path in tmp_path.iterdir()} == expected
    assert artifacts.normalized_data.height == 252
    assert summary["primary_discoveries"]["ha_bort_vs_ha_chow"] == 47
