from typing import TypeAlias, cast

import numpy as np
import scipy.special
import scipy.stats
from numpy.typing import NDArray

from nrf1_proteomics.moderated import fit_moderated_model

FloatArray: TypeAlias = NDArray[np.float64]


def _synthetic_model(
    *,
    protein_count: int = 60,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    rng = np.random.default_rng(20260728)
    design: FloatArray = np.asarray(
        np.column_stack(
            (
                np.repeat([1.0, 0.0], 6),
                np.repeat([0.0, 1.0], 6),
            )
        ),
        dtype=np.float64,
    )
    sample_count = cast(int, design.shape[0])
    expression: FloatArray = np.empty(
        (protein_count, sample_count),
        dtype=np.float64,
    )
    for protein in range(protein_count):
        group_means: FloatArray = np.array(
            [8.0 + protein / 50.0, 8.4 + protein / 50.0],
            dtype=np.float64,
        )
        expression[protein] = design @ group_means + rng.normal(
            0.0,
            0.2 + protein / 1000.0,
            sample_count,
        )
    contrasts: FloatArray = np.array([[-1.0], [1.0]], dtype=np.float64)
    return expression, design, contrasts


def test_standard_moderation_matches_closed_form_statistics() -> None:
    expression, design, contrasts = _synthetic_model()

    fit = fit_moderated_model(
        expression,
        design,
        contrasts,
        trend=False,
        robust=False,
        treat_lfc=0.5,
    )

    covariance = np.linalg.inv(design.T @ design)
    base_coefficients = expression @ design @ covariance
    residuals = expression - base_coefficients @ design.T
    residual_df = float(cast(int, design.shape[0]) - cast(int, design.shape[1]))
    residual_variance: FloatArray = np.sum(residuals**2, axis=1) / residual_df
    adjusted_log_variance = cast(
        FloatArray,
        (
            np.log(residual_variance)
            - float(scipy.special.digamma(residual_df / 2))
            + np.log(residual_df / 2)
        ),
    )
    excess_dispersion = float(np.var(adjusted_log_variance, ddof=1)) - float(
        scipy.special.polygamma(
            1,
            residual_df / 2,
        )
    )

    assert np.allclose(fit.coefficients, base_coefficients @ contrasts, atol=1e-12)
    assert np.allclose(fit.residual_variance, residual_variance, atol=1e-12)
    assert bool((fit.prior_degrees_of_freedom > 0).all())
    assert np.isclose(
        float(
            scipy.special.polygamma(
                1,
                fit.prior_degrees_of_freedom.item(0) / 2,
            )
        ),
        excess_dispersion,
        atol=1e-12,
    )
    assert np.allclose(fit.moderated_f, fit.moderated_t[:, 0] ** 2, atol=1e-12)


def test_treat_uses_interval_null_p_value() -> None:
    expression, design, contrasts = _synthetic_model()
    threshold = 0.5

    fit = fit_moderated_model(
        expression,
        design,
        contrasts,
        trend=False,
        robust=False,
        treat_lfc=threshold,
    )

    standard_error = fit.stdev_unscaled * np.sqrt(fit.posterior_variance)[:, np.newaxis]
    lower_t = (np.abs(fit.coefficients) - threshold) / standard_error
    upper_t = (np.abs(fit.coefficients) + threshold) / standard_error
    expected_p = scipy.stats.t.sf(
        lower_t,
        fit.total_degrees_of_freedom[:, np.newaxis],
    ) + scipy.stats.t.sf(
        upper_t,
        fit.total_degrees_of_freedom[:, np.newaxis],
    )

    assert np.allclose(fit.treat_p_value, expected_p, atol=1e-12)
    within_threshold = np.abs(fit.coefficients) <= threshold
    within_threshold_t = np.extract(within_threshold, fit.treat_t)
    assert np.count_nonzero(within_threshold_t) == 0
    assert bool((fit.treat_p_value <= 1).all())


def test_robust_moderation_limits_hypervariable_protein_shrinkage() -> None:
    expression, design, contrasts = _synthetic_model(protein_count=100)
    expression[0] += np.tile([5.0, -5.0], 6)

    fit = fit_moderated_model(
        expression,
        design,
        contrasts,
        trend=True,
        robust=True,
    )

    assert fit.residual_variance[0] == np.max(fit.residual_variance)
    assert fit.prior_degrees_of_freedom[0] < np.median(fit.prior_degrees_of_freedom)
    assert bool(np.isfinite(fit.posterior_variance).all())
    assert bool((fit.posterior_variance > 0).all())
