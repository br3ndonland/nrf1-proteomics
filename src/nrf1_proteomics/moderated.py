"""Protein-wise linear models with empirical Bayes variance moderation.

This module implements linear models useful for proteomic analysis:

- ordinary least-squares contrasts
- Smyth-style scaled-F variance shrinkage
- an abundance-dependent variance trend
- robust protection for hypervariable proteins
- moderated omnibus F tests
- TREAT interval-null tests.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TypeAlias, cast

import numpy as np
import scipy.linalg
import scipy.optimize
import scipy.special
import scipy.stats
from numpy.typing import NDArray
from scipy.interpolate import BSpline

_WINSOR_LOWER = 0.05
_WINSOR_UPPER = 0.10
_LOWESS_FRACTION = 0.40
FloatArray: TypeAlias = NDArray[np.float64]
_GAUSS_NODES, _GAUSS_WEIGHTS = cast(
    tuple[FloatArray, FloatArray],
    np.polynomial.legendre.leggauss(128),
)


@dataclass(frozen=True)
class ModeratedModelFit:
    """Statistics from one protein-wise model and its contrast matrix."""

    coefficients: FloatArray
    stdev_unscaled: FloatArray
    residual_variance: FloatArray
    prior_variance: FloatArray
    prior_degrees_of_freedom: FloatArray
    posterior_variance: FloatArray
    total_degrees_of_freedom: FloatArray
    moderated_t: FloatArray
    p_value: FloatArray
    moderated_f: FloatArray
    f_p_value: FloatArray
    treat_t: FloatArray
    treat_p_value: FloatArray


def fit_moderated_model(
    expression: FloatArray,
    design: FloatArray,
    contrasts: FloatArray,
    *,
    trend: bool = True,
    robust: bool = True,
    treat_lfc: float = 0.0,
) -> ModeratedModelFit:
    """Fit blocked OLS models and return moderated contrast statistics."""

    expression = np.asarray(expression, dtype=np.float64)
    design = np.asarray(design, dtype=np.float64)
    contrasts = np.asarray(contrasts, dtype=np.float64)
    _validate_model_inputs(expression, design, contrasts)

    protein_count = cast(int, expression.shape[0])
    sample_count = cast(int, expression.shape[1])
    design_rank = cast(int, np.linalg.matrix_rank(design))
    residual_df = float(sample_count - design_rank)
    coefficient_covariance = np.linalg.inv(design.T @ design)
    base_coefficients = expression @ design @ coefficient_covariance
    residuals = expression - base_coefficients @ design.T
    residual_variance = np.sum(residuals**2, axis=1) / residual_df

    coefficients = base_coefficients @ contrasts
    contrast_covariance = contrasts.T @ coefficient_covariance @ contrasts
    unscaled = np.sqrt(np.diag(contrast_covariance))
    stdev_unscaled = np.broadcast_to(
        unscaled,
        (protein_count, contrasts.shape[1]),
    ).copy()

    covariate = np.mean(expression, axis=1) if trend else None
    prior_variance, prior_df = _estimate_variance_prior(
        residual_variance,
        residual_df,
        covariate=covariate,
        robust=robust,
    )
    posterior_variance = _posterior_variance(
        residual_variance,
        residual_df,
        prior_variance,
        prior_df,
    )
    pooled_df = residual_df * protein_count
    total_df = np.minimum(residual_df + prior_df, pooled_df)
    standard_error = stdev_unscaled * np.sqrt(posterior_variance)[:, np.newaxis]
    moderated_t = coefficients / standard_error
    p_value = np.asarray(
        2.0
        * scipy.stats.t.sf(
            np.abs(moderated_t),
            total_df[:, np.newaxis],
        ),
        dtype=np.float64,
    )

    contrast_rank = cast(int, np.linalg.matrix_rank(contrasts))
    contrast_precision = np.linalg.pinv(contrast_covariance)
    quadratic = cast(
        FloatArray,
        np.einsum(
            "gi,ij,gj->g",
            coefficients,
            contrast_precision,
            coefficients,
        ),
    )
    moderated_f = quadratic / (contrast_rank * posterior_variance)
    f_p_value = np.asarray(
        scipy.stats.f.sf(moderated_f, contrast_rank, total_df),
        dtype=np.float64,
    )

    absolute_coefficients = np.abs(coefficients)
    lower_t = (absolute_coefficients - treat_lfc) / standard_error
    upper_t = (absolute_coefficients + treat_lfc) / standard_error
    treat_t = np.sign(coefficients) * np.maximum(lower_t, 0.0)
    treat_p_value = np.asarray(
        np.minimum(
            scipy.stats.t.sf(lower_t, total_df[:, np.newaxis])
            + scipy.stats.t.sf(upper_t, total_df[:, np.newaxis]),
            1.0,
        ),
        dtype=np.float64,
    )
    return ModeratedModelFit(
        coefficients=coefficients,
        stdev_unscaled=stdev_unscaled,
        residual_variance=residual_variance,
        prior_variance=prior_variance,
        prior_degrees_of_freedom=prior_df,
        posterior_variance=posterior_variance,
        total_degrees_of_freedom=total_df,
        moderated_t=moderated_t,
        p_value=p_value,
        moderated_f=moderated_f,
        f_p_value=f_p_value,
        treat_t=treat_t,
        treat_p_value=treat_p_value,
    )


def _validate_model_inputs(
    expression: FloatArray,
    design: FloatArray,
    contrasts: FloatArray,
) -> None:
    if expression.ndim != 2 or design.ndim != 2 or contrasts.ndim != 2:
        raise ValueError("Expression, design, and contrasts must be matrices.")
    if expression.shape[1] != design.shape[0]:
        raise ValueError("Expression samples must match design rows.")
    if design.shape[1] != contrasts.shape[0]:
        raise ValueError("Design coefficients must match contrast rows.")
    if design.shape[0] <= design.shape[1]:
        raise ValueError("The model must retain residual degrees of freedom.")
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("The design matrix must be full rank.")
    if np.linalg.matrix_rank(contrasts) != contrasts.shape[1]:
        raise ValueError("The contrast matrix must be full rank.")
    if not all(
        np.all(np.isfinite(matrix)) for matrix in (expression, design, contrasts)
    ):
        raise ValueError("Model inputs must contain only finite values.")


def _estimate_variance_prior(
    variances: FloatArray,
    residual_df: float,
    *,
    covariate: FloatArray | None,
    robust: bool,
) -> tuple[FloatArray, FloatArray]:
    """Estimate scaled-F prior variances and degrees of freedom."""
    adjusted_log_variance = (
        np.log(variances)
        - scipy.special.digamma(residual_df / 2.0)
        + math.log(residual_df / 2.0)
    )
    if covariate is None:
        fitted_log_variance = np.full(
            variances.size,
            float(np.mean(adjusted_log_variance)),
        )
        dispersion = float(np.var(adjusted_log_variance, ddof=1))
    else:
        fitted_log_variance, dispersion = _fit_abundance_trend(
            adjusted_log_variance,
            covariate,
        )

    excess_dispersion = dispersion - float(
        scipy.special.polygamma(1, residual_df / 2.0)
    )
    if excess_dispersion <= 0:
        baseline_df = math.inf
        baseline_variance = (
            np.full(variances.size, float(np.mean(variances)))
            if covariate is None
            else np.exp(fitted_log_variance)
        )
    else:
        baseline_df = 2.0 * _inverse_trigamma(excess_dispersion)
        log_scale_adjustment = scipy.special.digamma(baseline_df / 2.0) - math.log(
            baseline_df / 2.0
        )
        baseline_variance = np.exp(fitted_log_variance + log_scale_adjustment)

    if not robust or variances.size < 3:
        return (
            baseline_variance,
            np.full(variances.size, baseline_df),
        )
    return _robust_variance_prior(
        variances,
        residual_df,
        covariate,
        baseline_df,
        baseline_variance,
    )


def _fit_abundance_trend(
    adjusted_log_variance: FloatArray,
    covariate: FloatArray,
) -> tuple[FloatArray, float]:
    spline_df = (
        1
        + int(covariate.size >= 3)
        + int(covariate.size >= 6)
        + int(covariate.size >= 30)
    )
    spline_df = min(spline_df, np.unique(covariate).size)
    if spline_df < 2:
        fitted = np.full(
            adjusted_log_variance.size,
            float(np.mean(adjusted_log_variance)),
        )
        return fitted, float(np.var(adjusted_log_variance, ddof=1))

    spline = _natural_spline_design(covariate, spline_df)
    coefficients = np.linalg.lstsq(
        spline,
        adjusted_log_variance,
        rcond=None,
    )[0]
    fitted = spline @ coefficients
    rank = cast(int, np.linalg.matrix_rank(spline))
    residuals = adjusted_log_variance - fitted
    dispersion = float(np.sum(residuals**2) / (covariate.size - rank))
    return fitted, dispersion


def _natural_spline_design(values: FloatArray, degrees: int) -> FloatArray:
    values = np.asarray(values, dtype=np.float64)
    lower = float(np.min(values))
    upper = float(np.max(values))
    if lower == upper:
        return np.ones((values.size, 1))

    interior_count = max(degrees - 2, 0)
    probabilities = np.asarray(
        np.linspace(0.0, 1.0, interior_count + 2)[1:-1],
        dtype=np.float64,
    )
    interior = np.asarray(np.quantile(values, probabilities), dtype=np.float64)
    knots = np.asarray(
        np.concatenate(
            (
                np.repeat(lower, 4),
                interior,
                np.repeat(upper, 4),
            )
        ),
        dtype=np.float64,
    )
    basis_count = int(knots.size - 4)
    identity = np.asarray(np.eye(basis_count), dtype=np.float64)
    raw_basis = np.asarray(
        np.column_stack(
            [
                np.asarray(
                    BSpline(
                        knots,
                        cast(FloatArray, identity[index]),
                        3,
                    )(values),
                    dtype=np.float64,
                )
                for index in range(basis_count)
            ]
        ),
        dtype=np.float64,
    )
    constraints = np.asarray(
        np.column_stack(
            [
                np.asarray(
                    BSpline(
                        knots,
                        cast(FloatArray, identity[index]),
                        3,
                    ).derivative(2)([lower, upper]),
                    dtype=np.float64,
                )
                for index in range(basis_count)
            ]
        ),
        dtype=np.float64,
    )
    natural_space = np.asarray(
        scipy.linalg.null_space(constraints),
        dtype=np.float64,
    )
    return np.asarray(raw_basis @ natural_space, dtype=np.float64)


def _inverse_trigamma(value: float) -> float:
    if value <= 0:
        return math.inf
    lower = 1e-8
    upper = max(1.0, 2.0 / value)
    while float(scipy.special.polygamma(1, upper)) > value:
        upper *= 2.0
    return float(
        scipy.optimize.brentq(
            lambda candidate: float(scipy.special.polygamma(1, candidate)) - value,
            lower,
            upper,
        )
    )


def _robust_variance_prior(
    variances: FloatArray,
    residual_df: float,
    covariate: FloatArray | None,
    baseline_df: float,
    baseline_variance: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    """Winsorize log variances and reduce shrinkage for upper-tail outliers."""
    log_variances = np.log(variances)
    if covariate is None:
        sorted_values = np.sort(log_variances)
        trim = math.floor(_WINSOR_UPPER * variances.size)
        center = float(np.mean(sorted_values[trim : variances.size - trim]))
        log_trend = np.full(variances.size, center)
    else:
        log_trend = _lowess(log_variances, covariate)

    residuals = log_variances - log_trend
    quantiles: FloatArray = np.quantile(
        residuals,
        [_WINSOR_LOWER, 1.0 - _WINSOR_UPPER],
    )
    lower = quantiles.item(0)
    upper = quantiles.item(1)
    winsorized = np.clip(residuals, lower, upper)
    observed_mean = float(np.mean(winsorized))
    observed_variance = float(np.var(winsorized, ddof=1))

    robust_df = baseline_df
    if math.isfinite(baseline_df):
        baseline_moments = _winsorized_log_f_moments(
            residual_df,
            baseline_df,
        )
        infinite_moments = _winsorized_log_f_moments(
            residual_df,
            math.inf,
        )
        if observed_variance <= infinite_moments[1]:
            robust_df = math.inf
        elif observed_variance < baseline_moments[1]:
            robust_df = float(
                scipy.optimize.brentq(
                    lambda candidate: (
                        _winsorized_log_f_moments(
                            residual_df,
                            candidate,
                        )[1]
                        - observed_variance
                    ),
                    baseline_df,
                    1e8,
                )
            )

    theoretical_mean = _winsorized_log_f_moments(
        residual_df,
        robust_df,
    )[0]
    prior_variance = np.exp(log_trend + observed_mean - theoretical_mean)
    if not np.all(np.isfinite(prior_variance)):
        prior_variance = baseline_variance

    variance_ratio = variances / prior_variance
    if math.isfinite(robust_df):
        upper_tail = np.asarray(
            scipy.stats.f.sf(
                variance_ratio,
                residual_df,
                robust_df,
            ),
            dtype=np.float64,
        )
    else:
        upper_tail = np.asarray(
            scipy.stats.chi2.sf(
                variance_ratio * residual_df,
                residual_df,
            ),
            dtype=np.float64,
        )
    ranks = np.argsort(np.argsort(variance_ratio)[::-1]) + 1
    empirical_tail = (ranks - 0.5) / variances.size
    probability_regular = np.minimum(upper_tail / empirical_tail, 1.0)
    if np.all(probability_regular >= 1.0):
        return prior_variance, np.full(variances.size, robust_df)

    if math.isfinite(robust_df):
        outlier_df = _outlier_prior_df(
            variance_ratio,
            residual_df,
            robust_df,
            upper_tail,
        )
        prior_df = (
            probability_regular * robust_df + (1.0 - probability_regular) * outlier_df
        )
    else:
        pooled_df = residual_df * variances.size
        prior_df = probability_regular * pooled_df
    return prior_variance, _monotone_prior_df(prior_df, upper_tail)


def _winsorized_log_f_moments(
    numerator_df: float,
    denominator_df: float,
) -> tuple[float, float]:
    left = _WINSOR_LOWER
    right = 1.0 - _WINSOR_UPPER
    probabilities = (_GAUSS_NODES + 1.0) * (right - left) / 2.0 + left
    weights = _GAUSS_WEIGHTS * (right - left) / 2.0
    if math.isfinite(denominator_df):
        quantiles = np.asarray(
            scipy.stats.f.ppf(
                probabilities,
                numerator_df,
                denominator_df,
            ),
            dtype=np.float64,
        )
        lower = float(scipy.stats.f.ppf(left, numerator_df, denominator_df))
        upper = float(scipy.stats.f.ppf(right, numerator_df, denominator_df))
    else:
        quantiles = (
            np.asarray(
                scipy.stats.chi2.ppf(probabilities, numerator_df),
                dtype=np.float64,
            )
            / numerator_df
        )
        lower = float(scipy.stats.chi2.ppf(left, numerator_df) / numerator_df)
        upper = float(scipy.stats.chi2.ppf(right, numerator_df) / numerator_df)
    log_quantiles = np.log(quantiles)
    log_lower = math.log(lower)
    log_upper = math.log(upper)
    mean = (
        _WINSOR_LOWER * log_lower
        + float(np.sum(weights * log_quantiles))
        + _WINSOR_UPPER * log_upper
    )
    second_moment = (
        _WINSOR_LOWER * log_lower**2
        + float(np.sum(weights * log_quantiles**2))
        + _WINSOR_UPPER * log_upper**2
    )
    return mean, max(second_moment - mean**2, 0.0)


def _outlier_prior_df(
    variance_ratio: FloatArray,
    residual_df: float,
    baseline_df: float,
    upper_tail: FloatArray,
) -> float:
    most_extreme = int(np.argmax(variance_ratio))
    extreme_ratio = variance_ratio.item(most_extreme)
    log_tail = float(
        scipy.stats.f.logsf(
            extreme_ratio,
            residual_df,
            baseline_df,
        )
    )
    if not math.isfinite(log_tail):
        return 0.0
    candidate = math.log(0.5) / log_tail * baseline_df
    refined_log_tail = float(
        scipy.stats.f.logsf(
            extreme_ratio,
            residual_df,
            candidate,
        )
    )
    if math.isfinite(refined_log_tail) and refined_log_tail != 0:
        candidate *= math.log(0.5) / refined_log_tail
    if cast(float, np.min(upper_tail)) == 0.0:
        return 0.0
    return max(candidate, 0.0)


def _monotone_prior_df(
    prior_df: FloatArray,
    upper_tail: FloatArray,
) -> FloatArray:
    order = np.argsort(upper_tail)
    ordered = prior_df[order].copy()
    cumulative_means = np.cumsum(ordered) / np.arange(1, ordered.size + 1)
    minimum_index = int(np.argmin(cumulative_means))
    ordered[: minimum_index + 1] = cumulative_means[minimum_index]
    ordered = np.maximum.accumulate(ordered)
    result = np.empty_like(ordered)
    result[order] = ordered
    return result


def _lowess(
    response: FloatArray,
    covariate: FloatArray,
    *,
    fraction: float = _LOWESS_FRACTION,
    iterations: int = 3,
) -> FloatArray:
    order = np.asarray(np.argsort(covariate, kind="stable"), dtype=np.intp)
    x = np.asarray(covariate, dtype=np.float64)[order]
    y = np.asarray(response, dtype=np.float64)[order]
    point_count = int(x.size)
    neighborhood = max(2, math.ceil(fraction * point_count))
    robust_weights: FloatArray = np.ones(point_count, dtype=np.float64)
    fitted: FloatArray = np.empty(point_count, dtype=np.float64)

    for iteration in range(iterations + 1):
        for index in range(point_count):
            center = x.item(index)
            distances: FloatArray = np.abs(x - center)
            bandwidth = np.partition(distances, neighborhood - 1).item(neighborhood - 1)
            if bandwidth == 0:
                local_weights: FloatArray = np.asarray(
                    np.isclose(distances, 0.0, rtol=0.0, atol=0.0),
                    dtype=np.float64,
                )
            else:
                local_weights = np.asarray(
                    np.clip(
                        1.0 - (distances / bandwidth) ** 3,
                        0.0,
                        None,
                    )
                    ** 3,
                    dtype=np.float64,
                )
            local_weights *= robust_weights
            centered: FloatArray = x - center
            weight_sum = float(np.sum(local_weights))
            weighted_x = float(np.sum(local_weights * centered))
            weighted_xx = float(np.sum(local_weights * centered**2))
            weighted_y = float(np.sum(local_weights * y))
            weighted_xy = float(np.sum(local_weights * centered * y))
            denominator = weight_sum * weighted_xx - weighted_x**2
            fitted[index] = (
                weighted_y / weight_sum
                if denominator <= 1e-20
                else (weighted_xx * weighted_y - weighted_x * weighted_xy) / denominator
            )
        if iteration == iterations:
            break
        residuals = y - fitted
        scale = 6.0 * cast(float, np.median(np.abs(residuals)))
        robust_weights = (
            np.clip(1.0 - (residuals / scale) ** 2, 0.0, None) ** 2
            if scale > 0
            else np.ones(point_count)
        )

    result: FloatArray = np.empty(point_count, dtype=np.float64)
    result[order] = fitted
    return result


def _posterior_variance(
    residual_variance: FloatArray,
    residual_df: float,
    prior_variance: FloatArray,
    prior_df: FloatArray,
) -> FloatArray:
    finite = np.isfinite(prior_df)
    posterior = prior_variance.copy()
    posterior[finite] = (
        prior_df[finite] * prior_variance[finite]
        + residual_df * residual_variance[finite]
    ) / (prior_df[finite] + residual_df)
    return posterior
