pub(super) mod data;

use cyphus_interpolation::prelude::*;
use data::*;
use lazy_static::lazy_static;
use ndarray::prelude::*;

const SM_LOG_TEMP_MIN: f64 = -4.5;
const SM_LOG_TEMP_MAX: f64 = 4.0;

const SM_GSTAR_FRONT: f64 = 2.141289997868463;
const SM_GSTAR_BACK: f64 = 10.3359;

const SM_GEFF_FRONT: f64 = 3.3839699989395835;
const SM_GEFF_BACK: f64 = 106.83;

const SM_HEFF_FRONT: f64 = 3.9387999991430975;
const SM_HEFF_BACK: f64 = 106.83;

lazy_static! {
    static ref LOG_TEMP_DATA: Array1<f64> = Array::linspace(-4.5, 4.0, 341);
    static ref SM_SPLINE_SQRT_GSTAR: UnivariateSpline =
        UnivariateSplineBuilder::default(&(*LOG_TEMP_DATA), &(*GSTAR_DATA))
            .smoothing_factor(0.0)
            .degree(3)
            .extrapolation(3)
            .build()
            .unwrap();
    static ref SM_SPLINE_HEFF: UnivariateSpline = {
        UnivariateSplineBuilder::default(&(*LOG_TEMP_DATA), &(*HEFF_DATA))
            .smoothing_factor(0.0)
            .degree(3)
            .extrapolation(3)
            .build()
            .unwrap()
    };
    static ref SM_SPLINE_GEFF: UnivariateSpline = {
        UnivariateSplineBuilder::default(&(*LOG_TEMP_DATA), &(*GEFF_DATA))
            .smoothing_factor(0.0)
            .degree(3)
            .build()
            .unwrap()
    };
}

/// Compute the effective number of degrees of freedom in energy stored in
/// the standard model bath at a temperature `temp`.
pub fn sm_geff(temp: f64) -> f64 {
    let log_temp = temp.log10();
    if SM_LOG_TEMP_MIN < log_temp && log_temp < SM_LOG_TEMP_MAX {
        (*SM_SPLINE_GEFF).eval(log_temp)
    } else if log_temp < SM_LOG_TEMP_MIN {
        SM_GEFF_FRONT
    } else {
        SM_GEFF_BACK
    }
}

/// Compute the derivative w.r.t temperature of the effective number of
/// degrees of freedom in entropy stored in the standard model bath at a
/// temperature `temp`.
pub fn sm_geff_deriv(temp: f64) -> f64 {
    let log_temp = temp.log10();
    if SM_LOG_TEMP_MIN < log_temp && log_temp < SM_LOG_TEMP_MAX {
        (*SM_SPLINE_GEFF).derivative(1, log_temp) / temp
    } else {
        0.0
    }
}

/// Compute the effective number of degrees of freedom in entropy stored in
/// the standard model bath at a temperature `temp`.
pub fn sm_heff(temp: f64) -> f64 {
    let log_temp = temp.log10();
    if SM_LOG_TEMP_MIN < log_temp && log_temp < SM_LOG_TEMP_MAX {
        (*SM_SPLINE_HEFF).eval(log_temp)
    } else if log_temp < SM_LOG_TEMP_MIN {
        SM_HEFF_FRONT
    } else {
        SM_HEFF_BACK
    }
}

/// Compute the derivative w.r.t temperature of the effective number of
/// degrees of freedom in entropy stored in the standard model bath at a
/// temperature `temp`.
pub fn sm_heff_deriv(temp: f64) -> f64 {
    let log_temp = temp.log10();
    if SM_LOG_TEMP_MIN < log_temp && log_temp < SM_LOG_TEMP_MAX {
        (*SM_SPLINE_HEFF).derivative(1, log_temp) / temp
    } else {
        0.0
    }
}

/// Compute the square-root of g-star, which is given by:
///     g^{1/2}_*(T) = (1 + T/(3h) dh/dT) h / g^{1/2}
/// where `h` and `g` are the effective number of degrees of freedom stored in
/// entropy and energy and `T` is the temperature of the SM bath.
pub fn sm_sqrt_gstar(temp: f64) -> f64 {
    let log_temp = temp.log10();
    if SM_LOG_TEMP_MIN < log_temp && log_temp < SM_LOG_TEMP_MAX {
        (*SM_SPLINE_SQRT_GSTAR).eval(log_temp)
    } else if log_temp < SM_LOG_TEMP_MIN {
        SM_GSTAR_FRONT
    } else {
        SM_GSTAR_BACK
    }
}

/// Compute the derivative w.r.t. temperature of `sm_sqrt_gstar`.
pub fn sm_sqrt_gstar_deriv(temp: f64) -> f64 {
    let log_temp = temp.log10();
    if SM_LOG_TEMP_MIN < log_temp && log_temp < SM_LOG_TEMP_MAX {
        (*SM_SPLINE_SQRT_GSTAR).derivative(1, log_temp)
    } else {
        0.0
    }
}

/// Compute the energy density of the SM bath at a given temperature `temp`.
pub fn sm_energy_density(temp: f64) -> f64 {
    std::f64::consts::PI.powi(2) / 30.0 * sm_geff(temp) * temp.powi(4)
}

/// Compute the derivative w.r.t. temperature of `sm_energy_density`.
pub fn sm_energy_density_deriv(temp: f64) -> f64 {
    std::f64::consts::PI.powi(2) / 30.0
        * temp.powi(3)
        * (temp * sm_geff_deriv(temp) + 4.0 * sm_geff(temp))
}

/// Compute the entropy density of the SM bath at a given temperature `temp`.
pub fn sm_entropy_density(temp: f64) -> f64 {
    2.0 * std::f64::consts::PI.powi(2) / 45.0 * sm_heff(temp) * temp.powi(3)
}

/// Compute the derivative w.r.t. temperature of `sm_entropy_density`.
pub fn sm_entropy_density_deriv(temp: f64) -> f64 {
    2.0 * std::f64::consts::PI.powi(2) / 45.0
        * temp.powi(2)
        * (temp * sm_heff_deriv(temp) + 3.0 * sm_heff(temp))
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_sm_geff() {
        let temp = 1.0;
        dbg!(sm_geff(temp));
        assert!(false);
    }
}
