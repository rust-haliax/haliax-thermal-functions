use cyphus_integration::prelude::*;
use lazy_static::lazy_static;
//use cyphus_specfun::bessel::CylBesselK;
use std::f64::consts::PI;

lazy_static! {
    static ref GK: GaussKronrodIntegrator<f64> = GaussKronrodIntegratorBuilder::default()
        .epsabs(0.0)
        .epsrel(1e-3)
        .order(7)
        .build();
}

fn neq_scaled(x: f64, spin2: usize) -> f64 {
    let eta = if spin2 % 2 == 0 { 1.0 } else { -1.0 };
    let integrand = |z: f64| z * (z * z - x * x).sqrt() / (z.exp() - eta);
    (*GK).integrate(integrand, x, f64::INFINITY).val / (2.0 * PI * PI)
}

fn energy_density_scaled(x: f64, spin2: usize) -> f64 {
    let eta = if spin2 % 2 == 0 { 1.0 } else { -1.0 };
    let integrand = |z: f64| z * z * (z * z - x * x).sqrt() / (z.exp() - eta);
    (*GK).integrate(integrand, x, f64::INFINITY).val / (2.0 * PI * PI)
}

fn pressure_density_scaled(x: f64, spin2: usize) -> f64 {
    let eta = if spin2 % 2 == 0 { 1.0 } else { -1.0 };
    let integrand = |z: f64| (z * z - x * x).powf(1.5) / (z.exp() - eta);
    (*GK).integrate(integrand, x, f64::INFINITY).val / (6.0 * PI * PI)
}

fn entropy_density_scaled(x: f64, spin2: usize) -> f64 {
    let eta = if spin2 % 2 == 0 { 1.0 } else { -1.0 };
    let integrand = |z: f64| (4.0 * z * z - x * x) * (z * z - x * x).sqrt() / (z.exp() - eta);
    (*GK).integrate(integrand, x, f64::INFINITY).val / (6.0 * PI * PI)
}

/// Compute the equilibrium number density of a particle with a given
/// `temperature`, `mass`, internal degrees of freedom `g` and spin `spin2/2`.
pub fn neq(temperature: f64, mass: f64, g: f64, spin2: usize) -> f64 {
    g * temperature.powi(3) * neq_scaled(mass / temperature, spin2)
}

/// Compute the equilibrium energy density of a particle with a given
/// `temperature`, `mass`, internal degrees of freedom `g` and spin `spin2/2`.
pub fn energy_density(temperature: f64, mass: f64, g: f64, spin2: usize) -> f64 {
    g * temperature.powi(4) * energy_density_scaled(mass / temperature, spin2)
}

/// Compute the equilibrium pressure density of a particle with a given
/// `temperature`, `mass`, internal degrees of freedom `g` and spin `spin2/2`.
pub fn pressure_density(temperature: f64, mass: f64, g: f64, spin2: usize) -> f64 {
    g * temperature.powi(4) * pressure_density_scaled(mass / temperature, spin2)
}

/// Compute the equilibrium entropy density of a particle with a given
/// `temperature`, `mass`, internal degrees of freedom `g` and spin `spin2/2`.
pub fn entropy_density(temperature: f64, mass: f64, g: f64, spin2: usize) -> f64 {
    g * temperature.powi(3) * entropy_density_scaled(mass / temperature, spin2)
}

/// Compute the equilibrium number of degrees of freedom stored in energy of a
/// particle with a given `temperature`, `mass`, internal degrees of freedom
/// `g` and spin `spin2/2`.
pub fn geff(temperature: f64, mass: f64, g: f64, spin2: usize) -> f64 {
    30.0 / (PI * PI) * g * energy_density_scaled(mass / temperature, spin2)
}

/// Compute the equilibrium number of degrees of freedom stored in entropy of a
/// particle with a given `temperature`, `mass`, internal degrees of freedom
/// `g` and spin `spin2/2`.
pub fn heff(temperature: f64, mass: f64, g: f64, spin2: usize) -> f64 {
    45.0 / (2.0 * PI * PI) * g * entropy_density_scaled(mass / temperature, spin2)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_neq_scaled() {
        let x = 100.0;
        let spin2 = 2;
        println!("res = {:?}", neq_scaled(x, spin2));
    }
    #[test]
    fn test_energy_density_scaled() {
        let x = 10.0;
        let spin2 = 2;
        println!("res = {:?}", energy_density_scaled(x, spin2));
    }
}
