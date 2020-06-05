use cyphus_integration::GaussKronrodIntegratorBuilder;
use std::f64::consts::PI;

fn neq_scaled(x: f64, spin2: usize) -> f64 {
    let gk = GaussKronrodIntegratorBuilder::default()
        .epsrel(1e-8)
        .epsabs(0.0)
        .key(2)
        .build();

    let eta = if spin2 % 2 == 0 { 1.0 } else { -1.0 };
    let integrand = |z: f64| z * (z * z - x * x).sqrt() / (z.exp() - eta);
    gk.integrate(integrand, x, f64::INFINITY).val / (2.0 * PI * PI)
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
}
