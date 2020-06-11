pub struct ThermodynamicParticle {
    pub m: f64,
    pub g: f64,
    pub spin2: usize,
}

impl ThermodynamicParticle {
    pub fn neq(&self, temperature: f64) -> f64 {
        crate::thermal_functions::neq(temperature, self.m, self.g, self.spin2)
    }
    pub fn energy_density(&self, temperature: f64) -> f64 {
        crate::thermal_functions::energy_density(temperature, self.m, self.g, self.spin2)
    }
    pub fn pressure_density(&self, temperature: f64) -> f64 {
        crate::thermal_functions::pressure_density(temperature, self.m, self.g, self.spin2)
    }
    pub fn entropy_density(&self, temperature: f64) -> f64 {
        crate::thermal_functions::entropy_density(temperature, self.m, self.g, self.spin2)
    }
    pub fn geff(&self, temperature: f64) -> f64 {
        crate::thermal_functions::geff(temperature, self.m, self.g, self.spin2)
    }
    pub fn heff(&self, temperature: f64) -> f64 {
        crate::thermal_functions::geff(temperature, self.m, self.g, self.spin2)
    }
}
