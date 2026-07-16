//! Wolf damped-shifted-force (DSF) Coulomb pair potential for grid building.
//!
//! Reference: Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).

const INV_SQRT_PI: f64 = 0.564_189_583_547_756_3; // 1/sqrt(pi)

/// RASPA internal Coulomb conversion: (q_i q_j)/(4π ε₀) → Kelvin when q in e, r in Å.
pub const COULOMB_CONVERSION_K: f64 = 138_935.483_496_401_7;

/// Precomputed Wolf DSF constants for (alpha, cutoff).
pub struct WolfConstants {
    pub alpha: f64,
    pub alpha_sq: f64,
    pub cutoff: f64,
    pub erfc_a_r: f64,
    pub shift: f64,
    pub slope: f64,
}

impl WolfConstants {
    pub fn new(alpha: f64, cutoff: f64) -> Self {
        let alpha_sq = alpha * alpha;
        let er = libm::erfc(alpha * cutoff);
        let ex = (-alpha_sq * cutoff * cutoff).exp();
        let shift = er / cutoff;
        let slope = shift / cutoff + 2.0 * alpha * ex * INV_SQRT_PI / cutoff;
        Self {
            alpha,
            alpha_sq,
            cutoff,
            erfc_a_r: er,
            shift,
            slope,
        }
    }
}

/// Accumulate Wolf DSF Coulomb contribution into energy and radial derivatives.
///
/// Probe charge is +1 e; `q_j` is the framework atom charge in elementary charge.
/// All values accumulated in Kelvin (RASPA internal energy units).
#[inline]
pub fn accumulate_wolf_dsf(
    w: &WolfConstants,
    q_j: f64,
    r: f64,
    r2: f64,
    energy: &mut f64,
    f1: &mut f64,
    f2: &mut f64,
    f3: &mut f64,
) {
    if r >= w.cutoff {
        return;
    }

    let scale = COULOMB_CONVERSION_K * q_j;
    let erfc_ar = libm::erfc(w.alpha * r);
    let exp_ar2 = (-w.alpha_sq * r2).exp();

    *energy += scale * (erfc_ar / r - w.shift + w.slope * (r - w.cutoff));

    let inv_r2 = 1.0 / r2;
    let inv_r3 = inv_r2 / r;
    let inv_r4 = inv_r2 * inv_r2;
    let inv_r5 = inv_r4 / r;
    let inv_r6 = inv_r4 * inv_r2;
    let inv_r7 = inv_r6 / r;

    let two_a_exp = 2.0 * w.alpha * exp_ar2 * INV_SQRT_PI;

    let fc1 = scale * (-erfc_ar * inv_r3 - two_a_exp * inv_r2 + w.slope / r);
    let fc2 = scale
        * (3.0 * erfc_ar * inv_r5
            + 6.0 * w.alpha * exp_ar2 * INV_SQRT_PI * inv_r4
            + 4.0 * w.alpha * w.alpha_sq * exp_ar2 * INV_SQRT_PI * inv_r2
            - w.slope * inv_r3);
    let fc3 = scale
        * (-15.0 * erfc_ar * inv_r7
            - 30.0 * w.alpha * exp_ar2 * INV_SQRT_PI * inv_r6
            - 20.0 * w.alpha * w.alpha_sq * exp_ar2 * INV_SQRT_PI * inv_r4
            - 8.0 * w.alpha * w.alpha_sq * w.alpha_sq * exp_ar2 * INV_SQRT_PI * inv_r2
            + 3.0 * w.slope * inv_r5);

    *f1 += fc1;
    *f2 += fc2;
    *f3 += fc3;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wolf_force_zero_at_cutoff() {
        let w = WolfConstants::new(0.25, 12.0);
        let r = w.cutoff;
        let r2 = r * r;
        let mut e = 0.0;
        let (mut f1, mut f2, mut f3) = (0.0, 0.0, 0.0);
        accumulate_wolf_dsf(&w, 0.5, r, r2, &mut e, &mut f1, &mut f2, &mut f3);
        assert!(e.abs() < 1e-10);
        assert!(f1.abs() < 1e-8);
    }
}
