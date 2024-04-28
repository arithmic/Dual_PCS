use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use curve_traits::{AffinePoint, ProjectivePoint};

use crate::polynomial::Polynomial;

// Computes commitment to quotient polynomial q(X)=(phi(X)-y)/(X-z) where y is phi(z) i.e evaluation of the polynomial at point z.
// Commitment is done using the g2_powers of the setup
pub fn kzg_prover(
    mut polynomial: Polynomial,
    g1_powers: &Vec<ProjectivePoint<BlsCurve>>,
    z: Scalar,
) -> AffinePoint<BlsCurve> {
    // Computation of phi(z)
    let y = polynomial.evaluate(&z);
    polynomial.coeffs[0] -= y;
    // Computation of the polynomial q(X)=(phi(X)-y)/(X-z)
    let q_x = polynomial.divide_by_lin_factor(z);
    // Commitment of the polynomial q(X)=(phi(X)-y)/(X-z)
    ProjectivePoint::pippenger_msm(g1_powers, &q_x.coeffs).to_affine()
}
