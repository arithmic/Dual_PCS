use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use curve_traits::{AffinePoint, ProjectivePoint};
use fft::serial_fft::{get_inv_twiddles, interpolate_poly};

//  Commitment of the polynomial phi(x)= f_0 + x f_1 + x^2 f_2 + x^3 f_3 + x^4 f_4 +....x^d f_d
pub fn kzg_commit(
    mut polynomial: Vec<Scalar>,
    g1_powers: &Vec<ProjectivePoint<BlsCurve>>,
) -> (AffinePoint<BlsCurve>, Vec<Scalar>) {
    let inv_twiddles = get_inv_twiddles(polynomial.len() as u32);
    interpolate_poly(&mut polynomial, &inv_twiddles);
    // Computaion of the g^{f_0}+ (g^s)^{f_1} + (g^(s^2))^{f_2} +...+ (g^(s^d))^{f_d}
    (
        ProjectivePoint::pippenger_msm(&g1_powers, &polynomial).to_affine(),
        polynomial,
    )
}
