use bls381::scalar::Scalar;
use bls_curve::{
    bls::BlsCurve,
    fp2bls::{G2BlsCurve, G2ProjectivePoint},
};
use curve_traits::{AffinePoint, ProjectivePoint};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::Field;

// Function to compute the setup
// g is the generator of the group E(F_q) and h is the generator of the group E'(F_q^2).
// Computes g1_powers and g2_powers of setup.
// recieve1=(g,g^s,g^(s^2),...,g^(s^d)) and recieve2=(h,h^s,h^(s^2),...,h^(s^d)) where d=max_degree of polynomials to be committed and s is the secret key
pub fn kzg_setup(max_degree: usize) -> (Vec<ProjectivePoint<BlsCurve>>, AffinePoint<G2BlsCurve>) {
    let s = Scalar::random();
    let g1_powers = (0..max_degree)
        .into_par_iter()
        .map(|idx| {
            let s_power = s.power_by(&[idx as u64, 0, 0, 0]);
            // Computing g,g^s,g^(s^2),...,g^(s^d)
            ProjectivePoint::mul(ProjectivePoint::GENERATOR, s_power)
        })
        .collect();

    let g2_power_s = G2ProjectivePoint::mul(G2ProjectivePoint::GENERATOR, s);

    return (g1_powers, g2_power_s.to_affine());
}
