use bls381::scalar::Scalar;
use bls_curve::{bls::BlsCurve, fp2bls::G2BlsCurve};
use curve_traits::{AffinePoint, ProjectivePoint};
use pairing::pairing_traits::Pairing;

///Verify kzg proof
pub fn kzg_verify(
    commitment: AffinePoint<BlsCurve>,
    z: Scalar,
    evaluation: Scalar,
    proof: AffinePoint<BlsCurve>,
    verifier_key: &AffinePoint<G2BlsCurve>,
) -> bool {
    let g1 = ProjectivePoint::<BlsCurve>::GENERATOR;
    let g2 = ProjectivePoint::<G2BlsCurve>::GENERATOR;
    let commitment = commitment.to_projective().add(&g1.mul(-evaluation));

    // computing e(C-[b]_1, g_2)
    let lhs = BlsCurve::tate_pairing(commitment.to_affine(), g2.to_affine());

    let h_g2 = verifier_key.to_projective().add(&g2.mul(-z));
    // computing  e( [s-z]_1, pi)
    let rhs = BlsCurve::tate_pairing(proof, h_g2.to_affine());
    lhs == rhs
}
