use super::{bls_pairing::bls_compute_pairing, optimal::bls_optimal_ate_pairing};
use crate::pairing_traits::Pairing;
use bls_curve::gt::Gt;
use bls_curve::{bls::BlsCurve, fp2bls::G2BlsCurve};
use curve_traits::AffinePoint;
// Implementation of Pairing trait for BlsCurve
impl Pairing for BlsCurve {
    // Target field for the pairing
    type TargetField = Gt;

    //Base curve i.e G1 curve
    type G1Curve = BlsCurve;

    // Curve on extenion field i.e G2 curve
    type G2Curve = G2BlsCurve;
    //  Function to compute the tate pairing, e: G_1 * G_2 -> G_t
    fn tate_pairing(m: AffinePoint<BlsCurve>, n: AffinePoint<G2BlsCurve>) -> Self::TargetField {
        bls_curve::gt::Gt(bls_compute_pairing(m, n))
    }
    // Function to compute the optimal ate pairing, e: G_2 * G_1 -> G_t
    fn optimal_ate_pairing(
        m: AffinePoint<Self::G2Curve>,
        n: AffinePoint<Self::G1Curve>,
    ) -> Self::TargetField {
        bls_curve::gt::Gt(bls_optimal_ate_pairing(m, n))
    }
}
