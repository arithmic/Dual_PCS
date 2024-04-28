use super::{bn_pairing::bn_compute_pairing, optimal::bn_optimal_ate_pairing};
use crate::pairing_traits::Pairing;
use curve_traits::AffinePoint;
use BN_curve::{bncurve::BNCurve, fp2bn::G2BNCurve};
use BN_curve::gt::Gt;
// Implementation of Pairing trait for BNCurve
impl Pairing for BNCurve {
    // Target field for the pairing
    type TargetField = Gt;

    //Base curve i.e G1 curve
    type G1Curve = BNCurve;

    // Curve on extenion field i.e G2 curve
    type G2Curve = G2BNCurve;

    // Function to compute the tate pairing, e: G_1 * G_2 -> G_t
    fn tate_pairing(m: AffinePoint<BNCurve>, n: AffinePoint<G2BNCurve>) -> Self::TargetField {
        BN_curve::gt::Gt(bn_compute_pairing(m, n))
    }

    // Function to compute the optimal ate pairing, e: G_2 * G_1 -> G_t
    fn optimal_ate_pairing(
        m: AffinePoint<Self::G2Curve>,
        n: AffinePoint<Self::G1Curve>,
    ) -> Self::TargetField {
        BN_curve::gt::Gt(bn_optimal_ate_pairing(m, n))
    }
}
