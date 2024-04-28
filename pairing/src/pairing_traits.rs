use curve_traits::{AffinePoint, CurveParams};

//Trait defined for the Pairing
pub trait Pairing: CurveParams {

    // Target field for the pairing
    type TargetField;

    // Base curve i.e G1 curve
    type G1Curve: CurveParams;

    // Curve on extenion field i.e G2 curve
    type G2Curve: CurveParams;

    //  Function to compute the tate pairing, e: G_1 * G_2 -> G_t
    fn tate_pairing(m: AffinePoint<Self::G1Curve>, n: AffinePoint<Self::G2Curve>)
        -> Self::TargetField;

    // Function to compute the optimal ate pairing, e: G_2 * G_1 -> G_t
    fn optimal_ate_pairing(
        m: AffinePoint<Self::G2Curve>,
        n: AffinePoint<Self::G1Curve>,
    ) -> Self::TargetField;
}
