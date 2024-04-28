use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use curve_traits::{AffinePoint, ProjectivePoint};

pub fn kzg2commit(
    evaluations: &Vec<Scalar>,
    prover_key: &Vec<curve_traits::ProjectivePoint<BlsCurve>>,
) -> AffinePoint<BlsCurve> {
    assert_eq!(evaluations.len(), prover_key.len(), "length should be same");
    ProjectivePoint::pippenger_msm(prover_key, evaluations).to_affine()
}
