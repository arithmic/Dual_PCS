use bls381::scalar::Scalar;
use bls_curve::{
    bls::{AffinePoint, BlsCurve, ProjectivePoint},
    fp2bls::G2BlsCurve,
};
use common::KZGFFTEvaluationProof;
use pairing::pairing_traits::Pairing;

pub fn kzg2_verifier(
    verifier_key: curve_traits::AffinePoint<G2BlsCurve>,
    evaluation_proof: KZGFFTEvaluationProof,
    commitment_to_f: AffinePoint,
    z: Scalar,
) -> Option<Scalar> {
    //Generator of G1
    let g1 = <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;
    //Generator of G2
    let g2 = <G2BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;
    let h_g2 = verifier_key.to_projective() + g2.mul(-z);
    let pairing1 = BlsCurve::tate_pairing(
        (commitment_to_f.to_projective() + g1.mul(-evaluation_proof.evaluation_at_z)).to_affine(),
        g2.to_affine(),
    );
    let pairing2 = BlsCurve::tate_pairing(
        evaluation_proof.commitment_of_evaluations_of_q,
        h_g2.to_affine(),
    );
    if pairing1 == pairing2 {
        Some(evaluation_proof.evaluation_at_z)
    } else {
        None
    }
}

pub fn kzg2_batch_verify(
    evaluation_proof: KZGFFTEvaluationProof,
    verifier_key: curve_traits::AffinePoint<G2BlsCurve>,
    random_point: Scalar,
    commitments: Vec<AffinePoint>,
    evals: Vec<Scalar>,
    combiners: Vec<Scalar>,
) {
    let combined_eval =
        (0..evals.len()).fold(Scalar::ZERO, |acc, idx| acc + (combiners[idx] * evals[idx]));
    let combined_commit = (0..commitments.len())
        .fold(ProjectivePoint::IDENTITY, |acc, idx| {
            acc + (commitments[idx].to_projective().mul(combiners[idx]))
        })
        .to_affine();
    let evaluation = kzg2_verifier(
        verifier_key,
        evaluation_proof,
        combined_commit,
        random_point,
    )
    .unwrap();
    assert_eq!(evaluation, combined_eval, "Evaluations are incorrect");
}
