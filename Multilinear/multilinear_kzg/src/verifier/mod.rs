use curve_traits::{self, projectivepoint::Scalar, AffinePoint, ProjectivePoint};
use pairing::pairing_traits::Pairing;
use std::{fmt::Debug, ops::Add};

use crate::common::{MleCommit, MleEvalProof, VerificationKey};

pub fn verify<C: Pairing>(
    commit: &MleCommit<C>,
    proof: &MleEvalProof<C>,
    point: &Vec<Scalar<C::G1Curve>>,
    verification_key: &VerificationKey<C>,
) -> Scalar<C::G1Curve>
where
    C::TargetField: for<'a> Add<&'a C::TargetField, Output = C::TargetField> + PartialEq + Debug,
{
    let g_v: ProjectivePoint<<C as Pairing>::G1Curve> =
        commit.commitment - ProjectivePoint::<C::G1Curve>::GENERATOR.mul(proof.evaluation);
    let witness_msm =
        -ProjectivePoint::multi_exponentiation(proof.witnesses.clone(), point.to_vec());
    let lhs_check: C::TargetField =
        C::tate_pairing(g_v.into(), AffinePoint::<C::G2Curve>::GENERATOR);
    let rhs_check_1 = proof
        .witnesses
        .iter()
        .enumerate()
        .map(|(i, y_i)| {
            C::tate_pairing(
                y_i.to_affine(),
                verification_key.key[i + verification_key.key.len() - point.len()].to_affine(),
            )
        })
        .reduce(|acc, g| acc + &g)
        .unwrap();

    let rhs_check_2 = C::tate_pairing(witness_msm.to_affine(), AffinePoint::GENERATOR);

    let rhs_check = rhs_check_1 + &rhs_check_2;
    assert_eq!(lhs_check, rhs_check, "Check failed");
    proof.evaluation
}

pub fn batch_verify<C: Pairing>(
    commits: &Vec<MleCommit<C>>,
    evals: &Vec<Scalar<C::G1Curve>>,
    scalars: &Vec<Scalar<C::G1Curve>>,
    proof: &MleEvalProof<C>,
    point: &Vec<Scalar<C::G1Curve>>,
    verification_key: &VerificationKey<C>,
) -> Vec<Scalar<C::G1Curve>>
where
    C::TargetField: for<'a> Add<&'a C::TargetField, Output = C::TargetField> + PartialEq + Debug,
{
    let combined_eval = evals
        .iter()
        .zip(scalars.iter())
        .map(|(eval, scalar)| *eval * *scalar)
        .reduce(|acc, g| acc + g)
        .unwrap();

    assert_eq!(combined_eval, proof.evaluation);

    let combined_commit = MleCommit {
        commitment: commits
            .iter()
            .zip(scalars.iter())
            .map(|(commit, scalar)| commit.commitment.mul(*scalar))
            .reduce(|acc, g| acc + g)
            .unwrap(),
    };
    verify(&combined_commit, proof, point, verification_key);

    evals.to_vec()
}

pub fn batch_commits<C: Pairing>(
    commits: &Vec<MleCommit<C>>,
    evals: &Vec<Scalar<C::G1Curve>>,
    scalars: &Vec<Scalar<C::G1Curve>>,
    proof: &MleEvalProof<C>,
) -> (MleCommit<C>, Scalar<C::G1Curve>)
where
    C::TargetField: for<'a> Add<&'a C::TargetField, Output = C::TargetField> + PartialEq + Debug,
{
    let combined_eval = evals
        .iter()
        .zip(scalars.iter())
        .map(|(eval, scalar)| *eval * *scalar)
        .reduce(|acc, g| acc + g)
        .unwrap();

    assert_eq!(combined_eval, proof.evaluation);

    let combined_commit = MleCommit {
        commitment: commits
            .iter()
            .zip(scalars.iter())
            .map(|(commit, scalar)| commit.commitment.mul(*scalar))
            .reduce(|acc, g| acc + g)
            .unwrap(),
    };

    (combined_commit, combined_eval)
}

pub fn batch_verify_var_openings<C: Pairing>(
    commits: &Vec<MleCommit<C>>,
    evals: &Vec<Scalar<C::G1Curve>>,
    proofs: &Vec<MleEvalProof<C>>,
    points: &Vec<Vec<Scalar<C::G1Curve>>>,
    scalars: &Vec<Scalar<C::G1Curve>>,
    verification_key: &VerificationKey<C>,
) where
    C::TargetField: for<'a> Add<&'a C::TargetField, Output = C::TargetField> + PartialEq + Debug,
{
    //Computing lhs pairing.

    let mut batched_evals: ProjectivePoint<C::G1Curve> = ProjectivePoint::IDENTITY;
    let n_witnesses = points[0].len();
    for i in 0..commits.len() {
        batched_evals +=
            (commits[i].commitment - ProjectivePoint::GENERATOR.mul(evals[i])).mul(scalars[i]);
    }
    let lhs = C::tate_pairing(
        batched_evals.to_affine(),
        AffinePoint::<C::G2Curve>::GENERATOR,
    );

    let mut scaled_proofs = Vec::new();

    for i in 0..proofs.len() {
        let mut scaled_proof = Vec::new();

        for j in 0..proofs[i].witnesses.len() {
            scaled_proof.push(proofs[i].witnesses[j].mul(scalars[i]))
        }
        scaled_proofs.push(scaled_proof);
    }

    let mut pairing_checks = Vec::new();
    for i in 0..n_witnesses {
        let mut combined_witness = ProjectivePoint::IDENTITY;

        for j in 0..proofs.len() {
            combined_witness += scaled_proofs[j][i];
        }

        pairing_checks.push(C::tate_pairing(
            combined_witness.to_affine(),
            verification_key.key[i + verification_key.key.len() - n_witnesses].to_affine(),
        ));
    }

    let mut witness_msms: ProjectivePoint<C::G1Curve> = ProjectivePoint::IDENTITY;
    for i in 0..proofs.len() {
        witness_msms = witness_msms
            - ProjectivePoint::multi_exponentiation(scaled_proofs[i].clone(), points[i].clone());
    }

    pairing_checks.push(C::tate_pairing(
        witness_msms.to_affine(),
        AffinePoint::<C::G2Curve>::GENERATOR,
    ));

    let rhs = pairing_checks
        .into_iter()
        .reduce(|acc, g| acc + &g)
        .unwrap();

    assert_eq!(lhs, rhs, "Verification Failed");
}
