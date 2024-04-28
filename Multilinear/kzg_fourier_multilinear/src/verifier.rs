use bls381::scalar::Scalar;
use bls_curve::{
    bls::{BlsCurve, ProjectivePoint},
    fp2bls::G2BlsCurve,
};
use channel::Channel;
use common::{CommitmentsOfEvenOddPolys, EvaluationsAtRandomPoint, MultilinearKZG2Proof};
use curve_traits::AffinePoint;
use pairing::bls::bls_pairing::bls_compute_pairing;
use traits::traits::Field;

use crate::multilinear_common::reseed_with_even_odd_polys;
pub fn bench_verifier(
    proof: MultilinearKZG2Proof,
    random_points: Vec<Scalar>,
    verifier_key: AffinePoint<G2BlsCurve>,
    phi_k_commitments: Vec<AffinePoint<BlsCurve>>,
    commit_f: &AffinePoint<BlsCurve>,
) {
    let mut channel = Channel::initialize_with_affine_point([*commit_f].as_ref());

    verifier(
        proof,
        random_points,
        verifier_key,
        phi_k_commitments,
        commit_f,
        &mut channel,
    )
    .unwrap();
}

pub fn verifier(
    proof: MultilinearKZG2Proof,
    random_points: Vec<Scalar>,
    verifier_key: AffinePoint<G2BlsCurve>,
    phi_k_commitments: Vec<AffinePoint<BlsCurve>>,
    commit_f: &AffinePoint<BlsCurve>,
    channel: &mut Channel,
) -> Option<Scalar> {
    let num_var = random_points.len();

    let commitments_of_even_ood_polys = proof.get_commitments_of_even_odd_polys();

    //Reseed channel with commitments
    reseed_with_even_odd_polys(channel, commitments_of_even_ood_polys.clone());

    //Draw a random number using Fiat-shamir
    let z = channel.get_random_point();

    let evaluation = proof.get_evaluation();

    let evaluations_at_z = proof.get_evaluations_at_z();
    let phi_k_evaluation_at_z = proof.get_phi_k_at_z();
    let q_not = proof.get_q_not();
    first_verification_check(
        evaluations_at_z.clone(),
        evaluation,
        z,
        random_points,
        num_var,
        phi_k_evaluation_at_z,
        q_not,
    );

    //Draw random points
    let r_points = channel.get_random_points((5 * num_var - 4) as usize);

    let shi_commit = proof.get_shi_commit();
    channel.reseed_with_affine_point(&[shi_commit].to_vec());

    //Draw a random point using Fiat-shamir
    let s = channel.get_random_point();

    let evaluations_at_s = proof.get_evaluations_at_s();
    let phi_k_evaluation_at_s = proof.get_phi_k_at_s();

    let shi_at_s = evaluation_of_shi(
        evaluations_at_s,
        evaluations_at_z,
        s,
        z,
        &r_points,
        num_var,
        phi_k_evaluation_at_z,
        phi_k_evaluation_at_s,
    );

    let poly_combiners = channel.get_random_points(num_var * 3 - 2);

    let g_cap_evaluation = evaluation_of_g_cap(
        shi_at_s,
        evaluations_at_s,
        &poly_combiners,
        phi_k_evaluation_at_s,
    );

    let g_cap_commitment = compute_g_cap_commitment(
        shi_commit,
        commit_f,
        commitments_of_even_ood_polys,
        phi_k_commitments,
        &poly_combiners,
    );
    let g1_g_cap_evaluation =
        <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR
            .mul(g_cap_evaluation);

    //Generator of G2
    let g2 = <G2BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;
    let g2_s = g2.mul(-s);

    let commitment_r = proof.get_commitment_to_r();
    let lhs_pairing = bls_compute_pairing(
        (g_cap_commitment - g1_g_cap_evaluation).to_affine(),
        g2.to_affine(),
    );

    let rhs_pairing = bls_compute_pairing(
        commitment_r,
        (verifier_key.to_projective() + g2_s).to_affine(),
    );
    if lhs_pairing == rhs_pairing {
        Some(evaluation)
    } else {
        None
    }
}

pub(crate) fn first_verification_check(
    evaluations_at_z: EvaluationsAtRandomPoint,
    evaluation: Scalar,
    z: Scalar,
    mut random_points: Vec<Scalar>,
    num_var: usize,
    phi_k_evaluations: &Vec<Scalar>,
    q_zero: &Vec<Scalar>,
) {
    random_points.reverse();
    let even_polys_evaluations = evaluations_at_z.even_polys_evaluations;
    let odd_polys_evaluations = evaluations_at_z.odd_polys_evaluations;
    assert_eq!(
        even_polys_evaluations.len(),
        odd_polys_evaluations.len(),
        "length should match"
    );

    let (evaluations_of_g_at_z, evaluations_of_g_at_minus_z): (Vec<Scalar>, Vec<Scalar>) = (0
        ..even_polys_evaluations.len())
        .map(|idx| {
            let point = z.power_by([(1 << (num_var - (idx + 1) - 1)) as u64, 0, 0, 0]);
            let point_odd_evaluation = point * odd_polys_evaluations[idx];
            (
                even_polys_evaluations[idx] + point_odd_evaluation,
                even_polys_evaluations[idx] - point_odd_evaluation,
            )
        })
        .unzip();

    let mut evaluations_of_q = Vec::new();
    let mut evaluations_of_x_q = Vec::new();

    let inverse_of_minus_two = -Scalar::from(2_u8).invert().unwrap();
    let evaluations_of_x_qo = q_zero[0]
        * (z.power_by([(1 << (num_var - 1)) as u64, 0, 0, 0]) - Scalar::ONE)
        * inverse_of_minus_two;

    evaluations_of_q.push(q_zero[0]);
    evaluations_of_x_q.push(evaluations_of_x_qo);

    (0..evaluations_of_g_at_z.len()).for_each(|idx| {
        let evaluation = evaluations_of_g_at_z[idx] * phi_k_evaluations[idx];
        evaluations_of_x_q.push(evaluation);
        let z_power_minus_one = (z.power_by([1 << (num_var) as u64, 0, 0, 0])) - Scalar::ONE;
        let phi_k_positive = z_power_minus_one * (phi_k_evaluations[idx].invert().unwrap());
        let evaluation = evaluation + (evaluations_of_g_at_minus_z[idx] * phi_k_positive);
        evaluations_of_q.push(evaluation);
    });

    let lhs = evaluations_at_z.evaluation_of_f - evaluation;
    let rhs = (0..evaluations_of_q.len()).fold(Scalar::ZERO, |acc, idx| {
        acc + (evaluations_of_x_q[idx] - random_points[idx] * evaluations_of_q[idx])
    });
    assert_eq!(lhs, rhs, "Initial check failed.");
}

fn evaluation_of_shi(
    evaluations_at_s: &EvaluationsAtRandomPoint,
    evaluations_at_z: &EvaluationsAtRandomPoint,
    s: Scalar,
    z: Scalar,
    random_points: &[Scalar],
    num_var: usize,
    phi_k_at_z: &Vec<Scalar>,
    phi_k_at_s: &Vec<Scalar>,
) -> Scalar {
    let evaluation_of_f_at_s = evaluations_at_s.evaluation_of_f;
    let evaluation_of_f_at_z = evaluations_at_z.evaluation_of_f;
    let even_poly_evaluation_at_s = &evaluations_at_s.even_polys_evaluations;
    let even_poly_evaluation_at_z = &evaluations_at_z.even_polys_evaluations;
    let odd_poly_evaluation_at_s = &evaluations_at_s.odd_polys_evaluations;
    let odd_poly_evaluation_at_z = &evaluations_at_z.odd_polys_evaluations;
    let evaluation = (0..phi_k_at_z.len()).fold(Scalar::ZERO, |acc, idx| {
        let z_power_d_k = z.power_by([(1 << (num_var - (idx + 1))) as u64, 0, 0, 0]);
        let s_power_d_k = s.power_by([((1 << num_var) - (1 << idx)) as u64, 0, 0, 0]);
        let alpha = z_power_d_k.power_by([((1 << num_var) - (1 << idx)) as u64, 0, 0, 0]);
        acc + random_points[5 * idx + 1]
            * (phi_k_at_s[idx] - phi_k_at_z[idx])
            * (s.power_by([(1 << (num_var - (idx + 1) - 1)) as u64, 0, 0, 0])
                - z.power_by([(1 << (num_var - (idx + 1) - 1)) as u64, 0, 0, 0]))
            .invert()
            .unwrap()
            + ((random_points[5 * idx + 2]
                * (even_poly_evaluation_at_s[idx] - even_poly_evaluation_at_z[idx]))
                + (random_points[5 * idx + 3]
                    * (odd_poly_evaluation_at_s[idx] - odd_poly_evaluation_at_z[idx]))
                + (random_points[5 * idx + 4]
                    * (s_power_d_k * even_poly_evaluation_at_s[idx]
                        - even_poly_evaluation_at_z[idx] * alpha))
                + (random_points[5 * idx + 5]
                    * (s_power_d_k * odd_poly_evaluation_at_s[idx]
                        - odd_poly_evaluation_at_z[idx] * alpha)))
                * (s - z_power_d_k).invert().unwrap()
    }) + (random_points[0]
        * (evaluation_of_f_at_s - evaluation_of_f_at_z)
        * (s - z).invert().unwrap());

    evaluation
}
fn evaluation_of_g_cap(
    shi_at_s: Scalar,
    evaluations_at_s: &EvaluationsAtRandomPoint,
    random_points: &[Scalar],
    phi_k_evaluations_at_s: &Vec<Scalar>,
) -> Scalar {
    assert_eq!(
        evaluations_at_s.even_polys_evaluations.len(),
        evaluations_at_s.odd_polys_evaluations.len(),
        "no of even polys evaluations and odd polys evaluations should be same"
    );

    assert_eq!(
        evaluations_at_s.even_polys_evaluations.len(),
        phi_k_evaluations_at_s.len(),
        "length mismatched"
    );
    let f_s = evaluations_at_s.evaluation_of_f;
    let even_polys_evaluations = &evaluations_at_s.even_polys_evaluations;
    let odd_polys_evaluations = &evaluations_at_s.odd_polys_evaluations;

    (0..evaluations_at_s.even_polys_evaluations.len()).fold(Scalar::ZERO, |acc, iter| {
        acc + (random_points[3 * iter + 1] * even_polys_evaluations[iter])
            + (random_points[3 * iter + 2] * odd_polys_evaluations[iter])
            + (random_points[3 * iter + 3] * phi_k_evaluations_at_s[iter])
    }) + (random_points[0] * f_s)
        + shi_at_s
}
fn compute_g_cap_commitment(
    shi_commit: AffinePoint<BlsCurve>,
    commit_f: &AffinePoint<BlsCurve>,
    commitments_of_even_ood_polys: &CommitmentsOfEvenOddPolys,
    phi_k_commitments: Vec<AffinePoint<BlsCurve>>,
    random_points: &[Scalar],
) -> ProjectivePoint {
    assert_eq!(
        commitments_of_even_ood_polys.even_polys_commits.len(),
        commitments_of_even_ood_polys.odd_polys_commits.len(),
        "no of even polys commits and odd polys commits should be same"
    );

    assert_eq!(
        commitments_of_even_ood_polys.even_polys_commits.len(),
        phi_k_commitments.len(),
        "length mismatched 2"
    );
    let even_poly_commits = &commitments_of_even_ood_polys.even_polys_commits;
    let odd_polys_commits = &commitments_of_even_ood_polys.odd_polys_commits;

    let commitment_of_g = (0..commitments_of_even_ood_polys.even_polys_commits.len()).fold(
        ProjectivePoint::IDENTITY,
        |acc, iter| {
            acc + (even_poly_commits[iter]
                .to_projective()
                .mul(random_points[3 * iter + 1]))
                + (odd_polys_commits[iter]
                    .to_projective()
                    .mul(random_points[3 * iter + 2]))
                + (phi_k_commitments[iter]
                    .to_projective()
                    .mul(random_points[3 * iter + 3]))
        },
    ) + commit_f.to_projective().mul(random_points[0])
        + shi_commit.to_projective();

    commitment_of_g
}

pub fn batch_verify(
    proof: MultilinearKZG2Proof,
    commits: Vec<AffinePoint<BlsCurve>>,
    evals: Vec<Scalar>,
    combiners: Vec<Scalar>,
    random_points: Vec<Scalar>,
    verifier_key: AffinePoint<G2BlsCurve>,
    phi_k_commitments: Vec<AffinePoint<BlsCurve>>,
    channel: &mut Channel,
) -> Option<Vec<Scalar>> {
    let combined_eval =
        (0..evals.len()).fold(Scalar::ZERO, |acc, idx| acc + (combiners[idx] * evals[idx]));

    let combined_commit = (0..commits.len())
        .fold(ProjectivePoint::IDENTITY, |acc, idx| {
            acc + (commits[idx].to_projective().mul(combiners[idx]))
        })
        .to_affine();
    let evaluation = verifier(
        proof,
        random_points,
        verifier_key,
        phi_k_commitments,
        &combined_commit,
        channel,
    );
    match evaluation {
        Some(eval) => {
            assert_eq!(eval, combined_eval, "assertion failed");
            Some(evals)
        }
        None => None,
    }
}
