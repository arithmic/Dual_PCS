use crate::{multilinear_common::reseed_with_even_odd_polys, verifier::first_verification_check};
use bls381::scalar::Scalar;
use bls_curve::{
    bls::{BlsCurve, ProjectivePoint},
    fp2bls::G2BlsCurve,
};
use channel::Channel;
use common::{BatchMultilinearKZG2Proof, CommitmentsOfEvenOddPolys, EvaluationsAtRandomPoint};
use curve_traits::AffinePoint;
use pairing::bls::bls_pairing::bls_compute_pairing;
use traits::traits::Field;

pub fn batch_verifier(
    proof: BatchMultilinearKZG2Proof,
    random_points: Vec<Vec<Scalar>>,
    verifier_key: AffinePoint<G2BlsCurve>,
    phi_k_commitments: Vec<AffinePoint<BlsCurve>>,
    commit_f: &Vec<AffinePoint<BlsCurve>>,
    channel: &mut Channel,
) -> Option<Vec<Scalar>> {
    let num_var = random_points[0].len();
    let commitments_of_even_ood_polys = proof.commitments_of_even_ood_polys;

    //Reseed channel with commitments
    (0..commitments_of_even_ood_polys.len()).for_each(|idx| {
        reseed_with_even_odd_polys(channel, commitments_of_even_ood_polys[idx].clone())
    });

    //Draw a random number using Fiat-shamir
    let z = channel.get_random_point();

    let evaluations = proof.evaluations;

    let evaluations_at_z = proof.evaluations_at_z;
    let q_not_polys = proof.q_not_polys;
    let phi_k_evaluation_at_z = proof.phi_k_evaluation_at_z;
    (0..evaluations.len()).for_each(|idx| {
        first_verification_check(
            evaluations_at_z[idx].clone(),
            evaluations[idx],
            z,
            random_points[idx].clone(),
            num_var,
            &phi_k_evaluation_at_z,
            &q_not_polys[idx],
        )
    });
    //Draw random points
    let no_of_points = (4 * num_var as usize - 3) * random_points.len() + (num_var as usize - 1);
    let r_points = channel.get_random_points(no_of_points as usize);

    let shi_commit = proof.shi_commit;
    channel.reseed_with_affine_point(&[shi_commit].to_vec());

    //Draw a random point using Fiat-shamir
    let s = channel.get_random_point();

    let evaluations_at_s = proof.evaluations_at_s;
    let phi_k_evaluation_at_s = proof.phi_k_evaluation_at_s;

    let shi_at_s = evaluation_of_shi(
        &evaluations_at_s,
        &evaluations_at_z,
        s,
        z,
        &r_points,
        num_var,
        &phi_k_evaluation_at_z,
        &phi_k_evaluation_at_s,
    );

    //Draw random points to combine all the polys using Fiat-Shamir
    let poly_combiners = channel.get_random_points(
        (num_var as usize * 2 - 1) * random_points.len() + (num_var as usize - 1),
    );

    let g_cap_evaluation = evaluation_of_g_cap(
        shi_at_s,
        &evaluations_at_s,
        &poly_combiners,
        phi_k_evaluation_at_s,
    );

    let g_cap_commitment = compute_g_cap_commitment(
        shi_commit,
        commit_f,
        &commitments_of_even_ood_polys,
        phi_k_commitments,
        &poly_combiners,
    );
    let g1_g_cap_evaluation =
        <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR
            .mul(g_cap_evaluation);

    //Generator of G2
    let g2 = <G2BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;
    let g2_s = g2.mul(-s);

    let commitment_r = proof.commitment_to_r;
    let lhs_pairing = bls_compute_pairing(
        (g_cap_commitment - g1_g_cap_evaluation).to_affine(),
        g2.to_affine(),
    );

    let rhs_pairing = bls_compute_pairing(
        commitment_r,
        (verifier_key.to_projective() + g2_s).to_affine(),
    );
    if lhs_pairing == rhs_pairing {
        Some(evaluations)
    } else {
        None
    }
}

fn evaluation_of_shi(
    evaluations_at_s: &Vec<EvaluationsAtRandomPoint>,
    evaluations_at_z: &Vec<EvaluationsAtRandomPoint>,
    s: Scalar,
    z: Scalar,
    random_points: &[Scalar],
    num_var: usize,
    phi_k_at_z: &Vec<Scalar>,
    phi_k_at_s: &Vec<Scalar>,
) -> Scalar {
    let num_of_polys = evaluations_at_s.len();
    let evaluation_of_f_at_s = &(0..evaluations_at_s.len())
        .map(|idx| evaluations_at_s[idx].evaluation_of_f)
        .collect::<Vec<Scalar>>();
    let evaluation_of_f_at_z = &(0..evaluations_at_s.len())
        .map(|idx| evaluations_at_z[idx].evaluation_of_f)
        .collect::<Vec<Scalar>>();
    let even_poly_evaluation_at_s = &(0..evaluations_at_s.len())
        .map(|idx| evaluations_at_s[idx].clone().even_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();
    let even_poly_evaluation_at_z = &(0..evaluations_at_s.len())
        .map(|idx| evaluations_at_z[idx].clone().even_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();
    let odd_poly_evaluation_at_s = &(0..evaluations_at_s.len())
        .map(|idx| evaluations_at_s[idx].clone().odd_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();
    let odd_poly_evaluation_at_z = &(0..evaluations_at_s.len())
        .map(|idx| evaluations_at_z[idx].clone().odd_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();

    let evaluation = (0..phi_k_at_z.len()).fold(Scalar::ZERO, |acc, outer_idx| {
        let z_power_d_k = z.power_by([(1 << (num_var - (outer_idx + 1))) as u64, 0, 0, 0]);
        let s_power_d_k = s.power_by([((1 << num_var) - (1 << outer_idx)) as u64, 0, 0, 0]);
        let alpha = z_power_d_k.power_by([((1 << num_var) - (1 << outer_idx)) as u64, 0, 0, 0]);
        acc + random_points[num_of_polys + outer_idx]
            * (phi_k_at_s[outer_idx] - phi_k_at_z[outer_idx])
            * ((s.power_by([(1 << (num_var - (outer_idx + 1) - 1)) as u64, 0, 0, 0])
                - z.power_by([(1 << (num_var - (outer_idx + 1) - 1)) as u64, 0, 0, 0]))
            .invert()
            .unwrap())
            + ((0..num_of_polys).fold(Scalar::ZERO, |acc1, inner_idx| {
                acc1 + (random_points[4 * outer_idx + phi_k_at_z.len() + num_of_polys + inner_idx]
                    * (even_poly_evaluation_at_s[inner_idx][outer_idx]
                        - even_poly_evaluation_at_z[inner_idx][outer_idx]))
                    + (random_points
                        [4 * outer_idx + phi_k_at_z.len() + num_of_polys + 1 + inner_idx]
                        * (odd_poly_evaluation_at_s[inner_idx][outer_idx]
                            - odd_poly_evaluation_at_z[inner_idx][outer_idx]))
                    + (random_points
                        [4 * outer_idx + phi_k_at_z.len() + num_of_polys + 2 + inner_idx]
                        * (s_power_d_k * even_poly_evaluation_at_s[inner_idx][outer_idx]
                            - even_poly_evaluation_at_z[inner_idx][outer_idx] * alpha))
                    + (random_points
                        [4 * outer_idx + phi_k_at_z.len() + num_of_polys + 3 + inner_idx]
                        * (s_power_d_k * odd_poly_evaluation_at_s[inner_idx][outer_idx]
                            - odd_poly_evaluation_at_z[inner_idx][outer_idx] * alpha))
            })) * (s - z_power_d_k).invert().unwrap()
    }) + (0..num_of_polys).fold(Scalar::ZERO, |acc, idx| {
        acc + random_points[idx] * (evaluation_of_f_at_s[idx] - evaluation_of_f_at_z[idx])
    }) * (s - z).invert().unwrap();

    evaluation
}
fn compute_g_cap_commitment(
    shi_commit: AffinePoint<BlsCurve>,
    commit_f: &Vec<AffinePoint<BlsCurve>>,
    commitments_of_even_ood_polys: &Vec<CommitmentsOfEvenOddPolys>,
    phi_k_commitments: Vec<AffinePoint<BlsCurve>>,
    random_points: &[Scalar],
) -> ProjectivePoint {
    let num_polys = commit_f.len();
    assert_eq!(
        commitments_of_even_ood_polys[0].even_polys_commits.len(),
        commitments_of_even_ood_polys[0].odd_polys_commits.len(),
        "no of even polys commits and odd polys commits should be same"
    );

    assert_eq!(
        commitments_of_even_ood_polys[0].even_polys_commits.len(),
        phi_k_commitments.len(),
        "length mismatched 2"
    );
    let even_poly_commits = &(0..commitments_of_even_ood_polys.len())
        .map(|idx| {
            commitments_of_even_ood_polys[idx]
                .clone()
                .even_polys_commits
        })
        .collect::<Vec<Vec<AffinePoint<BlsCurve>>>>();
    let odd_polys_commits = &(0..commitments_of_even_ood_polys.len())
        .map(|idx| commitments_of_even_ood_polys[idx].clone().odd_polys_commits)
        .collect::<Vec<Vec<AffinePoint<BlsCurve>>>>();
    let num_var_minus_1 = commitments_of_even_ood_polys[0].even_polys_commits.len();
    let commitment_of_g = (0..num_var_minus_1).fold(ProjectivePoint::IDENTITY, |acc, iter| {
        acc + (0..num_polys).fold(ProjectivePoint::IDENTITY, |acc1, inner_idx| {
            acc1 + (even_poly_commits[inner_idx][iter]
                .to_projective()
                .mul(random_points[2 * iter + num_polys + num_var_minus_1 + inner_idx]))
                + (odd_polys_commits[inner_idx][iter]
                    .to_projective()
                    .mul(random_points[2 * iter + num_polys + num_var_minus_1 + 1 + inner_idx]))
        }) + (phi_k_commitments[iter]
            .to_projective()
            .mul(random_points[num_polys + iter]))
    }) + (0..num_polys).fold(ProjectivePoint::IDENTITY, |acc, idx| {
        acc + commit_f[idx].to_projective().mul(random_points[idx])
    }) + shi_commit.to_projective();

    commitment_of_g
}
fn evaluation_of_g_cap(
    shi_at_s: Scalar,
    evaluations_at_s: &Vec<EvaluationsAtRandomPoint>,
    random_points: &[Scalar],
    phi_k_evaluation_at_s: Vec<Scalar>,
) -> Scalar {
    let num_polys = evaluations_at_s.len();
    assert_eq!(
        evaluations_at_s[0].even_polys_evaluations.len(),
        evaluations_at_s[0].odd_polys_evaluations.len(),
        "no of even polys evaluations and odd polys evaluations should be same"
    );

    assert_eq!(
        evaluations_at_s[0].even_polys_evaluations.len(),
        phi_k_evaluation_at_s.len(),
        "length mismatched"
    );
    let f_s = (0..evaluations_at_s.len())
        .map(|idx| evaluations_at_s[idx].evaluation_of_f)
        .collect::<Vec<Scalar>>();
    let even_polys_evaluations = &(0..evaluations_at_s.len())
        .map(|idx| evaluations_at_s[idx].clone().even_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();
    let odd_polys_evaluations = &(0..evaluations_at_s.len())
        .map(|idx| evaluations_at_s[idx].clone().odd_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();

    (0..phi_k_evaluation_at_s.len()).fold(Scalar::ZERO, |acc, iter| {
        acc + (0..num_polys).fold(Scalar::ZERO, |acc1, inner_idx| {
            acc1 + (random_points[2 * iter + num_polys + phi_k_evaluation_at_s.len() + inner_idx]
                * even_polys_evaluations[inner_idx][iter])
                + (random_points
                    [2 * iter + num_polys + phi_k_evaluation_at_s.len() + 1 + inner_idx]
                    * odd_polys_evaluations[inner_idx][iter])
        }) + (random_points[num_polys + iter] * phi_k_evaluation_at_s[iter])
    }) + (0..num_polys).fold(Scalar::ZERO, |acc, idx| {
        acc + (random_points[idx] * f_s[idx])
    }) + shi_at_s
}
