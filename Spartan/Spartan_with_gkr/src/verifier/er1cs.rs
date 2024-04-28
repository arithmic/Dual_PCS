#![allow(unused)]
use crate::{spartan_common, verifier::batch_opening::batch_opening};
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use channel::Channel;
use multilinear_kzg::{
    common::VerificationKey,
    verifier::{batch_verify_var_openings, verify},
};
use polynomial::MultPolynomial;
use spartan_common::{
    eR1CSCommitments, eR1CStranscript, InitialSumCheckTranscript, ParSumCheckTranscript,
};
use traits::traits::Field;
#[allow(non_snake_case)]
pub fn verify_sat(
    transcript: eR1CStranscript,
    commitments: eR1CSCommitments,
    u: Scalar,
    PI: MultPolynomial,
    pi_indices: Vec<usize>,
    ver_key: &VerificationKey<BlsCurve>,
    channel: &mut Channel,
) {
    let rx = transcript.first_sum_check_transcript.random_points.clone();
    let ry = transcript.par_sum_check_transcript.random_points.clone();

    let Az_claimed_val = transcript.first_sum_check_transcript.Az_claimed_val;
    let Bz_claimed_val = transcript.first_sum_check_transcript.Bz_claimed_val;
    let Cz_claimed_val = transcript.first_sum_check_transcript.Cz_claimed_val;
    let E_final_eval = transcript.E_eval_proof.evaluation;

    let first_sum_check_final_eval =
        (Az_claimed_val * Bz_claimed_val) - (u * Cz_claimed_val + E_final_eval);

    initial_sum_check_verification(
        &transcript.first_sum_check_transcript,
        first_sum_check_final_eval,
        channel,
    );

    let random_coeffs = channel.get_random_points(3);

    let eval = random_coeffs[0] * Az_claimed_val
        + random_coeffs[1] * Bz_claimed_val
        + random_coeffs[2] * Cz_claimed_val;

    let PI_eval = evaluate_PI(pi_indices, PI, &ry, ry.len());
    let eval_point = [rx, ry.clone()].concat();
    let W_eval = transcript.W_eval_proof.evaluation;

    let Z_final_eval = PI_eval + W_eval;

    let A_claimed_val = transcript.par_sum_check_transcript.A_claimed_val;
    let B_claimed_val = transcript.par_sum_check_transcript.B_claimed_val;
    let C_claimed_val = transcript.par_sum_check_transcript.C_claimed_val;

    let par_sum_check_final_eval = Z_final_eval
        * (random_coeffs[0] * A_claimed_val
            + random_coeffs[1] * B_claimed_val
            + random_coeffs[2] * C_claimed_val);

    par_sum_check_verification(
        transcript.par_sum_check_transcript,
        eval,
        par_sum_check_final_eval,
        channel,
    );

    // Verifying the claimed values.
    let (
        sum_check_batch_commit,
        _sum_check_batch_eval,
        sum_check_eval_proof,
        sum_check_random_points,
        gkr_commit1,
        gkr_commit2,
        gkr_eval1,
        _gkr_eval2,
        gkr_batch_eval_proof1,
        gkr_batch_eval_proof2,
        gkr_final_layer_point1,
        gkr_final_layer_point2,
    ) = batch_opening(
        vec![A_claimed_val, B_claimed_val, C_claimed_val],
        transcript.BatchProof,
        vec![commitments.A, commitments.B, commitments.C],
        &eval_point,
        channel,
    );

    let mut commits = Vec::new();
    let mut points = Vec::new();
    let mut proofs = Vec::new();
    let mut evals = Vec::new();

    commits.push(commitments.E);
    points.push(transcript.first_sum_check_transcript.random_points);
    proofs.push(transcript.E_eval_proof);
    evals.push(E_final_eval);

    commits.push(commitments.W);
    points.push(ry);
    proofs.push(transcript.W_eval_proof);
    evals.push(W_eval);
    commits.push(gkr_commit1);
    points.push(gkr_final_layer_point1);
    proofs.push(gkr_batch_eval_proof1);
    evals.push(gkr_eval1);

    let scalars = channel.get_random_points(4);
    batch_verify_var_openings(&commits, &evals, &proofs, &points, &scalars, ver_key);

    verify(
        &sum_check_batch_commit,
        &sum_check_eval_proof,
        &sum_check_random_points,
        ver_key,
    );

    verify(
        &gkr_commit2,
        &gkr_batch_eval_proof2,
        &gkr_final_layer_point2,
        ver_key,
    );
}

pub fn initial_sum_check_verification(
    transcript: &InitialSumCheckTranscript,
    final_evaluation: Scalar,
    channel: &mut Channel,
) {
    let polynomials = &transcript.polynomials;
    let rounds = polynomials.len();
    let tau = channel.get_random_points(rounds);

    let mut current_sum = Scalar::ZERO;

    let mut r = vec![Scalar::ZERO; rounds];
    for i in 0..rounds {
        let poly = polynomials[i].get_coefficients();
        assert_eq!(
            current_sum,
            eval(poly, Scalar::ZERO) + eval(poly, Scalar::ONE),
            "f(0) + f(1) did not match binding at round {:?} in the eR1CS initial sum check",
            i
        );

        channel.reseed_with_scalars(poly);
        let r_i = channel.get_random_point();
        r[i] = r_i;
        current_sum = eval(poly, r_i)
    }
    let eq = evaluate_eq(tau, r);
    assert_eq!(
        current_sum,
        eq * final_evaluation,
        "Final assertion in eR1CS initial sum check failed"
    )
}

pub fn par_sum_check_verification(
    transcript: ParSumCheckTranscript,
    initial_evaluation: Scalar,
    final_evaluation: Scalar,
    channel: &mut Channel,
) {
    let polynomials = &transcript.polynomials;
    let rounds = polynomials.len();

    let mut current_sum = initial_evaluation;

    let mut r = vec![Scalar::ZERO; rounds];
    for i in 0..rounds {
        let poly = polynomials[i].get_coefficients();
        assert_eq!(
            current_sum,
            eval(poly, Scalar::ZERO) + eval(poly, Scalar::ONE),
            "f(0) + f(1) did not match binding at round {:?} in the second eR1CS sum check",
            i
        );

        channel.reseed_with_scalars(poly);
        let r_i = channel.get_random_point();
        r[i] = r_i;
        current_sum = eval(poly, r_i)
    }

    assert_eq!(
        current_sum, final_evaluation,
        "Final assertion in eR1CS second sum check failed"
    )
}
#[allow(non_snake_case)]
pub fn evaluate_PI(
    pi_indices: Vec<usize>,
    PI: MultPolynomial,
    point: &Vec<Scalar>,
    bits: usize,
) -> Scalar {
    let mut eval = Scalar::ZERO;

    for j in 0..PI.len() {
        let mut basis_eval = PI.get_coeff(j);
        for i in (0..bits).rev() {
            if (pi_indices[j] >> (i)) & 1 == 1 {
                basis_eval *= point[bits - i - 1]
            } else {
                basis_eval *= Scalar::ONE - point[bits - i - 1]
            }
        }
        eval += basis_eval;
    }
    eval
}

pub fn evaluate_eq(basis_point: Vec<Scalar>, evaluation_point: Vec<Scalar>) -> Scalar {
    let mut res = Scalar::ONE;
    for (x, y) in basis_point.iter().zip(evaluation_point.iter()) {
        res *= Scalar::ONE - *x - *y + (*x * *y).double()
    }
    res
}

//...........
// CODE  for evaluating polynomial at points
//.............
pub fn eval(p: &[Scalar], x: Scalar) -> Scalar {
    // Horner evaluation
    p.iter()
        .rev()
        .fold(Scalar::ZERO, |acc, &coeff| acc * x + coeff)
}
