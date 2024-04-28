use bls381::scalar::Scalar;
use bls_curve::{bls::AffinePoint, fp2bls::G2AffinePoint};
use channel::Channel;
use kzg_fourier_multilinear::{
    batch_verifier::batch_verifier,
    common::{BatchMultilinearKZG2Proof, KZGFourierDegreeBoundSetup},
};
use polynomial::MultPolynomial;
use spartan_common::{
    eR1CSCommitments, eR1CStranscript, eval, InitialSumCheckTranscript, ParSumCheckTranscript,
};
use traits::traits::Field;

use crate::{spartan_common, verifier::batch_opening::batch_opening};
#[allow(non_snake_case)]
pub fn verify_sat(
    transcript: eR1CStranscript,
    u: Scalar,
    PI: MultPolynomial,
    pi_indices: Vec<usize>,
    setup: &KZGFourierDegreeBoundSetup,
    eR1CSCommitments: eR1CSCommitments,
    channel: &mut Channel,
) {
    let rx = transcript.first_sum_check_transcript.random_points.clone();
    let ry = transcript.par_sum_check_transcript.random_points.clone();

    let Az_claimed_val = transcript.first_sum_check_transcript.Az_claimed_val;
    let Bz_claimed_val = transcript.first_sum_check_transcript.Bz_claimed_val;
    let Cz_claimed_val = transcript.first_sum_check_transcript.Cz_claimed_val;

    let rounds = transcript.first_sum_check_transcript.polynomials.len();
    let tau = channel.get_random_points(rounds);
    let mut initial_sum_check_random_points = Vec::new();
    for idx in 0..rounds {
        channel.reseed_with_scalars(
            &transcript.first_sum_check_transcript.polynomials[idx].get_coefficients(),
        );
        initial_sum_check_random_points.push(channel.get_random_point())
    }

    let random_coeffs = channel.get_random_points(3);
    let eval = random_coeffs[0] * Az_claimed_val
        + random_coeffs[1] * Bz_claimed_val
        + random_coeffs[2] * Cz_claimed_val;

    let PI_eval = evaluate_PI(pi_indices, PI, &ry, ry.len());

    let rounds = transcript.par_sum_check_transcript.polynomials.len();
    let mut par_sum_check_random_points = Vec::new();
    for idx in 0..rounds {
        channel.reseed_with_scalars(
            transcript.par_sum_check_transcript.polynomials[idx].get_coefficients(),
        );
        par_sum_check_random_points.push(channel.get_random_point());
    }

    let setup_for_specific_degree = setup.get_setup(1 << rx.len());
    let commits = vec![
        eR1CSCommitments.E_Commit,
        eR1CSCommitments.W_Commit,
        transcript.rx_basis_commit,
        transcript.ry_basis_commit,
    ];
    let random_points = channel.get_random_points(initial_sum_check_random_points.len());
    let combiners = channel.get_random_points(2);
    let evaluations = multilinear_batch_open(
        transcript.e_w_eval_proof,
        [
            initial_sum_check_random_points.clone(),
            par_sum_check_random_points.clone(),
            random_points.clone(),
        ]
        .to_vec(),
        setup_for_specific_degree.setup.verifier_key,
        setup_for_specific_degree.phi_commitment,
        commits,
        channel,
        &combiners,
    );
    let eq_rx = evaluate_eq(&rx, &random_points);
    let eq_ry = evaluate_eq(&ry, &random_points);
    assert_eq!(evaluations[2], combiners[0] * eq_rx + combiners[1] * eq_ry);
    let E_final_eval = evaluations[0];
    let W_eval = evaluations[1];
    let first_sum_check_final_eval =
        (Az_claimed_val * Bz_claimed_val) - (u * Cz_claimed_val + E_final_eval);

    initial_sum_check_verification(
        &transcript.first_sum_check_transcript,
        first_sum_check_final_eval,
        initial_sum_check_random_points,
        tau,
    );

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
        &par_sum_check_random_points,
    );
    // Verifying the claimed values.
    batch_opening(
        vec![A_claimed_val, B_claimed_val, C_claimed_val],
        transcript.BatchProof,
        setup,
        eR1CSCommitments,
        transcript.rx_basis_commit,
        transcript.ry_basis_commit,
        channel,
    );
}
fn multilinear_batch_open(
    proof: BatchMultilinearKZG2Proof,
    random_points: Vec<Vec<Scalar>>,
    verifier_key: G2AffinePoint,
    phi_k_commitments: Vec<AffinePoint>,
    commits: Vec<AffinePoint>,
    channel: &mut Channel,
    combiners: &Vec<Scalar>,
) -> Vec<Scalar> {
    let combined_commit =
        commits[2].to_projective().mul(combiners[0]) + commits[3].to_projective().mul(combiners[1]);
    batch_verifier(
        proof,
        random_points,
        verifier_key,
        phi_k_commitments,
        &[commits[0], commits[1], combined_commit.to_affine()].to_vec(),
        channel,
    )
    .unwrap()
}
pub fn initial_sum_check_verification(
    transcript: &InitialSumCheckTranscript,
    final_evaluation: Scalar,
    random_points: Vec<Scalar>,
    tau: Vec<Scalar>,
) {
    let polynomials = &transcript.polynomials;
    let rounds = polynomials.len();
    let mut current_sum = Scalar::ZERO;

    for i in 0..rounds {
        let poly = polynomials[i].get_coefficients();
        assert_eq!(
            current_sum,
            eval(poly, Scalar::ZERO) + eval(poly, Scalar::ONE),
            "f(0) + f(1) did not match binding at round {:?} in the eR1CS initial sum check",
            i
        );

        current_sum = eval(poly, random_points[i])
    }
    let eq = evaluate_eq(&tau, &random_points);
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
    random_points: &Vec<Scalar>,
) {
    let polynomials = &transcript.polynomials;
    let rounds = polynomials.len();

    let mut current_sum = initial_evaluation;

    for i in 0..rounds {
        let poly = polynomials[i].get_coefficients();
        assert_eq!(
            current_sum,
            eval(poly, Scalar::ZERO) + eval(poly, Scalar::ONE),
            "f(0) + f(1) did not match binding at round {:?} in the second eR1CS sum check",
            i
        );
        current_sum = eval(&poly, random_points[i])
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

pub fn evaluate_eq(basis_point: &Vec<Scalar>, evaluation_point: &Vec<Scalar>) -> Scalar {
    let mut res = Scalar::ONE;

    for (x, y) in basis_point.iter().zip(evaluation_point.iter()) {
        res *= Scalar::ONE - *x - *y + (*x * *y).double()
    }

    res
}
