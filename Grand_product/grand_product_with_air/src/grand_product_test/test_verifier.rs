use crate::{grand_product_common, verifier::gkr_verifier};
use bls381::scalar::Scalar;
use bls_curve::bls::{AffinePoint, ProjectivePoint};
use channel::Channel;
use grand_product_common::{Commitments, Evaluations, GkrTranscript};
use kzg_fft::{
    commit::kzg2commit,
    common::{KZGFFTDegreeBoundSetup, KZGFFTEvaluationProof},
    verifier::kzg2_batch_verify,
};
use traits::traits::{Field, PrimeField};

use super::common::reseed_with_commits;

pub fn verifier(
    gkr_transcript: GkrTranscript,
    commitments: Commitments,
    gkr_evaluations: Evaluations,
    setup: &KZGFFTDegreeBoundSetup,
    proof1: KZGFFTEvaluationProof,
    proof2: KZGFFTEvaluationProof,
) {
    let n_circuits = gkr_transcript.get_final_layers().len();
    let mut channel = Channel::initialize_with_affine_point(commitments.get_A_commits());
    reseed_with_commits(&commitments, &mut channel);

    let gamma_tau = channel.get_random_points(2);

    channel.reseed_with_affine_point(gkr_transcript.get_trace_commits());

    let constraint_comp_coeffs = channel.get_random_points(6 * n_circuits);

    channel.reseed_with_affine_point(&[gkr_transcript.composition_poly_commit()].to_vec());

    //Draw a random point using Fiat-Shamir
    let z = channel.get_random_point();
    let mut leaf_layer_1_at_z = vec![Scalar::ZERO; n_circuits];
    let mut leaf_layer_1_at_gz = vec![Scalar::ZERO; n_circuits];

    let gamma_square = gamma_tau[0].square();
    let trace_length = gkr_transcript.get_trace_length();
    let g_trace = Scalar::get_root_of_unity(trace_length.trailing_zeros());
    for c in 0..n_circuits {
        let term1 = gamma_square * gkr_evaluations.A_at_z[c];
        let term2 = gamma_tau[0] * gkr_evaluations.B_at_z[c];
        leaf_layer_1_at_z[c] = term1 + term2 + gkr_evaluations.C_at_z[c] - gamma_tau[1];

        let term1 = gamma_square * gkr_evaluations.A_at_gz[c];
        let term2 = gamma_tau[0] * gkr_evaluations.B_at_gz[c];
        leaf_layer_1_at_gz[c] = term1 + term2 + gkr_evaluations.C_at_gz[c] - gamma_tau[1];
    }

    gkr_transcript
        .get_odd_frame()
        .reseed_with_ood_frame(&mut channel);

    let poly_combiners = channel.get_random_points(n_circuits + 1);
    let A_commits = commitments.get_A_commits();
    let B_commits = commitments.get_B_commits();
    let C_commits = commitments.get_C_commits();
    let mut next_polys_commits = Vec::new();
    A_commits
        .iter()
        .for_each(|commit| next_polys_commits.push(*commit));
    B_commits
        .iter()
        .for_each(|commit| next_polys_commits.push(*commit));
    C_commits
        .iter()
        .for_each(|commit| next_polys_commits.push(*commit));
    let trace_commits = gkr_transcript.get_trace_commits();
    let composition_commit = gkr_transcript.composition_poly_commit();
    trace_commits
        .iter()
        .for_each(|commit| next_polys_commits.push(*commit));
    let mut current_polys_commits = next_polys_commits.clone();
    current_polys_commits.push(composition_commit);

    let current_frame = gkr_transcript.get_odd_frame().get_current_frame();
    let next_frame = gkr_transcript.get_odd_frame().get_next_frame();
    let compositon_frame = gkr_transcript.get_odd_frame().get_composition_frame();

    let mut current_evals = Vec::new();
    gkr_evaluations
        .A_at_z
        .iter()
        .for_each(|value| current_evals.push(*value));
    gkr_evaluations
        .B_at_z
        .iter()
        .for_each(|value| current_evals.push(*value));
    gkr_evaluations
        .C_at_z
        .iter()
        .for_each(|value| current_evals.push(*value));
    current_frame
        .iter()
        .for_each(|value| current_evals.push(*value));
    current_evals.push(compositon_frame);

    let mut next_evals = Vec::new();
    gkr_evaluations
        .A_at_gz
        .iter()
        .for_each(|value| next_evals.push(*value));
    gkr_evaluations
        .B_at_gz
        .iter()
        .for_each(|value| next_evals.push(*value));
    gkr_evaluations
        .C_at_gz
        .iter()
        .for_each(|value| next_evals.push(*value));
    next_frame.iter().for_each(|value| next_evals.push(*value));

    let combiners = channel.get_random_points(current_evals.len());
    let verifier_key = setup.get_setup(trace_length).verifier_key;
    kzg2_batch_verify(
        proof1,
        verifier_key,
        z,
        current_polys_commits,
        current_evals,
        combiners.clone(),
    );

    kzg2_batch_verify(
        proof2,
        verifier_key,
        g_trace * z,
        next_polys_commits,
        next_evals,
        combiners[0..combiners.len() - 1].to_vec(),
    );
    gkr_verifier(
        gkr_transcript,
        constraint_comp_coeffs,
        z,
        poly_combiners,
        n_circuits,
        trace_length,
        setup.get_setup(trace_length).verifier_key,
        leaf_layer_1_at_z,
        leaf_layer_1_at_gz,
    );
}
#[allow(dead_code)]
fn commit_idx_vec(prover_key: &Vec<ProjectivePoint>, trace_length1: usize) -> AffinePoint {
    let idx_vec = (0..trace_length1)
        .map(|idx| Scalar::from(idx as u32))
        .collect::<Vec<Scalar>>();
    kzg2commit(&idx_vec, prover_key)
}
#[allow(dead_code)]
fn evaluate_idx_vec(trace_length1: usize, evaluation_point: Scalar) -> (Scalar, Scalar) {
    let degree = trace_length1;
    let log2_degree = degree.trailing_zeros();
    let omega = Scalar::get_root_of_unity(log2_degree);

    let g_z1 = evaluation_point * omega;

    let mut omega_powers = vec![Scalar::ONE; degree];
    for idx in 1..degree {
        let temp = omega_powers[idx - 1] * omega;
        omega_powers[idx] = temp;
    }
    //Compute alpha_i
    let degree_inverse = (Scalar::from(degree as u64)).invert().unwrap();
    let beta_z1 = (0..degree)
        .map(|i| {
            let mut temp = Scalar::ONE;
            for j in 0..log2_degree {
                temp *= Scalar::ONE
                    + (omega_powers[(degree - i) % degree] * evaluation_point).power_by([
                        1_u64 << j,
                        0,
                        0,
                        0,
                    ])
            }
            temp * degree_inverse
        })
        .collect::<Vec<Scalar>>();

    let beta_g_z1 = (0..degree)
        .map(|i| {
            let mut temp = Scalar::ONE;
            for j in 0..log2_degree {
                temp *= Scalar::ONE
                    + (omega_powers[(degree - i) % degree] * g_z1).power_by([1_u64 << j, 0, 0, 0])
            }
            temp * degree_inverse
        })
        .collect::<Vec<Scalar>>();

    let current = (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
        acc + (Scalar::from(idx as u64) * beta_z1[idx])
    });
    let next = (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
        acc + (Scalar::from(idx as u64) * beta_g_z1[idx])
    });
    (current, next)
}
