use bls381::scalar::Scalar;
use bls_curve::bls::AffinePoint;
use channel::Channel;
use grand_product_with_air::verifier::gkr_verifier;
use kzg_fft::verifier::kzg2_batch_verify;
use kzg_fourier_multilinear::{common::KZGFourierDegreeBoundSetup, verifier::batch_verify};
use spartan_common::{
    eR1CSCommitments, eval, BatchSparseEvalProof, BatchSpartanSumCheckTranscript,
};
use traits::traits::{Field, PrimeField};

use crate::spartan_common;

#[allow(non_snake_case)]
pub fn batch_opening(
    evaluations: Vec<Scalar>,
    eval_proof: BatchSparseEvalProof,
    setup: &KZGFourierDegreeBoundSetup,
    eR1CSCommitments: eR1CSCommitments,
    rx_basis_commit: AffinePoint,
    ry_basis_commit: AffinePoint,
    channel: &mut Channel,
) {
    let e_rx_evals_sum_check = eval_proof.e_rx_evals_sum_check.clone();
    let e_ry_evals_sum_check = eval_proof.e_ry_evals_sum_check.clone();
    let val_evals_sum_check = eval_proof.val_evals_sum_check.clone();

    let random_coeffs = channel.get_random_points(3);

    let evaluation = random_coeffs
        .iter()
        .zip(evaluations.iter())
        .map(|(coeff, eval)| *coeff * *eval)
        .reduce(|acc, g| acc + g)
        .unwrap();

    let val_commits = eR1CSCommitments.preprocess_commits.val_commits.clone();

    let poly_evaluations = BatchSumCheckEvals::new(
        e_rx_evals_sum_check.clone(),
        e_ry_evals_sum_check.clone(),
        val_evals_sum_check.clone(),
    );

    batch_spartan_sum_check_verification(
        evaluation,
        &eval_proof.sum_check_transcript,
        poly_evaluations,
        random_coeffs,
        channel,
    );

    let gkr_transcript1 = eval_proof.clone().gkr_transcript1;
    let gkr_transcript2 = eval_proof.clone().gkr_transcript2;
    let gkr_evaluations = eval_proof.gkr_evaluations;
    let gkr_commits = eR1CSCommitments.preprocess_commits;

    let n_circuits = gkr_transcript1.get_final_layers().len();
    let transcript1_final_layers = gkr_transcript1.get_final_layers();
    let transcript2_final_layers = gkr_transcript2.get_final_layers();
    for idx in 0..n_circuits / 2 {
        let w_init = transcript1_final_layers[idx];
        let s = transcript1_final_layers[n_circuits / 2 + idx];
        let w_update = transcript2_final_layers[idx];
        let r = transcript2_final_layers[n_circuits / 2 + idx];

        assert_eq!(
            w_init * w_update,
            r * s,
            "assertion failed for circuit {}",
            idx
        )
    }

    let gamma_tau = channel.get_random_points(2);

    channel.reseed_with_affine_point(gkr_transcript1.get_trace_commits());

    let constraint_comp_coeffs1 = channel.get_random_points(6 * n_circuits);

    channel.reseed_with_affine_point(&[gkr_transcript1.composition_poly_commit()].to_vec());

    //Draw a random point using Fiat-Shamir
    let z1 = channel.get_random_point();
    let trace_length1 = gkr_transcript1.get_trace_length();
    let trace_length2 = gkr_transcript2.get_trace_length();

    let g_trace1 = Scalar::get_root_of_unity(trace_length1.trailing_zeros());
    let g_trace2 = Scalar::get_root_of_unity(trace_length2.trailing_zeros());
    let g_z1 = g_trace1 * z1;
    let idx_at_z1 = gkr_evaluations.idx_at_z1;
    let idx_at_g_z1 = gkr_evaluations.idx_at_g_z1;
    let rx_basis_evals_evaluation_at_z1 = gkr_evaluations.rx_basis_evals_evaluation_at_z1;
    let rx_basis_evals_evaluation_at_g_z1 = gkr_evaluations.rx_basis_evals_evaluation_at_g_z1;
    let ry_basis_evals_evaluation_at_z1 = gkr_evaluations.ry_basis_evals_evaluation_at_z1;
    let ry_basis_evals_evaluation_at_g_z1 = gkr_evaluations.ry_basis_evals_evaluation_at_g_z1;
    let final_ts_for_rows_evaluations_at_z1 = gkr_evaluations.final_ts_for_rows_evaluations_at_z1;
    let final_ts_for_cols_evaluations_at_z1 = gkr_evaluations.final_ts_for_cols_evaluations_at_z1;
    let final_ts_for_rows_evaluations_at_g_z1 =
        gkr_evaluations.final_ts_for_rows_evaluations_at_g_z1;
    let final_ts_for_cols_evaluations_at_g_z1 =
        gkr_evaluations.final_ts_for_cols_evaluations_at_g_z1;

    let gamma_square = gamma_tau[0].square();
    let term1 = (gamma_square * idx_at_z1) - gamma_tau[1];
    let term2 = gamma_tau[0] * rx_basis_evals_evaluation_at_z1;
    let term3 = term1 + term2;

    let term4 = gamma_tau[0] * ry_basis_evals_evaluation_at_z1;
    let term5 = term1 + term4;

    let w_init_at_z = [vec![term3; n_circuits / 4], vec![term5; n_circuits / 4]].concat();

    let s_at_z = [
        (0..final_ts_for_rows_evaluations_at_z1.len())
            .map(|idx| term3 + final_ts_for_rows_evaluations_at_z1[idx])
            .collect::<Vec<Scalar>>(),
        (0..final_ts_for_cols_evaluations_at_z1.len())
            .map(|idx| term5 + final_ts_for_cols_evaluations_at_z1[idx])
            .collect::<Vec<Scalar>>(),
    ]
    .concat();

    let term1 = (gamma_square * idx_at_g_z1) - gamma_tau[1];
    let term2 = gamma_tau[0] * rx_basis_evals_evaluation_at_g_z1;
    let term3 = term1 + term2;

    let term4 = gamma_tau[0] * ry_basis_evals_evaluation_at_g_z1;
    let term5 = term1 + term4;

    let w_init_at_g_z = [vec![term3; n_circuits / 4], vec![term5; n_circuits / 4]].concat();

    let s_at_g_z = [
        (0..final_ts_for_rows_evaluations_at_g_z1.len())
            .map(|idx| term3 + final_ts_for_rows_evaluations_at_g_z1[idx])
            .collect::<Vec<Scalar>>(),
        (0..final_ts_for_cols_evaluations_at_g_z1.len())
            .map(|idx| term5 + final_ts_for_cols_evaluations_at_g_z1[idx])
            .collect::<Vec<Scalar>>(),
    ]
    .concat();

    let final_ts_for_rows_commits = gkr_commits.final_ts_for_rows_commits;
    let final_ts_for_cols_commits = gkr_commits.final_ts_for_cols_commits;

    eval_proof
        .gkr_transcript1
        .get_odd_frame()
        .reseed_with_ood_frame(channel);

    let poly_combiners1 = channel.get_random_points(n_circuits + 1);

    channel.reseed_with_affine_point(gkr_transcript2.get_trace_commits());

    let constraint_comp_coeffs2 = channel.get_random_points(6 * n_circuits);
    channel.reseed_with_affine_point(&[gkr_transcript2.composition_poly_commit()].to_vec());

    //Draw a random point using Fiat-Shamir
    let z2 = channel.get_random_point();
    let g_z2 = g_trace2 * z2;

    let rows_evaluations_at_z2 = gkr_evaluations.rows_evaluations_at_z2;
    let rows_evaluations_at_g_z2 = gkr_evaluations.rows_evaluations_at_g_z2;

    let cols_evaluations_at_z2 = gkr_evaluations.cols_evaluations_at_z2;
    let cols_evaluations_at_g_z2 = gkr_evaluations.cols_evaluations_at_g_z2;

    let e_rx_evaluations_at_z2 = gkr_evaluations.e_rx_evaluations_at_z2;
    let e_rx_evaluations_at_g_z2 = gkr_evaluations.e_rx_evaluations_at_g_z2;

    let e_ry_evaluations_at_z2 = gkr_evaluations.e_ry_evaluations_at_z2;
    let e_ry_evaluations_at_g_z2 = gkr_evaluations.e_ry_evaluations_at_g_z2;

    let read_ts_for_rows_evaluations_at_z2 = gkr_evaluations.read_ts_for_rows_evaluations_at_z2;
    let read_ts_for_rows_evaluations_at_g_z2 = gkr_evaluations.read_ts_for_rows_evaluations_at_g_z2;

    let read_ts_for_cols_evaluations_at_z2 = gkr_evaluations.read_ts_for_cols_evaluations_at_z2;
    let read_ts_for_cols_evaluations_at_g_z2 = gkr_evaluations.read_ts_for_cols_evaluations_at_g_z2;

    let result1: (Vec<Scalar>, Vec<Scalar>) = (0..rows_evaluations_at_z2.len())
        .map(|idx| {
            let term1 = gamma_square * rows_evaluations_at_z2[idx];
            let term2 = gamma_tau[0] * e_rx_evaluations_at_z2[idx];
            let term3 = term1 + term2 + read_ts_for_rows_evaluations_at_z2[idx];
            let term4 = term3 - gamma_tau[1];
            (term4, term4 + Scalar::ONE)
        })
        .unzip();
    let result2: (Vec<Scalar>, Vec<Scalar>) = (0..cols_evaluations_at_z2.len())
        .map(|idx| {
            let term1 = gamma_square * cols_evaluations_at_z2[idx];
            let term2 = gamma_tau[0] * e_ry_evaluations_at_z2[idx];
            let term3 = term1 + term2 + read_ts_for_cols_evaluations_at_z2[idx];
            let term4 = term3 - gamma_tau[1];
            (term4, term4 + Scalar::ONE)
        })
        .unzip();

    let (r_at_z, w_update_at_z): (Vec<Scalar>, Vec<Scalar>) = (
        [result1.0, result2.0].concat(),
        [result1.1, result2.1].concat(),
    );

    let result1: (Vec<Scalar>, Vec<Scalar>) = (0..rows_evaluations_at_g_z2.len())
        .map(|idx| {
            let term1 = gamma_square * rows_evaluations_at_g_z2[idx];
            let term2 = gamma_tau[0] * e_rx_evaluations_at_g_z2[idx];
            let term3 = term1 + term2 + read_ts_for_rows_evaluations_at_g_z2[idx];
            let term4 = term3 - gamma_tau[1];
            (term4, term4 + Scalar::ONE)
        })
        .unzip();
    let result2: (Vec<Scalar>, Vec<Scalar>) = (0..cols_evaluations_at_g_z2.len())
        .map(|idx| {
            let term1 = gamma_square * cols_evaluations_at_g_z2[idx];
            let term2 = gamma_tau[0] * e_ry_evaluations_at_g_z2[idx];
            let term3 = term1 + term2 + read_ts_for_cols_evaluations_at_g_z2[idx];
            let term4 = term3 - gamma_tau[1];
            (term4, term4 + Scalar::ONE)
        })
        .unzip();

    let (r_at_g_z, w_update_at_g_z): (Vec<Scalar>, Vec<Scalar>) = (
        [result1.0, result2.0].concat(),
        [result1.1, result2.1].concat(),
    );

    let rows_commits = gkr_commits.rows_commits;
    let cols_commits = gkr_commits.cols_commits;
    let e_rx_commits = eval_proof.e_rx_commits;
    let e_ry_commits = eval_proof.e_ry_commits;
    let read_ts_for_rows_commits = gkr_commits.read_ts_for_rows_commits;
    let read_ts_for_cols_commits = gkr_commits.read_ts_for_cols_commits;

    eval_proof
        .gkr_transcript2
        .get_odd_frame()
        .reseed_with_ood_frame(channel);

    let poly_combiners2 = channel.get_random_points(n_circuits + 1);

    let mut commits1_next = Vec::new();
    (0..final_ts_for_rows_commits.len())
        .for_each(|idx| commits1_next.push(final_ts_for_rows_commits[idx]));
    (0..final_ts_for_cols_commits.len())
        .for_each(|idx| commits1_next.push(final_ts_for_cols_commits[idx]));
    let trace_commits = gkr_transcript1.get_trace_commits();
    let composition_commits = gkr_transcript1.composition_poly_commit();
    (0..trace_commits.len()).for_each(|idx| commits1_next.push(trace_commits[idx]));
    commits1_next.push(gkr_commits.index_commit);
    commits1_next.push(rx_basis_commit);
    commits1_next.push(ry_basis_commit);
    let mut commits1_current = commits1_next.clone();
    commits1_current.push(composition_commits);

    let mut eval1_at_z1 = Vec::new();
    (0..final_ts_for_rows_evaluations_at_z1.len())
        .for_each(|idx| eval1_at_z1.push(final_ts_for_rows_evaluations_at_z1[idx]));
    (0..final_ts_for_cols_evaluations_at_z1.len())
        .for_each(|idx| eval1_at_z1.push(final_ts_for_cols_evaluations_at_z1[idx]));
    let trace_frame = gkr_transcript1.get_odd_frame().get_current_frame();
    let composition_frame = gkr_transcript1.get_odd_frame().get_composition_frame();
    (0..trace_frame.len()).for_each(|idx| eval1_at_z1.push(trace_frame[idx]));
    eval1_at_z1.push(idx_at_z1);
    eval1_at_z1.push(rx_basis_evals_evaluation_at_z1);
    eval1_at_z1.push(ry_basis_evals_evaluation_at_z1);
    eval1_at_z1.push(composition_frame);

    let mut eval1_at_g_z1 = Vec::new();
    (0..final_ts_for_rows_evaluations_at_g_z1.len())
        .for_each(|idx| eval1_at_g_z1.push(final_ts_for_rows_evaluations_at_g_z1[idx]));
    (0..final_ts_for_cols_evaluations_at_g_z1.len())
        .for_each(|idx| eval1_at_g_z1.push(final_ts_for_cols_evaluations_at_g_z1[idx]));
    let next_frame = gkr_transcript1.get_odd_frame().get_next_frame();
    (0..next_frame.len()).for_each(|idx| eval1_at_g_z1.push(next_frame[idx]));
    eval1_at_g_z1.push(idx_at_g_z1);
    eval1_at_g_z1.push(rx_basis_evals_evaluation_at_g_z1);
    eval1_at_g_z1.push(ry_basis_evals_evaluation_at_g_z1);
    let combiners1 = channel.get_random_points(eval1_at_z1.len());

    ////
    let mut commits2_next = Vec::new();
    (0..rows_commits.len()).for_each(|idx| commits2_next.push(rows_commits[idx]));
    (0..cols_commits.len()).for_each(|idx| commits2_next.push(cols_commits[idx]));
    (0..e_rx_commits.len()).for_each(|idx| commits2_next.push(e_rx_commits[idx]));
    (0..e_ry_commits.len()).for_each(|idx| commits2_next.push(e_ry_commits[idx]));
    (0..read_ts_for_rows_commits.len())
        .for_each(|idx| commits2_next.push(read_ts_for_rows_commits[idx]));
    (0..read_ts_for_cols_commits.len())
        .for_each(|idx| commits2_next.push(read_ts_for_cols_commits[idx]));
    let trace_commits = gkr_transcript2.get_trace_commits();
    (0..trace_commits.len()).for_each(|idx| commits2_next.push(trace_commits[idx]));
    let mut commits2_current = commits2_next.clone();
    commits2_current.push(gkr_transcript2.composition_poly_commit());

    let mut eval1_at_z2 = Vec::new();
    (0..rows_evaluations_at_z2.len()).for_each(|idx| eval1_at_z2.push(rows_evaluations_at_z2[idx]));
    (0..cols_evaluations_at_z2.len()).for_each(|idx| eval1_at_z2.push(cols_evaluations_at_z2[idx]));
    (0..e_rx_evaluations_at_z2.len()).for_each(|idx| eval1_at_z2.push(e_rx_evaluations_at_z2[idx]));
    (0..e_ry_evaluations_at_z2.len()).for_each(|idx| eval1_at_z2.push(e_ry_evaluations_at_z2[idx]));
    (0..read_ts_for_rows_evaluations_at_z2.len())
        .for_each(|idx| eval1_at_z2.push(read_ts_for_rows_evaluations_at_z2[idx]));
    (0..read_ts_for_cols_evaluations_at_z2.len())
        .for_each(|idx| eval1_at_z2.push(read_ts_for_cols_evaluations_at_z2[idx]));

    let trace_frame = gkr_transcript2.get_odd_frame().get_current_frame();
    let composition_frame = gkr_transcript2.get_odd_frame().get_composition_frame();
    (0..trace_frame.len()).for_each(|idx| eval1_at_z2.push(trace_frame[idx]));
    eval1_at_z2.push(composition_frame);

    let mut eval1_at_g_z2 = Vec::new();
    (0..rows_evaluations_at_g_z2.len())
        .for_each(|idx| eval1_at_g_z2.push(rows_evaluations_at_g_z2[idx]));
    (0..cols_evaluations_at_g_z2.len())
        .for_each(|idx| eval1_at_g_z2.push(cols_evaluations_at_g_z2[idx]));
    (0..e_rx_evaluations_at_g_z2.len())
        .for_each(|idx| eval1_at_g_z2.push(e_rx_evaluations_at_g_z2[idx]));
    (0..e_ry_evaluations_at_g_z2.len())
        .for_each(|idx| eval1_at_g_z2.push(e_ry_evaluations_at_g_z2[idx]));
    (0..read_ts_for_rows_evaluations_at_g_z2.len())
        .for_each(|idx| eval1_at_g_z2.push(read_ts_for_rows_evaluations_at_g_z2[idx]));
    (0..read_ts_for_cols_evaluations_at_g_z2.len())
        .for_each(|idx| eval1_at_g_z2.push(read_ts_for_cols_evaluations_at_g_z2[idx]));

    let next_frame = gkr_transcript2.get_odd_frame().get_next_frame();
    (0..next_frame.len()).for_each(|idx| eval1_at_g_z2.push(next_frame[idx]));

    kzg2_batch_verify(
        eval_proof.proof1,
        setup.get_setup(trace_length1).setup.verifier_key,
        z1,
        commits1_current,
        eval1_at_z1,
        combiners1.clone(),
    );

    kzg2_batch_verify(
        eval_proof.proof2,
        setup.get_setup(trace_length1).setup.verifier_key,
        g_z1,
        commits1_next,
        eval1_at_g_z1,
        combiners1[0..combiners1.len() - 1].to_vec(),
    );
    let combiners = channel.get_random_points(eval1_at_z2.len());

    kzg2_batch_verify(
        eval_proof.proof3,
        setup.get_setup(trace_length2).setup.verifier_key,
        z2,
        commits2_current,
        eval1_at_z2,
        combiners.clone(),
    );

    kzg2_batch_verify(
        eval_proof.proof4,
        setup.get_setup(trace_length2).setup.verifier_key,
        g_z2,
        commits2_next,
        eval1_at_g_z2,
        combiners[0..combiners.len() - 1].to_vec(),
    );
    gkr_verifier(
        gkr_transcript1,
        constraint_comp_coeffs1,
        z1,
        poly_combiners1,
        n_circuits,
        trace_length1,
        setup.get_setup(trace_length1).setup.verifier_key,
        [w_init_at_z, s_at_z].concat(),
        [w_init_at_g_z, s_at_g_z].concat(),
    );
    gkr_verifier(
        gkr_transcript2,
        constraint_comp_coeffs2,
        z2,
        poly_combiners2,
        n_circuits,
        trace_length2,
        setup.get_setup(trace_length2).setup.verifier_key,
        [w_update_at_z, r_at_z].concat(),
        [w_update_at_g_z, r_at_g_z].concat(),
    );

    let batch_eval_coeffs_sum_check = channel.get_random_points(3 * 3);
    let random_points = eval_proof.sum_check_transcript.clone().random_points;
    let degree = 1 << random_points.len();
    batch_verify(
        eval_proof.sum_check_eval_proof,
        [e_rx_commits, e_ry_commits, val_commits].concat(),
        [
            e_rx_evals_sum_check,
            e_ry_evals_sum_check,
            val_evals_sum_check,
        ]
        .concat(),
        batch_eval_coeffs_sum_check,
        random_points,
        setup.setup.get_setup(degree).verifier_key,
        setup.get_setup(degree).phi_commitment.clone(),
        channel,
    );
}

pub fn batch_spartan_sum_check_verification(
    evaluation: Scalar,
    transcript: &BatchSpartanSumCheckTranscript,
    poly_evaluations: BatchSumCheckEvals,
    random_coeffs: Vec<Scalar>,
    channel: &mut Channel,
) {
    let polynomials = &transcript.polynomials;
    let rounds = polynomials.len();

    let mut current_sum = evaluation;

    let mut r = vec![Scalar::ZERO; rounds];
    for i in 0..rounds {
        let poly = polynomials[i].get_coefficients();
        assert_eq!(
            current_sum,
            eval(poly, Scalar::ZERO) + eval(poly, Scalar::ONE),
            "f(0) + f(1) did not match binding at round {:?}",
            i
        );

        channel.reseed_with_scalars(poly);
        let r_i = channel.get_random_point();
        r[i] = r_i;
        current_sum = eval(poly, r_i)
    }

    let mut final_eval = Scalar::ZERO;

    for i in 0..3 {
        final_eval += random_coeffs[i]
            * (poly_evaluations.e_rx_evals[i]
                * poly_evaluations.e_ry_evals[i]
                * poly_evaluations.val_evals[i])
    }
    assert_eq!(
        current_sum, final_eval,
        "Final sum check verification failed in Spartan"
    )
}

pub struct BatchSumCheckEvals {
    pub e_rx_evals: Vec<Scalar>,
    pub e_ry_evals: Vec<Scalar>,
    pub val_evals: Vec<Scalar>,
}

impl BatchSumCheckEvals {
    pub fn new(
        e_rx_evals: Vec<Scalar>,
        e_ry_evals: Vec<Scalar>,
        val_evals: Vec<Scalar>,
    ) -> BatchSumCheckEvals {
        BatchSumCheckEvals {
            e_rx_evals,
            e_ry_evals,
            val_evals,
        }
    }
}
