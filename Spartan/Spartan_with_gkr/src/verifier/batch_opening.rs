#![allow(non_snake_case)]
#![allow(unused)]
use crate::{preprocessing, spartan_common, verifier::er1cs::eval};
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use channel::Channel;
use grand_product_with_gkr::verifier::gkr_verifier;
use multilinear_kzg::{
    common::{MleCommit, MleEvalProof},
    verifier::batch_commits,
};
use preprocessing::{compute_coeff, SparseCommit};
use spartan_common::{BatchSparseEvalProof, BatchSpartanSumCheckTranscript};
use traits::traits::Field;
pub fn batch_opening(
    evaluations: Vec<Scalar>,
    eval_proof: BatchSparseEvalProof,
    commitments: Vec<SparseCommit>,
    eval_point: &Vec<Scalar>,
    channel: &mut Channel,
) -> (
    MleCommit<BlsCurve>,
    Scalar,
    MleEvalProof<BlsCurve>,
    Vec<Scalar>,
    MleCommit<BlsCurve>,
    MleCommit<BlsCurve>,
    Scalar,
    Scalar,
    MleEvalProof<BlsCurve>,
    MleEvalProof<BlsCurve>,
    Vec<Scalar>,
    Vec<Scalar>,
) {
    let e_rx_evals_sum_check = eval_proof.e_rx_evals_sum_check.clone();
    let e_ry_evals_sum_check = eval_proof.e_ry_evals_sum_check.clone();
    let val_evals_sum_check = eval_proof.val_evals_sum_check.clone();

    let (rx, ry) = eval_point.split_at(eval_point.len() / 2);
    let rx_basis_evals = compute_coeff(&rx.to_vec());
    let ry_basis_evals = compute_coeff(&ry.to_vec());

    let poly_evaluations = BatchSumCheckEvals::new(
        e_rx_evals_sum_check,
        e_ry_evals_sum_check,
        val_evals_sum_check,
    );

    let random_coeffs = channel.get_random_points(3);
    let evaluation = random_coeffs
        .iter()
        .zip(evaluations.iter())
        .map(|(coeff, eval)| *coeff * *eval)
        .reduce(|acc, g| acc + g)
        .unwrap();

    batch_spartan_sum_check_verification(
        evaluation,
        &eval_proof.sum_check_transcript,
        poly_evaluations,
        random_coeffs,
        channel,
    );

    let e_rx_commits = eval_proof.e_rx_commits.clone();
    let e_ry_commits = eval_proof.e_ry_commits.clone();
    let val_commits = commitments.iter().map(|Commit| Commit.val_commit).collect();

    let sum_check_e_rx_evals = eval_proof.e_rx_evals_sum_check;
    let sum_check_e_ry_evals = eval_proof.e_ry_evals_sum_check;
    let val_evals = eval_proof.val_evals_sum_check;

    let final_point = &eval_proof.gkr_transcript1.final_layer_point;
    let basis_evals_for_binding = compute_coeff(final_point);

    let eq_rx_eval = basis_evals_for_binding
        .iter()
        .zip(rx_basis_evals.iter())
        .fold(Scalar::ZERO, |acc, (eq_val, basis_val)| {
            acc + (*eq_val * *basis_val)
        });

    let eq_ry_eval = basis_evals_for_binding
        .iter()
        .zip(ry_basis_evals.iter())
        .fold(Scalar::ZERO, |acc, (eq_val, basis_val)| {
            acc + (*eq_val * *basis_val)
        });

    let gamma_tau = channel.get_random_points(2);

    let indices_eval = evaluate_indicies(final_point);

    let final_ts_evals_row_mem_check = eval_proof.final_ts_evals_row_mem_check;
    let final_ts_evals_col_mem_check = eval_proof.final_ts_evals_col_mem_check;
    let row_evals_mem_check = eval_proof.row_evals_mem_check;
    let col_evals_mem_check = eval_proof.col_evals_mem_check;
    let read_ts_evals_row_mem_check = eval_proof.read_ts_evals_row_mem_check;
    let read_ts_evals_col_mem_check = eval_proof.read_ts_evals_col_mem_check;
    let e_rx_evals = eval_proof.e_rx_evals_mem_check;
    let e_ry_evals = eval_proof.e_ry_evals_mem_check;

    let mut w_init_evaluations = vec![Scalar::ZERO; 6];
    let mut s_evaluations = vec![Scalar::ZERO; 6];
    let mut r_evaluations = vec![Scalar::ZERO; 6];
    let mut w_update_evaluations = vec![Scalar::ZERO; 6];

    for c in 0..3 {
        w_init_evaluations[c] =
            gamma_tau[0].square() * indices_eval + gamma_tau[0] * eq_rx_eval - gamma_tau[1];

        s_evaluations[c] = w_init_evaluations[c] + final_ts_evals_row_mem_check[c];

        r_evaluations[c] = gamma_tau[0].square() * row_evals_mem_check[c]
            + gamma_tau[0] * e_rx_evals[c]
            + read_ts_evals_row_mem_check[c]
            - gamma_tau[1];

        w_update_evaluations[c] = r_evaluations[c] + Scalar::ONE;

        w_init_evaluations[c + 3] =
            gamma_tau[0].square() * indices_eval + gamma_tau[0] * eq_ry_eval - gamma_tau[1];

        s_evaluations[c + 3] = w_init_evaluations[c + 3] + final_ts_evals_col_mem_check[c];

        r_evaluations[c + 3] = gamma_tau[0].square() * col_evals_mem_check[c]
            + gamma_tau[0] * e_ry_evals[c]
            + read_ts_evals_col_mem_check[c]
            - gamma_tau[1];

        w_update_evaluations[c + 3] = r_evaluations[c + 3] + Scalar::ONE;
    }

    let depth_1 = rx.len();

    let n_circuits = w_init_evaluations.len() + s_evaluations.len();
    gkr_verifier(
        &eval_proof.gkr_transcript1,
        depth_1,
        channel,
        [w_init_evaluations, s_evaluations].concat(),
        n_circuits,
    );
    let depth_2 = eval_proof.circuit2_depth;
    let n_circuits = w_update_evaluations.len() + r_evaluations.len();
    gkr_verifier(
        &eval_proof.gkr_transcript2,
        depth_2,
        channel,
        [w_update_evaluations, r_evaluations].concat(),
        n_circuits,
    );
    let batch_eval_coeffs_sum_check = channel.get_random_points(3 * 3);

    let (sum_check_batch_commit, sum_check_batch_eval) = batch_commits(
        [e_rx_commits, e_ry_commits, val_commits].concat().as_ref(),
        [sum_check_e_rx_evals, sum_check_e_ry_evals, val_evals]
            .concat()
            .as_ref(),
        &batch_eval_coeffs_sum_check,
        &eval_proof.sum_check_eval_proof,
    );

    let batch_eval_coeffs_gkr = channel.get_random_points(3 * 8);

    let final_ts_row_commits: Vec<_> = commitments
        .iter()
        .map(|commit| commit.final_ts_row_commit)
        .collect();
    let final_ts_col_commits: Vec<_> = commitments
        .iter()
        .map(|commit| commit.final_ts_col_commit)
        .collect();
    let row_commits = commitments.iter().map(|commit| commit.row_commit).collect();
    let col_commits = commitments.iter().map(|commit| commit.col_commit).collect();
    let read_ts_row_commits = commitments
        .iter()
        .map(|commit| commit.read_ts_row_commit)
        .collect();
    let read_ts_col_commits = commitments
        .iter()
        .map(|commit| commit.read_ts_col_commit)
        .collect();
    let e_rx_commits = eval_proof.e_rx_commits;
    let e_ry_commits = eval_proof.e_ry_commits;
    let (gkr_commit1, gkr_eval1) = batch_commits(
        &[final_ts_row_commits, final_ts_col_commits].concat(),
        &[final_ts_evals_row_mem_check, final_ts_evals_col_mem_check].concat(),
        &batch_eval_coeffs_gkr,
        &eval_proof.gkr_batch_eval_proof1,
    );
    let (gkr_commit2, gkr_eval2) = batch_commits(
        &[
            row_commits,
            col_commits,
            read_ts_row_commits,
            read_ts_col_commits,
            e_rx_commits,
            e_ry_commits,
        ]
        .concat(),
        &[
            row_evals_mem_check,
            col_evals_mem_check,
            read_ts_evals_row_mem_check,
            read_ts_evals_col_mem_check,
            e_rx_evals,
            e_ry_evals,
        ]
        .concat(),
        &batch_eval_coeffs_gkr,
        &eval_proof.gkr_batch_eval_proof2,
    );

    (
        sum_check_batch_commit,
        sum_check_batch_eval,
        eval_proof.sum_check_eval_proof,
        eval_proof.sum_check_transcript.random_points,
        gkr_commit1,
        gkr_commit2,
        gkr_eval1,
        gkr_eval2,
        eval_proof.gkr_batch_eval_proof1,
        eval_proof.gkr_batch_eval_proof2,
        eval_proof.gkr_transcript1.final_layer_point,
        eval_proof.gkr_transcript2.final_layer_point,
    )
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

fn evaluate_indicies(random_values: &Vec<Scalar>) -> Scalar {
    let mut evaluation = Scalar::ZERO;
    let bits = random_values.len() - 1;
    for i in 0..random_values.len() {
        evaluation += Scalar::from(1u64 << (bits - i)) * random_values[i];
    }
    evaluation
}
