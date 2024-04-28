#![allow(non_snake_case)]
#![allow(unused)]
use bls381::scalar::Scalar;
use bls_curve::bls::{AffinePoint, BlsCurve};
use channel::Channel;
use kzg_fft::{commit::kzg2commit, common::KZGFFTDegreeBoundSetup};
use kzg_fourier_multilinear::{common::KZGFourierDegreeBoundSetup, setup::multilinearkzg2setup};
use polynomial::MultPolynomial;
use preprocessing::{compute_coeff, ColumnData, SparseRep};
use prover::er1cs::{prove_sat, sparse_matrix_multiply};
use rand::Rng;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
    IntoParallelRefMutIterator, ParallelIterator,
};
use spartan_common::{eR1CSCommitments, eR1CSmetadata, PreprocessCommits};
use std::{collections::HashMap, time::Instant};
use traits::traits::Field;
use verifier::er1cs::verify_sat;

use crate::{preprocessing, prover, spartan_common, verifier};
#[test]
pub fn er1cs_test() {
    let setup = multilinearkzg2setup(1 << 7);
    let num_const = 1 << 7;
    let num_inputs = 10;
    let num_var = num_const - 1;
    let sparsity = 1;
    let (A, B, C, z, E, W, u, PI) =
        construct_matrices(sparsity as usize, num_const, num_var as usize, num_inputs);
    let (Z, E, W, er1cs_metadata, er1cs_commitments) = er1cs_commit(
        A.clone(),
        B.clone(),
        C.clone(),
        z,
        E,
        W,
        u,
        &setup,
        sparsity as usize,
    );

    let mut channel = Channel::initialize_with_scalar(&[Scalar::ONE]);
    let er1cs_transcript = prove_sat(
        A,
        B,
        C,
        &u,
        &Z,
        &E,
        &W,
        er1cs_metadata,
        &mut channel,
        &setup,
    );
    let mut channel = Channel::initialize_with_scalar(&[Scalar::ONE]);
    let pi_indices: Vec<usize> = (0..1 << 5).collect();

    verify_sat(
        er1cs_transcript,
        u,
        MultPolynomial(PI),
        pi_indices,
        &setup,
        er1cs_commitments,
        &mut channel,
    );
}
pub fn construct_matrices(
    sparsity: usize,
    num_const: usize,
    num_var: usize,
    num_pi: usize,
) -> (
    SparseRep,
    SparseRep,
    SparseRep,
    Vec<Scalar>,
    Vec<Scalar>,
    Vec<Scalar>,
    Scalar,
    Vec<Scalar>,
) {
    let W = vec![Scalar::random(); num_const];
    let u = Scalar::random();
    let PI = vec![Scalar::random(); num_pi];
    let mut E = vec![Scalar::ZERO; num_const];

    let mut Z = W.clone();
    Z.par_iter_mut()
        .enumerate()
        .take(PI.len())
        .for_each(|(i, W)| *W += PI[i]);
    let mut A: HashMap<usize, Vec<ColumnData>> = HashMap::new();
    let mut B: HashMap<usize, Vec<ColumnData>> = HashMap::new();
    let mut C: HashMap<usize, Vec<ColumnData>> = HashMap::new();
    let Z = MultPolynomial::new(Z);
    for i in 0..num_const - 1 {
        let mut rng = rand::thread_rng();
        let A_row: Vec<ColumnData> = (0..sparsity)
            .map(|_| ColumnData::new((rng.gen_range(0..num_var - 1)) as usize, Scalar::random()))
            .collect();
        let B_row: Vec<ColumnData> = (0..sparsity)
            .map(|_| ColumnData::new((rng.gen_range(0..num_var - 1)) as usize, Scalar::random()))
            .collect();
        let C_row: Vec<ColumnData> = (0..sparsity)
            .map(|_| ColumnData::new((rng.gen_range(0..num_var - 1)) as usize, Scalar::random()))
            .collect();
        A.insert(i, A_row);
        B.insert(i, B_row);
        C.insert(i, C_row);
    }

    let A = SparseRep::new(A);
    let B = SparseRep::new(B);
    let C = SparseRep::new(C);
    let Az = sparse_matrix_multiply(&A, &Z);
    let Bz = sparse_matrix_multiply(&B, &Z);
    let Cz = sparse_matrix_multiply(&C, &Z);
    E.par_iter_mut()
        .enumerate()
        .for_each(|(i, E)| *E = Az.get_coeff(i) * Bz.get_coeff(i) - u * Cz.get_coeff(i));

    assert_eq!(
        Az.get_coeff(0) * Bz.get_coeff(0),
        u * Cz.get_coeff(0) + E[0]
    );

    let z = Z.0;
    (A, B, C, z, E, W, u, PI)
}
pub fn er1cs_commit(
    A: SparseRep,
    B: SparseRep,
    C: SparseRep,
    Z: Vec<Scalar>,
    E: Vec<Scalar>,
    W: Vec<Scalar>,
    u: Scalar,
    setup: &KZGFourierDegreeBoundSetup,
    num_cols: usize,
) -> (
    MultPolynomial,
    MultPolynomial,
    MultPolynomial,
    eR1CSmetadata,
    eR1CSCommitments,
) {
    let A_metadata = A.get_metadata(num_cols);
    let B_metadata = B.get_metadata(num_cols);
    let C_metadata = C.get_metadata(num_cols);

    let sparse_metadata = [A_metadata.clone(), B_metadata.clone(), C_metadata.clone()];
    let er1cs_metadata = eR1CSmetadata::new(A_metadata, B_metadata, C_metadata);
    let rows: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.row.as_coeffs())
        .collect();
    let cols: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.col.as_coeffs())
        .collect();
    let read_ts_for_rows: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.timestamps.read_ts_row.as_coeffs())
        .collect();
    let read_ts_for_cols: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.timestamps.read_ts_col.as_coeffs())
        .collect();
    let final_ts_for_rows: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.timestamps.final_ts_row.as_coeffs())
        .collect();
    let final_ts_for_cols: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.timestamps.final_ts_col.as_coeffs())
        .collect();
    let val: Vec<Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.val.as_coeffs().clone())
        .collect();
    let index_vec = (0..final_ts_for_rows[0].len())
        .map(|idx| Scalar::from(idx as u64))
        .collect::<Vec<Scalar>>();
    let start_time = Instant::now();
    let E_commit = kzg2commit(&E, &setup.setup.get_setup(E.len()).prover_key);
    let W_commit = kzg2commit(&W, &setup.setup.get_setup(W.len()).prover_key);

    let preprocess_commit = commit(
        &setup.setup,
        &final_ts_for_rows,
        &final_ts_for_cols,
        &rows,
        &cols,
        &read_ts_for_rows,
        &read_ts_for_cols,
        &val,
        index_vec,
    );
    let eR1CSCommitments = eR1CSCommitments::new(preprocess_commit, E_commit, W_commit);
    let commit_size = eR1CSCommitments.to_bytes().len() as f64 / 1024f64;

    let Z = MultPolynomial::new(Z);
    let E = MultPolynomial::new(E);
    let W = MultPolynomial::new(W);
    (Z, E, W, er1cs_metadata, eR1CSCommitments)
}
pub(crate) fn commit(
    setup: &KZGFFTDegreeBoundSetup,
    final_ts_for_rows: &Vec<&Vec<Scalar>>,
    final_ts_for_cols: &Vec<&Vec<Scalar>>,
    rows: &Vec<&Vec<Scalar>>,
    cols: &Vec<&Vec<Scalar>>,
    read_ts_for_rows: &Vec<&Vec<Scalar>>,
    read_ts_for_cols: &Vec<&Vec<Scalar>>,
    val: &Vec<Vec<Scalar>>,
    index_vec: Vec<Scalar>,
) -> PreprocessCommits {
    let prover_key_1 = &setup.get_setup(final_ts_for_rows[0].len()).prover_key;
    let prover_key_2 = &setup.get_setup(rows[0].len()).prover_key;

    let final_ts_for_rows_commits = (0..final_ts_for_rows.len())
        .map(|idx| kzg2commit(final_ts_for_rows[idx], prover_key_1))
        .collect::<Vec<AffinePoint>>();

    let final_ts_for_cols_commits = (0..final_ts_for_cols.len())
        .map(|idx| kzg2commit(final_ts_for_cols[idx], prover_key_1))
        .collect::<Vec<AffinePoint>>();

    let rows_commits = (0..rows.len())
        .map(|idx| kzg2commit(rows[idx], prover_key_2))
        .collect::<Vec<AffinePoint>>();
    let cols_commits = (0..cols.len())
        .map(|idx| kzg2commit(cols[idx], prover_key_2))
        .collect::<Vec<AffinePoint>>();

    let read_ts_for_rows_commits = (0..read_ts_for_rows.len())
        .map(|idx| kzg2commit(read_ts_for_rows[idx], prover_key_2))
        .collect::<Vec<AffinePoint>>();
    let read_ts_for_cols_commits = (0..read_ts_for_cols.len())
        .map(|idx| kzg2commit(read_ts_for_cols[idx], prover_key_2))
        .collect::<Vec<AffinePoint>>();

    let val_commits = (0..val.len())
        .map(|idx| kzg2commit(&val[idx], &setup.get_setup(val[idx].len()).prover_key))
        .collect::<Vec<AffinePoint>>();
    let index_commit = kzg2commit(&index_vec, &setup.get_setup(index_vec.len()).prover_key);
    PreprocessCommits {
        final_ts_for_rows_commits,
        final_ts_for_cols_commits,
        rows_commits,
        cols_commits,
        read_ts_for_rows_commits,
        read_ts_for_cols_commits,
        val_commits,
        index_commit,
    }
}
