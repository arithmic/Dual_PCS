#![allow(non_snake_case)]
#![allow(unused)]
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use channel::Channel;
use multilinear_kzg::{
    common::{setup, VerificationKey, SRS},
    prover::commit,
};
use polynomial::MultPolynomial;
use preprocessing::{compute_coeff, ColumnData, SparseCommit, SparseRep};
use prover::er1cs::{prove_sat, sparse_matrix_multiply};
use rand::Rng;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
    IntoParallelRefMutIterator, ParallelIterator,
};
use spartan_common::{eR1CSCommitments, eR1CSmetadata};
use std::{collections::HashMap, time::Instant};
use traits::traits::Field;
use verifier::er1cs::verify_sat;

use crate::{preprocessing, prover, spartan_common, verifier};

#[test]
pub fn er1cs_test() {
    let num_const = 1 << 7;
    let num_inputs = 10;
    let num_var = num_const - 1;
    let sparsity = 1;
    let toxic_waste = (0..((num_const as u32).trailing_zeros() as usize + sparsity))
        .into_par_iter()
        .map(|_| Scalar::random())
        .collect();
    let (srs, ver_key) = setup::<BlsCurve>(toxic_waste);
    let (A, B, C, z, E, W, u, PI) =
        construct_matrices(sparsity as usize, num_const, num_var as usize, num_inputs);
    let (er1cs_metadata, er1cs_commitments) = er1cs_commit(
        A.clone(),
        B.clone(),
        C.clone(),
        E.clone(),
        W.clone(),
        &srs,
        sparsity,
    );
    let mut channel = Channel::initialize_with_affine_point(
        [
            er1cs_commitments.E.commitment.to_affine(),
            er1cs_commitments.W.commitment.to_affine(),
        ]
        .as_ref(),
    );

    let time = Instant::now();
    let er1cs_transcript = prove_sat(
        A,
        B,
        C,
        &u,
        &MultPolynomial::new(z),
        &MultPolynomial::new(E),
        &MultPolynomial::new(W),
        er1cs_metadata,
        &srs,
        &mut channel,
    );
    println!("Time to generate er1cs proof is {:?}", time.elapsed());
    println!(
        "Proof size {:?}",
        er1cs_transcript.to_bytes().len() as f64 / 1024f64
    );
    let mut channel = Channel::initialize_with_affine_point(
        [
            er1cs_commitments.E.commitment.to_affine(),
            er1cs_commitments.W.commitment.to_affine(),
        ]
        .as_ref(),
    );

    let pi_indices: Vec<usize> = (0..1 << 5).collect();
    let time = Instant::now();
    verify_sat(
        er1cs_transcript,
        er1cs_commitments,
        u,
        MultPolynomial::new(PI),
        pi_indices,
        &ver_key,
        &mut channel,
    );
    println!("Time to verify er1cs proof is {:?}", time.elapsed());
}
#[allow(unused)]
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
    E: Vec<Scalar>,
    W: Vec<Scalar>,
    srs: &SRS<BlsCurve>,
    sparsity: usize,
) -> (eR1CSmetadata, eR1CSCommitments) {
    let A_metadata = A.get_metadata(sparsity);
    let B_metadata = B.get_metadata(sparsity);
    let C_metadata = C.get_metadata(sparsity);
    let start_time = Instant::now();
    let E_commit = commit(&E, srs);
    let W_commit = commit(&W, srs);
    let A_commit = SparseCommit::new(&A_metadata, srs);
    let B_commit = SparseCommit::new(&B_metadata, srs);
    let C_commit = SparseCommit::new(&C_metadata, srs);
    println!("er1cs commit time {:?}", start_time.elapsed());
    let er1cs_commitments = eR1CSCommitments::new(A_commit, B_commit, C_commit, E_commit, W_commit);
    let commit_size = er1cs_commitments.to_bytes();
    println!("Commit size {:?}", commit_size);

    let er1cs_metadata = eR1CSmetadata::new(A_metadata, B_metadata, C_metadata);
    (er1cs_metadata, er1cs_commitments)
}
pub fn main(
    A: SparseRep,
    B: SparseRep,
    C: SparseRep,
    Z: Vec<Scalar>,
    E: Vec<Scalar>,
    W: Vec<Scalar>,
    u: Scalar,
    PI: Vec<Scalar>,
    srs: &SRS<BlsCurve>,
    ver_key: &VerificationKey<BlsCurve>,
    num_cols: usize,
) {
    let (er1cs_metadata, er1cs_commitments) = er1cs_commit(
        A.clone(),
        B.clone(),
        C.clone(),
        E.clone(),
        W.clone(),
        srs,
        num_cols,
    );
    let mut channel = Channel::initialize_with_affine_point(
        [
            er1cs_commitments.E.commitment.to_affine(),
            er1cs_commitments.W.commitment.to_affine(),
        ]
        .as_ref(),
    );

    let time = Instant::now();
    let er1cs_transcript = prove_sat(
        A,
        B,
        C,
        &u,
        &MultPolynomial::new(Z),
        &MultPolynomial::new(E),
        &MultPolynomial::new(W),
        er1cs_metadata,
        srs,
        &mut channel,
    );
    println!("Time to generate er1cs proof is {:?}", time.elapsed());
    println!(
        "Proof size {:?}",
        er1cs_transcript.to_bytes().len() as f64 / 1024f64
    );
    let mut channel = Channel::initialize_with_affine_point(
        [
            er1cs_commitments.E.commitment.to_affine(),
            er1cs_commitments.W.commitment.to_affine(),
        ]
        .as_ref(),
    );

    let pi_indices: Vec<usize> = (0..1 << 5).collect();
    let time = Instant::now();
    verify_sat(
        er1cs_transcript,
        er1cs_commitments,
        u,
        MultPolynomial::new(PI),
        pi_indices,
        ver_key,
        &mut channel,
    );
    println!("Time to verify er1cs proof is {:?}", time.elapsed());
}
