#![allow(unused)]
#![allow(non_snake_case)]
pub mod test_prover;
pub mod test_verifier;
use std::time::Instant;
pub mod common;
use bls381::scalar::Scalar;
use kzg_fft::setup::kzg2_setup;
use rand::{thread_rng, Rng};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use traits::traits::Field;

use crate::grand_product_test::{test_prover::prover, test_verifier::verifier};
#[test]
fn test_grand_product_with_air() {
    let total_circuits = 1;
    let input_length = 1 << 10;
    let setup = kzg2_setup(input_length);
    let (A, B, C) = construct_input(total_circuits, input_length);

    let start_time = Instant::now();
    let (commitments, gkr_transcript, evaluations, proof1, proof2) =
        prover(&A, &B, &C, total_circuits, &setup);
    println!("Prover time {:?}", start_time.elapsed());

    let commitment_size = commitments.to_bytes();
    println!("Commitment_size {:?}", commitment_size);

    let proof_size = (gkr_transcript.to_bytes().len() as f64 / 1024f64)
        + (evaluations.to_bytes())
        + (proof1.to_bytes().len() as f64 / 1024f64)
        + (proof2.to_bytes().len() as f64 / 1024f64);
    println!("Proof size {:?}", proof_size);

    let start_time = Instant::now();
    verifier(
        gkr_transcript,
        commitments,
        evaluations,
        &setup,
        proof1,
        proof2,
    );
    println!("Verify time {:?}", start_time.elapsed());
}
pub fn construct_input(
    n_circuits: usize,
    input_length: usize,
) -> (Vec<Vec<Scalar>>, Vec<Vec<Scalar>>, Vec<Vec<Scalar>>) {
    let A: Vec<Vec<Scalar>> = (0..n_circuits)
        .map(|_| {
            (0..input_length)
                .into_par_iter()
                .map(|_| Scalar::random())
                .collect()
        })
        .collect();

    let B: Vec<Vec<Scalar>> = (0..n_circuits)
        .map(|_| {
            (0..input_length)
                .into_par_iter()
                .map(|_| Scalar::random())
                .collect()
        })
        .collect();
    let C: Vec<Vec<Scalar>> = (0..n_circuits)
        .map(|_| {
            (0..input_length)
                .into_par_iter()
                .map(|_| Scalar::random())
                .collect()
        })
        .collect();
    (A, B, C)
}
