#![allow(unused_imports)]
#![allow(non_snake_case)]
use crate::gkr_common::mle_eval_proof_to_bytes;
use crate::grand_product_test::test_prover::prover;
use crate::grand_product_test::test_verifier::verifier;
use crate::prover::{compute_fourier_bases, gkr_prover, reed_solomon};
use bls381::scalar::Scalar;
use bls_curve::bls::AffinePoint;
use bls_curve::bls::BlsCurve;
use channel::Channel;
use fft::serial_fft::log2;
use multilinear_kzg::common::{setup, SRS};
use rand::{thread_rng, Rng};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use reed_solomon::reed_solomon;
use std::time::Instant;
use traits::traits::Field;
pub mod common;
pub mod test_prover;
pub mod test_verifier;
#[test]
fn test_gkr() {
    let total_circuits = 1;
    let input_length = 1 << 10;
    let toxic_waste = (0..10).into_par_iter().map(|_| Scalar::random()).collect();
    let (srs, ver_key) = setup(toxic_waste);

    let (A, B, C) = construct_input(total_circuits, input_length);

    let start_time = Instant::now();
    let (commitments, gkr_transcript, evaluation_proof, evaluations) =
        prover(&A, &B, &C, total_circuits, &srs);
    println!("Prover time {:?}", start_time.elapsed());

    let commitment_size = commitments.to_bytes();
    println!("Commitment size {:?}", commitment_size);

    let evaluation_proof_size =
        mle_eval_proof_to_bytes(&[evaluation_proof.clone()].to_vec()).len() as f64 / 1024f64;
    let proof_size = (gkr_transcript.to_bytes().len() as f64 / 1024f64)
        + evaluation_proof_size
        + evaluations.to_bytes();
    println!("Proof size {:?}", proof_size);

    let start_time = Instant::now();
    verifier(
        gkr_transcript,
        commitments,
        evaluation_proof,
        evaluations,
        ver_key,
        total_circuits,
    );
    println!("Verify time {:?}", start_time.elapsed());

    println!("----------------------------")
}
#[allow(dead_code)]
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
