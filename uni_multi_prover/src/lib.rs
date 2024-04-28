use std::time::Instant;

use bls381::scalar::Scalar;
use helper::{EvaluationProof, Setup, WitnessCommitment};
use multilinear::multilinear_evaluation_prover;
use univariate::univariate_evaluation_prover;

pub mod multilinear;
pub mod univariate;
extern crate bls381;
extern crate bls_curve;
extern crate channel;
extern crate curve_traits;
extern crate dory;
extern crate fft;
extern crate helper;
extern crate pairing;
extern crate rayon;
extern crate setup;
extern crate traits;
///Execute Dory for Univariate and multilinear polynomials
pub fn evaluation_prover(
    witness: Vec<Scalar>,
    witness_commitments: WitnessCommitment,
    setup: &Setup,
) -> (EvaluationProof, EvaluationProof) {
    let start_time = Instant::now();
    let univariate_evaluation_proof = univariate_evaluation_prover(
        setup,
        setup.tau_1.clone(),
        [
            witness_commitments.commitment_to_univariate,
            setup.clone().tau_ipp[0],
            witness_commitments.commitment_to_univariate,
        ]
        .to_vec(),
        witness_commitments.g2_power_poly,
        witness_commitments.univariate_polynomial,
    );
    println!("Univariate Evaluation time {:?},", start_time.elapsed());

    let start_time = Instant::now();
    let multilinear_evaluation_proof = multilinear_evaluation_prover(
        &setup,
        setup.tau_1.clone(),
        [
            witness_commitments.commitment_to_witness,
            setup.tau_ipp[0],
            witness_commitments.commitment_to_witness,
        ]
        .to_vec(),
        witness_commitments.g2_power_witness,
        witness,
    );
    println!("Multilinear Evaluation time {:?},", start_time.elapsed());

    (univariate_evaluation_proof, multilinear_evaluation_proof)
}
