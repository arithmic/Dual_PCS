use bls381::scalar::Scalar;
use helper::{Setup, WitnessCommitment};
use linking_prover::{commit::commit_witness, prover::linking_prover};
use linking_verifier::verifier::linking_verifier;
use setup::setup;
use std::time::Instant;
use traits::traits::Field;
use uni_multi_prover::evaluation_prover;
use uni_multi_verifier::evaluation_verifier;
extern crate bls381;
extern crate helper;
extern crate linking_prover;
extern crate linking_verifier;
extern crate setup;
extern crate traits;
fn main() {
    let degree_bound_setup = setup(1 << 8);
    let setup_size = degree_bound_setup.to_bytes();
    println!("Setup size {:?}KB", setup_size);

    let degree = 1 << 5;
    let setup = degree_bound_setup.setup_for_specific_degree(degree);
    let witness = (0..degree)
        .map(|_| Scalar::random())
        .collect::<Vec<Scalar>>();
    let witness_commitments = commit_witness(&witness, &setup);
    println!("----------");

    dory_evaluation(witness_commitments.clone(), witness, &setup);
    println!("--------------------------------------");
    linking_proof_system(witness_commitments, setup);
    println!("--------------------------------------");
    println!("--------------------------------------");
}

fn linking_proof_system(witness_commitments: WitnessCommitment, setup: Setup) {
    let start_time = Instant::now();
    let linking_proof = linking_prover(witness_commitments, setup.clone());
    println!("Prover time {:?},", start_time.elapsed());

    let proof_size = linking_proof.to_bytes();
    println!("Proof size {:?}KB,", proof_size);

    let start_time = Instant::now();
    linking_verifier(linking_proof, setup);
    println!("Verify time {:?},", start_time.elapsed());
    println!("------------------------")
}

fn dory_evaluation(witness_commitment: WitnessCommitment, witness: Vec<Scalar>, setup: &Setup) {
    let evaluation_proof = evaluation_prover(witness, witness_commitment, setup);

    let univariate_proof_size = evaluation_proof.0.to_bytes();
    println!("Univariate proof size {:?}KB,", univariate_proof_size);

    let multivariate_proof_size = evaluation_proof.1.to_bytes();
    println!("Multilinear proof size {:?}KB,", multivariate_proof_size);

    evaluation_verifier(evaluation_proof, setup);
    println!("------------------------")
}
