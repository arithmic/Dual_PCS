#![allow(unused)]
use crate::{prover::kzg2_batch_prover, verifier::kzg2_batch_verify};
use bls381::scalar::Scalar;
use bls_curve::bls::AffinePoint;
use channel::Channel;
use commit::kzg2commit;
use fft::{
    par_fft::{get_inv_twiddles, par_interpolate_poly},
    serial_fft::eval,
};
use prover::kzg2_prover;
use setup::kzg2_setup;
use std::time::Instant;
use traits::traits::Field;
use verifier::kzg2_verifier;
#[test]
fn test_uni_poly() {
    let degree_bound: usize = 1 << 8;
    let degree_bound_setup = kzg2_setup(degree_bound);
    print!(
        "Setup size is {:?}KB",
        degree_bound_setup.to_bytes().len() as f64 / 1024f64
    );

    let degree = 1 << 7;

    let evaluations = (0..degree)
        .map(|_| Scalar::random())
        .collect::<Vec<Scalar>>();
    let setup = degree_bound_setup.get_setup(degree);
    let start_time = Instant::now();
    let commitment_to_evaluations_of_f = kzg2commit(&evaluations, &setup.prover_key);

    let mut channel =
        Channel::initialize_with_affine_point(&[commitment_to_evaluations_of_f].as_ref());
    let z = channel.get_random_point();

    let evaluation_proof = kzg2_prover(&evaluations, z, &setup.prover_key);

    let mut channel =
        Channel::initialize_with_affine_point(&[commitment_to_evaluations_of_f].as_ref());
    let z = channel.get_random_point();
    let start_time = Instant::now();
    let verify = kzg2_verifier(
        setup.verifier_key,
        evaluation_proof,
        commitment_to_evaluations_of_f,
        z,
    )
    .expect("verification failed");
}

#[test]
fn test_batch() {
    let degree_bound_setup = kzg2_setup(1 << 10);
    let degree = 1 << 6;
    let setup = degree_bound_setup.get_setup(degree);
    let mut evaluations = (0..2)
        .map(|idx| {
            (0..degree)
                .map(|_| Scalar::random())
                .collect::<Vec<Scalar>>()
        })
        .collect::<Vec<Vec<Scalar>>>();
    let inv_twiddles = get_inv_twiddles(degree as u32);

    let start_time = Instant::now();
    let commitment_to_evaluations_of_f = (0..2)
        .map(|idx| kzg2commit(&evaluations[idx], &setup.prover_key))
        .collect::<Vec<AffinePoint>>();

    let mut channel = Channel::initialize_with_affine_point(&commitment_to_evaluations_of_f);
    let z = channel.get_random_point();
    let mut polys = evaluations.clone();
    let evals = (0..2)
        .map(|idx| {
            par_interpolate_poly(&mut polys[idx], inv_twiddles.clone());
            eval(&polys[idx], z)
        })
        .collect::<Vec<Scalar>>();
    let combiners = [Scalar::random(), Scalar::random()].to_vec();
    let evaluation_proof =
        kzg2_batch_prover(&evaluations, z, &evals, &combiners, &setup.prover_key);

    let mut channel =
        Channel::initialize_with_affine_point(commitment_to_evaluations_of_f.as_ref());
    let z = channel.get_random_point();
    let start_time = Instant::now();
    kzg2_batch_verify(
        evaluation_proof,
        setup.verifier_key,
        z,
        commitment_to_evaluations_of_f,
        evals,
        combiners,
    );
}
