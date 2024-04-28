#![allow(unused_imports)]
use crate::batch_prover::batch_eval_multiple_evaluation_points;
use batch_verifier::batch_verifier;
use bls381::fp;
use bls381::scalar::Scalar;
use bls_curve::bls::{AffinePoint, BlsCurve};
use channel::Channel;
use kzg_fft::commit;
use kzg_fft::commit::kzg2commit;
use prover::prover;
use prover::{batch_eval, evaluate_multilinear_poly};
use rayon::result;
use setup::multilinearkzg2setup;
use std::io::Write;
use std::time::Instant;
use traits::traits::Field;
use verifier::batch_verify;
use verifier::verifier;

#[test]
fn test_multilinear_commitment_scheme() {
    let degree_bound = 1 << 10;
    let setup = multilinearkzg2setup(degree_bound);

    let degree = 1 << 9;

    let poly = (0..degree)
        .map(|_| Scalar::random())
        .collect::<Vec<Scalar>>();

    //Commit polynomial
    let commit_f = commit::kzg2commit(&poly, &setup.setup.get_setup(poly.len()).prover_key);

    //Initialize channel
    let mut channel = Channel::initialize_with_affine_point(&[commit_f].as_ref());
    let random_points = (0..9).map(|_| Scalar::random()).collect::<Vec<Scalar>>();

    let start_time = Instant::now();
    let proof = prover(&poly, random_points.clone(), &setup, &mut channel);
    println!("Proof time: {:?} ", start_time.elapsed());

    let start_time = Instant::now();
    let mut channel = Channel::initialize_with_affine_point(&[commit_f].as_ref());
    verifier(
        proof,
        random_points,
        setup.setup.get_setup(degree).verifier_key,
        setup.get_setup(degree).phi_commitment.clone(),
        &commit_f,
        &mut channel,
    )
    .unwrap();

    println!("Verifier time {:?} ", start_time.elapsed());
}

#[test]
fn test_batched_version() {
    let n_poly = 3;
    let degree_bound = 1 << 10;
    let setup = multilinearkzg2setup(degree_bound);
    for idx in 6..degree_bound.trailing_zeros() as usize {
        let degree = 1 << idx;
        let polys = (0..n_poly)
            .map(|_| {
                (0..degree)
                    .map(|_| Scalar::random())
                    .collect::<Vec<Scalar>>()
            })
            .collect::<Vec<Vec<Scalar>>>();
        let prover_key = &setup.setup.get_setup(degree).prover_key;
        let polys_commit = (0..n_poly)
            .map(|idx| kzg2commit(&polys[idx], prover_key))
            .collect::<Vec<AffinePoint>>();
        let mut channel = Channel::initialize_with_affine_point(&polys_commit);
        let random_points = (0..idx).map(|_| Scalar::random()).collect::<Vec<Scalar>>();
        let evals = (0..n_poly)
            .map(|idx| evaluate_multilinear_poly(&polys[idx], &random_points))
            .collect::<Vec<Scalar>>();
        let combiners = (0..n_poly)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let proof = batch_eval(
            &polys,
            &combiners,
            evals.clone(),
            random_points.clone(),
            &setup,
            &mut channel,
        );
        let mut channel = Channel::initialize_with_affine_point(&polys_commit);
        let _ = batch_verify(
            proof,
            polys_commit,
            evals,
            combiners,
            random_points,
            setup.setup.get_setup(degree).verifier_key,
            setup.get_setup(degree).phi_commitment.clone(),
            &mut channel,
        )
        .unwrap();
    }
}

#[test]
fn test_batched_version2() {
    let setup = multilinearkzg2setup(1 << 8);

    let degree = 1 << 5;

    let poly = (0..2)
        .map(|_| {
            (0..degree)
                .map(|_| Scalar::random())
                .collect::<Vec<Scalar>>()
        })
        .collect::<Vec<Vec<Scalar>>>();

    //Commit polynomial
    let commit_f = (0..2)
        .map(|iter| {
            commit::kzg2commit(
                &poly[iter],
                &setup.setup.get_setup(poly[iter].len()).prover_key,
            )
        })
        .collect::<Vec<AffinePoint>>();

    //Initialize channel
    let mut channel = Channel::initialize_with_affine_point(&commit_f);
    let random_points = (0..2)
        .map(|_| (0..5).map(|_| Scalar::random()).collect::<Vec<Scalar>>())
        .collect::<Vec<Vec<Scalar>>>();
    let evals = (0..2)
        .map(|iter| evaluate_multilinear_poly(&poly[iter], &random_points[iter]))
        .collect::<Vec<Scalar>>();
    let start_time = Instant::now();
    let proof = batch_eval_multiple_evaluation_points(
        &poly,
        evals,
        random_points.clone(),
        &setup,
        &mut channel,
    );
    println!("time {:?}", start_time.elapsed());
    let mut channel = Channel::initialize_with_affine_point(&commit_f);

    batch_verifier(
        proof,
        random_points,
        setup.setup.get_setup(degree).verifier_key,
        setup.get_setup(degree).phi_commitment.clone(),
        &commit_f,
        &mut channel,
    )
    .unwrap();
}
