#![allow(unused_imports)]
use bls381::{
    fp::{self, Fp},
    fp12::Fp12,
    scalar::Scalar,
};
use bls_curve::{
    bls::BlsCurve,
    fp2bls::{G2AffinePoint, G2ProjectivePoint},
};
use curve_traits::{
    self, projectivepoint::Scalar as BNSCALAR, CurveArithmetic, CurveParams, ProjectivePoint,
};
use pairing::{self, pairing_traits::Pairing};
use rayon::{
    iter::IntoParallelIterator,
    prelude::{IntoParallelRefIterator, ParallelIterator},
};
use std::{
    mem,
    time::{Duration, Instant},
};
use traits::traits::Field;

use crate::{
    common::{compute_fourier_bases, setup, MleCommit, MleEvalProof},
    prover::{batch_eval, commit, evaluate},
    verifier::{batch_verify, batch_verify_var_openings, verify},
};

#[test]
pub fn kzg_test() {
    let toxic_waste: Vec<_> = (0..10).map(|_| <Scalar as Field>::random()).collect();

    let start_time = Instant::now();
    let (srs, verifier_key) = setup::<BlsCurve>(toxic_waste);
    println!("Setup time {:?}", start_time.elapsed());

    let srs_size = srs.to_bytes();
    let verification_key_size = verifier_key.to_bytes();
    println!("Setup size {:?} KB", srs_size + verification_key_size);

    let poly: Vec<Scalar> = (0..1 << 10).map(|_| <Scalar as Field>::random()).collect();
    let point: Vec<_> = (0..10).map(|_| <Scalar as Field>::random()).collect();

    let start_time = Instant::now();
    let commitment: MleCommit<BlsCurve> = commit(&poly, &srs);
    println!("Commit time: {:?} ", start_time.elapsed());

    let start_time = Instant::now();
    let proof: MleEvalProof<BlsCurve> = evaluate(&poly, &point, &srs);
    println!("Proof time: {:?} ", start_time.elapsed());

    let proof_size = proof.to_bytes();
    println!("proof size {:?} KB", proof_size);

    let start_time = Instant::now();
    verify(&commitment, &proof, &point, &verifier_key);
    println!("Verify time: {:?} ", start_time.elapsed());

    println!("------------------------------------------------")
}

#[test]
fn batch_proof_test() {
    println!("Starting set up \n \n");
    let n_polys = 10;
    let toxic_waste: Vec<_> = (0..10).map(|_| <Scalar as Field>::random()).collect();

    let time = Instant::now();
    let (srs, verifier_key) = setup::<BlsCurve>(toxic_waste);
    println!("Setup took time {:?} \n \n", time.elapsed());
    println!("Vector size =  2^{:?} \n \n", 4);

    let polys: Vec<Vec<Scalar>> = (0..n_polys)
        .into_par_iter()
        .map(|_| (0..1 << 4).map(|_| <Scalar as Field>::random()).collect())
        .collect();
    let polys: Vec<&Vec<Scalar>> = polys.iter().map(|poly| poly).collect();
    let point: Vec<_> = (0..4).map(|_| <Scalar as Field>::random()).collect();
    let bases_eval: Vec<Scalar> = compute_fourier_bases::<BlsCurve>(&point);
    let evals: Vec<Scalar> = (0..n_polys)
        .into_par_iter()
        .map(|j| {
            bases_eval
                .iter()
                .zip(polys[j].iter())
                .fold(<Scalar as Field>::ZERO, |acc, (basis, coefficient)| {
                    acc + (*basis * *coefficient)
                })
        })
        .collect();

    println!("Starting Commit \n");

    let time_3 = Instant::now();
    let commitments: Vec<MleCommit<BlsCurve>> =
        polys.iter().map(|poly| commit(&poly, &srs)).collect();
    println!("Commit time: {:?} \n", time_3.elapsed());

    let random_scalars: Vec<Scalar> = (0..n_polys).map(|_| Scalar::random()).collect();
    println!(
        " Commitment size is: {:?} \n \n",
        mem::size_of_val(&commitments[0])
    );

    let time_4 = Instant::now();
    let proof: MleEvalProof<BlsCurve> = batch_eval(&polys, &evals, &point, &random_scalars, &srs);
    println!("Proof time: {:?} \n", time_4.elapsed());
    println!(
        " proof size is: {:?} \n \n",
        mem::size_of_val(&*proof.witnesses)
    );

    let time_5 = Instant::now();
    batch_verify(
        &commitments,
        &evals,
        &random_scalars,
        &proof,
        &point,
        &verifier_key,
    );
    println!("Verification time: {:?} \n", time_5.elapsed());

    println!("------------------------------------------------\n\n")
}

#[test]

fn var_opening_test() {
    let n_polys = 5;
    let toxic_waste: Vec<_> = (0..10).map(|_| <Scalar as Field>::random()).collect();

    let time = Instant::now();
    let (srs, verifier_key) = setup::<BlsCurve>(toxic_waste);
    println!("Setup took time {:?} \n \n", time.elapsed());
    println!("Vector size =  2^{:?} \n \n", 4);

    let polys: Vec<Vec<Scalar>> = (0..n_polys)
        .into_par_iter()
        .map(|_| (0..1 << 4).map(|_| <Scalar as Field>::random()).collect())
        .collect();

    let polys: Vec<&Vec<Scalar>> = polys.iter().map(|poly| poly).collect();

    let points: Vec<_> = (0..n_polys)
        .into_par_iter()
        .map(|_| (0..4).map(|_| <Scalar as Field>::random()).collect())
        .collect();

    let evals: Vec<Scalar> = (0..n_polys)
        .into_par_iter()
        .map(|i| {
            let basis_eval = compute_fourier_bases::<BlsCurve>(&points[i]);
            basis_eval
                .iter()
                .zip(polys[i].iter())
                .fold(<Scalar as Field>::ZERO, |acc, (basis, coefficient)| {
                    acc + (*basis * *coefficient)
                })
        })
        .collect();

    println!("Starting Commit \n");

    let time_3 = Instant::now();
    let commitments: Vec<MleCommit<BlsCurve>> =
        polys.iter().map(|poly| commit(&poly, &srs)).collect();
    println!("Commit time: {:?} \n", time_3.elapsed());

    let random_scalars: Vec<Scalar> = (0..n_polys).map(|_| Scalar::ONE).collect();
    println!(
        " Commitment size is: {:?} \n \n",
        mem::size_of_val(&commitments[0])
    );

    let time_4 = Instant::now();
    let mut proofs: Vec<MleEvalProof<BlsCurve>> = Vec::new();
    for i in 0..n_polys {
        proofs.push(evaluate(polys[i], &points[i], &srs));
    }
    println!("Proof time: {:?} \n", time_4.elapsed());
    println!(
        " proof size is: {:?} \n \n",
        mem::size_of_val(&*proofs[0].witnesses)
    );

    let time_5 = Instant::now();
    batch_verify_var_openings(
        &commitments,
        &evals,
        &proofs,
        &points,
        &random_scalars,
        &verifier_key,
    );
    println!("Verification time: {:?} \n", time_5.elapsed());

    println!("------------------------------------------------\n\n")
}
