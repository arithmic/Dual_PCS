#![allow(unused_imports)]
use std::time::Instant;

use bls381::scalar::Scalar;
use fft::{
    par_fft::{get_inv_twiddles, par_interpolate_poly},
    serial_fft::interpolate_poly,
};
use traits::traits::{Field, PrimeField};

use crate::{
    commit::kzg_commit,
    polynomial::Polynomial,
    prover::kzg_prover,
    setup::kzg_setup,
    verifier::{self, kzg_verify},
};

#[test]
fn test_kzg() {
    let degree = 1 << 5;
    let start_time = Instant::now();
    let (g1_powers, verifier_key) = kzg_setup(degree);
    println!("Seetup time is {:?} ", start_time.elapsed());

    let evaluations: Vec<_> = (0..degree).map(|_| <Scalar as Field>::random()).collect();

    let start_time = Instant::now();
    let (commitment, coeff) = kzg_commit(evaluations.clone(), &g1_powers);
    println!("Commit time {:?}", start_time.elapsed());

    let polynomial = Polynomial { coeffs: coeff };
    let z = Scalar::random();

    let eval = polynomial.evaluate(&z);

    let start_time = Instant::now();
    let proof = kzg_prover(polynomial, &g1_powers, z);
    println!("Prover time {:?}", start_time.elapsed());

    let start_time = Instant::now();
    let v = verifier::kzg_verify(commitment, z, eval, proof, &verifier_key);
    println!("Verify time {:?}", start_time.elapsed());
    assert_eq!(v, true);
}
