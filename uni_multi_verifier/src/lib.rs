use std::time::Instant;

use helper::{EvaluationProof, Setup};
use multilinear::multilinear_evaluation_verifier;
use univariate::univariate_evaluation_verifier;

pub mod multilinear;
pub mod univariate;
extern crate bls381;
extern crate bls_curve;
extern crate channel;
extern crate crossbeam;
extern crate curve_traits;
extern crate dory;
extern crate fft;
extern crate helper;
extern crate pairing;
extern crate rayon;
extern crate setup;
extern crate traits;
#[allow(unused)]
pub fn evaluation_verifier(evaluation_proof: (EvaluationProof, EvaluationProof), setup: &Setup) {
    let start_time = Instant::now();
    let univariate_evaluation = univariate_evaluation_verifier(evaluation_proof.0, setup).unwrap();
    println!("Univariate verifier time {:?},", start_time.elapsed());

    let start_time = Instant::now();
    let multivariate_evaluation =
        multilinear_evaluation_verifier(evaluation_proof.1, setup).unwrap();
    println!("Multivariate verifier time {:?},", start_time.elapsed());
}
