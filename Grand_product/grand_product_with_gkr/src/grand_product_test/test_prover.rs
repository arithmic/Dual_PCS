#![allow(non_snake_case)]
#![allow(unused)]
use crate::gkr_common::{Commitments, Evaluations, GkrTranscript};
use crate::prover::gkr_prover;
use crate::prover::reed_solomon::reed_solomon;
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use channel::Channel;
use multilinear_kzg::common::{MleCommit, MleEvalProof, SRS};
use multilinear_kzg::prover::{batch_eval, commit};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use super::common::reseed_commits;
pub fn prover(
    A: &Vec<Vec<Scalar>>,
    B: &Vec<Vec<Scalar>>,
    C: &Vec<Vec<Scalar>>,
    n_circuits: usize,
    srs: &SRS<BlsCurve>,
) -> (
    Commitments,
    GkrTranscript,
    MleEvalProof<BlsCurve>,
    Evaluations,
) {
    let commitments = commits(A, B, C, srs);

    let mut channel =
        Channel::initialize_with_affine_point(&[commitments.A_commits[0].commitment.to_affine()]);
    reseed_commits(commitments.clone(), &mut channel);

    let circuits = reed_solomon(&A, &B, &C, n_circuits, &mut channel);
    let transcript = gkr_prover(&circuits, &mut channel);

    let point = transcript.final_layer_point.clone();
    let (evaluation_proof, evaluations) =
        evaluations(A, B, C, n_circuits, &point, srs, &mut channel);
    (commitments, transcript, evaluation_proof, evaluations)
}

fn commits(
    A: &Vec<Vec<Scalar>>,
    B: &Vec<Vec<Scalar>>,
    C: &Vec<Vec<Scalar>>,
    srs: &SRS<BlsCurve>,
) -> Commitments {
    let A_commis = (0..A.len())
        .map(|idx| commit(&A[idx], &srs))
        .collect::<Vec<MleCommit<BlsCurve>>>();
    let B_commits = (0..B.len())
        .map(|idx| commit(&B[idx], &srs))
        .collect::<Vec<MleCommit<BlsCurve>>>();
    let C_commits = (0..C.len())
        .map(|idx| commit(&C[idx], &srs))
        .collect::<Vec<MleCommit<BlsCurve>>>();

    Commitments::new(A_commis, B_commits, C_commits)
}

fn evaluations(
    A: &Vec<Vec<Scalar>>,
    B: &Vec<Vec<Scalar>>,
    C: &Vec<Vec<Scalar>>,
    n_circuits: usize,
    point: &Vec<Scalar>,
    srs: &SRS<BlsCurve>,
    channel: &mut Channel,
) -> (MleEvalProof<BlsCurve>, Evaluations) {
    let bases = compute_coeff(point);
    let A_evals: Vec<Scalar> = A
        .into_iter()
        .map(|a| {
            a.par_iter()
                .zip(bases.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();
    let B_evals: Vec<Scalar> = B
        .into_iter()
        .map(|a| {
            a.par_iter()
                .zip(bases.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();
    let C_evals: Vec<Scalar> = C
        .into_iter()
        .map(|a| {
            a.par_iter()
                .zip(bases.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();
    let scalars = channel.get_random_points(3 * n_circuits);
    let A: Vec<&Vec<Scalar>> = A.iter().map(|A| A).collect();
    let B: Vec<&Vec<Scalar>> = B.iter().map(|B| B).collect();
    let C: Vec<&Vec<Scalar>> = C.iter().map(|C| C).collect();
    let proof = batch_eval(
        &[A, B, C].concat(),
        &[A_evals.clone(), B_evals.clone(), C_evals.clone()].concat(),
        &point,
        &scalars,
        srs,
    );
    (proof, Evaluations::new(A_evals, B_evals, C_evals))
}
pub fn compute_coeff(r: &Vec<Scalar>) -> Vec<Scalar> {
    //Initialize fc_eq with (1- r[0]) and r[0]
    let mut fc = [Scalar::ONE - r[0], r[0]].to_vec();
    //Iterate over the length of the r vector
    for k in 1..r.len() {
        let temp = fc;
        fc = vec![Scalar::ZERO; temp.len() * 2];
        for iter in 0..temp.len() {
            fc[2 * iter] = temp[iter] * (Scalar::ONE - r[k as usize]);
            fc[2 * iter + 1] = temp[iter] * r[k as usize];
        }
    }
    fc
}
