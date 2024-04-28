#![allow(unused)]
#![allow(non_snake_case)]
use super::common::reseed_commits;
use crate::{
    gkr_common::{Commitments, Evaluations, GkrTranscript},
    verifier::gkr_verifier,
};
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use channel::Channel;
use multilinear_kzg::{
    common::{MleEvalProof, VerificationKey},
    verifier::batch_verify,
};
use traits::traits::Field;

pub fn verifier(
    gkr_transcript1: GkrTranscript,
    commitments: Commitments,
    evaluation_proof: MleEvalProof<BlsCurve>,
    evaluations: Evaluations,
    ver_key: VerificationKey<BlsCurve>,
    n_circuits: usize,
) {
    let mut channel =
        Channel::initialize_with_affine_point(&[commitments.A_commits[0].commitment.to_affine()]);
    reseed_commits(commitments.clone(), &mut channel);
    let point1 = gkr_transcript1.clone().final_layer_point;

    let gamma_tau = channel.get_random_points(2);
    let mut circuit_evals = vec![Scalar::ZERO; n_circuits];
    let A_evals = evaluations.A_evals;
    let B_evals = evaluations.B_evals;
    let C_evals = evaluations.C_evals;
    for c in 0..n_circuits {
        circuit_evals[c] =
            gamma_tau[0].square() * A_evals[c] + gamma_tau[0] * B_evals[c] + C_evals[c]
                - gamma_tau[1];
    }

    gkr_verifier(
        &gkr_transcript1,
        point1.len(),
        &mut channel,
        circuit_evals,
        n_circuits,
    );
    let scalars = channel.get_random_points(3 * n_circuits);
    batch_verify(
        &[
            commitments.A_commits,
            commitments.B_commits,
            commitments.C_commits,
        ]
        .concat(),
        &[A_evals, B_evals, C_evals].concat(),
        &scalars,
        &evaluation_proof,
        &gkr_transcript1.final_layer_point,
        &ver_key,
    );
}
