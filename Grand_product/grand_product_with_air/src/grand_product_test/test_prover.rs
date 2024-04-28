use bls381::scalar::Scalar;
use bls_curve::bls::AffinePoint;
use channel::Channel;
use kzg_fft::{
    commit::kzg2commit,
    common::{KZGFFTDegreeBoundSetup, KZGFFTEvaluationProof},
    prover::kzg2_batch_prover,
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::{Field, PrimeField};

use crate::{
    grand_product_common::{Commitments, Evaluations, GkrTranscript},
    prover::{gkr_prover, grand_product::grand_product},
};

use super::common::{self, reseed_with_commits};

pub fn prover(
    A: &Vec<Vec<Scalar>>,
    B: &Vec<Vec<Scalar>>,
    C: &Vec<Vec<Scalar>>,
    n_circuits: usize,
    setup: &KZGFFTDegreeBoundSetup,
) -> (
    Commitments,
    GkrTranscript,
    Evaluations,
    KZGFFTEvaluationProof,
    KZGFFTEvaluationProof,
) {
    let commitments = commit(&A, &B, &C, setup);
    let mut channel = Channel::initialize_with_affine_point(commitments.get_A_commits());
    reseed_with_commits(&commitments, &mut channel);

    let (leaf_layers, final_layers) = grand_product(&A, &B, &C, n_circuits, &mut channel);

    let (gkr_transcript, evaluation_point, trace, constraint_evaluations) = gkr_prover(
        &leaf_layers,
        final_layers,
        &mut channel,
        &setup.get_setup(leaf_layers[0].len()).prover_key,
    );

    let gkr_evaluations = evaluations(&A, &B, &C, evaluation_point);
    let mut next_polys = Vec::new();
    (0..A.len()).for_each(|idx| next_polys.push(A[idx].clone()));
    (0..B.len()).for_each(|idx| next_polys.push(B[idx].clone()));
    (0..C.len()).for_each(|idx| next_polys.push(C[idx].clone()));

    (0..trace.num_of_columns()).for_each(|idx| {
        if idx % 2 != 0 {
            next_polys.push(trace.column_at(idx))
        }
    });
    let mut current_polys = next_polys.clone();
    current_polys.push(constraint_evaluations);

    let current_frame = gkr_transcript.get_odd_frame().get_current_frame();
    let next_frame = gkr_transcript.get_odd_frame().get_next_frame();
    let compositon_frame = gkr_transcript.get_odd_frame().get_composition_frame();

    let mut current_evals = Vec::new();
    gkr_evaluations
        .A_at_z
        .iter()
        .for_each(|value| current_evals.push(*value));
    gkr_evaluations
        .B_at_z
        .iter()
        .for_each(|value| current_evals.push(*value));
    gkr_evaluations
        .C_at_z
        .iter()
        .for_each(|value| current_evals.push(*value));
    current_frame
        .iter()
        .for_each(|value| current_evals.push(*value));
    current_evals.push(compositon_frame);

    let mut next_evals = Vec::new();
    gkr_evaluations
        .A_at_gz
        .iter()
        .for_each(|value| next_evals.push(*value));
    gkr_evaluations
        .B_at_gz
        .iter()
        .for_each(|value| next_evals.push(*value));
    gkr_evaluations
        .C_at_gz
        .iter()
        .for_each(|value| next_evals.push(*value));
    next_frame.iter().for_each(|value| next_evals.push(*value));

    let prover_key = setup.get_setup(A[0].len()).prover_key;
    let combiners = channel.get_random_points(current_evals.len());

    let proof1 = kzg2_batch_prover(
        &current_polys,
        evaluation_point,
        &current_evals,
        &combiners,
        &prover_key,
    );
    let trace_length = A[0].len();
    let g_trace = Scalar::get_root_of_unity(trace_length.trailing_zeros());
    let proof2 = kzg2_batch_prover(
        &next_polys,
        g_trace * evaluation_point,
        &next_evals,
        &combiners[0..&combiners.len() - 1].to_vec(),
        &prover_key,
    );
    (commitments, gkr_transcript, gkr_evaluations, proof1, proof2)
}

#[allow(dead_code)]
pub fn commit(
    A: &Vec<Vec<Scalar>>,
    B: &Vec<Vec<Scalar>>,
    C: &Vec<Vec<Scalar>>,
    setup: &KZGFFTDegreeBoundSetup,
) -> Commitments {
    let prover_key = &setup.get_setup(A[0].len()).prover_key;
    let A_commit = (0..B.len())
        .map(|idx| kzg2commit(&A[idx], prover_key))
        .collect::<Vec<AffinePoint>>();
    let B_commits = (0..B.len())
        .map(|idx| kzg2commit(&B[idx], prover_key))
        .collect::<Vec<AffinePoint>>();
    let C_commits = (0..C.len())
        .map(|idx| kzg2commit(&C[idx], prover_key))
        .collect::<Vec<AffinePoint>>();
    Commitments::new(A_commit, B_commits, C_commits)
}
#[allow(dead_code)]
pub fn evaluations(
    A: &Vec<Vec<Scalar>>,
    B: &Vec<Vec<Scalar>>,
    C: &Vec<Vec<Scalar>>,
    z: Scalar,
) -> Evaluations {
    let trace_length1 = A[0].len();
    assert!(
        trace_length1.is_power_of_two(),
        " program code lenght should be power of 2"
    );

    let g_trace = Scalar::get_root_of_unity(trace_length1.trailing_zeros());
    let g_z = g_trace * z;

    let degree = trace_length1;
    let log2_degree = degree.trailing_zeros();
    let omega = g_trace;
    let mut omega_powers = vec![Scalar::ONE; degree];
    for idx in 1..degree {
        let temp = omega_powers[idx - 1] * omega;
        omega_powers[idx] = temp;
    }
    //Compute alpha_i
    let degree_inverse = (Scalar::from(degree as u64)).invert().unwrap();
    let beta_z = (0..degree)
        .into_par_iter()
        .map(|i| {
            let mut temp = Scalar::ONE;
            for j in 0..log2_degree {
                temp *= Scalar::ONE
                    + (omega_powers[(degree - i) % degree] * z).power_by([1_u64 << j, 0, 0, 0])
            }
            temp * degree_inverse
        })
        .collect::<Vec<Scalar>>();

    let beta_gz = (0..degree)
        .into_par_iter()
        .map(|i| {
            let mut temp = Scalar::ONE;
            for j in 0..log2_degree {
                temp *= Scalar::ONE
                    + (omega_powers[(degree - i) % degree] * g_z).power_by([1_u64 << j, 0, 0, 0])
            }
            temp * degree_inverse
        })
        .collect::<Vec<Scalar>>();

    let A_at_z = (0..A.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + A[outer_idx][idx] * beta_z[idx]
            })
        })
        .collect();
    let A_at_gz = (0..A.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + A[outer_idx][idx] * beta_gz[idx]
            })
        })
        .collect();

    let B_at_z = (0..B.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + B[outer_idx][idx] * beta_z[idx]
            })
        })
        .collect::<Vec<Scalar>>();
    let B_at_gz = (0..B.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + B[outer_idx][idx] * beta_gz[idx]
            })
        })
        .collect::<Vec<Scalar>>();
    let C_at_z = (0..C.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + C[outer_idx][idx] * beta_z[idx]
            })
        })
        .collect::<Vec<Scalar>>();
    let C_at_gz = (0..C.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + C[outer_idx][idx] * beta_gz[idx]
            })
        })
        .collect::<Vec<Scalar>>();
    Evaluations::new(A_at_z, A_at_gz, B_at_z, B_at_gz, C_at_z, C_at_gz)
}
