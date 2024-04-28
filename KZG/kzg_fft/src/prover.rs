use crate::commit;
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use common::KZGFFTEvaluationProof;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::{Field, PrimeField};

pub fn kzg2_prover(
    evaluations: &Vec<Scalar>,
    z: Scalar,
    prover_key: &Vec<curve_traits::ProjectivePoint<BlsCurve>>,
) -> KZGFFTEvaluationProof {
    let degree = evaluations.len();
    let log2_degree = degree.trailing_zeros();
    let omega = Scalar::get_root_of_unity(log2_degree);
    let mut omega_powers = vec![Scalar::ONE; degree];
    for idx in 1..degree {
        let temp = omega_powers[idx - 1] * omega;
        omega_powers[idx] = temp;
    }
    //Compute alpha_i
    let degree_inverse = (Scalar::from(degree as u64)).invert().unwrap();
    let beta = (0..degree)
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

    let evaluation_at_z =
        (0..degree).fold(Scalar::ZERO, |acc, idx| acc + evaluations[idx] * beta[idx]);

    let evaluation_of_q = (0..degree)
        .into_par_iter()
        .map(|idx| (evaluations[idx] - evaluation_at_z) * (omega_powers[idx] - z).invert().unwrap())
        .collect::<Vec<Scalar>>();

    let commitment_of_evaluations_of_q = commit::kzg2commit(&evaluation_of_q, prover_key);
    KZGFFTEvaluationProof {
        commitment_of_evaluations_of_q,
        evaluation_at_z,
    }
}

pub fn kzg2_batch_prover(
    evaluations: &Vec<Vec<Scalar>>,
    random_point: Scalar,
    evals: &Vec<Scalar>,
    combiners: &Vec<Scalar>,
    prover_key: &Vec<curve_traits::ProjectivePoint<BlsCurve>>,
) -> KZGFFTEvaluationProof {
    let n_evals = evaluations.len();
    let combined_evaluations = (0..evaluations[0].len())
        .map(|outer_idx| {
            (0..n_evals).fold(Scalar::ZERO, |acc, inner_idx| {
                acc + (combiners[inner_idx] * evaluations[inner_idx][outer_idx])
            })
        })
        .collect::<Vec<Scalar>>();

    let combined_eval =
        (0..evals.len()).fold(Scalar::ZERO, |acc, idx| acc + (combiners[idx] * evals[idx]));
    let batched_proof = kzg2_prover(&combined_evaluations, random_point, prover_key);
    assert_eq!(batched_proof.evaluation_at_z, combined_eval);
    batched_proof
}
