use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use common::KZGFourierDegreeBoundSetup;
use curve_traits::AffinePoint;
use kzg_fft::{commit, setup};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::PrimeField;
pub fn multilinearkzg2setup(degree_bound: usize) -> KZGFourierDegreeBoundSetup {
    let setup = setup::kzg2_setup(degree_bound);
    let phi_commitments = (0..degree_bound.trailing_zeros())
        .map(|idx| {
            let degree = 1 << idx;
            let degree_setup = setup.get_setup(1 << (idx + 1));
            let evaluations = &compute_evaluations(degree);
            commit::kzg2commit(
                &evaluations[degree..].to_vec(),
                &degree_setup.prover_key[degree..].to_vec(),
            )
        })
        .collect::<Vec<AffinePoint<BlsCurve>>>();
    KZGFourierDegreeBoundSetup {
        setup,
        phi_commitments,
        degree_bound,
    }
}
pub fn compute_evaluations(degree: usize) -> Vec<Scalar> {
    let w_2_k_plus_one = Scalar::get_root_of_unity(degree.trailing_zeros() + 1);
    let mut w_2_k_plus_one_powers = vec![Scalar::ONE; 2 * degree];

    for idx in 1..w_2_k_plus_one_powers.len() {
        w_2_k_plus_one_powers[idx] = w_2_k_plus_one_powers[idx - 1] * w_2_k_plus_one
    }
    [
        vec![Scalar::ZERO; degree],
        (degree..2 * degree)
            .into_par_iter()
            .map(|outer_idx| {
                let mut eval = Scalar::ONE;
                for inner_idx in 0..degree {
                    eval *= w_2_k_plus_one_powers[outer_idx] - w_2_k_plus_one_powers[inner_idx]
                }
                eval
            })
            .collect::<Vec<Scalar>>(),
    ]
    .concat()
}
