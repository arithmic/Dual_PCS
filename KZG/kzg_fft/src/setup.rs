use bls381::scalar::Scalar;
use bls_curve::{bls::BlsCurve, fp2bls::G2BlsCurve};
use curve_traits::{AffinePoint, ProjectivePoint};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::{Field, PrimeField};

use crate::common::KZGFFTDegreeBoundSetup;

pub fn kzg2_setup(degree_bound: usize) -> KZGFFTDegreeBoundSetup {
    let log2_degree_bound = degree_bound.trailing_zeros();
    //Generator of G1
    let g1 = <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;
    //Generator of G2
    let g2 = <G2BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;
    let omega = Scalar::get_root_of_unity(log2_degree_bound);
    let mut omega_powers = vec![Scalar::ONE; degree_bound];
    for idx in 1..degree_bound {
        omega_powers[idx] = omega_powers[idx - 1] * omega;
    }
    let r = Scalar::random();

    //Compute alpha_i
    let degree_bound_inverse = Scalar::from(degree_bound as u64).invert().unwrap();
    let alpha = (0..degree_bound)
        .into_par_iter()
        .map(|i| {
            let mut temp = Scalar::ONE;
            for j in 0..log2_degree_bound {
                temp *= Scalar::ONE
                    + (omega_powers[(degree_bound - i) % degree_bound] * r).power_by([
                        1_u64 << j,
                        0,
                        0,
                        0,
                    ])
            }
            temp * degree_bound_inverse
        })
        .collect::<Vec<Scalar>>();

    let mut prover_keys = Vec::new();

    let prover_key = (0..degree_bound)
        .into_par_iter()
        .map(|idx| g1.mul(alpha[idx]))
        .collect::<Vec<ProjectivePoint<BlsCurve>>>();

    prover_keys.push(prover_key.clone());

    let mut degree = degree_bound / 2;
    let mut folded_prover_key = prover_key.clone();
    while degree > 0 {
        for j in 0..folded_prover_key.len() / 2 {
            folded_prover_key[j] =
                folded_prover_key[j] + folded_prover_key[j + folded_prover_key.len() / 2];
        }
        folded_prover_key.truncate(folded_prover_key.len() / 2);
        prover_keys.push(folded_prover_key.clone());
        degree /= 2;
    }

    let mut verifier_key =
        vec![AffinePoint::<G2BlsCurve>::IDENTITY; (log2_degree_bound + 1) as usize];
    let mut r_square = r;

    for v_k in verifier_key
        .iter_mut()
        .take((log2_degree_bound + 1) as usize)
    {
        *v_k = g2.mul(r_square).to_affine();
        r_square = r_square.square();
    }
    KZGFFTDegreeBoundSetup {
        prover_keys,
        verifier_key,
    }
}
