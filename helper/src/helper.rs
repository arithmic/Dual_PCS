use bls381::scalar::Scalar;
use bls_curve::{bls::BlsCurve, fp2bls::G2BlsCurve, gt::Gt};
use pairing::pairing_traits::Pairing;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator, ParallelIterator,
};

/// Compute the Inner Pairing Product of two vectors IPP = /summation_i{i, n} e(w_i, g_i).
/// Tate pairing is used to compute pairing of two elements.
///
/// The size of the both vectors is assumed to be equal.
///
///Crossbeam and rayon crates are used for multithreading.
///
pub fn compute_ipp(
    left_input: &Vec<curve_traits::ProjectivePoint<BlsCurve>>,
    right_input: &Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
) -> Gt {
    assert_eq!(
        left_input.len(),
        right_input.len(),
        "Size of both inputs should be same for pairing"
    );

    let pairing_output = left_input
        .par_iter()
        .zip(right_input.par_iter())
        .fold_with(Gt::ONE, |acc, (value1, value2)| {
            acc + &(BlsCurve::tate_pairing(value1.to_affine(), value2.to_affine()))
        })
        .reduce(|| Gt::ONE, |acc, part_sum| acc + &part_sum);

    pairing_output
}

pub fn compute_dot_product(left_input: &Vec<Scalar>, right_input: &Vec<Scalar>) -> Scalar {
    assert_eq!(
        left_input.len(),
        right_input.len(),
        "Size of both inputs should be same for pairing"
    );
    let dot_product = left_input
        .par_iter()
        .zip_eq(right_input.par_iter())
        .fold_with(Scalar::ZERO, |acc, (value1, value2)| {
            acc + (*value1 * *value2)
        })
        .reduce(|| Scalar::ZERO, |acc, part_sum| acc + part_sum);
    dot_product
}

pub fn compute_y(r: &Vec<Scalar>) -> Vec<Scalar> {
    //Initialize fc_eq with (1- r[0]) and r[0]
    let mut y = [Scalar::ONE - r[r.len() - 1], r[r.len() - 1]].to_vec();
    //Iterate over the length of the r vector
    for k in (0..r.len() - 1).rev() {
        let temp = y;
        //initialize fc_eq of double size with zero
        y = vec![Scalar::ZERO; temp.len() * 2];
        for iter in 0..temp.len() {
            y[2 * iter] = temp[iter] * (Scalar::ONE - r[k]);
            y[(2 * iter) + 1] = temp[iter] * (r[k]);
        }
    }
    y
}

pub fn compute_fourier_bases(r: &Vec<Scalar>) -> Vec<Scalar> {
    let mut bases: Vec<Scalar> = vec![Scalar::ZERO; 1 << r.len()];
    let mut size = 1;
    bases[0] = Scalar::ONE;

    for r in r.iter().rev() {
        let (bases_left, bases_right) = bases.split_at_mut(size);
        let (bases_right, _) = bases_right.split_at_mut(size);

        bases_left
            .par_iter_mut()
            .zip_eq(bases_right.par_iter_mut())
            .for_each(|(b_l, b_r)| {
                *b_r = *b_l * *r;
                *b_l -= *b_r;
            });

        size *= 2;
    }

    bases
}
