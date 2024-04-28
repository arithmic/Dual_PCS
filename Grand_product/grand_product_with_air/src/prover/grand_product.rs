use bls381::scalar::Scalar;
use channel::Channel;
use fft::serial_fft::log2;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use traits::traits::Field;
#[allow(non_snake_case)]
pub fn grand_product(
    A: &Vec<Vec<Scalar>>,
    B: &Vec<Vec<Scalar>>,
    C: &Vec<Vec<Scalar>>,
    n_circuits: usize,
    channel: &mut Channel,
) -> (Vec<Vec<Scalar>>, Vec<Scalar>) {
    let input_length = A[0].len();
    let depth = log2(input_length) as usize;
    let gamma_tau = channel.get_random_points(2);
    let mut leaf_layer_1 = vec![vec![Scalar::ONE; 1 << depth]; n_circuits];

    let mut final_layer_1 = Vec::new();

    for c in 0..n_circuits {
        for i in 0..input_length {
            leaf_layer_1[c][i] =
                gamma_tau[0].square() * A[c][i] + gamma_tau[0] * B[c][i] + C[c][i] - gamma_tau[1];
        }

        final_layer_1.push(
            leaf_layer_1[c]
                .par_iter()
                .fold(|| Scalar::ONE, |acc, value| acc * *value)
                .reduce(|| Scalar::ONE, |acc, x| acc * x),
        );
    }

    (leaf_layer_1, final_layer_1)
}
