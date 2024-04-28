#![allow(non_snake_case)]
use bls381::scalar::Scalar;
use channel::Channel;
use fft::serial_fft::log2;
use gkr_common::{CircuitBinaryTree, CircuitLayer};
use traits::traits::Field;

use crate::gkr_common;

pub fn reed_solomon(
    A: &Vec<Vec<Scalar>>,
    B: &Vec<Vec<Scalar>>,
    C: &Vec<Vec<Scalar>>,
    n_circuits: usize,
    channel: &mut Channel,
) -> Vec<CircuitBinaryTree> {
    let input_length = A[0].len();
    let depth_1 = log2(input_length) as usize;

    let mut circuits: Vec<CircuitBinaryTree> = Vec::new();

    let gamma_tau = channel.get_random_points(2);

    for c in 0..n_circuits {
        let mut circuit_layers = Vec::<CircuitLayer>::new();
        let mut w_init_leaf = vec![Scalar::ONE; 1 << depth_1];
        for i in 0..input_length {
            w_init_leaf[i] =
                gamma_tau[0].square() * A[c][i] + gamma_tau[0] * B[c][i] + C[c][i] - gamma_tau[1];
        }
        circuit_layers.push(CircuitLayer::new(w_init_leaf));

        for k in 1..depth_1 + 1 {
            let layer_size = 1 << (depth_1 - k);
            let mut temp_w_init = vec![Scalar::ZERO; layer_size];
            for i in 0..layer_size {
                temp_w_init[i] = circuit_layers[k - 1].get_position(2 * i)
                    * circuit_layers[k - 1].get_position(2 * i + 1);
            }
            circuit_layers.push(CircuitLayer::new(temp_w_init));
        }

        circuits.push(CircuitBinaryTree::new(circuit_layers));
    }
    circuits
}
