use bls381::scalar::Scalar;
use channel::Channel;
use crypto_bigint::{Limb, Uint};
use gkr_common::{CircuitBinaryTree, GkrTranscript, MleLayerEvaluation};
use polynomial::{LayerPolynomials, MultPolynomial, Polynomial};
use rayon::{
    iter::{
        IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator,
        ParallelIterator,
    },
    slice::ParallelSliceMut,
};
use traits::traits::Field;

use crate::gkr_common;

pub mod reed_solomon;

//Prover for the sub-circuit corresponding to the leaf layer inputs of length of the table.
pub fn gkr_prover(circuits: &Vec<CircuitBinaryTree>, channel: &mut Channel) -> GkrTranscript {
    let depth = circuits[0].get_depth();
    circuits
        .par_iter()
        .for_each(|circuit| assert_eq!(depth, circuit.get_depth(), "Circuits do not have same"));

    let n_circuits = circuits.len();
    //This vector contains the values of the circuits at depth 1 i.e. the layer below the output layer.
    let mut final_evaluations = Vec::new();

    for c in 0..n_circuits {
        let circuit_layer_1 = circuits[c].get_layer_at_depth(depth - 1);
        final_evaluations.push(vec![
            circuit_layer_1.get_position(0),
            circuit_layer_1.get_position(1),
        ]);
    }
    let final_values: Vec<Scalar> = final_evaluations
        .iter()
        .flat_map(|value| value.iter())
        .cloned()
        .collect();
    channel.reseed_with_scalars(&final_values);

    let mut initial_random_point = vec![channel.get_random_point()];

    channel.reseed_with_scalars(&initial_random_point);

    //We are verifying the circuit evaluations in a batched manner by taking a linear combination of
    //the gate evaluation MLEs for each circuit, at each layer.
    let random_coeff = channel.get_random_points(n_circuits);

    channel.reseed_with_scalars(&random_coeff);

    //This vector contains the values the prover claims for the MLE evaluation at each layer.
    let mut claimed_values = Vec::new();

    //This vector contains the polynomials the prover sends for the sum check instance at each layer.
    let mut sum_check_polynomials: Vec<LayerPolynomials> = Vec::new();

    for layer in (1..depth).rev() {
        //The polynomials for the sum-check instance at the current layer.
        let mut polynomials_current_layer = LayerPolynomials::new(Vec::new());

        //The dense representation of the lagrange basis functions evaluated at the random point.
        let mut lagrange_bases_eval =
            MultPolynomial::new(compute_fourier_bases(initial_random_point.clone()));

        let current_depth = depth - layer;
        let layer_size = 1 << current_depth;

        //Random points generated over the sum check instance.
        let mut sum_check_random_points: Vec<Scalar> = vec![Scalar::ONE; current_depth + 1];

        //Dense representation of the next depth's left child gate and right child gate mle's respectively W_{d-1}(x;0), W_{d-1}(x:1)
        let mut child_left_extension =
            vec![MultPolynomial::new(vec![Scalar::ONE; layer_size]); n_circuits];
        let mut child_right_extension =
            vec![MultPolynomial::new(vec![Scalar::ONE; layer_size]); n_circuits];

        for i in 0..layer_size {
            for c in 0..n_circuits {
                child_left_extension[c].0[i] = circuits[c]
                    .get_layer_at_depth(layer - 1)
                    .get_position(2 * i);
                child_right_extension[c].0[i] = circuits[c]
                    .get_layer_at_depth(layer - 1)
                    .get_position(2 * i + 1)
            }
        }

        //This is the sum check instance for this layer.
        for i in 0..current_depth {
            //This contains the circuit-wise values for univariate polynomial to be sent evaluated at 0,1,-1 and 2.
            let mut eval = vec![[Scalar::ZERO; 4]; n_circuits];

            //This contains the coefficients of the linear combination of the circuit-wise evaluations.
            let mut combined_polynomial = vec![Scalar::ZERO; 4];
            let halfsize = layer_size >> (i + 1);

            //We evaluate the dense representation of the multilinear polynomial that is the linear coefficient of lagrange_bases_eval
            //when viewed as a polynomial in one variable.
            let mut lagrange_bases_lin_coeff = MultPolynomial(vec![Scalar::ZERO; halfsize]);

            lagrange_bases_lin_coeff
                .as_coeffs_mut()
                .par_iter_mut()
                .enumerate()
                .for_each(|(j, fc_coeff)| {
                    *fc_coeff = lagrange_bases_eval.get_coeff((j << 1) + 1)
                        - lagrange_bases_eval.get_coeff(j << 1)
                });

            //We compute the sum over all binary strings for the remaining variables, parallelised over the circuits.
            eval.par_iter_mut().enumerate().for_each(|(c, eval_c)| {
                for j in 0..halfsize {
                    //We use the fact that for any multilinear polynomial W in variables, x_1, ..., x_d+1,
                    //W(x_1, ..., x_d+1) = (1-x_1).W(x_1, ... , x_d,0) + x_1.W(x_1,...,x_d,1).

                    let child_left_temp = child_left_extension[c].get_coeff((j << 1) + 1)
                        - child_left_extension[c].get_coeff(j << 1);
                    let child_right_temp = child_right_extension[c].get_coeff((j << 1) + 1)
                        - child_right_extension[c].get_coeff(j << 1);

                    //Evaluation at 0
                    eval_c[0] += lagrange_bases_eval.get_coeff(j << 1)
                        * child_left_extension[c].get_coeff(j << 1)
                        * child_right_extension[c].get_coeff(j << 1);

                    //Evaluation at 1
                    eval_c[1] += lagrange_bases_eval.get_coeff((j << 1) + 1)
                        * child_left_extension[c].get_coeff((j << 1) + 1)
                        * child_right_extension[c].get_coeff((j << 1) + 1);

                    //Evaluation at -1
                    eval_c[2] += (lagrange_bases_eval.get_coeff(j << 1)
                        - lagrange_bases_lin_coeff.get_coeff(j))
                        * (child_left_extension[c].get_coeff(j << 1) - child_left_temp)
                        * (child_right_extension[c].get_coeff(j << 1) - child_right_temp);

                    //Evaluation at 2
                    eval_c[3] += (lagrange_bases_eval.get_coeff((j << 1) + 1)
                        + lagrange_bases_lin_coeff.get_coeff(j))
                        * (child_left_extension[c].get_coeff((j << 1) + 1) + child_left_temp)
                        * (child_right_extension[c].get_coeff((j << 1) + 1) + child_right_temp);
                }
            });

            //We conduct the inverse linear transform fromthe evaluations to get the coefficients of the circuit-wise polynomials.
            eval.par_iter_mut()
                .for_each(|eval_c| len_4_interpolate(eval_c));

            //We add the linear combination of the coefficients to get the batched polynomial for the sum check verifieer to verify.
            for c in 0..n_circuits {
                for k in 0..4 {
                    combined_polynomial[k] += random_coeff[c] * eval[c][k];
                }
            }

            //Here we push the current round's batched polynomial to the vector of polynomials for the current layer.
            polynomials_current_layer.add_polynomial(Polynomial::new(combined_polynomial.clone()));

            //The next round's random point is obtained by reseeding with this round's batched polynomial as a vector of scalars
            // and then drawing from the channel.
            channel.reseed_with_scalars(&combined_polynomial);

            let random_point = channel.get_random_point();
            sum_check_random_points[i + 1] = random_point;

            //Now we fix the leading variable of the the multilinear polynomials lagrange_bases_eval, child_left and child_right to get multilinear polynomials
            //in one less variable for the next sum check round

            //Since lagrange_bases_eval is a common computation accross circuits, we fold it in parallel separately.

            if i < 8 {
                lagrange_bases_eval = lagrange_bases_eval.fold_by_lsb(random_point)
            } else {
                lagrange_bases_eval = lagrange_bases_eval.par_fold_by_lsb(random_point)
            }

            //We fix the current leading variable in the sum check to the current random point of the protocol for the extension for the
            //child_left and child_right MLEs
            //
            child_left_extension
                .par_iter_mut()
                .zip(child_right_extension.par_iter_mut())
                .for_each(|(child_left_extension_c, child_right_extension_c)| {
                    *child_left_extension_c = child_left_extension_c.fold_by_lsb(random_point);
                    *child_right_extension_c = child_right_extension_c.fold_by_lsb(random_point);
                });
        }

        //We push the current layer's sum check polynomials to the vector that holds all the layer wise
        // sum check polynomials.
        sum_check_polynomials.push(polynomials_current_layer);

        //These are the variables that will contain the values for the appropriate circuit-wise linear-
        // combination for the claimed values, i.e. the scalar obtained after binding all but the last
        // variable to the random values obtained over the sum-check for the current layer. The last variable
        // is understood to be fixed to 0 and 1 respectively to obtain the MLEs of child_left and child_right
        // respectively.

        let mut mle_layer_evaluation = Vec::new();
        for c in 0..n_circuits {
            let mle_eval = MleLayerEvaluation::new(
                child_left_extension[c].get_coeff(0),
                child_right_extension[c].get_coeff(0),
            );
            mle_layer_evaluation.push(mle_eval);
        }

        //We push the current layer's MLE evaluations to the
        claimed_values.push(mle_layer_evaluation);

        //The line is of the form L(t) = (r_d_i;t), thus q(t)= W_{d-1}( L(t) ) is of degree 1 as W is linear in each variable.
        channel.reseed_with_scalars(&sum_check_random_points);
        let r = channel.get_random_point();
        sum_check_random_points[0] = r;
        initial_random_point = sum_check_random_points
    }
    initial_random_point.reverse();
    GkrTranscript::new(
        final_evaluations,
        claimed_values,
        sum_check_polynomials,
        initial_random_point,
    )
}

pub fn compute_fourier_bases(r: Vec<Scalar>) -> Vec<Scalar> {
    //Initialize fc_eq with (1- r[0]) and r[0]
    let mut fc_eq = [Scalar::ONE - r[r.len() - 1], r[r.len() - 1]].to_vec();
    //Iterate over the length of the r vector
    for k in (0..r.len() - 1).rev() {
        let temp = fc_eq;
        //initialize fc_eq of double size with zero
        fc_eq = vec![Scalar::ZERO; temp.len() * 2];

        if k < 8 {
            for iter in 0..temp.len() {
                fc_eq[2 * iter + 1] = temp[iter] * r[k];
                fc_eq[2 * iter] = temp[iter] - fc_eq[2 * iter + 1];
            }
        } else {
            fc_eq
                .par_chunks_mut(2)
                .zip(temp)
                .for_each(|(fc_eq_pair, temp)| {
                    fc_eq_pair[1] = temp * (r[k as usize]);
                    fc_eq_pair[0] = temp - fc_eq_pair[1];
                })
        }
    }
    fc_eq
}

const SIX_INV: Scalar = Scalar(Uint {
    limbs: [
        Limb(15372286724512153601),
        Limb(1954008828163476649),
        Limb(9224930440102993583),
        Limb(6961264049553707793),
    ],
});
const TWO_INV: Scalar = Scalar(Uint {
    limbs: [
        Limb(9223372034707292161),
        Limb(12240451741123816959),
        Limb(1845609449319885826),
        Limb(4176758429732224676),
    ],
});

pub fn len_4_interpolate(evaluations: &mut [Scalar; 4]) {
    let t0 = TWO_INV * (evaluations[1] + evaluations[2] - evaluations[0].double());
    let t1 = evaluations[1] - evaluations[2] + evaluations[0] + t0.double().double();
    let t2 = SIX_INV * (evaluations[3] - t1);
    *evaluations = [
        evaluations[0],
        evaluations[1] - (evaluations[0] + t0 + t2),
        t0,
        t2,
    ]
}
