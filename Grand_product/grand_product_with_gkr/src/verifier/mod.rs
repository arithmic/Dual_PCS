use bls381::scalar::Scalar;
use channel::Channel;
use fft::serial_fft::eval;
use gkr_common::GkrTranscript;
use traits::traits::Field;

use crate::gkr_common;

pub fn gkr_verifier(
    transcript: &GkrTranscript,
    depth: usize,
    channel: &mut Channel,
    circuit_evals: Vec<Scalar>,
    n_circuits: usize,
) {
    let polynomials = &transcript.polynomials;
    let final_evaluations = &transcript.final_evaluations;
    let claimed_values = &transcript.claimed_values;
    let final_point = &transcript.final_layer_point;
    let final_values: Vec<Scalar> = final_evaluations
        .iter()
        .flat_map(|value| value.iter())
        .cloned()
        .collect();
    channel.reseed_with_scalars(&final_values);

    let mut initial_random_point = vec![channel.get_random_point()];

    channel.reseed_with_scalars(&initial_random_point);

    //Verifier obtains the random coefficients the prover uses.
    let random_coeff = channel.get_random_points(n_circuits);

    let mut binding_per_layer = Scalar::ZERO;

    //The value for the claim of the first round of the protocol.
    for c in 0..n_circuits {
        binding_per_layer += random_coeff[c]
            * ((Scalar::ONE - initial_random_point[0]) * final_evaluations[c][0]
                + initial_random_point[0] * final_evaluations[c][1])
    }

    channel.reseed_with_scalars(&random_coeff);

    for d in 0..depth - 1 {
        let polynomials_for_layer = &polynomials[d];
        let rounds = polynomials_for_layer.number_of_rounds();
        let claimed_values_d = &claimed_values[d];

        let mut current_sum = binding_per_layer;
        let mut sum_check_random_points = vec![Scalar::ONE; rounds + 1];

        for i in 0..rounds {
            let poly = polynomials_for_layer
                .polynomial_for_round(i)
                .get_coefficients();
            assert_eq!(
                current_sum,
                poly[0].double() + poly[1] + poly[2] + poly[3],
                "Sum check failed on round {i} at depth {d}"
            );
            channel.reseed_with_scalars(&poly);

            let r = channel.get_random_point();
            current_sum = eval(poly, r);
            sum_check_random_points[i + 1] = r;
        }

        channel.reseed_with_scalars(&sum_check_random_points);
        let r = channel.get_random_point();

        sum_check_random_points[0] = r;

        let eq = evaluate_eq(
            initial_random_point,
            sum_check_random_points[1..].to_vec().clone(),
        );
        let mut temp = Scalar::ZERO;
        for c in 0..claimed_values_d.len() {
            temp += random_coeff[c] * (claimed_values_d[c].left * claimed_values_d[c].right)
        }

        assert_eq!(current_sum, eq * temp, "assertion failed at layer {d}");
        initial_random_point = sum_check_random_points;

        //After sum check for the layer completes successfully, the verifier can compute the challenge
        //for the next layer using the claimed values sent by the prover. Which correspond to,
        //W(x_1, x_2, ..., x_d, 0),W(x_1,x_2,...,x_d,1) respectively.
        //For any multilinear polynomial W in variables, x_1, ..., x_d+1.
        //W(x_1, ..., x_d+1) = (1-x_1).W(x_1, ... , x_d,0) + x_1.W(x_1,...,x_d,1)
        let mut next_layer_claimed_values = Scalar::ZERO;
        for c in 0..claimed_values_d.len() {
            next_layer_claimed_values += random_coeff[c]
                * ((Scalar::ONE - r) * claimed_values_d[c].left + r * claimed_values_d[c].right)
        }
        binding_per_layer = next_layer_claimed_values;
    }
    initial_random_point.reverse();
    assert_eq!(*final_point, initial_random_point);

    let mut final_claimed_values = Scalar::ZERO;
    for c in 0..n_circuits {
        final_claimed_values += random_coeff[c] * circuit_evals[c]
    }
    assert_eq!(
        binding_per_layer, final_claimed_values,
        "Final depth check failed"
    )
}
pub fn evaluate_eq(r_x: Vec<Scalar>, r_y: Vec<Scalar>) -> Scalar {
    let mut temp = Scalar::ONE;
    assert_eq!(r_x.len(), r_y.len());
    for k in 0..r_y.len() {
        temp = temp * ((r_x[k] * r_y[k]) + ((Scalar::ONE - r_x[k]) * (Scalar::ONE - r_y[k])));
    }
    temp
}
