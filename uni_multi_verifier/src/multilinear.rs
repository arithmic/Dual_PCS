use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use channel::Channel;
use helper::{EvaluationProof, Setup, GT};
use pairing::pairing_traits::Pairing;
use traits::traits::Field;

pub fn multilinear_evaluation_verifier(
    evaluation_proof: EvaluationProof,
    setup: &Setup,
) -> Option<Scalar> {
    let intermediate_layer_commitments = evaluation_proof
        .clone()
        .dory_proof
        .intermediate_layer_commitments;
    let mut channel = Channel::initialize_with_gt(&intermediate_layer_commitments[0]);
    let degree = setup.tau_1.len();
    let log2_degree = degree.trailing_zeros() as usize;
    let r = channel.get_random_points(log2_degree);
    let tau_ipp = &setup.tau_ipp;
    let delta_1_left_commitment = &setup.delta_1_left_commitment;
    let delta_1_right_commitment = &setup.delta_1_right_commitment;
    let delta_2_right_commitment = &setup.delta_2_right_commitment;
    let final_values = evaluation_proof.clone().dory_proof.final_values;

    let g1 = <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;

    let tau_1 = &setup.tau_1;
    let tau_2 = &setup.tau_2;
    let mut alpha_vec = Vec::new();
    let mut alpha_inv_vec = Vec::new();

    let mut v = GT.mul(evaluation_proof.evaluation);
    let mut flag = false;
    for round in 0..log2_degree {
        let intermediate_d = &evaluation_proof.dory_proof.intermediate_d[round];
        //Reseed channel with intermediate_d
        channel.reseed_with_gt(intermediate_d);

        //Generate beta using Fiat-shamir
        let beta = channel.get_random_point();
        //Compute inverse of beta
        let beta_inv = beta.invert().unwrap();

        let intermediate_cross_products =
            &evaluation_proof.dory_proof.intermediate_cross_products[round];
        //Reseed channel with intermediate_cross_products
        channel.reseed_with_gt(intermediate_cross_products);

        //Generate alpa using Fiat-shamir
        let alpha = channel.get_random_point();
        alpha_vec.push(alpha);

        let alpha_square = alpha.square();

        //Compute inverse of alpha
        let alpha_inv = alpha.invert().unwrap();
        alpha_inv_vec.push(alpha_inv);
        let alpha_inv_square = alpha_inv.square();
        let c1 = intermediate_layer_commitments[round][0]
            + &(intermediate_layer_commitments[round][1] * &beta_inv)
            + &(intermediate_layer_commitments[round][2] * &beta)
            + &tau_ipp[round]
            + &(intermediate_cross_products[0] * &alpha_square)
            + &(intermediate_cross_products[1] * &alpha_inv_square);
        let c2 = (intermediate_d[0] * &alpha)
            + &(delta_1_left_commitment[round] * &(alpha * beta))
            + &(intermediate_d[1] * &alpha_inv)
            + &(delta_1_right_commitment[round] * &(alpha_inv * beta));
        let c3 = (intermediate_d[2] * &alpha_inv)
            + &(delta_1_left_commitment[round] * &(alpha_inv * beta_inv))
            + &(intermediate_d[3] * &alpha)
            + &(delta_2_right_commitment[round] * &(alpha * beta_inv));

        v = v
            + &((intermediate_d[4] + &(intermediate_d[5] * &beta_inv)) * &alpha_square)
            + &((intermediate_d[6] + &(intermediate_d[7] * &beta_inv)) * &alpha_inv_square)
            + &(intermediate_d[8] * &beta_inv);

        //Reseed channel with intermediate layer commitments
        channel.reseed_with_gt(&[c1, c2, c3].to_vec());
        if round < log2_degree - 1 {
            assert_eq!(
                c1,
                intermediate_layer_commitments[round + 1][0],
                "assertion failed for c1"
            );
            assert_eq!(
                c2,
                intermediate_layer_commitments[round + 1][1],
                "assertion failed for c2"
            );
            assert_eq!(
                c3,
                intermediate_layer_commitments[round + 1][2],
                "assertion failed for c3"
            );
        } else {
            let expected_c1 = BlsCurve::tate_pairing(final_values.0, final_values.1);
            let expected_c2 = BlsCurve::tate_pairing(final_values.0, tau_2[0].to_affine());
            let expected_c3 = BlsCurve::tate_pairing(tau_1[0].to_affine(), final_values.1);
            let y_logn = (0..log2_degree).fold(Scalar::ONE, |acc, idx| {
                acc * ((alpha_vec[idx] * (Scalar::ONE - r[log2_degree - (idx + 1)]))
                    + (alpha_inv_vec[idx] * r[log2_degree - (idx + 1)]))
            });
            let h_logn = g1.mul(y_logn);
            let expected_v = BlsCurve::tate_pairing(h_logn.to_affine(), final_values.1);
            assert_eq!(v, expected_v, "assertion failed");
            if c1 == expected_c1 && c2 == expected_c2 && c3 == expected_c3 && v == expected_v {
                flag = true
            }
        }
    }
    if flag {
        return Some(evaluation_proof.evaluation);
    }
    None
}
