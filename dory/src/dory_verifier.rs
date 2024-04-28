use bls_curve::bls::BlsCurve;
use channel::Channel;
use helper::{DoryProof, Setup};
use pairing::pairing_traits::Pairing;
use traits::traits::Field;

pub fn dory_verifier(proof: DoryProof, setup: Setup) {
    let intermediate_layer_commitments = proof.intermediate_layer_commitments;
    let mut channel = Channel::initialize_with_gt(&intermediate_layer_commitments[0]);
    let degree = setup.tau_1.len();
    let log2_degree = degree.trailing_zeros() as usize;
    let tau_ipp = setup.tau_ipp;
    let delta_1_left_commitment = setup.delta_1_left_commitment;
    let delta_1_right_commitment = setup.delta_1_right_commitment;
    let delta_2_right_commitment = setup.delta_2_right_commitment;
    let final_values = proof.final_values;
    let tau_1 = setup.tau_1;
    let tau_2 = setup.tau_2;
    for round in 0..log2_degree {
        let intermediate_d = &proof.intermediate_d[round];
        //Reseed channel with intermediate_d
        channel.reseed_with_gt(intermediate_d);

        //Generate beta using Fiat-shamir
        let beta = channel.get_random_point();
        //Compute inverse of beta
        let beta_inv = beta.invert().unwrap();

        let intermediate_cross_products = &proof.intermediate_cross_products[round];
        //Reseed channel with intermediate_cross_products
        channel.reseed_with_gt(intermediate_cross_products);

        //Generate alpa using Fiat-shamir
        let alpha = channel.get_random_point();
        //Compute inverse of alpha
        let alpha_inv = alpha.invert().unwrap();

        let c1 = intermediate_layer_commitments[round][0]
            + &(intermediate_layer_commitments[round][1] * &beta_inv)
            + &(intermediate_layer_commitments[round][2] * &beta)
            + &tau_ipp[round]
            + &(intermediate_cross_products[0] * &alpha.square())
            + &(intermediate_cross_products[1] * &alpha_inv.square());
        let c2 = (intermediate_d[0] * &alpha)
            + &(delta_1_left_commitment[round] * &(alpha * beta))
            + &(intermediate_d[1] * &alpha_inv)
            + &(delta_1_right_commitment[round] * &(alpha_inv * beta));
        let c3 = (intermediate_d[2] * &alpha_inv)
            + &(delta_1_left_commitment[round] * &(alpha_inv * beta_inv))
            + &(intermediate_d[3] * &alpha)
            + &(delta_2_right_commitment[round] * &(alpha * beta_inv));

        //Reseed channel with intermediate layer commitments
        channel.reseed_with_gt(&[c1, c2, c3].to_vec());

        if round < log2_degree - 1 {
            assert_eq!(
                c1,
                intermediate_layer_commitments[round + 1][0],
                "assertion failed for c1 at round {}",
                round
            );
            assert_eq!(
                c2,
                intermediate_layer_commitments[round + 1][1],
                "assertion failed for c2 at round {}",
                round
            );
            assert_eq!(
                c3,
                intermediate_layer_commitments[round + 1][2],
                "assertion failed for c3 at round {}",
                round
            );
        } else {
            let expected_c1 = BlsCurve::tate_pairing(final_values.0, final_values.1);
            let expected_c2 = BlsCurve::tate_pairing(final_values.0, tau_2[0].to_affine());
            let expected_c3 = BlsCurve::tate_pairing(tau_1[0].to_affine(), final_values.1);
            assert_eq!(c1, expected_c1, "final assertion failed for c1");
            assert_eq!(c2, expected_c2, "final assertion failed for c2");
            assert_eq!(c3, expected_c3, "final assertion failed for c3");
        }
    }
}
