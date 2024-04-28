use bls381::scalar::Scalar;
use bls_curve::{bls::BlsCurve, fp2bls::G2BlsCurve, gt::Gt};
use channel::Channel;
use curve_traits::AffinePoint;
use helper::{helper::compute_ipp, DoryProof, Setup};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::Field;

pub fn dory_prover(
    setup: Setup,
    mut afgho_commitment: Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
    initial_commitments: Vec<Gt>,
    mut g: Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
) -> DoryProof {
    let mut intermediate_layer_commitments = Vec::new();
    let mut intermediate_d = Vec::new();
    let mut intermediate_cross_products = Vec::new();
    let mut final_values = (
        AffinePoint::<BlsCurve>::IDENTITY,
        AffinePoint::<G2BlsCurve>::IDENTITY,
    );
    //tau1 is the powers of generator of G1 group
    let tau_1 = setup.tau_1;
    //tau2 is the powers of generator of G2 group
    let tau_2 = setup.tau_2;
    let delta_2_right_commitment = setup.delta_2_right_commitment;
    let delta_1_right_commitment = setup.delta_1_right_commitment;
    //degree of the polynomial
    let degree = afgho_commitment.len();
    //Compute log2(degree)
    let log2_degree = degree.trailing_zeros() as usize;
    //Store
    intermediate_layer_commitments.push(initial_commitments);
    let mut channel = Channel::initialize_with_gt(&intermediate_layer_commitments[0]);
    for round in 0..log2_degree {
        let u_left_half = afgho_commitment[0..afgho_commitment.len() / 2].to_vec();
        let u_right_half = afgho_commitment[(afgho_commitment.len() / 2)..].to_vec();
        let g_left_half = g[0..(g.len() / 2)].to_vec();
        let g_right_half = g[(g.len() / 2)..].to_vec();
        let tau_1_left_half = tau_1[0..(degree / (1 << (round + 1)))].to_vec();
        let tau_1_right_half =
            tau_1[(degree / (1 << (round + 1)))..(degree / (1 << round))].to_vec();
        let tau_2_left_half = tau_2[0..(degree / (1 << (round + 1)))].to_vec();
        let tau_2_right_half =
            tau_2[(degree / (1 << (round + 1)))..(degree / (1 << round))].to_vec();

        //Compute  d_1_l, d_1_r, d_2_l, d_2_r
        let intermediate_values = [
            (&u_left_half, &tau_2_left_half),
            (&u_right_half, &tau_2_left_half),
            (&tau_1_left_half, &g_left_half),
            (&tau_1_left_half, &g_right_half),
        ]
        .into_par_iter()
        .map(|input| compute_ipp(input.0, input.1))
        .collect::<Vec<Gt>>();

        //reseed channel with d_1_l, d_1_r, d_2_l, d_2_r
        channel.reseed_with_gt(&intermediate_values);

        //Store  d_1_l, d_1_r, d_2_l, d_2_r
        intermediate_d.push(intermediate_values.clone());

        //Generate beta using fiat-shamir
        let beta = channel.get_random_point();
        //Compute inverse of beta
        let beta_inv = beta.invert().unwrap();

        //Compute V_l = <W_1_l, W_2_r> and V_r = <W_1_r, W_2_l>
        let ipp = [
            (&u_left_half, &g_right_half),
            (&u_left_half, &tau_2_right_half),
            (&u_right_half, &g_left_half),
            (&tau_1_right_half, &g_left_half),
        ]
        .into_par_iter()
        .map(|input| compute_ipp(input.0, input.1))
        .collect::<Vec<Gt>>();

        //<u_l, g_r> + <u_l, tau2_r> * beta^-1 + <tau1_l, g_r> * beta + delta_2_right_commitment
        let v_l = ipp[0]
            + &ipp[1].mul(beta_inv)
            + &intermediate_values[3].mul(beta)
            + &delta_2_right_commitment[round];

        //<u_r, g_l> + <u_r, tau2_l> * beta^-1 + <tau1_r, g_l> * beta + delta_1_right_commitment
        let v_r = ipp[2]
            + &intermediate_values[1].mul(beta_inv)
            + &ipp[3].mul(beta)
            + &delta_1_right_commitment[round];

        //Reseed channel with v_l, v_r
        channel.reseed_with_gt(&[v_l, v_r].to_vec());

        //Store v_l, v_r]
        intermediate_cross_products.push([v_l, v_r].to_vec());

        //Generate alpha using fiat-shamir
        let alpha = channel.get_random_point();

        //Compute u' and g'
        (afgho_commitment, g) = half_length_vector(
            &afgho_commitment,
            &g,
            &tau_1[0..(degree / (1 << (round)))].to_vec(),
            &tau_2[0..(degree / (1 << (round)))].to_vec(),
            alpha,
            beta,
            beta_inv,
        );
        //Compute next layer commitments
        if round < log2_degree - 1 {
            let layer_commitments = [
                (&afgho_commitment, &g),
                (&afgho_commitment, &tau_2_left_half),
                (&tau_1_left_half, &g),
            ]
            .into_par_iter()
            .map(|input| compute_ipp(input.0, input.1))
            .collect::<Vec<Gt>>();

            channel.reseed_with_gt(&layer_commitments);
            intermediate_layer_commitments.push(layer_commitments);
        }
        //For last layer prover sends u and g.
        else {
            final_values = (
                afgho_commitment[0].to_affine().clone(),
                g[0].to_affine().clone(),
            );
        }
    }
    DoryProof {
        intermediate_layer_commitments,
        intermediate_d,
        intermediate_cross_products,
        final_values,
    }
}

pub(crate) fn half_length_vector(
    data1: &Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
    data2: &Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    tau1: &Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
    tau2: &Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    alpha: Scalar,
    beta: Scalar,
    beta_inv: Scalar,
) -> (
    Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
    Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
) {
    let left_half_data1 = data1[0..data1.len() / 2].to_vec();
    let right_half_data1 = data1[data1.len() / 2..].to_vec();
    let left_half_data2 = data2[0..data2.len() / 2].to_vec();
    let right_half_data2 = data2[data2.len() / 2..].to_vec();
    let left_half_tau1 = tau1[0..tau1.len() / 2].to_vec();
    let right_half_tau1 = tau1[tau1.len() / 2..].to_vec();
    let left_half_tau2 = tau2[0..tau2.len() / 2].to_vec();
    let right_half_tau2 = tau2[tau2.len() / 2..].to_vec();
    let alpha_inv = alpha.invert().unwrap();
    let (u, g) = (0..left_half_data1.len())
        .into_par_iter()
        .map(|idx| {
            let u = ((left_half_data1[idx] + &(left_half_tau1[idx].mul(beta))).mul(alpha))
                + &((right_half_data1[idx] + &(right_half_tau1[idx].mul(beta))).mul(alpha_inv));

            let g = ((left_half_data2[idx] + &(left_half_tau2[idx].mul(beta_inv))).mul(alpha_inv))
                + &((right_half_data2[idx] + &(right_half_tau2[idx].mul(beta_inv))).mul(alpha));
            (u, g)
        })
        .unzip();
    (u, g)
}
