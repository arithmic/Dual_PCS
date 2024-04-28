use bls381::scalar::Scalar;
use bls_curve::{bls::BlsCurve, fp2bls::G2BlsCurve, gt::Gt};
use channel::Channel;
use curve_traits::{AffinePoint, ProjectivePoint};
use fft::serial_fft::eval;
use helper::{helper::compute_ipp, DoryProof, EvaluationProof, Setup};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::Field;

pub fn univariate_evaluation_prover(
    setup: &Setup,
    mut afgho_commitment: Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
    initial_commitments: Vec<Gt>,
    mut g: Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    polynomial: Vec<Scalar>,
) -> EvaluationProof {
    let mut dory_commitments = Vec::new();
    let mut intermediate_d = Vec::new();
    let mut intermediate_cross_products = Vec::new();
    let mut final_values = (
        AffinePoint::<BlsCurve>::IDENTITY,
        AffinePoint::<G2BlsCurve>::IDENTITY,
    );

    //Generator of G1 group
    let g1 = <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;

    //tau1 is the powers of generator of G1 group
    let tau_1 = &setup.tau_1;

    //tau2 is the powers of generator of G2 group
    let tau_2 = &setup.tau_2;
    let delta_2_right_commitment = &setup.delta_2_right_commitment;
    let delta_1_right_commitment = &setup.delta_1_right_commitment;

    //degree of the polynomial
    let degree = afgho_commitment.len();

    //Compute log2(degree)
    let log2_degree = degree.trailing_zeros() as usize;

    let mut channel = Channel::initialize_with_gt(&initial_commitments);
    dory_commitments.push(initial_commitments);

    let random_point = channel.get_random_point();
    let evaluation = eval(&polynomial, random_point);

    let y = (0..degree)
        .into_par_iter()
        .map(|idx| random_point.power_by(&[idx as u64, 0, 0, 0]))
        .collect::<Vec<Scalar>>();

    let mut h = (0..y.len())
        .into_par_iter()
        .map(|idx| g1.mul(y[idx]))
        .collect::<Vec<ProjectivePoint<BlsCurve>>>();

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
        let h_left_half = h[0..h.len() / 2].to_vec();
        let h_right_half = h[h.len() / 2..].to_vec();
        let intermediate_d_values = [
            (&u_left_half, &tau_2_left_half),
            (&u_right_half, &tau_2_left_half),
            (&tau_1_left_half, &g_left_half),
            (&tau_1_left_half, &g_right_half),
            (&h_left_half, &g_right_half),
            (&h_left_half, &tau_2_right_half),
            (&h_right_half, &g_left_half),
            (&h_right_half, &tau_2_left_half),
            (&h, &tau_2[0..(degree / (1 << (round)))].to_vec()),
        ]
        .into_par_iter()
        .map(|input| compute_ipp(input.0, input.1))
        .collect::<Vec<Gt>>();

        //reseed channel with d_1_l, d_1_r, d_2_l, d_2_r
        channel.reseed_with_gt(&intermediate_d_values);

        //Store  d_1_l, d_1_r, d_2_l, d_2_r
        intermediate_d.push(intermediate_d_values.clone());

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

        let v_l = ipp[0]
            + &ipp[1].mul(beta_inv)
            + &intermediate_d_values[3].mul(beta)
            + &delta_2_right_commitment[round];

        let v_r = ipp[2]
            + &intermediate_d_values[1].mul(beta_inv)
            + &ipp[3].mul(beta)
            + &delta_1_right_commitment[round];

        //Reseed channel with v_l, v_r
        channel.reseed_with_gt(&[v_l, v_r].to_vec());

        //Store v_l, v_r]
        intermediate_cross_products.push([v_l, v_r].to_vec());

        //Generate alpha using fiat-shamir
        let alpha = channel.get_random_point();

        //Compute u' and g'
        (afgho_commitment, g, h) = half_length_vector(
            &afgho_commitment,
            &g,
            &h,
            &tau_1[0..(degree / (1 << (round)))].to_vec(),
            &tau_2[0..(degree / (1 << (round)))].to_vec(),
            alpha,
            beta,
        );
        //Compute next layer commitments
        if round < log2_degree - 1 {
            let intermediate_commits = [
                (&afgho_commitment, &g),
                (&afgho_commitment, &tau_2_left_half),
                (&tau_1_left_half, &g),
            ]
            .into_par_iter()
            .map(|input| compute_ipp(input.0, input.1))
            .collect::<Vec<Gt>>();

            channel.reseed_with_gt(&intermediate_commits);
            dory_commitments.push(intermediate_commits);
        }
        //For last layer prover sends u and g.
        else {
            final_values = (
                afgho_commitment[0].to_affine().clone(),
                g[0].to_affine().clone(),
            );
        }
    }
    EvaluationProof {
        dory_proof: DoryProof {
            intermediate_layer_commitments: dory_commitments,
            intermediate_d,
            intermediate_cross_products,
            final_values,
        },
        evaluation,
    }
}
/// Compute
pub(crate) fn half_length_vector(
    data1: &Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
    data2: &Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    data3: &Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
    tau1: &Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
    tau2: &Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    alpha: Scalar,
    beta: Scalar,
) -> (
    Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
    Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
    Vec<curve_traits::ProjectivePoint<bls_curve::bls::BlsCurve>>,
) {
    let left_half_data1 = data1[0..data1.len() / 2].to_vec();
    let right_half_data1 = data1[data1.len() / 2..].to_vec();
    let left_half_data2 = data2[0..data2.len() / 2].to_vec();
    let right_half_data2 = data2[data2.len() / 2..].to_vec();
    let left_half_data3 = data3[0..data1.len() / 2].to_vec();
    let right_half_data3 = data3[data1.len() / 2..].to_vec();
    let left_half_tau1 = tau1[0..tau1.len() / 2].to_vec();
    let right_half_tau1 = tau1[tau1.len() / 2..].to_vec();
    let left_half_tau2 = tau2[0..tau2.len() / 2].to_vec();
    let right_half_tau2 = tau2[tau2.len() / 2..].to_vec();
    let alpha_inv = alpha.invert().unwrap();
    let beta_inv = beta.invert().unwrap();
    let (uh, g): (
        (
            Vec<ProjectivePoint<BlsCurve>>,
            Vec<ProjectivePoint<BlsCurve>>,
        ),
        Vec<ProjectivePoint<G2BlsCurve>>,
    ) = (0..left_half_data1.len())
        .into_par_iter()
        .map(|idx| {
            let u = ((left_half_data1[idx] + &(left_half_tau1[idx].mul(beta))).mul(alpha))
                + &((right_half_data1[idx] + &(right_half_tau1[idx].mul(beta))).mul(alpha_inv));
            let g = ((left_half_data2[idx] + &(left_half_tau2[idx].mul(beta_inv))).mul(alpha_inv))
                + &((right_half_data2[idx] + &(right_half_tau2[idx].mul(beta_inv))).mul(alpha));
            let h = left_half_data3[idx].mul(alpha) + right_half_data3[idx].mul(alpha_inv);
            ((u, h), g)
        })
        .unzip();
    (uh.0, g, uh.1)
}
