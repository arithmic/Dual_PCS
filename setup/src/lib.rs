use bls381::scalar::Scalar;
use bls_curve::fp2bls::G2BlsCurve;
use bls_curve::{bls::BlsCurve, gt::Gt};
use curve_traits::ProjectivePoint;
use fft::serial_fft::log2;
use helper::helper::compute_ipp;
use helper::{DegreeBoundSetup, MatrixCommitment};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::PrimeField;
extern crate bls381;
extern crate bls_curve;
extern crate curve_traits;
extern crate fft;
extern crate helper;
extern crate pairing;
extern crate rayon;
extern crate traits;
pub fn setup(degree: usize) -> DegreeBoundSetup {
    //Generator of G1 group
    let g1 = <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;
    //Generator of G2 group
    let g2: curve_traits::ProjectivePoint<G2BlsCurve> =
        <G2BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;

    //Compute the power of G1
    let tau_1: Vec<curve_traits::ProjectivePoint<BlsCurve>> = (1..degree + 1)
        .map(|i| g1.mul(Scalar::from(i as u8)))
        .collect::<Vec<<BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint>>();
    //Compute the power of G2
    let tau_2: Vec<curve_traits::ProjectivePoint<G2BlsCurve>> = (1..degree + 1)
        .map(|i| g2.mul(Scalar::from(i as u8)))
        .collect::<Vec<<G2BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint>>();

    //Compute log2 of degree
    let log2_degree = degree.trailing_zeros() as usize;
    //Compute Commitment to matrix
    let matrix_commitments = matrix_commitment(degree, &tau_1, &tau_2);

    let ipp_commits = (0..log2_degree)
        .into_par_iter()
        .map(|iter| {
            //tau1_l^{i-1} = tau1^{i}
            let tau1_left = tau_1[0..(degree / (1 << (iter + 1)))].to_vec();
            //tau1_r^{i-1}
            let tau1_right = tau_1[(degree / (1 << (iter + 1)))..(degree / (1 << iter))].to_vec();
            //tau2^{i} = tau2_l^{i-1}
            let tau2_left = tau_2[0..(degree / (1 << (iter + 1)))].to_vec();
            //tau2_r^{i-1}
            let tau2_right = tau_2[(degree / (1 << (iter + 1)))..(degree / (1 << iter))].to_vec();

            [
                (&tau1_left, &tau2_left),
                (&tau1_right, &tau2_left),
                (&tau1_left, &tau2_right),
                (
                    &tau_1[0..(degree / (1 << (iter)))].to_vec(),
                    &tau_2[0..(degree / (1 << (iter)))].to_vec(),
                ),
            ]
            .into_par_iter()
            .map(|data| compute_ipp(data.0, &data.1))
            .collect::<Vec<Gt>>()
        })
        .collect::<Vec<Vec<Gt>>>();

    //Initialize Delta 1 left commitment with Zero of Gt
    let mut delta_1_left_commitment = vec![Gt::ZERO; log2_degree];

    //Initialize Delta 1 right commitment with Zero of Gt
    let mut delta_1_right_commitment = vec![Gt::ZERO; log2_degree];

    // //Initialize Delta 2 right commitment with Zero of Gt
    let mut delta_2_right_commitment = vec![Gt::ZERO; log2_degree];

    //<tau_1, tau_2>
    let mut tau_ipp = vec![Gt::ZERO; log2_degree];

    //Collect data
    for (iter, ipp) in ipp_commits.iter().enumerate() {
        delta_1_left_commitment[iter] = ipp[0];
        delta_1_right_commitment[iter] = ipp[1];
        delta_2_right_commitment[iter] = ipp[2];
        tau_ipp[iter] = ipp[3];
    }

    DegreeBoundSetup {
        tau_1,
        tau_2,
        delta_1_left_commitment,
        delta_1_right_commitment,
        delta_2_right_commitment,
        tau_ipp,
        matrix_commitments,
    }
}

pub fn matrix_commitment(
    degree: usize,
    tau_1: &Vec<curve_traits::ProjectivePoint<BlsCurve>>,
    tau_2: &Vec<curve_traits::ProjectivePoint<G2BlsCurve>>,
) -> Vec<MatrixCommitment> {
    let mut matrix_commitments = Vec::new();
    let degree_log2 = log2(degree);
    for index in 1..degree_log2 + 1 {
        let degree = 1 << index;
        let omega = Scalar::get_root_of_unity(log2(degree));
        let mut omega_powers = vec![Scalar::ONE; degree];
        for idx in 1..degree {
            omega_powers[idx] = omega_powers[idx - 1] * omega;
        }
        let matrix_row_commits = (0..degree)
            .into_par_iter()
            .map(|outer_idx| {
                let powers = (0..degree)
                    .map(|inner_idx| omega_powers[(inner_idx * outer_idx) % degree])
                    .collect::<Vec<Scalar>>();

                //Compute B_{d,j} = Prod(Bd[i, j] ∗ hi) where i belong {0..degree}
                let row_commit =
                    <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::pippenger_msm(
                        &tau_1[0..degree],
                        &powers,
                    );
                row_commit
            })
            .collect::<Vec<ProjectivePoint<BlsCurve>>>();
        let matrix_commitment = compute_ipp(&matrix_row_commits, &tau_2[0..degree].to_vec());
        matrix_commitments.push(MatrixCommitment {
            matrix_row_commits,
            matrix_commitment,
        });
    }
    matrix_commitments
}
