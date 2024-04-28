use std::time::Instant;

use bls381::scalar::Scalar;
use bls_curve::{bls::BlsCurve, fp2bls::G2BlsCurve};
use curve_traits::ProjectivePoint;
use fft::par_fft::{get_inv_twiddles, par_interpolate_poly};
use helper::{Setup, WitnessCommitment};
use pairing::pairing_traits::Pairing;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

pub fn commit_witness(witness: &Vec<Scalar>, setup: &Setup) -> WitnessCommitment {
    //length of the witness vector
    let degree = witness.len();
    let tau1 = &setup.tau_1;
    let g2: curve_traits::ProjectivePoint<G2BlsCurve> =
        <G2BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;

    let start_time = Instant::now();
    //Compute the afgho commitment of the witness.
    let pippenger_commit_to_witness =
        <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::pippenger_msm(tau1, &witness);

    //Compute IPP(afgho_commitment_to_witness, tau2)
    let commitment_to_witness =
        BlsCurve::tate_pairing(pippenger_commit_to_witness.to_affine(), g2.to_affine());
    println!("Witness Commit time {:?}", start_time.elapsed());

    let start_time = Instant::now();
    //Compute inverse twiddles of size degree to interpolate witness.
    let inv_twiddles = get_inv_twiddles(degree as u32);
    //Copy witness into polynomial
    let mut polynomial = witness.clone();
    //Interpolate witness.
    par_interpolate_poly(&mut polynomial, inv_twiddles);

    //Compute the MSM of univariate polynomial.
    let polynomial_pipenger_commit =
        <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::pippenger_msm(
            tau1,
            &polynomial,
        );

    //Compute tate pairing of polynomial pipenger commit commit with g2
    let commitment_to_univariate =
        BlsCurve::tate_pairing(polynomial_pipenger_commit.to_affine(), g2.to_affine());

    println!("Univariate Commit time {:?}", start_time.elapsed());

    //Compute g2^{f_i}
    let g2_power_poly = (0..degree)
        .into_par_iter()
        .map(|idx| g2.mul(polynomial[idx as usize]))
        .collect::<Vec<ProjectivePoint<G2BlsCurve>>>();
    //Compute g2^{w_i}
    let g2_power_witness = (0..degree)
        .into_par_iter()
        .map(|idx| g2.mul(witness[idx as usize]))
        .collect::<Vec<ProjectivePoint<G2BlsCurve>>>();

    //Return witness commitment
    WitnessCommitment {
        commitment_to_witness,
        commitment_to_univariate,
        g2_power_poly,
        univariate_polynomial: polynomial,
        g2_power_witness,
    }
}
