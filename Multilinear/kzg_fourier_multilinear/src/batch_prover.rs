use crate::{
    multilinear_common::reseed_with_even_odd_polys,
    prover::{
        commit_even_odd_polys, compute_even_odd_polys, compute_quotient_polys,
        evaluate_phi_at_random_point, evaluate_polynomials_at_s, evaluations_at_random_point,
    },
    setup::compute_evaluations,
};
use bls381::scalar::Scalar;
use channel::Channel;
use common::{
    BatchMultilinearKZG2Proof, CommitmentsOfEvenOddPolys, EvaluationOfEvenOddPolys,
    EvaluationsAtRandomPoint, EvenOddPolys, KZGFourierDegreeBoundSetup,
};
use fft::{
    par_fft::{get_inv_twiddles, par_interpolate_poly},
    serial_fft::eval,
};
use kzg_fft::commit::{self};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::Field;
use traits::traits::PrimeField;

pub fn batch_eval_multiple_evaluation_points(
    polys: &Vec<Vec<Scalar>>,
    evaluations: Vec<Scalar>,
    random_points: Vec<Vec<Scalar>>,
    setup: &KZGFourierDegreeBoundSetup,
    channel: &mut Channel,
) -> BatchMultilinearKZG2Proof {
    let num_var = polys[0].len().trailing_zeros();

    //Compute q_k s.t. f - f(u) = sum(X_k - u_k) * q_k , k = {0, 1, .., n-1}
    let quotient_polys = (0..polys.len())
        .map(|idx| {
            compute_quotient_polys(&polys[idx], random_points[idx].clone(), evaluations[idx])
        })
        .collect::<Vec<Vec<Vec<Scalar>>>>();

    //Evaluate g_{q_k,e}^{k} and g_{q_k,o}^{k}
    let (even_odd_polys, evaluations_of_even_odd_polys): (
        Vec<EvenOddPolys>,
        Vec<EvaluationOfEvenOddPolys>,
    ) = (0..polys.len())
        .map(|idx| compute_even_odd_polys(quotient_polys[idx].clone()))
        .unzip();

    //Commit evaluations of even and odd polys
    let commitments_of_even_ood_polys = (0..polys.len())
        .map(|idx| commit_even_odd_polys(&evaluations_of_even_odd_polys[idx], &setup))
        .collect::<Vec<CommitmentsOfEvenOddPolys>>();

    //Reseed channel with commitments
    (0..polys.len()).for_each(|idx| {
        reseed_with_even_odd_polys(channel, commitments_of_even_ood_polys[idx].clone())
    });

    //Draw a random point using Fiat-shamir
    let z = channel.get_random_point();
    let q_not_polys = (0..polys.len())
        .map(|idx| even_odd_polys[idx].const_poly.clone())
        .collect::<Vec<Vec<Scalar>>>();
    let (evaluations_at_z, interpolated_poly): (Vec<EvaluationsAtRandomPoint>, Vec<Vec<Scalar>>) =
        (0..polys.len())
            .map(|idx| evaluations_at_random_point(z, &even_odd_polys[idx], polys[idx].clone()))
            .unzip();

    let phi_k_evaluation_at_z = evaluate_phi_at_random_point(num_var as usize, z);
    //Draw random points
    let no_of_points = (4 * num_var as usize - 3) * polys.len() + (num_var as usize - 1);
    let r_points = channel.get_random_points(no_of_points as usize);

    let (shi, extended_phi_k_polys_evals) = compute_shi(
        polys,
        &evaluations_at_z,
        &r_points,
        &z,
        &evaluations_of_even_odd_polys,
        &phi_k_evaluation_at_z,
    );

    let shi_commit = commit::kzg2commit(&shi, &setup.setup.get_setup(shi.len()).prover_key);

    channel.reseed_with_affine_point(&[shi_commit].to_vec());

    //Draw a random point using Fiat-shamir
    let s = channel.get_random_point();

    //Draw random points to combine all the polys using Fiat-Shamir
    let poly_combiners = channel
        .get_random_points((num_var as usize * 2 - 1) * polys.len() + (num_var as usize - 1));

    let evaluations_at_s = (0..polys.len())
        .map(|idx| {
            evaluate_polynomials_at_s(s, &even_odd_polys[idx], interpolated_poly[idx].clone())
        })
        .collect::<Vec<EvaluationsAtRandomPoint>>();

    let phi_k_evaluation_at_s = evaluate_phi_at_random_point(num_var as usize, s);

    let r = combine_polynomials(
        shi,
        polys,
        &evaluations_of_even_odd_polys,
        &poly_combiners,
        s,
        &evaluations_at_s,
        extended_phi_k_polys_evals,
        &phi_k_evaluation_at_s,
    );
    let commitment_to_r = commit::kzg2commit(&r, &setup.setup.get_setup(r.len()).prover_key);
    BatchMultilinearKZG2Proof {
        evaluations,
        commitments_of_even_ood_polys,
        evaluations_at_z,
        evaluations_at_s,
        commitment_to_r,
        shi_commit,
        q_not_polys,
        phi_k_evaluation_at_z,
        phi_k_evaluation_at_s,
    }
}

fn combine_polynomials(
    mut shi: Vec<Scalar>,
    polys: &Vec<Vec<Scalar>>,
    even_odd_polys_evals: &Vec<EvaluationOfEvenOddPolys>,
    random_points: &[Scalar],
    s: Scalar,
    evaluations_at_s: &Vec<EvaluationsAtRandomPoint>,
    extended_phi_k_polys_evals: Vec<Vec<Scalar>>,
    phi_k_evaluation_at_s: &Vec<Scalar>,
) -> Vec<Scalar> {
    let num_polys = polys.len();
    let num_var = polys[0].len().trailing_zeros();

    let even_poly_evaluations = &(0..polys.len())
        .map(|idx| even_odd_polys_evals[idx].clone().even_polys_evaluations)
        .collect::<Vec<Vec<Vec<Scalar>>>>();
    let odd_poly_evaluations = &(0..polys.len())
        .map(|idx| even_odd_polys_evals[idx].clone().odd_polys_evaluations)
        .collect::<Vec<Vec<Vec<Scalar>>>>();
    assert_eq!(
        even_poly_evaluations.len(),
        odd_poly_evaluations.len(),
        "no of even polys and no of odd polys should be equal"
    );
    let g_num_var = Scalar::get_root_of_unity(num_var);
    let evaluations_of_g = (0..1 << num_var)
        .into_par_iter()
        .map(|outer_idx| {
            let mut temp = Scalar::ZERO;
            for inner_idx in 0..even_odd_polys_evals[0].odd_polys_evaluations.len() {
                for idx in 0..polys.len() {
                    temp += (random_points[2 * inner_idx + num_polys + num_var as usize - 1 + idx]
                        * even_poly_evaluations[idx][inner_idx][outer_idx])
                        + (random_points[2 * inner_idx + num_polys + num_var as usize + idx]
                            * odd_poly_evaluations[idx][inner_idx][outer_idx])
                }
                temp += random_points[num_polys + inner_idx]
                    * extended_phi_k_polys_evals[inner_idx][outer_idx]
            }
            temp + ((0..polys.len()).fold(Scalar::ZERO, |acc, iter| {
                acc + random_points[iter] * polys[iter][outer_idx]
            })) + shi[outer_idx]
        })
        .collect::<Vec<Scalar>>();

    let even_polys_evaluations = &(0..polys.len())
        .map(|idx| evaluations_at_s[idx].clone().even_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();
    let odd_polys_evaluations = &(0..polys.len())
        .map(|idx| evaluations_at_s[idx].clone().odd_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();

    let f_s = &(0..polys.len())
        .map(|idx| evaluations_at_s[idx].evaluation_of_f)
        .collect::<Vec<Scalar>>();

    let inv_twiddles = get_inv_twiddles(shi.len() as u32);
    par_interpolate_poly(&mut shi, inv_twiddles);
    let shi_at_s = eval(&shi, s);

    let g_at_s = (0..phi_k_evaluation_at_s.len())
        .into_par_iter()
        .map(|iter| {
            (0..polys.len()).fold(Scalar::ZERO, |acc, idx| {
                acc + (random_points[2 * iter + num_polys + num_var as usize - 1 + idx]
                    * even_polys_evaluations[idx][iter])
                    + (random_points[2 * iter + num_polys + num_var as usize + idx]
                        * odd_polys_evaluations[idx][iter])
            }) + (random_points[num_polys + iter] * phi_k_evaluation_at_s[iter])
        })
        .collect::<Vec<Scalar>>()
        .into_iter()
        .fold(Scalar::ZERO, |acc, value| acc + value)
        + (0..polys.len()).fold(Scalar::ZERO, |acc, idx| {
            acc + (random_points[idx] * f_s[idx])
        })
        + shi_at_s;

    let mut omega_powers = vec![Scalar::ONE; 1 << num_var];
    for idx in 1..omega_powers.len() {
        omega_powers[idx] = omega_powers[idx - 1] * g_num_var;
    }
    let quotient_poly = (0..1 << num_var)
        .into_par_iter()
        .map(|idx| ((evaluations_of_g[idx] - g_at_s) * (omega_powers[idx] - s).invert().unwrap()))
        .collect::<Vec<Scalar>>();

    quotient_poly
}
fn compute_shi(
    polys: &Vec<Vec<Scalar>>,
    evaluations_at_z: &Vec<EvaluationsAtRandomPoint>,
    random_points: &Vec<Scalar>,
    z: &Scalar,
    evaluations_of_even_odd_polys: &Vec<EvaluationOfEvenOddPolys>,
    phi_k_evaluation_at_z: &Vec<Scalar>,
) -> (Vec<Scalar>, Vec<Vec<Scalar>>) {
    let num_var = polys[0].len().trailing_zeros();

    let omega = Scalar::get_root_of_unity(polys[0].len().trailing_zeros());

    let mut omega_power = vec![Scalar::ONE; polys[0].len()];
    for idx in 1..omega_power.len() {
        omega_power[idx] = omega * omega_power[idx - 1]
    }
    let evaluations = (0..polys.len())
        .map(|idx| evaluations_at_z[idx].evaluation_of_f)
        .collect::<Vec<Scalar>>();

    let phi_k_evaluations = (0..num_var - 1)
        .map(|idx| compute_evaluations(1 << (idx + 1)))
        .collect::<Vec<Vec<Scalar>>>();

    let even_poly_evaluations = &(0..polys.len())
        .map(|idx| {
            evaluations_of_even_odd_polys[idx]
                .clone()
                .even_polys_evaluations
        })
        .collect::<Vec<Vec<Vec<Scalar>>>>();

    let odd_poly_evaluations = &(0..polys.len())
        .map(|idx| {
            evaluations_of_even_odd_polys[idx]
                .clone()
                .odd_polys_evaluations
        })
        .collect::<Vec<Vec<Vec<Scalar>>>>();

    let even_poly_evaluation_at_z = &(0..polys.len())
        .map(|idx| evaluations_at_z[idx].clone().even_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();
    let odd_poly_evaluation_at_z = &(0..polys.len())
        .map(|idx| evaluations_at_z[idx].clone().odd_polys_evaluations)
        .collect::<Vec<Vec<Scalar>>>();

    let mut extended_phi_k_polys_evals = Vec::new();
    (0..num_var as usize - 1).for_each(|idx| {
        let mut phi_poly = phi_k_evaluations[idx].clone();
        for _ in 0..(1 << num_var) / phi_k_evaluations[idx].len() - 1 {
            phi_poly = [phi_poly, phi_k_evaluations[idx].clone()].concat()
        }
        assert_eq!(
            phi_poly.len(),
            (1 << num_var),
            "length mismatch for phi k polys"
        );
        extended_phi_k_polys_evals.push(phi_poly);
    });

    // Compute z^(1<<(d - (k +1)))
    let z_power_d_k = (0..num_var - 1)
        .into_par_iter()
        .map(|idx| z.power_by([(1 << (num_var - (idx + 1))) as u64, 0, 0, 0]))
        .collect::<Vec<Scalar>>();

    //Compute (z^(1<<(d - (k +1))))^(1<<d - 1<<k)
    let alpha = (0..num_var - 1)
        .into_par_iter()
        .map(|idx| {
            z_power_d_k[idx as usize].power_by([((1 << num_var) - (1 << idx)) as u64, 0, 0, 0])
        })
        .collect::<Vec<Scalar>>();

    let shi = (0..polys[0].len())
        .into_par_iter()
        .map(|outer_idx| {
            (0..num_var as usize - 1).fold(Scalar::ZERO, |acc1, inner_idx| {
                let omega_power_d_k = omega_power[outer_idx].power_by([
                    ((1 << num_var) - (1 << inner_idx)) as u64,
                    0,
                    0,
                    0,
                ]);

                acc1 + (random_points[polys.len() + inner_idx]
                    * (extended_phi_k_polys_evals[inner_idx][outer_idx]
                        - phi_k_evaluation_at_z[inner_idx])
                    * ((omega_power[outer_idx].power_by([
                        (1 << (num_var as usize - (inner_idx + 1) - 1)) as u64,
                        0,
                        0,
                        0,
                    ]) - z.power_by([
                        (1 << (num_var as usize - (inner_idx + 1) - 1)) as u64,
                        0,
                        0,
                        0,
                    ]))
                    .invert()
                    .unwrap()))
                    + ((0..polys.len()).fold(Scalar::ZERO, |acc, idx| {
                        acc + (random_points
                            [4 * inner_idx + polys.len() + num_var as usize - 1 + idx]
                            * (even_poly_evaluations[idx][inner_idx][outer_idx]
                                - even_poly_evaluation_at_z[idx][inner_idx]))
                            + (random_points[4 * inner_idx + polys.len() + num_var as usize + idx]
                                * (odd_poly_evaluations[idx][inner_idx][outer_idx]
                                    - odd_poly_evaluation_at_z[idx][inner_idx]))
                            + (random_points
                                [4 * inner_idx + polys.len() + num_var as usize + 1 + idx]
                                * (omega_power_d_k
                                    * even_poly_evaluations[idx][inner_idx][outer_idx]
                                    - even_poly_evaluation_at_z[idx][inner_idx] * alpha[inner_idx]))
                            + (random_points
                                [4 * inner_idx + polys.len() + num_var as usize + 2 + idx]
                                * (omega_power_d_k
                                    * odd_poly_evaluations[idx][inner_idx][outer_idx]
                                    - odd_poly_evaluation_at_z[idx][inner_idx] * alpha[inner_idx]))
                    }) * (omega_power[outer_idx] - z_power_d_k[inner_idx])
                        .invert()
                        .unwrap())
            }) + ((0..polys.len()).fold(Scalar::ZERO, |acc, idx| {
                acc + random_points[idx] * (polys[idx][outer_idx] - evaluations[idx])
            }) * (omega_power[outer_idx] - *z).invert().unwrap())
        })
        .collect::<Vec<Scalar>>();
    (shi, extended_phi_k_polys_evals)
}
