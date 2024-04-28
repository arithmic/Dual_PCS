use crate::{multilinear_common::reseed_with_even_odd_polys, setup::compute_evaluations};
use bls381::scalar::Scalar;
use channel::Channel;
use common::{
    CommitmentsOfEvenOddPolys, EvaluationOfEvenOddPolys, EvaluationsAtRandomPoint, EvenOddPolys,
    KZGFourierDegreeBoundSetup, MultilinearKZG2Proof, OFFSET,
};

use fft::{
    par_fft::{
        get_inv_twiddles, get_twiddles, par_eval_poly_with_offset, par_interpolate_poly,
        par_interpolate_poly_with_offset,
    },
    serial_fft::{eval, log2},
};
use helper::helper::{compute_dot_product, compute_fourier_bases};
use kzg_fft::commit;
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
    IntoParallelRefMutIterator, ParallelIterator,
};
use traits::traits::{Field, PrimeField};
pub fn prover(
    poly: &Vec<Scalar>,
    random_points: Vec<Scalar>,
    setup: &KZGFourierDegreeBoundSetup,
    channel: &mut Channel,
) -> MultilinearKZG2Proof {
    let num_var = poly.len().trailing_zeros();
    //Compute bases
    let evaluation = evaluate_multilinear_poly(poly, &random_points);

    //Compute q_k s.t. f - f(u) = sum(X_k - u_k) * q_k , k = {0, 1, .., n-1}
    let quotient_polys = compute_quotient_polys(poly, random_points, evaluation);

    //Evaluate g_{q_k,e}^{k} and g_{q_k,o}^{k}
    let (even_odd_polys, evaluations_of_even_odd_polys) = compute_even_odd_polys(quotient_polys);

    //Commit evaluations of even and odd polys
    let commitments_of_even_ood_polys =
        commit_even_odd_polys(&evaluations_of_even_odd_polys, &setup);

    //Reseed channel with commitments
    reseed_with_even_odd_polys(channel, commitments_of_even_ood_polys.clone());

    //Draw a random point using Fiat-shamir
    let z = channel.get_random_point();

    //Evaluate all the polys at r
    let q_not_poly = even_odd_polys.const_poly.clone();
    let (evaluations_at_z, interpolated_poly) =
        evaluations_at_random_point(z, &even_odd_polys, poly.clone());

    let phi_k_evaluation_at_z = evaluate_phi_at_random_point(num_var as usize, z);

    //Draw random points
    let r_points = channel.get_random_points((5 * num_var - 4) as usize);

    let (shi, extended_phi_k_polys_evals) = compute_shi(
        poly,
        &evaluations_at_z,
        &r_points,
        &z,
        &evaluations_of_even_odd_polys,
        &phi_k_evaluation_at_z,
    );

    //Commit shi
    let shi_commit = commit::kzg2commit(&shi, &setup.setup.get_setup(shi.len()).prover_key);
    channel.reseed_with_affine_point(&[shi_commit].to_vec());

    //Draw a random point using Fiat-shamir
    let s = channel.get_random_point();

    //Draw random points to combine all the polys using Fiat-Shamir
    let poly_combiners = channel.get_random_points((num_var * 3 - 2) as usize);

    //Evaluate all the polys at r
    let evaluations_at_s = evaluate_polynomials_at_s(s, &even_odd_polys, interpolated_poly);
    let phi_k_evaluation_at_s = evaluate_phi_at_random_point(num_var as usize, s);

    let r = combine_polynomials(
        shi,
        poly,
        &evaluations_of_even_odd_polys,
        &poly_combiners,
        s,
        &evaluations_at_s,
        extended_phi_k_polys_evals,
        &phi_k_evaluation_at_s,
    );

    let commitment_to_r = commit::kzg2commit(&r, &setup.setup.get_setup(r.len()).prover_key);
    MultilinearKZG2Proof::new(
        evaluation,
        commitments_of_even_ood_polys,
        evaluations_at_z,
        commitment_to_r,
        evaluations_at_s,
        shi_commit,
        phi_k_evaluation_at_z,
        phi_k_evaluation_at_s,
        q_not_poly,
    )
}
pub fn evaluate_multilinear_poly(poly: &Vec<Scalar>, random_points: &Vec<Scalar>) -> Scalar {
    let fourier_bases = compute_fourier_bases(random_points);
    compute_dot_product(poly, &fourier_bases)
}

pub(crate) fn compute_quotient_polys(
    poly: &Vec<Scalar>,
    point: Vec<Scalar>,
    eval: Scalar,
) -> Vec<Vec<Scalar>> {
    let n_vars = poly.len().trailing_zeros() as usize;
    let mut rec_step: Vec<Scalar> = poly.par_iter().map(|coeff| *coeff - eval).collect();
    let mut witnesses: Vec<Vec<Scalar>> = vec![Vec::<Scalar>::new(); n_vars];

    //Computing witness polynomials
    for variable in 0..n_vars {
        let halfsize = rec_step.len() >> 1;
        let mut r = vec![Scalar::ZERO; halfsize];
        let mut witness_i = vec![Scalar::ZERO; halfsize];
        let current_point = point[variable];

        let point_inv = (-current_point).invert().unwrap();

        r.par_iter_mut().enumerate().for_each(|(i, r_i)| {
            *r_i = rec_step[i] + (rec_step[i + halfsize] - rec_step[i]) * current_point;
        });

        witness_i
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, witness_variable_i)| {
                *witness_variable_i = (rec_step[i] - r[i]) * point_inv;
            });

        witnesses[variable] = witness_i;
        rec_step = r;
    }
    witnesses.reverse();
    witnesses
}
pub(crate) fn compute_even_odd_polys(
    quotient_polys: Vec<Vec<Scalar>>,
) -> (EvenOddPolys, EvaluationOfEvenOddPolys) {
    let num_vars = quotient_polys.len();

    //Shift the polynomials by 2^k
    let q_not = quotient_polys[0].clone();
    let shifted_quotient_polys = (1..num_vars)
        .into_par_iter()
        .map(|k| {
            let length = quotient_polys[k].len();
            [vec![Scalar::ZERO; length], quotient_polys[k].clone()].concat()
        })
        .collect::<Vec<Vec<Scalar>>>();

    //Interpolate Polynomials over the domain 2^k+1
    let shifted_quotient_polys = (0..shifted_quotient_polys.len())
        .into_par_iter()
        .map(|idx| {
            let inv_twiddles = get_inv_twiddles(shifted_quotient_polys[idx].len() as u32);
            let mut poly = shifted_quotient_polys[idx].clone();
            par_interpolate_poly(&mut poly, inv_twiddles);
            poly
        })
        .collect::<Vec<Vec<Scalar>>>();

    let shifted_polys_evaluations = (0..shifted_quotient_polys.len())
        .into_par_iter()
        .map(|outer_idx| {
            let domain = shifted_quotient_polys[outer_idx].len();
            let half_domain = domain / 2;
            let g_half_domain = Scalar::get_root_of_unity(log2(half_domain));
            let g_domain = Scalar::get_root_of_unity(log2(domain));
            let mut offset = vec![OFFSET; half_domain];
            for idx in 1..half_domain {
                offset[idx] = offset[idx - 1] * g_half_domain
            }
            let evaluations = (0..half_domain)
                .into_par_iter()
                .map(|iter| {
                    let mut denominator = Scalar::ONE;
                    let mut g_domain_power = Scalar::ONE;
                    for _ in 0..half_domain {
                        denominator *= offset[iter] - g_domain_power;
                        g_domain_power *= g_domain
                    }

                    eval(&shifted_quotient_polys[outer_idx], offset[iter])
                        * (denominator.invert().unwrap())
                })
                .collect::<Vec<Scalar>>();
            evaluations
        })
        .collect::<Vec<Vec<Scalar>>>();

    //Interpolate polynomials.
    let shifted_polys_evaluations = (0..shifted_polys_evaluations.len())
        .into_par_iter()
        .map(|idx| {
            let inv_twiddles = get_inv_twiddles(shifted_polys_evaluations[idx].len() as u32);
            let mut evaluations = shifted_polys_evaluations[idx].clone();
            par_interpolate_poly_with_offset(&mut evaluations, inv_twiddles, OFFSET);
            evaluations
        })
        .collect::<Vec<Vec<Scalar>>>();

    //Break the polynomial into even polynomial and odd polynomial i.e. f(X) = f_e(X^2) + Xf_o(X^2)
    let (even_polys, odd_polys) = (0..num_vars - 1)
        .into_par_iter()
        .map(|outer_idx| {
            let mut even_poly = Vec::new();
            let mut odd_poly = Vec::new();
            (0..shifted_polys_evaluations[outer_idx].len()).for_each(|inner_idx| {
                if inner_idx % 2 == 0 {
                    even_poly.push(shifted_polys_evaluations[outer_idx][inner_idx])
                } else {
                    odd_poly.push(shifted_polys_evaluations[outer_idx][inner_idx])
                }
            });
            (even_poly, odd_poly)
        })
        .collect::<(Vec<Vec<Scalar>>, Vec<Vec<Scalar>>)>();

    let even_odd_polys = EvenOddPolys {
        const_poly: q_not.clone(),
        even_polys: even_polys.clone(),
        odd_polys: odd_polys.clone(),
    };

    // Evaluate even and odd poly over the domain D
    let (even_polys_evaluations, odd_polys_evaluations) = (0..num_vars - 1)
        .into_par_iter()
        .map(|outer_idx| {
            let offset = (1 << num_vars) / even_polys[outer_idx].len();
            let twiddles = get_twiddles(even_polys[outer_idx].len() as u32);
            (
                par_eval_poly_with_offset(&even_polys[outer_idx], &twiddles, Scalar::ONE, offset),
                par_eval_poly_with_offset(&odd_polys[outer_idx], &twiddles, Scalar::ONE, offset),
            )
        })
        .collect::<(Vec<Vec<Scalar>>, Vec<Vec<Scalar>>)>();

    (
        even_odd_polys,
        EvaluationOfEvenOddPolys {
            const_poly: q_not,
            even_polys_evaluations,
            odd_polys_evaluations,
        },
    )
}

pub(crate) fn commit_even_odd_polys(
    evaluations: &EvaluationOfEvenOddPolys,
    setup: &KZGFourierDegreeBoundSetup,
) -> CommitmentsOfEvenOddPolys {
    assert_eq!(
        evaluations.even_polys_evaluations.len(),
        evaluations.odd_polys_evaluations.len(),
        "length of even polys should be equal to length of odd polys"
    );
    let len = evaluations.even_polys_evaluations.len();

    let (even_polys_commits, odd_polys_commits) = (0..len)
        .into_par_iter()
        .map(|idx| {
            let prover_key = &setup
                .setup
                .get_setup(evaluations.even_polys_evaluations[idx].len())
                .prover_key;
            let commit1 = commit::kzg2commit(&evaluations.even_polys_evaluations[idx], prover_key);
            let commit2 = commit::kzg2commit(&evaluations.odd_polys_evaluations[idx], prover_key);
            (commit1, commit2)
        })
        .unzip();

    CommitmentsOfEvenOddPolys::new(even_polys_commits, odd_polys_commits)
}
pub(crate) fn evaluate_phi_at_random_point(num_var: usize, random_point: Scalar) -> Vec<Scalar> {
    let phi_k_evaluations = (0..num_var - 1)
        .into_par_iter()
        .map(|idx| {
            let k = idx + 1;
            let mut phi_k_negative = Scalar::ONE;
            let omega = Scalar::get_root_of_unity((k + 1) as u32);
            let eval_point = random_point.power_by([(1 << (num_var - k - 1)) as u64, 0, 0, 0]);
            let mut omega_power = Scalar::ONE;
            for _ in 0..(1 << k) {
                phi_k_negative *= eval_point - omega_power;
                omega_power *= omega;
            }
            phi_k_negative
        })
        .collect::<Vec<Scalar>>();
    phi_k_evaluations
}
pub(crate) fn evaluations_at_random_point(
    z: Scalar,
    even_odd_polys: &EvenOddPolys,
    mut poly: Vec<Scalar>,
) -> (EvaluationsAtRandomPoint, Vec<Scalar>) {
    let num_var = poly.len().trailing_zeros() as usize;
    assert_eq!(
        even_odd_polys.even_polys.len(),
        even_odd_polys.odd_polys.len(),
        "length should be equal"
    );

    let (even_polys_evaluations, odd_polys_evaluations) = (0..even_odd_polys.even_polys.len())
        .into_par_iter()
        .map(|idx| {
            let eval_point = z.power_by([(1 << (num_var - (idx + 1))) as u64, 0, 0, 0]);
            let even_poly_eval = eval(&even_odd_polys.even_polys[idx], eval_point);
            let odd_poly_eval = eval(&even_odd_polys.odd_polys[idx], eval_point);
            (even_poly_eval, odd_poly_eval)
        })
        .collect::<(Vec<Scalar>, Vec<Scalar>)>();

    let inv_twiddles = get_inv_twiddles(poly.len() as u32);
    par_interpolate_poly(&mut poly, inv_twiddles);

    let evaluation_of_f = eval(&poly, z);

    (
        EvaluationsAtRandomPoint {
            evaluation_of_f,
            even_polys_evaluations,
            odd_polys_evaluations,
            // phi_k_evaluations,
            // q_zero,
        },
        poly,
    )
}

pub(crate) fn evaluate_polynomials_at_s(
    s: Scalar,
    even_odd_polys: &EvenOddPolys,
    interpolated_poly: Vec<Scalar>,
) -> EvaluationsAtRandomPoint {
    let even_polys = &even_odd_polys.even_polys;
    let odd_polys = &even_odd_polys.odd_polys;

    assert_eq!(even_polys.len(), odd_polys.len(), "length should be equal");

    let (even_polys_evaluations, odd_polys_evaluations) = (0..odd_polys.len())
        .into_par_iter()
        .map(|idx| {
            let even_poly_eval = eval(&even_polys[idx], s);
            let odd_poly_eval = eval(&odd_polys[idx], s);
            (even_poly_eval, odd_poly_eval)
        })
        .collect::<(Vec<Scalar>, Vec<Scalar>)>();

    let evaluation_of_f = eval(&interpolated_poly, s);

    EvaluationsAtRandomPoint {
        evaluation_of_f,
        even_polys_evaluations,
        odd_polys_evaluations,
    }
}
fn compute_shi(
    poly: &Vec<Scalar>,
    evaluations_at_z: &EvaluationsAtRandomPoint,
    random_points: &Vec<Scalar>,
    z: &Scalar,
    evaluations_of_even_odd_polys: &EvaluationOfEvenOddPolys,
    phi_k_evaluation_at_z: &Vec<Scalar>,
) -> (Vec<Scalar>, Vec<Vec<Scalar>>) {
    let num_var = poly.len().trailing_zeros();
    let omega = Scalar::get_root_of_unity(poly.len().trailing_zeros());
    let mut omega_power = vec![Scalar::ONE; poly.len()];
    for idx in 1..omega_power.len() {
        omega_power[idx] = omega * omega_power[idx - 1]
    }
    let evaluation = evaluations_at_z.evaluation_of_f;

    let phi_k_evaluations = (0..num_var - 1)
        .map(|idx| compute_evaluations(1 << (idx + 1)))
        .collect::<Vec<Vec<Scalar>>>();

    let even_poly_evaluations = &evaluations_of_even_odd_polys.even_polys_evaluations;
    let odd_poly_evaluations = &evaluations_of_even_odd_polys.odd_polys_evaluations;
    let even_poly_evaluation_at_z = &evaluations_at_z.even_polys_evaluations;
    let odd_poly_evaluation_at_z = &evaluations_at_z.odd_polys_evaluations;

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

    let shi = (0..poly.len())
        .into_par_iter()
        .map(|outer_idx| {
            (0..num_var as usize - 1).fold(Scalar::ZERO, |acc, inner_idx| {
                let omega_power_d_k = omega_power[outer_idx].power_by([
                    ((1 << num_var) - (1 << inner_idx)) as u64,
                    0,
                    0,
                    0,
                ]);

                acc + random_points[5 * inner_idx + 1]
                    * (extended_phi_k_polys_evals[inner_idx][outer_idx]
                        - phi_k_evaluation_at_z[inner_idx])
                    * (omega_power[outer_idx].power_by([
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
                    .unwrap()
                    + ((random_points[5 * inner_idx + 2]
                        * (even_poly_evaluations[inner_idx][outer_idx]
                            - even_poly_evaluation_at_z[inner_idx]))
                        + (random_points[5 * inner_idx + 3]
                            * (odd_poly_evaluations[inner_idx][outer_idx]
                                - odd_poly_evaluation_at_z[inner_idx]))
                        + (random_points[5 * inner_idx + 4]
                            * (omega_power_d_k * even_poly_evaluations[inner_idx][outer_idx]
                                - even_poly_evaluation_at_z[inner_idx] * alpha[inner_idx]))
                        + (random_points[5 * inner_idx + 5]
                            * (omega_power_d_k * odd_poly_evaluations[inner_idx][outer_idx]
                                - odd_poly_evaluation_at_z[inner_idx] * alpha[inner_idx])))
                        * (omega_power[outer_idx] - z_power_d_k[inner_idx])
                            .invert()
                            .unwrap()
            }) + (random_points[0]
                * ((poly[outer_idx] - evaluation)
                    * (omega_power[outer_idx] - *z).invert().unwrap()))
        })
        .collect::<Vec<Scalar>>();
    (shi, extended_phi_k_polys_evals)
}

fn combine_polynomials(
    mut shi: Vec<Scalar>,
    poly: &Vec<Scalar>,
    even_odd_polys_evals: &EvaluationOfEvenOddPolys,
    random_points: &[Scalar],
    s: Scalar,
    evaluations_at_s: &EvaluationsAtRandomPoint,
    extended_phi_k_polys_evals: Vec<Vec<Scalar>>,
    phi_k_evaluation_at_s: &Vec<Scalar>,
) -> Vec<Scalar> {
    let num_var = poly.len().trailing_zeros();
    let even_poly_evaluations = &even_odd_polys_evals.even_polys_evaluations;
    let odd_poly_evaluations = &even_odd_polys_evals.odd_polys_evaluations;
    assert_eq!(
        even_poly_evaluations.len(),
        odd_poly_evaluations.len(),
        "no of even polys and no of odd polys should be equal"
    );
    let g_num_var = Scalar::get_root_of_unity(num_var);

    let evaluations_of_g = (0..1 << num_var)
        .into_par_iter()
        .map(|idx| {
            let mut temp = Scalar::ZERO;
            for inner_idx in 0..even_odd_polys_evals.odd_polys_evaluations.len() {
                temp += (random_points[3 * inner_idx + 1] * even_poly_evaluations[inner_idx][idx])
                    + (random_points[3 * inner_idx + 2] * odd_poly_evaluations[inner_idx][idx])
                    + (random_points[3 * inner_idx + 3]
                        * extended_phi_k_polys_evals[inner_idx][idx])
            }
            temp + (random_points[0] * poly[idx]) + shi[idx]
        })
        .collect::<Vec<Scalar>>();

    let even_polys_evaluations = &evaluations_at_s.even_polys_evaluations;
    let odd_polys_evaluations = &evaluations_at_s.odd_polys_evaluations;

    let f_s = evaluations_at_s.evaluation_of_f;

    let inv_twiddles = get_inv_twiddles(shi.len() as u32);
    par_interpolate_poly(&mut shi, inv_twiddles);
    let shi_at_s = eval(&shi, s);

    let g_at_s = (0..even_polys_evaluations.len())
        .into_par_iter()
        .map(|iter| {
            (random_points[3 * iter + 1] * even_polys_evaluations[iter])
                + (random_points[3 * iter + 2] * odd_polys_evaluations[iter])
                + (random_points[3 * iter + 3] * phi_k_evaluation_at_s[iter])
        })
        .collect::<Vec<Scalar>>()
        .into_iter()
        .fold(Scalar::ZERO, |acc, value| acc + value)
        + (random_points[0] * f_s)
        + shi_at_s;

    let mut omega_powers = vec![Scalar::ONE; 1 << num_var];
    for idx in 1..omega_powers.len() {
        let temp = omega_powers[idx - 1] * g_num_var;
        omega_powers[idx] = temp;
    }
    let quotient_poly = (0..1 << num_var)
        .into_par_iter()
        .map(|idx| ((evaluations_of_g[idx] - g_at_s) * (omega_powers[idx] - s).invert().unwrap()))
        .collect::<Vec<Scalar>>();

    quotient_poly
}

pub fn batch_eval(
    polys: &Vec<Vec<Scalar>>,
    combiners: &Vec<Scalar>,
    evals: Vec<Scalar>,
    random_points: Vec<Scalar>,
    setup: &KZGFourierDegreeBoundSetup,
    channel: &mut Channel,
) -> MultilinearKZG2Proof {
    let n_polys = polys.len();
    let combined_poly = (0..polys[0].len())
        .map(|outer_idx| {
            (0..n_polys).fold(Scalar::ZERO, |acc, inner_idx| {
                acc + (combiners[inner_idx] * polys[inner_idx][outer_idx])
            })
        })
        .collect::<Vec<Scalar>>();
    let combined_eval =
        (0..evals.len()).fold(Scalar::ZERO, |acc, idx| acc + (combiners[idx] * evals[idx]));
    let proof = prover(&combined_poly, random_points, setup, channel);
    assert_eq!(proof.get_evaluation(), combined_eval);
    proof
}
