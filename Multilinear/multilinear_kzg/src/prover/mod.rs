use curve_traits::{self, projectivepoint::Scalar, CurveArithmetic, CurveParams, ProjectivePoint};
use pairing::pairing_traits::Pairing;
use rayon::{
    self,
    prelude::{
        IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
        IntoParallelRefMutIterator, ParallelIterator,
    },
};
use std::vec;
use traits::traits::Field;

use crate::common::{MleCommit, MleEvalProof, SRS};

pub fn commit<C: Pairing>(poly: &Vec<Scalar<C::G1Curve>>, srs: &SRS<C>) -> MleCommit<C> {
    let mut commit: ProjectivePoint<C::G1Curve>;

    if poly.len() < srs.g_1.len() {
        let halfsize = srs.g_1.len() >> 1;
        let mut srs_g1: Vec<_> = (0..halfsize)
            .into_par_iter()
            .map(|j| srs.g_1[j] + srs.g_1[j + halfsize])
            .collect();

        while srs_g1.len() > poly.len() {
            let halfsize = srs_g1.len() >> 1;
            srs_g1 = (0..halfsize)
                .into_par_iter()
                .map(|j| srs_g1[j] + srs_g1[j + halfsize])
                .collect();
        }
        commit = ProjectivePoint::pippenger_msm(&srs_g1, &poly);
    } else {
        commit = ProjectivePoint::pippenger_msm(&srs.g_1, &poly);
    }

    MleCommit { commitment: commit }
}

pub fn evaluate<C: Pairing>(
    poly: &Vec<Scalar<C::G1Curve>>,
    point: &Vec<Scalar<C::G1Curve>>,
    srs: &SRS<C>,
) -> MleEvalProof<C> {
    let n_vars = poly.len().trailing_zeros() as usize;

    debug_assert_eq!(n_vars, point.len());

    let mut srs_g1: Vec<ProjectivePoint<C::G1Curve>> =
        vec![ProjectivePoint::IDENTITY; srs.g_1.len()];
    if poly.len() < srs.g_1.len() {
        let halfsize = srs.g_1.len() >> 1;
        srs_g1 = (0..halfsize)
            .into_par_iter()
            .map(|j| srs.g_1[j] + srs.g_1[j + halfsize])
            .collect();
        while srs_g1.len() > poly.len() {
            let halfsize = srs_g1.len() >> 1;
            srs_g1 = (0..halfsize)
                .into_par_iter()
                .map(|j| srs_g1[j] + srs_g1[j + halfsize])
                .collect();
        }
    } else {
        srs_g1.copy_from_slice(&srs.g_1);
    }

    let mut witnesses: Vec<Vec<Scalar<C::G1Curve>>> =
        vec![Vec::<Scalar<C::G1Curve>>::new(); n_vars];

    let bases_eval: Vec<Scalar<C::G1Curve>> = compute_fourier_bases::<C::G1Curve>(&point);

    let eval = bases_eval
        .par_iter()
        .zip(poly.par_iter())
        .fold(
            || <Scalar<C::G1Curve> as Field>::ZERO,
            |acc, (basis, coefficient)| acc + (*basis * *coefficient),
        )
        .reduce(
            || <Scalar<C::G1Curve> as Field>::ZERO,
            |acc, part_sum| acc + part_sum,
        );

    let mut rec_step: Vec<Scalar<C::G1Curve>> =
        poly.par_iter().map(|coeff| *coeff - eval).collect();

    //Computing witness polynomials
    for variable in 0..n_vars {
        let halfsize = rec_step.len() >> 1;
        let mut r = vec![<Scalar<C::G1Curve> as Field>::ZERO; halfsize];
        let mut witness_i = vec![<Scalar<C::G1Curve> as Field>::ZERO; halfsize];
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

    let mut witness_commits = vec![ProjectivePoint::IDENTITY; n_vars];
    let halfsize = poly.len() >> 1;
    let mut batching_srs: Vec<ProjectivePoint<C::G1Curve>> = (0..halfsize)
        .into_par_iter()
        .map(|i| srs_g1[i] + srs_g1[i + halfsize])
        .collect();

    //Computing commit
    for i in 0..n_vars {
        witness_commits[i] =
            ProjectivePoint::pippenger_msm(&batching_srs, &witnesses[i][..batching_srs.len()]);

        let next_length = batching_srs.len() >> 1;
        batching_srs = (0..next_length)
            .into_par_iter()
            .map(|i| batching_srs[i] + batching_srs[i + next_length])
            .collect();
    }

    MleEvalProof {
        evaluation: eval,
        witnesses: witness_commits,
    }
}

pub fn batch_eval<C: Pairing>(
    polys: &Vec<&Vec<Scalar<C::G1Curve>>>,
    evals: &Vec<Scalar<C::G1Curve>>,
    point: &Vec<Scalar<C::G1Curve>>,
    scalars: &Vec<Scalar<C::G1Curve>>,
    srs: &SRS<C>,
) -> MleEvalProof<C> {
    let mut poly: Vec<Scalar<C::G1Curve>> = vec![Scalar::<C::G1Curve>::ZERO; 1 << point.len()];
    let n_polys = polys.len();

    for i in 0..poly.len() {
        for j in 0..n_polys {
            poly[i] += polys[j][i] * scalars[j]
        }
    }

    let linear_combination_eval = evals
        .iter()
        .zip(scalars.iter())
        .map(|(eval, scalar)| *eval * *scalar)
        .reduce(|acc, g| acc + g)
        .unwrap();

    let n_vars = poly.len().trailing_zeros() as usize;

    debug_assert_eq!(n_vars, point.len());

    let mut srs_g1: Vec<ProjectivePoint<C::G1Curve>> =
        vec![ProjectivePoint::IDENTITY; srs.g_1.len()];
    if poly.len() < srs.g_1.len() {
        let halfsize = srs.g_1.len() >> 1;
        srs_g1 = (0..halfsize)
            .into_iter()
            .map(|j| srs.g_1[j] + srs.g_1[j + halfsize])
            .collect();
        while srs_g1.len() > poly.len() {
            let halfsize = srs_g1.len() >> 1;
            srs_g1 = (0..halfsize)
                .into_iter()
                .map(|j| srs_g1[j] + srs_g1[j + halfsize])
                .collect();
        }
    } else {
        srs_g1.copy_from_slice(&srs.g_1);
    }
    let mut witnesses: Vec<Vec<Scalar<C::G1Curve>>> =
        vec![Vec::<Scalar<C::G1Curve>>::new(); n_vars];

    let mut rec_step: Vec<Scalar<C::G1Curve>> = poly
        .iter()
        .map(|coeff| *coeff - linear_combination_eval)
        .collect();

    //Computing witness polynomials
    for variable in 0..n_vars {
        let halfsize = rec_step.len() >> 1;
        let mut r = vec![<Scalar<C::G1Curve> as Field>::ZERO; halfsize];
        let mut witness_i = vec![<Scalar<C::G1Curve> as Field>::ZERO; halfsize];
        let current_point = point[variable];

        let point_inv = (-current_point).invert().unwrap();

        r.iter_mut().enumerate().for_each(|(i, r_i)| {
            *r_i = rec_step[i] + (rec_step[i + halfsize] - rec_step[i]) * current_point;
        });

        witness_i
            .iter_mut()
            .enumerate()
            .for_each(|(i, witness_variable_i)| {
                *witness_variable_i = (rec_step[i] - r[i]) * point_inv;
            });

        witnesses[variable] = witness_i;
        rec_step = r;
    }

    let mut witness_commits = vec![ProjectivePoint::IDENTITY; n_vars];
    let halfsize = poly.len() >> 1;
    let mut batching_srs: Vec<ProjectivePoint<C::G1Curve>> = (0..halfsize)
        .into_iter()
        .map(|i| srs_g1[i] + srs_g1[i + halfsize])
        .collect();

    //Computing commit
    for i in 0..n_vars {
        witness_commits[i] =
            ProjectivePoint::pippenger_msm(&batching_srs, &witnesses[i][..batching_srs.len()]);

        let next_length = batching_srs.len() >> 1;
        batching_srs = (0..next_length)
            .into_iter()
            .map(|i| batching_srs[i] + batching_srs[i + next_length])
            .collect();
    }

    MleEvalProof {
        evaluation: linear_combination_eval,
        witnesses: witness_commits,
    }
}

pub fn compute_fourier_bases<C: CurveArithmetic + CurveParams>(
    r: &Vec<Scalar<C>>,
) -> Vec<Scalar<C>> {
    //Initialize fc_eq with (1- r[0]) and r[0]
    let mut fc_eq = [<Scalar<C> as Field>::ONE - r[0], r[0]].to_vec();
    //Iterate over the length of the r vector
    for k in 1..r.len() {
        let temp = fc_eq;
        //initialize fc_eq of double size with zero
        fc_eq = vec![<Scalar<C> as Field>::ZERO; temp.len() * 2];
        for iter in 0..temp.len() {
            fc_eq[2 * iter] = temp[iter] * (<Scalar<C> as Field>::ONE - r[k]);
            fc_eq[(2 * iter) + 1] = temp[iter] * (r[k]);
        }
    }
    fc_eq
}
