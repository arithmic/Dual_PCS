#![allow(unused)]
use crate::{
    preprocessing,
    prover::batch_eval::{batch_eval_proof, len_4_interpolate},
    spartan_common,
};
use bls381::scalar::Scalar;
use channel::Channel;
use kzg_fft::commit::kzg2commit;
use kzg_fourier_multilinear::{
    batch_prover::batch_eval_multiple_evaluation_points,
    common::{BatchMultilinearKZG2Proof, KZGFourierDegreeBoundSetup},
    prover::evaluate_multilinear_poly,
};
use polynomial::{MultPolynomial, Polynomial};
use preprocessing::{compute_coeff, SparseRep};
use spartan_common::{
    eR1CSmetadata, eR1CStranscript, interpolate, InitialSumCheckTranscript, ParSumCheckTranscript,
};
use traits::traits::Field;
#[allow(non_snake_case)]
pub fn prove_sat(
    A: SparseRep,
    B: SparseRep,
    C: SparseRep,
    u: &Scalar,
    z: &MultPolynomial,
    E: &MultPolynomial,
    W: &MultPolynomial,
    metadatas: eR1CSmetadata,
    channel: &mut Channel,
    setup: &KZGFourierDegreeBoundSetup,
) -> eR1CStranscript {
    let first_sum_check_transcript = first_layer_sum_check(&A, &B, &C, u, z, E, channel);

    let rx = first_sum_check_transcript.random_points.clone();

    let rx_basis_evals = compute_coeff(&first_sum_check_transcript.random_points);

    let par_sum_check_transcript = parallel_sum_checks(&A, &B, &C, z, rx_basis_evals, channel);

    let ry = par_sum_check_transcript.random_points.clone();
    let E_evals = evaluate_multilinear_poly(E.as_coeffs(), &rx);
    let W_evals = evaluate_multilinear_poly(W.as_coeffs(), &ry);

    let rx_basis_evals = compute_coeff(&rx.to_vec());
    let ry_basis_evals = compute_coeff(&ry.to_vec());
    let rx_basis_commit = kzg2commit(
        &rx_basis_evals,
        &setup.get_setup(rx_basis_evals.len()).setup.prover_key,
    );
    let ry_basis_commit = kzg2commit(
        &ry_basis_evals,
        &setup.get_setup(ry_basis_evals.len()).setup.prover_key,
    );
    let random_points = channel.get_random_points(rx.len());
    let e_w_eval_proof = batch_multilinear_polys(
        E.as_coeffs().to_vec(),
        W.as_coeffs().to_vec(),
        E_evals,
        W_evals,
        rx,
        ry,
        rx_basis_evals.clone(),
        ry_basis_evals.clone(),
        random_points,
        setup,
        channel,
    );
    let metadatas = vec![metadatas.A, metadatas.B, metadatas.C];
    let eval_proofs = batch_eval_proof(metadatas, channel, setup, rx_basis_evals, ry_basis_evals);

    eR1CStranscript::new(
        first_sum_check_transcript,
        par_sum_check_transcript,
        eval_proofs,
        e_w_eval_proof,
        rx_basis_commit,
        ry_basis_commit,
    )
}
#[allow(non_snake_case)]
fn batch_multilinear_polys(
    E: Vec<Scalar>,
    W: Vec<Scalar>,
    E_evals: Scalar,
    W_evals: Scalar,
    rx: Vec<Scalar>,
    ry: Vec<Scalar>,
    rx_basis_evals: Vec<Scalar>,
    ry_basis_evals: Vec<Scalar>,
    random_points: Vec<Scalar>,
    setup: &KZGFourierDegreeBoundSetup,
    channel: &mut Channel,
) -> BatchMultilinearKZG2Proof {
    let rx_basis_evaluation = evaluate_multilinear_poly(&rx_basis_evals, &random_points);
    let ry_basis_evaluation = evaluate_multilinear_poly(&ry_basis_evals, &random_points);
    let combiner = channel.get_random_points(2);
    let combined_eval = combiner[0] * rx_basis_evaluation + combiner[1] * ry_basis_evaluation;

    let mut combined_poly = vec![Scalar::ZERO; rx_basis_evals.len()];
    for idx in 0..rx_basis_evals.len() {
        combined_poly[idx] = combiner[0] * rx_basis_evals[idx] + combiner[1] * ry_basis_evals[idx]
    }
    let e_w_eval_proof = batch_eval_multiple_evaluation_points(
        &vec![E, W, combined_poly],
        vec![E_evals, W_evals, combined_eval],
        vec![rx.clone(), ry.clone(), random_points],
        setup,
        channel,
    );
    e_w_eval_proof
}
#[allow(non_snake_case)]
pub fn first_layer_sum_check(
    A: &SparseRep,
    B: &SparseRep,
    C: &SparseRep,
    u: &Scalar,
    z: &MultPolynomial,
    E: &MultPolynomial,
    channel: &mut Channel,
) -> InitialSumCheckTranscript {
    let mut E = MultPolynomial::new(E.as_coeffs().clone().to_vec());

    let mut Az = sparse_matrix_multiply(A, z);
    let mut Bz = sparse_matrix_multiply(B, z);
    let mut Cz = sparse_matrix_multiply(C, z);

    let sum_check_rounds = z.len().trailing_zeros() as usize;

    let tau = channel.get_random_points(sum_check_rounds);

    let mut fourcoeffs = MultPolynomial::new(compute_coeff(&tau));

    let mut polynomials = Vec::new();
    let mut random_points = vec![Scalar::ZERO; sum_check_rounds];

    for round in 0..sum_check_rounds {
        let mut eval = [Scalar::ZERO; 4];
        let halfsize = 1 << (sum_check_rounds - round - 1);
        for k in 0..halfsize {
            eval[0] += fourcoeffs.get_coeff(k)
                * (Az.get_coeff(k) * Bz.get_coeff(k) - (*u * Cz.get_coeff(k) + E.get_coeff(k)));

            eval[1] += fourcoeffs.get_coeff(k + halfsize)
                * (Az.get_coeff(k + halfsize) * Bz.get_coeff(k + halfsize)
                    - (*u * Cz.get_coeff(k + halfsize) + E.get_coeff(k + halfsize)));

            eval[2] += (fourcoeffs.get_coeff(k).double() - fourcoeffs.get_coeff(k + halfsize))
                * ((Az.get_coeff(k).double() - Az.get_coeff(k + halfsize))
                    * (Bz.get_coeff(k).double() - Bz.get_coeff(k + halfsize))
                    - (*u * (Cz.get_coeff(k).double() - Cz.get_coeff(k + halfsize))
                        + (E.get_coeff(k).double() - E.get_coeff(k + halfsize))));

            eval[3] += (fourcoeffs.get_coeff(k + halfsize).double() - fourcoeffs.get_coeff(k))
                * ((Az.get_coeff(k + halfsize).double() - Az.get_coeff(k))
                    * (Bz.get_coeff(k + halfsize).double() - Bz.get_coeff(k))
                    - (*u * (Cz.get_coeff(k + halfsize).double() - Cz.get_coeff(k))
                        + (E.get_coeff(k + halfsize).double() - E.get_coeff(k))));
        }

        len_4_interpolate(&mut eval);

        channel.reseed_with_scalars(&eval);

        polynomials.push(Polynomial::new(eval.to_vec()));

        let r_i = channel.get_random_point();

        random_points[round] = r_i;

        fourcoeffs = fourcoeffs.par_fold_by_msb(r_i);
        Az = Az.par_fold_by_msb(r_i);
        Bz = Bz.par_fold_by_msb(r_i);
        Cz = Cz.par_fold_by_msb(r_i);
        E = E.par_fold_by_msb(r_i);
    }
    let Az_claimed_val = Az.get_coeff(0);
    let Bz_claimed_val = Bz.get_coeff(0);
    let Cz_claimed_val = Cz.get_coeff(0);
    let E_claimed_val = E.get_coeff(0);

    InitialSumCheckTranscript {
        polynomials,
        random_points,
        Az_claimed_val,
        Bz_claimed_val,
        Cz_claimed_val,
        E_claimed_val,
    }
}
#[allow(non_snake_case)]
pub fn parallel_sum_checks(
    A: &SparseRep,
    B: &SparseRep,
    C: &SparseRep,
    z: &MultPolynomial,
    rx_basis_evals: Vec<Scalar>,
    channel: &mut Channel,
) -> ParSumCheckTranscript {
    let sum_check_rounds = z.len().trailing_zeros() as usize;
    let A_rx: Vec<Scalar> = A.bind_row_variable(&rx_basis_evals, z.len());

    let B_rx: Vec<Scalar> = B.bind_row_variable(&rx_basis_evals, z.len());

    let C_rx: Vec<Scalar> = C.bind_row_variable(&rx_basis_evals, z.len());

    let A_rx = MultPolynomial::new(A_rx);
    let B_rx = MultPolynomial::new(B_rx);
    let C_rx = MultPolynomial::new(C_rx);

    let mut Z = z.clone();

    let mut batch = vec![A_rx, B_rx, C_rx];
    let random_coeffs = channel.get_random_points(3);

    let mut par_sum_check_polys = Vec::new();
    let mut par_sum_check_random_points = Vec::new();
    for round in 0..sum_check_rounds {
        let mut eval = vec![vec![Scalar::ZERO; 3]; 3];
        let halfsize = 1 << (sum_check_rounds - round - 1);

        let mut comb_poly = vec![Scalar::ZERO; 3];
        for k in 0..halfsize {
            let temp = Z.get_coeff(k + halfsize).double() - Z.get_coeff(k);
            for p in 0..3 {
                eval[p][0] += batch[p].get_coeff(k) * Z.get_coeff(k);
                eval[p][1] += batch[p].get_coeff(k + halfsize) * Z.get_coeff(k + halfsize);

                eval[p][2] +=
                    (batch[p].get_coeff(k + halfsize).double() - batch[p].get_coeff(k)) * temp;
            }
        }

        for p in 0..3 {
            interpolate(&mut eval[p])
        }
        for p in 0..3 {
            for c in 0..3 {
                comb_poly[c] += random_coeffs[p] * eval[p][c];
            }
        }

        par_sum_check_polys.push(Polynomial::new(comb_poly.clone()));

        channel.reseed_with_scalars(&comb_poly);

        let r_i = channel.get_random_point();

        par_sum_check_random_points.push(r_i);

        Z = Z.par_fold_by_msb(r_i);
        for p in 0..3 {
            batch[p] = batch[p].par_fold_by_msb(r_i);
        }
    }
    ParSumCheckTranscript::new(
        par_sum_check_polys,
        par_sum_check_random_points,
        batch[0].get_coeff(0),
        batch[1].get_coeff(0),
        batch[2].get_coeff(0),
        Z.get_coeff(0),
    )
}
pub fn sparse_matrix_multiply(mat: &SparseRep, z: &MultPolynomial) -> MultPolynomial {
    let number_of_rows = mat.fourcoeffs.len();
    let mut a_z_vec: Vec<Scalar> = vec![Scalar::ZERO; z.len()];
    for i in 0..number_of_rows {
        match mat.fourcoeffs.get(&i) {
            Some(a_column_data_vec) => {
                let mut a_z = Scalar::ZERO;
                for column_data in a_column_data_vec.iter() {
                    let variable = z.get_coeff(column_data.column);
                    a_z += variable * column_data.value;
                }
                a_z_vec[i] = a_z;
            }
            None => {
                println!("Not found");
                panic!();
            }
        }
    }
    MultPolynomial(a_z_vec)
}
