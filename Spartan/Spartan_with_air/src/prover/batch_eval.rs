#![allow(unused)]
use bls381::scalar::Scalar;
use bls_curve::bls::AffinePoint;
use channel::Channel;
use crypto_bigint::{Limb, Uint};
use grand_product_with_air::{
    grand_product_common::{FinalLayers, GkrTranscript, LeafLayers},
    prover::{gkr_prover, trace::Trace},
};
use kzg_fft::{commit::kzg2commit, common::KZGFFTDegreeBoundSetup, prover::kzg2_batch_prover};
use kzg_fourier_multilinear::{common::KZGFourierDegreeBoundSetup, prover::batch_eval};
use polynomial::{MultPolynomial, Polynomial};
use preprocessing::SparseMetaData;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
use spartan_common::{BatchSparseEvalProof, BatchSpartanSumCheckTranscript, GKREvaluations};
use std::time::Instant;
use traits::traits::{Field, PrimeField};

use crate::{preprocessing, spartan_common};

pub fn batch_eval_sum_check(
    sparse_metadata: &Vec<SparseMetaData>,
    mut e_rx: Vec<MultPolynomial>,
    mut e_ry: Vec<MultPolynomial>,
    channel: &mut Channel,
) -> BatchSpartanSumCheckTranscript {
    let mut val: Vec<MultPolynomial> = sparse_metadata
        .iter()
        .map(|metadata| metadata.val.clone())
        .collect();

    let sum_check_rounds = val[0].len().trailing_zeros() as usize;

    let random_coeffs = channel.get_random_points(3);
    let mut sum_check_polynomials: Vec<Polynomial> = Vec::new();
    let mut sum_check_random_points: Vec<Scalar> = vec![Scalar::ZERO; sum_check_rounds];

    for i in 0..sum_check_rounds {
        let halfsize = 1usize << (sum_check_rounds - 1usize - i);

        let mut combined_poly = [Scalar::ZERO; 4];
        let mut eval = vec![[Scalar::ZERO; 4]; 3];

        for c in 0..sparse_metadata.len() {
            eval[c] = (0..halfsize).fold([Scalar::ZERO; 4], |mut acc, k| {
                acc[0] +=
                    val[c].get_coeff(2 * k) * e_rx[c].get_coeff(2 * k) * e_ry[c].get_coeff(2 * k);

                acc[1] += val[c].get_coeff(2 * k + 1)
                    * e_rx[c].get_coeff(2 * k + 1)
                    * e_ry[c].get_coeff(2 * k + 1);

                acc[2] += (val[c].get_coeff(2 * k).double() - val[c].get_coeff(2 * k + 1))
                    * (e_rx[c].get_coeff(2 * k).double() - e_rx[c].get_coeff(2 * k + 1))
                    * (e_ry[c].get_coeff(2 * k).double() - e_ry[c].get_coeff(2 * k + 1));

                acc[3] += (val[c].get_coeff(2 * k + 1).double() - val[c].get_coeff(2 * k))
                    * (e_rx[c].get_coeff(2 * k + 1).double() - e_rx[c].get_coeff(2 * k))
                    * (e_ry[c].get_coeff(2 * k + 1).double() - e_ry[c].get_coeff(2 * k));

                acc
            });
            len_4_interpolate(&mut eval[c])
        }

        for c in 0..sparse_metadata.len() {
            for k in 0..4 {
                combined_poly[k] += random_coeffs[c] * eval[c][k]
            }
        }

        channel.reseed_with_scalars(&combined_poly);

        let r_i = channel.get_random_point();

        sum_check_random_points[sum_check_rounds - 1 - i] = r_i;

        sum_check_polynomials.push(Polynomial::new(combined_poly.to_vec()));

        for k in 0..sparse_metadata.len() {
            e_rx[k] = e_rx[k].fold_by_lsb(r_i);
            e_ry[k] = e_ry[k].fold_by_lsb(r_i);
            val[k] = val[k].fold_by_lsb(r_i);
        }
    }

    BatchSpartanSumCheckTranscript::new(
        sum_check_polynomials,
        sum_check_random_points,
        e_rx,
        e_ry,
        val,
    )
}

pub fn batch_offline_mem_check(
    rows: &Vec<&Vec<Scalar>>,
    cols: &Vec<&Vec<Scalar>>,
    read_ts_for_rows: &Vec<&Vec<Scalar>>,
    read_ts_for_cols: &Vec<&Vec<Scalar>>,
    final_ts_for_rows: &Vec<&Vec<Scalar>>,
    final_ts_for_cols: &Vec<&Vec<Scalar>>,
    rx_basis_evals: &Vec<Scalar>,
    ry_basis_evals: &Vec<Scalar>,
    e_rx: &Vec<&Vec<Scalar>>,
    e_ry: &Vec<&Vec<Scalar>>,
    channel: &mut Channel,
    setup: &KZGFFTDegreeBoundSetup,
) -> (
    GkrTranscript,
    GkrTranscript,
    Scalar,
    Scalar,
    Trace,
    Vec<Scalar>,
    Trace,
    Vec<Scalar>,
) {
    let (leaf_layers, final_layers) = get_batch_grand_product_circuits(
        rx_basis_evals,
        ry_basis_evals,
        rows,
        cols,
        e_rx,
        e_ry,
        read_ts_for_rows,
        read_ts_for_cols,
        final_ts_for_rows,
        final_ts_for_cols,
        channel,
    );

    let leaf_layer_1 = leaf_layers.leaf_layer_1();
    let final_layer_1 = final_layers.final_layer_1();

    let (gkr_transcript1, evaluation_point1, trace1, composition_evaluations1) = gkr_prover(
        leaf_layer_1,
        final_layer_1.to_vec(),
        channel,
        &setup.get_setup(leaf_layer_1[0].len()).prover_key,
    );
    let leaf_layer_2 = leaf_layers.leaf_layer_2();
    let final_layer_2 = final_layers.final_layer_2();

    let (gkr_transcript2, evaluation_point2, trace2, composition_evaluations2) = gkr_prover(
        leaf_layer_2,
        final_layer_2.to_vec(),
        channel,
        &setup.get_setup(leaf_layer_2[0].len()).prover_key,
    );

    (
        gkr_transcript1,
        gkr_transcript2,
        evaluation_point1,
        evaluation_point2,
        trace1,
        composition_evaluations1,
        trace2,
        composition_evaluations2,
    )
}

pub fn batch_eval_proof(
    sparse_metadata: Vec<SparseMetaData>,
    channel: &mut Channel,
    setup: &KZGFourierDegreeBoundSetup,
    rx_basis_evals: Vec<Scalar>,
    ry_basis_evals: Vec<Scalar>,
) -> BatchSparseEvalProof {
    let e_rx: Vec<Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| {
            metadata
                .row
                .as_coeffs()
                .par_iter()
                .map(|row_idx| rx_basis_evals[row_idx.0.as_words()[0] as usize])
                .collect()
        })
        .collect();
    let e_ry: Vec<Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| {
            metadata
                .col
                .as_coeffs()
                .par_iter()
                .map(|col_idx| ry_basis_evals[col_idx.0.as_words()[0] as usize])
                .collect()
        })
        .collect();

    let e_rx_polys: Vec<MultPolynomial> = e_rx
        .iter()
        .map(|rx| MultPolynomial::new(rx.to_vec()))
        .collect();
    let e_ry_polys: Vec<MultPolynomial> = e_ry
        .iter()
        .map(|ry| MultPolynomial::new(ry.to_vec()))
        .collect();

    let val: Vec<Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.val.as_coeffs().clone())
        .collect();

    let sum_check_transcript =
        batch_eval_sum_check(&sparse_metadata, e_rx_polys, e_ry_polys, channel);

    let e_rx_evals_sum_check: Vec<Scalar> = sum_check_transcript
        .e_rx
        .iter()
        .map(|poly| poly.get_coeff(0))
        .collect();
    let e_ry_evals_sum_check: Vec<Scalar> = sum_check_transcript
        .e_ry
        .iter()
        .map(|poly| poly.get_coeff(0))
        .collect();
    let val_evals_sum_check: Vec<Scalar> = sum_check_transcript
        .val
        .iter()
        .map(|poly| poly.get_coeff(0))
        .collect();

    let e_rx_polys: Vec<&Vec<Scalar>> = e_rx.iter().collect();
    let e_ry_polys: Vec<&Vec<Scalar>> = e_ry.iter().collect();

    let e_rx_commits = (0..e_rx.len())
        .map(|idx| {
            kzg2commit(
                &e_rx[idx],
                &setup.get_setup(e_rx[idx].len()).setup.prover_key,
            )
        })
        .collect::<Vec<AffinePoint>>();
    let e_ry_commits = (0..e_ry.len())
        .map(|idx| {
            kzg2commit(
                &e_ry[idx],
                &setup.get_setup(e_ry[idx].len()).setup.prover_key,
            )
        })
        .collect::<Vec<AffinePoint>>();

    let rows: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.row.as_coeffs())
        .collect();
    let cols: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.col.as_coeffs())
        .collect();
    let read_ts_for_rows: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.timestamps.read_ts_row.as_coeffs())
        .collect();
    let read_ts_for_cols: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.timestamps.read_ts_col.as_coeffs())
        .collect();
    let final_ts_for_rows: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.timestamps.final_ts_row.as_coeffs())
        .collect();
    let final_ts_for_cols: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.timestamps.final_ts_col.as_coeffs())
        .collect();
    let (
        gkr_transcript1,
        gkr_transcript2,
        z1,
        z2,
        trace1,
        composition_evaluations1,
        trace2,
        composition_evaluations2,
    ) = batch_offline_mem_check(
        &rows,
        &cols,
        &read_ts_for_rows,
        &read_ts_for_cols,
        &final_ts_for_rows,
        &final_ts_for_cols,
        &rx_basis_evals,
        &ry_basis_evals,
        &e_rx_polys,
        &e_ry_polys,
        channel,
        &setup.setup,
    );
    let gkr_evaluations = evaluate(
        z1,
        z2,
        &final_ts_for_rows,
        &final_ts_for_cols,
        &rows,
        &cols,
        &e_rx_polys,
        &e_ry_polys,
        &read_ts_for_rows,
        &read_ts_for_cols,
        &rx_basis_evals,
        &ry_basis_evals,
    );
    let idx_vec = (0..trace1.num_of_rows())
        .map(|idx| Scalar::from(idx as u64))
        .collect::<Vec<Scalar>>();
    let g_trace1 = Scalar::get_root_of_unity(trace1.num_of_rows().trailing_zeros());
    let g_trace2 = Scalar::get_root_of_unity(trace2.num_of_rows().trailing_zeros());
    let mut poly1_next = Vec::new();
    (0..final_ts_for_rows.len()).for_each(|idx| poly1_next.push(final_ts_for_rows[idx].to_vec()));
    (0..final_ts_for_cols.len()).for_each(|idx| poly1_next.push(final_ts_for_cols[idx].to_vec()));
    (0..trace1.num_of_columns()).for_each(|idx| {
        if idx % 2 != 0 {
            poly1_next.push(trace1.column_at(idx))
        }
    });
    poly1_next.push(idx_vec);
    poly1_next.push(rx_basis_evals);
    poly1_next.push(ry_basis_evals);
    let mut poly1_current = poly1_next.clone();
    poly1_current.push(composition_evaluations1);

    let mut eval1_at_z1 = Vec::new();
    (0..final_ts_for_rows.len())
        .for_each(|idx| eval1_at_z1.push(gkr_evaluations.final_ts_for_rows_evaluations_at_z1[idx]));
    (0..final_ts_for_cols.len())
        .for_each(|idx| eval1_at_z1.push(gkr_evaluations.final_ts_for_cols_evaluations_at_z1[idx]));
    let trace_frame = gkr_transcript1.get_odd_frame().get_current_frame();
    let composition_frame = gkr_transcript1.get_odd_frame().get_composition_frame();
    (0..trace_frame.len()).for_each(|idx| eval1_at_z1.push(trace_frame[idx]));
    eval1_at_z1.push(gkr_evaluations.idx_at_z1);
    eval1_at_z1.push(gkr_evaluations.rx_basis_evals_evaluation_at_z1);
    eval1_at_z1.push(gkr_evaluations.ry_basis_evals_evaluation_at_z1);
    eval1_at_z1.push(composition_frame);

    let mut eval1_at_g_z1 = Vec::new();
    (0..final_ts_for_rows.len()).for_each(|idx| {
        eval1_at_g_z1.push(gkr_evaluations.final_ts_for_rows_evaluations_at_g_z1[idx])
    });
    (0..final_ts_for_cols.len()).for_each(|idx| {
        eval1_at_g_z1.push(gkr_evaluations.final_ts_for_cols_evaluations_at_g_z1[idx])
    });
    let next_frame = gkr_transcript1.get_odd_frame().get_next_frame();
    (0..next_frame.len()).for_each(|idx| eval1_at_g_z1.push(next_frame[idx]));
    eval1_at_g_z1.push(gkr_evaluations.idx_at_g_z1);
    eval1_at_g_z1.push(gkr_evaluations.rx_basis_evals_evaluation_at_g_z1);
    eval1_at_g_z1.push(gkr_evaluations.ry_basis_evals_evaluation_at_g_z1);
    ///////
    let mut poly2_next = Vec::new();
    (0..rows.len()).for_each(|idx| poly2_next.push(rows[idx].to_vec()));
    (0..cols.len()).for_each(|idx| poly2_next.push(cols[idx].to_vec()));
    (0..cols.len()).for_each(|idx| poly2_next.push(e_rx[idx].to_vec()));
    (0..cols.len()).for_each(|idx| poly2_next.push(e_ry[idx].to_vec()));
    (0..cols.len()).for_each(|idx| poly2_next.push(read_ts_for_rows[idx].to_vec()));
    (0..cols.len()).for_each(|idx| poly2_next.push(read_ts_for_cols[idx].to_vec()));
    (0..trace2.num_of_columns()).for_each(|idx| {
        if idx % 2 != 0 {
            poly2_next.push(trace2.column_at(idx))
        }
    });
    let mut poly2_current = poly2_next.clone();
    poly2_current.push(composition_evaluations2);

    let mut eval1_at_z2 = Vec::new();
    (0..final_ts_for_rows.len())
        .for_each(|idx| eval1_at_z2.push(gkr_evaluations.rows_evaluations_at_z2[idx]));
    (0..final_ts_for_cols.len())
        .for_each(|idx| eval1_at_z2.push(gkr_evaluations.cols_evaluations_at_z2[idx]));
    (0..final_ts_for_cols.len())
        .for_each(|idx| eval1_at_z2.push(gkr_evaluations.e_rx_evaluations_at_z2[idx]));
    (0..final_ts_for_cols.len())
        .for_each(|idx| eval1_at_z2.push(gkr_evaluations.e_ry_evaluations_at_z2[idx]));
    (0..final_ts_for_cols.len())
        .for_each(|idx| eval1_at_z2.push(gkr_evaluations.read_ts_for_rows_evaluations_at_z2[idx]));
    (0..final_ts_for_cols.len())
        .for_each(|idx| eval1_at_z2.push(gkr_evaluations.read_ts_for_cols_evaluations_at_z2[idx]));
    let trace_frame = gkr_transcript2.get_odd_frame().get_current_frame();
    let composition_frame = gkr_transcript2.get_odd_frame().get_composition_frame();
    (0..trace_frame.len()).for_each(|idx| eval1_at_z2.push(trace_frame[idx]));
    eval1_at_z2.push(composition_frame);

    let mut eval1_at_g_z2 = Vec::new();
    (0..final_ts_for_rows.len())
        .for_each(|idx| eval1_at_g_z2.push(gkr_evaluations.rows_evaluations_at_g_z2[idx]));
    (0..final_ts_for_cols.len())
        .for_each(|idx| eval1_at_g_z2.push(gkr_evaluations.cols_evaluations_at_g_z2[idx]));
    (0..final_ts_for_rows.len())
        .for_each(|idx| eval1_at_g_z2.push(gkr_evaluations.e_rx_evaluations_at_g_z2[idx]));
    (0..final_ts_for_cols.len())
        .for_each(|idx| eval1_at_g_z2.push(gkr_evaluations.e_ry_evaluations_at_g_z2[idx]));
    (0..final_ts_for_rows.len()).for_each(|idx| {
        eval1_at_g_z2.push(gkr_evaluations.read_ts_for_rows_evaluations_at_g_z2[idx])
    });
    (0..final_ts_for_cols.len()).for_each(|idx| {
        eval1_at_g_z2.push(gkr_evaluations.read_ts_for_cols_evaluations_at_g_z2[idx])
    });

    let next_frame = gkr_transcript2.get_odd_frame().get_next_frame();
    (0..next_frame.len()).for_each(|idx| eval1_at_g_z2.push(next_frame[idx]));

    let combiners = channel.get_random_points(poly1_current.len());
    assert_eq!(eval1_at_z1.len(), poly1_current.len());

    let proof1 = kzg2_batch_prover(
        &poly1_current,
        z1,
        &eval1_at_z1,
        &combiners,
        &setup.setup.get_setup(poly1_current[0].len()).prover_key,
    );
    //
    let proof2 = kzg2_batch_prover(
        &poly1_next,
        g_trace1 * z1,
        &eval1_at_g_z1,
        &combiners[0..&combiners.len() - 1].to_vec(),
        &setup.setup.get_setup(poly1_next[0].len()).prover_key,
    );

    let combiners = channel.get_random_points(poly2_current.len());
    let proof3 = kzg2_batch_prover(
        &poly2_current,
        z2,
        &eval1_at_z2,
        &combiners,
        &setup.setup.get_setup(poly2_current[0].len()).prover_key,
    );

    let proof4 = kzg2_batch_prover(
        &poly2_next,
        g_trace2 * z2,
        &eval1_at_g_z2,
        &combiners[0..&combiners.len() - 1].to_vec(),
        &setup.setup.get_setup(poly2_next[0].len()).prover_key,
    );

    /////
    let start_time = Instant::now();
    let batch_eval_coeffs_sum_check = channel.get_random_points(3 * sparse_metadata.len());
    let sum_check_eval_proof = batch_eval(
        &[e_rx, e_ry, val].concat(),
        &batch_eval_coeffs_sum_check,
        [
            e_rx_evals_sum_check.clone(),
            e_ry_evals_sum_check.clone(),
            val_evals_sum_check.clone(),
        ]
        .concat(),
        sum_check_transcript.random_points.clone(),
        setup,
        channel,
    );
    BatchSparseEvalProof::new(
        sum_check_transcript,
        gkr_transcript1,
        gkr_transcript2,
        e_rx_evals_sum_check,
        e_ry_evals_sum_check,
        val_evals_sum_check,
        sum_check_eval_proof,
        e_rx_commits,
        e_ry_commits,
        gkr_evaluations,
        proof1,
        proof2,
        proof3,
        proof4,
    )
}

pub(crate) fn evaluate(
    z1: Scalar,
    z2: Scalar,
    final_ts_for_rows: &Vec<&Vec<Scalar>>,
    final_ts_for_cols: &Vec<&Vec<Scalar>>,
    rows: &Vec<&Vec<Scalar>>,
    cols: &Vec<&Vec<Scalar>>,
    e_rx: &Vec<&Vec<Scalar>>,
    e_ry: &Vec<&Vec<Scalar>>,
    read_ts_for_rows: &Vec<&Vec<Scalar>>,
    read_ts_for_cols: &Vec<&Vec<Scalar>>,
    rx_basis_evals: &Vec<Scalar>,
    ry_basis_evals: &Vec<Scalar>,
) -> GKREvaluations {
    let trace_length1 = final_ts_for_rows[0].len();
    let trace_length2 = rows[0].len();

    assert!(
        trace_length1.is_power_of_two(),
        " program code length should be power of 2"
    );
    assert!(
        trace_length2.is_power_of_two(),
        " length of rows for each circuit must be power of 2"
    );
    let g_trace1 = Scalar::get_root_of_unity(trace_length1.trailing_zeros());
    let g_trace2 = Scalar::get_root_of_unity(trace_length2.trailing_zeros());
    let g_z1 = g_trace1 * z1;
    let g_z2 = g_trace2 * z2;

    let log2_degree = trace_length1.trailing_zeros();
    let omega = g_trace1;
    let mut omega_powers = vec![Scalar::ONE; trace_length1];
    for idx in 1..trace_length1 {
        let temp = omega_powers[idx - 1] * omega;
        omega_powers[idx] = temp;
    }
    //Compute alpha_i
    let degree_inverse = (Scalar::from(trace_length1 as u64)).invert().unwrap();
    let beta_z1 = (0..trace_length1)
        .into_par_iter()
        .map(|i| {
            let mut temp = Scalar::ONE;
            for j in 0..log2_degree {
                temp *= Scalar::ONE
                    + (omega_powers[(trace_length1 - i) % trace_length1] * z1).power_by([
                        1_u64 << j,
                        0,
                        0,
                        0,
                    ])
            }
            temp * degree_inverse
        })
        .collect::<Vec<Scalar>>();

    let beta_g_z1 = (0..trace_length1)
        .into_par_iter()
        .map(|i| {
            let mut temp = Scalar::ONE;
            for j in 0..log2_degree {
                temp *= Scalar::ONE
                    + (omega_powers[(trace_length1 - i) % trace_length1] * g_z1).power_by([
                        1_u64 << j,
                        0,
                        0,
                        0,
                    ])
            }
            temp * degree_inverse
        })
        .collect::<Vec<Scalar>>();

    let log2_degree = trace_length2.trailing_zeros();
    let degree_inverse = (Scalar::from(trace_length2 as u64)).invert().unwrap();

    let omega = g_trace2;
    let mut omega_powers = vec![Scalar::ONE; trace_length2];
    for idx in 1..trace_length2 {
        let temp = omega_powers[idx - 1] * omega;
        omega_powers[idx] = temp;
    }
    let beta_z2 = (0..trace_length2)
        .into_par_iter()
        .map(|i| {
            let mut temp = Scalar::ONE;
            for j in 0..log2_degree {
                temp *= Scalar::ONE
                    + (omega_powers[(trace_length2 - i) % trace_length2] * z2).power_by([
                        1_u64 << j,
                        0,
                        0,
                        0,
                    ])
            }
            temp * degree_inverse
        })
        .collect::<Vec<Scalar>>();

    let beta_g_z2 = (0..trace_length2)
        .into_par_iter()
        .map(|i| {
            let mut temp = Scalar::ONE;
            for j in 0..log2_degree {
                temp *= Scalar::ONE
                    + (omega_powers[(trace_length2 - i) % trace_length2] * g_z2).power_by([
                        1_u64 << j,
                        0,
                        0,
                        0,
                    ])
            }
            temp * degree_inverse
        })
        .collect::<Vec<Scalar>>();

    let final_ts_for_rows_evaluations_at_z1 = (0..final_ts_for_rows.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + final_ts_for_rows[outer_idx][idx] * beta_z1[idx]
            })
        })
        .collect::<Vec<Scalar>>();
    let final_ts_for_rows_evaluations_at_g_z1 = (0..final_ts_for_rows.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + final_ts_for_rows[outer_idx][idx] * beta_g_z1[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let final_ts_for_cols_evaluations_at_z1 = (0..final_ts_for_cols.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + final_ts_for_cols[outer_idx][idx] * beta_z1[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let final_ts_for_cols_evaluations_at_g_z1 = (0..final_ts_for_cols.len())
        .map(|outer_idx| {
            (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
                acc + final_ts_for_cols[outer_idx][idx] * beta_g_z1[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let rows_evaluations_at_z2 = (0..rows.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + rows[outer_idx][idx] * beta_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let rows_evaluations_at_g_z2 = (0..rows.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + rows[outer_idx][idx] * beta_g_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let cols_evaluations_at_z2 = (0..cols.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + cols[outer_idx][idx] * beta_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let cols_evaluations_at_g_z2 = (0..cols.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + cols[outer_idx][idx] * beta_g_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let e_rx_evaluations_at_z2 = (0..e_rx.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + e_rx[outer_idx][idx] * beta_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let e_rx_evaluations_at_g_z2 = (0..e_rx.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + e_rx[outer_idx][idx] * beta_g_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let e_ry_evaluations_at_z2 = (0..e_ry.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + e_ry[outer_idx][idx] * beta_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let e_ry_evaluations_at_g_z2 = (0..e_ry.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + e_ry[outer_idx][idx] * beta_g_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let read_ts_for_rows_evaluations_at_z2 = (0..read_ts_for_rows.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + read_ts_for_rows[outer_idx][idx] * beta_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let read_ts_for_rows_evaluations_at_g_z2 = (0..read_ts_for_rows.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + read_ts_for_rows[outer_idx][idx] * beta_g_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let read_ts_for_cols_evaluations_at_z2 = (0..read_ts_for_cols.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + read_ts_for_cols[outer_idx][idx] * beta_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();

    let read_ts_for_cols_evaluations_at_g_z2 = (0..read_ts_for_cols.len())
        .map(|outer_idx| {
            (0..trace_length2).fold(Scalar::ZERO, |acc, idx| {
                acc + read_ts_for_cols[outer_idx][idx] * beta_g_z2[idx]
            })
        })
        .collect::<Vec<Scalar>>();
    let rx_basis_evals_evaluation_at_z1 = (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
        acc + rx_basis_evals[idx] * beta_z1[idx]
    });
    let rx_basis_evals_evaluation_at_g_z1 = (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
        acc + rx_basis_evals[idx] * beta_g_z1[idx]
    });
    let ry_basis_evals_evaluation_at_z1 = (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
        acc + ry_basis_evals[idx] * beta_z1[idx]
    });
    let ry_basis_evals_evaluation_at_g_z1 = (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
        acc + ry_basis_evals[idx] * beta_g_z1[idx]
    });
    let idx_at_z1 = (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
        acc + Scalar::from(idx as u64) * beta_z1[idx]
    });
    let idx_at_g_z1 = (0..trace_length1).fold(Scalar::ZERO, |acc, idx| {
        acc + Scalar::from(idx as u64) * beta_g_z1[idx]
    });
    GKREvaluations::new(
        rx_basis_evals_evaluation_at_z1,
        ry_basis_evals_evaluation_at_z1,
        rx_basis_evals_evaluation_at_g_z1,
        ry_basis_evals_evaluation_at_g_z1,
        final_ts_for_rows_evaluations_at_z1,
        final_ts_for_cols_evaluations_at_z1,
        final_ts_for_rows_evaluations_at_g_z1,
        final_ts_for_cols_evaluations_at_g_z1,
        rows_evaluations_at_z2,
        rows_evaluations_at_g_z2,
        cols_evaluations_at_z2,
        cols_evaluations_at_g_z2,
        e_rx_evaluations_at_z2,
        e_rx_evaluations_at_g_z2,
        e_ry_evaluations_at_z2,
        e_ry_evaluations_at_g_z2,
        read_ts_for_rows_evaluations_at_z2,
        read_ts_for_rows_evaluations_at_g_z2,
        read_ts_for_cols_evaluations_at_z2,
        read_ts_for_cols_evaluations_at_g_z2,
        idx_at_z1,
        idx_at_g_z1,
    )
}
pub fn get_batch_grand_product_circuits(
    rx_basis_evals: &Vec<Scalar>,
    ry_basis_evals: &Vec<Scalar>,
    rows: &Vec<&Vec<Scalar>>,
    cols: &Vec<&Vec<Scalar>>,
    e_rx: &Vec<&Vec<Scalar>>,
    e_ry: &Vec<&Vec<Scalar>>,
    read_ts_for_rows: &Vec<&Vec<Scalar>>,
    read_ts_for_cols: &Vec<&Vec<Scalar>>,
    final_ts_for_rows: &Vec<&Vec<Scalar>>,
    final_ts_for_cols: &Vec<&Vec<Scalar>>,
    channel: &mut Channel,
) -> (LeafLayers, FinalLayers) {
    let reads_length = rows[0].len();
    let init_length = rx_basis_evals.len();
    let n_circuits = rows.len() + cols.len();

    //The depth for the circuit of ws and s will be log_2(program_length)
    let depth_1 = init_length.trailing_zeros() as usize;

    //The depth for the circuit of ws_updates and rs will be log_2(instructions_length)
    let depth_2 = reads_length.trailing_zeros() as usize;

    let gamma_tau = channel.get_random_points(2);

    let mut w_init_leaf_layers = vec![vec![Scalar::ONE; 1 << depth_1]; n_circuits];
    let mut s_leaf_layers = vec![vec![Scalar::ONE; 1 << depth_1]; n_circuits];
    let mut w_update_leaf_layers = vec![vec![Scalar::ONE; 1 << depth_2]; n_circuits];
    let mut r_leaf_layers = vec![vec![Scalar::ONE; 1 << depth_2]; n_circuits];

    let mut w_init_final_layers = Vec::new();
    let mut s_final_layers = Vec::new();
    let mut w_update_final_layers = Vec::new();
    let mut r_final_layers = Vec::new();

    //The fingerprint for an entry is a polynomial expression of the form, gamma^2 * index + gamma * value + timestamp - tau.
    for c in 0..n_circuits / 2 {
        for i in 0..init_length {
            w_init_leaf_layers[c][i] = gamma_tau[0].square() * Scalar::from(i as u32)
                + gamma_tau[0] * rx_basis_evals[i]
                - gamma_tau[1];

            s_leaf_layers[c][i] = w_init_leaf_layers[c][i] + final_ts_for_rows[c][i]
        }

        for i in 0..reads_length {
            r_leaf_layers[c][i] = gamma_tau[0].square() * rows[c][i]
                + gamma_tau[0] * e_rx[c][i]
                + read_ts_for_rows[c][i]
                - gamma_tau[1];
            w_update_leaf_layers[c][i] = r_leaf_layers[c][i] + Scalar::ONE;
        }
        w_init_final_layers.push(
            w_init_leaf_layers[c]
                .par_iter()
                .fold(|| Scalar::ONE, |acc, value| acc * *value)
                .reduce(|| Scalar::ONE, |acc, x| acc * x),
        );
        s_final_layers.push(
            s_leaf_layers[c]
                .par_iter()
                .fold(|| Scalar::ONE, |acc, value| acc * *value)
                .reduce(|| Scalar::ONE, |acc, x| acc * x),
        );
        w_update_final_layers.push(
            w_update_leaf_layers[c]
                .par_iter()
                .fold(|| Scalar::ONE, |acc, value| acc * *value)
                .reduce(|| Scalar::ONE, |acc, x| acc * x),
        );
        r_final_layers.push(
            r_leaf_layers[c]
                .par_iter()
                .fold(|| Scalar::ONE, |acc, value| acc * *value)
                .reduce(|| Scalar::ONE, |acc, x| acc * x),
        );
    }
    for c in n_circuits / 2..n_circuits {
        for i in 0..init_length {
            w_init_leaf_layers[c][i] = gamma_tau[0].square() * Scalar::from(i as u32)
                + gamma_tau[0] * ry_basis_evals[i]
                - gamma_tau[1];

            s_leaf_layers[c][i] = w_init_leaf_layers[c][i] + final_ts_for_cols[c - 3][i]
        }

        for i in 0..reads_length {
            r_leaf_layers[c][i] = gamma_tau[0].square() * cols[c - 3][i]
                + gamma_tau[0] * e_ry[c - 3][i]
                + read_ts_for_cols[c - 3][i]
                - gamma_tau[1];
            w_update_leaf_layers[c][i] = r_leaf_layers[c][i] + Scalar::ONE
        }
        w_init_final_layers.push(
            w_init_leaf_layers[c]
                .par_iter()
                .fold(|| Scalar::ONE, |acc, value| acc * *value)
                .reduce(|| Scalar::ONE, |acc, x| acc * x),
        );
        s_final_layers.push(
            s_leaf_layers[c]
                .par_iter()
                .fold(|| Scalar::ONE, |acc, value| acc * *value)
                .reduce(|| Scalar::ONE, |acc, x| acc * x),
        );
        w_update_final_layers.push(
            w_update_leaf_layers[c]
                .par_iter()
                .fold(|| Scalar::ONE, |acc, value| acc * *value)
                .reduce(|| Scalar::ONE, |acc, x| acc * x),
        );
        r_final_layers.push(
            r_leaf_layers[c]
                .par_iter()
                .fold(|| Scalar::ONE, |acc, value| acc * *value)
                .reduce(|| Scalar::ONE, |acc, x| acc * x),
        );
    }

    let leaf_layer = LeafLayers::new(
        [w_init_leaf_layers, s_leaf_layers].concat(),
        [w_update_leaf_layers, r_leaf_layers].concat(),
    );

    let final_layer = FinalLayers::new(
        [w_init_final_layers, s_final_layers].concat(),
        [w_update_final_layers, r_final_layers].concat(),
    );
    (leaf_layer, final_layer)
}

const SIX_INV: Scalar = Scalar(Uint {
    limbs: [
        Limb(15372286724512153601),
        Limb(1954008828163476649),
        Limb(9224930440102993583),
        Limb(6961264049553707793),
    ],
});
const TWO_INV: Scalar = Scalar(Uint {
    limbs: [
        Limb(9223372034707292161),
        Limb(12240451741123816959),
        Limb(1845609449319885826),
        Limb(4176758429732224676),
    ],
});

pub fn len_4_interpolate(evaluations: &mut [Scalar; 4]) {
    let t0 = TWO_INV * (evaluations[1] + evaluations[2] - evaluations[0].double());
    let t1 = evaluations[1] - evaluations[2] + evaluations[0] + t0.double().double();
    let t2 = SIX_INV * (evaluations[3] - t1);
    *evaluations = [
        evaluations[0],
        evaluations[1] - (evaluations[0] + t0 + t2),
        t0,
        t2,
    ]
}
