#![allow(unused)]
use crate::preprocessing;
use crate::spartan_common::{
    BatchGrandProductCircuits, BatchSparseEvalProof, BatchSpartanSumCheckTranscript,
};
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use channel::Channel;
use grand_product_with_gkr::gkr_common::{CircuitBinaryTree, CircuitLayer, GkrTranscript};
use grand_product_with_gkr::prover::{gkr_prover, len_4_interpolate};
use multilinear_kzg::common::SRS;
use multilinear_kzg::prover::{batch_eval, commit};
use polynomial::{MultPolynomial, Polynomial};
use preprocessing::{compute_coeff, SparseMetaData};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use traits::traits::Field;
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
            eval[c] = (0..halfsize)
                .into_iter()
                .fold([Scalar::ZERO; 4], |mut acc, k| {
                    acc[0] += val[c].get_coeff(2 * k)
                        * e_rx[c].get_coeff(2 * k)
                        * e_ry[c].get_coeff(2 * k);

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
) -> (GkrTranscript, GkrTranscript, usize) {
    let mut circuits = get_batch_grand_product_circuits(
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
    let circuit2_depth = rows[0].len().trailing_zeros() as usize;
    let mut circuit = Vec::new();
    circuit.append(&mut circuits.w_init_circuit);
    circuit.append(&mut circuits.s_circuit);

    let transcript1 = gkr_prover(&mut circuit, channel);
    let mut circuit = Vec::new();
    circuit.append(&mut circuits.w_update_circuit);
    circuit.append(&mut circuits.r_circuit);
    let transcript2 = gkr_prover(&mut circuit, channel);
    (transcript1, transcript2, circuit2_depth)
}

pub fn batch_eval_proof(
    sparse_metadata: Vec<SparseMetaData>,
    eval_point: &Vec<Scalar>,
    srs: &SRS<BlsCurve>,
    channel: &mut Channel,
) -> BatchSparseEvalProof {
    let (rx, ry) = eval_point.split_at(eval_point.len() / 2);
    let rx_basis_evals = compute_coeff(&rx.to_vec());
    let ry_basis_evals = compute_coeff(&ry.to_vec());

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

    let e_rx_refs: Vec<&Vec<Scalar>> = e_rx.iter().collect();
    let e_ry_refs: Vec<&Vec<Scalar>> = e_ry.iter().collect();

    let e_rx_commits = e_rx_polys
        .iter()
        .map(|rx_poly| commit(rx_poly.as_coeffs(), srs))
        .collect();
    let e_ry_commits = e_ry_polys
        .iter()
        .map(|ry_poly| commit(ry_poly.as_coeffs(), srs))
        .collect();

    let val: Vec<&Vec<Scalar>> = sparse_metadata
        .iter()
        .map(|metadata| metadata.val.as_coeffs())
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

    let rx_basis_evals = MultPolynomial::new(rx_basis_evals);
    let ry_basis_evals = MultPolynomial::new(ry_basis_evals);

    let e_rx_polys: Vec<&Vec<Scalar>> = e_rx.iter().collect();
    let e_ry_polys: Vec<&Vec<Scalar>> = e_ry.iter().collect();

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

    let (gkr_transcript1, gkr_transcript2, circuit2_depth) = batch_offline_mem_check(
        &rows,
        &cols,
        &read_ts_for_rows,
        &read_ts_for_cols,
        &final_ts_for_rows,
        &final_ts_for_cols,
        rx_basis_evals.as_coeffs(),
        ry_basis_evals.as_coeffs(),
        &e_rx_polys,
        &e_ry_polys,
        channel,
    );

    let final_point_basis_evals = compute_coeff(&gkr_transcript1.final_layer_point);

    let final_ts_evals_row_mem_check: Vec<Scalar> = final_ts_for_rows
        .into_iter()
        .map(|final_ts| {
            final_ts
                .par_iter()
                .zip(final_point_basis_evals.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();

    let final_ts_evals_col_mem_check: Vec<Scalar> = final_ts_for_cols
        .into_iter()
        .map(|final_ts| {
            final_ts
                .par_iter()
                .zip(final_point_basis_evals.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();

    let final_point_basis_evals = compute_coeff(&gkr_transcript2.final_layer_point);

    let row_evals_mem_check: Vec<_> = rows
        .iter()
        .map(|row| {
            row.par_iter()
                .zip(final_point_basis_evals.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();
    let col_evals_mem_check: Vec<_> = cols
        .iter()
        .map(|col| {
            col.par_iter()
                .zip(final_point_basis_evals.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();
    let read_ts_evals_row_mem_check: Vec<_> = read_ts_for_rows
        .iter()
        .map(|read_ts_row| {
            read_ts_row
                .par_iter()
                .zip(final_point_basis_evals.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();
    let read_ts_evals_col_mem_check: Vec<_> = read_ts_for_cols
        .iter()
        .map(|read_ts_col| {
            read_ts_col
                .par_iter()
                .zip(final_point_basis_evals.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();
    let e_rx_evals_mem_check: Vec<_> = e_rx_polys
        .iter()
        .map(|e_rx| {
            e_rx.par_iter()
                .zip(final_point_basis_evals.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();
    let e_ry_evals_mem_check: Vec<_> = e_ry_polys
        .iter()
        .map(|e_ry| {
            e_ry.par_iter()
                .zip(final_point_basis_evals.par_iter())
                .map(|(coeff, basis)| *coeff * *basis)
                .reduce(|| Scalar::ZERO, |acc, g| acc + g)
        })
        .collect();
    let batch_eval_coeffs_sum_check = channel.get_random_points(3 * sparse_metadata.len());

    let sum_check_eval_proof = batch_eval(
        &[e_rx_refs, e_ry_refs, val].concat(),
        &[
            e_rx_evals_sum_check.clone(),
            e_ry_evals_sum_check.clone(),
            val_evals_sum_check.clone(),
        ]
        .concat(),
        &sum_check_transcript.random_points,
        &batch_eval_coeffs_sum_check,
        srs,
    );

    let batch_eval_coeffs_gkr = channel.get_random_points(8 * sparse_metadata.len());

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

    let e_rx_refs: Vec<&Vec<Scalar>> = e_rx.iter().map(|rx| rx).collect();
    let e_ry_refs: Vec<&Vec<Scalar>> = e_ry.iter().map(|ry| ry).collect();

    let gkr_batch_eval_proof1 = batch_eval(
        &[final_ts_for_rows, final_ts_for_cols].concat(),
        &[
            final_ts_evals_row_mem_check.clone(),
            final_ts_evals_col_mem_check.clone(),
        ]
        .concat(),
        &gkr_transcript1.final_layer_point,
        &batch_eval_coeffs_gkr,
        srs,
    );
    let gkr_batch_eval_proof2 = batch_eval(
        &[
            rows,
            cols,
            read_ts_for_rows,
            read_ts_for_cols,
            e_rx_refs,
            e_ry_refs,
        ]
        .concat(),
        &[
            row_evals_mem_check.clone(),
            col_evals_mem_check.clone(),
            read_ts_evals_row_mem_check.clone(),
            read_ts_evals_col_mem_check.clone(),
            e_rx_evals_mem_check.clone(),
            e_ry_evals_mem_check.clone(),
        ]
        .concat(),
        &gkr_transcript2.final_layer_point,
        &batch_eval_coeffs_gkr,
        srs,
    );

    BatchSparseEvalProof::new(
        sum_check_transcript.clone(),
        gkr_transcript1,
        gkr_transcript2,
        e_rx_evals_sum_check,
        e_ry_evals_sum_check,
        val_evals_sum_check,
        sum_check_eval_proof,
        e_rx_commits,
        e_ry_commits,
        final_ts_evals_row_mem_check,
        final_ts_evals_col_mem_check,
        row_evals_mem_check,
        col_evals_mem_check,
        read_ts_evals_row_mem_check,
        read_ts_evals_col_mem_check,
        e_rx_evals_mem_check,
        e_ry_evals_mem_check,
        gkr_batch_eval_proof1,
        gkr_batch_eval_proof2,
        circuit2_depth,
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
) -> BatchGrandProductCircuits {
    let reads_length = rows[0].len();
    let init_length = rx_basis_evals.len();
    let n_circuits = rows.len() + cols.len();

    //The depth for the circuit of ws and s will be log_2(program_length)
    let depth_1 = init_length.trailing_zeros() as usize;

    //The depth for the circuit of ws_updates and rs will be log_2(instructions_length)
    let depth_2 = reads_length.trailing_zeros() as usize;

    let mut w_init_circuit_layers = vec![Vec::new(); n_circuits];
    let mut s_circuit_layers = vec![Vec::new(); n_circuits];
    let mut r_circuit_layers = vec![Vec::new(); n_circuits];
    let mut w_update_circuit_layers = vec![Vec::new(); n_circuits];

    let mut w_init_circuits = Vec::new();
    let mut s_circuits = Vec::new();
    let mut w_update_circuits = Vec::new();
    let mut r_circuits = Vec::new();

    //Draw of 2 random integers for the Reed-Solomon fingerprint
    let gamma_tau = channel.get_random_points(2);

    //The fingerprint for an entry is a polynomial expression of the form, gamma^2 * index + gamma * value + timestamp - tau.
    for c in 0..n_circuits / 2 {
        let mut w_init_leaves = vec![Scalar::ONE; 1 << depth_1];
        let mut s_leaves = vec![Scalar::ONE; 1 << depth_1];
        let mut w_update_leaves = vec![Scalar::ONE; 1 << depth_2];
        let mut r_leaves = vec![Scalar::ONE; 1 << depth_2];

        for i in 0..init_length {
            w_init_leaves[i] = gamma_tau[0].square() * Scalar::from(i as u32)
                + gamma_tau[0] * rx_basis_evals[i]
                - gamma_tau[1];

            s_leaves[i] = gamma_tau[0].square() * Scalar::from(i as u32)
                + gamma_tau[0] * rx_basis_evals[i]
                + final_ts_for_rows[c][i]
                - gamma_tau[1];
        }

        for i in 0..reads_length {
            w_update_leaves[i] = gamma_tau[0].square() * rows[c][i]
                + gamma_tau[0] * e_rx[c][i]
                + (read_ts_for_rows[c][i] + Scalar::ONE)
                - gamma_tau[1];

            r_leaves[i] = gamma_tau[0].square() * rows[c][i]
                + gamma_tau[0] * e_rx[c][i]
                + read_ts_for_rows[c][i]
                - gamma_tau[1];
        }

        w_init_circuit_layers[c].push(CircuitLayer::new(w_init_leaves));
        s_circuit_layers[c].push(CircuitLayer::new(s_leaves));
        w_update_circuit_layers[c].push(CircuitLayer::new(w_update_leaves));
        r_circuit_layers[c].push(CircuitLayer::new(r_leaves));

        //Here we recurse over the layers and compute the circuit description for the grand product argument,
        //which we will verify with the GKR protocol.
        for k in 1..depth_1 + 1 {
            let layer_size = 1 << (depth_1 - k);
            let mut temp_w_init = vec![Scalar::ZERO; layer_size];
            let mut temp_s = vec![Scalar::ZERO; layer_size];

            for i in 0..layer_size {
                temp_w_init[i] = w_init_circuit_layers[c][k - 1].get_position(2 * i)
                    * w_init_circuit_layers[c][k - 1].get_position(2 * i + 1);
                temp_s[i] = s_circuit_layers[c][k - 1].get_position(2 * i)
                    * s_circuit_layers[c][k - 1].get_position(2 * i + 1);
            }

            w_init_circuit_layers[c].push(CircuitLayer::new(temp_w_init));
            s_circuit_layers[c].push(CircuitLayer::new(temp_s));
        }

        for k in 1..depth_2 + 1 {
            let layer_size = 1 << (depth_2 - k);
            let mut temp_w_update = vec![Scalar::ZERO; layer_size];
            let mut temp_r = vec![Scalar::ZERO; layer_size];
            for i in 0..layer_size {
                temp_w_update[i] = w_update_circuit_layers[c][k - 1].get_position(2 * i)
                    * w_update_circuit_layers[c][k - 1].get_position(2 * i + 1);
                temp_r[i] = r_circuit_layers[c][k - 1].get_position(2 * i)
                    * r_circuit_layers[c][k - 1].get_position(2 * i + 1);
            }

            w_update_circuit_layers[c].push(CircuitLayer::new(temp_w_update));
            r_circuit_layers[c].push(CircuitLayer::new(temp_r));
        }

        w_init_circuits.push(CircuitBinaryTree::new(w_init_circuit_layers[c].clone()));
        s_circuits.push(CircuitBinaryTree::new(s_circuit_layers[c].clone()));
        w_update_circuits.push(CircuitBinaryTree::new(w_update_circuit_layers[c].clone()));
        r_circuits.push(CircuitBinaryTree::new(r_circuit_layers[c].clone()));
    }
    for c in n_circuits / 2..n_circuits {
        let mut w_init_leaves = vec![Scalar::ONE; 1 << depth_1];
        let mut s_leaves = vec![Scalar::ONE; 1 << depth_1];
        let mut w_update_leaves = vec![Scalar::ONE; 1 << depth_2];
        let mut r_leaves = vec![Scalar::ONE; 1 << depth_2];

        for i in 0..init_length {
            w_init_leaves[i] = gamma_tau[0].square() * Scalar::from(i as u32)
                + gamma_tau[0] * ry_basis_evals[i]
                - gamma_tau[1];

            s_leaves[i] = gamma_tau[0].square() * Scalar::from(i as u32)
                + gamma_tau[0] * ry_basis_evals[i]
                + final_ts_for_cols[c - 3][i]
                - gamma_tau[1];
        }

        for i in 0..reads_length {
            w_update_leaves[i] = gamma_tau[0].square() * cols[c - 3][i]
                + gamma_tau[0] * e_ry[c - 3][i]
                + (read_ts_for_cols[c - 3][i] + Scalar::ONE)
                - gamma_tau[1];

            r_leaves[i] = gamma_tau[0].square() * cols[c - 3][i]
                + gamma_tau[0] * e_ry[c - 3][i]
                + read_ts_for_cols[c - 3][i]
                - gamma_tau[1];
        }

        w_init_circuit_layers[c].push(CircuitLayer::new(w_init_leaves));
        s_circuit_layers[c].push(CircuitLayer::new(s_leaves));
        w_update_circuit_layers[c].push(CircuitLayer::new(w_update_leaves));
        r_circuit_layers[c].push(CircuitLayer::new(r_leaves));

        //Here we recurse over the layers and compute the circuit description for the grand product argument,
        //which we will verify with the GKR protocol.
        for k in 1..depth_1 + 1 {
            let layer_size = 1 << (depth_1 - k);
            let mut temp_w_init = vec![Scalar::ZERO; layer_size];
            let mut temp_s = vec![Scalar::ZERO; layer_size];

            for i in 0..layer_size {
                temp_w_init[i] = w_init_circuit_layers[c][k - 1].get_position(2 * i)
                    * w_init_circuit_layers[c][k - 1].get_position(2 * i + 1);
                temp_s[i] = s_circuit_layers[c][k - 1].get_position(2 * i)
                    * s_circuit_layers[c][k - 1].get_position(2 * i + 1);
            }

            w_init_circuit_layers[c].push(CircuitLayer::new(temp_w_init));
            s_circuit_layers[c].push(CircuitLayer::new(temp_s));
        }

        for k in 1..depth_2 + 1 {
            let layer_size = 1 << (depth_2 - k);
            let mut temp_w_update = vec![Scalar::ZERO; layer_size];
            let mut temp_r = vec![Scalar::ZERO; layer_size];
            for i in 0..layer_size {
                temp_w_update[i] = w_update_circuit_layers[c][k - 1].get_position(2 * i)
                    * w_update_circuit_layers[c][k - 1].get_position(2 * i + 1);
                temp_r[i] = r_circuit_layers[c][k - 1].get_position(2 * i)
                    * r_circuit_layers[c][k - 1].get_position(2 * i + 1);
            }

            w_update_circuit_layers[c].push(CircuitLayer::new(temp_w_update));
            r_circuit_layers[c].push(CircuitLayer::new(temp_r));
        }

        w_init_circuits.push(CircuitBinaryTree::new(w_init_circuit_layers[c].clone()));
        s_circuits.push(CircuitBinaryTree::new(s_circuit_layers[c].clone()));
        w_update_circuits.push(CircuitBinaryTree::new(w_update_circuit_layers[c].clone()));
        r_circuits.push(CircuitBinaryTree::new(r_circuit_layers[c].clone()));
    }
    BatchGrandProductCircuits::new(w_init_circuits, s_circuits, w_update_circuits, r_circuits)
}
