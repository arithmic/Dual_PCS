use crate::grand_product_common::{GkrTranscript, OODFrame, OFFSET};
use bls381::scalar::Scalar;
use bls_curve::bls::{AffinePoint, ProjectivePoint};
use channel::Channel;
use evaluations::Evaluations;
use fft::{
    par_fft::{get_inv_twiddles, get_twiddles, par_eval_poly, par_interpolate_poly_with_offset},
    serial_fft::eval,
};
use kzg_fft::commit::kzg2commit;
use poly::PolyTable;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use trace::Trace;
use traits::traits::{Field, PrimeField};
mod evaluations;
pub mod grand_product;
mod poly;
pub mod trace;

pub fn gkr_prover(
    leaf_layers: &Vec<Vec<Scalar>>,
    final_layers: Vec<Scalar>,
    channel: &mut Channel,
    prover_key: &Vec<ProjectivePoint>,
) -> (GkrTranscript, Scalar, Trace, Vec<Scalar>) {
    let n_circuits = leaf_layers.len();
    let trace_length = leaf_layers[0].len();

    let trace = Trace::new(leaf_layers);
    let poly_table = PolyTable::new(&trace);

    let trace_extension = Evaluations::new(&poly_table);

    let trace_commits = trace.commit_columns(prover_key);
    channel.reseed_with_affine_point(&trace_commits.clone());

    let constraint_comp_coeff = channel.get_random_points(6 * n_circuits);

    let mut constraint_evaluations =
        trace_extension.evaluate_and_merge_constraints(&final_layers, &constraint_comp_coeff);

    let inv_twiddles = get_inv_twiddles(constraint_evaluations.len() as u32);
    par_interpolate_poly_with_offset(&mut constraint_evaluations, inv_twiddles, OFFSET);

    let composition_poly = constraint_evaluations.clone();

    let twiddles = get_twiddles(constraint_evaluations.len() as u32);
    par_eval_poly(&mut constraint_evaluations, &twiddles);

    let composition_poly_commit = kzg2commit(&constraint_evaluations, prover_key);
    channel.reseed_with_affine_point(&[composition_poly_commit].to_vec());

    //Draw a random point using Fiat-Shamir
    let z = channel.get_random_point();
    let g_trace = Scalar::get_root_of_unity(trace_length.trailing_zeros());
    let g_trace_z = g_trace * z;

    let current_frame = poly_table.evaluate_poly_table(z);

    let next_frame = poly_table.evaluate_poly_table(g_trace_z);

    let composition_frame = eval(&composition_poly, z);

    let ood_frame = OODFrame::new(current_frame, next_frame, composition_frame);
    ood_frame.reseed_with_ood_frame(channel);

    let poly_combiners = channel.get_random_points(n_circuits + 1);

    let (quotient_poly_commit_1, quotient_poly_commit_2) = commit_quotient_polys(
        &trace,
        &constraint_evaluations,
        poly_combiners,
        ood_frame.get_current_frame(),
        ood_frame.get_next_frame(),
        composition_frame,
        z,
        g_trace_z,
        prover_key,
    );
    (
        GkrTranscript::new(
            ood_frame,
            quotient_poly_commit_1,
            quotient_poly_commit_2,
            trace_commits,
            composition_poly_commit,
            trace_length,
            final_layers,
        ),
        z,
        trace,
        constraint_evaluations,
    )
}

fn commit_quotient_polys(
    trace: &Trace,
    composition_poly_evaluations: &Vec<Scalar>,
    poly_combiners: Vec<Scalar>,
    evaluation_at_z: &Vec<Scalar>,
    evaluation_at_gz: &Vec<Scalar>,
    composition_frame: Scalar,
    z: Scalar,
    g_trace_z: Scalar,
    prover_key: &Vec<ProjectivePoint>,
) -> (AffinePoint, AffinePoint) {
    let domain = trace.num_of_rows();
    let (poly1, poly2): (Vec<Scalar>, Vec<Scalar>) = (0..domain)
        .into_par_iter()
        .map(|idx| {
            let mut next = Scalar::ZERO;
            for inner_idx in 0..trace.num_of_columns() {
                if inner_idx % 2 != 0 {
                    next += poly_combiners[inner_idx / 2 + 1] * trace.cell_at(idx, inner_idx)
                }
            }
            let current = next + (poly_combiners[0] * composition_poly_evaluations[idx]);
            (current, next)
        })
        .unzip();

    let mut poly1_evaluation = Scalar::ZERO;
    let mut poly2_evaluation = Scalar::ZERO;
    for idx in 0..trace.num_of_columns() / 2 {
        poly1_evaluation += poly_combiners[idx + 1] * evaluation_at_z[idx];
        poly2_evaluation += poly_combiners[idx + 1] * evaluation_at_gz[idx];
    }
    poly1_evaluation += poly_combiners[0] * composition_frame;

    let g_domain = Scalar::get_root_of_unity(domain.trailing_zeros());
    let mut omega_powers = vec![Scalar::ONE; domain];
    for idx in 1..omega_powers.len() {
        omega_powers[idx] = omega_powers[idx - 1] * g_domain;
    }
    let (quotient_poly1, quotient_poly2): (Vec<Scalar>, Vec<Scalar>) = (0..domain)
        .into_par_iter()
        .map(|idx| {
            (
                (poly1[idx] - poly1_evaluation) * (omega_powers[idx] - z).invert().unwrap(),
                (poly2[idx] - poly2_evaluation) * (omega_powers[idx] - g_trace_z).invert().unwrap(),
            )
        })
        .unzip();
    (
        kzg2commit(&quotient_poly1, prover_key),
        kzg2commit(&quotient_poly2, prover_key),
    )
}
