use crate::grand_product_common::{evaluate_constraints, merge_constraints, GkrTranscript};
use bls381::scalar::Scalar;
use bls_curve::{
    bls::{BlsCurve, ProjectivePoint},
    fp2bls::G2BlsCurve,
};
use curve_traits::AffinePoint;
use pairing::pairing_traits::Pairing;
use traits::traits::PrimeField;

pub fn gkr_verifier(
    gkr_transcript: GkrTranscript,
    constraint_comp_coeffs: Vec<Scalar>,
    z: Scalar,
    poly_combiners: Vec<Scalar>,
    n_circuits: usize,
    trace_length: usize,
    verifier_key: AffinePoint<G2BlsCurve>,
    leaf_layer1_current: Vec<Scalar>,
    leaf_layer1_next: Vec<Scalar>,
) {
    let final_layers_1 = gkr_transcript.get_final_layers();
    let ood_frame = gkr_transcript.get_odd_frame();
    let current_frame = ood_frame.get_current_frame();
    let next_frame = ood_frame.get_next_frame();
    let mut current = vec![Scalar::ZERO; 2 * n_circuits];
    for idx in 0..n_circuits {
        current[2 * idx] = leaf_layer1_current[idx];
        current[(2 * idx) + 1] = current_frame[idx];
    }
    let mut next = vec![Scalar::ZERO; 2 * n_circuits];
    for idx in 0..n_circuits {
        next[2 * idx] = leaf_layer1_next[idx];
        next[(2 * idx) + 1] = next_frame[idx];
    }
    let mut evaluations = vec![Scalar::ZERO; 3 * n_circuits];
    evaluate_constraints(
        &current,
        &next,
        n_circuits,
        final_layers_1,
        &mut evaluations,
    );
    let g_trace = Scalar::get_root_of_unity(trace_length.trailing_zeros());
    let g_trace_z = g_trace * z;
    let evaluation = merge_constraints(
        &evaluations,
        z,
        constraint_comp_coeffs,
        trace_length,
        g_trace,
        n_circuits,
    );
    let composition_frame = gkr_transcript.get_odd_frame().get_composition_frame();
    assert_eq!(evaluation, composition_frame, "OOD check failed");
    let trace_commits = gkr_transcript.get_trace_commits();
    let composition_poly_commit = gkr_transcript.composition_poly_commit();

    //Generator of G1
    let g1 = <BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;
    //Generator of G2
    let g2 = <G2BlsCurve as curve_traits::CurveArithmetic>::ProjectivePoint::GENERATOR;
    let h_g2 = verifier_key.to_projective() + g2.mul(-z);
    let h_g2_z = verifier_key.to_projective() + g2.mul(-g_trace_z);

    let evaluation_at_z = (0..current_frame.len()).fold(Scalar::ZERO, |acc, idx| {
        acc + (poly_combiners[idx + 1] * current_frame[idx])
    }) + (poly_combiners[0] * composition_frame);
    let evaluation_at_g_z = (0..next_frame.len()).fold(Scalar::ZERO, |acc, idx| {
        acc + (poly_combiners[idx + 1] * next_frame[idx])
    });
    let commitment2 = (0..trace_commits.len()).fold(ProjectivePoint::IDENTITY, |acc, idx| {
        acc + trace_commits[idx]
            .to_projective()
            .mul(poly_combiners[idx + 1])
    });
    let commitment1 = commitment2
        + composition_poly_commit
            .to_projective()
            .mul(poly_combiners[0]);
    let (quotient_poly_commit_1, quotient_poly_commit_2) =
        gkr_transcript.get_quotient_poly_commits();

    let pairing1 = BlsCurve::tate_pairing(
        (commitment1 + g1.mul(-evaluation_at_z)).to_affine(),
        g2.to_affine(),
    );
    let pairing2 = BlsCurve::tate_pairing(quotient_poly_commit_1, h_g2.to_affine());
    assert_eq!(pairing1, pairing2, "pairing check 1 failed");

    let pairing1 = BlsCurve::tate_pairing(
        (commitment2 + g1.mul(-evaluation_at_g_z)).to_affine(),
        g2.to_affine(),
    );
    let pairing2 = BlsCurve::tate_pairing(quotient_poly_commit_2, h_g2_z.to_affine());
    assert_eq!(pairing1, pairing2, "pairing check 2 failed");
}
