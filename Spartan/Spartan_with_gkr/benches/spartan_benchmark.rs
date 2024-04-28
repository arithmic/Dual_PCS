#![allow(non_snake_case)]
use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use channel::Channel;
use criterion::*;
use multilinear_kzg::common::{setup, VerificationKey};
use polynomial::MultPolynomial;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use traits::traits::Field;
use Spartan_with_gkr::{
    prover::er1cs::prove_sat,
    spartan_common::{eR1CSCommitments, eR1CStranscript},
    test::{construct_matrices, er1cs_commit},
    verifier::er1cs::verify_sat,
};
extern crate criterion;
fn spartan_prover_benchmark(c: &mut Criterion) {
    let toxic_waste = (0..16).into_par_iter().map(|_| Scalar::random()).collect();
    let (srs, _ver_key) = setup::<BlsCurve>(toxic_waste);
    for &idx in [10, 12, 16].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("spartan_prover_benchmark");
        group.plot_config(plot_config);
        let srs = srs.clone();
        let num_const = 1 << idx;
        let num_inputs = 10;
        let num_var = num_const - 1;
        let sparsity = 1;
        let (A, B, C, z, E, W, u, _PI) =
            construct_matrices(sparsity as usize, num_const, num_var as usize, num_inputs);

        let (er1cs_metadata, er1cs_commitments) = er1cs_commit(
            A.clone(),
            B.clone(),
            C.clone(),
            E.clone(),
            W.clone(),
            &srs,
            sparsity,
        );

        let mut channel = Channel::initialize_with_affine_point(
            &[
                er1cs_commitments.E.commitment.to_affine(),
                er1cs_commitments.W.commitment.to_affine(),
            ]
            .to_vec(),
        );

        let name = format!("Spartan Prover");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                prove_sat(
                    black_box(A.clone()),
                    black_box(B.clone()),
                    black_box(C.clone()),
                    &u,
                    black_box(&MultPolynomial::new(z.clone())),
                    black_box(&MultPolynomial::new(E.clone())),
                    black_box(&MultPolynomial::new(W.clone())),
                    black_box(er1cs_metadata.clone()),
                    black_box(&srs),
                    black_box(&mut channel),
                )
            });
        });
        group.finish();
    }
}
fn spartan_verifier_benchmark(c: &mut Criterion) {
    let toxic_waste = (0..16).into_par_iter().map(|_| Scalar::random()).collect();
    let (srs, ver_key) = setup::<BlsCurve>(toxic_waste);
    for &idx in [10, 12, 16].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("spartan_verifier_benchmark");
        group.plot_config(plot_config);
        let srs = srs.clone();
        let num_const = 1 << idx;
        let num_inputs = 10;
        let num_var = num_const - 1;
        let sparsity = 1;
        let (A, B, C, z, E, W, u, PI) =
            construct_matrices(sparsity as usize, num_const, num_var as usize, num_inputs);

        let (er1cs_metadata, er1cs_commitments) = er1cs_commit(
            A.clone(),
            B.clone(),
            C.clone(),
            E.clone(),
            W.clone(),
            &srs,
            sparsity,
        );

        let mut channel = Channel::initialize_with_affine_point(
            &[
                er1cs_commitments.E.commitment.to_affine(),
                er1cs_commitments.W.commitment.to_affine(),
            ]
            .to_vec(),
        );
        let er1cs_transcript = prove_sat(
            A,
            B,
            C,
            &u,
            &MultPolynomial::new(z),
            &MultPolynomial::new(E),
            &MultPolynomial::new(W),
            er1cs_metadata,
            &srs,
            &mut channel,
        );
        group.sample_size(10);

        let ver_key = ver_key.clone();
        let pi_indices: Vec<usize> = (0..1 << 5).collect();
        let name = format!("Spartan Verifier");
        group.bench_function(&name, move |b| {
            b.iter(|| {
                bench_verify_sat(
                    black_box(er1cs_transcript.clone()),
                    black_box(er1cs_commitments.clone()),
                    black_box(u),
                    black_box(MultPolynomial::new(PI.clone())),
                    black_box(pi_indices.clone()),
                    black_box(&ver_key),
                );
            });
        });
        group.finish();
    }
}
criterion_group!(
    benches,
    spartan_prover_benchmark,
    spartan_verifier_benchmark
);
criterion_main!(benches);

//--------Helper functions ----------
pub fn bench_verify_sat(
    er1cs_transcript: eR1CStranscript,
    er1cs_commitments: eR1CSCommitments,
    u: Scalar,
    PI: MultPolynomial,
    pi_indices: Vec<usize>,
    ver_key: &VerificationKey<BlsCurve>,
) {
    let mut channel = Channel::initialize_with_affine_point(
        [
            er1cs_commitments.E.commitment.to_affine(),
            er1cs_commitments.W.commitment.to_affine(),
        ]
        .as_ref(),
    );
    verify_sat(
        er1cs_transcript,
        er1cs_commitments,
        u,
        PI,
        pi_indices,
        &ver_key,
        &mut channel,
    );
}
