#![allow(non_snake_case)]
use bls381::scalar::Scalar;
use channel::Channel;
use criterion::*;
use kzg_fourier_multilinear::{common::KZGFourierDegreeBoundSetup, setup::multilinearkzg2setup};
use polynomial::MultPolynomial;
use Spartan_with_air::{
    preprocessing::SparseRep,
    prover::er1cs::prove_sat,
    spartan_common::{eR1CSCommitments, eR1CSmetadata, eR1CStranscript},
    test::{construct_matrices, er1cs_commit},
    verifier::er1cs::verify_sat,
};
extern crate criterion;
fn spartan_prover_benchmark(c: &mut Criterion) {
    let setup = multilinearkzg2setup(1 << 15);
    for &idx in [9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("spartan_prover_benchmark");
        group.plot_config(plot_config);

        let num_const = 1 << idx;
        let num_inputs = 10;
        let num_var = num_const - 1;
        let sparsity = 1;
        let (A, B, C, z, E, W, u, _PI) =
            construct_matrices(sparsity as usize, num_const, num_var as usize, num_inputs);
        let (Z, E, W, er1cs_metadata, _er1cs_commitments) = er1cs_commit(
            A.clone(),
            B.clone(),
            C.clone(),
            z,
            E,
            W,
            u,
            &setup,
            sparsity as usize,
        );

        let setup = setup.clone();
        let name = format!("Spartan Prover");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                bench_prove_sat(
                    black_box(A.clone()),
                    black_box(B.clone()),
                    black_box(C.clone()),
                    black_box(&u),
                    black_box(&Z),
                    black_box(&E),
                    black_box(&W),
                    black_box(er1cs_metadata.clone()),
                    black_box(&setup),
                )
            });
        });
        group.finish();
    }
}

fn spartan_verifier_benchmark(c: &mut Criterion) {
    let setup = multilinearkzg2setup(1 << 15);
    for &idx in [9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("spartan_verifier_benchmark");
        group.plot_config(plot_config);

        let num_const = 1 << idx;
        let num_inputs = 10;
        let num_var = num_const - 1;
        let sparsity = 1;
        let (A, B, C, z, E, W, u, PI) =
            construct_matrices(sparsity as usize, num_const, num_var as usize, num_inputs);

        let (Z, E, W, er1cs_metadata, er1cs_commitments) = er1cs_commit(
            A.clone(),
            B.clone(),
            C.clone(),
            z,
            E,
            W,
            u,
            &setup,
            sparsity as usize,
        );

        let mut channel = Channel::initialize_with_scalar(&[Scalar::ONE]);
        let setup = setup.clone();
        let er1cs_transcript = prove_sat(
            A,
            B,
            C,
            &u,
            &Z,
            &E,
            &W,
            er1cs_metadata,
            &mut channel,
            &setup,
        );

        let pi_indices: Vec<usize> = (0..1 << 5).collect();
        let name = format!("Spartan Verifier");
        group.bench_function(&name, move |b| {
            b.iter(|| {
                bench_verify_sat(
                    black_box(er1cs_transcript.clone()),
                    black_box(u),
                    black_box(MultPolynomial::new(PI.clone())),
                    black_box(pi_indices.clone()),
                    black_box(&setup),
                    black_box(er1cs_commitments.clone()),
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
//--------------------------Helper functions--------------------------
fn bench_prove_sat(
    A: SparseRep,
    B: SparseRep,
    C: SparseRep,
    u: &Scalar,
    z: &MultPolynomial,
    E: &MultPolynomial,
    W: &MultPolynomial,
    metadatas: eR1CSmetadata,
    setup: &KZGFourierDegreeBoundSetup,
) {
    let mut channel = Channel::initialize_with_scalar(&[Scalar::ONE]);
    prove_sat(A, B, C, u, z, E, W, metadatas, &mut channel, setup);
}
fn bench_verify_sat(
    er1cs_transcript: eR1CStranscript,
    u: Scalar,
    PI: MultPolynomial,
    pi_indices: Vec<usize>,
    setup: &KZGFourierDegreeBoundSetup,
    eR1CSCommitments: eR1CSCommitments,
) {
    let mut channel = Channel::initialize_with_scalar(&[Scalar::ONE]);
    verify_sat(
        er1cs_transcript,
        u,
        PI,
        pi_indices,
        &setup,
        eR1CSCommitments,
        &mut channel,
    );
}
