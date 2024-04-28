#![allow(non_snake_case)]
use bls381::scalar::Scalar;
use criterion::*;
use grand_product_with_gkr::grand_product_test::{
    construct_input, test_prover::prover, test_verifier::verifier,
};
use multilinear_kzg::common::setup;
use traits::traits::Field;
extern crate criterion;
fn gp_gkr_prove_benchmark(c: &mut Criterion) {
    let toxic_waste = (0..20).map(|_| Scalar::random()).collect::<Vec<Scalar>>();
    let (srs, _ver_key) = setup(toxic_waste);

    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("gp_gkr_prove_benchmark");
        group.plot_config(plot_config);
        let srs = srs.clone();
        let total_circuits = 1;
        let input_length = 1 << idx;
        let (A, B, C) = construct_input(total_circuits, input_length);

        let name = format!("GP_GKR_Prover");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                prover(
                    black_box(&A),
                    black_box(&B),
                    black_box(&C),
                    black_box(total_circuits),
                    black_box(&srs),
                );
            });
        });
        group.finish();
    }
}
fn gp_gkr_verify_benchmark(c: &mut Criterion) {
    let toxic_waste = (0..20).map(|_| Scalar::random()).collect::<Vec<Scalar>>();
    let (srs, ver_key) = setup(toxic_waste);

    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("gp_gkr_verify_benchmark");
        group.plot_config(plot_config);
        let total_circuits = 1;
        let input_length = 1 << idx;
        let (A, B, C) = construct_input(total_circuits, input_length);
        let ver_key = ver_key.clone();

        let name = format!("GP_GKR_Verifier");
        let (commitments, gkr_transcript, evaluation_proof, evaluations) =
            prover(&A, &B, &C, total_circuits, &srs);

        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                verifier(
                    black_box(gkr_transcript.clone()),
                    black_box(commitments.clone()),
                    black_box(evaluation_proof.clone()),
                    black_box(evaluations.clone()),
                    black_box(ver_key.clone()),
                    black_box(total_circuits),
                );
            });
        });
        group.finish();
    }
}

criterion_group!(benches, gp_gkr_prove_benchmark, gp_gkr_verify_benchmark);
criterion_main!(benches);
