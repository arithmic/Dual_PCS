#![allow(non_snake_case)]
use criterion::*;
use grand_product_with_air::grand_product_test::{
    construct_input, test_prover::prover, test_verifier::verifier,
};
use kzg_fft::setup::kzg2_setup;
extern crate criterion;

fn gp_air_prove_benchmark(c: &mut Criterion) {
    let setup = kzg2_setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("gp_air_prove_benchmark");
        group.plot_config(plot_config);
        let setup = setup.clone();
        let total_circuits = 1;
        let input_length = 1 << idx;
        let (A, B, C) = construct_input(total_circuits, input_length);

        let name = format!("GP_AIR_Prover");
        group.bench_function(&name, move |b| {
            b.iter(|| {
                prover(
                    black_box(&A),
                    black_box(&B),
                    black_box(&C),
                    black_box(total_circuits),
                    black_box(&setup),
                );
            });
        });
        group.finish();
    }
}
fn gp_air_verify_benchmark(c: &mut Criterion) {
    let setup = kzg2_setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("gp_air_verify_benchmark");
        group.plot_config(plot_config);
        let total_circuits = 1;
        let input_length = 1 << idx;
        let setup = setup.clone();

        let (A, B, C) = construct_input(total_circuits, input_length);
        let (commitments, gkr_transcript, evaluations, proof1, proof2) =
            prover(&A, &B, &C, total_circuits, &setup);

        let name = format!("GP_AIR_Verifier");
        group.bench_function(&name, move |b| {
            b.iter(|| {
                verifier(
                    black_box(gkr_transcript.clone()),
                    black_box(commitments.clone()),
                    black_box(evaluations.clone()),
                    black_box(&setup),
                    black_box(proof1.clone()),
                    black_box(proof2.clone()),
                );
            });
        });
        group.finish();
    }
}

criterion_group!(benches, gp_air_prove_benchmark, gp_air_verify_benchmark);
criterion_main!(benches);
