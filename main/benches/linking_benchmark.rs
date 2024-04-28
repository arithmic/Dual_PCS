use bls381::scalar::Scalar;
use criterion::*;
use linking_prover::{commit::commit_witness, prover::linking_prover};
use linking_verifier::verifier::linking_verifier;
use setup::setup;
use traits::traits::Field;
extern crate bls381;
extern crate criterion;
extern crate traits;
fn linking_commitment_benchmark(c: &mut Criterion) {
    let degree_bound_setup = setup(1 << 15);
    for &iter in [7, 8, 9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("linking_commitment_benchmark");
        group.plot_config(plot_config);
        let degree = 1 << iter;
        let setup = degree_bound_setup.setup_for_specific_degree(degree);
        let witness = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();

        let name = format!("Linking Commitment");

        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| commit_witness(black_box(&witness), black_box(&setup)));
        });
        group.finish();
    }
}
fn linking_prover_benchmark(c: &mut Criterion) {
    let degree_bound_setup = setup(1 << 15);
    for &iter in [7, 8, 9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("linking_prover_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << iter;
        let setup = degree_bound_setup.setup_for_specific_degree(degree);
        let witness = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let witness_commitments = commit_witness(&witness, &setup);

        let name = format!("Linking Prover");

        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                linking_prover(
                    black_box(witness_commitments.clone()),
                    black_box(setup.clone()),
                )
            });
        });
        group.finish();
    }
}
fn linking_verifier_benchmark(c: &mut Criterion) {
    let degree_bound_setup = setup(1 << 16);
    for &iter in [7, 8, 9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("linking_verifier_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << iter;
        let setup = degree_bound_setup.setup_for_specific_degree(degree);
        let witness = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let witness_commitments = commit_witness(&witness, &setup);
        let linking_proof = linking_prover(witness_commitments, setup.clone());
        let name = format!("Linking Verifier");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| linking_verifier(black_box(linking_proof.clone()), black_box(setup.clone())));
        });
        group.finish();
    }
}

criterion_group!(
    benches,
    linking_commitment_benchmark,
    linking_prover_benchmark,
    linking_verifier_benchmark
);
criterion_main!(benches);
