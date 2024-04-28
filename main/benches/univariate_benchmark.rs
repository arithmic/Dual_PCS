use bls381::scalar::Scalar;
use criterion::*;
use linking_prover::commit::commit_witness;
use setup::setup;
use traits::traits::Field;

use uni_multi_prover::univariate::univariate_evaluation_prover;
use uni_multi_verifier::univariate::univariate_evaluation_verifier;
extern crate bls381;
extern crate criterion;
extern crate traits;

fn univariate_prover_benchmark(c: &mut Criterion) {
    let degree_bound_setup = setup(1 << 15);
    for &iter in [7, 8, 9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("univariate_prover_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << iter;
        let setup = degree_bound_setup.setup_for_specific_degree(degree);
        let witness = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let witness_commitments = commit_witness(&witness, &setup);

        let name = format!("Univariate Prover");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                univariate_evaluation_prover(
                    black_box(&setup),
                    black_box(setup.tau_1.clone()),
                    black_box(
                        [
                            witness_commitments.commitment_to_univariate,
                            setup.tau_ipp[0],
                            witness_commitments.commitment_to_univariate,
                        ]
                        .to_vec(),
                    ),
                    black_box(witness_commitments.g2_power_poly.clone()),
                    black_box(witness_commitments.univariate_polynomial.clone()),
                )
            });
        });
        group.finish();
    }
}
fn univariate_verifier_benchmark(c: &mut Criterion) {
    let degree_bound_setup = setup(1 << 15);
    for &iter in [7, 8, 9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("univariate_verifier_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << iter;
        let setup = degree_bound_setup.setup_for_specific_degree(degree);
        let witness = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let witness_commitments = commit_witness(&witness, &setup);
        let univariate_evaluation_proof = univariate_evaluation_prover(
            &setup,
            setup.tau_1.clone(),
            [
                witness_commitments.commitment_to_univariate,
                setup.tau_ipp[0],
                witness_commitments.commitment_to_univariate,
            ]
            .to_vec(),
            witness_commitments.g2_power_poly,
            witness_commitments.univariate_polynomial,
        );

        let name = format!("Univariate Verifier");
        group.sample_size(10);

        group.bench_function(&name, move |b| {
            b.iter(|| {
                univariate_evaluation_verifier(
                    black_box(univariate_evaluation_proof.clone()),
                    black_box(&setup),
                )
            });
        });
        group.finish();
    }
}

criterion_group!(
    benches,
    univariate_prover_benchmark,
    univariate_verifier_benchmark
);
criterion_main!(benches);
