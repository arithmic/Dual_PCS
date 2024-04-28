use bls381::scalar::Scalar;
use criterion::*;
use linking_prover::commit::commit_witness;
use setup::setup;
use traits::traits::Field;

use uni_multi_prover::multilinear::multilinear_evaluation_prover;
use uni_multi_verifier::multilinear::multilinear_evaluation_verifier;
extern crate bls381;
extern crate criterion;
extern crate traits;

fn multilinear_prover_benchmark(c: &mut Criterion) {
    let degree_bound_setup = setup(1 << 15);
    for &iter in [7, 8, 9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("multilinear_prover_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << iter;
        let setup = degree_bound_setup.setup_for_specific_degree(degree);
        let witness = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let witness_commitments = commit_witness(&witness, &setup);

        let name = format!("multilinear Prover");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                multilinear_evaluation_prover(
                    black_box(&setup),
                    black_box(setup.tau_1.clone()),
                    black_box(
                        [
                            witness_commitments.commitment_to_witness,
                            setup.tau_ipp[0],
                            witness_commitments.commitment_to_witness,
                        ]
                        .to_vec(),
                    ),
                    black_box(witness_commitments.g2_power_witness.clone()),
                    black_box(witness.clone()),
                )
            });
        });
        group.finish();
    }
}
fn multilinear_verifier_benchmark(c: &mut Criterion) {
    let degree_bound_setup = setup(1 << 15);
    for &iter in [7, 8, 9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("multilinear_verifier_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << iter;
        let setup = degree_bound_setup.setup_for_specific_degree(degree);
        let witness = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let witness_commitments = commit_witness(&witness, &setup);
        let multilinear_evaluation_proof = multilinear_evaluation_prover(
            &setup,
            setup.tau_1.clone(),
            [
                witness_commitments.commitment_to_witness,
                setup.tau_ipp[0],
                witness_commitments.commitment_to_witness,
            ]
            .to_vec(),
            witness_commitments.g2_power_witness,
            witness,
        );
        let name = format!("multilinear Verifier");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                multilinear_evaluation_verifier(
                    black_box(multilinear_evaluation_proof.clone()),
                    black_box(&setup),
                )
            });
        });
        group.finish();
    }
}

criterion_group!(
    benches,
    multilinear_prover_benchmark,
    multilinear_verifier_benchmark
);
criterion_main!(benches);
