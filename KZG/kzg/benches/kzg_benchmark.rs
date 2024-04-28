use bls381::scalar::Scalar;
use criterion::*;
use kzg::{
    commit::kzg_commit, polynomial::Polynomial, prover::kzg_prover, setup::kzg_setup,
    verifier::kzg_verify,
};
use traits::traits::Field;
extern crate criterion;
fn kzg_setup_benchmark(c: &mut Criterion) {
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_setup_benchmark");
        group.plot_config(plot_config);

        let name = format!("KZG Setup");
        group.bench_function(&name, move |b| {
            b.iter(|| kzg_setup(1 << idx));
        });
        group.finish();
    }
}
fn kzg_commit_benchmark(c: &mut Criterion) {
    let setup = kzg_setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_commit_benchmark");
        group.plot_config(plot_config);
        let g1_powers = setup.0.clone()[0..1 << idx].to_vec();
        let evaluations = (0..1 << idx)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let name = format!("KZG Commit");
        group.sample_size(5);
        group.bench_function(&name, move |b| {
            b.iter(|| kzg_commit(black_box(evaluations.clone()), black_box(&g1_powers)));
        });
        group.finish();
    }
}
fn kzg_prover_benchmark(c: &mut Criterion) {
    let setup = kzg_setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_prover_benchmark");
        group.plot_config(plot_config);

        let polynomial = Polynomial::random(1 << idx);
        let z = Scalar::random();
        let g1_powers = setup.0[0..1 << idx].to_vec();
        let name = format!("KZG Open");
        group.sample_size(5);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                kzg_prover(
                    black_box(polynomial.clone()),
                    black_box(&g1_powers),
                    black_box(z),
                )
            });
        });
        group.finish();
    }
}
fn kzg_verify_benchmark(c: &mut Criterion) {
    let setup = kzg_setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_verify_benchmark");
        group.plot_config(plot_config);

        let evaluations: Vec<_> = (0..(1 << idx))
            .map(|_| <Scalar as Field>::random())
            .collect();
        let g1_powers = setup.0[0..1 << idx].to_vec();
        let (commitment, coeff) = kzg_commit(evaluations.clone(), &g1_powers);

        let polynomial = Polynomial { coeffs: coeff };
        let z = Scalar::random();
        let eval = polynomial.evaluate(&z);
        let verifier_key = setup.1;
        let proof = kzg_prover(polynomial.clone(), &g1_powers, z);

        group.sample_size(5);
        let name = format!("KZG Verify");
        group.bench_function(&name, move |b| {
            b.iter(|| {
                kzg_verify(
                    black_box(commitment),
                    black_box(z),
                    black_box(eval),
                    black_box(proof),
                    black_box(&verifier_key),
                )
            });
        });
        group.finish();
    }
}
criterion_group!(
    benches,
    kzg_setup_benchmark,
    kzg_commit_benchmark,
    kzg_prover_benchmark,
    kzg_verify_benchmark
);
criterion_main!(benches);
