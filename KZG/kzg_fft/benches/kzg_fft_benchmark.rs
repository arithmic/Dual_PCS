use bls381::scalar::Scalar;
use channel::Channel;
use criterion::*;
use kzg_fft::{
    commit::kzg2commit, prover::kzg2_prover, setup::kzg2_setup, verifier::kzg2_verifier,
};
use traits::traits::Field;
extern crate bls381;
extern crate channel;
extern crate criterion;
extern crate helper;
extern crate kzg_fft;
extern crate traits;
fn kzg_fft_setup_benchmark(c: &mut Criterion) {
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_fft_setup_benchmark");
        group.plot_config(plot_config);

        let name = format!("KZG-FFT Setup");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| kzg2_setup(1 << idx));
        });
        group.finish();
    }
}
fn kzg_fft_commit_benchmark(c: &mut Criterion) {
    let setup = kzg2_setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_fft_commit_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << idx;

        let evaluations = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let prover_key = setup.get_setup(degree).prover_key;

        let name = format!("KZG-FFT Prover");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| kzg2commit(black_box(&evaluations), black_box(&prover_key)));
        });
        group.finish();
    }
}
fn kzg_fft_prover_benchmark(c: &mut Criterion) {
    let degree_bound_setup = kzg2_setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_fft_prover_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << idx;

        let evaluations = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();

        //Commit polynomial
        let setup = degree_bound_setup.get_setup(degree);

        let commitment_to_evaluations_of_f = kzg2commit(&evaluations, &setup.prover_key);

        //Initialize channel
        let mut channel =
            Channel::initialize_with_affine_point(&[commitment_to_evaluations_of_f].as_ref());

        let z = channel.get_random_point();

        group.sample_size(10);
        let name = format!("KZG-FFT Prover");
        group.bench_function(&name, move |b| {
            b.iter(|| {
                kzg2_prover(
                    black_box(&evaluations),
                    black_box(z),
                    black_box(&setup.prover_key),
                )
            });
        });
        group.finish();
    }
}
fn kzg_fft_verifier_benchmark(c: &mut Criterion) {
    let degree_bound_setup = kzg2_setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_fft_verifier_benchmark");
        group.plot_config(plot_config);
        let degree = 1 << idx;

        let evaluations = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();

        //Commit polynomial
        let setup = degree_bound_setup.get_setup(degree);

        let commitment_to_evaluations_of_f = kzg2commit(&evaluations, &setup.prover_key);

        //Initialize channel
        let mut channel =
            Channel::initialize_with_affine_point(&[commitment_to_evaluations_of_f].as_ref());

        let z = channel.get_random_point();

        let evaluation_proof = kzg2_prover(&evaluations, z, &setup.prover_key);

        let mut channel =
            Channel::initialize_with_affine_point(&[commitment_to_evaluations_of_f].as_ref());
        let z = channel.get_random_point();

        let name = format!("KZG-FFT Verifier");
        group.sample_size(10);

        group.bench_function(&name, move |b| {
            b.iter(|| {
                kzg2_verifier(
                    black_box(setup.verifier_key),
                    black_box(evaluation_proof.clone()),
                    black_box(commitment_to_evaluations_of_f),
                    black_box(z),
                )
            });
        });
        group.finish();
    }
}
criterion_group!(
    benches,
    kzg_fft_setup_benchmark,
    kzg_fft_commit_benchmark,
    kzg_fft_prover_benchmark,
    kzg_fft_verifier_benchmark
);
criterion_main!(benches);
