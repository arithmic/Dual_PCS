use bls381::scalar::Scalar;
use channel::Channel;
use criterion::*;
use kzg_fft::commit;

use kzg_fourier_multilinear::{
    prover::prover, setup::multilinearkzg2setup, verifier::bench_verifier,
};
use traits::traits::Field;
extern crate bls381;
extern crate channel;
extern crate criterion;
extern crate helper;
extern crate kzg_fft;
extern crate kzg_fourier_multilinear;
extern crate traits;
fn kzg_fourier_setup_benchmark(c: &mut Criterion) {
    for &s in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_fourier_setup_benchmark");
        group.plot_config(plot_config);
        let name = format!("KZG Fourier Setup");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| multilinearkzg2setup(1 << s));
        });
        group.finish();
    }
}
fn kzg_fourier_prover_benchmark(c: &mut Criterion) {
    let setup = multilinearkzg2setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_fourier_prover_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << idx;

        let poly = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let setup = setup.clone();
        let prover_key = setup.get_setup(poly.len()).setup.prover_key;
        //Commit polynomial
        let commit_f = commit::kzg2commit(&poly, &prover_key);

        //Initialize channel
        let mut channel = Channel::initialize_with_affine_point(&[commit_f].as_ref());
        let random_points = (0..idx).map(|_| Scalar::random()).collect::<Vec<Scalar>>();

        group.sample_size(10);
        let name = format!("KZG Fourier Prover");
        group.bench_function(&name, move |b| {
            b.iter(|| {
                prover(
                    black_box(&poly),
                    black_box(random_points.clone()),
                    black_box(&setup.clone()),
                    black_box(&mut channel),
                )
            });
        });
        group.finish();
    }
}

fn kzg_fourier_verifier_benchmark(c: &mut Criterion) {
    let setup = multilinearkzg2setup(1 << 20);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("kzg_fourier_verifier_benchmark");
        group.plot_config(plot_config);
        let degree = 1 << idx;

        let poly = (0..degree)
            .map(|_| Scalar::random())
            .collect::<Vec<Scalar>>();
        let setup = setup.clone();
        //Commit polynomial
        let commit_f = commit::kzg2commit(&poly, &setup.setup.get_setup(poly.len()).prover_key);

        //Initialize channel
        let mut channel = Channel::initialize_with_affine_point(&[commit_f].as_ref());
        let random_points = (0..idx).map(|_| Scalar::random()).collect::<Vec<Scalar>>();

        let proof = prover(&poly, random_points.clone(), &setup, &mut channel);
        let name = format!("KZG Fourier Verifier");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                bench_verifier(
                    black_box(proof.clone()),
                    black_box(random_points.clone()),
                    black_box(setup.setup.get_setup(degree).verifier_key.clone()),
                    black_box(setup.get_setup(degree).phi_commitment.clone()),
                    black_box(&commit_f),
                )
            });
        });
        group.finish();
    }
}
criterion_group!(
    benches,
    kzg_fourier_setup_benchmark,
    kzg_fourier_prover_benchmark,
    kzg_fourier_verifier_benchmark
);
criterion_main!(benches);
