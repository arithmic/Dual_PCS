use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use criterion::{
    black_box, criterion_group, criterion_main, AxisScale, Criterion, PlotConfiguration,
};
use multilinear_kzg::{
    common::setup,
    prover::{commit, evaluate},
    verifier::verify,
};
use traits::traits::Field;

fn multilinear_kzg_setup_benchmark(c: &mut Criterion) {
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("multilinear_kzg_setup_benchmark");
        group.plot_config(plot_config);

        let toxic_waste: Vec<_> = (0..idx).map(|_| Scalar::random()).collect();

        let name = format!("Multilinear KZG Setup");
        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| setup::<BlsCurve>(black_box(toxic_waste.clone())));
        });
        group.finish();
    }
}
fn multilinear_kzg_commit_benchmark(c: &mut Criterion) {
    let toxic_waste: Vec<_> = (0..20).map(|_| Scalar::random()).collect();
    let setup = setup::<BlsCurve>(toxic_waste);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("multilinear_kzg_commit_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << idx;

        let (srs, _verification_key) = setup.clone();

        let poly: Vec<Scalar> = (0..degree).map(|_| <Scalar as Field>::random()).collect();
        let point: Vec<_> = (0..idx).map(|_| <Scalar as Field>::random()).collect();

        group.sample_size(10);
        let name = format!("Multilinear KZG  Prover");
        group.bench_function(&name, move |b| {
            b.iter(|| {
                evaluate(black_box(&poly), black_box(&point), black_box(&srs));
            });
        });
        group.finish();
    }
}

fn multilinear_kzg_prover_benchmark(c: &mut Criterion) {
    let toxic_waste: Vec<_> = (0..20).map(|_| Scalar::random()).collect();
    let setup = setup::<BlsCurve>(toxic_waste);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("multilinear_kzg_commit_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << idx;

        let (srs, _verification_key) = setup.clone();

        let poly: Vec<Scalar> = (0..degree).map(|_| <Scalar as Field>::random()).collect();

        group.sample_size(10);
        let name = format!("Multilinear KZG Commit");
        group.bench_function(&name, move |b| {
            b.iter(|| {
                commit(&black_box(poly.clone()), black_box(&srs));
            });
        });
        group.finish();
    }
}

fn multilinear_kzg_verifier_benchmark(c: &mut Criterion) {
    let toxic_waste: Vec<_> = (0..20).map(|_| Scalar::random()).collect();
    let setup = setup::<BlsCurve>(toxic_waste);
    for &idx in [15, 16, 17, 18, 19, 20].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("multilinear_kzg_verifier_benchmark");
        group.plot_config(plot_config);

        let degree = 1 << idx;

        let (srs, verification_key) = setup.clone();

        let poly: Vec<Scalar> = (0..degree).map(|_| <Scalar as Field>::random()).collect();
        let point: Vec<_> = (0..idx).map(|_| <Scalar as Field>::random()).collect();

        //Commit polynomial
        let commitment = commit(&poly, &srs);

        let proof = evaluate(&poly, &point, &srs);
        let name = format!("KZG Fourier Verifier");

        group.sample_size(10);
        group.bench_function(&name, move |b| {
            b.iter(|| {
                verify(
                    black_box(&commitment),
                    black_box(&proof),
                    black_box(&point),
                    black_box(&verification_key),
                );
            });
        });
        group.finish();
    }
}
criterion_group!(
    benches,
    multilinear_kzg_setup_benchmark,
    multilinear_kzg_commit_benchmark,
    multilinear_kzg_prover_benchmark,
    multilinear_kzg_verifier_benchmark,
);
criterion_main!(benches);
