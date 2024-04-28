use criterion::*;
use setup::setup;
extern crate criterion;
extern crate setup;
fn linking_setup_benchmark(c: &mut Criterion) {
    for &s in [7, 8, 9, 10, 11, 12, 13, 14, 15].iter() {
        let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
        let mut group = c.benchmark_group("linking_setup_benchmark");
        group.plot_config(plot_config);
        group.sample_size(10);
        let name = format!("Linking Setup");
        group.bench_function(&name, move |b| {
            b.iter(|| setup(1 << s));
        });
        group.finish();
    }
}
criterion_group!(benches, linking_setup_benchmark);
criterion_main!(benches);
