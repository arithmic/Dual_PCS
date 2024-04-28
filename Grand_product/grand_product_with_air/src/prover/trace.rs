use crate::grand_product_common::evaluate_constraints;
use bls381::scalar::Scalar;
use bls_curve::bls::{AffinePoint, ProjectivePoint};
use kzg_fft::commit;
use table::Table;

pub struct Trace {
    trace: Table<Scalar>,
}
#[allow(dead_code)]
impl Trace {
    pub fn new(leaf_layers: &Vec<Vec<Scalar>>) -> Self {
        (0..leaf_layers.len()).for_each(|idx| {
            assert!(
                leaf_layers[idx].len().is_power_of_two(),
                "Length of each leaf layer should be power of 2"
            );
        });
        let n_circuits = leaf_layers.len();
        let no_columns = 2 * n_circuits;
        let no_rows = leaf_layers[0].len();
        let mut trace = Table::new(vec![vec![Scalar::ZERO; no_columns]; no_rows]);
        for j in 0..n_circuits {
            trace.update_column_at(j * 2, &leaf_layers[j]);
            trace.set_cell(0, (j * 2) + 1, leaf_layers[j][0]);
        }
        for row in 1..no_rows {
            for column in 0..n_circuits {
                let value =
                    trace.cell_at(row - 1, (2 * column) + 1) * trace.cell_at(row, 2 * column);
                trace.set_cell(row, (2 * column) + 1, value);
            }
        }
        Self { trace }
    }
    pub fn get_trace_table(&self) -> &Table<Scalar> {
        &self.trace
    }
    pub fn num_of_columns(&self) -> usize {
        self.trace.num_of_columns()
    }
    pub fn num_of_rows(&self) -> usize {
        self.trace.num_of_rows()
    }
    pub fn column_at(&self, idx: usize) -> Vec<Scalar> {
        self.trace.column_at(idx)
    }
    pub fn row_at(&self, idx: usize) -> Vec<Scalar> {
        self.trace.row_at(idx)
    }
    pub fn cell_at(&self, row_index: usize, column_index: usize) -> Scalar {
        self.trace.cell_at(row_index, column_index)
    }
    ///We commit only even number index columns of the trace
    /// since verifier can compute commits of odd number indexed columns
    pub fn commit_columns(&self, prover_key: &Vec<ProjectivePoint>) -> Vec<AffinePoint> {
        (0..self.num_of_columns())
            .filter_map(|idx| {
                if idx % 2 != 0 {
                    Some(commit::kzg2commit(&self.column_at(idx), prover_key))
                } else {
                    None
                }
            })
            .collect::<Vec<AffinePoint>>()
    }
    pub fn debug_trace(&self, final_layers: &Vec<Scalar>, n_circuits: usize) {
        for idx in 0..self.num_of_rows() {
            let current = &self.row_at(idx);
            let next = &self.row_at((idx + 1) % self.num_of_rows());
            let mut evaluations = vec![Scalar::ONE; 6 * n_circuits];
            evaluate_constraints(current, next, n_circuits, final_layers, &mut evaluations);
            if idx == 0 {
                assert_eq!(
                    evaluations[0..n_circuits],
                    vec![Scalar::ZERO; n_circuits],
                    "assertion failed at idx {}",
                    idx
                );
            }
            if idx != self.num_of_rows() - 1 {
                assert_eq!(
                    evaluations[n_circuits..2 * n_circuits],
                    vec![Scalar::ZERO; n_circuits],
                    "assertion failed at idx {}",
                    idx
                )
            }
            if idx == self.num_of_rows() - 1 {
                assert_eq!(
                    evaluations[2 * n_circuits..3 * n_circuits],
                    vec![Scalar::ZERO; n_circuits],
                    "assertion failed at idx {}",
                    idx
                )
            }
        }
    }
}
