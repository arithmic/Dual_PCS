use super::trace::Trace;
use bls381::scalar::Scalar;
use fft::{
    par_fft::{get_inv_twiddles, par_interpolate_poly},
    serial_fft::eval,
};
use table::Table;

pub(crate) struct PolyTable {
    poly_table: Table<Scalar>,
}
#[allow(dead_code)]
impl PolyTable {
    pub fn new(trace: &Trace) -> Self {
        let no_columns = trace.num_of_columns();
        let no_rows = trace.num_of_rows();
        let mut poly_table = Table::new(vec![vec![Scalar::ZERO; no_columns]; no_rows]);
        let inv_twiddles = get_inv_twiddles(no_rows as u32);
        (0..no_columns).for_each(|col| {
            let mut column = trace.column_at(col);
            par_interpolate_poly(&mut column, inv_twiddles.clone());
            poly_table.update_column_at(col, &column)
        });
        Self { poly_table }
    }
    pub fn get_poly_table(&self) -> &Table<Scalar> {
        &self.poly_table
    }
    pub fn num_of_columns(&self) -> usize {
        self.poly_table.num_of_columns()
    }
    pub fn num_of_rows(&self) -> usize {
        self.poly_table.num_of_rows()
    }
    pub fn column_at(&self, idx: usize) -> Vec<Scalar> {
        self.poly_table.column_at(idx)
    }
    pub fn row_at(&self, idx: usize) -> Vec<Scalar> {
        self.poly_table.row_at(idx)
    }
    pub fn cell_at(&self, row_index: usize, column_index: usize) -> Scalar {
        self.poly_table.cell_at(row_index, column_index)
    }
    pub fn evaluate_poly_table(&self, point: Scalar) -> Vec<Scalar> {
        let mut evaluations = Vec::new();
        (0..self.num_of_columns()).for_each(|idx| {
            if idx % 2 != 0 {
                evaluations.push(eval(&self.column_at(idx), point))
            }
        });
        evaluations
    }
}
