use super::poly::PolyTable;
use crate::grand_product_common::{evaluate_constraints, merge_constraints, OFFSET};
use bls381::scalar::Scalar;
use fft::par_fft::{get_twiddles, par_eval_poly_with_offset};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use table::Table;
use traits::traits::PrimeField;
pub(crate) struct Evaluations {
    evaluations: Table<Scalar>,
}
#[allow(dead_code)]
impl Evaluations {
    pub fn new(poly_table: &PolyTable) -> Self {
        let no_columns = poly_table.num_of_columns();
        let no_rows = poly_table.num_of_rows();
        let mut evaluations = Table::new(vec![vec![Scalar::ZERO; no_columns]; no_rows]);
        let twiddles = get_twiddles(no_rows as u32);
        (0..no_columns).for_each(|col| {
            let column =
                par_eval_poly_with_offset(&poly_table.column_at(col), &twiddles, OFFSET, 1);
            evaluations.update_column_at(col, &column)
        });
        Self { evaluations }
    }
    pub fn get_evaluations(&self) -> &Table<Scalar> {
        &self.evaluations
    }
    pub fn num_of_columns(&self) -> usize {
        self.evaluations.num_of_columns()
    }
    pub fn num_of_rows(&self) -> usize {
        self.evaluations.num_of_rows()
    }
    pub fn column_at(&self, idx: usize) -> Vec<Scalar> {
        self.evaluations.column_at(idx)
    }
    pub fn row_at(&self, idx: usize) -> Vec<Scalar> {
        self.evaluations.row_at(idx)
    }
    pub fn cell_at(&self, row_index: usize, column_index: usize) -> Scalar {
        self.evaluations.cell_at(row_index, column_index)
    }
    pub fn evaluate_and_merge_constraints(
        &self,
        final_layers_1: &Vec<Scalar>,
        constraint_comp_coeffs: &Vec<Scalar>,
    ) -> Vec<Scalar> {
        let n_circuits = final_layers_1.len();
        let domain_size = self.num_of_rows();
        let g_domain = Scalar::get_root_of_unity(domain_size.trailing_zeros());

        let mut x_arr = vec![Scalar::ZERO; domain_size];
        x_arr[0] = OFFSET;
        for i in 1..domain_size {
            x_arr[i] = x_arr[i - 1] * g_domain;
        }
        let evaluations = (0..domain_size)
            .into_par_iter()
            .map(|row| {
                let current = self.row_at(row);
                let next = self.row_at((row + 1) % domain_size);
                let mut evaluations = vec![Scalar::ZERO; 3 * n_circuits];
                evaluate_constraints(
                    &current,
                    &next,
                    n_circuits,
                    final_layers_1,
                    &mut evaluations,
                );
                merge_constraints(
                    &evaluations,
                    x_arr[row],
                    constraint_comp_coeffs.to_vec(),
                    domain_size,
                    g_domain,
                    n_circuits,
                )
            })
            .collect::<Vec<Scalar>>();
        evaluations
    }
}
