use crate::bit_reverse::{par_permute, permute, permute_index};
use crypto_bigint::U256;
use rayon::{
    self,
    prelude::{IntoParallelIterator, ParallelIterator}, iter::IntoParallelRefIterator,
};
use traits::traits::{ExtensionOf, Field, PrimeField};
use utils::uninit_vector;

pub fn get_twiddles<F: PrimeField>(N: u32) -> Vec<F> {
    let bits = 32 - N.leading_zeros() - 1;
    let omega_n = F::get_root_of_unity(bits);
    let mut omega: Vec<F> = (0..N >> 1)
        .into_par_iter()
        .map(|i| omega_n.power_by([i as u64]))
        .collect();
    par_permute(&mut omega);
    omega
}

pub fn get_inv_twiddles<F: PrimeField>(N: u32) -> Vec<F> {
    let bits = 32 - N.leading_zeros() - 1;
    let omega_n = F::invert(F::get_root_of_unity(bits)).unwrap();
    let mut omega: Vec<F> = (0..N >> 1)
        .into_par_iter()
        .map(|i| omega_n.power_by([i as u64]))
        .collect();
    par_permute(&mut omega);
    omega
}

// Code for serial FFT

// POLYNOMIAL EVALUATION
// ================================================================================================
/// Evaluates polynomial `p` in-place over the domain of length `p.len()` in the field specified
/// by `B` using the FFT algorithm.
pub fn evaluate_poly<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField>(
    p: &mut Vec<F>,
    twiddles: &[P],
) {
    fft_in_place(p, twiddles, 1, 1, 0);
    par_permute(p);
}
/// Evaluates polynomial `p` over the domain of length `p.len()` * `blowup_factor` shifted by
/// `domain_offset` in the field specified `B` using the FFT algorithm and returns the result.
pub fn evaluate_poly_with_offset<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField>(
    p: &[F],
    twiddles: &[P],
    domain_offset: P,
    blowup_factor: usize,
) -> Vec<F> {
    let domain_size = p.len() * blowup_factor;
    let g = P::get_root_of_unity(log2(domain_size));

    let mut result = unsafe { uninit_vector(domain_size) };
    result
        .as_mut_slice()
        .chunks_mut(p.len())
        .enumerate()
        .for_each(|(i, chunk)| {
            let idx = permute_index(blowup_factor, i) as u64;
            let offset = g.power_by(&U256::from(idx).to_words()) * domain_offset;
            let mut factor = P::ONE;
            for (d, c) in chunk.iter_mut().zip(p.iter()) {
                *d = (*c).mul_base(factor);
                factor = factor * offset;
            }
            fft_in_place(chunk, twiddles, 1, 1, 0);
        });
    permute(&mut result);
    result
}
pub(super) fn fft_in_place<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField>(
    values: &mut [F],
    twiddles: &[P],
    count: usize,
    stride: usize,
    offset: usize,
) {
    let size = values.len() / stride;
    debug_assert!(size.is_power_of_two());
    debug_assert!(offset < stride);
    debug_assert_eq!(values.len() % size, 0);
    // Keep recursing until size is 2

    const MAX_LOOP: usize = 256;
    if size > 2 {
        if stride == count && count < MAX_LOOP {
            fft_in_place(values, twiddles, 2 * count, 2 * stride, offset);
        } else {
            fft_in_place(values, twiddles, count, 2 * stride, offset);
            fft_in_place(values, twiddles, count, 2 * stride, offset + stride);
        }
    }

    for offset in offset..(offset + count) {
        butterfly(values, offset, stride);
    }
    let last_offset = offset + size * stride;
    for (i, offset) in (offset..last_offset)
        .step_by(2 * stride)
        .enumerate()
        .skip(1)
    {
        for j in offset..(offset + count) {
            butterfly_twiddle(values, twiddles[i], j, stride);
        }
    }
}
// HELPER FUNCTIONS
// ================================================================================================
#[inline(always)]
fn butterfly<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField>(
    values: &mut [F],
    offset: usize,
    stride: usize,
) {
    let i = offset;
    let j = offset + stride;
    let temp = values[i];
    values[i] = temp + values[j];
    values[j] = temp - values[j];
}

#[inline(always)]
fn butterfly_twiddle<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField>(
    values: &mut [F],
    twiddle: P,
    offset: usize,
    stride: usize,
) {
    let i = offset;
    let j = offset + stride;
    let temp = values[i];
    values[j] = values[j].mul_base(twiddle);
    values[i] = temp + values[j];
    values[j] = temp - values[j];
}
pub fn log2(n: usize) -> u32 {
    assert!(n.is_power_of_two(), "n must be a power of two");
    n.trailing_zeros()
}
// POLYNOMIAL INTERPOLATION
// ================================================================================================
/// Interpolates `evaluations` over a domain of length `evaluations.len()` in the field specified
/// `B` into a polynomial in coefficient form using the FFT algorithm.
pub fn interpolate_poly<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField>(
    evaluations: &mut [F],
    inv_twiddles: &[P],
) {
    fft_in_place(evaluations, inv_twiddles, 1, 1, 0);
    let inv_length = P::invert(P::from(evaluations.len() as u64)).unwrap();
    for e in evaluations.iter_mut() {
        *e = e.mul_base(inv_length);
    }
    permute(evaluations);
}

/// Interpolates `evaluations` over a domain of length `evaluations.len()` and shifted by
/// `domain_offset` in the field specified by `B` into a polynomial in coefficient form using
/// the FFT algorithm.
pub fn interpolate_poly_with_offset<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField>(
    evaluations: &mut [F],
    inv_twiddles: &[P],
    domain_offset: P,
) {
    fft_in_place(evaluations, inv_twiddles, 1, 1, 0);
    permute(evaluations);
    let domain_offset = domain_offset.invert().unwrap();
    let mut offset = P::invert(P::from(evaluations.len() as u64)).unwrap();
    for coeff in evaluations.iter_mut() {
        *coeff = coeff.mul_base(offset);
        offset *= domain_offset;
    }
}

/// Returns the degree of the provided polynomial.
///
/// If the size of the provided slice is much larger than the degree of the polynomial (i.e.,
/// a large number of leading coefficients is ZERO), this operation can be quite inefficient.
pub fn degree_of<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField>(poly: &[F]) -> usize {
    for i in (0..poly.len()).rev() {
        if poly[i] != F::ZERO {
            return i;
        }
    }
    0
}

//...........
// CODE  for evaluating polynomial at points
//.............
pub fn eval<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField >(p: &[F], x: P) -> F

{
    // Horner evaluation
    p.iter()
        .rev()
        .fold(F::ZERO, |acc, &coeff| F::mul_base(acc, x) + coeff)
}

pub fn eval_many<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField >(p: &[F], xs: &[P]) -> Vec<F>
{
    xs.par_iter().map(|x| eval(p, *x)).collect()
}