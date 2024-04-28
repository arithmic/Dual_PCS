use crypto_bigint::{U256, U64};
use rayon::{iter::{IntoParallelRefMutIterator, IndexedParallelIterator, IntoParallelIterator, ParallelIterator}, slice::{ParallelSliceMut, ParallelSlice}};
use stark252::field::Fp;
use traits::traits::{Field, ExtensionOf, PrimeField};
use utils::uninit_vector;
use crate::{bit_reverse::{par_permute, permute, permute_index}, serial_fft::{fft_in_place, log2}};

//assuming a is an array whose length is a power of 2 and at most 2^20.

// pub fn parallel_fft<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField >(poly: &mut [F], omega: &Vec<P>, threads:usize){
//     let N = poly.len();
//     let bits = N.trailing_zeros();
//     //We iterate over the bits of binary representation of the array length.
 
//     for k in 0..(bits>>1){
//         //iterate over number of problems which double every round of subidivison
//         let problems = 2<<k;
//         let problemsize = N/problems;
        
//         println!("Size of problem: {:?} \n \n", problemsize);
        
//         //We if we are in a problem of size n, we only iterate over n/2 of them
//         //for butterfly computations

//         let halfsize = problemsize>>1;
//         let breakpoints = [problemsize>>2, problemsize>>1, 3*(problemsize>>2), problemsize];

//         poly.par_chunks_mut(problemsize).for_each(|chunk| {
//             unsafe{
//                 let first_quarter = (&mut *(&mut chunk[..breakpoints[0]] as *mut [F])).par_iter_mut();
//                 let second_quarter = (&mut *(&mut chunk[breakpoints[0]..breakpoints[1]] as *mut [F])).par_iter_mut();
//                 let third_quarter = (&mut *(&mut chunk[breakpoints[1]..breakpoints[2]] as *mut [F])).par_iter_mut();
//                 let fourth_quarter =  (&mut *(&mut chunk[breakpoints[2]..] as *mut [F])).par_iter_mut();
                
//                 let double_1 = first_quarter.zip(third_quarter);
//                 let double_2 = second_quarter.zip(fourth_quarter);
//                 double_1.zip(double_2).enumerate().for_each(|(i,((x_1,x_3),(x_2,x_4)))|{
//                   let k_0 = (*x_1 - *x_3).mul_base(omega[i*problems]);
//                   let k_1 = (*x_2 - *x_4).mul_base(omega[i*problems]);
//                   let k_2 = *x_1 + *x_3;
//                   let k_3 = *x_2 + *x_4;

//                   let t_0 = k_2 + k_3;
//                   let t_1 = k_2 - k_3).mul_base(omega[2*i*problems]);
//                   let t_2 = (k_0 + k_1);
//                   let t_3 = (k_0 - k_1).mul_base(omega[2*i*problems]);
                  
//                   *x_1 = t_0;
//                   *x_2 = t_1;
//                   *x_3 = t_2;
//                   *x_4 = t_3;
//                 })
//               }
//         })
//         //We sort the results according to their index.
//     }
//     //The loop returns the array in bit-reversed order, we fix the indexes.
//    par_permute(poly)
    
// }

pub fn par_eval_poly<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField >(poly: &mut [F], twiddles: &Vec<P>){
    //Empty call to check if root of unity is implemented in the field.
    assert!(poly.len().is_power_of_two());
    assert!(twiddles.len()==(poly.len()/2));
    let threads = rayon::current_num_threads();
    fork_join_fft(poly, twiddles, 1, 1, threads);
    par_permute(poly)
}

pub fn fork_join_fft<F: Field<BaseField = P>+ExtensionOf<P>, P: PrimeField >(poly: &mut [F], omega: &Vec<P>, step:usize, branches:usize,threads:usize){
    let N = poly.len();
    //We if we are in a problem of size n, we only iterate over n/2 of them
     //for butterfly computations
    if N>1 {
        //We create threads on subdivisions until our branches are equal to our threads.
        let halfsize: usize = (N>>1);
        let poly_ref = unsafe { &mut *(poly as *mut [F])};
        if branches < threads{
            let poly_2 =  poly_ref.par_iter_mut().skip(halfsize);

            poly.par_iter_mut().take(halfsize).zip(poly_2).enumerate().into_par_iter().for_each(|(i, x)| {
                let t1 = *x.0 + *x.1;
                let t0 = (*x.0 - *x.1).mul_base(omega[i*step]);
                *x.0=t1;
                *x.1=t0;
        });
        
        }

        //If the branches are equal to the threads in the current recursion level, the algorithm continues sequentially
        //in each thread without spawning new threads.
        else {
        let poly_2 =  poly_ref.iter_mut().skip(halfsize);
        
        poly.iter_mut().take(halfsize).zip(poly_2).enumerate().into_iter().for_each(|(i, x)| {
            let t1 = *x.0 + *x.1;
            let t0 = (*x.0 - *x.1).mul_base(omega[i*step]);
            *x.0=t1;
            *x.1=t0;
    })
    }
    rayon::join(||fork_join_fft(&mut poly[..halfsize], omega, step<<1, branches<<1,threads),
    || fork_join_fft(&mut poly_ref[halfsize..], &omega, step<<1, branches<<1,threads));
}
}



pub fn len_4_interpolate_offset<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField >(evaluations:&mut Vec<F>, twiddles:&Vec<P>, domain_offset: F){
    let t_0 = evaluations[0] + evaluations[2];
    let t_1 = evaluations[1] + evaluations[3];
    let t_2 = evaluations[0] - evaluations[2];
    let t_3 = (evaluations[1] - evaluations[3]).mul_base(twiddles[1]);

    evaluations[0] = t_0 + t_1;
    evaluations[1] = t_2 + t_3;
    evaluations[2] = t_0 - t_1;
    evaluations[3] =  t_2 - t_3;

    let domain_offset = domain_offset.invert().unwrap();
    let mut offset = F::from(evaluations.len() as u64)
        .invert()
        .unwrap();
    for coeff in evaluations.iter_mut() {
        *coeff *= offset;
        offset *= domain_offset;
    }
}

pub fn interpolate_poly_with_offset<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField >(
    evaluations: &mut [F],
    inv_twiddles: &[P],
    domain_offset: F,
) {
    fft_in_place(evaluations, inv_twiddles, 1, 1, 0);
    permute(evaluations);
    let domain_offset = domain_offset.invert().unwrap();
    let mut offset = F::from(evaluations.len() as u64)
        .invert()
        .unwrap();
    for coeff in evaluations.iter_mut() {
        *coeff *= offset;
        offset *= domain_offset;
    }
}

pub fn inplace_butterfly<F:Field>(a:&mut [F], i: usize, j: usize, omega: &F){
    let a_0 = a[i] + a[j];
    let a_1 = (a[i]- a[j])**omega;
    a[i]=a_0;
    a[j]=a_1;
}


//Given an array of coefficients does the butterfly computation on chosen indices
//and returns a double with the result and appropriate position.
pub fn indexed_butterfly(a:&[Fp], i: usize, j: usize, omega: &Fp)->[(Fp,usize);2]{
    [(a[i] + a[j],i), ((a[i] - a[j])*omega,j)]
}


pub fn get_twiddles<F:PrimeField>(N:u32)->Vec<F>{
    let bits = 32 - N.leading_zeros()-1;
    let omega_n = F::get_root_of_unity(bits);
    let omega:Vec<F> = (0..N>>1).into_par_iter().map(|i| omega_n.power_by([i as u64])).collect();
    omega
}

pub fn get_inv_twiddles<F:PrimeField>(N:u32)->Vec<F>{
    let bits = 32 - N.leading_zeros()-1;
    let omega_n = F::invert(F::get_root_of_unity(bits)).unwrap();
    let mut omega:Vec<F> = (0..N>>1).into_par_iter().map(|i| omega_n.power_by([i as u64])).collect();
    omega
}


// POLYNOMIAL INTERPOLATION
// ================================================================================================
/// Interpolates `evaluations` over a domain of length `evaluations.len()` in the field specified
/// `B` into a polynomial in coefficient form using the FFT algorithm.

pub fn par_eval_poly_with_offset<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField >(
    p: &[F],
    twiddles: &Vec<P>,
    domain_offset: P,
    blowup_factor: usize
) -> Vec<F> {
    debug_assert!(blowup_factor.is_power_of_two(),"blow_up factor may not be power of two");
    let domain_size = p.len() * blowup_factor;
    let g = P::get_root_of_unity(log2(domain_size));
    let mut result = unsafe { uninit_vector(domain_size) };
    result
        .as_mut_slice()
        .par_chunks_mut(p.len())
        .enumerate()
        .for_each(|(i, chunk)| {
            let idx = permute_index(blowup_factor, i) as u64;
            let offset = g.power_by(&U256::from(idx).to_words()) * domain_offset;
            let mut factor = P::ONE;
            for (d, c) in chunk.iter_mut().zip(p.iter()) {
                *d = (*c).mul_base(factor);
                factor = factor * offset;
            }
            //fork_join_fft(chunk, twiddles, 1, 1, 8);
            let threads = rayon::current_num_threads();
            fork_join_fft(chunk, twiddles, 1, 1, threads);
        });
    par_permute(&mut result);
    result
}

pub fn par_interpolate_poly<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField >(evaluations: &mut [F], inv_twiddles: Vec<P>){
    
    let threads = rayon::current_num_threads();
    fork_join_fft(evaluations, &inv_twiddles, 1, 1, threads);

    //fork_join_fft(evaluations, &inv_twiddles, 1, 1, 8);
    let inv_length = P::invert(P::from(evaluations.len() as u64)).unwrap();
    evaluations.par_iter_mut().for_each(|e|*e = e.mul_base(inv_length));
    par_permute(evaluations);
}
/// Interpolates `evaluations` over a domain of length `evaluations.len()` and shifted by
/// `domain_offset` in the field specified by `B` into a polynomial in coefficient form using
/// the FFT algorithm.
pub fn par_interpolate_poly_with_offset<F: Field<BaseField = P> + ExtensionOf<P>, P: PrimeField >(
    evaluations: &mut [F],
    inv_twiddles: Vec<P>,
    domain_offset: P
    
) {
    par_eval_poly(evaluations, &inv_twiddles);
    
    let domain_offset = domain_offset.invert().unwrap();
    let mut offset = P::invert(P::from(evaluations.len() as u64)).unwrap();
    for coeff in evaluations.iter_mut() {
        //*coeff *= offset;
        *coeff=coeff.mul_base(offset);
        offset *= domain_offset;
    }
}



#[test]
fn fork_join(){
    use std::time::{Duration, Instant};
    let mut total1 = Duration::ZERO;
    let mut total2 = Duration::ZERO;
    for i in 1..1000{
    let mut poly:Vec<Fp> = (0..4).map(|x: u8| Fp::random()).collect();
    let mut poly2 = poly.clone();
    let omega: Vec<_> = vec![Fp::ONE, Fp::get_root_of_unity(Fp::TWO_ADDICITY - 4)];
    let domain_offset = Fp::random();
    let time1=Instant::now();
    interpolate_poly_with_offset(&mut poly, &omega, domain_offset);
    total1+=time1.elapsed();

    let time2 = Instant::now();
    len_4_interpolate_offset(&mut poly2, &omega, domain_offset);
    total2+=time2.elapsed();
    assert_eq!(poly,poly2)
    }
    println!("time for naiive impl: {:?}, time for spec. len 4 implementation: {:?}", total1, total2)
}

#[test]
fn impl_test(){
    
    let mut poly:Vec<Fp> = (1..5).map(|x: u8| Fp::from(x)).collect();
    let omega = poly.clone();
    let domain_offset = Fp::from(3u8);
    interpolate_poly_with_offset(&mut poly, &omega, domain_offset);
    for p in poly.iter() {
        println!("{:?}", p.0.to_string());
    }

    interpolate_poly_with_offset(&mut poly, &omega, domain_offset);
    for p in poly {
        println!("{:?}", p.0.to_string());
    }
}




////////Winterfell FFT.....................
/// ///..............
/// ...................

pub fn winter_evaluate_poly<F: Field<BaseField = P>+ExtensionOf<P>, P: PrimeField >(p: &mut [F], twiddles: &[P]) {
    split_radix_fft(p, twiddles);
    par_permute(p);
}

pub(super) fn split_radix_fft<F: Field<BaseField = P>+ExtensionOf<P>, P: PrimeField >(
    values: &mut [F],
    twiddles: &[P],
) {
    // generator of the domain should be in the middle of twiddles
    let n = values.len();
    let g = twiddles[twiddles.len() / 2];
    //debug_assert_eq!(g.power_by(U64::from(n as u32).to_words()), B::ONE);

    let inner_len = 1_usize << (n.ilog2() / 2);
    let outer_len = n / inner_len;
    let stretch = outer_len / inner_len;
    debug_assert!(outer_len == inner_len || outer_len == 2 * inner_len);
    debug_assert_eq!(outer_len * inner_len, n);

    // transpose inner x inner x stretch square matrix
    transpose_square_stretch(values, inner_len, stretch);

    // apply inner FFTs
    values
        .par_chunks_mut(outer_len)
        .for_each(|row| fft_in_place(row, &twiddles, stretch, stretch, 0));

    // transpose inner x inner x stretch square matrix
    transpose_square_stretch(values, inner_len, stretch);

    // apply outer FFTs
    values
        .par_chunks_mut(outer_len)
        .enumerate()
        .for_each(|(i, row)| {
            if i > 0 {
                let i = permute_index(inner_len, i);
                let inner_twiddle = g.power_by(U64::from(i as u32).to_words());
                let mut outer_twiddle = inner_twiddle;
                for element in row.iter_mut().skip(1) {
                    *element = (*element).mul_base(outer_twiddle);
                    outer_twiddle = outer_twiddle * inner_twiddle;
                }
            }
            fft_in_place( row,&twiddles, stretch, stretch, 0);
        });
}

// TRANSPOSING
// ================================================================================================

fn transpose_square_stretch<T>(matrix: &mut [T], size: usize, stretch: usize) {
    assert_eq!(matrix.len(), size * size * stretch);
    match stretch {
        1 => transpose_square_1(matrix, size),
        2 => transpose_square_2(matrix, size),
        _ => unimplemented!("only stretch sizes 1 and 2 are supported"),
    }
}

fn transpose_square_1<T>(matrix: &mut [T], size: usize) {
    debug_assert_eq!(matrix.len(), size * size);
    if size % 2 != 0 {
        unimplemented!("odd sizes are not supported");
    }

    // iterate over upper-left triangle, working in 2x2 blocks
    for row in (0..size).step_by(2) {
        let i = row * size + row;
        matrix.swap(i + 1, i + size);
        for col in (row..size).step_by(2).skip(1) {
            let i = row * size + col;
            let j = col * size + row;
            matrix.swap(i, j);
            matrix.swap(i + 1, j + size);
            matrix.swap(i + size, j + 1);
            matrix.swap(i + size + 1, j + size + 1);
        }
    }
}

fn transpose_square_2<T>(matrix: &mut [T], size: usize) {
    debug_assert_eq!(matrix.len(), 2 * size * size);

    // iterate over upper-left triangle, working in 1x2 blocks
    for row in 0..size {
        for col in (row..size).skip(1) {
            let i = (row * size + col) * 2;
            let j = (col * size + row) * 2;
            matrix.swap(i, j);
            matrix.swap(i + 1, j + 1);
        }
    }
}

// HELPER FUNCTIONS
// ================================================================================================

fn clone_and_shift<E: PrimeField>(source: &[E], destination: &mut [E], offset: E) {
    let batch_size = source.len() / rayon::current_num_threads().next_power_of_two();
    source
        .par_chunks(batch_size)
        .zip(destination.par_chunks_mut(batch_size))
        .enumerate()
        .for_each(|(i, (source, destination))| {
            let mut factor = offset.power_by((U64::from((i * batch_size) as u32).to_words()));
            for (s, d) in source.iter().zip(destination.iter_mut()) {
                *d = (*s).mul_base(factor);
                factor = factor * offset;
            }
        });
    }



pub fn cache_fft<F: Field<BaseField = P>+ExtensionOf<P>, P: PrimeField >(
    values: &mut [F],
    twiddles: &[P],
){
    par_permute(values);
    let threshold = (log2(rayon::current_num_threads().next_power_of_two()) -1) as usize;
    DL_fft(values, twiddles, 1, 0,threshold);
}



pub fn DL_fft<F: Field<BaseField = P>+ExtensionOf<P>, P: PrimeField >(
    values: &mut [F],
    twiddles: &[P],
    step: usize,
    depth: usize,
    threshold:usize,
){
    let N = values.len();
    
    if N>1{
        let halfsize = N>>1;

        let values_ref = unsafe {
         &mut *(values as *mut [F])
        };
      
        if depth<threshold{        
            rayon::join(||DL_fft(&mut values_ref[..halfsize], twiddles, step<<1,depth+1,threshold),
            || DL_fft(&mut values[halfsize..], &twiddles, step<<1, depth+1,threshold));
            
        values_ref[..halfsize].par_iter_mut().zip(values[halfsize..].par_iter_mut()).enumerate().for_each(|(i,(x,y))|{
            let temp= *x;
            let q = y.mul_base(twiddles[i * step]);
            *x = temp + q;
            *y = temp - q;
        });
        }
        else {
            rayon::join(||DL_fft(&mut values_ref[..halfsize], twiddles, step<<1,depth+1,threshold),
            || DL_fft(&mut values[halfsize..], &twiddles, step<<1, depth+1,threshold));
            
            for i in 0..halfsize {
                let temp = values[i];
                let q = values[i+halfsize].mul_base(twiddles[i * step]);
                values[i] = temp + q;
                values[i+ halfsize] = temp -q;
            }
        }

    }
}
