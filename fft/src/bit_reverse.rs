use std::cmp::max;
use traits::traits::Field;

pub fn permute_index(size: usize, index: usize) -> usize {
    const USIZE_BITS: usize = 0_usize.count_zeros() as usize;
    debug_assert!(index < size);
    if size == 1 {
        0
    } else {
        debug_assert!(size.is_power_of_two());
        let bits = size.trailing_zeros() as usize;
        index.reverse_bits() >> (USIZE_BITS - bits)
    }
}
/// Permute an array of FFT results.
pub fn permute<T>(v: &mut [T]) {
    let n = v.len();
    for i in 0..n {
        let j = permute_index(n, i);
        if j > i {
            v.swap(i, j);
        }
    }
}
pub fn par_permute<E: Field>(v: &mut [E]) {
    let n = v.len();
    if n<= 8{ permute(v);}
    else{
    let num_subproblem = rayon::current_num_threads().next_power_of_two();
    let subproblem_size = max(n / num_subproblem, 1);
    rayon::scope(|s| {
        for i in 0..num_subproblem {
            // create another mutable reference to the slice of values to use in a new thread; this
            // is OK because we never write the same positions in the slice from different threads
            let values = unsafe { &mut *(&mut v[..] as *mut [E]) };
            s.spawn(move |_| {
                let first = i * subproblem_size;
                let last = first + subproblem_size;
                for j in first..last {
                    let k = permute_index(n, j);
                    if k > j {
                        values.swap(j, k);
                    }
                }
            });
        }
    });
   }
}
// pub fn par_permute<E: Field>(mut v: &mut [E]) {

//     let n = v.len();
//    // let b:Vec<_> =(0..n).into_par_iter().map(|i| permute_index(n,i)).collect();
//     a.par_iter_mut().enumerate().for_each(|(i,x)|);