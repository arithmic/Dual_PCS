# Fast Fourier Transforms

This repository  provides an implementation of the Fast Fourier Transform (FFT) algorithm in Rust. The implementation is based on the Cooley-Tukey radix-2 decimation-in-time (DIT) algorithm, which is commonly used for efficient FFT calculations. Fast Fourier Transform describes an algorithm to compute polynomial multiplication in O(n log n) as opposed to O(n^2) when computed naiively. The Fast Fourier Transform is an algorithm to compute the Discrete Fourier transform of a degree n polynomial in O(n log n) as opposed to O (n^2), as after computing the Discrete Fourier Transform, polynomial multiplication can be done in O(n) time, thus bringing the complexity to O(n log n). The Discrete Fourier Transform uses properties specific to primitive roots of unity to characterise polynomials in terms of their evaluations, such that the correspondence is efficient to compute. 

**FFT crate** inculdes the following functionalities :
* `evaluate_poly()`
* `evaluate_poly_with_offset()`
* `interpolate_poly()`
* `interpolate_poly_with_offset()`
* `get_twiddles()`
* `get_inv_twiddles()`

These function are implemented in both serial(std) mode and the concurrent mode(operations will be executed in multiple threads and number of threads can be configured via RAYON_NUM_THREADS environment variable).

