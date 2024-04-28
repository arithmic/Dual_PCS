# Dual Polynomial Commitment Scheme and Applications

## Installation

Clone this repo:

`git@github.com:arithmic/Dual_PCS.git`

Clone dependencies

Field
`git@github.com:arithmic/Field_open.git`

Elliptic Curve
`git@github.com:arithmic/ECC_Open.git`

Pairing
`git@github.com:arithmic/Pairing_Open.git`

Hashes
`git@github.com:arithmic/non-algebraic-hashes.git`

FFT
`git@github.com:arithmic/FFT_Open.git`

Random Coin
`git@github.com:arithmic/Random_Coin_Open.git`

Table
`git@github.com:arithmic/Table_Open.git`

To run tests:

```text
RUSTFLAGS="-C target_cpu=native" cargo test
```

## Performance

### benchmarks

`main` includes three benches: `benches/linking_benchmark.rs` and `benches/univariate_benchmark.rs` and `benches/multilinear_benchmark.rs`. to bench the `Linking Proof using Argument System`

`setup` includes 1 bench: `benches/linking_setup_benchmark.rs` to bench the setup of `Linking Proof using Argument System`

`Multilinear` includes two benches: `Multilinear/kzg_fourier_multilinear/benches/kzg_fourier_benchmark.rs` and `Multilinear/multilinear_kzg/benches/multilinear_kzg_benchmark.rs` to bench the `Multilinear KZG Fourier` and `Multilinear KZG`

`Spartan` includes two benches: `Spartan/Spartan_with_air/benches/spartan_benchmark.rs` and `Spartan/Spartan_with_gkr/benches/spartan_benchmark.rs` to bench the `Spartan with AIR` and `Spartan with GKR`

`KZG` includes two benches: `KZG/kzg/benches/kzg_benchmark.rs` and `KZG/kzg_fft/benches/kzg_fft_benchmark.rs` to bench the `KZG` and `KZG-FFT`

`Grand_product` includes two benches: `Grand_product/grand_product_with_air/benches/grand_product_benchmark.rs` and `Grand_product/grand_product_with_gkr/benches/grand_product_benchmark.rs` to bench the `Grand Product with AIR` and `Grand Product with GKR`

To run end-to-end benchmarks:

```text
RUSTFLAGS="-C target_cpu=native" cargo bench
```
