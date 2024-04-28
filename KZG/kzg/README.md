# KZG Commitment Scheme

KZG stands for Kate, Zaverucha and Goldberg. 
This scheme is used for the commitmnet of the polynomial.
The KZG commitment scheme proceeds in following steps:
1. Trusted setup 
2. Commitment of the polynomial
3. Evaluation proof of the polynomial
4. Verification

The KZG commitment scheme is implemented for 
  * BN254 curve
  * BLS12-381 curve


### Example
Here is an example of KZG commitment scheme using the BLS12-381 curve.
```       
            // trusted setup
            let (g1_powers, g2_powers) = BlsCurve::setup(1000);
            //Generating the random polynomial 'p'
            let p = Polynomial::<Scalar<BlsCurve>>::random(1000);
            let z = Scalar::<BlsCurve>::random();
            // Evaluation of the polynomial 'p' at z
            let y = p.evaluate(&z);
            // Commitment of the polynomial
            let c = BlsCurve::commit(&p, &g1_powers);
            // Evaluation proof of the polynomial
            let pi = BlsCurve::open(p, &g2_powers, z);
            // Verification
            let v = BlsCurve::verify(c, z, y, pi, &g1_powers, &g2_powers[0]);
            assert_eq!(v, true); 
```

### Dependencies :
1. crypto_bigint [https://github.com/arithmic/crypto_bigint].
2. field traits repository [https://github.com/arithmic/Field_Open.git].
3. ecc_curve repository [https://github.com/arithmic/ECC_Open.git].
