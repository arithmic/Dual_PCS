
# Pairing

This repository contains pairing which is used in signature scheme of pairing-friendly elliptic curves. An elliptic curve pairing is a nondegenerate bilinear map e: $\ E[ r ] * E[ r] -> u_r$ , where $\ E[ r ]$ is r-torison subgroup of elliptic curve group and $u_r$ is a multiplicative group of r-th roots of unity in the field $\ F_p^k$(where p is the characteristic of the field over which the elliptic curve is constructed and k is the embedding degree).
The  Pairing crate provides the following traits:
* **Pairing  trait**: it describes the interface for pairing-friendly elliptic curves by defining the functions for the tate pairing and the optimal ate pairing.

The Pairing trait is implementated for:
1. BLS12-381 curve
2. BN254 curve


#### Example
Example defined below uses the BN254 curve, we can also use BLS12-381 instaed of BN254 curve.
``` 
        let g1_generator = AffinePoint::GENERATOR;
        let g2_generator = G2AffinePoint::GENERATOR;
        let tate_pairing = BNCurve::tate_pairing(g1_generator, g2_generator);
        let optimal_ate_pairing = BNCurve::optimal_ate_pairing(g2_generator, g1_generator);

```

#### How to use pairing crate in your project:
Bring the pairing crate into your project by using : <br>
pairing = { git = "ssh://git@github.com/arithmic/digital_signatures.git", branch ="main" }.

#### Dependencies :
1. crypto_bigint [https://github.com/arithmic/crypto_bigint].
2. field traits repository [https://github.com/arithmic/field_trait.git].
3. ecc_curve repository [https://github.com/arithmic/ecc_curve.git].
