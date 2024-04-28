use bn254::{fp::Fp, fp12::Fp12, fp2::Fp2, fp6::Fp6, scalar::SCALAR_MODULUS};
use crypto_bigint::U256;
use std::ops::{Add, Mul, Not};
use traits::traits::Field;
use BN_curve::{
    bncurve::{AffinePoint, BNCurve, ProjectivePoint},
    fp12bns::FP12AffinePoint,
    fp2bn::G2AffinePoint,
};

use crate::line::Line;

// The BN parameter x for BN12-254 is 4965661367192848881
pub const BN_X: u64 = 4965661367192848881;
pub const BN_X_IS_NEGATIVE: bool = false;

// Function to map a point of E[F_p^2]  to E[F_p^{12}]
pub fn untwist(point: G2AffinePoint) -> FP12AffinePoint {
    let res = point;
    let x1 = res.x;
    let y1 = res.y;
    let root = Fp6::<Fp>::new(Fp2::ZERO, Fp2::ONE, Fp2::ZERO);
    let mut x = Fp12::<Fp>::new(
        Fp6 {
            c0: x1,
            c1: Fp2::ZERO,
            c2: Fp2::ZERO,
        },
        Fp6::ZERO,
    );
    x *= Fp12::new(root, Fp6::ZERO);
    let mut y = Fp12::new(Fp6::new(y1, Fp2::ZERO, Fp2::ZERO), Fp6::ZERO);
    y *= Fp12::new(Fp6::ZERO, root);
    FP12AffinePoint { x, y, infinity: 0 }
}

// Returns the multipliction of F_{p^12} element and the element obtained when line is evaluated at E[F_(p^12)] point
pub fn line_eval_fp12(f: Fp12<Fp>, l: Line<BNCurve>, n1: &FP12AffinePoint) -> Fp12<Fp> {
    let x = n1.x;
    let y = n1.y;
    let z = Fp12::<Fp> {
        c0: Fp6::<Fp> {
            c0: Fp2 {
                c0: l.c2,
                c1: Fp::ZERO,
            },
            c1: Fp2::ZERO,
            c2: Fp2::ZERO,
        },
        c1: Fp6::ZERO,
    };
    let a = x.mul_scalar(l.c0);
    let b = y.mul_scalar(l.c1);
    let add = (a.add(b)).add(z);
    let res = f.mul(add);
    res
}

// Function to compute the miller loop for pairing
fn miller_loop(m: ProjectivePoint, n: FP12AffinePoint, r: U256) -> Fp12<Fp> {
    let mut t = m;
    let mut f = Fp12::ONE;

    if (!m.to_affine().is_identity()).into() && n.is_infinity().not() {
        let k = r.bits() - 1;
        for i in (0..k).rev() {
            f = f.square();
            let mut l = Line::line_fn(t, t);

            f = line_eval_fp12(f, l, &n);
            t = t.double();

            if r.bit(i) == 1 {
                l = Line::line_fn(t, m);
                f = line_eval_fp12(f, l, &n);
                t = t.add(&m);
            }
        }
    }
    f
}

// Function to compute the final exponentiation i.e  point^{p^12-1/r}
pub fn final_exponentiation(r: &Fp12<Fp>) -> Fp12<Fp> {
    let mut f1 = *r;
    f1 = f1.conjugate(); //(f)^{p^6}=conjugate(f) when f is in cyclotomic subgroup of Fp12
                         // f1.frobenius_map(6);
    let mut f2 = r.invert().unwrap();
    let mut r = f1;
    r *= f2; //r=f^(p^6-1)
    f2 = r;
    r.frobenius_map(2); //r={f^(p^6-1)}^p^2
    r *= f2; //r=f^(p^6-1)(p^2+1)  : Easy Part
    fn exp_by_x(value: &mut Fp12<Fp>, x: u64) {
        *value = value.power_by(&[x]);
    }
    //Hard Part
    let f = r;
    let x = BN_X;
    // let mut a1 = f.power_by(U256::from_u128(6).to_words());
    let mut a1 = f * f * f * f * f * f;
    exp_by_x(&mut a1, x);
    // a1 = (a1 * f.power_by([5 as u64, 0, 0, 0])).invert().unwrap();
    a1 = (a1 * f * f * f * f * f).invert().unwrap();
    let a = a1;
    a1.frobenius_map(1);
    let mut b = a1;
    b = a * b;

    let mut f_p = f;
    f_p.frobenius_map(1);

    let mut f_p2 = f;
    f_p2.frobenius_map(2);

    let mut f_p3 = f;
    f_p3.frobenius_map(3);
    let mut y1 = b * f_p.square() * f_p2;
    let value = y1;
    // y1 = y1.power_by(U256::from_u128(6).to_words());
    y1 = y1 * y1 * y1 * y1 * y1 * y1;
    exp_by_x(&mut y1, x);
    exp_by_x(&mut y1, x);
    y1 = y1 * value;
    // let y2 = b * (f_p * *f).power_by([9 as u64, 0, 0, 0]);
    let mut value = f_p * f;
    value = value * value * value * value * value * value * value * value * value;
    let y2 = b * (value);
    // let y3 = a * f.power_by(&[4 as u64, 0, 0, 0]);
    let y3 = a * (f * f * f * f);
    f_p3 * y1 * y2 * y3
}

// Computes tate pairing
pub fn bn_compute_pairing(m: AffinePoint, n: G2AffinePoint) -> Fp12<Fp> {
    // Convert Fp2Affine Point to Fp12Affine Point
    let n_fp12 = untwist(n);

    // Perform miller loop
    let result = miller_loop(m.to_projective(), n_fp12, SCALAR_MODULUS);

    // Apply final exponentiation to the result
    return final_exponentiation(&result);
}
