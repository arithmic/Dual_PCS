use crate::line::Line;
use bls381::{fp::Fp, fp12::Fp12, fp2::Fp2, fp6::Fp6, scalar::SCALAR_MODULUS};
use bls_curve::{
    bls::{BlsCurve, ProjectivePoint, AffinePoint},
    fp12bls::FP12AffinePoint,
    fp2bls::G2AffinePoint,
};
use crypto_bigint::U256;
use std::ops::{Add, Mul, Not};
use traits::traits::Field;

// The BLS parameter x for BLS12-381 is -0xd201000000010000
pub const BLS_X: u64 = 0xd201000000010000;
pub const BLS_X_IS_NEGATIVE: bool = true;

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
    x *= Fp12::new(root, Fp6::ZERO).invert().unwrap();
    let mut y = Fp12::new(Fp6::new(y1, Fp2::ZERO, Fp2::ZERO), Fp6::ZERO);
    y *= Fp12::new(Fp6::ZERO, root).invert().unwrap();
    FP12AffinePoint { x, y, infinity: 0 }
}

// Returns the multipliction of F_{p^12} element and the element obtained when line is evaluated at E[F_(p^12)] point
pub fn line_eval_fp12(f: Fp12<Fp>, l: Line<BlsCurve>, n1: &FP12AffinePoint) -> Fp12<Fp> {
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

// Function to compute the miller loop
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

//Algorithm2 in Technical Specifications Document
// Function to compute the final exponentiation i.e  point^{p^12-1/r}
pub fn final_exponentiation(r: &Fp12<Fp>) -> Fp12<Fp> {
    let mut f1 = *r;
    f1 = f1.conjugate(); //(f)^{p^6}=conjugate(f) when f is in cyclotomic subgroup of Fp12
    let mut f2 = r.invert().unwrap();
    let mut r = f1;
    r *= f2; //r=f^(p^6-1)
    f2 = r;
    r.frobenius_map(2); //r={f^(p^6-1)}^p^2
    r *= f2; //r=f^(p^6-1)(p^2+1)  : Easy Part
    fn exp_by_x(f: &mut Fp12<Fp>, x: u64) {
        *f = f.power_by(&[x]);
        if BLS_X_IS_NEGATIVE {
            *f = f.conjugate();
        }
    }
    //Hard Part
    let mut x = BLS_X;
    let mut y0 = r;
    y0 = y0.square(); //line 3
    let mut y1 = y0;
    exp_by_x(&mut y1, x); //line 5, y1=r^(2x)
    x >>= 1;
    let mut y2 = y1;
    exp_by_x(&mut y2, x); //y2=r^(x^2)
    x <<= 1;
    let mut y3 = r;
    y3 = y3.conjugate();
    y1 *= y3; //y1=r^(2x).r^(-1)
    y1 = y1.conjugate(); //y1=r^(-2x).r
    y1 *= y2; //y1=r^(x^2-2x).r
    y2 = y1;
    exp_by_x(&mut y2, x); //y2=r^(x^3-2x^2).r^x
    y3 = y2;
    exp_by_x(&mut y3, x); //y3=r^(x^4-2x^3).r^(x^2)
    y1 = y1.conjugate();
    y3 *= y1; //y3=r^(x^4-2x^3+2x).r^(-1)
    y1 = y1.conjugate();
    y1.frobenius_map(3); //y1=[r^(x^2-2x).r]^(p^3)       a
    y2.frobenius_map(2); //y2=[r^(x^3-2x^2).r^x]^(p^2)   b
    y1 *= y2; //y1=ab
    y2 = y3;
    exp_by_x(&mut y2, x); //y2=r^(x^5-2x^4+2x^2).r^(-x)
    y2 *= y0;
    y2 *= r; //y2=r^(x^5-2x^4+2x^2).(r^(x-2))^-1.r  c
    y1 *= y2; //y1=abc
    y2 = y3;
    y2.frobenius_map(1); //y2=[r^(x^4-2x^3+2x).r^(-1)]^(p)     d
    y1 *= y2; //y1=abcd
    y1
}

// Computes tate pairing for the bls curve
pub fn bls_compute_pairing(m: AffinePoint, n: G2AffinePoint) -> Fp12<Fp> {
     // Convert Fp2Affine Point to Fp12Affine Point
    let n_fp12 = untwist(n);
        // Perform miller loop 
    let res_fp12 = miller_loop(m.to_projective(), n_fp12, SCALAR_MODULUS);
    return final_exponentiation(&res_fp12);
}
