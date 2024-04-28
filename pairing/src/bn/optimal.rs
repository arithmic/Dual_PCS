use super::bn_pairing::{final_exponentiation, untwist};
use bn254::{fp::Fp, fp12::Fp12};
use crypto_bigint::U256;
use std::ops::{Mul, Neg, Not};
use traits::traits::Field;
use BN_curve::{
    bncurve::AffinePoint,
    fp12bns::FP12AffinePoint,
    fp2bn::{G2AffinePoint, G2ProjectivePoint},
};

pub const SIX_X_2: u128 = 29793968203157093288;

// Function for the line evaluation used in the miller loop of optimal ate pairing
pub fn optimal_line_eval(
    f: Fp12<Fp>,
    p: G2AffinePoint,
    q: G2AffinePoint,
    n: AffinePoint,
) -> Fp12<Fp> {
    let mut res = Fp12::one().mul_scalar(n.y);
    let p_fp12 = untwist(p);

    let q_fp12 = untwist(q);
    let mut slope = Fp12::one();
    let mut v = Fp12::one();
    if p == q {
        slope = p_fp12.x.square().double() + p_fp12.x.square();
        slope = slope * p_fp12.y.double().invert().unwrap();
        v = p_fp12.y - (slope * p_fp12.x);
    } else {
        slope = p_fp12.y - q_fp12.y;
        slope = slope * (p_fp12.x - q_fp12.x).invert().unwrap();
        v = (p_fp12.x * q_fp12.y) - (p_fp12.y * q_fp12.x);
        v *= (p_fp12.x - q_fp12.x).invert().unwrap();
    }
    res = res - slope.mul_scalar(n.x) - v;
    f.mul(res)
}

// Miller loop for the optimal ate pairing
fn optimal_miller_loop(m: G2ProjectivePoint, n: AffinePoint, r: U256) -> Fp12<Fp> {
    let mut t = m;
    let mut f = Fp12::ONE;
    let m_a = m.to_affine();
    if m.is_identity().not() && n.is_identity().not().into() {
        let k = r.bits() - 1;
        let mut _count = 0;
        for i in (0..k).rev() {
            f = f.square();
            f = optimal_line_eval(f, t.to_affine(), t.to_affine(), n);
            t = t.double();

            if r.bit(i) == 1 {
                _count += 1;
                f = optimal_line_eval(f, t.to_affine(), m_a, n);
                t = t.add(&m);
            }
        }
        let value = untwist(m.into());
        let q1 = endomorohism(value, 1);
        let t = untwist(t.into());
        f = optimal_line_eval_2(f, t, q1, n);
        let q1 = q1.to_projective();
        let mut t = t.to_projective();
        t = t + q1;
        let q2 = endomorohism(value, 2);
        let q2 = q2.to_projective().neg();
        let t = t.to_affine();
        let q2 = q2.to_affine();
        f = optimal_line_eval_2(f, t, q2, n);
    }
    f
}

pub fn bn_optimal_ate_pairing(m: G2AffinePoint, n: AffinePoint) -> Fp12<Fp> {
    // Miller loop with optimal ate pairing
    let res_fp12 = optimal_miller_loop(m.to_projective(), n, U256::from_u128(SIX_X_2));
    // Apply final exponentiation to the result
    return final_exponentiation(&res_fp12);
}

// Function for the line evaluation used in the miller loop of optimal ate pairing
// Point P and Q are untwist points
pub fn optimal_line_eval_2(
    f: Fp12<Fp>,
    p: FP12AffinePoint,
    q: FP12AffinePoint,
    n: AffinePoint,
) -> Fp12<Fp> {
    let mut res = Fp12::one().mul_scalar(n.y);
    let p_fp12 = p;
    let q_fp12 = q;
    let mut slope = Fp12::one();
    let mut v = Fp12::one();
    if p_fp12 == q_fp12 {
        slope = p_fp12.x.square().double() + p_fp12.x.square();
        slope = slope * p_fp12.y.double().invert().unwrap();
        v = p_fp12.y - (slope * p_fp12.x);
    } else {
        slope = p_fp12.y - q_fp12.y;
        slope = slope * (p_fp12.x - q_fp12.x).invert().unwrap();
        v = (p_fp12.x * q_fp12.y) - (p_fp12.y * q_fp12.x);
        v *= (p_fp12.x - q_fp12.x).invert().unwrap();
    }
    res = res - slope.mul_scalar(n.x) - v;
    f.mul(res)
}

// Function used to compute a -> a^{p}
pub fn endomorohism(value: FP12AffinePoint, power: usize) -> FP12AffinePoint {
    let mut a = value.x;
    Fp12::<Fp>::frobenius_map(&mut a, power);
    let mut b = value.y;
    Fp12::<Fp>::frobenius_map(&mut b, power);
    let q = FP12AffinePoint {
        x: a,
        y: b,
        infinity: 0,
    };
   // assert!(<F12BNCurve as OnCurve<F12BNCurve>>::is_on_curve(q));
    q
}
