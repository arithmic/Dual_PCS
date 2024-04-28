use std::ops::{Mul, Neg, Not};
use bls381::{fp12::Fp12, fp::Fp};
use bls_curve::{fp2bls::{G2AffinePoint, G2ProjectivePoint}, bls::AffinePoint};
use crypto_bigint::U256;
use traits::traits::Field;
use super::bls_pairing::{untwist, BLS_X_IS_NEGATIVE, final_exponentiation, BLS_X};

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
        slope *= p_fp12.y.double().invert().unwrap();
        v = p_fp12.y - (slope * p_fp12.x);
    } else {
        if p_fp12.x == q_fp12.x && p_fp12.y == q_fp12.y.neg() {
            return f.mul(Fp12::one().mul_scalar(n.x) - p_fp12.x);
        }
        slope = p_fp12.y - q_fp12.y;
        slope *= (p_fp12.x - q_fp12.x).invert().unwrap();
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
    }
    if BLS_X_IS_NEGATIVE {
        f = f.conjugate();
    }

    f
}


// Function to define optimal ate pairing for bls curve 
pub fn bls_optimal_ate_pairing(m: G2AffinePoint, n: AffinePoint) -> Fp12<Fp> {
    // Miller loop with optimal ate pairing
    let res_fp12 = optimal_miller_loop(m.to_projective(), n, U256::from_u64(BLS_X));
     // Apply final exponentiation to the result
    return final_exponentiation(&res_fp12);

}


