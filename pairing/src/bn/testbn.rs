#[cfg(test)]
mod tests {
    use crate::{
        bn::{
            bn_pairing::{final_exponentiation, untwist},
            optimal::optimal_line_eval,
        },
        pairing_traits::Pairing,
    };
    use bn254::{
        fp::Fp,
        fp12::Fp12,
        fp2::Fp2,
        fp6::Fp6,
        scalar::{Scalar, SCALAR_MODULUS},
    };
    use std::ops::Mul;
    use traits::traits::Field;
    use BN_curve::{
        bncurve::{AffinePoint, BNCurve, ProjectivePoint},
        fp2bn::{G2AffinePoint, G2ProjectivePoint}, gt::Gt,
    };

    // Test case to check the tate_pairing function with property e(Summation a[i],n) = Product {e(a[i],n)}
    #[test]
    pub fn prop1_first_arg() {
        for _ in 0..10 {
            let mut list  = [ProjectivePoint::IDENTITY; 50];
            let mut sum = ProjectivePoint::IDENTITY;
            for i in 0..50{
                list[i] = ProjectivePoint::GENERATOR.mul(Scalar::random());
                sum = sum + list[i]
            }
            // let a = ProjectivePoint::GENERATOR.mul(Scalar::random());
            let n = G2ProjectivePoint::GENERATOR
                .mul(Scalar::random())
                .to_affine();
            let pair1: BN_curve::gt::Gt = BNCurve::tate_pairing(sum.to_affine(), n);
            let mut pair;
            let mut pair2 = Gt::ONE;
            for i in 0..50{
                pair = BNCurve::tate_pairing(list[i].to_affine(), n);
                pair2 = pair + &pair2;
            }
            assert_eq!(pair1, pair2);
        }
    }
    // Test case to check the tate_pairing function with property e(a, Summation b[i]) = Product {e(a, b[i])}
    #[test]
    pub fn prop2_second_arg() {
        for _ in 0..10 {
            let mut list  = [G2ProjectivePoint::IDENTITY; 50];
            let mut sum = G2ProjectivePoint::IDENTITY;
            for i in 0..50{
                list[i] = G2ProjectivePoint::GENERATOR.mul(Scalar::random());
                sum = sum + list[i]
            }
            // let a = ProjectivePoint::GENERATOR.mul(Scalar::random());
            let a = ProjectivePoint::GENERATOR
                .mul(Scalar::random())
                .to_affine();
            let pair1: BN_curve::gt::Gt = BNCurve::tate_pairing(a, sum.to_affine());
            let mut pair;
            let mut pair2 = Gt::ONE;
            for i in 0..50{
                pair = BNCurve::tate_pairing(a, list[i].to_affine());
                pair2 = pair + &pair2;
            }
            assert_eq!(pair1, pair2);
        }
    }

    #[test]
    fn test_final_exponentiation() {
        let a = Fp12::<Fp>::random();
        let pow = a.power_by(&SCALAR_MODULUS.to_words());
        let result = final_exponentiation(&pow);
        assert_eq!(result, Fp12::ONE);
    }

    #[test]
    fn test_twist() {
        let a = G2AffinePoint::GENERATOR;
        let b = untwist(a);
        let lhs = b.y.square();
        let mut rhs = b.x.square().mul(b.x);
        let three = Fp12::<Fp>::new(
            Fp6::<Fp>::new(
                Fp2::<Fp>::new(Fp::new((3 as u128).into()), Fp::ZERO),
                Fp2::ZERO,
                Fp2::ZERO,
            ),
            Fp6::ZERO,
        );
        rhs += three;
        assert_eq!(lhs, rhs);
    }

    // Test case to check the optimal_ate_pairing function with property e(m1,2*n1)  =e(2*m1,n1)
    #[test]
    pub fn optimal_test_pairing() {
        let m1 = G2ProjectivePoint::GENERATOR;
        let m2 = m1.mul(Scalar([2, 0, 0, 0].into()));
        let n1 = AffinePoint::GENERATOR;
        let n2 = n1.to_projective().double().to_affine();
        let pair1 = BNCurve::optimal_ate_pairing(m1.to_affine(), n2);
        let pair2 = BNCurve::optimal_ate_pairing(m2.to_affine(), n1);
        assert_eq!(pair1, pair2);
    }

    // Test case to check the optimal_ate_pairing function with property e(a,n) * e(b,n) =e(a + b,n)
    #[test]
    pub fn optimal_ate_pairing_check() {
        for _ in 0..100 {
            let a = G2ProjectivePoint::GENERATOR.mul(Scalar::random());
            let b = G2ProjectivePoint::GENERATOR.mul(Scalar::random());
            let n = AffinePoint::GENERATOR
                .to_projective()
                .mul(Scalar::random())
                .to_affine();
            let pair1 = BNCurve::optimal_ate_pairing(a.to_affine(), n);
            let pair2 = BNCurve::optimal_ate_pairing(b.to_affine(), n);
            let sum = a + b;
            let pair3 = BNCurve::optimal_ate_pairing(sum.to_affine(), n);
            assert_eq!(pair3, pair1 + &pair2);
        }
    }

    #[test]
    pub fn line_test() {
        for _ in 0..100 {
            let f = Fp12::<Fp>::random();
            let p = G2ProjectivePoint::GENERATOR;
            let q = p;
            let n = ProjectivePoint::GENERATOR;
            let l1 = optimal_line_eval(f, p.to_affine(), q.to_affine(), n.to_affine());
            let point2 = -(p + q);
            let l2 = optimal_line_eval(f, p.to_affine(), point2.to_affine(), n.to_affine());
            assert_eq!(l1, l2);
        }
    }
}

// // Test case to check the tate_pairing function with property e(a,n)*e(2a,n) = e(3a,n)
// #[test]
// pub fn pairing_test() {
//     let a = AffinePoint::GENERATOR;
//     let n = G2AffinePoint::GENERATOR;
//     let pair1 = BNCurve::tate_pairing(a, n);
//     let pair2 = BNCurve::tate_pairing(a.to_projective().double().into(), n);
//     let first =
//         BNCurve::tate_pairing(a.to_projective().add(&a.to_projective().double()).into(), n);
//     let second = pair1.add(&pair2);
//     assert_eq!(first, second);
// }

// // Test case to check the tate_pairing function with property e(a,n)*e(b,n) = e(a+b,n)
// #[test]
// pub fn pairing_check() {
//     for _ in 0..100 {
//         let a = ProjectivePoint::GENERATOR.mul(Scalar::random());
//         let b = ProjectivePoint::GENERATOR.mul(Scalar::random());
//         let n = G2AffinePoint::GENERATOR
//             .to_projective()
//             .mul(Scalar::random())
//             .to_affine();
//         let pair1 = BNCurve::tate_pairing(a.to_affine(), n);
//         let pair2 = BNCurve::tate_pairing(b.to_affine(), n);
//         let sum = a + b;
//         let pair3 = BNCurve::tate_pairing(sum.to_affine(), n);
//         assert_eq!(pair3, pair1 + &pair2);
//     }
// }