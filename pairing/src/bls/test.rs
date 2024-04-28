#[cfg(test)]
mod tests {
    use crate::pairing_traits::Pairing;
    use crate::{bls::bls_pairing::final_exponentiation, line::Line};
    use bls381::{
        fp::Fp,
        fp12::Fp12,
        scalar::{Scalar, SCALAR_MODULUS},
    };
    use bls_curve::gt::Gt;
    use bls_curve::{
        bls::{AffinePoint, BlsCurve, ProjectivePoint},
        fp2bls::{G2AffinePoint, G2ProjectivePoint},
    };
    use traits::traits::Field;

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
            let pair1: bls_curve::gt::Gt = BlsCurve::tate_pairing(sum.to_affine(), n);
            let mut pair;
            let mut pair2 = Gt::ONE;
            for i in 0..50{
                pair = BlsCurve::tate_pairing(list[i].to_affine(), n);
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
            let pair1: bls_curve::gt::Gt = BlsCurve::tate_pairing(a, sum.to_affine());
            let mut pair;
            let mut pair2 = Gt::ONE;
            for i in 0..50{
                pair = BlsCurve::tate_pairing(a, list[i].to_affine());
                pair2 = pair + &pair2;
            }
            assert_eq!(pair1, pair2);
        }
    }
    
    // Test to check  whether the function final exponentiation works properly.
    #[test]
    fn test_final_exponentiation() {
        let a = Fp12::<Fp>::random();
        let pow = a.power_by(&SCALAR_MODULUS.to_words());
        let result = final_exponentiation(&pow);
        assert_eq!(result, Fp12::ONE);
    }

    // Test case to check the optimal_ate_pairing function with property e(m1,2*n1)  =e(2*m1,n1)
    #[test]
    pub fn optimal_test_pairing() {
        let m1 = G2ProjectivePoint::GENERATOR;
        let m2 = m1.mul(Scalar([2, 0, 0, 0].into()));
        let n1 = AffinePoint::GENERATOR;
        let n2 = n1.to_projective().double().to_affine();
        let pair1 = BlsCurve::optimal_ate_pairing(m1.to_affine(), n2);
        let pair2 = BlsCurve::optimal_ate_pairing(m2.to_affine(), n1);
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
            let pair1 = BlsCurve::optimal_ate_pairing(a.to_affine(), n);
            let pair2 = BlsCurve::optimal_ate_pairing(b.to_affine(), n);
            let sum = a + b;
            let pair3 = BlsCurve::optimal_ate_pairing(sum.to_affine(), n);
            assert_eq!(pair3, pair1 + &pair2);
        }
    }

    // Test case to check  whether the line function works properly
    #[test]
    pub fn line_test() {
        let m1 = ProjectivePoint::GENERATOR.double();
        let m2 = ProjectivePoint::GENERATOR;
        let l1 = Line::<BlsCurve>::line_fn(m1, m2);
        let m3 = (m1 + m2).neg();
        let l2 = Line::<BlsCurve>::line_fn(m2, m3);
        assert_eq!(
            l2.c0 * l1.c0.invert().unwrap(),
            l2.c1 * l1.c1.invert().unwrap()
        );
        assert_eq!(
            l2.c2 * l1.c2.invert().unwrap(),
            l2.c1 * l1.c1.invert().unwrap()
        );
    }
}
// // Test case to check the tate_pairing function with property e(a,n)*e(2a,n) = e(3a,n)
    // #[test]
    // pub fn pairing_1() {
    //     let a = ProjectivePoint::GENERATOR;
    //     let n = G2AffinePoint::GENERATOR;
    //     let pair1 = BlsCurve::tate_pairing(a.to_affine(), n);
    //     let pair2 = BlsCurve::tate_pairing(a.double().to_affine(), n);
    //     let first = BlsCurve::tate_pairing(a.add(&a.double()).to_affine(), n);
    //     let second = pair1.add(&pair2);
    //     assert_eq!(first, second);
    // }

    // // Test case to check the tate_pairing function with property e(a,n)*e(b,n) = e(a+b,n)
    // #[test]
    // pub fn pairing_2() {
    //     for _ in 0..100 {
    //         let a = ProjectivePoint::GENERATOR.mul(Scalar::random());
    //         let b = ProjectivePoint::GENERATOR.mul(Scalar::random());
    //         let n = G2AffinePoint::GENERATOR
    //             .to_projective()
    //             .mul(Scalar::random())
    //             .to_affine();
    //         let pair1 = BlsCurve::tate_pairing(a.to_affine(), n);
    //         let pair2 = BlsCurve::tate_pairing(b.to_affine(), n);
    //         let sum = a + b;
    //         let pair3 = BlsCurve::tate_pairing(sum.to_affine(), n);
    //         assert_eq!(pair3, pair1 + &pair2);
    //     }
    // }
