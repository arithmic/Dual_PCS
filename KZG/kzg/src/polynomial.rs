use bls381::scalar::Scalar;
use std::ops::{Add, Neg, Sub};
use traits::traits::Field;

// struct defination of the Polynomial
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Polynomial {
    pub coeffs: Vec<Scalar>,
}

impl Polynomial {
    /// Returns the zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }
    /// Checks if the given polynomial is zero.
    fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || self.coeffs.iter().all(|coeff| coeff.is_zero())
    }
    /// Outputs a univariate polynomial of degree `d` where
    /// each coefficient is sampled uniformly at random.
    pub fn random(d: usize) -> Self {
        let mut random_coeffs = Vec::new();
        for _ in 0..d {
            random_coeffs.push(Scalar::random());
        }
        Self::from_coefficients_vec(random_coeffs)
    }
    // While there are zeros at the end of the coefficient vector, pop them off.
    fn truncate_leading_zeros(&mut self) {
        while self.coeffs.last().map_or(false, |c| c.is_zero()) {
            self.coeffs.pop();
        }
    }
    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_vec(coeffs: Vec<Scalar>) -> Self {
        let mut result = Self { coeffs };
        // While there are zeros at the end of the coefficient vector, pop them off.
        result.truncate_leading_zeros();
        // Check that either the coefficients vec is empty or that the last coeff is
        // non-zero.
        assert!(result.coeffs.last().map_or(true, |coeff| !coeff.is_zero()));
        result
    }
    /// Returns the total degree of the polynomial
    pub fn degree(&self) -> usize {
        if self.is_zero() {
            0
        } else {
            assert!(self.coeffs.last().map_or(false, |coeff| !coeff.is_zero()));
            self.coeffs.len() - 1
        }
    }
    /// Evaluates polynomial at the given `point`
    pub fn evaluate(&self, point: &Scalar) -> Scalar {
        if self.is_zero() {
            return Scalar::ZERO;
        } else if point.is_zero() {
            return self.coeffs[0];
        }
        Self::horner_evaluate(&self.coeffs, &point)
    }
    // Horner's method for polynomial evaluation
    pub fn horner_evaluate(poly_coeffs: &[Scalar], point: &Scalar) -> Scalar {
        poly_coeffs
            .iter()
            .rfold(Scalar::ZERO, move |result, coeff| result * *point + *coeff)
    }
    // compute q(x) = phi(x)/(x-z)
    pub fn divide_by_lin_factor(&self, z: Scalar) -> Self {
        let z_inv = z.invert().unwrap();
        let mut quotient_coeffs: Vec<Scalar> = Vec::new();
        let mut temp = Scalar::ZERO;
        for i in 0..self.coeffs.len() {
            let a = ((self.coeffs[i] - temp) * z_inv).neg();
            quotient_coeffs.push(a);
            temp = a;
        }
        Polynomial {
            coeffs: quotient_coeffs,
        }
    }
}

impl Add for Polynomial {
    type Output = Polynomial;

    fn add(self, other: Polynomial) -> Self {
        &self + &other
    }
}

impl Add<&Polynomial> for &Polynomial {
    type Output = Polynomial;
    // Addition of the two polynomials f(x) and g(x) i.e  f(x)+g(x)
    fn add(self, other: &Polynomial) -> Polynomial {
        let mut result = if self.is_zero() {
            other.clone()
        } else if other.is_zero() {
            self.clone()
        } else if self.degree() >= other.degree() {
            let mut result = self.clone();
            result
                .coeffs
                .iter_mut()
                .zip(&other.coeffs)
                .for_each(|(a, b)| {
                    *a += *b;
                });
            result
        } else {
            let mut result = other.clone();
            result
                .coeffs
                .iter_mut()
                .zip(&self.coeffs)
                .for_each(|(a, b)| {
                    *a += *b;
                });
            result
        };
        result.truncate_leading_zeros();
        result
    }
}

impl Neg for Polynomial {
    type Output = Polynomial;

    #[inline]
    fn neg(mut self) -> Polynomial {
        self.coeffs.iter_mut().for_each(|coeff| {
            *coeff = -*coeff;
        });
        self
    }
}

impl Sub<&Polynomial> for &Polynomial {
    type Output = Polynomial;

    // Subtraction of the two polynomials f(x) and g(x) i.e f(x)-g(x)
    #[inline]
    fn sub(self, other: &Polynomial) -> Polynomial {
        let mut result = if self.is_zero() {
            let mut result = other.clone();
            result.coeffs.iter_mut().for_each(|c| *c = -(*c));
            result
        } else if other.is_zero() {
            self.clone()
        } else if self.degree() >= other.degree() {
            let mut result = self.clone();
            result
                .coeffs
                .iter_mut()
                .zip(&other.coeffs)
                .for_each(|(a, b)| *a -= *b);
            result
        } else {
            let mut result = self.clone();
            result.coeffs.resize(other.coeffs.len(), Scalar::ZERO);
            result
                .coeffs
                .iter_mut()
                .zip(&other.coeffs)
                .for_each(|(a, b)| *a -= *b);
            result
        };
        result.truncate_leading_zeros();
        result
    }
}
