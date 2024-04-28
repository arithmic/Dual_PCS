use bls381::scalar::Scalar;
use bls_curve::bls::BlsCurve;
use crypto_bigint::Encoding;
use multilinear_kzg::{
    common::{MleCommit, SRS},
    prover::commit,
};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
use std::{fmt::Display, io::Write};
extern crate crypto_bigint;
#[derive(Clone)]
pub struct Polynomial(pub Vec<Scalar>);
impl Display for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for idx in 0..self.get_coefficients().len() - 1 {
            write!(f, "{},", self.get_coefficients()[idx]).expect("failed to writer");
        }
        write!(
            f,
            "{}",
            self.get_coefficients()[self.get_coefficients().len() - 1]
        )
    }
}
impl Polynomial {
    pub fn new(coeffs: Vec<Scalar>) -> Polynomial {
        Polynomial(coeffs)
    }
    pub fn get_coefficients(&self) -> &Vec<Scalar> {
        &self.0
    }
    pub fn get_coefficient_of_degree(&self, degree: usize) -> Scalar {
        self.0[degree]
    }
    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();

        //Write polynomials
        for value in &self.0 {
            result
                .write(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        result
    }
}
#[derive(Clone)]
pub struct LayerPolynomials(Vec<Polynomial>);

impl LayerPolynomials {
    pub fn new(polynomials: Vec<Polynomial>) -> LayerPolynomials {
        LayerPolynomials(polynomials)
    }
    pub fn add_polynomial(&mut self, polynomial: Polynomial) {
        self.0.push(polynomial)
    }
    pub fn polynomial_for_round(&self, round: usize) -> &Polynomial {
        &self.0[round]
    }
    pub fn number_of_rounds(&self) -> usize {
        self.0.len()
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();

        //Write polynomials
        for value in &self.0 {
            result.write(&value.to_bytes()).expect("failed to write");
        }
        result
    }
}

#[derive(Clone)]
pub struct MultPolynomial(pub Vec<Scalar>);

impl Display for MultPolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for idx in 0..self.len() - 1 {
            write!(f, "{},", self.0[idx]).expect("failed to writer");
        }
        write!(f, "{}", self.0[self.len() - 1])
    }
}
impl MultPolynomial {
    pub fn new(fourcoeffs: Vec<Scalar>) -> MultPolynomial {
        MultPolynomial(fourcoeffs)
    }
    pub fn get_coeff(&self, idx: usize) -> Scalar {
        self.0[idx]
    }

    pub fn as_coeffs(&self) -> &Vec<Scalar> {
        &self.0
    }

    pub fn set_coeff(&mut self, idx: usize, value: Scalar) {
        self.0[idx] = value
    }

    pub fn as_coeffs_mut(&mut self) -> &mut Vec<Scalar> {
        &mut self.0
    }

    pub fn par_fold_by_msb(&self, point: Scalar) -> MultPolynomial {
        let halfsize = self.0.len() >> 1;
        let mut res = vec![Scalar::ZERO; halfsize];
        res.par_iter_mut().enumerate().for_each(|(j, res_j)| {
            *res_j = self.0[j] + (self.0[j + halfsize] - self.0[j]) * point;
        });

        MultPolynomial(res)
    }

    pub fn fold_by_msb(&self, point: Scalar) -> MultPolynomial {
        let halfsize = self.0.len() >> 1;
        let mut res = vec![Scalar::ZERO; halfsize];
        for j in 0..halfsize {
            res[j] = self.0[j] + (self.0[j + halfsize] - self.0[j]) * point;
        }

        MultPolynomial(res)
    }

    pub fn fold_by_lsb(&self, point: Scalar) -> MultPolynomial {
        let halfsize = self.0.len() >> 1;
        let mut res = vec![Scalar::ZERO; halfsize];
        for k in 0..halfsize {
            res[k] = (self.0[2 * k] * (Scalar::ONE - point)) + (self.0[2 * k + 1] * point)
        }
        MultPolynomial(res)
    }

    pub fn par_fold_by_lsb(&self, point: Scalar) -> MultPolynomial {
        let halfsize = self.0.len() >> 1;
        let mut res = vec![Scalar::ZERO; halfsize];
        res.par_iter_mut().enumerate().for_each(|(k, res_k)| {
            *res_k = (self.0[2 * k] * (Scalar::ONE - point)) + (self.0[2 * k + 1] * point)
        });
        MultPolynomial(res)
    }

    pub fn commit(&self, srs: &SRS<BlsCurve>) -> MleCommit<BlsCurve> {
        commit(&self.0, srs)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();

        //Write polynomials
        for value in &self.0 {
            result
                .write(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        result
    }
}
