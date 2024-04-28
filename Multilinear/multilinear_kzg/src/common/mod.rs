#![allow(non_snake_case)]
use crypto_bigint::U256;
use curve_traits::{projectivepoint::Scalar, CurveArithmetic, CurveParams, ProjectivePoint};
use pairing::pairing_traits::Pairing;
use rayon::{
    iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator},
    slice::ParallelSlice,
};
use std::{fmt::Display, io::Write};
use traits::traits::Field;

#[derive(Debug, Clone, Copy)]
pub struct MleCommit<C: Pairing> {
    pub commitment: ProjectivePoint<C::G1Curve>,
}

impl<C: Pairing> Display for MleCommit<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} {} {}",
            self.commitment.x, self.commitment.y, self.commitment.z
        )
    }
}

impl<C: Pairing> MleCommit<C> {
    pub fn to_string(&self) -> String {
        format!(
            "{} {} {}",
            self.commitment.x, self.commitment.y, self.commitment.z
        )
    }

    pub fn from_string(str: &str) -> Option<MleCommit<C>> {
        let mut split = str.split(" ");
        let state1 = split.next()?;
        let state2 = split.next()?;
        let state3 = split.next()?;
        Some(MleCommit::<C> {
            commitment: ProjectivePoint::<C::G1Curve> {
                x: <<C as Pairing>::G1Curve as CurveParams>::FieldElement::from(U256::from_be_hex(
                    state1,
                )),
                y: <<C as Pairing>::G1Curve as CurveParams>::FieldElement::from(U256::from_be_hex(
                    state2,
                )),
                z: <<C as Pairing>::G1Curve as CurveParams>::FieldElement::from(U256::from_be_hex(
                    state3,
                )),
            },
        })
    }
}

#[derive(Debug, Clone)]
pub struct MleEvalProof<C: Pairing> {
    pub evaluation: Scalar<C::G1Curve>,
    pub witnesses: Vec<ProjectivePoint<C::G1Curve>>,
}
impl<C: Pairing> MleEvalProof<C> {
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        result
            .write_all(&self.evaluation.to_curve_bytes())
            .expect("failed to write");
        for witness in &self.witnesses {
            let affine_point = witness.to_affine();
            result
                .write_all(&affine_point.x.to_curve_bytes())
                .expect("failed to write");
            result
                .write_all(&affine_point.y.to_curve_bytes())
                .expect("failed to write");
        }
        result.len() as f64 / 1024f64
    }
}

pub fn setup<C: Pairing>(toxic_waste: Vec<Scalar<C::G1Curve>>) -> (SRS<C>, VerificationKey<C>) {
    let g_1 = ProjectivePoint::<C::G1Curve>::GENERATOR;
    let g_2 = ProjectivePoint::<C::G2Curve>::GENERATOR;

    let g1_exponents: Vec<<<C as Pairing>::G1Curve as CurveArithmetic>::Scalar> =
        compute_fourier_bases::<C::G1Curve>(&toxic_waste);
    let srs_G1: Vec<_> = g1_exponents.par_iter().map(|X| g_1.mul(*X)).collect();
    let ver_key: Vec<ProjectivePoint<<C as Pairing>::G2Curve>> = toxic_waste
        .par_iter()
        .map(|r_i| g_2.mul(Scalar::<C::G2Curve>::from_words(&r_i.to_words())))
        .collect();
    (
        SRS {
            g_1: srs_G1,
            variables: toxic_waste.len(),
        },
        VerificationKey { key: ver_key },
    )
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SRS<C: Pairing> {
    pub g_1: Vec<ProjectivePoint<C::G1Curve>>,
    pub variables: usize,
}

impl<C: Pairing> SRS<C> {
    pub fn from_string(str: &str) -> SRS<C> {
        let mut split = str.split(" ");
        let field_elements: Vec<<<C as Pairing>::G1Curve as CurveParams>::FieldElement> = split
            .par_bridge()
            .map(|x| {
                <<C as Pairing>::G1Curve as CurveParams>::FieldElement::from(U256::from_be_hex(x))
            })
            .collect();
        let projective_points: Vec<ProjectivePoint<C::G1Curve>> = field_elements
            .par_chunks(3)
            .map(|points| ProjectivePoint::<C::G1Curve> {
                x: points[0],
                y: points[1],
                z: points[2],
            })
            .collect();

        let variables = projective_points.len().trailing_zeros() as usize;
        SRS {
            g_1: projective_points,
            variables: variables,
        }
    }
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        for value in &self.g_1 {
            let affine_point = value.to_affine();
            result
                .write_all(&affine_point.x.to_curve_bytes())
                .expect("failed to write");
            result
                .write_all(&affine_point.y.to_curve_bytes())
                .expect("failed to write");
        }
        result
            .write_all(&self.variables.to_be_bytes())
            .expect("failed to write");
        result.len() as f64 / 1024f64
    }
}

impl<C: Pairing> Display for SRS<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Ok(for p in self.g_1.iter() {
            write!(f, "{}  {}  {} ", p.x, p.y, p.z);
        })
    }
}
#[derive(Clone)]
pub struct VerificationKey<C: Pairing> {
    pub key: Vec<ProjectivePoint<C::G2Curve>>,
}
impl<C: Pairing> VerificationKey<C> {
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        for value in &self.key {
            let affine_value = value.to_affine();
            result
                .write_all(&affine_value.x.to_curve_bytes())
                .expect("failed to write");
            result
                .write_all(&affine_value.x.to_curve_bytes())
                .expect("failed to write");
        }
        result.len() as f64 / 1024f64
    }
}
pub struct BatchEvalProof {}

pub fn compute_fourier_bases<C: CurveArithmetic + CurveParams>(
    r: &Vec<Scalar<C>>,
) -> Vec<Scalar<C>> {
    //Initialize fc_eq with (1- r[0]) and r[0]
    let mut fc_eq = [<Scalar<C> as Field>::ONE - r[0], r[0]].to_vec();
    //Iterate over the length of the r vector
    for k in 1..r.len() {
        let temp = fc_eq;
        //initialize fc_eq of double size with zero
        fc_eq = vec![<Scalar<C> as Field>::ZERO; temp.len() * 2];
        // if temp.len() < (2 << 8) {
        for iter in 0..temp.len() {
            fc_eq[2 * iter] = temp[iter] * (<Scalar<C> as Field>::ONE - r[k]);
            fc_eq[(2 * iter) + 1] = temp[iter] * (r[k]);
        }
    }
    fc_eq
}
