use std::{fmt::Display, io::Write};

use bls381::{fp, scalar::Scalar};
use bls_curve::bls::BlsCurve;
use crypto_bigint::Encoding;
use curve_traits::AffinePoint;
use kzg_fft::common::{KZGFFTDegreeBoundSetup, KZGFFTSetup};
use traits::traits::PrimeField;
pub const OFFSET: Scalar = Scalar::GENERATOR;

#[derive(Clone)]
pub struct KZGFourierDegreeBoundSetup {
    pub setup: KZGFFTDegreeBoundSetup,
    pub phi_commitments: Vec<AffinePoint<BlsCurve>>,
    pub degree_bound: usize,
}

impl KZGFourierDegreeBoundSetup {
    pub fn get_setup(&self, degree: usize) -> KZGFourierSetup {
        let setup = self.setup.get_setup(degree);
        let phi_commitment = self.phi_commitments[1..degree.trailing_zeros() as usize].to_vec();
        KZGFourierSetup {
            setup,
            phi_commitment,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write setup
        result
            .write_all(&self.setup.to_bytes())
            .expect("failed to write bytes");
        println!("size of Kzg2setup is {:?}", result.len() as f64 / 1024f64);
        //Write phi commitments
        for commit in &self.phi_commitments {
            result
                .write_all(&fp::Fp::to_bytes(&commit.x))
                .expect("failed to write bytes");
            result
                .write_all(&fp::Fp::to_bytes(&commit.y))
                .expect("failed to write bytes");
        }
        //write degree bound
        result
            .write_all(&self.degree_bound.to_be_bytes())
            .expect("failed to write");
        result
    }
}
impl Display for KZGFourierDegreeBoundSetup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //Write setup
        write!(f, "{}*", self.setup).expect("failed to write");
        for idx in 0..self.phi_commitments.len() - 1 {
            write!(
                f,
                "{} {},",
                self.phi_commitments[idx].x, self.phi_commitments[idx].y
            )
            .expect("failed to write");
        }

        //Write phi_commitments
        write!(
            f,
            "{} {}*",
            self.phi_commitments[self.phi_commitments.len() - 1].x,
            self.phi_commitments[self.phi_commitments.len() - 1].y
        )
        .expect("failed to write");

        //Write degree_bound
        write!(f, "{}", Scalar::from(self.degree_bound as u64))
    }
}
#[derive(Clone)]
pub struct KZGFourierSetup {
    pub setup: KZGFFTSetup,
    pub phi_commitment: Vec<AffinePoint<BlsCurve>>,
}
#[derive(Clone)]
pub struct EvaluationOfEvenOddPolys {
    pub const_poly: Vec<Scalar>,
    pub even_polys_evaluations: Vec<Vec<Scalar>>,
    pub odd_polys_evaluations: Vec<Vec<Scalar>>,
}
impl EvaluationOfEvenOddPolys {
    pub fn new(
        const_poly: Vec<Scalar>,
        even_polys_evaluations: Vec<Vec<Scalar>>,
        odd_polys_evaluations: Vec<Vec<Scalar>>,
    ) -> Self {
        Self {
            const_poly,
            even_polys_evaluations,
            odd_polys_evaluations,
        }
    }
}
#[derive(Clone)]
pub struct EvenOddPolys {
    pub const_poly: Vec<Scalar>,
    pub even_polys: Vec<Vec<Scalar>>,
    pub odd_polys: Vec<Vec<Scalar>>,
}
impl EvenOddPolys {
    pub fn new(
        const_poly: Vec<Scalar>,
        even_polys_evaluations: Vec<Vec<Scalar>>,
        odd_polys_evaluations: Vec<Vec<Scalar>>,
    ) -> Self {
        Self {
            const_poly,
            even_polys: even_polys_evaluations,
            odd_polys: odd_polys_evaluations,
        }
    }
}
#[derive(Clone)]
pub struct CommitmentsOfEvenOddPolys {
    pub even_polys_commits: Vec<AffinePoint<BlsCurve>>,
    pub odd_polys_commits: Vec<AffinePoint<BlsCurve>>,
}
impl CommitmentsOfEvenOddPolys {
    pub fn new(
        even_polys_commitment: Vec<AffinePoint<BlsCurve>>,
        odd_polys_commitment: Vec<AffinePoint<BlsCurve>>,
    ) -> Self {
        Self {
            even_polys_commits: even_polys_commitment,
            odd_polys_commits: odd_polys_commitment,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write even polys commits
        for commit in &self.even_polys_commits {
            result
                .write_all(&fp::Fp::to_bytes(&commit.x))
                .expect("failed to write bytes");
            result
                .write_all(&fp::Fp::to_bytes(&commit.y))
                .expect("failed to write bytes");
        }
        //Write odd polys commits
        for commit in &self.odd_polys_commits {
            result
                .write_all(&fp::Fp::to_bytes(&commit.x))
                .expect("failed to write bytes");
            result
                .write_all(&fp::Fp::to_bytes(&commit.y))
                .expect("failed to write bytes");
        }
        result
    }
}

#[derive(Clone)]
pub struct EvaluationsAtRandomPoint {
    pub evaluation_of_f: Scalar,
    pub even_polys_evaluations: Vec<Scalar>,
    pub odd_polys_evaluations: Vec<Scalar>,
}
impl EvaluationsAtRandomPoint {
    pub fn new(
        evaluation_of_f: Scalar,
        even_polys_evaluations: Vec<Scalar>,
        odd_polys_evaluations: Vec<Scalar>,
    ) -> EvaluationsAtRandomPoint {
        Self {
            evaluation_of_f,
            even_polys_evaluations,
            odd_polys_evaluations,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write evaluation of f
        result
            .write_all(&self.evaluation_of_f.0.to_be_bytes())
            .expect("failed to write");

        //Write even_polys_evaluations
        for eval in &self.even_polys_evaluations {
            result
                .write_all(&eval.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write odd_polys_evaluations
        for eval in &self.odd_polys_evaluations {
            result
                .write_all(&eval.0.to_be_bytes())
                .expect("failed to write")
        }

        result
    }
}

#[derive(Clone)]
pub struct BatchMultilinearKZG2Proof {
    pub evaluations: Vec<Scalar>,
    pub commitments_of_even_ood_polys: Vec<CommitmentsOfEvenOddPolys>,
    pub evaluations_at_z: Vec<EvaluationsAtRandomPoint>,
    pub evaluations_at_s: Vec<EvaluationsAtRandomPoint>,
    pub commitment_to_r: AffinePoint<BlsCurve>,
    pub shi_commit: AffinePoint<BlsCurve>,
    pub q_not_polys: Vec<Vec<Scalar>>,
    pub phi_k_evaluation_at_z: Vec<Scalar>,
    pub phi_k_evaluation_at_s: Vec<Scalar>,
}
impl BatchMultilinearKZG2Proof {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write evaluation
        for eval in &self.evaluations {
            result
                .write_all(&eval.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write commitments_of_even_ood_polys
        for commitments in &self.commitments_of_even_ood_polys {
            result
                .write_all(&commitments.to_bytes())
                .expect("failed to write");
        }

        //Write evaluations_at_z
        for evaluation in &self.evaluations_at_z {
            result
                .write(&evaluation.to_bytes())
                .expect("failed to write");
        }

        //Write evaluations_at_s
        for evaluation in &self.evaluations_at_s {
            result
                .write(&evaluation.to_bytes())
                .expect("failed to write");
        }
        //Write commitment_to_r
        result
            .write_all(&fp::Fp::to_bytes(&self.commitment_to_r.x))
            .expect("failed to write bytes");
        result
            .write_all(&fp::Fp::to_bytes(&self.commitment_to_r.y))
            .expect("failed to write bytes");

        //Write shi_commit
        result
            .write_all(&fp::Fp::to_bytes(&self.shi_commit.x))
            .expect("failed to write bytes");
        result
            .write_all(&fp::Fp::to_bytes(&self.shi_commit.y))
            .expect("failed to write bytes");

        //Write phi_k_evaluation_at_z
        for eval in &self.phi_k_evaluation_at_z {
            result
                .write_all(&eval.0.to_be_bytes())
                .expect("failed to write")
        }
        //Write phi_k_evaluation_at_s
        for eval in &self.phi_k_evaluation_at_s {
            result
                .write_all(&eval.0.to_be_bytes())
                .expect("failed to write")
        }
        //Write q_zero
        for poly in &self.q_not_polys {
            for eval in poly {
                result
                    .write_all(&eval.0.to_be_bytes())
                    .expect("failed to write")
            }
        }
        result
    }
}
#[derive(Clone)]
pub struct MultilinearKZG2Proof {
    evaluation: Scalar,
    commitments_of_even_ood_polys: CommitmentsOfEvenOddPolys,
    evaluations_at_z: EvaluationsAtRandomPoint,
    evaluations_at_s: EvaluationsAtRandomPoint,
    commitment_to_r: AffinePoint<BlsCurve>,
    shi_commit: AffinePoint<BlsCurve>,
    phi_k_evaluation_at_z: Vec<Scalar>,
    phi_k_evaluation_at_s: Vec<Scalar>,
    q_not_poly: Vec<Scalar>,
}
impl MultilinearKZG2Proof {
    pub fn new(
        evaluation: Scalar,
        commitments_of_even_ood_polys: CommitmentsOfEvenOddPolys,
        evaluations_at_z: EvaluationsAtRandomPoint,
        commitment_to_r: AffinePoint<BlsCurve>,
        evaluations_at_s: EvaluationsAtRandomPoint,
        shi_commit: AffinePoint<BlsCurve>,
        phi_k_evaluation_at_z: Vec<Scalar>,
        phi_k_evaluation_at_s: Vec<Scalar>,
        q_not_poly: Vec<Scalar>,
    ) -> Self {
        Self {
            evaluation,
            commitments_of_even_ood_polys,
            evaluations_at_z,
            commitment_to_r,
            evaluations_at_s,
            shi_commit,
            phi_k_evaluation_at_z,
            phi_k_evaluation_at_s,
            q_not_poly,
        }
    }
    pub fn get_evaluation(&self) -> Scalar {
        self.evaluation
    }
    pub fn get_commitments_of_even_odd_polys(&self) -> &CommitmentsOfEvenOddPolys {
        &self.commitments_of_even_ood_polys
    }
    pub fn get_evaluations_at_z(&self) -> &EvaluationsAtRandomPoint {
        &self.evaluations_at_z
    }
    pub fn get_commitment_to_r(&self) -> AffinePoint<BlsCurve> {
        self.commitment_to_r
    }
    pub fn get_evaluations_at_s(&self) -> &EvaluationsAtRandomPoint {
        &self.evaluations_at_s
    }
    pub fn get_shi_commit(&self) -> AffinePoint<BlsCurve> {
        self.shi_commit
    }
    pub fn get_phi_k_at_z(&self) -> &Vec<Scalar> {
        &self.phi_k_evaluation_at_z
    }
    pub fn get_phi_k_at_s(&self) -> &Vec<Scalar> {
        &self.phi_k_evaluation_at_s
    }
    pub fn get_q_not(&self) -> &Vec<Scalar> {
        &self.q_not_poly
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write evaluation
        result
            .write_all(&self.evaluation.0.to_be_bytes())
            .expect("failed to write");

        //Write commitments_of_even_ood_polys
        result
            .write_all(&self.commitments_of_even_ood_polys.to_bytes())
            .expect("failed to write");

        //Write evaluations_at_z
        result
            .write(&self.evaluations_at_z.to_bytes())
            .expect("failed to write");

        //Write evaluations_at_s
        result
            .write(&self.evaluations_at_s.to_bytes())
            .expect("failed to write");

        //Write commitment_to_r
        result
            .write_all(&fp::Fp::to_bytes(&self.commitment_to_r.x))
            .expect("failed to write bytes");
        result
            .write_all(&fp::Fp::to_bytes(&self.commitment_to_r.y))
            .expect("failed to write bytes");

        //Write shi_commit
        result
            .write_all(&fp::Fp::to_bytes(&self.shi_commit.x))
            .expect("failed to write bytes");
        result
            .write_all(&fp::Fp::to_bytes(&self.shi_commit.y))
            .expect("failed to write bytes");

        //Write phi_k_evaluation_at_ss
        for eval in &self.phi_k_evaluation_at_z {
            result
                .write_all(&eval.0.to_be_bytes())
                .expect("failed to write")
        }
        //Write phi_k_evaluation_at_s
        for eval in &self.phi_k_evaluation_at_s {
            result
                .write_all(&eval.0.to_be_bytes())
                .expect("failed to write")
        }
        //Write q_zero
        for eval in &self.q_not_poly {
            result
                .write_all(&eval.0.to_be_bytes())
                .expect("failed to write")
        }
        result
    }
}
