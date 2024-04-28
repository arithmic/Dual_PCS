use std::{fmt::Display, io::Write};

use bls381::{
    fp::{self, Fp},
    fp2::Fp2,
    scalar::Scalar,
};
use bls_curve::{bls::BlsCurve, fp2bls::G2BlsCurve};
use curve_traits::AffinePoint;

#[derive(Clone)]
pub struct KZGFFTSetup {
    pub prover_key: Vec<curve_traits::ProjectivePoint<BlsCurve>>,
    pub verifier_key: curve_traits::AffinePoint<G2BlsCurve>,
}
impl KZGFFTSetup {
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();
        //Write prover key
        for p_k in self.prover_key.clone() {
            let affine_pk = p_k.to_affine();
            result
                .write_all(&fp::Fp::to_bytes(&affine_pk.x))
                .expect("failed to write bytes");
            result
                .write_all(&fp::Fp::to_bytes(&affine_pk.y))
                .expect("failed to write bytes");
        }

        //Write Verifier key
        result
            .write_all(&Fp2::<Fp>::to_bytes(&self.verifier_key.x))
            .expect("failed to write bytes");
        result
            .write_all(&Fp2::<Fp>::to_bytes(&self.verifier_key.y))
            .expect("failed to write bytes");

        result.len() as f64 / 1024f64
    }
}
#[derive(Clone)]
pub struct KZGFFTEvaluationProof {
    pub commitment_of_evaluations_of_q: AffinePoint<BlsCurve>,
    pub evaluation_at_z: Scalar,
}
impl KZGFFTEvaluationProof {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write commitment_of_evaluations_of_q
        result
            .write_all(&fp::Fp::to_bytes(&self.commitment_of_evaluations_of_q.x))
            .expect("failed to write bytes");
        result
            .write_all(&fp::Fp::to_bytes(&self.commitment_of_evaluations_of_q.y))
            .expect("failed to write bytes");

        //Write evaluation_at_z
        result
            .write_all(&Scalar::to_bytes(&self.evaluation_at_z))
            .expect("failed to write bytes");

        result
    }
}
#[derive(Clone)]
pub struct KZGFFTDegreeBoundSetup {
    pub prover_keys: Vec<Vec<curve_traits::ProjectivePoint<BlsCurve>>>,
    pub verifier_key: Vec<curve_traits::AffinePoint<G2BlsCurve>>,
}
impl KZGFFTDegreeBoundSetup {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write prover key
        for prover_key in self.prover_keys.clone() {
            for p_k in prover_key {
                let affine_pk = p_k.to_affine();
                result
                    .write_all(&fp::Fp::to_bytes(&affine_pk.x))
                    .expect("failed to write bytes");
                result
                    .write_all(&fp::Fp::to_bytes(&affine_pk.y))
                    .expect("failed to write bytes");
            }
        }

        //Write Verifier key
        for v_k in self.verifier_key.clone() {
            result
                .write_all(&Fp2::<Fp>::to_bytes(&v_k.x))
                .expect("failed to write bytes");
            result
                .write_all(&Fp2::<Fp>::to_bytes(&v_k.y))
                .expect("failed to write bytes");
        }
        result
    }
}
impl Display for KZGFFTDegreeBoundSetup {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //Write prover key
        for outer_idx in 0..self.prover_keys.len() - 1 {
            for inner_idx in 0..self.prover_keys[outer_idx].len() - 1 {
                write!(
                    f,
                    "{} {} {},",
                    self.prover_keys[outer_idx][inner_idx].x,
                    self.prover_keys[outer_idx][inner_idx].y,
                    self.prover_keys[outer_idx][inner_idx].z,
                )
                .expect("failed to write");
            }
            write!(
                f,
                "{} {} {}$",
                self.prover_keys[outer_idx][self.prover_keys[outer_idx].len() - 1].x,
                self.prover_keys[outer_idx][self.prover_keys[outer_idx].len() - 1].y,
                self.prover_keys[outer_idx][self.prover_keys[outer_idx].len() - 1].z,
            )
            .expect("failed to write");
        }
        for inner_idx in 0..self.prover_keys[self.prover_keys.len() - 1].len() - 1 {
            write!(
                f,
                "{} {} {},",
                self.prover_keys[self.prover_keys.len() - 1][inner_idx].x,
                self.prover_keys[self.prover_keys.len() - 1][inner_idx].y,
                self.prover_keys[self.prover_keys.len() - 1][inner_idx].z,
            )
            .expect("failed to write");
        }
        write!(
            f,
            "{} {} {}#",
            self.prover_keys[self.prover_keys.len() - 1]
                [self.prover_keys[self.prover_keys.len() - 1].len() - 1]
                .x,
            self.prover_keys[self.prover_keys.len() - 1]
                [self.prover_keys[self.prover_keys.len() - 1].len() - 1]
                .y,
            self.prover_keys[self.prover_keys.len() - 1]
                [self.prover_keys[self.prover_keys.len() - 1].len() - 1]
                .z,
        )
        .expect("failed to write");
        //Write verifier key
        for idx in 0..self.verifier_key.len() - 1 {
            write!(
                f,
                "{}&{} {}&{},",
                self.verifier_key[idx].x.c0,
                self.verifier_key[idx].x.c1,
                self.verifier_key[idx].y.c0,
                self.verifier_key[idx].y.c1
            )
            .expect("failed to write");
        }
        write!(
            f,
            "{}&{} {}&{}",
            self.verifier_key[self.verifier_key.len() - 1].x.c0,
            self.verifier_key[self.verifier_key.len() - 1].x.c1,
            self.verifier_key[self.verifier_key.len() - 1].y.c0,
            self.verifier_key[self.verifier_key.len() - 1].y.c1,
        )
    }
}
impl KZGFFTDegreeBoundSetup {
    pub fn get_setup(&self, degree: usize) -> KZGFFTSetup {
        let folded_prover_key = self.prover_keys
            [(self.prover_keys[0].len().trailing_zeros() - degree.trailing_zeros()) as usize]
            .clone();

        let verifier_key = self.verifier_key
            [(self.prover_keys[0].len().trailing_zeros() - degree.trailing_zeros()) as usize];
        assert_eq!(folded_prover_key.len(), degree);
        KZGFFTSetup {
            prover_key: folded_prover_key,
            verifier_key,
        }
    }
}
