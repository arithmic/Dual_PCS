#![allow(unused)]
use bls381::scalar::Scalar;
use bls_curve::bls::AffinePoint;
use channel::Channel;
use crypto_bigint::{Encoding, Limb, Uint};
use grand_product_with_air::grand_product_common::GkrTranscript;
use kzg_fft::common::KZGFFTEvaluationProof;
use kzg_fourier_multilinear::common::{BatchMultilinearKZG2Proof, MultilinearKZG2Proof};
use polynomial::{MultPolynomial, Polynomial};
use preprocessing::SparseMetaData;
use std::io::Write;
use traits::traits::Field;

use crate::preprocessing;
extern crate bls381;
extern crate bls_curve;
extern crate channel;
extern crate crypto_bigint;
extern crate polynomial;
extern crate traits;

#[allow(unused)]
const SIX_INV: Scalar = Scalar(Uint {
    limbs: [
        Limb(15372286724512153601),
        Limb(1954008828163476649),
        Limb(9224930440102993583),
        Limb(6961264049553707793),
    ],
});
const TWO_INV: Scalar = Scalar(Uint {
    limbs: [
        Limb(9223372034707292161),
        Limb(12240451741123816959),
        Limb(1845609449319885826),
        Limb(4176758429732224676),
    ],
});
#[derive(Clone)]
pub struct GKREvaluations {
    pub rx_basis_evals_evaluation_at_z1: Scalar,
    pub ry_basis_evals_evaluation_at_z1: Scalar,
    pub rx_basis_evals_evaluation_at_g_z1: Scalar,
    pub ry_basis_evals_evaluation_at_g_z1: Scalar,
    pub final_ts_for_rows_evaluations_at_z1: Vec<Scalar>,
    pub final_ts_for_cols_evaluations_at_z1: Vec<Scalar>,
    pub final_ts_for_rows_evaluations_at_g_z1: Vec<Scalar>,
    pub final_ts_for_cols_evaluations_at_g_z1: Vec<Scalar>,
    pub rows_evaluations_at_z2: Vec<Scalar>,
    pub rows_evaluations_at_g_z2: Vec<Scalar>,
    pub cols_evaluations_at_z2: Vec<Scalar>,
    pub cols_evaluations_at_g_z2: Vec<Scalar>,
    pub e_rx_evaluations_at_z2: Vec<Scalar>,
    pub e_rx_evaluations_at_g_z2: Vec<Scalar>,
    pub e_ry_evaluations_at_z2: Vec<Scalar>,
    pub e_ry_evaluations_at_g_z2: Vec<Scalar>,
    pub read_ts_for_rows_evaluations_at_z2: Vec<Scalar>,
    pub read_ts_for_rows_evaluations_at_g_z2: Vec<Scalar>,
    pub read_ts_for_cols_evaluations_at_z2: Vec<Scalar>,
    pub read_ts_for_cols_evaluations_at_g_z2: Vec<Scalar>,
    pub idx_at_z1: Scalar,
    pub idx_at_g_z1: Scalar,
}
impl GKREvaluations {
    pub fn new(
        rx_basis_evals_evaluation_at_z1: Scalar,
        ry_basis_evals_evaluation_at_z1: Scalar,
        rx_basis_evals_evaluation_at_g_z1: Scalar,
        ry_basis_evals_evaluation_at_g_z1: Scalar,
        final_ts_for_rows_evaluations_at_z1: Vec<Scalar>,
        final_ts_for_cols_evaluations_at_z1: Vec<Scalar>,
        final_ts_for_rows_evaluations_at_g_z1: Vec<Scalar>,
        final_ts_for_cols_evaluations_at_g_z1: Vec<Scalar>,
        rows_evaluations_at_z2: Vec<Scalar>,
        rows_evaluations_at_g_z2: Vec<Scalar>,
        cols_evaluations_at_z2: Vec<Scalar>,
        cols_evaluations_at_g_z2: Vec<Scalar>,
        e_rx_evaluations_at_z2: Vec<Scalar>,
        e_rx_evaluations_at_g_z2: Vec<Scalar>,
        e_ry_evaluations_at_z2: Vec<Scalar>,
        e_ry_evaluations_at_g_z2: Vec<Scalar>,
        read_ts_for_rows_evaluations_at_z2: Vec<Scalar>,
        read_ts_for_rows_evaluations_at_g_z2: Vec<Scalar>,
        read_ts_for_cols_evaluations_at_z2: Vec<Scalar>,
        read_ts_for_cols_evaluations_at_g_z2: Vec<Scalar>,
        idx_at_z1: Scalar,
        idx_at_g_z1: Scalar,
    ) -> Self {
        Self {
            rx_basis_evals_evaluation_at_z1,
            ry_basis_evals_evaluation_at_z1,
            rx_basis_evals_evaluation_at_g_z1,
            ry_basis_evals_evaluation_at_g_z1,
            final_ts_for_rows_evaluations_at_z1,
            final_ts_for_cols_evaluations_at_z1,
            final_ts_for_rows_evaluations_at_g_z1,
            final_ts_for_cols_evaluations_at_g_z1,
            rows_evaluations_at_z2,
            rows_evaluations_at_g_z2,
            cols_evaluations_at_z2,
            cols_evaluations_at_g_z2,
            e_rx_evaluations_at_z2,
            e_rx_evaluations_at_g_z2,
            e_ry_evaluations_at_z2,
            e_ry_evaluations_at_g_z2,
            read_ts_for_rows_evaluations_at_z2,
            read_ts_for_rows_evaluations_at_g_z2,
            read_ts_for_cols_evaluations_at_z2,
            read_ts_for_cols_evaluations_at_g_z2,
            idx_at_z1,
            idx_at_g_z1,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write rx_basis_evals_evaluation_at_z1
        result
            .write_all(&self.rx_basis_evals_evaluation_at_z1.0.to_be_bytes())
            .expect("failed to write");

        //Write ry_basis_evals_evaluation_at_z1
        result
            .write_all(&self.ry_basis_evals_evaluation_at_z1.0.to_be_bytes())
            .expect("failed to write");

        //Write rx_basis_evals_evaluation_at_g_z1
        result
            .write_all(&self.rx_basis_evals_evaluation_at_g_z1.0.to_be_bytes())
            .expect("failed to write");

        //Write ry_basis_evals_evaluation_at_g_z1
        result
            .write_all(&self.ry_basis_evals_evaluation_at_g_z1.0.to_be_bytes())
            .expect("failed to write");

        //Write final_ts_for_rows_evaluations_at_z1
        for value in &self.final_ts_for_rows_evaluations_at_z1 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write final_ts_for_cols_evaluations_at_z1
        for value in &self.final_ts_for_cols_evaluations_at_z1 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write final_ts_for_rows_evaluations_at_g_z1
        for value in &self.final_ts_for_rows_evaluations_at_g_z1 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write final_ts_for_cols_evaluations_at_g_z1
        for value in &self.final_ts_for_cols_evaluations_at_g_z1 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write rows_evaluations_at_z2
        for value in &self.rows_evaluations_at_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write rows_evaluations_at_g_z2
        for value in &self.rows_evaluations_at_g_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write cols_evaluations_at_z2
        for value in &self.cols_evaluations_at_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write cols_evaluations_at_g_z2
        for value in &self.cols_evaluations_at_g_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write e_rx_evaluations_at_z2
        for value in &self.e_rx_evaluations_at_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write e_rx_evaluations_at_g_z2
        for value in &self.e_rx_evaluations_at_g_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write e_ry_evaluations_at_z2
        for value in &self.e_ry_evaluations_at_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write e_ry_evaluations_at_g_z2
        for value in &self.e_ry_evaluations_at_g_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write read_ts_for_rows_evaluations_at_z2
        for value in &self.read_ts_for_rows_evaluations_at_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write read_ts_for_rows_evaluations_at_g_z2
        for value in &self.read_ts_for_rows_evaluations_at_g_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write read_ts_for_cols_evaluations_at_z2
        for value in &self.read_ts_for_cols_evaluations_at_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write read_ts_for_cols_evaluations_at_g_z2
        for value in &self.read_ts_for_cols_evaluations_at_g_z2 {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write");
        }
        result
    }
}
#[derive(Clone)]
pub struct PreprocessCommits {
    pub final_ts_for_rows_commits: Vec<AffinePoint>,
    pub final_ts_for_cols_commits: Vec<AffinePoint>,
    pub rows_commits: Vec<AffinePoint>,
    pub cols_commits: Vec<AffinePoint>,
    pub read_ts_for_rows_commits: Vec<AffinePoint>,
    pub read_ts_for_cols_commits: Vec<AffinePoint>,
    pub val_commits: Vec<AffinePoint>,
    pub index_commit: AffinePoint,
}
impl PreprocessCommits {
    pub fn new(
        final_ts_for_rows_commits: Vec<AffinePoint>,
        final_ts_for_cols_commits: Vec<AffinePoint>,
        rows_commits: Vec<AffinePoint>,
        cols_commits: Vec<AffinePoint>,
        read_ts_for_rows_commits: Vec<AffinePoint>,
        read_ts_for_cols_commits: Vec<AffinePoint>,
        val_commits: Vec<AffinePoint>,
        index_commit: AffinePoint,
    ) -> Self {
        Self {
            final_ts_for_rows_commits,
            final_ts_for_cols_commits,
            rows_commits,
            cols_commits,
            read_ts_for_rows_commits,
            read_ts_for_cols_commits,
            val_commits,
            index_commit,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();

        //Write final_ts_for_rows_commits
        for commit in &self.final_ts_for_rows_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write bytes");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write bytes");
        }

        //Write final_ts_for_cols_commits
        for commit in &self.final_ts_for_cols_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write bytes");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write bytes");
        }

        //Write rows_commits
        for commit in &self.rows_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write bytes");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write bytes");
        }

        //Write cols_commits
        for commit in &self.cols_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write bytes");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write bytes");
        }

        //Write read_ts_for_rows_commits
        for commit in &self.read_ts_for_rows_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write bytes");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write bytes");
        }

        //Write read_ts_for_cols_commits
        for commit in &self.read_ts_for_cols_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write bytes");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write bytes");
        }
        //Write index_commit
        result
            .write_all(&self.index_commit.x.0.to_be_bytes())
            .expect("failed to write bytes");
        result
            .write_all(&self.index_commit.y.0.to_be_bytes())
            .expect("failed to write bytes");

        result
    }
}

#[derive(Clone)]
pub struct BatchSpartanSumCheckTranscript {
    pub polynomials: Vec<Polynomial>,
    pub random_points: Vec<Scalar>,
    pub e_rx: Vec<MultPolynomial>,
    pub e_ry: Vec<MultPolynomial>,
    pub val: Vec<MultPolynomial>,
}

impl BatchSpartanSumCheckTranscript {
    pub fn new(
        polynomials: Vec<Polynomial>,
        random_points: Vec<Scalar>,
        e_rx: Vec<MultPolynomial>,
        e_ry: Vec<MultPolynomial>,
        val: Vec<MultPolynomial>,
    ) -> BatchSpartanSumCheckTranscript {
        BatchSpartanSumCheckTranscript {
            polynomials,
            random_points,
            e_rx,
            e_ry,
            val,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write polynomials
        for poly in &self.polynomials {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write e_rx
        for poly in &self.e_rx {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write e_ry
        for poly in &self.e_ry {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write val
        for poly in &self.val {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }
        result
    }
}
#[derive(Clone)]
pub struct BatchSparseEvalProof {
    pub sum_check_transcript: BatchSpartanSumCheckTranscript,
    pub gkr_transcript1: GkrTranscript,
    pub gkr_transcript2: GkrTranscript,
    pub e_rx_evals_sum_check: Vec<Scalar>,
    pub e_ry_evals_sum_check: Vec<Scalar>,
    pub val_evals_sum_check: Vec<Scalar>,
    pub sum_check_eval_proof: MultilinearKZG2Proof,
    pub e_rx_commits: Vec<AffinePoint>,
    pub e_ry_commits: Vec<AffinePoint>,
    pub gkr_evaluations: GKREvaluations,
    pub proof1: KZGFFTEvaluationProof,
    pub proof2: KZGFFTEvaluationProof,
    pub proof3: KZGFFTEvaluationProof,
    pub proof4: KZGFFTEvaluationProof,
}

impl BatchSparseEvalProof {
    pub fn new(
        sum_check_transcript: BatchSpartanSumCheckTranscript,
        gkr_transcript1: GkrTranscript,
        gkr_transcript2: GkrTranscript,
        e_rx_evals_sum_check: Vec<Scalar>,
        e_ry_evals_sum_check: Vec<Scalar>,
        val_evals_sum_check: Vec<Scalar>,
        sum_check_eval_proof: MultilinearKZG2Proof,
        e_rx_commits: Vec<AffinePoint>,
        e_ry_commits: Vec<AffinePoint>,
        gkr_evaluations: GKREvaluations,
        proof1: KZGFFTEvaluationProof,
        proof2: KZGFFTEvaluationProof,
        proof3: KZGFFTEvaluationProof,
        proof4: KZGFFTEvaluationProof,
    ) -> BatchSparseEvalProof {
        Self {
            sum_check_transcript,
            gkr_transcript1,
            gkr_transcript2,
            e_rx_evals_sum_check,
            e_ry_evals_sum_check,
            val_evals_sum_check,
            sum_check_eval_proof,
            e_rx_commits,
            e_ry_commits,
            gkr_evaluations,
            proof1,
            proof2,
            proof3,
            proof4,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write sum_check_transcripts
        result
            .write_all(&self.sum_check_transcript.to_bytes())
            .expect("failed to write");

        //Write gkr_transcript1
        result
            .write_all(&self.gkr_transcript1.to_bytes())
            .expect("failed to write");

        //Write gkr_transcript2
        result
            .write_all(&self.gkr_transcript2.to_bytes())
            .expect("failed to write");

        //Write e_rx_evals_sum_check
        for value in &self.e_rx_evals_sum_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write e_ry_evals_sum_check
        for value in &self.e_ry_evals_sum_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write val_evals_sum_check
        for value in &self.val_evals_sum_check {
            result
                .write_all(&value.0.to_be_bytes())
                .expect("failed to write")
        }

        //Write sum_check_eval_proof
        result
            .write_all(&self.sum_check_eval_proof.to_bytes())
            .expect("failed to write");

        //Write e_rx_commits
        for commit in &self.e_rx_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write e_ry_commits
        for commit in &self.e_ry_commits {
            result
                .write_all(&commit.x.0.to_be_bytes())
                .expect("failed to write");
            result
                .write_all(&commit.y.0.to_be_bytes())
                .expect("failed to write");
        }

        //Write gkr_evaluations
        result
            .write_all(&self.gkr_evaluations.to_bytes())
            .expect("failed to write");

        //Write proof1
        result
            .write_all(&self.proof1.to_bytes())
            .expect("failed to write");
        result
            .write_all(&self.proof2.to_bytes())
            .expect("failed to write");
        result
            .write_all(&self.proof3.to_bytes())
            .expect("failed to write");
        result
            .write_all(&self.proof4.to_bytes())
            .expect("failed to write");
        result
    }
}
pub fn interpolate(eval: &mut Vec<Scalar>) {
    let t0 = eval[0] - eval[1].double();
    eval[1] = (-(eval[0] + eval[2]) - t0.double()) * TWO_INV;
    eval[2] = (t0 + eval[2]) * TWO_INV;
}
#[derive(Clone)]
#[allow(non_snake_case)]
pub struct InitialSumCheckTranscript {
    pub polynomials: Vec<Polynomial>,
    pub random_points: Vec<Scalar>,
    pub Az_claimed_val: Scalar,
    pub Bz_claimed_val: Scalar,
    pub Cz_claimed_val: Scalar,
    pub E_claimed_val: Scalar,
}
#[allow(non_snake_case)]
impl InitialSumCheckTranscript {
    pub fn new(
        polynomials: Vec<Polynomial>,
        random_points: Vec<Scalar>,
        Az_claimed_val: Scalar,
        Bz_claimed_val: Scalar,
        Cz_claimed_val: Scalar,
        E_claimed_val: Scalar,
    ) -> InitialSumCheckTranscript {
        InitialSumCheckTranscript {
            polynomials,
            random_points,
            Az_claimed_val,
            Bz_claimed_val,
            Cz_claimed_val,
            E_claimed_val,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write poplynomials
        for poly in &self.polynomials {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write Az_claimed_val
        result
            .write_all(&self.Az_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write Bz_claimed_val
        result
            .write_all(&self.Bz_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write Cz_claimed_val
        result
            .write_all(&self.Cz_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write E_claimed_val
        result
            .write_all(&self.E_claimed_val.0.to_be_bytes())
            .expect("failed to write");
        result
    }
}

#[derive(Clone)]
#[allow(non_snake_case)]
pub struct ParSumCheckTranscript {
    pub polynomials: Vec<Polynomial>,
    pub random_points: Vec<Scalar>,
    pub A_claimed_val: Scalar,
    pub B_claimed_val: Scalar,
    pub C_claimed_val: Scalar,
    pub Z_claimed_val: Scalar,
}

#[allow(non_snake_case)]
impl ParSumCheckTranscript {
    pub fn new(
        polynomials: Vec<Polynomial>,
        random_points: Vec<Scalar>,
        A_claimed_val: Scalar,
        B_claimed_val: Scalar,
        C_claimed_val: Scalar,
        Z_claimed_val: Scalar,
    ) -> ParSumCheckTranscript {
        ParSumCheckTranscript {
            polynomials,
            random_points,
            A_claimed_val,
            B_claimed_val,
            C_claimed_val,
            Z_claimed_val,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        //Write poplynomials
        for poly in &self.polynomials {
            result.write_all(&poly.to_bytes()).expect("failed to write");
        }

        //Write A_claimed_val
        result
            .write_all(&self.A_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write B_claimed_val
        result
            .write_all(&self.B_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write C_claimed_val
        result
            .write_all(&self.C_claimed_val.0.to_be_bytes())
            .expect("failed to write");

        //Write Z_claimed_val
        result
            .write_all(&self.Z_claimed_val.0.to_be_bytes())
            .expect("failed to write");
        result
    }
}

#[derive(Clone)]
#[allow(non_snake_case, non_camel_case_types)]
pub struct eR1CStranscript {
    pub first_sum_check_transcript: InitialSumCheckTranscript,
    pub par_sum_check_transcript: ParSumCheckTranscript,
    pub BatchProof: BatchSparseEvalProof,
    pub e_w_eval_proof: BatchMultilinearKZG2Proof,
    pub rx_basis_commit: AffinePoint,
    pub ry_basis_commit: AffinePoint,
}

#[allow(non_snake_case)]
impl eR1CStranscript {
    pub fn new(
        first_sum_check_transcript: InitialSumCheckTranscript,
        par_sum_check_transcript: ParSumCheckTranscript,
        BatchProof: BatchSparseEvalProof,
        e_w_eval_proof: BatchMultilinearKZG2Proof,
        rx_basis_commit: AffinePoint,
        ry_basis_commit: AffinePoint,
    ) -> eR1CStranscript {
        eR1CStranscript {
            first_sum_check_transcript,
            par_sum_check_transcript,
            BatchProof,
            e_w_eval_proof,
            rx_basis_commit,
            ry_basis_commit,
        }
    }
    pub fn to_bytes(&self) -> f64 {
        let mut result = Vec::new();

        //Write first_sum_check_transcript
        result
            .write_all(&self.first_sum_check_transcript.to_bytes())
            .expect("failed to write");

        //Write par_sum_check_transcript
        result
            .write_all(&self.par_sum_check_transcript.to_bytes())
            .expect("failed to write");

        //Write BatchProof
        result
            .write_all(&self.BatchProof.to_bytes())
            .expect("failed to write");

        //Write e_we_eval_proof
        result
            .write_all(&self.e_w_eval_proof.to_bytes())
            .expect("failed to write");

        result.len() as f64 / 1024f64
    }
}

#[allow(non_snake_case, non_camel_case_types)]
#[derive(Clone)]
pub struct eR1CSmetadata {
    pub A: SparseMetaData,
    pub B: SparseMetaData,
    pub C: SparseMetaData,
}

#[allow(non_snake_case)]
impl eR1CSmetadata {
    pub fn new(A: SparseMetaData, B: SparseMetaData, C: SparseMetaData) -> eR1CSmetadata {
        eR1CSmetadata { A, B, C }
    }
}
//...........
// CODE  for evaluating polynomial at points
//.............
pub fn eval(p: &[Scalar], x: Scalar) -> Scalar {
    // Horner evaluation
    p.iter()
        .rev()
        .fold(Scalar::ZERO, |acc, &coeff| acc * x + coeff)
}
//TODO:- Reseed with cpmmits
pub fn reseed_with_commits(gkr_commits: &PreprocessCommits, channel: &mut Channel) {
    channel.reseed_with_affine_point(&gkr_commits.final_ts_for_rows_commits);
    channel.reseed_with_affine_point(&gkr_commits.final_ts_for_cols_commits);
    channel.reseed_with_affine_point(&gkr_commits.rows_commits);
    channel.reseed_with_affine_point(&gkr_commits.cols_commits);
    channel.reseed_with_affine_point(&gkr_commits.read_ts_for_rows_commits);
    channel.reseed_with_affine_point(&gkr_commits.read_ts_for_cols_commits);
}
#[allow(non_snake_case, non_camel_case_types)]
#[derive(Clone)]
pub struct eR1CSCommitments {
    pub preprocess_commits: PreprocessCommits,
    pub E_Commit: AffinePoint,
    pub W_Commit: AffinePoint,
}

#[allow(non_snake_case)]
impl eR1CSCommitments {
    pub fn new(
        preprocess_commits: PreprocessCommits,
        E_Commit: AffinePoint,
        W_Commit: AffinePoint,
    ) -> eR1CSCommitments {
        eR1CSCommitments {
            preprocess_commits,
            E_Commit,
            W_Commit,
        }
    }
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut result = Vec::new();
        result
            .write_all(&self.preprocess_commits.to_bytes())
            .expect("failed to write");

        // Write E_commit
        result
            .write_all(&self.E_Commit.x.0.to_be_bytes())
            .expect("failed to write");
        result
            .write_all(&self.E_Commit.y.0.to_be_bytes())
            .expect("failed to write");

        // Write W_commit
        result
            .write_all(&self.W_Commit.x.0.to_be_bytes())
            .expect("failed to write");
        result
            .write_all(&self.W_Commit.y.0.to_be_bytes())
            .expect("failed to write");
        result
    }
}
