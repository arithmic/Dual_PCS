use channel::Channel;
use common::CommitmentsOfEvenOddPolys;

pub(crate) fn reseed_with_even_odd_polys(
    channel: &mut Channel,
    commitments: CommitmentsOfEvenOddPolys,
) {
    commitments
        .even_polys_commits
        .into_iter()
        .for_each(|commit| channel.reseed_with_affine_point(&[commit].to_vec()));

    commitments.odd_polys_commits.into_iter().for_each(
        |commit: curve_traits::AffinePoint<bls_curve::bls::BlsCurve>| {
            channel.reseed_with_affine_point(&[commit].to_vec())
        },
    );
}
