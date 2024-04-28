#![allow(unused)]
use crate::grand_product_common;
use channel::Channel;
use grand_product_common::Commitments;

pub(crate) fn reseed_with_commits(commitments: &Commitments, channel: &mut Channel) {
    channel.reseed_with_affine_point(commitments.get_B_commits());
    channel.reseed_with_affine_point(commitments.get_C_commits());
}
