use channel::Channel;

use crate::gkr_common::Commitments;

pub(crate) fn reseed_commits(commitments: Commitments, channel: &mut Channel) {
    channel.reseed_mle_commit(commitments.A_commits[1..].to_vec());
    channel.reseed_mle_commit(commitments.B_commits);
    channel.reseed_mle_commit(commitments.C_commits);
}
