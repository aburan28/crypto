#![recursion_limit = "512"]
//! Public fixture interface for the independently complete rank and shared-log implementations.
#[path = "koblitz_rank_fixture/rank.rs"]
mod rank;
#[path = "koblitz_rank_fixture/shared.rs"]
mod shared;

fn main() {
    let arguments: Vec<_> = std::env::args().collect();
    let mode = arguments.get(8).map(String::as_str).unwrap_or("pointwise");
    let shared_mode = mode.starts_with("pair_pair_dual_")
        || mode.starts_with("pair_pair_guided_")
        || std::env::var("KIC_SHARED_FACTOR_LOG_PRECOMPUTATION").as_deref() == Ok("1")
        || std::env::var("KIC_SHARED_PUBLIC_FIXTURE_DOMAIN").as_deref() == Ok("1")
        || std::env::var_os("KIC_X_FILTER_BITS_PER_EXPECTED").is_some()
        || std::env::var_os("KIC_X_FILTER_HASHES").is_some();
    if shared_mode {
        assert!(arguments.len() <= 10, "the shared-log interface takes its public fixture domain through KIC_SHARED_PUBLIC_FIXTURE_DOMAIN");
        shared::run();
    } else {
        rank::run();
    }
}
