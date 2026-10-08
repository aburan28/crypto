//! Layout receipt for research/f6_n83_compact_pairs_20261005.
use crypto_lib::binary_ecc::BinaryPoint;

#[allow(dead_code)]
struct Before {
    point: BinaryPoint,
    pair: (usize, usize),
    neg_pair: (usize, usize),
}

#[allow(dead_code)]
struct After {
    point: Option<(u128, u128)>,
    pair: (u32, u32),
    neg_pair: (u32, u32),
}

fn main() {
    println!(
        "old_entry_bytes={} new_entry_bytes={}",
        std::mem::size_of::<Before>(),
        std::mem::size_of::<After>()
    );
}
