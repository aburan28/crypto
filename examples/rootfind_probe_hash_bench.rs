//! Probe-path microbench: shared `pair_filter_hash` vs hashing twice per miss.
//!
//! Models the walking `m = 2` descent reject on **misses** (the common case).
//! Not a security claim. Public Koblitz fixtures only.

use std::time::Instant;

use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_frobenius_factor_base, KoblitzCurve, PairSumTable,
};

fn mix(key: u64) -> u64 {
    let mut h = key.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h.wrapping_mul(0xc4ce_b9fe_1a85_ec53)
}

fn main() {
    let n: u32 = std::env::var("PROBE_N")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(13);
    let probes: u64 = std::env::var("PROBE_COUNT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(5_000_000);
    let kc = KoblitzCurve::new(0, n).expect("curve");
    let fb = build_frobenius_factor_base(&kc, 0).expect("fb");
    let pair = PairSumTable::build(&kc, &fb).expect("table");

    // Synthetic packed keys that are almost never pair sums: high bit set,
    // odd abscissa lane, so the presence filter rejects without recovery.
    let mut keys = Vec::with_capacity(4096);
    for k in 0..4096u64 {
        keys.push((k.wrapping_mul(0x9e37_79b9_7f4a_7c15) << 1) | 1 | (1u64 << 63));
    }

    let t0 = Instant::now();
    let mut hits = 0u64;
    for i in 0..probes {
        let key = keys[(i as usize) & 4095];
        let h = mix(key);
        pair.prefetch_key_hash(h);
        if pair.contains_key_hash(key, h) {
            hits += 1;
        }
    }
    let shared = t0.elapsed();

    let t1 = Instant::now();
    let mut hits2 = 0u64;
    for i in 0..probes {
        let key = keys[(i as usize) & 4095];
        pair.prefetch_key(key);
        if pair.contains_key(key) {
            hits2 += 1;
        }
    }
    let doubled = t1.elapsed();

    assert_eq!(hits, hits2, "probe answers must agree");
    let ns_shared = shared.as_secs_f64() * 1e9 / probes as f64;
    let ns_doubled = doubled.as_secs_f64() * 1e9 / probes as f64;
    println!(
        "n={n} |F|={} probes={probes} hits={hits} (miss-dominated)",
        fb.points.len()
    );
    println!("shared hash (prefetch+reject): {ns_shared:.2} ns/probe");
    println!("doubled hash (pre cycle-1):    {ns_doubled:.2} ns/probe");
    println!(
        "speedup (doubled/shared):       {:.3}×",
        ns_doubled / ns_shared
    );
}
