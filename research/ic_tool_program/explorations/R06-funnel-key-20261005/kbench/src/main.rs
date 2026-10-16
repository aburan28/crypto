//! How often does a rotation-invariant hash of a random n-bit word match the
//! invariant of some stored orbit while the orbit itself is not stored?
//! Stored: N random words (their orbits). Probes: M random words.
use std::collections::HashSet;

fn rotl(c: u64, t: u32, n: u32, mask: u64) -> u64 {
    let t = t % n;
    if t == 0 { c } else { ((c << t) | (c >> (n - t))) & mask }
}

fn least_rotation(c: u64, n: u32, mask: u64) -> u64 {
    let mut best = c;
    for t in 1..n {
        let v = rotl(c, t, n, mask);
        if v < best { best = v; }
    }
    best
}

fn inv(c: u64, n: u32, mask: u64, ks: &[u32]) -> u64 {
    let mut acc = c.count_ones() as u64;
    for &k in ks {
        let a = (c & rotl(c, k, n, mask)).count_ones() as u64;
        acc = acc.rotate_left(6) ^ a;
    }
    acc
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: u32 = args[1].parse().unwrap();
    let stored: usize = args[2].parse().unwrap();
    let probes: usize = args[3].parse().unwrap();
    let mask = (1u64 << n) - 1;
    let variants: Vec<(&str, Vec<u32>)> = vec![
        ("k1-4", (1..=4).collect()),
        ("k1-6", (1..=6).collect()),
        ("k1-8", (1..=8).collect()),
        ("k1-10", (1..=10).collect()),
        ("spread8", vec![1, 2, 3, 5, 8, 13, 21, 29]),
    ];
    let mut s = 0x9E37_79B9_7F4A_7C15u64 ^ (n as u64);
    let mut next = move || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; s };
    let mut canon = HashSet::with_capacity(stored);
    let mut words = Vec::with_capacity(stored);
    while canon.len() < stored {
        let c = next() & mask;
        if canon.insert(least_rotation(c, n, mask)) { words.push(c); }
    }
    let probe_words: Vec<u64> = (0..probes).map(|_| next() & mask).collect();
    for (name, ks) in &variants {
        let set: HashSet<u64> = words.iter().map(|&c| inv(c, n, mask, ks)).collect();
        let mut false_match = 0usize;
        let mut true_match = 0usize;
        for &c in &probe_words {
            if set.contains(&inv(c, n, mask, ks)) {
                if canon.contains(&least_rotation(c, n, mask)) { true_match += 1; } else { false_match += 1; }
            }
        }
        println!("n={n} stored={stored} distinct_inv={} probes={probes} {name}: false-match rate {:.3e} (true {true_match})",
            set.len(), false_match as f64 / probes as f64);
    }
}
