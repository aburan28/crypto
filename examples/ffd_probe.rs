//! First-fall-degree probe for the Crossbred AutoLab X5 beat.
//!
//! ```bash
//! cargo run --release --example ffd_probe -- 9:3 9:4 15:3 15:4
//! FFD_TRIALS=4 cargo run --release --example ffd_probe -- 9:4 15:4
//! ```
//!
//! Uses `ffd_summary` with factor index 0, `d_max = 3`, seed `0x5EED`.
//! That is the scaling-target convention, not a per-rung choice.

use crypto_lib::cryptanalysis::koblitz_bench::{ffd_summary, format_ffd_table};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let cases: Vec<(u32, usize)> = if args.is_empty() {
        vec![(9, 3), (9, 4), (15, 3), (15, 4)]
    } else {
        args.iter()
            .filter_map(|a| {
                let (n, m) = a.split_once(':')?;
                Some((n.parse().ok()?, m.parse().ok()?))
            })
            .collect()
    };
    let trials: usize = std::env::var("FFD_TRIALS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4);
    let d_max = 3u32;
    let seed = 0x5EED_u64;
    let mut rows = Vec::new();
    for (n, m) in cases {
        match ffd_summary(n, 0, m, d_max, seed, trials) {
            Some(s) => rows.push(s),
            None => eprintln!("skip n={n} m={m}: no system under 64 vars"),
        }
    }
    print!("{}", format_ffd_table(&rows));
}
