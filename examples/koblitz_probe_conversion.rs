//! The conversion factor a tier ladder needs to be priced in one unit.
//!
//! AGENTS.md §2 asks for one unit, and §6 rules out wall-clock as the
//! headline.  The pair-table thread has been reporting seconds, which
//! makes its tiers incomparable with everything else on the scoreboard:
//! the ledger's unit is `S = group-addition equivalents / sqrt(r)`, and
//! a pair-table probe is not a group addition.
//!
//! So measure the conversion, on this host, in this process, on the
//! base the ladder actually runs: how many group additions a single
//! table probe costs, per tier.  Everything else in the pipeline is
//! already counted natively — the build does one addition per stored
//! pair, rho does one per step — so this one factor is what stands
//! between the tier ladder and the ledger's axis.
//!
//! The factor is per tier on purpose.  A probe is a hash and a short
//! scan either way; what differs is whether the table it reads is in
//! cache, and at 12,688 points the folded table is 16 MiB against the
//! full table's 1,592.  That difference is a property of this host and
//! is recorded as such, not smuggled into an operation count.
//!
//!     cargo run --release --example koblitz_probe_conversion -- 61 12000

use crypto_lib::cryptanalysis::koblitz_fast::{BatchScratch, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::*;
use std::time::Instant;

/// Enough probes that the timer's resolution is not the measurement.
const PROBES: usize = 200_000;
/// Group additions timed against the same clock.
const ADDS: usize = 2_000_000;

fn main() {
    let degree: u32 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(61);
    let points: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(12000);
    let kc = KoblitzCurve::new(0, degree).expect("curve");
    let fb = build_subgroup_orbit_factor_base(&kc, 1, points).expect("factor base");
    let fc = FastCurve::new(&kc.curve).expect("fast curve");
    let n = fb.points.len();
    let orbits = fb.signed_orbits.len();
    let r = kc.subgroup_order.to_u64_digits()[0];

    // The unit: one affine group addition on this curve, on this host.
    // Measured twice, because the pipeline performs both kinds.
    //
    // A *chained* addition pays its own field inversion; a *batched*
    // one shares a single inversion across the slice by Montgomery's
    // trick, which is what the rho walk does across its parallel walks
    // and what the pair-table build does down a row.  The batched
    // addition is therefore the honest unit here: it is the cheaper of
    // the two, so quoting probes in batched additions charges the
    // probe the most, which is the direction an accounting correction
    // should err in.
    let g = fc.lift(kc.generator());
    let mut acc = fc.lift(&fb.points[0]);
    let t = Instant::now();
    for _ in 0..ADDS {
        acc = fc.add(acc, g);
    }
    let chained_ns = t.elapsed().as_secs_f64() * 1e9 / ADDS as f64;
    std::hint::black_box(acc);

    let batch: Vec<FastPoint> = (0..1024).map(|i| fc.lift(&fb.points[i % n])).collect();
    let mut scratch = BatchScratch::default();
    let mut out: Vec<FastPoint> = Vec::with_capacity(batch.len());
    let rounds = ADDS / batch.len();
    let t = Instant::now();
    for _ in 0..rounds {
        out.clear();
        fc.add_many(g, &batch, &mut out, &mut scratch);
    }
    let add_ns = t.elapsed().as_secs_f64() * 1e9 / (rounds * batch.len()) as f64;
    std::hint::black_box(out.len());

    println!("base {n} points, {orbits} signed orbits, degree {degree}, r = 2^{:.1}", (r as f64).log2());
    println!("one group addition, chained (own inversion):  {chained_ns:.2} ns");
    println!("one group addition, batched (shared inverse): {add_ns:.2} ns   <- the unit");

    // The same probes against each representation of the same base, so
    // the only thing that differs between the rows is the table.
    let mut rng: u64 = 0x9e37_79b9_7f4a_7c15;
    let mut targets = Vec::with_capacity(PROBES);
    for _ in 0..PROBES {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        targets.push(fc.mul_u64(g, 1 + rng % (r - 1)).pack());
    }

    let tiers: [(&str, Option<PairSumTable>); 3] = [
        ("folded", PairSumTable::build_folded_within(&kc, &fb, PairSumTable::DEFAULT_BYTE_BUDGET)),
        ("compact", PairSumTable::build_compact_within(&kc, &fb, PairSumTable::DEFAULT_BYTE_BUDGET)),
        ("full", PairSumTable::build_full_within(&kc, &fb, PairSumTable::DEFAULT_BYTE_BUDGET)),
    ];
    println!("\n{:<8} {:>14} {:>12} {:>12} {:>16}", "tier", "stored pairs", "build s", "probe ns", "adds per probe");
    for (name, table) in tiers {
        let Some(table) = table else {
            println!("{name:<8} {:>14}", "did not fit");
            continue;
        };
        // Rebuild under the clock: the build is a counted phase too, and
        // its native count is one addition per stored pair.
        let t = Instant::now();
        let rebuilt = match name {
            "folded" => PairSumTable::build_folded_within(&kc, &fb, PairSumTable::DEFAULT_BYTE_BUDGET),
            "compact" => PairSumTable::build_compact_within(&kc, &fb, PairSumTable::DEFAULT_BYTE_BUDGET),
            _ => PairSumTable::build_full_within(&kc, &fb, PairSumTable::DEFAULT_BYTE_BUDGET),
        }
        .expect("it fit a moment ago");
        let build_s = t.elapsed().as_secs_f64();
        drop(rebuilt);

        let t = Instant::now();
        let mut hits = 0usize;
        for &target in &targets {
            if table.contains_pair(target) {
                hits += 1;
            }
        }
        let probe_ns = t.elapsed().as_secs_f64() * 1e9 / PROBES as f64;
        std::hint::black_box(hits);
        println!(
            "{name:<8} {:>14} {build_s:>12.2} {probe_ns:>12.1} {:>16.2}",
            table.len(),
            probe_ns / add_ns
        );
    }
    println!(
        "\nrecord both numbers with the run: the adds-per-probe factor is a\n\
         property of this host and this base, and a ledger that quotes S\n\
         without it has converted a foreign unit by assumption."
    );
}
