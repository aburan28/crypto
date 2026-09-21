//! **What the `m = 4` arm actually pays per `(k, l)`.**
//!
//! `RESEARCH_KOBLITZ_INDEX_CALCULUS.md` names the `m = 4` arm as the
//! one place that still probes a folded table one target at a time,
//! and prices that at "the 265 rather than the 145".  That is the
//! smaller half of what the loop does one at a time.  Its inner step
//! is
//!
//! ```ignore
//! let rest = self.curve.add(target, self.curve.neg(pair));
//! self.pairs_for(rest, &mut pairs);
//! ```
//!
//! and `FastCurve::add` calls `Gf2::inv`, a Fermat inversion — `n − 1`
//! squarings and as many multiplications — where the `m = 3` arm's
//! `add_many` amortises one inversion over a whole block by
//! Montgomery's trick.  This prices the two against each other.
use std::time::Instant;

use crypto_lib::cryptanalysis::koblitz_fast::{BatchScratch, FastCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::*;

fn main() {
    let degree: u32 = 61;
    let kc = KoblitzCurve::new(0, degree).expect("curve");
    let fc = FastCurve::new(&kc.curve).expect("fast curve");
    let g = fc.lift(kc.generator());
    let n = 200_000usize;

    let pts: Vec<FastPoint> = (1u64..=n as u64)
        .map(|t| fc.mul_u64(g, t * 1_000_003 + 5))
        .filter(|p| !p.infinity)
        .collect();
    let target = fc.mul_u64(g, 7_700_017);
    let xs: Vec<u64> = pts.iter().map(|p| p.x).collect();

    // The Fermat inversion on its own.
    let t0 = Instant::now();
    let mut acc = 0u64;
    for &x in &xs {
        acc ^= fc.field.inv(x);
    }
    let inv_ns = t0.elapsed().as_secs_f64() * 1e9 / xs.len() as f64;
    std::hint::black_box(acc);

    // What the `m = 4` inner step does: one `add` per pair, which is
    // one inversion per pair.
    let negs: Vec<FastPoint> = pts.iter().map(|&p| fc.neg(p)).collect();
    let t0 = Instant::now();
    let mut sink = 0u64;
    for &q in &negs {
        sink ^= fc.add(target, q).x;
    }
    let add_ns = t0.elapsed().as_secs_f64() * 1e9 / negs.len() as f64;
    std::hint::black_box(sink);

    // What the `m = 3` arm does instead: one inversion a block.
    let mut scratch = BatchScratch::default();
    let mut out = Vec::new();
    let t0 = Instant::now();
    for block in negs.chunks(1024) {
        out.clear();
        fc.add_many(target, block, &mut out, &mut scratch);
        sink ^= out[0].x;
    }
    let add_many_ns = t0.elapsed().as_secs_f64() * 1e9 / negs.len() as f64;
    std::hint::black_box(sink);

    println!("at n = {degree}, {} points:", negs.len());
    println!("  Gf2::inv alone                        {inv_ns:8.1} ns");
    println!("  FastCurve::add      (m = 4's step)    {add_ns:8.1} ns  <- one inversion each");
    println!("  add_many, 1024 wide (m = 3's step)    {add_many_ns:8.1} ns");
    println!();
    println!(
        "  the unbatched inversion costs {:.1} ns a pair, {:.0}x the batched step",
        add_ns - add_many_ns,
        add_ns / add_many_ns
    );
    println!(
        "  for comparison, the lone-vs-blocked PROBE difference is about 104 ns,\n  \
         so the inversion is {:.1}x the cost the note named",
        (add_ns - add_many_ns) / 104.0
    );
}
