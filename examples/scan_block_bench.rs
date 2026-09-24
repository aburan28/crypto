//! Times the two halves of the m = 3 decomposition scan on one block:
//! the batched subtraction `R − P_k` (`FastCurve::add_many_lazy`) and the
//! Frobenius-orbit key of every rest (`FrobeniusCanon::canon`).
//!
//! `cargo run --release --example scan_block_bench`
use crypto_lib::cryptanalysis::koblitz_fast::{BatchScratch, FastCurve, FastPoint, FrobeniusCanon};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use std::hint::black_box;
use std::time::Instant;

fn main() {
    let kc = KoblitzCurve::new(0, 53).expect("k0n53");
    let fc = FastCurve::new(&kc.curve).unwrap();
    let canon = FrobeniusCanon::new(&fc.field, 53).unwrap();
    let g = fc.lift(kc.generator());
    let mut pts: Vec<FastPoint> = Vec::new();
    let mut p = g;
    let step = fc.mul(g, &BigUint::from(0x9e37_79b9_7f4au64));
    for _ in 0..1024 {
        pts.push(fc.neg(p));
        p = fc.add(p, step);
    }
    let target = fc.mul(g, &BigUint::from(123_456_789u64));
    let (mut out, mut lambdas, mut scratch) = (Vec::new(), Vec::new(), BatchScratch::default());
    let iters = 2000;
    let t = Instant::now();
    let mut chk = 0u64;
    for _ in 0..iters {
        out.clear();
        lambdas.clear();
        fc.add_many_lazy(
            black_box(target),
            &pts,
            &mut out,
            &mut lambdas,
            &mut scratch,
        );
        chk ^= out[17].x;
    }
    let sub_ns = t.elapsed().as_nanos() as f64 / (iters * 1024) as f64;
    let t = Instant::now();
    for _ in 0..iters {
        for q in &out {
            chk ^= canon.canon(black_box(q.x));
        }
    }
    let key_ns = t.elapsed().as_nanos() as f64 / (iters * 1024) as f64;
    println!("per point: subtraction {sub_ns:.2} ns, key {key_ns:.2} ns   (chk {chk:x})");
}
