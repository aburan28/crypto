//! A conditional cost model: what a bounded last fall degree would mean for
//! index calculus on binary curves, against rho.
//!
//! **Model, not a measurement.** It assumes the measured toy-size value
//! `d_F = 4` holds at every size (which the HKYY conditional in
//! `research/notes/index-calculus/RESEARCH_FALL_DEGREE_BOUNDS.md` §2.2 says
//! cannot be true for growing `m` unless ECDLP over `F_{2^n}` is
//! subexponential), and that a decomposition is decided and solved by linear
//! algebra in `V_{F,4}` (HKY Prop 2.8) at dense cost `M(4)^ω`,
//! `M(4) = Σ_{k≤4} C(N, k)`, with `N = (m−2)n + mℓ` unknowns (chained `S₃`).
//!
//! For each `n` and `m`: `ℓ = ⌈(n + log₂ m!)/m⌉` so a random target
//! decomposes with probability about 1; relations needed `≈ 2^ℓ` (factor-base
//! size; Frobenius folding would divide by `n`, a constant ignored here); total
//! `log₂(2^ℓ · M(4)^ω + 2^{2ℓ})` (relation search + sparse linear algebra);
//! rho `n/2` bits. Prints the best `m` per `n` and the rho margin.
//!
//! Usage: `fall_degree_crossover [omega]` (default 2.81).
fn log2_binom_sum(n: f64, d: u32) -> f64 {
    let mut s = 0.0f64;
    let mut c = 1.0f64;
    for k in 0..=d {
        if k > 0 {
            c *= (n - (k as f64 - 1.0)) / k as f64;
        }
        s += c;
    }
    s.log2()
}

fn log2_add(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (1.0 + (lo - hi).exp2()).log2()
}

fn main() {
    let omega: f64 = std::env::args().nth(1).map_or(2.81, |s| s.parse().unwrap());
    println!("| n | rho (bits) | best m | ℓ | N | log₂ M(4) | log₂ per-target | log₂ total (model) | margin vs rho |");
    println!("|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    for n in [131u32, 163, 233, 283, 409, 571, 800, 1000] {
        let nf = n as f64;
        let mut best: Option<(f64, u32, u32, f64, f64)> = None;
        for m in 3..=40u32 {
            let lf: f64 = (2..=m).map(|k| (k as f64).log2()).sum();
            let ell = ((nf + lf) / m as f64).ceil() as u32;
            let nv = (m as f64 - 2.0) * nf + m as f64 * ell as f64;
            let lm = log2_binom_sum(nv, 4);
            let per = omega * lm;
            let total = log2_add(ell as f64 + per, 2.0 * ell as f64);
            if best.is_none_or(|b| total < b.0) {
                best = Some((total, m, ell, nv, lm));
            }
        }
        let (total, m, ell, nv, lm) = best.unwrap();
        println!(
            "| {n} | {:.1} | {m} | {ell} | {nv:.0} | {lm:.1} | {:.1} | {total:.1} | {:+.1} |",
            nf / 2.0,
            omega * lm,
            total - nf / 2.0
        );
    }
}
