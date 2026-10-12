//! The Koblitz curve `K_0 : y² + xy = x³ + 1` (`a = 0`, `b = 1`) over a field
//! with the [`Fld`] interface.
//!
//! Group law (`a_1 = 1`, `a_2 = a = 0`, `a_6 = b = 1`):
//!
//! ```text
//!   −(x, y)    = (x, y + x)
//!   P₁ ≠ ±P₂:  λ = (y₁ + y₂)/(x₁ + x₂),  x₃ = λ² + λ + x₁ + x₂ + a,  y₃ = λ(x₁ + x₃) + x₃ + y₁
//!   2P:        λ = x + y/x,              x₃ = λ² + λ + a,            y₃ = x² + (λ + 1)x₃
//! ```
//!
//! with `a = 0`: there is no `+1` in `x₃`, unlike the `a = b = 1` curve of the
//! PS-1 formulas.  The only 2-torsion point is `T = (0, 1)` and translation by
//! it is `x ↦ 1/x` (checked in `checks.rs`).

use crate::field::Fld;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Pt {
    Inf,
    Aff(u64, u64),
}

pub fn neg(p: Pt) -> Pt {
    match p {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => Pt::Aff(x, y ^ x),
    }
}

/// Exact group law, every special case, one inversion per call.
pub fn add<F: Fld>(f: &F, p: Pt, q: Pt) -> Pt {
    match (p, q) {
        (Pt::Inf, r) | (r, Pt::Inf) => r,
        (Pt::Aff(x1, y1), Pt::Aff(x2, y2)) => {
            if x1 == x2 {
                if y1 ^ y2 == x1 {
                    return Pt::Inf;
                }
                // y1 == y2 and x1 != 0: doubling.
                let lam = x1 ^ f.mul(y1, f.inv(x1));
                let x3 = f.sqr(lam) ^ lam;
                let y3 = f.sqr(x1) ^ f.mul(lam ^ f.one(), x3);
                return Pt::Aff(x3, y3);
            }
            let lam = f.mul(y1 ^ y2, f.inv(x1 ^ x2));
            let x3 = f.sqr(lam) ^ lam ^ x1 ^ x2;
            let y3 = f.mul(lam, x1 ^ x3) ^ x3 ^ y1;
            Pt::Aff(x3, y3)
        }
    }
}

pub fn on_curve<F: Fld>(f: &F, x: u64, y: u64) -> bool {
    let lhs = f.sqr(y) ^ f.mul(x, y);
    let rhs = f.mul(f.sqr(x), x) ^ f.one();
    lhs == rhs
}

/// `Σ_{i=0}^{(n−1)/2} c^{4^i}`, a solution of `z² + z = c` when `Tr(c) = 0` (`n` odd).
pub fn half_trace<F: Fld>(f: &F, c: u64) -> u64 {
    let mut out = 0u64;
    let mut term = c;
    for _ in 0..=f.n() / 2 {
        out ^= term;
        term = f.sqr(f.sqr(term));
    }
    out
}

/// A rational `y` over `x`, the library's first root (`y = x·z`), or `None`.
pub fn lift<F: Fld>(f: &F, x: u64) -> Option<u64> {
    if x == 0 {
        return Some(f.one());
    }
    // y = xz with z² + z = x + 1/x²  (from y² + xy = x³ + 1).
    let c = x ^ f.sqr(f.inv(x));
    if c.count_ones() & 1 == 1 {
        return None;
    }
    let z = half_trace(f, c);
    Some(f.mul(x, z))
}

/// The four candidate third abscissae of the pair `(P₁, P₂)` for the target
/// `R`: `x(R + P₁ + P₂), x(R + P₁ − P₂), x(R − P₁ + P₂), x(R − P₁ − P₂)`
/// (`None` where that point is the point at infinity).  Exact, slow.
pub fn thirds<F: Fld>(f: &F, p1: (u64, u64), p2: (u64, u64), r: (u64, u64)) -> [Option<u64>; 4] {
    let (p1, p2) = (Pt::Aff(p1.0, p1.1), Pt::Aff(p2.0, p2.1));
    let r = Pt::Aff(r.0, r.1);
    let mut out = [None; 4];
    for (k, z) in [add(f, p1, p2), add(f, p1, neg(p2))]
        .into_iter()
        .enumerate()
    {
        for (m, rr) in [r, neg(r)].into_iter().enumerate() {
            if let Pt::Aff(x, _) = add(f, z, rr) {
                out[2 * k + m] = Some(x);
            }
        }
    }
    out
}

pub fn sorted3(a: u64, b: u64, c: u64) -> [u64; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}
