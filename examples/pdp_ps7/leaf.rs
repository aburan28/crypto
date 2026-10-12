//! The leaf: given the first two summands `(a, b)` of a decomposition, the
//! third abscissae `t ∈ V` with `S₄(a, b, t, x_R) = 0`.
//!
//! * [`GroupLeaf`] is PS-2's leaf (`examples/pdp_speedup/kernel.rs`): the four
//!   group-law candidates `x(P₁ ± P₂ ± R)`, a hit when a candidate is below
//!   `2^ℓ`.  Called for one pair at a time (a CDCL theory sees one pair), so
//!   the two denominators of the pair share **one inversion** (Montgomery's
//!   trick on two elements) instead of the row-batched inversion of the plain
//!   loop; the row points `U± = P₁ ± R` are cached per `x₁`.
//! * [`QuarticLeaf`] is PS-0's independent root enumerator (`corpus::pair_roots`:
//!   the quartic in `x₃`, `gcd` with `L_V`): the algebraic question over all of
//!   `V`, used where the group law does not apply (non-point targets).
//!
//! Every field operation of a [`GroupLeaf`] is counted; [`Gates`] prices the
//! counts in PS-1's gate model so that no arm gets the leaf for free.

use crate::common::Rung;
use crate::corpus::pair_roots;
use crate::group::{add, neg, thirds, Pt};
use crypto_lib::cryptanalysis::semaev_decomp::{subspace_poly, Gf2};

/// Field-operation counts of leaf work (not of shell work).
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct LeafCounts {
    pub mul: u64,
    pub sqr: u64,
    pub inv: u64,
    /// Additions of two `n`-bit elements (`n` gates each).
    pub add_n: u64,
    /// Flips of one bit (the constant `+1`), one gate each.
    pub flip1: u64,
    /// Zero tests of an `n`-bit element (denominators).
    pub zero_tests: u64,
    /// Membership tests `t ∈ V` (the top `n − ℓ` bits are zero).
    pub member_tests: u64,
    /// Pairs sent to the exact reference (zero denominator); priced as four exact additions.
    pub exact_fallbacks: u64,
    /// Leaf evaluations (pairs).
    pub evals: u64,
}

impl LeafCounts {
    pub fn plus(&self, o: &LeafCounts) -> LeafCounts {
        LeafCounts {
            mul: self.mul + o.mul,
            sqr: self.sqr + o.sqr,
            inv: self.inv + o.inv,
            add_n: self.add_n + o.add_n,
            flip1: self.flip1 + o.flip1,
            zero_tests: self.zero_tests + o.zero_tests,
            member_tests: self.member_tests + o.member_tests,
            exact_fallbacks: self.exact_fallbacks + o.exact_fallbacks,
            evals: self.evals + o.evals,
        }
    }

    /// Cost of one generic exact group addition (`group::add`, `x₁ ≠ x₂`):
    /// one inversion, two products, one squaring, seven `n`-bit additions, one bit flip.
    pub fn exact_add() -> LeafCounts {
        LeafCounts {
            mul: 2,
            sqr: 1,
            inv: 1,
            add_n: 7,
            flip1: 1,
            ..Default::default()
        }
    }

    pub fn scaled(&self, k: u64) -> LeafCounts {
        LeafCounts {
            mul: self.mul * k,
            sqr: self.sqr * k,
            inv: self.inv * k,
            add_n: self.add_n * k,
            flip1: self.flip1 * k,
            zero_tests: self.zero_tests * k,
            member_tests: self.member_tests * k,
            exact_fallbacks: self.exact_fallbacks * k,
            evals: self.evals * k,
        }
    }
}

/// PS-1's gate counts (fused ternary-logic issues, `ps1_algebra/raw/gates.jsonl`).
#[derive(Clone, Copy, Debug)]
pub struct Gates {
    pub n: u64,
    pub mul: u64,
    pub sqr: u64,
    pub inv: u64,
    pub zero_mask: u64,
    pub in_subspace_mask: u64,
}

impl Gates {
    /// Read the row for `(n, l)` from PS-1's `gates.jsonl`.
    pub fn from_ps1(text: &str, n: u32, l: u32) -> Option<Gates> {
        for line in text.lines() {
            let v: serde_json::Value = serde_json::from_str(line).ok()?;
            if v["n"].as_u64()? == n as u64 && v["l"].as_u64()? == l as u64 {
                return Some(Gates {
                    n: n as u64,
                    mul: v["mul"].as_u64()?,
                    sqr: v["sqr"].as_u64()?,
                    inv: v["inv"].as_u64()?,
                    zero_mask: v["zero_mask"].as_u64()?,
                    in_subspace_mask: v["in_subspace_mask"].as_u64()?,
                });
            }
        }
        None
    }

    /// Counted bit operations (gates) of the leaf work in `c`.
    pub fn price(&self, c: &LeafCounts) -> u64 {
        c.mul * self.mul
            + c.sqr * self.sqr
            + c.inv * self.inv
            + c.add_n * self.n
            + c.flip1
            + c.zero_tests * self.zero_mask
            + c.member_tests * self.in_subspace_mask
            + c.exact_fallbacks * 4 * self.price_exact_add()
    }

    fn price_exact_add(&self) -> u64 {
        let e = LeafCounts::exact_add();
        e.mul * self.mul + e.sqr * self.sqr + e.inv * self.inv + e.add_n * self.n + e.flip1
    }
}

pub trait Leaf {
    /// Sorted candidate abscissae `S`; the pair space is `{(a, b) : a ≤ b, a, b ∈ S}`.
    fn space(&self) -> &[u64];
    /// Third abscissae `t` (all of them, no duplicates) for the pair `(a, b)`, `a ≤ b`, both in `S`.
    fn eval(&mut self, a: u64, b: u64, out: &mut Vec<u64>);
    fn counts(&self) -> LeafCounts;
    fn name(&self) -> &'static str;
}

// ── PS-2's leaf, one pair at a time ────────────────────────────────────

#[derive(Clone, Copy)]
struct RowPts {
    xu: u64,
    yu: u64,
    xv: u64,
    yv: u64,
}

/// Rational base points of one rung, shared by every target of the rung.
pub struct BaseTable {
    pub base: Vec<(u64, u64)>,
    pub xs: Vec<u64>,
    /// `y` of each base abscissa, indexed by `x`; `u64::MAX` for a non-liftable `x`.
    pub ytab: Vec<u64>,
}

impl BaseTable {
    pub fn new(base: Vec<(u64, u64)>, l: u32) -> Self {
        let mut ytab = vec![u64::MAX; 1usize << l];
        for &(x, y) in &base {
            ytab[x as usize] = y;
        }
        let xs = base.iter().map(|p| p.0).collect();
        Self { base, xs, ytab }
    }
}

pub struct GroupLeaf<'a> {
    gf: &'a Gf2,
    l: u32,
    tab: &'a BaseTable,
    r: (u64, u64),
    rows: Vec<Option<RowPts>>,
    cnt: LeafCounts,
}

impl<'a> GroupLeaf<'a> {
    pub fn new(gf: &'a Gf2, l: u32, tab: &'a BaseTable, r: (u64, u64)) -> Self {
        Self {
            gf,
            l,
            tab,
            r,
            rows: vec![None; 1usize << l],
            cnt: LeafCounts::default(),
        }
    }

    fn row(&mut self, a: u64) -> Option<RowPts> {
        if let Some(r) = self.rows[a as usize] {
            return Some(r);
        }
        let p1 = Pt::Aff(a, self.tab.ytab[a as usize]);
        let rp = Pt::Aff(self.r.0, self.r.1);
        let (Pt::Aff(xu, yu), Pt::Aff(xv, yv)) =
            (add(self.gf, p1, rp), add(self.gf, p1, neg(rp)))
        else {
            return None;
        };
        // Two exact additions of the generic case.
        self.cnt = self.cnt.plus(&LeafCounts::exact_add().scaled(2));
        let r = RowPts { xu, yu, xv, yv };
        self.rows[a as usize] = Some(r);
        Some(r)
    }
}

impl Leaf for GroupLeaf<'_> {
    fn space(&self) -> &[u64] {
        &self.tab.xs
    }

    fn eval(&mut self, a: u64, b: u64, out: &mut Vec<u64>) {
        out.clear();
        let gf = self.gf;
        let span = 1u64 << self.l;
        self.cnt.evals += 1;
        let (x2, y2) = (b, self.tab.ytab[b as usize]);
        let Some(RowPts { xu, yu, xv, yv }) = self.row(a) else {
            // P₁ = ±R: only when x_R ∈ V, which the caller excludes.
            self.cnt.exact_fallbacks += 1;
            for t in thirds(gf, (a, self.tab.ytab[a as usize]), (x2, y2), self.r)
                .into_iter()
                .flatten()
            {
                if t < span && !out.contains(&t) {
                    out.push(t);
                }
            }
            return;
        };
        let d = xu ^ x2;
        let dd = xv ^ x2;
        self.cnt.add_n += 2;
        self.cnt.zero_tests += 2;
        if d == 0 || dd == 0 {
            self.cnt.exact_fallbacks += 1;
            for t in thirds(gf, (a, self.tab.ytab[a as usize]), (x2, y2), self.r)
                .into_iter()
                .flatten()
            {
                if t < span && !out.contains(&t) {
                    out.push(t);
                }
            }
            return;
        }
        // One inversion for both denominators: inv(d·dd), then 1/d and 1/dd by one product each.
        let inv = gf.inv(gf.mul(d, dd));
        let (j, jj) = (gf.mul(inv, dd), gf.mul(inv, d));
        self.cnt.mul += 3;
        self.cnt.inv += 1;
        let lam = gf.mul(yu ^ y2, j);
        let t1 = gf.sqr(lam) ^ lam ^ d ^ 1;
        let lam2 = lam ^ gf.mul(x2, j);
        let t2 = gf.sqr(lam2) ^ lam2 ^ d ^ 1;
        let mu = gf.mul(yv ^ y2, jj);
        let t3 = gf.sqr(mu) ^ mu ^ dd ^ 1;
        let mu2 = mu ^ gf.mul(x2, jj);
        let t4 = gf.sqr(mu2) ^ mu2 ^ dd ^ 1;
        self.cnt.mul += 4;
        self.cnt.sqr += 4;
        // yu^y2, yv^y2; per candidate chain: sqr^lam^d (2 adds); lam2 and mu2 (1 add each).
        self.cnt.add_n += 2 + 4 * 2 + 2;
        self.cnt.flip1 += 4;
        self.cnt.member_tests += 4;
        for t in [t1, t2, t3, t4] {
            if t < span && !out.contains(&t) {
                out.push(t);
            }
        }
    }

    fn counts(&self) -> LeafCounts {
        self.cnt
    }

    fn name(&self) -> &'static str {
        "group-law (PS-2), one pair at a time"
    }
}

/// Field operations that PS-2's row-batched plain loop does for a run that
/// entered `rows` rows of a base of `m` points and evaluated `pairs` pairs
/// (derived from the loop's structure in `kernel.rs`, not instrumented):
/// per row two exact additions (`U±`) and one `batch_inv` of `2·cnt` elements
/// (`3(k − 1)` products and one inversion for `k` elements); per pair four
/// products, four squarings and the add/test work of the candidate chains.
pub fn plain_loop_counts(m: u64, rows: u64, pairs: u64) -> LeafCounts {
    let mut c = LeafCounts::default();
    for i in 0..rows {
        let k = 2 * (m - i);
        c = c.plus(&LeafCounts::exact_add().scaled(2));
        c.mul += 3 * (k - 1);
        c.inv += 1;
    }
    c.mul += 4 * pairs;
    c.sqr += 4 * pairs;
    // d, dd; yu^y2, yv^y2; four chains of two adds; lam2, mu2.
    c.add_n += (2 + 2 + 8 + 2) * pairs;
    c.flip1 += 4 * pairs;
    c.zero_tests += 2 * pairs;
    c.member_tests += 4 * pairs;
    c.evals = pairs;
    c
}

/// Rows entered by the plain loop after `pairs` pairs of a base of `m` points
/// (a row is entered, and its batched inversion paid, when its first pair is visited).
pub fn rows_entered(m: u64, pairs: u64) -> u64 {
    let (mut rows, mut done) = (0u64, 0u64);
    while done < pairs {
        done += m - rows;
        rows += 1;
    }
    rows
}

// ── PS-0's algebraic leaf (quartic and gcd), all of V ──────────────────

pub struct QuartInfo {
    pub lv: Vec<u64>,
}

pub struct QuarticLeaf<'a> {
    rung: &'a Rung,
    xr: u64,
    lv: Vec<u64>,
    space: Vec<u64>,
    cnt: LeafCounts,
    scratch: Vec<[u64; 3]>,
}

impl<'a> QuarticLeaf<'a> {
    pub fn new(rung: &'a Rung, xr: u64) -> Self {
        Self {
            rung,
            xr,
            lv: subspace_poly(rung.l, &rung.gf),
            space: (0..rung.span()).collect(),
            cnt: LeafCounts::default(),
            scratch: Vec::new(),
        }
    }
}

impl Leaf for QuarticLeaf<'_> {
    fn space(&self) -> &[u64] {
        &self.space
    }

    fn eval(&mut self, a: u64, b: u64, out: &mut Vec<u64>) {
        out.clear();
        self.cnt.evals += 1;
        self.scratch.clear();
        pair_roots(self.rung, &self.lv, self.xr, a, b, &mut self.scratch);
        for tr in &self.scratch {
            // The triple is sorted(a, b, t): remove one a and one b, the rest is t.
            let mut rest = tr.to_vec();
            let ia = rest.iter().position(|&v| v == a).expect("a in triple");
            rest.remove(ia);
            let ib = rest.iter().position(|&v| v == b).expect("b in triple");
            rest.remove(ib);
            let t = rest[0];
            if !out.contains(&t) {
                out.push(t);
            }
        }
    }

    fn counts(&self) -> LeafCounts {
        self.cnt
    }

    fn name(&self) -> &'static str {
        "quartic and gcd with L_V (PS-0 enumerator)"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rows_entered_matches_a_direct_walk() {
        let m = 7u64;
        let total = m * (m + 1) / 2;
        let mut walk = Vec::new();
        for i in 0..m {
            for _ in i..m {
                walk.push(i);
            }
        }
        for p in 1..=total {
            assert_eq!(rows_entered(m, p), walk[(p - 1) as usize] + 1, "pairs {p}");
        }
    }
}
