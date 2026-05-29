//! **Elliptic divisibility sequences, elliptic nets, and the EDS-Residue
//! handle on the ECDLP** (Shipsey 2000, Stange 2007, Lauter–Stange 2008).
//!
//! This module is the runnable companion to `RESEARCH_EDS_RESIDUE.md`.  It
//! builds the elliptic divisibility sequence (EDS) `W_{E,P}(n)` attached to
//! a point `P` on a short-Weierstrass curve over `F_p`, and measures the
//! structures that the Lauter–Stange equivalence chain hangs on:
//!
//! 1. **Rank of apparition.**  `W(n) = 0  ⟺  [n]P = O`, so the first zero
//!    of the sequence is exactly `ord(P)`.  This is the anchor that ties
//!    the sequence to the group, and we cross-check it against honest
//!    point arithmetic.
//! 2. **Periodicity.**  Over `F_p` the sequence `W(n) mod p` is purely
//!    periodic; we measure its period `π_W`.
//! 3. **The quadratic-residuosity / Legendre sequence** `χ(n) = (W(n) | p)`.
//!    The *EDS-Residue problem* (Lauter–Stange) is: given the curve, `P`,
//!    and `Q = [k]P` (but **not** `k`), decide whether `W_{E,P}(k)` is a
//!    quadratic residue.  We measure the period `π_χ`, the QR / NQR / zero
//!    balance, and whether the χ-pattern localises `k`.
//! 4. **The zero-lattice of the 2-D elliptic net.**  For the pair `(P, Q)`
//!    the net value `W(a, b)` vanishes iff `[a]P + [b]Q = O`.  When
//!    `Q = [k]P` the zero set is the sublattice `{(a,b) : a + bk ≡ 0
//!    (mod m)}`, whose slope *is* the discrete log.  We exhibit that
//!    lattice and recover `k` from it.
//!
//! ## The recurrence
//!
//! An EDS is a sequence with `W(0)=0`, `W(1)=1` satisfying Ward's relation
//!
//! ```text
//!   W(m+n)·W(m−n) = W(m+1)·W(m−1)·W(n)² − W(n+1)·W(n−1)·W(m)².
//! ```
//!
//! For the sequence attached to a curve, `W(n) = ψ_n(P)` where `ψ_n` is the
//! `n`-th division polynomial, so `W(2)=2y`, `W(3)=ψ_3(P)`, `W(4)=ψ_4(P)`.
//! We extend with the standard duplication formulas (valid for `W(1)=1`):
//!
//! ```text
//!   W(2n+1) = W(n+2)·W(n)³ − W(n−1)·W(n+1)³
//!   W(2n)·W(2) = W(n)·( W(n+2)·W(n−1)² − W(n−2)·W(n+1)² )
//! ```
//!
//! Both right-hand sides reference only strictly smaller indices, so the
//! array fills left-to-right in `O(n)` field operations.
//!
//! ## What this is and is not
//!
//! Lauter–Stange prove EDS-Residue is solvable in sub-exponential time iff
//! the ECDLP is — so **this is not a generic-beating attack**, and the
//! module does not claim one.  It is an instrument for the *structural*
//! question the equivalence leaves open: does the Legendre sequence carry
//! exploitable bias on special curves, and how cheaply does the QR pattern
//! of the net localise the zero-lattice in practice?  See the research doc
//! for the honest scorecard.

use crate::cryptanalysis::ec_index_calculus::sqrt_mod_p;
use crate::ecc::field::FieldElement;
use crate::ecc::point::Point;
use num_bigint::BigUint;
use num_traits::{One, Zero};

// ── small modular helpers (BigUint, fixed modulus p) ───────────────────────

#[inline]
fn m_mul(x: &BigUint, y: &BigUint, p: &BigUint) -> BigUint {
    (x * y) % p
}

#[inline]
fn m_sub(x: &BigUint, y: &BigUint, p: &BigUint) -> BigUint {
    // (x − y) mod p without signed arithmetic.
    ((x % p) + p - (y % p)) % p
}

#[inline]
fn m_pow(x: &BigUint, e: u32, p: &BigUint) -> BigUint {
    x.modpow(&BigUint::from(e), p)
}

/// Modular inverse via Fermat (p prime).  `None` for zero.
fn m_inv(x: &BigUint, p: &BigUint) -> Option<BigUint> {
    if (x % p).is_zero() {
        return None;
    }
    Some(x.modpow(&(p - BigUint::from(2u32)), p))
}

/// Legendre symbol `(w | p)` returned as `+1` (QR), `-1` (non-residue) or
/// `0` (w ≡ 0).  `p` is an odd prime.
pub fn legendre(w: &BigUint, p: &BigUint) -> i8 {
    let r = (w % p).modpow(&((p - BigUint::one()) >> 1), p);
    if r.is_zero() {
        0
    } else if r.is_one() {
        1
    } else {
        // r == p − 1
        -1
    }
}

// ── the elliptic divisibility sequence ─────────────────────────────────────

/// Build `W_{E,P}(0..count)` over `F_p` for the curve `y² = x³ + ax + b`
/// and the affine point `P = (x, y)` with `y ≠ 0`.
///
/// Returns the residues `W(n) mod p`.  Requires `count ≥ 5`.
pub fn eds_sequence(
    p: &BigUint,
    a: &BigUint,
    b: &BigUint,
    x: &BigUint,
    y: &BigUint,
    count: usize,
) -> Vec<BigUint> {
    assert!(count >= 5, "need at least the four seed terms");
    let x = x % p;
    let y = y % p;
    let a = a % p;
    let b = b % p;

    let two = BigUint::from(2u32);
    let three = BigUint::from(3u32);
    let four = BigUint::from(4u32);
    let five = BigUint::from(5u32);
    let six = BigUint::from(6u32);
    let twelve = BigUint::from(12u32);
    let twenty = BigUint::from(20u32);
    let eight = BigUint::from(8u32);

    // W(0)=0, W(1)=1, W(2)=2y.
    let mut w: Vec<BigUint> = Vec::with_capacity(count);
    w.push(BigUint::zero());
    w.push(BigUint::one());
    w.push(m_mul(&two, &y, p));

    // W(3) = 3x⁴ + 6a x² + 12b x − a².
    let x2 = m_pow(&x, 2, p);
    let x4 = m_pow(&x, 4, p);
    let w3a = m_mul(&three, &x4, p);
    let w3b = m_mul(&m_mul(&six, &a, p), &x2, p);
    let w3c = m_mul(&m_mul(&twelve, &b, p), &x, p);
    let w3d = m_pow(&a, 2, p);
    let w3 = m_sub(&((w3a + w3b + w3c) % p), &w3d, p);
    w.push(w3);

    // W(4) = 4y( x⁶ + 5a x⁴ + 20b x³ − 5a² x² − 4ab x − 8b² − a³ ).
    let x3 = m_pow(&x, 3, p);
    let x6 = m_pow(&x, 6, p);
    let a2 = m_pow(&a, 2, p);
    let t_x6 = x6;
    let t_5ax4 = m_mul(&m_mul(&five, &a, p), &x4, p);
    let t_20bx3 = m_mul(&m_mul(&twenty, &b, p), &x3, p);
    let t_5a2x2 = m_mul(&m_mul(&five, &a2, p), &x2, p);
    let t_4abx = m_mul(&m_mul(&m_mul(&four, &a, p), &b, p), &x, p);
    let t_8b2 = m_mul(&eight, &m_pow(&b, 2, p), p);
    let t_a3 = m_pow(&a, 3, p);
    let pos = (t_x6 + t_5ax4 + t_20bx3) % p;
    let mut inner = m_sub(&pos, &t_5a2x2, p);
    inner = m_sub(&inner, &t_4abx, p);
    inner = m_sub(&inner, &t_8b2, p);
    inner = m_sub(&inner, &t_a3, p);
    let w4 = m_mul(&m_mul(&four, &y, p), &inner, p);
    w.push(w4);

    let w2_inv = m_inv(&w[2], p).expect("W(2)=2y must be invertible (y≠0)");

    for k in 5..count {
        let val = if k % 2 == 1 {
            // W(2n+1) = W(n+2)·W(n)³ − W(n−1)·W(n+1)³
            let n = (k - 1) / 2;
            let lhs = m_mul(&w[n + 2], &m_pow(&w[n], 3, p), p);
            let rhs = m_mul(&w[n - 1], &m_pow(&w[n + 1], 3, p), p);
            m_sub(&lhs, &rhs, p)
        } else {
            // W(2n)·W(2) = W(n)·( W(n+2)·W(n−1)² − W(n−2)·W(n+1)² )
            let n = k / 2;
            let t1 = m_mul(&w[n + 2], &m_pow(&w[n - 1], 2, p), p);
            let t2 = m_mul(&w[n - 2], &m_pow(&w[n + 1], 2, p), p);
            let bracket = m_sub(&t1, &t2, p);
            let num = m_mul(&w[n], &bracket, p);
            m_mul(&num, &w2_inv, p)
        };
        w.push(val);
    }
    w
}

// ── honest point arithmetic (validation anchor) ────────────────────────────

fn fe(v: &BigUint, p: &BigUint) -> FieldElement {
    FieldElement::new(v.clone(), p.clone())
}

/// Order of `P`: smallest `m > 0` with `[m]P = O` (brute force, toy sizes).
pub fn point_order(x: &BigUint, y: &BigUint, a: &BigUint, p: &BigUint, cap: u64) -> Option<u64> {
    let base = Point::Affine {
        x: fe(x, p),
        y: fe(y, p),
    };
    let a_fe = fe(a, p);
    let mut acc = Point::Infinity; // [0]P
    for m in 1..=cap {
        acc = acc.add(&base, &a_fe); // acc = [m]P
        if let Point::Infinity = acc {
            return Some(m);
        }
    }
    None
}

/// Locate the first point on the curve with `x ≥ x_start`, `y ≠ 0`, whose
/// order is at least `min_order`.  Returns `(x, y, order)`.
pub fn find_toy_point(
    p: &BigUint,
    a: &BigUint,
    b: &BigUint,
    x_start: u64,
    min_order: u64,
) -> Option<(BigUint, BigUint, u64)> {
    let pu: u64 = p.try_into().ok()?;
    for xi in x_start..pu {
        let x = BigUint::from(xi);
        let rhs = (m_pow(&x, 3, p) + m_mul(a, &x, p) + b) % p;
        if rhs.is_zero() {
            continue; // y = 0 → 2-torsion, W(2)=0
        }
        if let Some(y) = sqrt_mod_p(&rhs, p) {
            if y.is_zero() {
                continue;
            }
            if let Some(ord) = point_order(&x, &y, a, p, pu + 2) {
                if ord >= min_order {
                    return Some((x, y, ord));
                }
            }
        }
    }
    None
}

// ── analysis report ────────────────────────────────────────────────────────

/// Result of analysing one `(E, P)` pair.
#[derive(Clone, Debug)]
pub struct EdsReport {
    pub p: BigUint,
    pub a: BigUint,
    pub b: BigUint,
    pub px: BigUint,
    pub py: BigUint,
    /// `ord(P)` from honest point arithmetic.
    pub order: u64,
    /// First `n>0` with `W(n) ≡ 0`; must equal `order`.
    pub rank_of_apparition: u64,
    /// Did the validation `W(n)=0 ⟺ [n]P=O` hold for all tested `n`?
    pub apparition_consistent: bool,
    /// The shift multiplier: over `F_p`, `W(n+r) = A · Bⁿ · W(n)` for all
    /// `n`, where `r = ord(P)`.  `(A, B)`.
    pub mult_a: BigUint,
    pub mult_b: BigUint,
    /// Did the multiplier law `W(n+r) = A·Bⁿ·W(n)` verify for all `n`?
    pub mult_law_holds: bool,
    /// Legendre characters of the multiplier: `(χ(A), χ(B))`.  These are the
    /// invariants the EDS-Residue period hangs on.
    pub chi_a: i8,
    pub chi_b: i8,
    /// Period of `W(n) mod p` = `r · j`, the smallest `j` closing the
    /// multiplier (`Bʲ=1` and `Aʲ·B^{r·j(j-1)/2}=1`).
    pub period_w: u64,
    /// Period of the Legendre sequence `χ(n)` = `r · j_χ`.
    pub period_chi: u64,
    /// `(#QR, #NQR, #zero)` over the apparition block `n ∈ [1, r]`.
    pub qr_balance: (u64, u64, u64),
}

/// `(±1)^e` from a sign `s ∈ {+1,-1}` and exponent parity.
#[inline]
fn pm_pow(s: i8, e_is_odd: bool) -> i8 {
    if s == 1 || !e_is_odd {
        1
    } else {
        -1
    }
}

/// Full analysis of the EDS attached to `P` on `y²=x³+ax+b / F_p`.
pub fn analyze(p: &BigUint, a: &BigUint, b: &BigUint, px: &BigUint, py: &BigUint) -> EdsReport {
    let pu: u64 = p.clone().try_into().expect("toy prime fits in u64");
    let order = point_order(px, py, a, p, pu + 2).expect("P has finite order");
    let r = order as usize;

    // Need ≥ 2 apparition blocks plus headroom to verify the multiplier law.
    let count = (2 * r + 4).max(8);
    let w = eds_sequence(p, a, b, px, py, count);

    // Rank of apparition: first zero.
    let rank = (1..count)
        .find(|&n| w[n].is_zero())
        .map(|n| n as u64)
        .unwrap_or(0);

    // Validate W(n)=0 ⟺ [n]P=O against honest point arithmetic.
    let a_fe = fe(a, p);
    let base = Point::Affine {
        x: fe(px, p),
        y: fe(py, p),
    };
    let mut consistent = true;
    let mut acc = Point::Infinity; // [0]P
    for wn in w.iter() {
        let is_identity = matches!(acc, Point::Infinity);
        if wn.is_zero() != is_identity {
            consistent = false;
            break;
        }
        acc = acc.add(&base, &a_fe);
    }

    // Shift multiplier (A, B) from W(n+r) = A·Bⁿ·W(n):
    //   B = W(r+2) / ( W(2)·W(r+1) ),   A = W(r+1) / B.
    let w2_inv = m_inv(&w[2], p).expect("W(2) invertible");
    let wr1_inv = m_inv(&w[r + 1], p).expect("W(r+1) nonzero (m ∤ r+1)");
    let mult_b = m_mul(&w[r + 2], &m_mul(&w2_inv, &wr1_inv, p), p);
    let b_inv = m_inv(&mult_b, p).expect("B nonzero");
    let mult_a = m_mul(&w[r + 1], &b_inv, p);

    // Verify the law across a full block.
    let mut mult_law_holds = true;
    let mut bn = BigUint::one(); // Bⁿ, starting n=0
    for n in 0..r {
        let predicted = m_mul(&mult_a, &m_mul(&bn, &w[n], p), p);
        if predicted != w[n + r] {
            mult_law_holds = false;
            break;
        }
        bn = m_mul(&bn, &mult_b, p);
    }

    let chi_a = legendre(&mult_a, p);
    let chi_b = legendre(&mult_b, p);

    // Period of W: smallest j with Bʲ=1 and Aʲ·B^{r·j(j-1)/2}=1.  π_W = r·j.
    let period_w = {
        let mut found = None;
        for j in 1u64..=pu {
            let bj = mult_b.modpow(&BigUint::from(j), p);
            if !bj.is_one() {
                continue;
            }
            // exponent r·j·(j-1)/2 (mod p−1 handled by modpow)
            let exp = BigUint::from(r as u64) * BigUint::from(j) * BigUint::from(j - 1) / 2u32;
            let lhs = m_mul(
                &mult_a.modpow(&BigUint::from(j), p),
                &mult_b.modpow(&exp, p),
                p,
            );
            if lhs.is_one() {
                found = Some(j);
                break;
            }
        }
        order * found.unwrap_or(0)
    };

    // Period of χ: smallest j with ε^j=1 and δ^j·ε^{r·j(j-1)/2}=1, χ-arithmetic.
    let period_chi = {
        let mut found = None;
        for j in 1u64..=8 {
            // coefficient of n: ε^j must be 1
            if pm_pow(chi_b, j % 2 == 1) != 1 {
                continue;
            }
            let exp_parity = {
                // parity of r·j·(j-1)/2
                let e = (r as u128) * (j as u128) * ((j as u128) - 1) / 2;
                e % 2 == 1
            };
            let lhs = pm_pow(chi_a, j % 2 == 1) * pm_pow(chi_b, exp_parity);
            if lhs == 1 {
                found = Some(j);
                break;
            }
        }
        order * found.unwrap_or(0)
    };

    // QR balance over the apparition block n ∈ [1, r].
    let (mut qr, mut nqr, mut zero) = (0u64, 0u64, 0u64);
    for n in 1..=r {
        match legendre(&w[n], p) {
            1 => qr += 1,
            -1 => nqr += 1,
            _ => zero += 1,
        }
    }

    EdsReport {
        p: p.clone(),
        a: a.clone(),
        b: b.clone(),
        px: px.clone(),
        py: py.clone(),
        order,
        rank_of_apparition: rank,
        apparition_consistent: consistent,
        mult_a,
        mult_b,
        mult_law_holds,
        chi_a,
        chi_b,
        period_w,
        period_chi,
        qr_balance: (qr, nqr, zero),
    }
}

// ── the 2-D net zero-lattice (DL recovery illustration) ────────────────────

/// For `P` of order `m` and `Q = [k]P`, the 2-D elliptic-net value
/// `W(a,b)` vanishes iff `[a]P + [b]Q = O`, i.e. iff `a + bk ≡ 0 (mod m)`.
/// Recover `k` from any net zero with `b` invertible mod `m`.
///
/// This uses honest point arithmetic to *locate* the vanishing locus (the
/// point of the net algorithm is that each `W(a,b)` is computable by the
/// recurrence without per-cell scalar multiplication; finding a zero is
/// still a search, which is why the whole construction is ECDLP-hard).
/// Returns `(k, the (a,b) zero used)`.
pub fn recover_dl_from_net_zeros(
    px: &BigUint,
    py: &BigUint,
    a: &BigUint,
    p: &BigUint,
    order: u64,
    k_true: u64,
    grid: u64,
) -> Option<(u64, (i64, i64))> {
    let a_fe = fe(a, p);
    let base = Point::Affine {
        x: fe(px, p),
        y: fe(py, p),
    };
    // Q = [k]P.
    let q = base.scalar_mul(&BigUint::from(k_true), &a_fe);

    for b_idx in 1..=(grid as i64) {
        for a_idx in 0..(grid as i64) {
            // [a]P + [b]Q
            let ap = if a_idx == 0 {
                Point::Infinity
            } else {
                base.scalar_mul(&BigUint::from(a_idx as u64), &a_fe)
            };
            let bq = match &q {
                Point::Infinity => Point::Infinity,
                _ => q.scalar_mul(&BigUint::from(b_idx as u64), &a_fe),
            };
            let sum = ap.add(&bq, &a_fe);
            if let Point::Infinity = sum {
                // a + b·k ≡ 0 (mod m)  ⇒  k ≡ −a · b⁻¹ (mod m)
                let m = BigUint::from(order);
                let b_big = BigUint::from(b_idx as u64) % &m;
                if let Some(b_inv) = m_inv(&b_big, &m) {
                    let neg_a = m_sub(&BigUint::zero(), &(BigUint::from(a_idx as u64) % &m), &m);
                    let k = (neg_a * b_inv) % &m;
                    let k_u: u64 = k.try_into().ok()?;
                    return Some((k_u, (a_idx, b_idx)));
                }
            }
        }
    }
    None
}

// ── markdown report ─────────────────────────────────────────────────────────

/// Render an [`EdsReport`] as a Markdown block for the research log.
pub fn format_report(r: &EdsReport) -> String {
    let (qr, nqr, zero) = r.qr_balance;
    let total = qr + nqr + zero;
    let bias = if (qr + nqr) > 0 {
        (qr as f64 - nqr as f64) / (qr + nqr) as f64
    } else {
        0.0
    };
    let sgn = |s: i8| match s {
        1 => "+1",
        -1 => "-1",
        _ => "0",
    };
    format!(
        "## EDS-Residue analysis\n\n\
         - Curve: `y² = x³ + {a}x + {b}  (mod {p})`\n\
         - Point `P = ({px}, {py})`, `ord(P) = {ord}`\n\
         - Rank of apparition (first `W(n)=0`): **{rank}**  {anchor}\n\
         - Apparition law `W(n)=0 ⟺ [n]P=O`: **{cons}**\n\
         - Shift multiplier `W(n+r) = A·Bⁿ·W(n)`: `A = {ma}`, `B = {mb}`  ({mlaw})\n\
         - Multiplier characters: `χ(A) = {ca}`, `χ(B) = {cb}`\n\
         - Period of `W(n) mod p`: `π_W = {pw}` = `{jw}·r`\n\
         - Period of Legendre seq `χ(n)`: `π_χ = {pc}` = `{jc}·r`\n\
         - QR balance over apparition block `[1,r]`: QR={qr}, NQR={nqr}, zero={zero} (n={total})\n\
         - Residue bias `(QR−NQR)/(QR+NQR)` = `{bias:+.4}`\n",
        a = r.a,
        b = r.b,
        p = r.p,
        px = r.px,
        py = r.py,
        ord = r.order,
        rank = r.rank_of_apparition,
        anchor = if r.rank_of_apparition == r.order {
            "✓ = ord(P)"
        } else {
            "✗ MISMATCH"
        },
        cons = if r.apparition_consistent { "✓ holds" } else { "✗ FAILS" },
        ma = r.mult_a,
        mb = r.mult_b,
        mlaw = if r.mult_law_holds { "✓ verified" } else { "✗ FAILS" },
        ca = sgn(r.chi_a),
        cb = sgn(r.chi_b),
        pw = r.period_w,
        jw = if r.order > 0 { r.period_w / r.order } else { 0 },
        pc = r.period_chi,
        jc = if r.order > 0 { r.period_chi / r.order } else { 0 },
        qr = qr,
        nqr = nqr,
        zero = zero,
        total = total,
        bias = bias,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A concrete toy curve: y² = x³ + 37x + 2 (mod 1009).
    fn toy() -> (BigUint, BigUint, BigUint) {
        (
            BigUint::from(1009u32),
            BigUint::from(37u32),
            BigUint::from(2u32),
        )
    }

    #[test]
    fn rank_of_apparition_equals_order() {
        let (p, a, b) = toy();
        let (px, py, ord) = find_toy_point(&p, &a, &b, 1, 7).expect("toy point exists");
        let r = analyze(&p, &a, &b, &px, &py);
        assert_eq!(r.order, ord);
        // The defining property: first zero of the EDS is exactly ord(P).
        assert_eq!(
            r.rank_of_apparition, r.order,
            "rank of apparition must equal ord(P)"
        );
    }

    #[test]
    fn apparition_law_holds() {
        let (p, a, b) = toy();
        let (px, py, _) = find_toy_point(&p, &a, &b, 1, 7).unwrap();
        let r = analyze(&p, &a, &b, &px, &py);
        // W(n) = 0 ⟺ [n]P = O, checked against honest point arithmetic.
        assert!(r.apparition_consistent, "apparition law W(n)=0 ⟺ [n]P=O must hold");
    }

    #[test]
    fn legendre_matches_euler() {
        let p = BigUint::from(1009u32);
        // 2 is a QR mod 1009? 1009 ≡ 1 (mod 8) ⇒ 2 is a QR.
        assert_eq!(legendre(&BigUint::from(0u32), &p), 0);
        // sanity: a square is always a QR.
        for v in 1u32..50 {
            let sq = (v * v) % 1009;
            assert_eq!(legendre(&BigUint::from(sq), &p), 1, "{} should be QR", sq);
        }
    }

    #[test]
    fn net_zero_lattice_recovers_discrete_log() {
        let (p, a, b) = toy();
        let (px, py, ord) = find_toy_point(&p, &a, &b, 1, 7).unwrap();
        // Choose a non-trivial true discrete log and recover it from the
        // vanishing locus of the 2-D net.
        let k_true = (ord / 2).max(2) % ord;
        let grid = ord + 1;
        let (k_rec, (za, zb)) =
            recover_dl_from_net_zeros(&px, &py, &a, &p, ord, k_true, grid).expect("found a net zero");
        assert_eq!(k_rec % ord, k_true % ord, "net zero must encode the DL");
        // Verify the recovered relation: za + zb·k ≡ 0 (mod ord).
        let lhs = (za as u64 + (zb as u64) * (k_true % ord)) % ord;
        assert_eq!(lhs, 0, "zero must lie on the slope-k lattice");
    }

    #[test]
    fn multiplier_law_and_periods() {
        let (p, a, b) = toy();
        let (px, py, _) = find_toy_point(&p, &a, &b, 1, 7).unwrap();
        let r = analyze(&p, &a, &b, &px, &py);
        // The Ward/Lauter–Stange shift law must verify exactly.
        assert!(r.mult_law_holds, "W(n+r) = A·Bⁿ·W(n) must hold");
        // χ period must divide W period (χ is a function of W).
        assert!(r.period_w > 0 && r.period_chi > 0);
        assert_eq!(r.period_w % r.period_chi, 0, "π_χ must divide π_W");
        // Both periods are integer multiples of the rank of apparition.
        assert_eq!(r.period_w % r.order, 0);
        assert_eq!(r.period_chi % r.order, 0);
        // QR balance accounts for the whole apparition block [1, r].
        let (qr, nqr, zero) = r.qr_balance;
        assert_eq!(qr + nqr + zero, r.order);
        // χ(A), χ(B) are genuine signs.
        assert!(r.chi_a != 0 && r.chi_b != 0);
    }

    #[test]
    fn report_renders() {
        let (p, a, b) = toy();
        let (px, py, _) = find_toy_point(&p, &a, &b, 1, 7).unwrap();
        let r = analyze(&p, &a, &b, &px, &py);
        let md = format_report(&r);
        assert!(md.contains("EDS-Residue analysis"));
        assert!(md.contains("Shift multiplier"));
    }
}
