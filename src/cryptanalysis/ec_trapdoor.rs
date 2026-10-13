//! # Teske's Elliptic-Curve Trapdoor System (J. Cryptology 19(1), 2006).
//!
//! Edlyn Teske showed how to construct an elliptic curve `E/F_{2^N}`
//! that **looks generic** to an outside observer but admits a low-genus
//! Weil-descent (GHS) attack known only to the trapdoor owner.  The
//! trapdoor is the *correct factorisation* `N = n·l` to descend through
//! (and, more concretely, which subfield contains the curve's `b`
//! parameter or its square root).
//!
//! This module ships the **complete** primitives behind the construction
//! and verification:
//!
//! 1. [`FieldTower`] — the field-extension lattice `F_2 ⊂ F_{2^l} ⊂
//!    F_{2^N}` with the Frobenius `σ: x ↦ x^{2^l}` and subfield
//!    membership tests.
//! 2. [`magic_number`] — the GHS magic number `m(b) = dim_{F_2}
//!    span{(1, σ^i(√b))}` of Maurer–Menezes–Teske (2002), eq. (5).
//! 3. [`ghs_genus`] — the genus of the GHS hyperelliptic descent
//!    curve, `g = 2^{m−1}`, or `2^{m−1} − 1` in the type-I case
//!    `(x − 1) ∤ Ord_b` ([`ghs_type_i`]).
//! 4. [`audit_curve`] — the **public auditor view**: for every
//!    factorisation `N = n_i · l_i`, compute the magic number and
//!    GHS-attack genus, returning the minimum.
//! 5. [`construct_trapdoor_curve`] — the **constructor view**: pick
//!    `(a, b)` so a *specific* factorisation yields a small genus.
//!
//! The actual descent (producing `C/F_{2^l}` and the homomorphism
//! `Φ: E(F_{2^N}) → Jac(C)(F_{2^l})`) is in
//! [`crate::cryptanalysis::ghs_descent`].  The HCDLP solver runs on
//! the Jacobian arithmetic provided by
//! [`crate::binary_ecc::hyperelliptic`].
//!
//! ## The GHS magic number — what it is and why it matters
//!
//! Let `K = F_{2^N}`, `k = F_{2^l}` with `N = n·l`, and
//! `σ ∈ Gal(K/k)` the Frobenius `x ↦ x^{2^l}` (order `n` on `K`).
//! Consider `E/K : y² + x·y = x³ + a·x² + b` (ordinary, so `b ≠ 0`).
//!
//! Substituting `y = x·t` and absorbing the resulting `(√b/x)²`
//! square via an Artin–Schreier shift `t ↦ t + √b/x` yields the
//! defining equation of an Artin–Schreier extension `K(C)/K(x)`:
//! ```text
//!     t² + t = x + a + (√b)/x
//! ```
//! Apply `σ` to get the conjugate functions `x + σ^i(a) + σ^i(√b)/x`
//! and take the compositum.  It is determined by their `F_2`-span modulo
//! `℘(K(x)) = {h² + h}`.  The constants `σ^i(a) − a` have absolute trace
//! zero, so they lie in `℘(K)` and drop out; a relation
//! `Σ c_i σ^i(f) ∈ ℘(K(x))` needs `Σ c_i = 0` (the `x` term) and
//! `Σ c_i σ^i(√b) = 0` (the `1/x` term).  The dimension of the span is
//! therefore (Maurer–Menezes–Teske, LMS J. Comput. Math. 5 (2002), eq. (5))
//! ```text
//!     m(b) := dim_{F_2} span_{F_2} { (1, σ^i(√b)) : i = 0, …, n−1 } ⊂ F_2 × K
//! ```
//! Write `Ord_b ∈ F_2[x]` for the `σ`-annihilator of `√b`.  MMT Theorem
//! 5(i): `m(b) = deg Ord_b + c`, with `c = 1` when `(x − 1) ∤ Ord_b` and
//! `c = 0` otherwise.  The plain orbit dimension `deg Ord_b`
//! ([`plain_orbit_dim`], this module's definition before 2026-10-13) is
//! one too small in the `c = 1` case.
//!
//! The genus of the compositum is the number of its quadratic subfields
//! `w² + w = u` of genus one, i.e. of the `(1, β)` in the span with
//! `β ≠ 0`:
//! ```text
//!     g(C) = 2^{m−1}      when (x − 1) | Ord_b   (type II)
//!     g(C) = 2^{m−1} − 1  when (x − 1) ∤ Ord_b   (type I: (1, 0) is in the span)
//! ```
//! (Galbraith–Hess–Smart 2002; Hess, LMS J. Comp. Math. 7 (2004).)  For
//! `n` odd, type I is exactly `Tr_{K/k}(b) = 0` with `b ∉ k`.  Whether GHS
//! applies at all when `Tr_{K/F_2}(a) = 1` is MMT Lemma 6, condition (8)
//! (`n` odd, or `Tr_{K/k}(b) ≠ 0`); this module does not check it.
//!
//! **For an attacker, small `m` is fatal:** the descent transports
//! ECDLP on `E(K)` (size `~ q^n`) to HCDLP on `Jac(C)(k)` (size
//! `~ q^g`, `q = 2^l`), and HCDLP for `g ≤ 4` is *much* easier than
//! ECDLP on the equivalent-size group.
//!
//! **For the trapdoor:** Teske picks `b` with `√b ∈ F_{2^{m·l}} \
//! F_{2^l}` for a small target `m`.  Then `σ` has orbit-size `m` on
//! `√b`; when that orbit is `F_2`-independent, `Ord_b = x^m − 1`, so the
//! magic number is `m` and the genus `2^{m−1}` is small.  An attacker
//! observing `E` over `K` cannot tell *which* factorisation
//! `N = n_i · l_i` is the "right" one without computing the magic
//! number for every divisor of `N` — and for highly composite `N`,
//! that's its own little discrete-log problem.
//!
//! ## What's in this module
//!
//! Field-arithmetic side:
//! - [`FieldTower::frobenius`], [`FieldTower::frobenius_iter`].
//! - [`FieldTower::sqrt`] — characteristic-2 unique square root,
//!   `√x = x^{2^{N−1}}`.
//! - [`FieldTower::is_in_subfield`] — `x ∈ F_{2^{sub_l}}` ⇔
//!   `σ^{sub_l/gcd}(x) = x`.
//!
//! Trapdoor analysis:
//! - [`f2_dim`] — `F_2`-rank of a list of elements of `K` via
//!   Gaussian elimination on packed bit-vectors.
//! - [`magic_number`], [`ghs_type_i`], [`plain_orbit_dim`].
//! - [`ghs_genus`], [`ghs_genus_with_type`].
//! - [`audit_curve`] — exhaustive over every factorisation of `N`.
//! - [`construct_trapdoor_curve`] — search loop producing `(a, b)`
//!   for a chosen genus target on a chosen factorisation.
//!
//! ## Honest scope
//!
//! - `m` depends on `b` only; [`magic_number_full`] takes `a` for API
//!   compatibility and equals [`magic_number`].  Its pre-2026-10-13 value,
//!   which added `(σ − 1)·σ^i(a)` as generators, is
//!   [`legacy_orbit_dim_with_a`].
//! - [`ghs_genus`] returns `2^{m−1}`; [`ghs_genus_with_type`] and
//!   [`audit_curve`] apply the type-I correction.
//! - The applicability condition for `Tr_{K/F_2}(a) = 1` (MMT Lemma 6)
//!   and the descent itself are not checked here.

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// A field tower `F_2 ⊂ F_{2^l} ⊂ F_{2^N}` with `N = n·l`.
///
/// We store **only the top field's irreducible polynomial**; subfield
/// membership is decided by the Frobenius fixed-point test
/// `x ∈ F_{2^l} ⇔ σ(x) = x`.  This sidesteps the need for an
/// explicit embedding of `F_{2^l}` into `F_{2^N}` (which would
/// require compatible irreducibles).
#[derive(Clone, Debug)]
pub struct FieldTower {
    /// `[K : F_2]`.
    pub big_n: u32,
    /// `[K : k]` — the order of `σ` on `K`.
    pub n: u32,
    /// `[k : F_2]`.  `N = n · l`.
    pub l: u32,
    /// Reduction polynomial of `K = F_{2^N}`.
    pub irr: IrreduciblePoly,
}

impl FieldTower {
    /// Construct from `(N, n, l, irr)` after checking `N = n · l`.
    pub fn new(big_n: u32, n: u32, l: u32, irr: IrreduciblePoly) -> Self {
        assert_eq!(big_n, n * l, "N must equal n·l");
        assert_eq!(irr.degree, big_n);
        Self { big_n, n, l, irr }
    }

    /// `σ(x) = x^{2^l}` — the Frobenius generator of `Gal(K/k)`.
    pub fn frobenius(&self, x: &F2mElement) -> F2mElement {
        x.square_k_times(self.l, &self.irr)
    }

    /// `σ^i(x) = x^{2^{i·l}}`.
    pub fn frobenius_iter(&self, x: &F2mElement, i: u32) -> F2mElement {
        x.square_k_times(i * self.l, &self.irr)
    }

    /// `√x` in `F_{2^N}` — the unique `c` with `c² = x`.
    ///
    /// Computed as `c = x^{2^{N−1}}` since the absolute Frobenius
    /// `φ(x) = x²` has order `N` on `F_{2^N}`, so its inverse is
    /// `φ^{N−1}`.
    pub fn sqrt(&self, x: &F2mElement) -> F2mElement {
        x.square_k_times(self.big_n - 1, &self.irr)
    }

    /// `true` iff `x ∈ F_{2^{sub_l}}`, where `sub_l` must divide
    /// `big_n`.  Equivalent to `x^{2^{sub_l}} = x`.
    pub fn is_in_subfield(&self, x: &F2mElement, sub_l: u32) -> bool {
        assert!(self.big_n.is_multiple_of(sub_l), "sub_l must divide N");
        x.square_k_times(sub_l, &self.irr) == *x
    }

    /// Order of the σ-orbit of `x`: the smallest positive `i ≤ n`
    /// with `σ^i(x) = x`, equivalently the smallest `i ≥ 1` such
    /// that `x ∈ F_{2^{i·l}}`.  Always divides `n`.
    pub fn frobenius_orbit_length(&self, x: &F2mElement) -> u32 {
        // Generate divisors of n in ascending order and test.
        let mut divisors: Vec<u32> = (1..=self.n).filter(|d| self.n.is_multiple_of(*d)).collect();
        divisors.sort_unstable();
        for d in divisors {
            // σ^d(x) = x  ⇔  x ∈ F_{2^{d·l}}.
            let xs = x.square_k_times(d * self.l, &self.irr);
            if &xs == x {
                return d;
            }
        }
        self.n
    }
}

/// Compute the `F_2`-dimension of the subspace of `F_{2^m}`
/// (viewed as `F_2^m`) spanned by `elements`.
///
/// Uses Gaussian elimination over `F_2` on packed `u64` bit-vectors.
pub fn f2_dim(elements: &[F2mElement], m: u32) -> u32 {
    if elements.is_empty() {
        return 0;
    }
    let n_words = m.div_ceil(64) as usize;
    // Copy bit-vectors into a mutable matrix.
    let mut rows: Vec<Vec<u64>> = elements
        .iter()
        .map(|e| {
            let mut row = vec![0u64; n_words];
            for (i, w) in e.raw_bits().iter().enumerate() {
                if i < n_words {
                    row[i] = *w;
                }
            }
            row
        })
        .collect();
    let mut rank = 0usize;
    for col in 0..m {
        let w = (col / 64) as usize;
        let b = col % 64;
        // Find a row at index ≥ rank with bit `col` set.
        let mut pivot = None;
        for r in rank..rows.len() {
            if (rows[r][w] >> b) & 1 == 1 {
                pivot = Some(r);
                break;
            }
        }
        if let Some(p) = pivot {
            rows.swap(rank, p);
            // Eliminate from all other rows.
            for r in 0..rows.len() {
                if r != rank && (rows[r][w] >> b) & 1 == 1 {
                    for i in 0..n_words {
                        rows[r][i] ^= rows[rank][i];
                    }
                }
            }
            rank += 1;
        }
        if rank >= rows.len() {
            break;
        }
    }
    rank as u32
}

/// `F_2`-dimension of the plain Frobenius orbit `span{σ^i(√b) : 0 ≤ i < n}`.
///
/// This is `deg Ord_b`, where `Ord_b ∈ F_2[x]` is the `σ`-annihilator of
/// `√b` (Maurer–Menezes–Teske's notation).  It is **not** the GHS magic
/// number: it omits the leading `1` of MMT eq. (5) and is one too small
/// exactly when `(x − 1) ∤ Ord_b`.  Until 2026-10-13 this was the body of
/// [`magic_number`]; it is kept under this name so older reports that
/// recorded it can be reproduced.
pub fn plain_orbit_dim(tower: &FieldTower, b: &F2mElement) -> u32 {
    f2_dim(&sqrt_b_orbit(tower, b), tower.big_n)
}

/// The value the pre-2026-10-13 [`magic_number_full`] returned: the plain
/// `√b` orbit plus the constants `(σ − 1)·σ^i(a)`, ranked in one space.
///
/// Those constants are always in `℘(K) = {α² + α}` (their absolute trace is
/// zero), so they never change the Artin–Schreier compositum; adding them to
/// the `√b` coefficients mixes the constant and `1/x` coefficients of the
/// Artin–Schreier functions.  Kept only to reproduce older reports.
pub fn legacy_orbit_dim_with_a(tower: &FieldTower, a: &F2mElement, b: &F2mElement) -> u32 {
    let mut gens = sqrt_b_orbit(tower, b);
    // (σ - 1)·σ^i(a) = σ^{i+1}(a) + σ^i(a) (char 2)
    let mut cur_a = a.clone();
    for _ in 0..tower.n {
        let next_a = tower.frobenius(&cur_a);
        gens.push(cur_a.add(&next_a));
        cur_a = next_a;
    }
    f2_dim(&gens, tower.big_n)
}

fn sqrt_b_orbit(tower: &FieldTower, b: &F2mElement) -> Vec<F2mElement> {
    let mut orbit = Vec::with_capacity(tower.n as usize);
    let mut cur = tower.sqrt(b);
    for _ in 0..tower.n {
        orbit.push(cur.clone());
        cur = tower.frobenius(&cur);
    }
    orbit
}

/// The GHS magic number of `E: y² + x·y = x³ + a·x² + b` relative to the
/// [`FieldTower`] (Maurer–Menezes–Teske, LMS J. Comput. Math. 5 (2002),
/// eq. (5)):
///
/// `m(b) := dim_{F_2} span_{F_2} { (1, σ^i(√b)) : 0 ≤ i < n } ⊂ F_2 × K`.
///
/// By MMT Theorem 5(i), `m(b) = deg Ord_b + c` with `c = 1` when
/// `(x − 1) ∤ Ord_b` and `c = 0` otherwise, i.e. [`plain_orbit_dim`] plus
/// [`ghs_type_i`].  The value depends on `b` only; `a` enters the GHS
/// attack through its applicability (MMT Lemma 6), not through `m`.
///
/// `b = 0` (not a curve) returns `0`.
pub fn magic_number(tower: &FieldTower, _a: &F2mElement, b: &F2mElement) -> u32 {
    mmt_ranks(tower, b).1
}

/// The GHS magic number with the curve's `a` passed through.  Equal to
/// [`magic_number`]: the Artin–Schreier functions are
/// `x + σ^i(a) + σ^i(√b)/x`, and their constant parts differ by elements of
/// `℘(K)`, so `a` does not change the compositum or `m`.  (Before
/// 2026-10-13 this added `(σ − 1)·σ^i(a)` as generators; that value is
/// [`legacy_orbit_dim_with_a`].)
pub fn magic_number_full(tower: &FieldTower, a: &F2mElement, b: &F2mElement) -> u32 {
    magic_number(tower, a, b)
}

/// Whether `(x − 1) ∤ Ord_b`, i.e. whether some odd-size subset of the
/// `σ`-orbit of `√b` sums to zero, i.e. whether `(1, 0)` lies in the MMT
/// span.  This is MMT Theorem 5(i)'s `c = 1` case, and it is exactly the
/// GHS genus-`2^{m−1} − 1` ("type I") case: the genus of the compositum is
/// the number of `(1, β)` in the span with `β ≠ 0`, which is `2^{m−1}`
/// less one when `(1, 0)` is present.  `b = 0` returns `false`.
pub fn ghs_type_i(tower: &FieldTower, b: &F2mElement) -> bool {
    let (plain, mmt) = mmt_ranks(tower, b);
    plain != 0 && mmt == plain + 1
}

/// `(plain, mmt)`: the ranks of `{σ^i(√b)}` and of `{(1, σ^i(√b))}`.  The
/// pairs are encoded by setting bit `big_n` (the `F_2` coordinate) on each
/// orbit element.  `mmt` is forced to `0` when `b = 0`.
fn mmt_ranks(tower: &FieldTower, b: &F2mElement) -> (u32, u32) {
    let orbit = sqrt_b_orbit(tower, b);
    let plain = f2_dim(&orbit, tower.big_n);
    if plain == 0 {
        return (0, 0);
    }
    let lifted: Vec<F2mElement> = orbit
        .iter()
        .map(|e| {
            let mut v = e.to_biguint();
            v.set_bit(u64::from(tower.big_n), true);
            F2mElement::from_biguint(&v, tower.big_n + 1)
        })
        .collect();
    (plain, f2_dim(&lifted, tower.big_n + 1))
}

/// GHS hyperelliptic-descent genus from the magic number.
///
/// `g = 2^{m−1}` in the type-II case `(x − 1) | Ord_b`; the type-I case
/// `(x − 1) ∤ Ord_b` has `2^{m−1} − 1` (see [`ghs_type_i`] and
/// [`ghs_genus_with_type`]).
pub fn ghs_genus(magic_m: u32) -> BigUint {
    if magic_m == 0 {
        BigUint::zero()
    } else {
        BigUint::one() << (magic_m - 1)
    }
}

/// GHS genus for a known type branch, without fixed-width truncation.
///
/// This is the common exact-arithmetic primitive used by
/// [`ghs_genus_with_type`] and [`audit_curve`].  Keeping the power of two as a
/// [`BigUint`] matters for ordinary cryptographic field degrees: a perfectly
/// ordinary screen can have `magic_m >= 33`, where the former `u32` shift
/// overflowed before the result could be reported.
pub fn ghs_genus_for_type(magic_m: u32, type_i: bool) -> BigUint {
    let generic = ghs_genus(magic_m);
    if type_i && magic_m != 0 {
        generic - BigUint::one()
    } else {
        generic
    }
}

/// Return `(genus, is_type_i)` for the MMT [`magic_number`] `m`.  Type I
/// (genus `2^{m−1} − 1`) is `(x − 1) ∤ Ord_b` ([`ghs_type_i`]); otherwise
/// the genus is `2^{m−1}`.
///
/// Before 2026-10-13 this paired the plain orbit dimension with the test
/// "orbit length == m", which inverts the type in both directions (for
/// `n` prime and `√b` of full orbit it reported `2^{n−1} − 1` where the
/// genus is `2^{n−1}`, and `2^{n−2}` where it is `2^{n−1} − 1`).
pub fn ghs_genus_with_type(tower: &FieldTower, b: &F2mElement) -> (BigUint, bool) {
    let (plain, m) = mmt_ranks(tower, b);
    if m == 0 {
        return (BigUint::zero(), false);
    }
    let type_i = m == plain + 1;
    (ghs_genus_for_type(m, type_i), type_i)
}

/// **Auditor view** of a curve `E/F_{2^N}`: try every non-trivial
/// factorisation `N = n·l` and report the descent magic number /
/// genus for each.  Returns the rows sorted by `(genus, l)` ascending
/// (smaller-genus attacks first).
pub fn audit_curve(
    big_n: u32,
    big_irr: &IrreduciblePoly,
    a: &F2mElement,
    b: &F2mElement,
) -> Vec<DescentRow> {
    let mut rows = Vec::new();
    for l in 1..big_n {
        if !big_n.is_multiple_of(l) {
            continue;
        }
        let n = big_n / l;
        if n < 2 {
            continue;
        }
        let tower = FieldTower::new(big_n, n, l, big_irr.clone());
        let m = magic_number_full(&tower, a, b);
        let (g, type_i) = ghs_genus_with_type(&tower, b);
        rows.push(DescentRow {
            n,
            l,
            magic_m: m,
            genus: g,
            type_i,
        });
    }
    // Sort: smallest genus first (most threatening descent); ties
    // broken by larger l (more "work" in the descent curve field).
    rows.sort_by(|x, y| x.genus.cmp(&y.genus).then(y.l.cmp(&x.l)));
    rows
}

/// One row of an [`audit_curve`] report.
#[derive(Clone, Debug, Default)]
pub struct DescentRow {
    /// `[K : k]` for this factorisation.
    pub n: u32,
    /// `[k : F_2]`.
    pub l: u32,
    /// MMT magic number `dim span{(1, σ^i(√b))}` for this factorisation.
    pub magic_m: u32,
    /// GHS genus `2^{m−1}` or `2^{m−1} − 1`.
    pub genus: BigUint,
    /// Type I, `(x − 1) ∤ Ord_b` (genus `2^{m−1} − 1`)?
    pub type_i: bool,
}

/// **Trapdoor constructor**: given a target factorisation `N = n·l`
/// and a target magic number `m_target` (small ⇒ small descent
/// genus), search for `b ∈ F_{2^N}` with `√b ∈ F_{2^{m_target·l}}`
/// whose σ-orbit has length exactly `m_target` and is `F_2`-independent.
/// Then `Ord_b = x^{m_target} − 1`, which `(x − 1)` divides, so the MMT
/// magic number is `m_target` and the curve is type II (genus
/// `2^{m_target−1}`).
///
/// Returns the curve coefficients `(a, b)` and the resulting
/// audit row.  We fix `a = 1` (so the `a`-contribution to the
/// magic number vanishes).
///
/// **Returns `None`** if the search budget is exhausted; this is
/// rare for small `N` but can happen for highly constrained towers.
///
/// The caller is responsible for checking the resulting curve is
/// not anomalous or otherwise weak via *other* attacks (e.g.
/// supersingular detection) — this routine only ensures the
/// GHS-descent trapdoor structure.
pub fn construct_trapdoor_curve(
    big_n: u32,
    n: u32,
    l: u32,
    m_target: u32,
    big_irr: &IrreduciblePoly,
    max_attempts: u32,
) -> Option<TrapdoorCurve> {
    assert!(m_target >= 1 && m_target <= n);
    let tower = FieldTower::new(big_n, n, l, big_irr.clone());
    let a = F2mElement::one(big_n);

    // Strategy: enumerate small candidates for √b directly.  We
    // want √b ∈ F_{2^{m_target·l}} but ∉ F_{2^{(m_target−1)·l}}
    // (or any smaller subfield containing it).  Brute-force search
    // over small bit-patterns suffices for the toy parameters this
    // module targets.
    let sub_l = m_target * l;
    let _strictly_below_l = (m_target - 1) * l;

    let mut tried = 0u32;
    let mut candidate_bits = 1u64;
    while tried < max_attempts {
        let sqrt_b = F2mElement::from_biguint(&num_bigint::BigUint::from(candidate_bits), big_n);
        candidate_bits = candidate_bits.wrapping_add(1);
        if candidate_bits == 0 {
            candidate_bits = 1;
        }
        tried += 1;
        if sqrt_b.is_zero() {
            continue;
        }
        // Must lie in F_{2^{m_target·l}}.
        if sub_l < big_n && !tower.is_in_subfield(&sqrt_b, sub_l) {
            continue;
        }
        // Must NOT lie in a strict subfield of F_{2^{m_target·l}}.
        // Equivalent: the σ-orbit must have length exactly `m_target`.
        let orbit = tower.frobenius_orbit_length(&sqrt_b);
        if orbit != m_target {
            continue;
        }
        // Also rule out: span has lower-than-expected F_2-rank
        // (rare but possible if there's an F_2-linear relation).
        let mut orbit_elts = Vec::with_capacity(m_target as usize);
        let mut cur = sqrt_b.clone();
        for _ in 0..m_target {
            orbit_elts.push(cur.clone());
            cur = tower.frobenius(&cur);
        }
        let rank = f2_dim(&orbit_elts, big_n);
        if rank != m_target {
            continue;
        }
        // Recover b = (√b)².
        let b = sqrt_b.square(big_irr);
        if b.is_zero() {
            continue;
        }
        // Build the audit row for the target factorisation.
        let m = magic_number_full(&tower, &a, &b);
        let (genus, type_i) = ghs_genus_with_type(&tower, &b);
        debug_assert_eq!(m, m_target);

        return Some(TrapdoorCurve {
            big_n,
            n,
            l,
            big_irr: big_irr.clone(),
            a,
            b,
            sqrt_b,
            row: DescentRow {
                n,
                l,
                magic_m: m,
                genus,
                type_i,
            },
            full_audit: None,
        });
    }
    None
}

/// The output of [`construct_trapdoor_curve`].
#[derive(Clone, Debug)]
pub struct TrapdoorCurve {
    pub big_n: u32,
    pub n: u32,
    pub l: u32,
    pub big_irr: IrreduciblePoly,
    /// Curve parameter `a` (we fix `a = 1` in the constructor).
    pub a: F2mElement,
    /// Curve parameter `b`.
    pub b: F2mElement,
    /// Precomputed `√b` — useful in the descent construction.
    pub sqrt_b: F2mElement,
    /// The target factorisation row (the trapdoor).
    pub row: DescentRow,
    /// Full auditor table, optionally populated by callers who want
    /// to verify the trapdoor is unique (the target row has the
    /// minimum genus across all factorisations).
    pub full_audit: Option<Vec<DescentRow>>,
}

impl TrapdoorCurve {
    /// Run [`audit_curve`] over every factorisation and attach the
    /// result.  Verifies the trapdoor is **non-trivial**: the target
    /// factorisation is the minimum-genus row.
    pub fn populate_audit(&mut self) {
        let rows = audit_curve(self.big_n, &self.big_irr, &self.a, &self.b);
        self.full_audit = Some(rows);
    }

    /// `true` iff the target factorisation `(self.n, self.l)` has
    /// the **smallest genus** across all factorisations — i.e., the
    /// trapdoor owner's descent really is the optimal attack.
    pub fn trapdoor_is_optimal(&self) -> Option<bool> {
        self.full_audit.as_ref().map(|rows| {
            // Rows are sorted ascending by genus; check the *first*
            // row matches our target factorisation.
            rows.first()
                .map(|r| r.n == self.n && r.l == self.l)
                .unwrap_or(false)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::IrreduciblePoly;
    use num_bigint::BigUint;

    /// `F_{2^6}` with `z⁶ + z + 1` (a known irreducible).
    fn f64_irr() -> IrreduciblePoly {
        IrreduciblePoly {
            degree: 6,
            low_terms: vec![0, 1],
        }
    }

    /// `F_{2^8}` with the AES irreducible `z⁸ + z⁴ + z³ + z + 1`.
    fn f256_irr() -> IrreduciblePoly {
        IrreduciblePoly::deg_8()
    }

    #[test]
    fn sqrt_squares_to_input() {
        let irr = f64_irr();
        let tower = FieldTower::new(6, 3, 2, irr.clone());
        for v in 1..64u64 {
            let x = F2mElement::from_biguint(&BigUint::from(v), 6);
            let s = tower.sqrt(&x);
            assert_eq!(s.square(&irr), x);
        }
    }

    #[test]
    fn subfield_detection() {
        let irr = f64_irr();
        let tower = FieldTower::new(6, 3, 2, irr.clone());
        // F_{2^1} ⊂ F_{2^2} ⊂ F_{2^6}.  In particular 0 and 1 are
        // in F_{2^1}.
        let one = F2mElement::one(6);
        assert!(tower.is_in_subfield(&one, 1));
        assert!(tower.is_in_subfield(&one, 2));
        // Element z = (bit 1) — has degree 1, in F_{2^6} but check
        // if F_{2^2}-membership is sensible.  z² + z + 1 = 0 in F_4
        // iff its Frobenius square fixes it.
        let z = F2mElement::from_bit_positions(&[1], 6);
        // z ∈ F_{2^2}  ⇔  z^4 = z.  Let's compute z^4 = z² · z² and
        // see if it equals z.  Almost certainly not for this z.
        let in_f4 = tower.is_in_subfield(&z, 2);
        // True or false — both are valid outcomes; just exercising
        // the API.  We expect z ∉ F_4 since z has degree 1 < 6 and
        // is not algebraic over F_4 in a way that lands it back.
        let _ = in_f4;
    }

    #[test]
    fn frobenius_orbit_lengths_divide_n() {
        let irr = f64_irr();
        let tower = FieldTower::new(6, 3, 2, irr.clone());
        // For each nonzero element of F_{2^6}, the σ-orbit length
        // must divide n = 3 — so it's 1 or 3.
        for v in 1u64..64 {
            let x = F2mElement::from_biguint(&BigUint::from(v), 6);
            let orbit = tower.frobenius_orbit_length(&x);
            assert!(
                orbit == 1 || orbit == 3,
                "orbit {} not 1 or 3 (v = {})",
                orbit,
                v
            );
            // Orbit 1 iff x ∈ F_4 (the fixed field of σ on K).
            let fixed = tower.is_in_subfield(&x, 2);
            assert_eq!(orbit == 1, fixed, "orbit / subfield mismatch for v = {}", v);
        }
    }

    #[test]
    fn f2_dim_singleton() {
        let m = 8;
        let x = F2mElement::from_biguint(&BigUint::from(0b10101u32), m);
        assert_eq!(f2_dim(&[x], m), 1);
        // Zero should give rank 0.
        assert_eq!(f2_dim(&[F2mElement::zero(m)], m), 0);
    }

    #[test]
    fn f2_dim_three_independent_vectors() {
        let m = 8;
        let a = F2mElement::from_biguint(&BigUint::from(0b0001u32), m);
        let b = F2mElement::from_biguint(&BigUint::from(0b0010u32), m);
        let c = F2mElement::from_biguint(&BigUint::from(0b0100u32), m);
        let d = F2mElement::from_biguint(&BigUint::from(0b0111u32), m); // a+b+c
                                                                        // {a, b, c} → rank 3.  {a, b, c, d} → still rank 3.
        assert_eq!(f2_dim(&[a.clone(), b.clone(), c.clone()], m), 3);
        assert_eq!(f2_dim(&[a, b, c, d], m), 3);
    }

    #[test]
    fn magic_number_b_in_subfield_is_one() {
        // If b ∈ F_{2^l}, then √b ∈ F_{2^l} too (char-2 subfields are
        // closed under sqrt), so the σ-orbit of √b is just {√b}.
        // Magic number = 1, genus = 1.
        let irr = f256_irr();
        let tower = FieldTower::new(8, 4, 2, irr.clone());
        // Pick an element of F_4: x s.t. x^4 = x.  Try x = 0x11 ∈ F_{2^8}:
        // probably not in F_4.  Easier: pick b ∈ F_2 ⊂ F_4.  b = 1
        // certainly works.
        let b = F2mElement::one(8);
        let a = F2mElement::zero(8);
        let m = magic_number(&tower, &a, &b);
        assert_eq!(m, 1);
        assert_eq!(ghs_genus(m), BigUint::one());
    }

    #[test]
    fn genus_is_exact_beyond_u32_shift_width() {
        assert_eq!(ghs_genus(32), BigUint::from(1u32) << 31usize);
        assert_eq!(ghs_genus(33), BigUint::from(1u32) << 32usize);
        assert_eq!(ghs_genus(131), BigUint::from(1u32) << 130usize);
        assert_eq!(
            ghs_genus_for_type(131, true),
            (BigUint::from(1u32) << 130usize) - BigUint::from(1u32)
        );
    }

    #[test]
    fn magic_number_constructor_targets_two() {
        // Try to construct a curve with magic number 2, descent
        // genus 2, on the (8, 4, 2) tower.
        let irr = f256_irr();
        let constructed = construct_trapdoor_curve(8, 4, 2, 2, &irr, 256);
        let tc = constructed.expect("found b with m=2");
        assert_eq!(tc.row.magic_m, 2);
        // Genus should be 2 (or 1 if type-I).  Either way ≤ 2.
        assert!(tc.row.genus <= BigUint::from(2u32));
    }

    #[test]
    fn audit_curve_includes_all_factorisations() {
        let irr = f256_irr();
        // Quick: pick b = 1 (so all factorisations give magic = 1).
        let a = F2mElement::zero(8);
        let b = F2mElement::one(8);
        let rows = audit_curve(8, &irr, &a, &b);
        // 8 = 2·4 = 4·2.  Both should appear.
        let factorisations: Vec<(u32, u32)> = rows.iter().map(|r| (r.n, r.l)).collect();
        assert!(factorisations.contains(&(4, 2)));
        assert!(factorisations.contains(&(2, 4)));
        // All magic numbers should be 1 (b is in F_2).
        for r in &rows {
            assert_eq!(r.magic_m, 1, "{:?}", r);
        }
    }

    /// `F_{2^12}` with `z^12 + z^3 + 1`.
    fn f4096_irr() -> IrreduciblePoly {
        IrreduciblePoly {
            degree: 12,
            low_terms: vec![0, 3],
        }
    }

    /// `Tr_{K/k}(x) = Σ_{i<n} σ^i(x)`.
    fn rel_trace(tower: &FieldTower, x: &F2mElement) -> F2mElement {
        let mut acc = F2mElement::zero(tower.big_n);
        let mut cur = x.clone();
        for _ in 0..tower.n {
            acc = acc.add(&cur);
            cur = tower.frobenius(&cur);
        }
        acc
    }

    /// MMT Theorem 5(i) over `N = 12`: for `n = 3`,
    /// `x^3 − 1 = (x − 1)(x^2 + x + 1)` gives `m = 1` on `k*` and `m = 3`
    /// elsewhere (`m = 2` is not attained); for `n = 2`, `m = 1` on `k*`
    /// and `m = 2` elsewhere.
    #[test]
    fn mmt_theorem5_counts_n12() {
        let irr = f4096_irr();
        for (n, l, expect) in [
            (3u32, 4u32, vec![0u64, 15, 0, 4096 - 16]),
            (2, 6, vec![0, 63, 4096 - 64]),
        ] {
            let tower = FieldTower::new(12, n, l, irr.clone());
            let mut h = vec![0u64; n as usize + 1];
            for v in 1u64..4096 {
                let b = F2mElement::from_biguint(&BigUint::from(v), 12);
                h[magic_number(&tower, &F2mElement::zero(12), &b) as usize] += 1;
            }
            assert_eq!(h, expect, "n = {n}");
        }
    }

    /// The plain orbit dimension is one short exactly on the `b ∉ k` with
    /// `Tr_{K/k}(b) = 0` (`N = 12`, `n = 3`), and those are exactly the
    /// type-I curves.  `magic_number_full` ignores `a`.
    #[test]
    fn plain_orbit_undercounts_on_relative_trace_zero() {
        let irr = f4096_irr();
        let tower = FieldTower::new(12, 3, 4, irr);
        let a = F2mElement::from_biguint(&BigUint::from(0b1011_0110_0101u32), 12);
        let mut extra = 0;
        for v in 1u64..4096 {
            let b = F2mElement::from_biguint(&BigUint::from(v), 12);
            let mmt = magic_number(&tower, &F2mElement::zero(12), &b);
            let plain = plain_orbit_dim(&tower, &b);
            let expect_extra = rel_trace(&tower, &b).is_zero() && !tower.is_in_subfield(&b, 4);
            assert_eq!(mmt, plain + u32::from(expect_extra), "v = {v}");
            assert_eq!(ghs_type_i(&tower, &b), expect_extra, "v = {v}");
            assert_eq!(magic_number_full(&tower, &a, &b), mmt, "v = {v}");
            extra += u32::from(expect_extra);
        }
        // Tr_{K/k} is onto k with kernel of size q^2 = 256, less the
        // 16 elements of k (which have Tr = 3b = b, so only b = 0).
        assert_eq!(extra, 256 - 1);
    }

    /// `n = 5` (`N = 10`, `l = 2`): `x^5 − 1 = (x − 1)(x^4 + x^3 + x^2 + x + 1)`.
    /// Generic `b` has `m = 5`, genus `16`; `b ∉ k` with `Tr_{K/k}(b) = 0` has
    /// `m = 5`, genus `15` (the 15/16 dichotomy of the `F_{2^155}`, `n = 5`
    /// GHS literature).  The old orbit-length rule gave 15 and 8.
    #[test]
    fn genus_type_n5() {
        let irr = IrreduciblePoly {
            degree: 10,
            low_terms: vec![0, 3],
        };
        let tower = FieldTower::new(10, 5, 2, irr);
        let mut seen = std::collections::BTreeMap::new();
        for v in 1u64..1024 {
            let b = F2mElement::from_biguint(&BigUint::from(v), 10);
            let (g, type_i) = ghs_genus_with_type(&tower, &b);
            let m = magic_number(&tower, &F2mElement::zero(10), &b);
            let in_k = tower.is_in_subfield(&b, 2);
            let tr0 = rel_trace(&tower, &b).is_zero();
            let class = if in_k {
                "k"
            } else if tr0 {
                "tr0"
            } else {
                "generic"
            };
            *seen
                .entry((class, m, g.to_string(), type_i))
                .or_insert(0u32) += 1;
        }
        let got: Vec<_> = seen.into_iter().collect();
        assert_eq!(
            got,
            vec![
                (("generic", 5, "16".to_string(), false), 1023 - 3 - 255),
                (("k", 1, "1".to_string(), false), 3),
                (("tr0", 5, "15".to_string(), true), 255),
            ]
        );
    }
}
