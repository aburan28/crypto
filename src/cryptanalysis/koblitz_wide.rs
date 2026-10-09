//! **Koblitz curves over two-word binary fields, counted, with the tuned
//! signed-Frobenius walk on them**: the group `ecbench` runs the m = 83
//! confidence gate (AGENTS.md §8a) in.
//!
//! The repository's packed binary curves ([`crate::cryptanalysis::koblitz_fast`])
//! hold a field element in one `u64` and stop at degree 62.  This module
//! holds one in a `u128` ([`GfWide`]) and reaches degree 126, with the
//! prime subgroup order, the scalars and the point keys in `u128` too.
//!
//! **The same walk, not a similar one.**  [`signed_frobenius_walk`] is
//! `ic_boundary::rho_walk_with` on `SignedFrobeniusClasses`, ported
//! statement for statement: the seed's draws, the jump table and its
//! coefficient width, the start stride, the jump index and the
//! distinguished-point test from `mix` of the packed point, the normal
//! element the canonicalisation is named in, the least rotation and its
//! tie rule, the sign rule, the look-ahead, the cycle window and its
//! escape, the cap, the verification charge, the counters.  On a curve
//! both can hold (`n ≤ 62`) the two walks take the same steps to the same
//! answer with the same ledger, and the tests check that they do, counter
//! for counter, so a figure the wide walk records at m = 83 is in the
//! unit the narrow walk's figures are in.  Only two things differ by
//! construction: a key's hash folds its high word in (`mix`), which is
//! the identity on a one-word key, and the scalar draws past `2^64` use
//! the RNG's `u128` range.
//!
//! **Scalars.**  The coefficients `a, b` of a walk's point `[a]G + [b]Q`
//! live modulo `r < 2^127`.  A step adds; a canonicalisation multiplies
//! by `±λ^t`, and those products go through a two-word Montgomery form
//! of the multipliers ([`Mont128`]), one reduction each.  The rare
//! products at a collision go through `BigUint`.
//!
//! **An instance is frozen parameters**, never a search:
//! [`WideInstance::explicit`] takes the degree, the modulus, `a`, the
//! subgroup order, the cofactor and a generator, checks them (Rabin's
//! irreducibility test, the Lucas order, primality, `G` on the curve,
//! `[r]G = O`) and derives `λ` with the function the narrow
//! [`crate::cryptanalysis::koblitz_index_calculus::KoblitzCurve`] uses, so
//! a wide instance built from a narrow curve's parameters is that curve.
//!
//! The field is Track B's (`research/ic_tool_program/track-b`, B3),
//! brought into the tree; the point, canon and instance code follows its
//! shape where the narrow code does not dictate one.

use std::collections::{BTreeMap, HashMap, VecDeque};
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::{One, ToPrimitive};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crate::cryptanalysis::curve_id::{self, CurveId};
use crate::cryptanalysis::ecc2k130_guard::is_probable_prime;
use crate::cryptanalysis::gf2_wide::GfWide;
use crate::cryptanalysis::ic_boundary::{generic_floor_ops, GroupOps};
use crate::cryptanalysis::koblitz_index_calculus::frobenius_eigenvalue_q;

/// The widest degree the walk takes: a point packs to `2(x + 1) + s`,
/// which must fit a `u128`, and the subgroup order must stay below
/// `2^127` for the Montgomery reduction's headroom.
pub const MAX_WIDE_DEGREE: u32 = 126;

/// The narrowest: below this the narrow code holds every curve.
pub const MIN_WIDE_DEGREE: u32 = 3;

// ── Points and the curve ───────────────────────────────────────────

/// An affine point of a [`WideCurve`], or `O`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct WidePoint {
    pub x: u128,
    pub y: u128,
    pub infinity: bool,
}

impl WidePoint {
    pub const INFINITY: Self = Self {
        x: 0,
        y: 0,
        infinity: true,
    };

    pub fn affine(x: u128, y: u128) -> Self {
        Self {
            x,
            y,
            infinity: false,
        }
    }

    /// `0` for `O`, otherwise `2(x + 1) + s` with the sign bit `s`
    /// separating `P` from `−P`: `FastPoint::pack` in two words, equal
    /// to it on a one-word point.
    #[inline]
    pub fn pack(self) -> u128 {
        if self.infinity {
            0
        } else {
            let sign = u128::from(self.y > (self.x ^ self.y));
            ((self.x + 1) << 1) | sign
        }
    }
}

/// `ic_boundary::mix` on a one-word key, and on both words once the high
/// word is used, so a walk on a one-word field draws the same jumps and
/// the same distinguished points as the narrow walk.
#[inline]
fn mix(v: u128) -> u64 {
    fn mix64(v: u64) -> u64 {
        let mut h = v ^ 0x2545_F491_4F6C_DD1D;
        h = h.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        h ^= h >> 31;
        h = h.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        h ^= h >> 29;
        h
    }
    let (lo, hi) = (v as u64, (v >> 64) as u64);
    if hi == 0 {
        mix64(lo)
    } else {
        mix64(lo ^ mix64(hi))
    }
}

/// `y² + xy = x³ + ax² + 1` over a [`GfWide`] field, `a ∈ {0, 1}`, in
/// affine coordinates, the formulas of `FastCurve`.
#[derive(Clone, Debug)]
pub struct WideCurve {
    pub field: GfWide,
    pub n: u32,
    pub a: u128,
}

impl WideCurve {
    pub fn is_on_curve(&self, p: WidePoint) -> bool {
        if p.infinity {
            return true;
        }
        let f = &self.field;
        let x2 = f.sqr(p.x);
        f.sqr(p.y) ^ f.mul(p.x, p.y) == f.mul(x2, p.x) ^ f.mul(self.a, x2) ^ 1
    }

    /// `−P = (x, x + y)`.
    #[inline]
    pub fn neg(&self, p: WidePoint) -> WidePoint {
        if p.infinity {
            p
        } else {
            WidePoint::affine(p.x, p.x ^ p.y)
        }
    }

    /// `[2]P`.
    pub fn double(&self, p: WidePoint) -> WidePoint {
        if p.infinity || p.x == 0 {
            return WidePoint::INFINITY;
        }
        let f = &self.field;
        let lambda = p.x ^ f.mul(p.y, f.inv(p.x));
        let x3 = f.sqr(lambda) ^ lambda ^ self.a;
        WidePoint::affine(x3, f.sqr(p.x) ^ f.mul(lambda ^ 1, x3))
    }

    #[inline]
    fn add_with_lambda(&self, p: WidePoint, q: WidePoint, lambda: u128) -> WidePoint {
        let f = &self.field;
        let x3 = f.sqr(lambda) ^ lambda ^ p.x ^ q.x ^ self.a;
        WidePoint::affine(x3, f.mul(lambda, p.x ^ x3) ^ x3 ^ p.y)
    }

    /// `P + Q`.
    pub fn add(&self, p: WidePoint, q: WidePoint) -> WidePoint {
        if p.infinity {
            return q;
        }
        if q.infinity {
            return p;
        }
        if p.x == q.x {
            return if p.y ^ q.y == p.x {
                WidePoint::INFINITY
            } else {
                self.double(p)
            };
        }
        let f = &self.field;
        let lambda = f.mul(p.y ^ q.y, f.inv(p.x ^ q.x));
        self.add_with_lambda(p, q, lambda)
    }

    /// `(x², y²)`: the Frobenius.
    pub fn frobenius(&self, p: WidePoint) -> WidePoint {
        if p.infinity {
            p
        } else {
            WidePoint::affine(self.field.sqr(p.x), self.field.sqr(p.y))
        }
    }

    /// `[k]P` by double-and-add, charged exactly as
    /// `CountedGroup::mul` charges it: one scalar multiplication,
    /// `bits(k) − 1` doublings and `popcount(k) − 1` additions.
    pub fn mul_counted(&self, ops: &mut GroupOps, p: WidePoint, k: u128) -> WidePoint {
        ops.scalar_mults += 1;
        if k == 0 || p.infinity {
            return WidePoint::INFINITY;
        }
        let bits = 128 - k.leading_zeros();
        let mut acc = p;
        for i in (0..bits - 1).rev() {
            ops.doubles += 1;
            acc = self.double(acc);
            if (k >> i) & 1 == 1 {
                ops.adds += 1;
                acc = self.add(acc, p);
            }
        }
        acc
    }

    /// `[k]P` on a ledger nobody reads.
    pub fn mul(&self, p: WidePoint, k: u128) -> WidePoint {
        self.mul_counted(&mut GroupOps::default(), p, k)
    }

    pub fn lift(&self, p: &BinaryPoint) -> WidePoint {
        match p {
            BinaryPoint::Infinity => WidePoint::INFINITY,
            BinaryPoint::Affine { x, y } => {
                WidePoint::affine(self.field.from_element(x), self.field.from_element(y))
            }
        }
    }

    pub fn lower(&self, p: WidePoint) -> BinaryPoint {
        if p.infinity {
            BinaryPoint::Infinity
        } else {
            BinaryPoint::Affine {
                x: self.field.to_element(p.x),
                y: self.field.to_element(p.y),
            }
        }
    }
}

// ── Artin–Schreier roots ───────────────────────────────────────────

/// `ic_boundary::ArtinSchreier` in two words: the echelon form of the
/// `F₂`-linear map `u ↦ u² + u`, so `points_with_x` lifts an abscissa to
/// the *same* root the narrow instance lifts it to (not the half-trace
/// root, which can be the other one), and a public target hashed on a
/// wide instance of a narrow curve is the narrow target.
#[derive(Clone, Debug)]
pub struct WideArtinSchreier {
    n: u32,
    by_lead: Vec<Option<(u128, u128)>>,
}

impl WideArtinSchreier {
    pub fn new(gf: &GfWide) -> Self {
        let n = gf.n;
        let mut by_lead: Vec<Option<(u128, u128)>> = vec![None; 128];
        for i in 0..n {
            let v = 1u128 << i;
            let mut img = gf.sqr(v) ^ v;
            let mut pre = v;
            while img != 0 {
                let lb = 127 - img.leading_zeros() as usize;
                match by_lead[lb] {
                    Some((pi, pp)) => {
                        img ^= pi;
                        pre ^= pp;
                    }
                    None => {
                        by_lead[lb] = Some((img, pre));
                        break;
                    }
                }
            }
        }
        Self { n, by_lead }
    }

    /// A solution of `u² + u = c`, or `None` when `Tr(c) = 1`.
    pub fn solve(&self, mut c: u128) -> Option<u128> {
        let mut u = 0u128;
        while c != 0 {
            let lb = 127 - c.leading_zeros() as usize;
            let (pi, pp) = self.by_lead[lb]?;
            c ^= pi;
            u ^= pp;
        }
        Some(u)
    }

    pub fn degree(&self) -> u32 {
        self.n
    }
}

// ── The canonical key ──────────────────────────────────────────────

/// The normal basis the walk names Frobenius orbits in, and the least
/// rotation: `koblitz_fast::FrobeniusCanon` in two words, from the
/// normal element its randomised search would find (the same candidate
/// sequence, two candidates to an element past 64 bits).
#[derive(Clone, Debug)]
pub struct WideCanon {
    n: u32,
    mask: u128,
    /// `to_normal[i][b]`: normal coordinates of the element whose `i`-th
    /// byte is `b`.
    to_normal: Vec<[u128; 256]>,
    /// `from_normal[i][b]`: the element whose normal coordinates have
    /// `i`-th byte `b`.
    from_normal: Vec<[u128; 256]>,
}

/// Columns → their matrix's inverse, rows as words.
fn invert_f2(columns: &[u128], n: u32) -> Option<Vec<u128>> {
    let mut a: Vec<u128> = (0..n)
        .map(|i| {
            (0..n).fold(0u128, |acc, k| {
                acc | (((columns[k as usize] >> i) & 1) << k)
            })
        })
        .collect();
    let mut inv: Vec<u128> = (0..n).map(|i| 1u128 << i).collect();
    for c in 0..n as usize {
        let pivot = (c..n as usize).find(|&r| (a[r] >> c) & 1 == 1)?;
        a.swap(c, pivot);
        inv.swap(c, pivot);
        for r in 0..n as usize {
            if r != c && (a[r] >> c) & 1 == 1 {
                a[r] ^= a[c];
                inv[r] ^= inv[c];
            }
        }
    }
    Some(inv)
}

fn byte_tables(n: u32, by_bit: &[u128]) -> Vec<[u128; 256]> {
    (0..n.div_ceil(8) as usize)
        .map(|bi| {
            let mut table = [0u128; 256];
            for (b, slot) in table.iter_mut().enumerate() {
                *slot = (0..8)
                    .filter(|t| (b >> t) & 1 == 1 && bi * 8 + t < n as usize)
                    .fold(0u128, |acc, t| acc ^ by_bit[bi * 8 + t]);
            }
            table
        })
        .collect()
}

impl WideCanon {
    /// The basis change, or `None` if no normal element turns up in
    /// 4096 candidates.
    pub fn new(field: &GfWide) -> Option<Self> {
        let n = field.n;
        let mask = field.mask();
        let mut candidate = 2u64;
        let mut next = || {
            candidate = candidate.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ (candidate >> 29);
            candidate
        };
        let mut inverse = None;
        let mut gamma_found = 0u128;
        for _ in 0..4096 {
            let low = next() as u128;
            let gamma = (if n <= 64 {
                low
            } else {
                low | ((next() as u128) << 64)
            }) & mask;
            if gamma == 0 {
                continue;
            }
            let mut column = gamma;
            let columns: Vec<u128> = (0..n)
                .map(|_| {
                    let c = column;
                    column = field.sqr(column);
                    c
                })
                .collect();
            if let Some(inv) = invert_f2(&columns, n) {
                inverse = Some(inv);
                gamma_found = gamma;
                break;
            }
        }
        let inverse = inverse?;
        // Normal coordinates: column j of the inverse per set bit j of x.
        let to_by_bit: Vec<u128> = (0..n)
            .map(|j| {
                (0..n).fold(0u128, |acc, i| {
                    acc | (((inverse[i as usize] >> j) & 1) << i)
                })
            })
            .collect();
        // Back: coordinate i stands for γ^{2^i}.
        let mut from_by_bit = Vec::with_capacity(n as usize);
        let mut column = gamma_found;
        for _ in 0..n {
            from_by_bit.push(column);
            column = field.sqr(column);
        }
        Some(Self {
            n,
            mask,
            to_normal: byte_tables(n, &to_by_bit),
            from_normal: byte_tables(n, &from_by_bit),
        })
    }

    #[inline(always)]
    fn apply(tables: &[[u128; 256]], x: u128) -> u128 {
        let bytes = x.to_le_bytes();
        let mut acc = 0u128;
        for (table, &b) in tables.iter().zip(bytes.iter()) {
            acc ^= table[usize::from(b)];
        }
        acc
    }

    /// The normal coordinates of a field element.
    #[inline(always)]
    pub fn coords(&self, x: u128) -> u128 {
        Self::apply(&self.to_normal, x)
    }

    /// The field element with these normal coordinates.
    #[inline(always)]
    pub fn element(&self, c: u128) -> u128 {
        Self::apply(&self.from_normal, c)
    }

    #[inline(always)]
    fn rotl(&self, v: u128, t: u32) -> u128 {
        if t == 0 {
            v
        } else {
            ((v << t) | (v >> (self.n - t))) & self.mask
        }
    }

    /// `rotl(v, K)` for a constant `0 < K < n`.
    #[inline(always)]
    fn rotl_by<const K: u32>(&self, v: u128) -> u128 {
        let wrapped = if self.n >= 64 + K {
            u128::from(((v >> 64) as u64) >> (self.n - 64 - K))
        } else {
            v >> (self.n - K)
        };
        ((v << K) & self.mask) | wrapped
    }

    /// The canonical name of `x`'s Frobenius orbit with the rotation
    /// that produced it: `FrobeniusCanon::canon_with_shift`.
    #[inline]
    pub fn canon_with_shift(&self, x: u128) -> (u128, u32) {
        self.least_rotation(self.coords(x))
    }

    /// The least of the `n` cyclic rotations of `c`, and the smallest
    /// `t` giving it: `FrobeniusCanon::least_rotation`'s answer.
    ///
    /// The least rotation starts at one of the longest cyclic zero runs
    /// of `c`.  `Z_k` has bit `p` set when bits `p, …, p − k + 1` are all
    /// zero and `Z_{a+b} = Z_a ∧ rotl(Z_b, a)`, so `Z_2, Z_4, Z_8, Z_16`
    /// come by doubling and the longest run below 16 by refinement; a run
    /// of 16 or more and the small degrees grow it one bit at a time.
    /// Only the rotations starting at those runs are compared.
    pub fn least_rotation(&self, c: u128) -> (u128, u32) {
        let n = self.n;
        let zeros = !c & self.mask;
        if zeros == 0 || c == 0 {
            return (c, 0);
        }
        let grow = |mut run: u128| {
            let top = 1u128 << (n - 1);
            loop {
                let next = run & (((run << 1) & self.mask) | u128::from(run & top != 0));
                if next == 0 {
                    break run;
                }
                run = next;
            }
        };
        let run = if n < 17 {
            grow(zeros)
        } else {
            let z2 = zeros & self.rotl_by::<1>(zeros);
            let z4 = z2 & self.rotl_by::<2>(z2);
            let z8 = z4 & self.rotl_by::<4>(z4);
            let z16 = z8 & self.rotl_by::<8>(z8);
            if z16 != 0 {
                grow(z16)
            } else {
                let keep = |cand: u128, cur: u128| if cand != 0 { cand } else { cur };
                let cur = keep(z8, self.mask);
                let cur = keep(z4 & self.rotl_by::<4>(cur), cur);
                let cur = keep(z2 & self.rotl_by::<2>(cur), cur);
                keep(zeros & self.rotl_by::<1>(cur), cur)
            }
        };
        if run & (run - 1) == 0 {
            let t = n - 1 - run.trailing_zeros();
            return (self.rotl(c, t), t);
        }
        let mut best = u128::MAX;
        let mut best_t = 0u32;
        let mut tops = run;
        while tops != 0 {
            let p = tops.trailing_zeros();
            tops &= tops - 1;
            let t = n - 1 - p;
            let v = self.rotl(c, t);
            if v < best || (v == best && t < best_t) {
                best = v;
                best_t = t;
            }
        }
        (best, best_t)
    }

    /// `x^{2^t}`: a rotation of the normal coordinates.
    #[inline]
    pub fn frobenius(&self, x: u128, t: u32) -> u128 {
        if t == 0 {
            x
        } else {
            self.element(self.rotl(self.coords(x), t))
        }
    }

    pub fn degree(&self) -> u32 {
        self.n
    }
}

// ── Scalars modulo r ───────────────────────────────────────────────

/// `a·b` as a 256-bit value, `(low, high)`.
#[inline]
fn mul256(a: u128, b: u128) -> (u128, u128) {
    let (a0, a1) = (a & u128::from(u64::MAX), a >> 64);
    let (b0, b1) = (b & u128::from(u64::MAX), b >> 64);
    let p00 = a0 * b0;
    let p01 = a0 * b1;
    let p10 = a1 * b0;
    let p11 = a1 * b1;
    let (mid, c1) = p01.overflowing_add(p10);
    let mid_lo = mid << 64;
    let mid_hi = (mid >> 64) | (u128::from(c1) << 64);
    let (lo, c2) = p00.overflowing_add(mid_lo);
    (lo, p11 + mid_hi + u128::from(c2))
}

/// Montgomery arithmetic modulo an odd `r < 2^127` with `R = 2^128`:
/// one reduction turns a plain coefficient times a multiplier kept in
/// Montgomery form into the plain product, which is all the walk needs.
#[derive(Clone, Copy, Debug)]
pub struct Mont128 {
    r: u128,
    /// `−r⁻¹ mod 2^128`.
    neg_inv: u128,
    /// `2^256 mod r`.
    r2: u128,
}

impl Mont128 {
    pub fn new(r: u128) -> Self {
        assert!(r & 1 == 1 && r > 1, "an odd modulus above 1");
        assert!(r < (1u128 << 127), "r < 2^127 for the reduction's headroom");
        // Newton: starting from r itself (r·r ≡ 1 mod 8) each step
        // doubles the correct low bits.
        let mut inv = r;
        for _ in 0..7 {
            inv = inv.wrapping_mul(2u128.wrapping_sub(r.wrapping_mul(inv)));
        }
        debug_assert_eq!(r.wrapping_mul(inv), 1);
        let r2 = ((BigUint::one() << 256usize) % BigUint::from(r))
            .to_u128()
            .expect("below r");
        Self {
            r,
            neg_inv: inv.wrapping_neg(),
            r2,
        }
    }

    pub fn modulus(&self) -> u128 {
        self.r
    }

    /// `T·2^{−128} mod r` for `T = lo + hi·2^128 < r·2^128`.
    #[inline]
    fn redc(&self, lo: u128, hi: u128) -> u128 {
        let m = lo.wrapping_mul(self.neg_inv);
        let (mr_lo, mr_hi) = mul256(m, self.r);
        let (_, carry) = lo.overflowing_add(mr_lo);
        let mut t = hi + mr_hi + u128::from(carry);
        if t >= self.r {
            t -= self.r;
        }
        t
    }

    /// `x·2^128 mod r`, for `x < r`.
    pub fn to_mont(&self, x: u128) -> u128 {
        let (lo, hi) = mul256(x, self.r2);
        self.redc(lo, hi)
    }

    /// `a·b mod r` for a plain `a < r` and `b` in Montgomery form.
    #[inline]
    pub fn mul_plain_by_mont(&self, a: u128, b_mont: u128) -> u128 {
        let (lo, hi) = mul256(a, b_mont);
        self.redc(lo, hi)
    }
}

#[inline]
fn addmod(a: u128, b: u128, r: u128) -> u128 {
    (a + b) % r
}

#[inline]
fn submod(a: u128, b: u128, r: u128) -> u128 {
    if a >= b {
        a - b
    } else {
        r - (b - a)
    }
}

fn mulmod_big(a: u128, b: u128, r: u128) -> u128 {
    ((BigUint::from(a) * BigUint::from(b)) % BigUint::from(r))
        .to_u128()
        .expect("below r")
}

/// `a⁻¹ mod r` for a prime `r`, or `None` for `0`.
fn invmod_big(a: u128, r: u128) -> Option<u128> {
    if a == 0 {
        return None;
    }
    let (a, r_big) = (BigUint::from(a), BigUint::from(r));
    let inv = a.modpow(&(&r_big - BigUint::from(2u8)), &r_big);
    ((&a * &inv) % &r_big == BigUint::one()).then(|| inv.to_u128().expect("below r"))
}

// ── Polynomials over GF(2) in a word pair ──────────────────────────

fn poly_degree(f: u128) -> u32 {
    127 - f.leading_zeros()
}

fn poly_gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        while a != 0 && poly_degree(a) >= poly_degree(b) {
            a ^= b << (poly_degree(a) - poly_degree(b));
        }
        std::mem::swap(&mut a, &mut b);
    }
    a
}

fn prime_factors(mut n: u32) -> Vec<u32> {
    let mut out = Vec::new();
    let mut d = 2;
    while d * d <= n {
        if n.is_multiple_of(d) {
            out.push(d);
            while n.is_multiple_of(d) {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        out.push(n);
    }
    out
}

/// Rabin's test (any degree the field holds, 2 to 127): `f` of degree
/// `n` is irreducible over `GF(2)` exactly
/// when `z^{2^n} ≡ z (mod f)` and `gcd(f, z^{2^{n/q}} − z) = 1` for
/// every prime `q | n`.  The powers are squarings in the quotient ring.
pub fn is_irreducible(irr: &IrreduciblePoly) -> bool {
    let n = irr.degree;
    if !(2..=crate::cryptanalysis::gf2_wide::MAX_WIDE_DEGREE).contains(&n)
        || !irr.low_terms.contains(&0)
        || irr.low_terms.iter().any(|&t| t >= n)
    {
        return false;
    }
    let Some(ring) = GfWide::new(irr) else {
        return false;
    };
    let f = irr
        .low_terms
        .iter()
        .fold(1u128 << n, |acc, &t| acc | (1u128 << t));
    let z = 2u128;
    if ring.sqr_k(z, n) != z {
        return false;
    }
    prime_factors(n).into_iter().all(|q| {
        let h = ring.sqr_k(z, n / q) ^ z;
        poly_gcd(f, h) == 1
    })
}

// ── The instance ───────────────────────────────────────────────────

/// A Koblitz curve over a two-word field with its prime-order subgroup,
/// a generator and `λ`, from frozen parameters.
#[derive(Clone, Debug)]
pub struct WideInstance {
    /// The ICV1 slug.
    pub name: String,
    pub n: u32,
    pub a: u8,
    pub irreducible: IrreduciblePoly,
    pub curve: WideCurve,
    pub group_order: u128,
    pub r: u128,
    pub cofactor: u128,
    pub generator: WidePoint,
    /// `λ mod r` with `π(P) = [λ]P` on the subgroup.
    pub lambda: u128,
    pub canon: WideCanon,
    artin: WideArtinSchreier,
}

/// `#E_a(F_{2^n}) = 2^n + 1 − V_n`, `V_0 = 2`, `V_1 = t`,
/// `V_{k+1} = t·V_k − 2·V_{k−1}`, with `t = −1` for `a = 0` and `+1`
/// for `a = 1`.
pub fn koblitz_group_order(a: u8, n: u32) -> u128 {
    let t: i128 = if a == 0 { -1 } else { 1 };
    let (mut v0, mut v1) = (2i128, t);
    for _ in 1..n {
        let next = t * v1 - 2 * v0;
        v0 = v1;
        v1 = next;
    }
    ((1i128 << n) + 1 - v1) as u128
}

impl WideInstance {
    /// `K_a` over `F_2[z]/(irreducible)` with its order-`r` subgroup,
    /// cofactor `h` and generator `(gx, gy)`, checked: the modulus is
    /// irreducible, `3 ≤ n ≤ 126`, `r` is an odd probable prime with
    /// `r·h = #E` by the Lucas sequence, `G` is on the curve, not `O`,
    /// and `[r]G = O`.  `λ` is derived as the narrow `KoblitzCurve`
    /// derives it and checked by `π(G) = [λ]G`.
    pub fn explicit(
        a: u8,
        irreducible: &IrreduciblePoly,
        r: u128,
        cofactor: u128,
        gx: u128,
        gy: u128,
    ) -> Result<Self, String> {
        let n = irreducible.degree;
        if a > 1 {
            return Err(format!("a Koblitz curve has a ∈ {{0, 1}}, not {a}"));
        }
        if !(MIN_WIDE_DEGREE..=MAX_WIDE_DEGREE).contains(&n) {
            return Err(format!(
                "the two-word group holds {MIN_WIDE_DEGREE} ≤ n ≤ {MAX_WIDE_DEGREE}, not {n}"
            ));
        }
        if !is_irreducible(irreducible) {
            return Err(format!(
                "the modulus 0x{:x} of degree {n} is not irreducible over GF(2)",
                curve_id::modulus_integer(irreducible)
            ));
        }
        if r < 5 || r & 1 == 0 || r >= (1u128 << 127) {
            return Err(format!(
                "r = {r} must be an odd prime, at least 5, below 2^127"
            ));
        }
        if !is_probable_prime(&BigUint::from(r)) {
            return Err(format!("r = {r} is not prime"));
        }
        let order = koblitz_group_order(a, n);
        if r.checked_mul(cofactor) != Some(order) {
            return Err(format!(
                "r · cofactor = {r} · {cofactor} is not #E = {order} (Lucas sequence for K_{a} / F_2^{n})"
            ));
        }
        let field = GfWide::new(irreducible).ok_or("no two-word field for this modulus")?;
        let curve = WideCurve {
            field,
            n,
            a: u128::from(a),
        };
        let generator = WidePoint::affine(gx & curve.field.mask(), gy & curve.field.mask());
        if gx != generator.x || gy != generator.y {
            return Err("the generator's coordinates exceed the field".into());
        }
        if !curve.is_on_curve(generator) {
            return Err("the generator is not on the curve".into());
        }
        if !curve.mul(generator, r).infinity {
            return Err("[r]G is not the identity".into());
        }
        // λ as KoblitzCurve::subfield derives it, on the general curve.
        let one = F2mElement::one(n);
        let general = BinaryCurve {
            m: n,
            irreducible: irreducible.clone(),
            a: if a == 1 {
                one.clone()
            } else {
                F2mElement::zero(n)
            },
            b: one,
            generator: curve.lower(generator),
            order: BigUint::from(r),
            cofactor: BigUint::from(cofactor),
        };
        let trace: i64 = if a == 0 { -1 } else { 1 };
        let lambda = frobenius_eigenvalue_q(&general, trace, 2, 1, &BigUint::from(r))
            .and_then(|l| l.to_u128())
            .ok_or("no Frobenius eigenvalue on the subgroup")?;
        if curve.mul(generator, lambda) != curve.frobenius(generator) {
            return Err("π(G) ≠ [λ]G: the derived eigenvalue does not act as the Frobenius".into());
        }
        let canon = WideCanon::new(&curve.field).ok_or("no normal element found")?;
        let artin = WideArtinSchreier::new(&curve.field);
        let mut inst = Self {
            name: String::new(),
            n,
            a,
            irreducible: irreducible.clone(),
            curve,
            group_order: order,
            r,
            cofactor,
            generator,
            lambda,
            canon,
            artin,
        };
        inst.name = inst.curve_id().slug;
        Ok(inst)
    }

    /// The curve's ICV1 identity, from the same facts the narrow
    /// instance hashes (`BinaryInstance::curve_id`).
    pub fn curve_id(&self) -> CurveId {
        curve_id::binary(
            self.n,
            &curve_id::modulus_integer(&self.irreducible),
            &BigUint::from(self.a),
            &BigUint::one(),
            &BigUint::from(self.group_order),
            Some(-7),
        )
        .expect("a checked Koblitz curve is non-singular and inside the Hasse interval")
    }

    /// The modulus as an integer, hex.
    pub fn modulus_hex(&self) -> String {
        format!("0x{:x}", curve_id::modulus_integer(&self.irreducible))
    }

    /// Both points with abscissa `x`, or none: `BinaryInstance::points_with_x`.
    pub fn points_with_x(&self, x: u128) -> Vec<WidePoint> {
        let f = &self.curve.field;
        if x == 0 {
            let y = f.sqr_k(1, self.n - 1);
            return vec![WidePoint::affine(0, y)];
        }
        let inv = f.inv(x);
        let c = x ^ self.curve.a ^ f.mul(1, f.sqr(inv));
        match self.artin.solve(c) {
            Some(u) => {
                let y = f.mul(x, u);
                let p = WidePoint::affine(x, y);
                vec![p, self.curve.neg(p)]
            }
            None => Vec::new(),
        }
    }

    /// `√(πr/2A)` with `A = 2n`: the walk's floor in steps.
    pub fn expected_steps(&self) -> f64 {
        generic_floor_ops(self.r as f64, 2.0 * f64::from(self.n))
    }
}

// ── The tuned signed-Frobenius walk ─────────────────────────────────

/// `SignedFrobeniusClasses` in two words: the representative of a
/// point's signed Frobenius orbit and the multiplier `σλ^t` that reaches
/// it, the latter also in Montgomery form for the step.
struct WideSignedFrobeniusClasses<'a> {
    inst: &'a WideInstance,
    mont: Mont128,
    /// `λ^t mod r`, plain.
    lambda_pow: Vec<u128>,
    /// `λ^t` and `−λ^t` in Montgomery form.
    pos_mont: Vec<u128>,
    neg_mont: Vec<u128>,
}

impl<'a> WideSignedFrobeniusClasses<'a> {
    fn new(inst: &'a WideInstance) -> Self {
        let r = inst.r;
        let mont = Mont128::new(r);
        let mut lambda_pow = Vec::with_capacity(inst.n as usize);
        let mut current = 1u128;
        for _ in 0..inst.n {
            lambda_pow.push(current);
            current = mulmod_big(current, inst.lambda, r);
        }
        let pos_mont = lambda_pow.iter().map(|&l| mont.to_mont(l)).collect();
        let neg_mont = lambda_pow
            .iter()
            .map(|&l| mont.to_mont(if l == 0 { 0 } else { r - l }))
            .collect();
        Self {
            inst,
            mont,
            lambda_pow,
            pos_mont,
            neg_mont,
        }
    }

    /// The representative `σ·φ^t(p)`, `t` and whether `σ = −1`:
    /// `SignedFrobeniusClasses::canon_parts`.
    #[inline]
    fn canon_parts(&self, p: WidePoint) -> (WidePoint, u32, bool) {
        if p.infinity {
            return (p, 0, false);
        }
        let canon = &self.inst.canon;
        let (_, t) = canon.canon_with_shift(p.x);
        let (x, y) = if t == 0 {
            (p.x, p.y)
        } else {
            (canon.frobenius(p.x, t), canon.frobenius(p.y, t))
        };
        // −(x, y) = (x, x + y).
        let negated = x ^ y;
        if negated < y {
            (WidePoint::affine(x, negated), t, true)
        } else {
            (WidePoint::affine(x, y), t, false)
        }
    }

    /// `rho_canon`: the representative with its coefficients, and
    /// whether the canonicalisation was exactly a negation.
    #[inline]
    fn canon(
        &self,
        p: WidePoint,
        a: u128,
        b: u128,
        canonicalisations: &mut u64,
    ) -> (WidePoint, u128, u128, bool) {
        *canonicalisations += 1;
        if p.infinity {
            return (p, a, b, false);
        }
        let (rep, t, negated) = self.canon_parts(p);
        let r = self.inst.r;
        let l = self.lambda_pow[t as usize];
        let mu = if negated {
            if l == 0 {
                0
            } else {
                r - l
            }
        } else {
            l
        };
        if mu == 1 {
            (rep, a, b, false)
        } else {
            let m = if negated {
                self.neg_mont[t as usize]
            } else {
                self.pos_mont[t as usize]
            };
            (
                rep,
                self.mont.mul_plain_by_mont(a, m),
                self.mont.mul_plain_by_mont(b, m),
                mu == r - 1,
            )
        }
    }

    fn method() -> &'static str {
        "signed-Frobenius r-adding walk: normal-basis canonicalisation, look-ahead, short-cycle escape by doubling, distinguished points, stride starts"
    }
}

/// `rho_jumps_for` by the bit length of `r`.
fn jumps_for_bits(bits: u32) -> usize {
    match bits {
        0..=15 => 4,
        16..=22 => 8,
        _ => 16,
    }
}

/// `jump_coefficient_bits`.
fn jump_coefficient_bits(jumps: usize) -> u32 {
    let log = usize::BITS - 1 - jumps.max(2).leading_zeros();
    (3 * log).clamp(6, 12)
}

/// `WalkShape::of(r, RhoWalk::negation())`.
struct Shape {
    jumps: usize,
    dp_bits: u32,
    dp_mask: u64,
    coefficient_bits: u32,
    walk_cap: u64,
}

impl Shape {
    fn of(r: u128) -> Self {
        let bits = 128 - r.leading_zeros();
        let jumps = jumps_for_bits(bits).max(1);
        let dp_bits = (bits / 4).min(40);
        Self {
            jumps,
            dp_bits,
            dp_mask: if dp_bits == 0 {
                0
            } else {
                (1u64 << dp_bits) - 1
            },
            coefficient_bits: jump_coefficient_bits(jumps),
            walk_cap: (1u64 << dp_bits) * 20 + 64,
        }
    }
}

/// `joint_mul`: `[a]P + [b]Q` by one joint double-and-add with `P + Q`
/// supplied; one scalar multiplication on the ledger.
fn joint_mul(
    g: &WideCurve,
    ops: &mut GroupOps,
    p: WidePoint,
    q: WidePoint,
    pq: WidePoint,
    a: u128,
    b: u128,
) -> WidePoint {
    ops.scalar_mults += 1;
    let top = 128 - (a | b).leading_zeros();
    if top == 0 {
        return WidePoint::INFINITY;
    }
    let pick = |i: u32| match ((a >> i) & 1, (b >> i) & 1) {
        (1, 0) => Some(p),
        (0, 1) => Some(q),
        (1, 1) => Some(pq),
        _ => None,
    };
    let mut acc = pick(top - 1).expect("the top bit of a | b is set");
    for i in (0..top - 1).rev() {
        ops.doubles += 1;
        acc = g.double(acc);
        if let Some(t) = pick(i) {
            ops.adds += 1;
            acc = g.add(acc, t);
        }
    }
    acc
}

#[derive(Clone, Copy, Default)]
struct Tally {
    steps: u64,
    look_ahead_adds: u64,
    escapes: u64,
    cycles_2: u64,
    cycles_4: u64,
    cycles_other: u64,
    abandoned: u64,
    capped: u64,
    canonicalisations: u64,
}

impl Tally {
    fn walk_ops(&self) -> u64 {
        self.steps + self.look_ahead_adds + self.escapes
    }
}

enum WalkEnd {
    Distinguished(WidePoint, u128, u128),
    Lost,
}

/// `TunedWalker` on the wide classes.
struct Walker<'a> {
    g: &'a WideCurve,
    classes: &'a WideSignedFrobeniusClasses<'a>,
    r: u128,
    jumps: Vec<(WidePoint, u128, u128)>,
    dp_mask: u64,
    walk_cap: u64,
    window: usize,
    recent: VecDeque<(u128, WidePoint, u128, u128)>,
}

impl Walker<'_> {
    #[inline]
    fn index(&self, x: &WidePoint) -> usize {
        ((mix(x.pack()) >> 24) % self.jumps.len() as u64) as usize
    }

    fn walk(
        &mut self,
        ops: &mut GroupOps,
        tally: &mut Tally,
        start: WidePoint,
        a: u128,
        b: u128,
        budget: u64,
    ) -> WalkEnd {
        let (g, r) = (self.g, self.r);
        let jump_count = self.jumps.len();
        let (mut x, mut a, mut b, _) =
            self.classes
                .canon(start, a, b, &mut tally.canonicalisations);
        self.recent.clear();
        let mut last_escape: Option<u128> = None;
        let mut len = 0u64;
        loop {
            // One step, with the look-ahead.
            let j0 = self.index(&x);
            let mut tried = 0usize;
            let (y, ya, yb) = loop {
                let j = (j0 + tried) % jump_count;
                let (m, aj, bj) = self.jumps[j];
                ops.adds += 1;
                let sum = g.add(x, m);
                let (y, ya, yb, flipped) = self.classes.canon(
                    sum,
                    addmod(a, aj, r),
                    addmod(b, bj, r),
                    &mut tally.canonicalisations,
                );
                tried += 1;
                if flipped && tried < jump_count && !y.infinity && self.index(&y) == j {
                    tally.look_ahead_adds += 1;
                    continue;
                }
                break (y, ya, yb);
            };
            x = y;
            a = ya;
            b = yb;
            tally.steps += 1;
            len += 1;
            if x.infinity {
                return WalkEnd::Lost;
            }
            if self.window > 0 {
                let key = x.pack();
                if let Some(pos) = self.recent.iter().position(|e| e.0 == key) {
                    let c = self.recent.len() - pos;
                    match c {
                        2 => tally.cycles_2 += 1,
                        4 => tally.cycles_4 += 1,
                        _ => tally.cycles_other += 1,
                    }
                    let (min_key, mp, ma, mb) = self
                        .recent
                        .iter()
                        .skip(pos)
                        .copied()
                        .min_by_key(|e| e.0)
                        .expect("a cycle has members");
                    if last_escape == Some(min_key) {
                        tally.abandoned += 1;
                        return WalkEnd::Lost;
                    }
                    last_escape = Some(min_key);
                    ops.doubles += 1;
                    let d2 = g.double(mp);
                    tally.escapes += 1;
                    let (e, ea, eb, _) = self.classes.canon(
                        d2,
                        addmod(ma, ma, r),
                        addmod(mb, mb, r),
                        &mut tally.canonicalisations,
                    );
                    x = e;
                    a = ea;
                    b = eb;
                    self.recent.clear();
                    if x.infinity {
                        return WalkEnd::Lost;
                    }
                }
                self.recent.push_back((x.pack(), x, a, b));
                if self.recent.len() > self.window {
                    self.recent.pop_front();
                }
            }
            if mix(x.pack()) & self.dp_mask == 0 {
                return WalkEnd::Distinguished(x, a, b);
            }
            if len > self.walk_cap {
                tally.capped += 1;
                return WalkEnd::Lost;
            }
            if tally.walk_ops() >= budget {
                return WalkEnd::Lost;
            }
        }
    }
}

/// A counted wide rho run: `ic_boundary::RhoResult` with a two-word
/// answer.
#[derive(Clone, Debug)]
pub struct WideRhoResult {
    pub method: String,
    pub automorphisms: u32,
    pub group_ops: GroupOps,
    pub steps: u64,
    pub walks: u64,
    pub distinguished_points: u64,
    pub gae: f64,
    pub s: f64,
    pub s_walk: f64,
    pub expected_steps: f64,
    pub steps_over_expected: f64,
    pub wall_ns: u64,
    pub recovered: Option<u128>,
    pub verified: bool,
    pub counters: BTreeMap<String, u64>,
}

/// `ic_boundary::signed_frobenius_rho_tuned` on a wide instance: the
/// tuned walk of `rho_walk_with` with `RhoWalk::negation()`'s shape on
/// the signed-Frobenius classes (`A = 2n`), one target, every operation
/// charged, until the logarithm is found and checked or `max_steps`
/// walk operations are spent.
pub fn signed_frobenius_walk(
    inst: &WideInstance,
    target: WidePoint,
    seed: u64,
    max_steps: u64,
) -> WideRhoResult {
    const CYCLE_WINDOW: usize = 16;
    let start = Instant::now();
    let g = &inst.curve;
    let r = inst.r;
    let classes = WideSignedFrobeniusClasses::new(inst);
    let mut ops = GroupOps::default();
    let mut rng = StdRng::seed_from_u64(seed);
    let shape = Shape::of(r);
    let short = r.min(1u128 << shape.coefficient_bits).max(2) as u64;
    let generator = inst.generator;
    // A draw in [1, r): the narrow range below 2^64, so the seeds draw
    // as the narrow walk draws them there.
    let draw = |rng: &mut StdRng| -> u128 {
        match u64::try_from(r) {
            Ok(r64) => u128::from(rng.gen_range(1..r64.max(2))),
            Err(_) => rng.gen_range(1..r),
        }
    };

    // Set-up: G + Q once, the table, the start stride.
    ops.adds += 1;
    let gq = g.add(generator, target);
    let jumps: Vec<(WidePoint, u128, u128)> = (0..shape.jumps)
        .map(|_| {
            let a = u128::from(rng.gen_range(1..short));
            let b = u128::from(rng.gen_range(1..short));
            (
                joint_mul(g, &mut ops, generator, target, gq, a, b),
                a % r,
                b % r,
            )
        })
        .collect();
    let (mut ta, mut tb) = (draw(&mut rng), draw(&mut rng));
    let mut stride = joint_mul(g, &mut ops, generator, target, gq, ta, tb);
    while stride.infinity {
        ta = draw(&mut rng);
        tb = draw(&mut rng);
        stride = joint_mul(g, &mut ops, generator, target, gq, ta, tb);
    }
    let setup = ops;

    let mut walker = Walker {
        g,
        classes: &classes,
        r,
        jumps,
        dp_mask: shape.dp_mask,
        walk_cap: shape.walk_cap,
        window: CYCLE_WINDOW,
        recent: VecDeque::with_capacity(CYCLE_WINDOW + 1),
    };
    let mut table: HashMap<u128, (u128, u128)> = HashMap::with_capacity(1024);
    let automorphisms = 2 * inst.n;
    let expected = generic_floor_ops(r as f64, f64::from(automorphisms));
    let mut tally = Tally::default();
    let (mut walks, mut dps, mut start_adds) = (0u64, 0u64, 0u64);
    let (mut useless, mut verification_ops) = (0u64, 0u64);
    let mut recovered = None;
    let (mut sx, mut sa, mut sb) = (stride, ta % r, tb % r);

    while tally.walk_ops() < max_steps && walks < max_steps {
        walks += 1;
        if walks > 1 {
            ops.adds += 1;
            sx = g.add(sx, stride);
            sa = addmod(sa, ta, r);
            sb = addmod(sb, tb, r);
            start_adds += 1;
        }
        if sx.infinity {
            continue;
        }
        let WalkEnd::Distinguished(x, a, b) =
            walker.walk(&mut ops, &mut tally, sx, sa, sb, max_steps)
        else {
            continue;
        };
        dps += 1;
        let key = x.pack();
        match table.get(&key) {
            Some(&(a2, b2)) => {
                if b2 != b {
                    // a + b d = a2 + b2 d  ⇒  d = (a − a2)/(b2 − b)
                    let num = submod(a, a2, r);
                    let den = submod(b2, b, r);
                    if let Some(inv) = invmod_big(den, r) {
                        let d = mulmod_big(num, inv, r);
                        let before = ops.gae();
                        let check = g.mul_counted(&mut ops, generator, d);
                        verification_ops += (ops.gae() - before) as u64;
                        if check == target {
                            recovered = Some(d);
                            break;
                        }
                    }
                }
                useless += 1;
            }
            None => {
                table.insert(key, (a, b));
            }
        }
    }
    let wall_ns = start.elapsed().as_nanos() as u64;
    let gae = ops.gae();
    let walked = tally.walk_ops();
    let mut counters = BTreeMap::new();
    for (k, v) in [
        ("jumps", shape.jumps as u64),
        ("jump_coefficient_bits", u64::from(shape.coefficient_bits)),
        ("distinguished_point_bits", u64::from(shape.dp_bits)),
        ("setup_additions", setup.adds),
        ("setup_doublings", setup.doubles),
        ("start_additions", start_adds),
        ("walk_operations", walked),
        ("verification_operations", verification_ops),
        ("look_ahead_additions", tally.look_ahead_adds),
        ("cycle_escape_doublings", tally.escapes),
        ("cycles_length_2", tally.cycles_2),
        ("cycles_length_4", tally.cycles_4),
        ("cycles_other_length", tally.cycles_other),
        ("walks_abandoned_in_a_cycle", tally.abandoned),
        ("walks_capped", tally.capped),
        ("useless_collisions", useless),
        ("canonicalisations_uncharged", tally.canonicalisations),
    ] {
        counters.insert(k.to_string(), v);
    }
    WideRhoResult {
        method: WideSignedFrobeniusClasses::method().into(),
        automorphisms,
        group_ops: ops,
        steps: tally.steps,
        walks,
        distinguished_points: dps,
        gae,
        s: gae / (r as f64).sqrt(),
        s_walk: walked as f64 / (r as f64).sqrt(),
        expected_steps: expected,
        steps_over_expected: walked as f64 / expected,
        wall_ns,
        verified: recovered.is_some(),
        recovered,
        counters,
    }
}

/// `rho_cap` for a two-word subgroup: `multiple` times the `A = 1`
/// floor, plus the margin.
pub fn wide_rho_cap(r: u128, multiple: f64) -> u64 {
    (generic_floor_ops(r as f64, 1.0) * multiple) as u64 + 4096
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::{
        koblitz_instance, signed_frobenius_rho_tuned, BinaryInstance,
    };
    use crate::cryptanalysis::koblitz_fast::FrobeniusCanon;

    /// The m = 83 gate's frozen parameters
    /// (`research/ic_tool_program/conformance/v2/params/gate-m83-T001.json`).
    pub(crate) fn gate_m83() -> WideInstance {
        WideInstance::explicit(
            0,
            &IrreduciblePoly {
                degree: 83,
                low_terms: vec![0, 1, 2, 45],
            },
            2_417_851_639_230_796_216_685_689,
            4,
            0x477f_7710_3dfa_d598_5080_0,
            0x2fa5_e737_d542_c4e4_fd5c_3,
        )
        .expect("the gate curve builds")
    }

    fn wide_of(inst: &BinaryInstance) -> WideInstance {
        WideInstance::explicit(
            inst.a as u8,
            &inst.irreducible,
            u128::from(inst.r),
            u128::from(inst.cofactor),
            u128::from(inst.generator.x),
            u128::from(inst.generator.y),
        )
        .expect("a narrow curve's parameters build a wide instance")
    }

    /// Every narrow Koblitz instance up to `max_n` the constructor
    /// builds, with its wide twin.
    fn pairs(max_n: u32) -> Vec<(BinaryInstance, WideInstance)> {
        let mut out = Vec::new();
        for n in [13u32, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61] {
            if n > max_n {
                break;
            }
            for a in [0u8, 1] {
                if let Some(inst) = koblitz_instance(a, n) {
                    let wide = wide_of(&inst);
                    out.push((inst, wide));
                }
            }
        }
        assert!(out.len() >= 6, "{} pairs", out.len());
        out
    }

    fn xorshift(seed: &mut u64) -> u64 {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        *seed
    }

    #[test]
    fn mont128_matches_biguint() {
        let mut seed = 0x1234_5678_9abc_def1u64;
        for r in [
            2_417_851_639_230_796_216_685_689u128,
            (1u128 << 126) + 1 + 2 * 15, // odd, near the top
            65587,
            0xffff_ffff_ffff_ffc5, // a 64-bit prime
            (1u128 << 100) + 277,
        ] {
            let m = Mont128::new(r);
            for _ in 0..200 {
                let a =
                    ((u128::from(xorshift(&mut seed)) << 64) | u128::from(xorshift(&mut seed))) % r;
                let b =
                    ((u128::from(xorshift(&mut seed)) << 64) | u128::from(xorshift(&mut seed))) % r;
                assert_eq!(
                    m.mul_plain_by_mont(a, m.to_mont(b)),
                    mulmod_big(a, b, r),
                    "r = {r}"
                );
            }
            assert_eq!(m.mul_plain_by_mont(1, m.to_mont(1)), 1 % r);
        }
    }

    #[test]
    fn wide_field_matches_the_narrow_field() {
        use crate::cryptanalysis::semaev_decomp::Gf2;
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        for (inst, wide) in pairs(61) {
            let narrow = Gf2::new(&inst.irreducible);
            let f = &wide.curve.field;
            for _ in 0..200 {
                let a = xorshift(&mut seed) & inst.gf.mask;
                let b = xorshift(&mut seed) & inst.gf.mask;
                assert_eq!(u128::from(narrow.mul(a, b)), f.mul(a.into(), b.into()));
                assert_eq!(u128::from(narrow.sqr(a)), f.sqr(a.into()));
                assert_eq!(u128::from(narrow.inv(a)), f.inv(a.into()));
            }
        }
    }

    #[test]
    fn wide_canon_matches_frobenius_canon() {
        let mut seed = 0x0bad_5eed_1234_5678u64;
        for (inst, wide) in pairs(61) {
            let narrow = FrobeniusCanon::new(&inst.gf, inst.n).unwrap();
            let n = inst.n;
            let mask = inst.gf.mask;
            // Random words, and the periodic patterns where rotations tie.
            let mut words: Vec<u64> = (0..300).map(|_| xorshift(&mut seed) & mask).collect();
            words.extend([
                0,
                1,
                mask,
                0x5555_5555_5555_5555 & mask,
                0xaaaa_aaaa_aaaa_aaaa & mask,
            ]);
            for x in words {
                let (c, t) = narrow.canon_with_shift(x);
                assert_eq!(
                    wide.canon.canon_with_shift(x.into()),
                    (u128::from(c), t),
                    "n = {n}, x = {x:#x}"
                );
            }
        }
    }

    #[test]
    fn wide_instance_of_a_narrow_curve_is_that_curve() {
        let mut seed = 0x5eed_0f_c0ffee_u64;
        for (inst, wide) in pairs(61) {
            let kc = inst.koblitz.as_ref().unwrap();
            assert_eq!(wide.name, inst.curve_id().slug, "slug");
            assert_eq!(wide.curve_id().icv1, inst.curve_id().icv1, "ICV1 string");
            assert_eq!(
                BigUint::from(wide.lambda),
                &kc.lambda % BigUint::from(inst.r),
                "λ mod r"
            );
            assert_eq!(wide.group_order, u128::from(inst.group_order));
            for _ in 0..50 {
                let x = xorshift(&mut seed) & inst.gf.mask;
                let narrow: Vec<(u128, u128)> = inst
                    .points_with_x(x)
                    .into_iter()
                    .map(|p| (p.x.into(), p.y.into()))
                    .collect();
                let got: Vec<(u128, u128)> = wide
                    .points_with_x(x.into())
                    .into_iter()
                    .map(|p| (p.x, p.y))
                    .collect();
                assert_eq!(got, narrow, "points above x = {x:#x} on n = {}", inst.n);
                let k = xorshift(&mut seed) % inst.r;
                let mut scratch = GroupOps::default();
                let p = crate::cryptanalysis::ic_boundary::CountedGroup::mul(
                    &crate::cryptanalysis::ic_boundary::BinaryGroup(&inst.fast),
                    &mut scratch,
                    inst.generator,
                    k,
                );
                let q = wide.curve.mul(wide.generator, k.into());
                assert_eq!((q.x, q.y, q.infinity), (p.x.into(), p.y.into(), p.infinity));
            }
        }
    }

    #[test]
    fn wide_walk_matches_the_narrow_walk_step_for_step() {
        // Several targets and seeds per curve: the same ledger, the same
        // steps, walks and distinguished points, the same counters, the
        // same answer.
        let mut seed = 0x7a1e_5eed_u64;
        let mut compared = 0;
        for (inst, wide) in pairs(41) {
            for round in 0..3u64 {
                let k = 1 + xorshift(&mut seed) % (inst.r - 1);
                let mut scratch = GroupOps::default();
                let target = crate::cryptanalysis::ic_boundary::CountedGroup::mul(
                    &crate::cryptanalysis::ic_boundary::BinaryGroup(&inst.fast),
                    &mut scratch,
                    inst.generator,
                    k,
                );
                let walk_seed = 1000 * u64::from(inst.n) + round;
                let cap = crate::cryptanalysis::ic_boundary::rho_cap(inst.r, 64.0);
                let narrow = signed_frobenius_rho_tuned(&inst, target, walk_seed, cap).unwrap();
                let got = signed_frobenius_walk(
                    &wide,
                    WidePoint::affine(target.x.into(), target.y.into()),
                    walk_seed,
                    cap,
                );
                let tag = format!("K_{} n = {} round {round}", inst.a, inst.n);
                assert_eq!(
                    got.recovered,
                    narrow.recovered.map(u128::from),
                    "{tag}: answer"
                );
                assert_eq!(got.recovered, Some(u128::from(k)), "{tag}: the planted k");
                assert_eq!(got.group_ops, narrow.group_ops, "{tag}: ledger");
                assert_eq!(
                    (got.steps, got.walks, got.distinguished_points),
                    (narrow.steps, narrow.walks, narrow.distinguished_points),
                    "{tag}: steps, walks, distinguished points"
                );
                assert_eq!(got.counters, narrow.counters, "{tag}: counters");
                assert_eq!(got.method, narrow.method, "{tag}: method");
                assert_eq!(got.automorphisms, narrow.automorphisms);
                assert_eq!(got.gae, narrow.gae);
                assert_eq!(got.expected_steps, narrow.expected_steps);
                compared += 1;
            }
        }
        assert!(compared >= 18, "{compared} runs compared");
    }

    #[test]
    fn irreducibility_is_decided() {
        let gate = IrreduciblePoly {
            degree: 83,
            low_terms: vec![0, 1, 2, 45],
        };
        assert!(is_irreducible(&gate));
        // z^83 + 1 has the factor z + 1.
        assert!(!is_irreducible(&IrreduciblePoly {
            degree: 83,
            low_terms: vec![0],
        }));
        // z^127 + z + 1 is the standard trinomial; z^64 + z + 1 is not
        // irreducible (no irreducible trinomial of degree 64 exists).
        assert!(is_irreducible(&IrreduciblePoly {
            degree: 127,
            low_terms: vec![0, 1],
        }));
        assert!(!is_irreducible(&IrreduciblePoly {
            degree: 64,
            low_terms: vec![0, 1],
        }));
        for (_, wide) in pairs(61) {
            assert!(is_irreducible(&wide.irreducible), "n = {}", wide.n);
        }
    }

    #[test]
    fn the_gate_curve_builds_with_its_known_facts() {
        let g = gate_m83();
        assert_eq!(g.name, "icv1-f2m83-tm6151469093347-debefd74");
        assert_eq!(g.group_order, 9_671_406_556_923_184_866_742_756);
        assert_eq!(g.lambda, 254_512_724_090_651_164_922_414);
        assert_eq!(g.modulus_hex(), "0x800000000200000000007");
        assert!(g.curve.is_on_curve(g.generator));
        assert!(g.curve.mul(g.generator, g.r).infinity);
        assert!(!g.curve.mul(g.generator, 2).infinity);
        // The floor in steps, √(πr/4n).
        assert!(
            (g.expected_steps() / 1.513e11 - 1.0).abs() < 0.01,
            "{}",
            g.expected_steps()
        );
        // The parameter checks refuse a wrong cofactor and a wrong generator.
        assert!(
            WideInstance::explicit(0, &g.irreducible, g.r, 2, g.generator.x, g.generator.y)
                .is_err()
        );
        assert!(WideInstance::explicit(
            0,
            &g.irreducible,
            g.r,
            4,
            g.generator.x,
            g.generator.y ^ 1
        )
        .is_err());
    }

    #[test]
    fn a_budgeted_walk_on_the_gate_is_deterministic_and_exhausts() {
        let g = gate_m83();
        // The frozen public target T001.
        let t = WidePoint::affine(0x355f_b5df_7a90_5f16_921e_b, 0x5900_a390_f42d_290f_1bbe);
        assert!(g.curve.is_on_curve(t));
        let a = signed_frobenius_walk(&g, t, 2_293_761, 20_000);
        let b = signed_frobenius_walk(&g, t, 2_293_761, 20_000);
        assert!(a.recovered.is_none());
        assert!(a.counters["walk_operations"] >= 20_000);
        assert_eq!(a.counters["distinguished_point_bits"], 20);
        assert_eq!(a.counters["jumps"], 16);
        assert_eq!(
            (a.group_ops, &a.counters, a.steps),
            (b.group_ops, &b.counters, b.steps)
        );
    }
}
