//! Shared pieces of the PDP-oracle speedup harness: the rung (field, curve,
//! subspace), the frozen RNG, and lifting an abscissa to curve points.
//! Part of `research/pdp_oracle_speedup_20261011/`.

use crypto_lib::binary_ecc::curve::{point_add, point_neg};
use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_semaev::solve_quadratic_f2m;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;

/// SplitMix64, the generator the protocols freeze.
pub struct SplitMix64(pub u64);

impl SplitMix64 {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

// ── GF(2)[x] on bitmasks, to verify the moduli (as semaev_decomp_bench) ──

fn poly_deg(a: u64) -> i32 {
    63 - a.leading_zeros() as i32
}

fn poly_mulmod(mut a: u64, mut b: u64, f: u64, n: u32) -> u64 {
    let mut r = 0u64;
    while b != 0 {
        if b & 1 == 1 {
            r ^= a;
        }
        b >>= 1;
        a <<= 1;
        if (a >> n) & 1 == 1 {
            a ^= f;
        }
    }
    r
}

fn poly_gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let db = poly_deg(b);
        while a != 0 && poly_deg(a) >= db {
            a ^= b << (poly_deg(a) - db);
        }
        std::mem::swap(&mut a, &mut b);
    }
    a
}

fn prime_factors(mut n: u32) -> Vec<u32> {
    let mut out = Vec::new();
    let mut p = 2;
    while p * p <= n {
        if n.is_multiple_of(p) {
            out.push(p);
            while n.is_multiple_of(p) {
                n /= p;
            }
        }
        p += 1;
    }
    if n > 1 {
        out.push(n);
    }
    out
}

/// Rabin's test.
fn is_irreducible(f: u64, n: u32) -> bool {
    let mut t = 2u64;
    for _ in 0..n {
        t = poly_mulmod(t, t, f, n);
    }
    if t != 2 {
        return false;
    }
    for p in prime_factors(n) {
        let mut t = 2u64;
        for _ in 0..(n / p) {
            t = poly_mulmod(t, t, f, n);
        }
        if poly_gcd(t ^ 2, f) != 1 {
            return false;
        }
    }
    true
}

/// Smallest-weight irreducible of degree `n` by the bench's own search:
/// trinomials first (smallest middle term), then pentanomials.  Returns the
/// low terms in `IrreduciblePoly` form (`[0, k]` or `[0, a, b, c]`).
pub fn find_irreducible(n: u32) -> Vec<u32> {
    for k in 1..n {
        let f = (1u64 << n) | (1u64 << k) | 1;
        if is_irreducible(f, n) {
            return vec![0, k];
        }
    }
    for c in 1..n {
        for b in 1..c {
            for a in 1..b {
                let f = (1u64 << n) | (1u64 << c) | (1u64 << b) | (1u64 << a) | 1;
                if is_irreducible(f, n) {
                    return vec![0, a, b, c];
                }
            }
        }
    }
    panic!("no low-weight irreducible of degree {n}");
}

/// One rung of the ladder: `n = 3ℓ`, the bench's modulus, the one-word
/// field, the Koblitz-family curve `y² + xy = x³ + x² + 1` and the
/// factor-base subspace `V = ⟨1, z, …, z^{ℓ−1}⟩` (words below `2^ℓ`).
pub struct Rung {
    pub l: u32,
    pub n: u32,
    pub irr: IrreduciblePoly,
    pub gf: Gf2,
    pub curve: BinaryCurve,
}

impl Rung {
    pub fn new(l: u32) -> Self {
        let n = 3 * l;
        let irr = IrreduciblePoly {
            degree: n,
            low_terms: find_irreducible(n),
        };
        let gf = Gf2::new(&irr);
        let one = F2mElement::one(n);
        let curve = BinaryCurve {
            m: n,
            irreducible: irr.clone(),
            a: one.clone(),
            b: one,
            generator: BinaryPoint::Infinity,
            order: BigUint::from(1u32),
            cofactor: BigUint::from(1u32),
        };
        Self {
            l,
            n,
            irr,
            gf,
            curve,
        }
    }

    pub fn span(&self) -> u64 {
        1u64 << self.l
    }

    pub fn mask(&self) -> u64 {
        (1u64 << self.n) - 1
    }

    /// The modulus as the integer whose bits are its coefficients.
    pub fn modulus(&self) -> u64 {
        self.irr
            .low_terms
            .iter()
            .fold(1u64 << self.n, |acc, &t| acc | (1u64 << t))
    }

    /// Pairs the baseline loop visits when it refutes: `span(span+1)/2`.
    pub fn pairs_total(&self) -> u64 {
        let s = self.span();
        s * (s + 1) / 2
    }

    /// The `y` with `y² + xy = x³ + x² + 1`, none, one or two of them.
    pub fn lifts(&self, x: u64) -> Vec<u64> {
        let irr = &self.irr;
        let xe = self.gf.to_element(x);
        let x2 = xe.square(irr);
        let x3 = x2.mul(&xe, irr);
        let one = F2mElement::one(self.n);
        let c = x3.add(&x2).add(&one);
        let mut ys: Vec<u64> = solve_quadratic_f2m(&one, &xe, &c, self.n, irr)
            .iter()
            .map(|y| self.gf.from_element(y))
            .collect();
        ys.sort_unstable();
        ys.dedup();
        for &y in &ys {
            debug_assert!(self.curve.is_on_curve(&self.point(x, y)));
        }
        ys
    }

    pub fn point(&self, x: u64, y: u64) -> BinaryPoint {
        BinaryPoint::Affine {
            x: self.gf.to_element(x),
            y: self.gf.to_element(y),
        }
    }

    pub fn add(&self, p: &BinaryPoint, q: &BinaryPoint) -> BinaryPoint {
        point_add(&self.curve, p, q)
    }

    pub fn neg(&self, p: &BinaryPoint) -> BinaryPoint {
        point_neg(p)
    }

    pub fn xy(&self, p: &BinaryPoint) -> Option<(u64, u64)> {
        match p {
            BinaryPoint::Infinity => None,
            BinaryPoint::Affine { x, y } => {
                Some((self.gf.from_element(x), self.gf.from_element(y)))
            }
        }
    }

    /// Re-add a root `(x₁, x₂, x₃)` of `S₄(·, ·, ·, x_R)` in the group: true
    /// iff every `xᵢ` lifts to a rational point and some choice of lifts
    /// (which covers the signs) sums to `±R`.  False when some `xᵢ` has no
    /// rational lift: the root is a root of the polynomial and not a
    /// decomposition over rational points.
    pub fn readd(&self, triple: [u64; 3], r: &BinaryPoint) -> Readd {
        let lifts: Vec<Vec<u64>> = triple.iter().map(|&x| self.lifts(x)).collect();
        if lifts.iter().any(|l| l.is_empty()) {
            return Readd::NotRational;
        }
        let neg_r = self.neg(r);
        for &y1 in &lifts[0] {
            for &y2 in &lifts[1] {
                for &y3 in &lifts[2] {
                    let s = self.add(
                        &self.add(&self.point(triple[0], y1), &self.point(triple[1], y2)),
                        &self.point(triple[2], y3),
                    );
                    if &s == r || s == neg_r {
                        return Readd::Sums;
                    }
                }
            }
        }
        Readd::RationalButNoSum
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Readd {
    /// Rational lifts exist and sum to ±R.
    Sums,
    /// Some abscissa has no rational lift: a root of the polynomial only.
    NotRational,
    /// All lift but no sign choice sums to ±R: would be a bug.
    RationalButNoSum,
}

pub fn sorted3(a: u64, b: u64, c: u64) -> [u64; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}
