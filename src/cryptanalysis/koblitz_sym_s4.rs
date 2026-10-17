//! **m = 3 decomposition through the symmetrised fourth Semaev polynomial.**
//!
//! The S₃ chain ([`crate::cryptanalysis::koblitz_groebner::build_decomposition_system`])
//! solves `x(P₁ + P₂ + P₃) = x(R)` with `3l + n` unknowns: one `n`-bit
//! intermediate point.  Here the three summands are replaced by their
//! elementary symmetric functions `e₁, e₂, e₃`, and `S₄(x₁, x₂, x₃, x_R)` is
//! written in them.  For `y² + xy = x³ + a₂x² + b` over `F_{2^n}`, rewriting
//! `S₄ = Res_T(S₃(x₁,x₂,T), S₃(x₃,x₄,T))` in `e₁, e₂, e₃` by the
//! leading-term algorithm gives (checked against the resultant by
//! `tests::symmetric_form_equals_the_resultant`):
//!
//! ```text
//! S₄ = x⁴b² + e₃x³b + e₃²(x⁴ + b) + e₃³x + e₃⁴ + e₂²x²b + e₂²e₃x³
//!    + e₂⁴x⁴ + e₁²e₃xb + e₁²e₃²x² + e₁⁴b²            (x = x_R)
//! ```
//!
//! `e₁` and `e₂` occur only squared, so after Weil descent every term has
//! Boolean degree ≤ 2: the system is **quadratic**, where the chain is cubic.
//!
//! **The subspaces.**  For a geometric factor base `V = θ·⟨1, g, …, g^{l−1}⟩`
//! (`koblitz_isogeny_cost::geometric_basis`; the monomial basis is θ = 1,
//! g = z) the symmetric functions of three elements of `V` lie in
//! `W_k = θ^k·⟨1, g, …, g^{k(l−1)}⟩` of dimension `k(l−1) + 1`.  So the
//! unknowns are `6l − 3` bits, with no intermediate point.  For a random `V`
//! the `W_k` are not small and this formulation does not apply.
//!
//! **Lifting.**  A root `(e₁, e₂, e₃)` is a decomposition only when
//! `t³ + e₁t² + e₂t + e₃` splits with all three roots (with multiplicity)
//! in `V`; [`split_in`] checks that, and the caller's point lift remains the
//! certificate.

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::koblitz_groebner::{FieldStructure, SymElement, MAX_VARS};
use crate::cryptanalysis::pq_groebner_f2::F2BoolPoly;

/// The subspaces `W₁, W₂, W₃` holding `e₁, e₂, e₃`, as bases.
#[derive(Clone, Debug)]
pub struct SymSpaces {
    pub w: [Vec<F2mElement>; 3],
}

impl SymSpaces {
    /// `W_k = θ^k·⟨1, g, …, g^{k(l−1)}⟩`, or the whole field when
    /// `k(l−1) + 1 ≥ n`.  Independent for prime `n` and `g ∉ F_2`.
    pub fn geometric(
        n: u32,
        l: u32,
        theta: &F2mElement,
        g: &F2mElement,
        irr: &IrreduciblePoly,
    ) -> Self {
        let mut theta_k = F2mElement::one(n);
        let w = std::array::from_fn(|k| {
            theta_k = theta_k.mul(theta, irr);
            let dim = (k as u32 + 1) * (l - 1) + 1;
            if dim >= n {
                (0..n)
                    .map(|i| F2mElement::from_bit_positions(&[i], n))
                    .collect()
            } else {
                let mut out = vec![theta_k.clone()];
                while out.len() < dim as usize {
                    let next = out.last().unwrap().mul(g, irr);
                    out.push(next);
                }
                out
            }
        });
        Self { w }
    }

    /// Unknowns: `Σ dim W_k`.
    pub fn n_vars(&self) -> usize {
        self.w.iter().map(Vec::len).sum()
    }

    fn offsets(&self) -> [usize; 3] {
        [0, self.w[0].len(), self.w[0].len() + self.w[1].len()]
    }

    /// The values `(e₁, e₂, e₃)` a Boolean assignment `root` gives.
    pub fn decode(&self, root: u64, n: u32) -> [F2mElement; 3] {
        let off = self.offsets();
        std::array::from_fn(|k| {
            let mut acc = F2mElement::zero(n);
            for (t, be) in self.w[k].iter().enumerate() {
                if (root >> (off[k] + t)) & 1 == 1 {
                    acc = acc.add(be);
                }
            }
            acc
        })
    }
}

/// The `n` Boolean equations of `S₄(e₁, e₂, e₃, x_R) = 0` over `spaces`.
/// `None` when the unknowns exceed [`MAX_VARS`].
pub fn build_sym_s4_system(
    spaces: &SymSpaces,
    x_r: &F2mElement,
    b: &F2mElement,
    st: &FieldStructure,
    irr: &IrreduciblePoly,
) -> Option<Vec<F2BoolPoly>> {
    let n = st.n;
    let n_vars = spaces.n_vars();
    if n_vars > MAX_VARS {
        return None;
    }
    let off = spaces.offsets();
    let e: [SymElement; 3] =
        std::array::from_fn(|k| SymElement::from_subspace_vars(&spaces.w[k], off[k], n, n_vars));

    // constants in x = x_R and b
    let x = x_r.clone();
    let x2 = x.mul(&x, irr);
    let x3 = x2.mul(&x, irr);
    let x4 = x2.mul(&x2, irr);
    let b2 = b.mul(b, irr);
    let c0 = x4.mul(&b2, irr);
    let c_e3 = x3.mul(b, irr);
    let c_e3sq = x4.add(b);
    let c_e3cu = x.clone();
    let c_e2sq = x2.mul(b, irr);
    let c_e2sq_e3 = x3.clone();
    let c_e2q = x4.clone();
    let c_e1sq_e3 = x.mul(b, irr);
    let c_e1sq_e3sq = x2.clone();
    let c_e1q = b2;

    let e1sq = e[0].square(st);
    let e1q = e1sq.square(st);
    let e2sq = e[1].square(st);
    let e2q = e2sq.square(st);
    let e3 = &e[2];
    let e3sq = e3.square(st);
    let e3cu = e3sq.mul(e3, st);
    let e3q = e3sq.square(st);

    let s = SymElement::constant(&c0, n, n_vars)
        .add(&e3.mul_constant(&c_e3, st))
        .add(&e3sq.mul_constant(&c_e3sq, st))
        .add(&e3cu.mul_constant(&c_e3cu, st))
        .add(&e3q)
        .add(&e2sq.mul_constant(&c_e2sq, st))
        .add(&e2sq.mul(e3, st).mul_constant(&c_e2sq_e3, st))
        .add(&e2q.mul_constant(&c_e2q, st))
        .add(&e1sq.mul(e3, st).mul_constant(&c_e1sq_e3, st))
        .add(&e1sq.mul(&e3sq, st).mul_constant(&c_e1sq_e3sq, st))
        .add(&e1q.mul_constant(&c_e1q, st));
    Some(s.coords.into_iter().filter(|p| !p.is_zero()).collect())
}

/// `S₄(e₁, e₂, e₃, x)` evaluated in the field, from the symmetric form.
pub fn s4_symmetric_eval(
    n: u32,
    e: &[F2mElement; 3],
    x: &F2mElement,
    b: &F2mElement,
    irr: &IrreduciblePoly,
) -> F2mElement {
    let m = |a: &F2mElement, c: &F2mElement| a.mul(c, irr);
    let (x2, b2) = (m(x, x), m(b, b));
    let (x3, x4) = (m(&x2, x), m(&x2, &x2));
    let e1sq = m(&e[0], &e[0]);
    let e2sq = m(&e[1], &e[1]);
    let e3sq = m(&e[2], &e[2]);
    [
        m(&x4, &b2),
        m(&e[2], &m(&x3, b)),
        m(&e3sq, &x4.add(b)),
        m(&m(&e3sq, &e[2]), x),
        m(&e3sq, &e3sq),
        m(&e2sq, &m(&x2, b)),
        m(&m(&e2sq, &e[2]), &x3),
        m(&m(&e2sq, &e2sq), &x4),
        m(&m(&e1sq, &e[2]), &m(x, b)),
        m(&m(&e1sq, &e3sq), &x2),
        m(&m(&e1sq, &e1sq), &b2),
    ]
    .iter()
    .fold(F2mElement::zero(n), |acc, t| acc.add(t))
}

/// `S₄(x₁, x₂, x₃, x₄)` as `Res_T(S₃(x₁,x₂,T), S₃(x₃,x₄,T))`, evaluated
/// directly — the reference the symmetric form is tested against.
pub fn s4_resultant_eval(
    xs: &[F2mElement; 4],
    b: &F2mElement,
    irr: &IrreduciblePoly,
) -> F2mElement {
    let m = |a: &F2mElement, c: &F2mElement| a.mul(c, irr);
    let s12 = xs[0].add(&xs[1]);
    let p12 = m(&xs[0], &xs[1]);
    let s34 = xs[2].add(&xs[3]);
    let p34 = m(&xs[2], &xs[3]);
    let (a, p, c) = (m(&s12, &s12), p12.clone(), m(&p12, &p12).add(b));
    let (d, e, f) = (m(&s34, &s34), p34.clone(), m(&p34, &p34).add(b));
    let t1 = m(&a, &f).add(&m(&c, &d));
    m(&t1, &t1).add(&m(&m(&a, &e).add(&m(&p, &d)), &m(&p, &f).add(&m(&c, &e))))
}

/// The roots of `t³ + e₁t² + e₂t + e₃`, with multiplicity, when all three
/// lie in `v_values`; `None` otherwise.
pub fn split_in(
    n: u32,
    e: &[F2mElement; 3],
    v_values: &[F2mElement],
    irr: &IrreduciblePoly,
) -> Option<[F2mElement; 3]> {
    // coefficients low → high: e₃, e₂, e₁, 1
    let mut poly = vec![e[2].clone(), e[1].clone(), e[0].clone(), F2mElement::one(n)];
    let mut roots = Vec::with_capacity(3);
    let eval = |p: &[F2mElement], t: &F2mElement| {
        p.iter()
            .rev()
            .fold(F2mElement::zero(n), |acc, c| acc.mul(t, irr).add(c))
    };
    for v in v_values {
        // divide out v as often as it divides
        while poly.len() > 1 && eval(&poly, v).is_zero() {
            // synthetic division by (t + v)
            let deg = poly.len() - 1;
            let mut q = vec![F2mElement::zero(n); deg];
            q[deg - 1] = poly[deg].clone();
            for i in (0..deg - 1).rev() {
                q[i] = poly[i + 1].add(&q[i + 1].mul(v, irr));
            }
            poly = q;
            roots.push(v.clone());
        }
        if roots.len() == 3 {
            return Some([roots[0].clone(), roots[1].clone(), roots[2].clone()]);
        }
    }
    None
}

/// One-word arithmetic in `F_{2^n}`, `n ≤ 32`.  With PCLMULQDQ the product
/// is one carry-less multiply and the reduction folds the high half by
/// `f − z^n` (a fixed number of folds for a sparse `f`); otherwise a
/// bit-serial loop whose cost grows with the operands' bit lengths.
/// `KOBLITZ_WORDFIELD_SERIAL=1` forces the serial path.
#[derive(Clone, Debug)]
pub struct WordField {
    n: u32,
    high: Vec<u64>,
    clmul: bool,
}

impl WordField {
    pub fn new(n: u32, irr: &IrreduciblePoly) -> Self {
        assert!(n <= 32, "WordField handles n ≤ 32");
        let z = F2mElement::from_bit_positions(&[1], n);
        let mut p = F2mElement::one(n);
        let mut high = Vec::with_capacity(2 * n as usize);
        for k in 0..2 * n {
            if k >= n {
                high.push(p.raw_bits().first().copied().unwrap_or(0));
            }
            p = p.mul(&z, irr);
        }
        let clmul =
            clmul::usable() && std::env::var("KOBLITZ_WORDFIELD_SERIAL").as_deref() != Ok("1");
        Self { n, high, clmul }
    }

    /// Whether products use the carry-less multiply instruction.
    pub fn uses_clmul(&self) -> bool {
        self.clmul
    }

    #[inline]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        if self.clmul {
            // SAFETY: `clmul` is set only when `clmul::usable()` held.
            return unsafe { self.mul_clmul(a, b) };
        }
        self.mul_serial(a, b)
    }

    /// Carry-less product, then fold `hi·z^n = hi·(f − z^n)` until it fits.
    ///
    /// # Safety
    /// PCLMULQDQ must be available.
    #[inline]
    unsafe fn mul_clmul(&self, a: u64, b: u64) -> u64 {
        let f_low = self.high[0];
        let low_mask = (1u64 << self.n) - 1;
        let mut prod = unsafe { clmul::mul(a, b) } as u64;
        loop {
            let hi = prod >> self.n;
            if hi == 0 {
                return prod;
            }
            prod = (prod & low_mask) ^ unsafe { clmul::mul(hi, f_low) } as u64;
        }
    }

    #[inline]
    fn mul_serial(&self, a: u64, b: u64) -> u64 {
        let mut prod = 0u64;
        let mut bb = b;
        let mut i = 0;
        while bb != 0 {
            if bb & 1 == 1 {
                prod ^= a << i;
            }
            bb >>= 1;
            i += 1;
        }
        let low_mask = (1u64 << self.n) - 1;
        let mut out = prod & low_mask;
        let mut hi = prod >> self.n;
        let mut k = 0;
        while hi != 0 {
            if hi & 1 == 1 {
                out ^= self.high[k];
            }
            hi >>= 1;
            k += 1;
        }
        out
    }

    #[inline]
    pub fn sq(&self, a: u64) -> u64 {
        self.mul(a, a)
    }
}

mod clmul {
    #[cfg(target_arch = "x86_64")]
    pub fn usable() -> bool {
        std::arch::is_x86_feature_detected!("pclmulqdq")
    }
    #[cfg(not(target_arch = "x86_64"))]
    pub fn usable() -> bool {
        false
    }

    /// # Safety
    /// PCLMULQDQ must be available.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "pclmulqdq")]
    pub unsafe fn mul(a: u64, b: u64) -> u128 {
        use std::arch::x86_64::{_mm_clmulepi64_si128, _mm_set_epi64x};
        let product =
            _mm_clmulepi64_si128(_mm_set_epi64x(0, a as i64), _mm_set_epi64x(0, b as i64), 0);
        // SAFETY: `__m128i` and `u128` are both 16 plain bytes; lane 0 is
        // the low half on little-endian x86-64.
        unsafe { std::mem::transmute::<std::arch::x86_64::__m128i, u128>(product) }
    }

    /// # Safety
    /// Never called: `usable()` is false off x86-64.
    #[cfg(not(target_arch = "x86_64"))]
    pub unsafe fn mul(_a: u64, _b: u64) -> u128 {
        unreachable!("no carry-less multiply on this target")
    }
}

/// Counters of one [`decompose_by_e3`] call.
#[derive(Clone, Copy, Debug, Default)]
pub struct E3Stats {
    /// `e₃` values visited.
    pub e3_visited: u64,
    /// `(e₁, e₂, e₃)` roots of `S₄` found (before the splitting test).
    pub roots: u64,
    /// Roots whose cubic split with all three roots in `V`.
    pub split: u64,
    /// `F_{2^n}` multiplications and squarings performed, setup, columns,
    /// right-hand sides and splitting tests included (counted at each call).
    pub field_muls: u64,
    /// One-word XORs in the per-`e₃` eliminations and solution sweeps.
    pub word_xors: u64,
}

/// **Exact m = 3 decomposition by enumerating `e₃`.**  With `e₃` fixed,
/// `S₄(e₁, e₂, e₃, x)` is `F_2`-affine in the bits of `e₁ ∈ W₁` and
/// `e₂ ∈ W₂` (they occur only as `e₁², e₁⁴, e₂², e₂⁴`), so each `e₃ ∈ W₃`
/// costs one `n × (dim W₁ + dim W₂)` elimination.  Every root is found;
/// `accept` gets each one whose cubic splits in `V` (values as one-word
/// bitmasks) and stops the search by returning `true`.
///
/// Cost `2^{dim W₃}·poly(n)` — exponent `3l`, the class of brute force over
/// `V³`, not a sub-exponential step.
pub fn decompose_by_e3(
    spaces: &SymSpaces,
    v_values: &[u64],
    x_r: u64,
    b: u64,
    f: &WordField,
    mut accept: impl FnMut([u64; 3]) -> bool,
) -> (bool, E3Stats) {
    let w = |k: usize| -> Vec<u64> {
        spaces.w[k]
            .iter()
            .map(|e| e.raw_bits().first().copied().unwrap_or(0))
            .collect()
    };
    let (w1, w2, w3) = (w(0), w(1), w(2));
    let (d1, d2) = (w1.len(), w2.len());
    let cols = d1 + d2;
    assert!(cols <= 63);
    let x = x_r;
    let mut muls = 0u64;
    let mut xors = 0u64;
    let (x2, b2) = (f.sq(x), f.sq(b));
    let (x3, x4) = (f.mul(x2, x), f.sq(x2));
    let xb = f.mul(x, b);
    let x2b = f.mul(x2, b);
    let x4b2 = f.mul(x4, b2);
    let x3b = f.mul(x3, b);
    let x4_b = x4 ^ b;
    // per-basis-vector squares and fourth powers
    let u2: Vec<u64> = w1.iter().map(|&u| f.sq(u)).collect();
    let u4b2: Vec<u64> = u2.iter().map(|&u| f.mul(f.sq(u), b2)).collect();
    let v2: Vec<u64> = w2.iter().map(|&v| f.sq(v)).collect();
    let v2x2b_v4x4: Vec<u64> = v2
        .iter()
        .map(|&v| f.mul(v, x2b) ^ f.mul(f.sq(v), x4))
        .collect();
    let v2x3: Vec<u64> = v2.iter().map(|&v| f.mul(v, x3)).collect();
    // The setup above: x², b², x³, x⁴, xb, x²b, x⁴b², x³b; per W₁ vector
    // u², u⁴, ·b²; per W₂ vector v², v²·x²b, v⁴, v⁴·x⁴, v²·x³.
    muls += 8 + 3 * d1 as u64 + 5 * d2 as u64;

    let mut stats = E3Stats::default();
    let d3 = w3.len();
    let mut e3 = 0u64;
    let mut column = vec![0u64; cols];
    for step in 0u64..(1u64 << d3) {
        if step > 0 {
            e3 ^= w3[step.trailing_zeros() as usize]; // Gray code
        }
        stats.e3_visited += 1;
        let e3sq = f.sq(e3);
        // constant part S4(0, 0, e3)
        let rhs =
            x4b2 ^ f.mul(e3, x3b) ^ f.mul(e3sq, x4_b) ^ f.mul(f.mul(e3sq, e3), x) ^ f.sq(e3sq);
        // e1 columns: u²(e3·x·b + e3²·x²) + u⁴b²
        let a1 = f.mul(e3, xb) ^ f.mul(e3sq, x2);
        // e3², the right-hand side's five products, a1's two, one per column.
        muls += 1 + 5 + 2 + (d1 + d2) as u64;
        for i in 0..d1 {
            column[i] = f.mul(u2[i], a1) ^ u4b2[i];
        }
        // e2 columns: v²(x²b + e3·x³) + v⁴x⁴
        for j in 0..d2 {
            column[d1 + j] = v2x2b_v4x4[j] ^ f.mul(v2x3[j], e3);
        }
        // solve Σ bit_c · column[c] = rhs over F_2^n
        // elimination on (column-bitmask, value) pairs
        let mut pivots: Vec<(u64, u64)> = Vec::new(); // (value, combo)
        let mut kernel: Vec<u64> = Vec::new();
        for (c, &val) in column.iter().enumerate() {
            let mut v = val;
            let mut combo = 1u64 << c;
            for &(pv, pc) in &pivots {
                if v ^ pv < v {
                    v ^= pv;
                    combo ^= pc;
                    xors += 2;
                }
            }
            if v == 0 {
                kernel.push(combo);
            } else {
                pivots.push((v, combo));
                pivots.sort_unstable_by_key(|a| std::cmp::Reverse(a.0));
            }
        }
        let mut r = rhs;
        let mut sol = 0u64;
        for &(pv, pc) in &pivots {
            if r ^ pv < r {
                r ^= pv;
                sol ^= pc;
                xors += 2;
            }
        }
        if r != 0 {
            continue;
        }
        for kmask in 0u64..(1u64 << kernel.len()) {
            let mut s = sol;
            for (t, &k) in kernel.iter().enumerate() {
                if (kmask >> t) & 1 == 1 {
                    s ^= k;
                }
            }
            // The kernel combination and the two span folds below.
            xors += u64::from(kmask.count_ones() + s.count_ones());
            stats.roots += 1;
            let span = |ws: &[u64], bits: u64| {
                ws.iter()
                    .enumerate()
                    .filter(|(t, _)| (bits >> t) & 1 == 1)
                    .fold(0, |a, (_, &u)| a ^ u)
            };
            let e1 = span(&w1, s & ((1u64 << d1) - 1));
            let e2 = span(&w2, s >> d1);
            if let Some(roots) = split_words_counted(f, [e1, e2, e3], v_values, &mut muls) {
                stats.split += 1;
                if accept(roots) {
                    stats.field_muls = muls;
                    stats.word_xors = xors;
                    return (true, stats);
                }
            }
        }
    }
    stats.field_muls = muls;
    stats.word_xors = xors;
    (false, stats)
}

/// **m = 2 decomposition by one linear solve.**  With `s = x₁ + x₂ ∈ W₁`
/// and `p = x₁x₂ ∈ W₂`, `S₃(x₁, x₂, X) = s²X² + pX + p² + b`, which is
/// `F_2`-affine in the bits of `(s, p)`: squaring is linear and `X` is a
/// constant.  So every pair `{x₁, x₂} ⊂ V` with `S₃(x₁, x₂, X) = 0` comes
/// from one `n × (dim W₁ + dim W₂)` elimination plus a splitting test of
/// `t² + st + p` over `V`.  `accept` gets each pair and stops the search by
/// returning `true`; the second value counts the `(s, p)` roots found.
pub fn decompose_pair_linear(
    w1: &[u64],
    w2: &[u64],
    x_t: u64,
    b: u64,
    f: &WordField,
    v_values: &[u64],
    mut accept: impl FnMut([u64; 2]) -> bool,
) -> (bool, u64) {
    let (d1, d2) = (w1.len(), w2.len());
    assert!(d1 + d2 <= 63);
    let x2 = f.sq(x_t);
    let columns: Vec<u64> = w1
        .iter()
        .map(|&u| f.mul(f.sq(u), x2))
        .chain(w2.iter().map(|&w| f.mul(w, x_t) ^ f.sq(w)))
        .collect();
    let mut pivots: Vec<(u64, u64)> = Vec::new();
    let mut kernel: Vec<u64> = Vec::new();
    for (c, &val) in columns.iter().enumerate() {
        let (mut v, mut combo) = (val, 1u64 << c);
        for &(pv, pc) in &pivots {
            if v ^ pv < v {
                v ^= pv;
                combo ^= pc;
            }
        }
        if v == 0 {
            kernel.push(combo);
        } else {
            pivots.push((v, combo));
            pivots.sort_unstable_by_key(|a| std::cmp::Reverse(a.0));
        }
    }
    let (mut r, mut sol) = (b, 0u64);
    for &(pv, pc) in &pivots {
        if r ^ pv < r {
            r ^= pv;
            sol ^= pc;
        }
    }
    if r != 0 {
        return (false, 0);
    }
    let span = |ws: &[u64], bits: u64| {
        ws.iter()
            .enumerate()
            .filter(|(t, _)| (bits >> t) & 1 == 1)
            .fold(0u64, |a, (_, &u)| a ^ u)
    };
    let mut roots = 0;
    for kmask in 0u64..(1u64 << kernel.len()) {
        let mut sp = sol;
        for (t, &k) in kernel.iter().enumerate() {
            if (kmask >> t) & 1 == 1 {
                sp ^= k;
            }
        }
        roots += 1;
        let sv = span(w1, sp & ((1u64 << d1) - 1));
        let pv = span(w2, sp >> d1);
        // t² + s·t + p splits in V iff one root lies in V (the other is s + t).
        if let Some(&x1) = v_values.iter().find(|&&t| f.sq(t) ^ f.mul(sv, t) ^ pv == 0) {
            if accept([x1, sv ^ x1]) {
                return (true, roots);
            }
        }
    }
    (false, roots)
}

/// [`split_in`] on one-word values.
pub fn split_words(f: &WordField, e: [u64; 3], v_values: &[u64]) -> Option<[u64; 3]> {
    let mut ignored = 0u64;
    split_words_counted(f, e, v_values, &mut ignored)
}

/// [`split_words`], counting its `F_{2^n}` multiplications on `muls`.
pub fn split_words_counted(
    f: &WordField,
    e: [u64; 3],
    v_values: &[u64],
    muls: &mut u64,
) -> Option<[u64; 3]> {
    let mut poly = vec![e[2], e[1], e[0], 1u64];
    let mut roots = [0u64; 3];
    let mut found = 0;
    let eval = |p: &[u64], t: u64, muls: &mut u64| {
        *muls += p.len() as u64;
        p.iter().rev().fold(0u64, |acc, &c| f.mul(acc, t) ^ c)
    };
    for &v in v_values {
        while poly.len() > 1 && eval(&poly, v, muls) == 0 {
            let deg = poly.len() - 1;
            let mut q = vec![0u64; deg];
            q[deg - 1] = poly[deg];
            for i in (0..deg - 1).rev() {
                q[i] = poly[i + 1] ^ f.mul(q[i + 1], v);
                *muls += 1;
            }
            poly = q;
            roots[found] = v;
            found += 1;
        }
        if found == 3 {
            return Some(roots);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::koblitz_groebner::{
        build_decomposition_system, solve_boolean_system, SolveOptions, SolverEngine, SplitRule,
    };
    use crate::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use std::collections::BTreeSet;

    fn el(v: u64, n: u32) -> F2mElement {
        let b: Vec<u32> = (0..n).filter(|i| (v >> i) & 1 == 1).collect();
        F2mElement::from_bit_positions(&b, n)
    }
    fn bits(e: &F2mElement) -> u64 {
        e.raw_bits().first().copied().unwrap_or(0)
    }
    fn elementary(x: &[F2mElement; 3], irr: &IrreduciblePoly) -> [F2mElement; 3] {
        let e1 = x[0].add(&x[1]).add(&x[2]);
        let e2 = x[0]
            .mul(&x[1], irr)
            .add(&x[0].mul(&x[2], irr))
            .add(&x[1].mul(&x[2], irr));
        let e3 = x[0].mul(&x[1], irr).mul(&x[2], irr);
        [e1, e2, e3]
    }
    fn opts() -> SolveOptions {
        SolveOptions {
            engine: SolverEngine::MatrixF4 { max_degree: 4 },
            max_solutions: usize::MAX,
            node_budget: 1_000_000,
            split_rule: SplitRule::default(),
        }
    }

    #[test]
    fn symmetric_form_equals_the_resultant() {
        let n = 17;
        let irr = find_irreducible_sparse(n).unwrap();
        let mut rng = StdRng::seed_from_u64(7);
        let mask = (1u64 << n) - 1;
        for _ in 0..500 {
            let x: [F2mElement; 4] = std::array::from_fn(|_| el(rng.gen::<u64>() & mask, n));
            let b = el(rng.gen::<u64>() & mask, n);
            let e = elementary(&[x[0].clone(), x[1].clone(), x[2].clone()], &irr);
            assert_eq!(
                s4_symmetric_eval(n, &e, &x[3], &b, &irr),
                s4_resultant_eval(&x, &b, &irr)
            );
        }
    }

    /// The linear m = 2 solve returns exactly the brute-force pairs.
    #[test]
    fn pair_linear_matches_brute_force() {
        let (n, l) = (13u32, 4u32);
        let irr = find_irreducible_sparse(n).unwrap();
        let f = WordField::new(n, &irr);
        let mut rng = StdRng::seed_from_u64(5);
        for (theta, g) in [(el(1, n), el(2, n)), (el(0x5a3, n), el(0x1c9, n))] {
            let spaces = SymSpaces::geometric(n, l, &theta, &g, &irr);
            let words = |k: usize| -> Vec<u64> { spaces.w[k].iter().map(bits).collect() };
            let (w1, w2) = (words(0), words(1));
            let v: Vec<u64> = (0..1u64 << l)
                .map(|mask| {
                    (0..l as usize)
                        .filter(|t| (mask >> t) & 1 == 1)
                        .fold(0, |a, t| a ^ w1[t])
                })
                .collect();
            let b = 1u64;
            let mut nonempty = 0;
            for _ in 0..60 {
                // a target that decomposes: a root X of S3(v_i, v_j, X)
                let (i, j) = (rng.gen_range(0..v.len()), rng.gen_range(0..v.len()));
                let (s0, p0) = (v[i] ^ v[j], f.mul(v[i], v[j]));
                let Some(x_t) = (0..1u64 << n)
                    .find(|&x| f.mul(f.sq(s0), f.sq(x)) ^ f.mul(p0, x) ^ f.sq(p0) ^ b == 0)
                else {
                    continue;
                };
                let mut brute = BTreeSet::new();
                for i in 0..v.len() {
                    for j in i..v.len() {
                        let s = v[i] ^ v[j];
                        let p = f.mul(v[i], v[j]);
                        if f.mul(f.sq(s), f.sq(x_t)) ^ f.mul(p, x_t) ^ f.sq(p) ^ b == 0 {
                            brute.insert([v[i].min(v[j]), v[i].max(v[j])]);
                        }
                    }
                }
                let mut lin = BTreeSet::new();
                decompose_pair_linear(&w1, &w2, x_t, b, &f, &v, |[a, c]| {
                    lin.insert([a.min(c), a.max(c)]);
                    false
                });
                assert_eq!(lin, brute);
                nonempty += !brute.is_empty() as usize;
            }
            assert!(nonempty > 0);
        }
    }

    #[test]
    fn word_field_matches_f2m() {
        let n = 23;
        let irr = find_irreducible_sparse(n).unwrap();
        let f = WordField::new(n, &irr);
        let mut rng = StdRng::seed_from_u64(3);
        for _ in 0..2000 {
            let (a, b) = (
                rng.gen::<u64>() & ((1 << n) - 1),
                rng.gen::<u64>() & ((1 << n) - 1),
            );
            let want = bits(&el(a, n).mul(&el(b, n), &irr));
            assert_eq!(f.mul(a, b), want);
            assert_eq!(f.mul_serial(a, b), want);
        }
    }

    #[test]
    fn split_in_recovers_roots_with_multiplicity() {
        let n = 13;
        let irr = find_irreducible_sparse(n).unwrap();
        let v: Vec<F2mElement> = (0..16).map(|i| el(i, n)).collect();
        for triple in [[3u64, 5, 9], [4, 4, 7], [6, 6, 6], [0, 2, 2]] {
            let x = triple.map(|t| el(t, n));
            let e = elementary(&x, &irr);
            let mut got: Vec<u64> = split_in(n, &e, &v, &irr)
                .unwrap()
                .iter()
                .map(bits)
                .collect();
            got.sort();
            let mut want = triple.to_vec();
            want.sort();
            assert_eq!(got, want);
        }
        // a root outside V is refused
        let e = elementary(&[el(3, n), el(5, n), el(1 << 9, n)], &irr);
        assert!(split_in(n, &e, &v, &irr).is_none());
    }

    /// Every multiset {x₁, x₂, x₃} ⊂ V with S₄(x₁, x₂, x₃, x_R) = 0 is found
    /// by the symmetric system, by the S₃ chain, and by brute force — the
    /// same set from all three, over monomial and geometric V.
    #[test]
    fn symmetric_solver_agrees_with_chain_and_brute_force() {
        let (n, l) = (11u32, 3u32);
        let irr = find_irreducible_sparse(n).unwrap();
        let st = FieldStructure::new(n, &irr);
        let b = el(1, n);
        let mut rng = StdRng::seed_from_u64(11);
        for (theta, g) in [(el(1, n), el(2, n)), (el(0x2a5, n), el(0x13b, n))] {
            let mut basis = vec![theta.clone()];
            while basis.len() < l as usize {
                let next = basis.last().unwrap().mul(&g, &irr);
                basis.push(next);
            }
            let v: Vec<F2mElement> = (0..1u64 << l)
                .map(|mask| {
                    (0..l as usize)
                        .filter(|t| (mask >> t) & 1 == 1)
                        .fold(F2mElement::zero(n), |a, t| a.add(&basis[t]))
                })
                .collect();
            let spaces = SymSpaces::geometric(n, l, &theta, &g, &irr);
            let mut nonempty = 0;
            for _ in 0..40 {
                let x_r = el(rng.gen::<u64>() & ((1 << n) - 1), n);
                // brute force over multisets
                let mut brute = BTreeSet::new();
                for i in 0..v.len() {
                    for j in i..v.len() {
                        for k in j..v.len() {
                            let xs = [v[i].clone(), v[j].clone(), v[k].clone(), x_r.clone()];
                            if s4_resultant_eval(&xs, &b, &irr).is_zero() {
                                let mut t = [bits(&v[i]), bits(&v[j]), bits(&v[k])];
                                t.sort();
                                brute.insert(t);
                            }
                        }
                    }
                }
                // symmetric system
                let eqs = build_sym_s4_system(&spaces, &x_r, &b, &st, &irr).unwrap();
                let (roots, stats) = solve_boolean_system(&eqs, spaces.n_vars(), &opts());
                assert!(!stats.exhausted);
                let mut sym = BTreeSet::new();
                for root in roots {
                    let e = spaces.decode(root, n);
                    assert!(s4_symmetric_eval(n, &e, &x_r, &b, &irr).is_zero());
                    if let Some(r) = split_in(n, &e, &v, &irr) {
                        let mut t = r.map(|x| bits(&x));
                        t.sort();
                        sym.insert(t);
                    }
                }
                // S3 chain (ordered triples, normalised)
                let sys = build_decomposition_system(&basis, &x_r, &b, 3, &st).unwrap();
                let (roots, _) = solve_boolean_system(&sys.equations, sys.n_vars, &opts());
                let mut chain = BTreeSet::new();
                for root in roots {
                    let mut t: [u64; 3] =
                        std::array::from_fn(|i| bits(&sys.summand_x(&basis, root, i, n)));
                    t.sort();
                    chain.insert(t);
                }
                // e3 enumeration
                let f = WordField::new(n, &irr);
                let vw: Vec<u64> = v.iter().map(bits).collect();
                let mut by_e3 = BTreeSet::new();
                let (_, st3) = decompose_by_e3(&spaces, &vw, bits(&x_r), bits(&b), &f, |mut r| {
                    r.sort();
                    by_e3.insert(r);
                    false
                });
                assert_eq!(st3.e3_visited, 1 << spaces.w[2].len());
                assert_eq!(by_e3, brute, "e3 enumeration vs brute force");
                assert_eq!(sym, brute, "symmetric vs brute force");
                assert_eq!(chain, brute, "chain vs brute force");
                nonempty += !brute.is_empty() as usize;
            }
            assert!(
                nonempty > 0,
                "no target decomposed; the test checked nothing"
            );
        }
    }
}
