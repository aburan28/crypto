//! # Index calculus for Koblitz (subfield) curves with Frobenius-invariant factor bases.
//!
//! Implementation of the algorithm of **Galbraith, Granger, Merz and
//! Petit**, *On index calculus algorithms for subfield curves*
//! (SAC 2020, ePrint 2020/1315).  A Koblitz curve
//!
//! ```text
//!     K_a :  y² + x y = x³ + a x² + 1,      a ∈ {0, 1}
//! ```
//!
//! is defined over `F_2` but used over `F_{2^n}`.  Because the curve
//! equation has `F_2`-coefficients, the `2`-power Frobenius
//!
//! ```text
//!     π(x, y) = (x², y²)
//! ```
//!
//! is an endomorphism of `E(F_{2^n})`, and on the large prime-order
//! subgroup `⟨P⟩` it acts as multiplication by a known scalar `λ`
//! (a root of `λ² − t λ + 2 ≡ 0 mod r`, `t = (−1)^{1−a}` the trace of
//! Frobenius over `F_2`).  Pollard rho exploits that with a `√n`
//! speed-up; the GGMP paper shows index calculus gets a **larger**
//! speed-up: `≈ n` in relation collection and `≈ n²` in the linear
//! algebra.
//!
//! ## Where the speed-up comes from
//!
//! Take a factor base `F ⊆ E(F_{2^n})` that is *Frobenius invariant*,
//! `π(F) = F`.  Then `F` splits into `π`-orbits, each of size dividing
//! `n`, and if `P' = π^k(Q)` with `Q` an orbit representative then
//!
//! ```text
//!     log_P P' = λ^k · log_P Q      (mod r).
//! ```
//!
//! So the `|F|` unknowns of a classical index calculus collapse to
//! `≈ |F| / n` unknowns — one per orbit.  Consequences:
//!
//! - **Relation collection:** we need `≈ |F| / n` relations instead of
//!   `≈ |F|`, i.e. a factor-`n` saving.
//! - **Linear algebra:** the relation matrix shrinks by a factor `n` in
//!   *both* dimensions, so a quadratic (sparse-linear-algebra) solve
//!   costs `n²` less.
//! - **Symmetry breaking:** a decomposition `R = P_1 + … + P_m` is
//!   found `m!` times if the summands are searched independently.
//!   Distributing the summands over the shifted factor bases
//!   `F_i = π^{i−1}(F)` — which for an invariant base all equal `F` —
//!   and requiring the tuple to be ordered removes that `m!`.  Each
//!   summand `P_i = π^{k_i}(Q_{j_i})` contributes `λ^{k_i}` to column
//!   `j_i` of the relation, which is exactly the "rewrite the relation
//!   so that all `P'_i ∈ F_1 = F`" step of the paper.
//!
//! ## Frobenius-invariant factor bases from linearised polynomials
//!
//! The paper's construction (their §4, and the "Examples of Frobenius
//! invariant factor bases" slide of the SAC talk):
//!
//! 1. Factor `x^n − 1 = (x − 1) f_1 f_2 ⋯ f_s` in `F_2[x]`.  When
//!    `gcd(n, 2) = 1` the `f_i` are distinct irreducible polynomials,
//!    all of degree `ℓ := ord_n(2)`.
//! 2. For `f_j = Σ_k f_{j,k} x^k` form the **linearised** polynomial
//!    `F_j(X) = Σ_k f_{j,k} X^{2^k}`.  Its root set in `F_{2^n}` is an
//!    `F_2`-subspace `V_j` of dimension `ℓ` (the `F_2[x] ≅ {linearised
//!    polynomials under composition}` isomorphism sends `x^n − 1` to
//!    `X^{2^n} − X`, so `F_j` divides `X^{2^n} − X`), and `V_j` is
//!    closed under squaring because the `f_{j,k}` lie in `F_2`.
//! 3. `F := { P ∈ E(F_{2^n}) : F_j(x(P)) = 0 }` is Frobenius invariant
//!    of size `≈ 2^ℓ`.
//!
//! `V_j` is computed here as the kernel of the `F_2`-linear map
//! `v ↦ F_j(v)` on `F_{2^n}`, which is both simpler and self-checking
//! (the kernel dimension must come out equal to `deg f_j`).
//!
//! ## Honest scope
//!
//! - **Two decomposition oracles, and the algebraic one is not yet the
//!   faster one.**  "Is `R` a sum of `m` factor-base points?" is
//!   answered either by
//!   [`DecompositionStrategy::Groebner`] — Semaev's `S₃` Weil-restricted
//!   to a low-degree Boolean system over the invariant subspace and
//!   solved with matrix-F4 (see
//!   [`crate::cryptanalysis::koblitz_groebner`]), which is the real
//!   algorithm's oracle — or by
//!   [`DecompositionStrategy::Enumerate`], a table-driven search over
//!   ordered tuples costing `|F|^{m−1}` group operations.  They are
//!   cross-checked against each other in the tests and always agree.
//!   At the sizes this module can reach the search is still 10–100×
//!   *faster* in wall-clock terms: `|F|` is a few dozen points, so
//!   `|F|^{m−1}` is nothing, while the Macaulay matrix already has
//!   thousands of columns.  The algebra earns its place on scaling
//!   rather than on these numbers — its cost tracks the degree of
//!   regularity of the system instead of `|F|`, and it *refutes* an
//!   undecomposable target with a certificate (the reduction yields the
//!   constant `1`) instead of merely failing to find one.
//! - **Toy parameters only.**  `n ≤ 24` or so: the factor base is
//!   materialised (`2^ℓ` elements) and the group order is found by
//!   trial division.  This does not threaten sect163k1 or any other
//!   deployed Koblitz curve — as the paper's own conclusion puts it,
//!   index calculus remains *worse* than rho for curves used in
//!   practice.
//!
//! ## References
//!
//! - S. D. Galbraith, R. Granger, S.-P. Merz, C. Petit, *On index
//!   calculus algorithms for subfield curves*, SAC 2020,
//!   ePrint 2020/1315.
//! - J.-M. Couveignes, R. Lercier, *Elliptic periods for finite
//!   fields*, Finite Fields Appl. 15 (2009) — the other source of
//!   Frobenius-invariant factor bases mentioned in the talk.
//! - I. Semaev, *Summation polynomials and the discrete logarithm
//!   problem on elliptic curves*, ePrint 2004/031.
//! - R. Gallant, R. Lambert, S. Vanstone, *Improving the parallelized
//!   Pollard lambda search on anomalous binary curves*, Math. Comp. 69
//!   (2000) — the `√n` rho speed-up we compare against.

use std::collections::HashMap;

use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use crate::binary_ecc::curve::{point_add, point_neg, scalar_mul};
use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crate::cryptanalysis::binary_semaev::solve_artin_schreier;
use crate::cryptanalysis::ec_index_calculus::{gaussian_eliminate_mod_n, sqrt_mod_p};
use crate::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, solve_boolean_system_filtered, FieldStructure, SolveOptions,
    SolveStats, SolverEngine,
};
use crate::cryptanalysis::sat::SolveResult;
use crate::cryptanalysis::semaev_sat::encode_boolean_system;
use crate::utils::mod_inverse;

/// Largest extension degree this module will build a curve for.  The
/// factor base and the point-counting/factoring helpers are all
/// materialised, so this is a deliberate guard rail, not a limit of
/// the mathematics.
pub const MAX_N: u32 = 24;

// ── F_2[x] helpers on `u64` bitmasks ───────────────────────────────
//
// Polynomials of degree < 64 over F_2, bit `i` = coefficient of `x^i`.
// Used for the irreducible-polynomial search and for factoring
// `x^n − 1`; the field arithmetic itself lives in `binary_ecc::f2m`.

/// Degree of a bitmask polynomial, or `None` for the zero polynomial.
fn poly_deg(f: u64) -> Option<u32> {
    if f == 0 {
        None
    } else {
        Some(63 - f.leading_zeros())
    }
}

/// Remainder of `a` modulo `f` in `F_2[x]`.
fn poly_rem(mut a: u64, f: u64) -> u64 {
    let df = match poly_deg(f) {
        Some(d) => d,
        None => return a,
    };
    while let Some(da) = poly_deg(a) {
        if da < df {
            break;
        }
        a ^= f << (da - df);
    }
    a
}

/// `a · b mod f` in `F_2[x]`.
fn poly_mulmod(a: u64, b: u64, f: u64) -> u64 {
    let mut acc = 0u64;
    let mut a = poly_rem(a, f);
    let mut b = b;
    while b != 0 {
        if b & 1 == 1 {
            acc ^= a;
        }
        b >>= 1;
        a <<= 1;
        a = poly_rem(a, f);
    }
    poly_rem(acc, f)
}

/// `gcd(a, b)` in `F_2[x]`, monic by construction (F_2 has one unit).
fn poly_gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = poly_rem(a, b);
        a = b;
        b = r;
    }
    a
}

/// `x^(2^k) mod f` by repeated squaring of the residue.
fn poly_x_pow_2k(k: u32, f: u64) -> u64 {
    let mut acc = poly_rem(0b10, f); // x
    for _ in 0..k {
        acc = poly_mulmod(acc, acc, f);
    }
    acc
}

/// Distinct prime divisors of `n` (trial division; `n` is tiny here).
fn prime_divisors(mut n: u32) -> Vec<u32> {
    let mut out = Vec::new();
    let mut d = 2u32;
    while d * d <= n {
        if n % d == 0 {
            out.push(d);
            while n % d == 0 {
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

/// **Rabin's irreducibility test** for a degree-`d` polynomial of
/// `F_2[x]`: `f` is irreducible iff `x^(2^d) ≡ x (mod f)` and
/// `gcd(x^(2^(d/p)) − x, f) = 1` for every prime `p | d`.
pub fn is_irreducible_f2(f: u64) -> bool {
    let d = match poly_deg(f) {
        Some(d) if d >= 1 => d,
        _ => return false,
    };
    if f & 1 == 0 && d > 1 {
        return false; // divisible by x
    }
    if poly_x_pow_2k(d, f) != poly_rem(0b10, f) {
        return false;
    }
    for p in prime_divisors(d) {
        let t = poly_x_pow_2k(d / p, f) ^ 0b10;
        if poly_deg(poly_gcd(poly_rem(t, f), f)) != Some(0) {
            return false;
        }
    }
    true
}

/// Smallest (as an integer bitmask) irreducible polynomial of degree
/// `n` over `F_2`, which for every `n` in range is a trinomial or
/// pentanomial.  Returned as an [`IrreduciblePoly`].
pub fn find_irreducible(n: u32) -> Option<IrreduciblePoly> {
    if n == 0 || n >= 64 {
        return None;
    }
    let hi = 1u64 << n;
    for low in 0..hi {
        let f = hi | low;
        if is_irreducible_f2(f) {
            let low_terms = (0..n).filter(|i| (low >> i) & 1 == 1).collect();
            return Some(IrreduciblePoly {
                degree: n,
                low_terms,
            });
        }
    }
    None
}

/// Multiplicative order of `2` modulo odd `n` — the degree `ℓ` of every
/// non-trivial irreducible factor of `x^n − 1` over `F_2`.
pub fn order_of_2_mod_n(n: u32) -> Option<u32> {
    if n < 3 || n % 2 == 0 {
        return None;
    }
    let mut acc = 1u64;
    for k in 1..=n {
        acc = (acc * 2) % n as u64;
        if acc == 1 {
            return Some(k);
        }
    }
    None
}

/// **Factor `x^n − 1` over `F_2`** into its non-trivial irreducible
/// factors (i.e. excluding `x − 1`), for odd `n`.
///
/// Every such factor has degree `ℓ = ord_n(2)`, so they are found by
/// scanning the `2^ℓ` monic degree-`ℓ` polynomials for the ones that
/// are irreducible and divide `x^n − 1`.  Returned as bitmasks.
pub fn factor_x_n_minus_1(n: u32) -> Vec<u64> {
    let ell = match order_of_2_mod_n(n) {
        Some(l) if l < 40 => l,
        _ => return Vec::new(),
    };
    let mut out = Vec::new();
    let hi = 1u64 << ell;
    for low in 0..hi {
        let f = hi | low;
        if !is_irreducible_f2(f) {
            continue;
        }
        // f | x^n − 1  ⟺  x^n ≡ 1 (mod f).
        let mut xn = poly_rem(0b10, f);
        let mut acc = 1u64;
        let mut e = n;
        while e > 0 {
            if e & 1 == 1 {
                acc = poly_mulmod(acc, xn, f);
            }
            xn = poly_mulmod(xn, xn, f);
            e >>= 1;
        }
        if acc == 1 {
            out.push(f);
        }
    }
    out
}

// ── Koblitz curve over F_{2^n} ─────────────────────────────────────

/// A Koblitz curve `K_a : y² + xy = x³ + a x² + 1` over `F_2`, taken
/// over the extension `F_{2^n}`, together with everything the index
/// calculus needs: the prime-order subgroup, a generator, and the
/// scalar `λ` by which Frobenius acts on that subgroup.
#[derive(Clone, Debug)]
pub struct KoblitzCurve {
    /// Curve parameter `a ∈ {0, 1}`.
    pub a: u8,
    /// Extension degree `n`.
    pub n: u32,
    /// Underlying binary curve (arithmetic + generator + order).
    pub curve: BinaryCurve,
    /// Trace of Frobenius over `F_2`: `t = (−1)^{1−a}`, i.e. `−1` for
    /// `a = 0` and `+1` for `a = 1`.
    pub trace: i64,
    /// `#E(F_{2^n})`.
    pub group_order: BigUint,
    /// Prime order `r` of the subgroup the DLP lives in.
    pub subgroup_order: BigUint,
    /// `h = #E / r`.
    pub cofactor: BigUint,
    /// The eigenvalue of `π` on `⟨generator⟩`: `π(Q) = [λ]Q`.
    pub lambda: BigUint,
}

/// `#E(K_a / F_{2^n}) = 2^n + 1 − s_n`, where `s_0 = 2`, `s_1 = t` and
/// `s_k = t·s_{k−1} − 2·s_{k−2}` (Koblitz's recurrence for the traces
/// of the Frobenius powers).
pub fn koblitz_point_count(a: u8, n: u32) -> BigUint {
    let t: i128 = if a == 0 { -1 } else { 1 };
    let (mut s_prev, mut s_cur) = (2i128, t);
    for _ in 1..n {
        let next = t * s_cur - 2 * s_prev;
        s_prev = s_cur;
        s_cur = next;
    }
    let count = (1i128 << n) + 1 - s_cur;
    BigUint::from(count as u128)
}

/// Trial-division factorisation into `(prime, exponent)` pairs.  Only
/// ever called on `#E` for toy `n`.
fn factorise(mut v: BigUint) -> Vec<(BigUint, u32)> {
    let mut out: Vec<(BigUint, u32)> = Vec::new();
    let mut d = BigUint::from(2u32);
    while &d * &d <= v {
        let mut e = 0;
        while (&v % &d).is_zero() {
            v /= &d;
            e += 1;
        }
        if e > 0 {
            out.push((d.clone(), e));
        }
        d += BigUint::one();
    }
    if v > BigUint::one() {
        out.push((v, 1));
    }
    out
}

impl KoblitzCurve {
    /// Build `K_a` over `F_{2^n}`: pick a defining irreducible
    /// polynomial, count points, split off the largest prime-order
    /// subgroup, find a generator of it, and determine `λ`.
    ///
    /// Returns `None` if `n` is out of range, if no irreducible
    /// polynomial is found, or if the largest prime factor of `#E` is
    /// too small for the subgroup to be usable (`r ≤ h`, or `r² | #E`).
    pub fn new(a: u8, n: u32) -> Option<Self> {
        if a > 1 || n < 3 || n > MAX_N {
            return None;
        }
        let irreducible = find_irreducible(n)?;
        let a_fe = if a == 0 {
            F2mElement::zero(n)
        } else {
            F2mElement::one(n)
        };
        let b_fe = F2mElement::one(n);

        // Arithmetic-only shell; generator/order are filled in below.
        let mut curve = BinaryCurve {
            m: n,
            irreducible,
            a: a_fe,
            b: b_fe,
            generator: BinaryPoint::Infinity,
            order: BigUint::zero(),
            cofactor: BigUint::one(),
        };

        let group_order = koblitz_point_count(a, n);
        let factors = factorise(group_order.clone());
        let (r, e) = factors.last()?.clone();
        if e != 1 {
            return None;
        }
        let cofactor = &group_order / &r;
        if r <= cofactor {
            return None;
        }

        // A generator of the order-r subgroup: kill the cofactor on
        // successive curve points until the result is non-trivial.
        let mut generator = BinaryPoint::Infinity;
        for raw in 0..(1u64 << n) {
            let x = F2mElement::from_biguint(&BigUint::from(raw), n);
            let pts = points_with_x(&curve, &x);
            let mut found = false;
            for p in pts {
                let cand = scalar_mul(&curve, &p, &cofactor);
                if cand != BinaryPoint::Infinity {
                    generator = cand;
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        if generator == BinaryPoint::Infinity {
            return None;
        }
        curve.generator = generator;
        curve.order = r.clone();
        curve.cofactor = cofactor.clone();

        let trace: i64 = if a == 0 { -1 } else { 1 };
        let lambda = frobenius_eigenvalue(&curve, trace, &r)?;

        Some(Self {
            a,
            n,
            curve,
            trace,
            group_order,
            subgroup_order: r,
            cofactor,
            lambda,
        })
    }

    /// The `2`-power Frobenius `π(x, y) = (x², y²)`.
    pub fn frobenius(&self, p: &BinaryPoint) -> BinaryPoint {
        match p {
            BinaryPoint::Infinity => BinaryPoint::Infinity,
            BinaryPoint::Affine { x, y } => BinaryPoint::Affine {
                x: x.square(&self.curve.irreducible),
                y: y.square(&self.curve.irreducible),
            },
        }
    }

    /// `[k]·P` on this curve.
    pub fn mul(&self, p: &BinaryPoint, k: &BigUint) -> BinaryPoint {
        scalar_mul(&self.curve, p, k)
    }

    /// `P + Q` on this curve.
    pub fn add(&self, p: &BinaryPoint, q: &BinaryPoint) -> BinaryPoint {
        point_add(&self.curve, p, q)
    }

    /// The generator of the prime-order subgroup.
    pub fn generator(&self) -> &BinaryPoint {
        &self.curve.generator
    }
}

/// Both points of `E` with the given `x`-coordinate, or an empty vector
/// when `x` is not the abscissa of any `F_{2^n}`-point.
///
/// For `x ≠ 0` the substitution `y = x·u` turns `y² + xy = x³ + ax² + b`
/// into the Artin–Schreier equation `u² + u = x + a + b/x²`, solvable
/// iff the right-hand side has absolute trace `0`.  For `x = 0` the
/// equation degenerates to `y² = b`, whose unique solution is `y = √b`.
pub fn points_with_x(curve: &BinaryCurve, x: &F2mElement) -> Vec<BinaryPoint> {
    let irr = &curve.irreducible;
    let m = curve.m;
    if x.is_zero() {
        // y² = b ⇒ y = b^(2^(m−1)) (squaring is a bijection).
        let y = curve.b.square_k_times(m - 1, irr);
        return vec![BinaryPoint::Affine { x: x.clone(), y }];
    }
    let x_inv = match x.flt_inverse(irr) {
        Some(v) => v,
        None => return Vec::new(),
    };
    let rhs = x.add(&curve.a).add(&curve.b.mul(&x_inv.square(irr), irr));
    let u = match solve_artin_schreier(&rhs, m, irr) {
        Some(u) => u,
        None => return Vec::new(),
    };
    let y = x.mul(&u, irr);
    let p = BinaryPoint::Affine { x: x.clone(), y };
    let np = point_neg(&p);
    if p == np {
        vec![p]
    } else {
        vec![p, np]
    }
}

/// The scalar `λ` with `π(Q) = [λ]Q` for all `Q ∈ ⟨G⟩`.
///
/// `π` satisfies its characteristic equation `π² − tπ + 2 = 0`, so on
/// the order-`r` subgroup `λ` is one of the two roots of
/// `λ² − tλ + 2 ≡ 0 (mod r)`; which one is settled by testing
/// `π(G) = [λ]G`.
pub fn frobenius_eigenvalue(curve: &BinaryCurve, trace: i64, r: &BigUint) -> Option<BigUint> {
    let t = if trace >= 0 {
        BigUint::from(trace as u64) % r
    } else {
        r - (BigUint::from((-trace) as u64) % r)
    };
    // disc = t² − 8
    let t_sq = (&t * &t) % r;
    let disc = (&t_sq + r - (BigUint::from(8u32) % r)) % r;
    let s = sqrt_mod_p(&disc, r)?;
    let inv2 = mod_inverse(&BigUint::from(2u32), r)?;

    let g = &curve.generator;
    let pi_g = match g {
        BinaryPoint::Infinity => return None,
        BinaryPoint::Affine { x, y } => BinaryPoint::Affine {
            x: x.square(&curve.irreducible),
            y: y.square(&curve.irreducible),
        },
    };
    for sign in [true, false] {
        let num = if sign {
            (&t + &s) % r
        } else {
            (&t + r - &s) % r
        };
        let lambda = (num * &inv2) % r;
        if scalar_mul(curve, g, &lambda) == pi_g {
            return Some(lambda);
        }
    }
    None
}

// ── Frobenius-invariant factor base ────────────────────────────────

/// A Frobenius-invariant factor base and its orbit structure.
#[derive(Clone, Debug)]
pub struct FrobeniusFactorBase {
    /// Degree `ℓ = ord_n(2)` of the chosen factor `f_j` of `x^n − 1`.
    pub ell: u32,
    /// The chosen `f_j`, as an `F_2[x]` bitmask.
    pub f_j: u64,
    /// Exponents `k` with `f_{j,k} = 1`, i.e. the linearised polynomial
    /// is `F_j(X) = Σ X^{2^k}` over these `k`.
    pub linearised_exponents: Vec<u32>,
    /// The `2^ℓ` roots of `F_j` in `F_{2^n}` — an `F_2`-subspace closed
    /// under squaring.
    pub subspace: Vec<F2mElement>,
    /// An `F_2`-basis of that subspace, `ℓ` elements.  The Gröbner
    /// decomposition writes its unknowns in this basis, which is what
    /// keeps the Semaev system quadratic.
    pub subspace_basis: Vec<F2mElement>,
    /// The factor base itself: every point whose abscissa is a root.
    pub points: Vec<BinaryPoint>,
    /// `π`-orbits, as lists of indices into `points`.  Orbit `o` is
    /// `[i_0, i_1, …]` with `points[i_{k}] = π^k(points[i_0])`.
    pub orbits: Vec<Vec<usize>>,
    /// For each point index: `(orbit, k)` with `point = π^k(representative)`.
    pub orbit_of: Vec<(usize, u32)>,
}

impl FrobeniusFactorBase {
    /// Number of unknowns the linear algebra actually carries — one per
    /// `π`-orbit, versus `points.len()` for a non-invariant base.
    pub fn unknowns(&self) -> usize {
        self.orbits.len()
    }

    /// Lookup table from point identity to factor-base index, as both
    /// decomposition oracles need.
    pub fn index_map(&self) -> HashMap<(BigUint, BigUint), usize> {
        self.points
            .iter()
            .enumerate()
            .map(|(i, p)| (point_key(p), i))
            .collect()
    }
}

/// An `F_2`-**basis** of the kernel of the linearised polynomial
/// `F(X) = Σ_{k ∈ exps} X^{2^k}` acting `F_2`-linearly on `F_{2^n}`.
///
/// The map is written as an `n × n` matrix over `F_2` in the polynomial
/// basis `1, z, …, z^{n−1}`; the kernel basis is read off by row
/// reducing `[image | preimage]`.
pub fn linearised_kernel_basis(exps: &[u32], n: u32, irr: &IrreduciblePoly) -> Vec<F2mElement> {
    // Column i = F(z^i), packed into the low n bits of a u64.
    let mut rows: Vec<(u64, u64)> = Vec::with_capacity(n as usize);
    for i in 0..n {
        let basis = F2mElement::from_bit_positions(&[i], n);
        let mut img = F2mElement::zero(n);
        for &k in exps {
            img = img.add(&basis.square_k_times(k, irr));
        }
        let img_bits = img.raw_bits().first().copied().unwrap_or(0);
        rows.push((img_bits, 1u64 << i));
    }

    // Gaussian elimination on the image halves; whatever reduces to
    // zero contributes its preimage to the kernel basis.
    let mut pivots: Vec<(u64, u64)> = Vec::new();
    let mut kernel_basis: Vec<u64> = Vec::new();
    for (mut img, mut pre) in rows {
        for &(pimg, ppre) in &pivots {
            let lead = 1u64 << (63 - pimg.leading_zeros());
            if img & lead != 0 {
                img ^= pimg;
                pre ^= ppre;
            }
        }
        if img == 0 {
            kernel_basis.push(pre);
        } else {
            pivots.push((img, pre));
            pivots.sort_by(|x, y| y.0.cmp(&x.0));
        }
    }

    kernel_basis
        .into_iter()
        .map(|v| F2mElement::from_biguint(&BigUint::from(v), n))
        .collect()
}

/// The kernel itself: every `F_2`-combination of
/// [`linearised_kernel_basis`], `2^dim` elements.
pub fn linearised_kernel(exps: &[u32], n: u32, irr: &IrreduciblePoly) -> Vec<F2mElement> {
    span_f2(&linearised_kernel_basis(exps, n, irr), n)
}

/// All `2^{|basis|}` `F_2`-combinations of `basis`.
pub fn span_f2(basis: &[F2mElement], n: u32) -> Vec<F2mElement> {
    let mut out = Vec::with_capacity(1usize << basis.len());
    for mask in 0..(1u64 << basis.len()) {
        let mut acc = F2mElement::zero(n);
        for (b, base) in basis.iter().enumerate() {
            if (mask >> b) & 1 == 1 {
                acc = acc.add(base);
            }
        }
        out.push(acc);
    }
    out
}

/// **Build a Frobenius-invariant factor base** for `curve` from the
/// `index`-th non-trivial irreducible factor of `x^n − 1`.
///
/// `index` selects which of the `(n − 1)/ℓ` factors to use; different
/// factors give different (and, as the paper notes, sometimes needed —
/// a single invariant base may not yield `n` independent relations)
/// factor bases of the same size.
pub fn build_frobenius_factor_base(kc: &KoblitzCurve, index: usize) -> Option<FrobeniusFactorBase> {
    let factors = factor_x_n_minus_1(kc.n);
    let f_j = *factors.get(index)?;
    let ell = poly_deg(f_j)?;
    let exps: Vec<u32> = (0..=ell).filter(|k| (f_j >> k) & 1 == 1).collect();

    let subspace_basis = linearised_kernel_basis(&exps, kc.n, &kc.curve.irreducible);
    let subspace = span_f2(&subspace_basis, kc.n);
    if subspace.len() != (1usize << ell) {
        // Kernel dimension must equal deg f_j; anything else means the
        // factor did not divide x^n − 1 after all.
        return None;
    }

    let mut points: Vec<BinaryPoint> = Vec::new();
    for x in &subspace {
        for p in points_with_x(&kc.curve, x) {
            points.push(p);
        }
    }

    // Index points for the orbit walk and for relation lookups.
    let mut index_of: HashMap<(BigUint, BigUint), usize> = HashMap::new();
    for (i, p) in points.iter().enumerate() {
        index_of.insert(point_key(p), i);
    }

    let mut orbit_of: Vec<(usize, u32)> = vec![(usize::MAX, 0); points.len()];
    let mut orbits: Vec<Vec<usize>> = Vec::new();
    for start in 0..points.len() {
        if orbit_of[start].0 != usize::MAX {
            continue;
        }
        let o = orbits.len();
        let mut cycle = Vec::new();
        let mut cur = points[start].clone();
        let mut k = 0u32;
        loop {
            let idx = *index_of.get(&point_key(&cur))?;
            if orbit_of[idx].0 != usize::MAX {
                break;
            }
            orbit_of[idx] = (o, k);
            cycle.push(idx);
            cur = kc.frobenius(&cur);
            k += 1;
        }
        orbits.push(cycle);
    }

    Some(FrobeniusFactorBase {
        ell,
        f_j,
        linearised_exponents: exps,
        subspace,
        subspace_basis,
        points,
        orbits,
        orbit_of,
    })
}

/// Hashable identity of a point — the key of
/// [`FrobeniusFactorBase::index_map`].
pub fn point_key(p: &BinaryPoint) -> (BigUint, BigUint) {
    match p {
        BinaryPoint::Infinity => (BigUint::zero(), BigUint::zero()),
        BinaryPoint::Affine { x, y } => {
            // Shift by one so that no affine point collides with the
            // all-zero key reserved for O.
            (x.to_biguint() + BigUint::one(), y.to_biguint())
        }
    }
}

// ── Relations ──────────────────────────────────────────────────────

/// One relation `[a]G + [b]Q = Σ_i P_i` rewritten over the orbit
/// representatives of a Frobenius-invariant factor base.
#[derive(Clone, Debug)]
pub struct KoblitzRelation {
    /// `a` in `R = [a]G + [b]Q`.
    pub coef_a: BigUint,
    /// `b` in `R = [a]G + [b]Q`.
    pub coef_b: BigUint,
    /// The decomposition, as `(orbit, k)` pairs: summand
    /// `P_i = π^{k_i}(rep_{o_i})`.
    pub summands: Vec<(usize, u32)>,
    /// Dense row over the orbit unknowns: entry `o` is `Σ_i λ^{k_i}`
    /// over the summands lying in orbit `o`, reduced mod `r`.
    pub row: Vec<BigUint>,
}

/// Search for a decomposition `R = P_1 + … + P_m` with all `P_i` in the
/// factor base, returning the summands as `(orbit, k)` pairs.
///
/// **Symmetry breaking.**  Only non-decreasing index tuples are
/// enumerated, so each decomposition is found once rather than `m!`
/// times — the enumeration analogue of the paper's trick of spreading
/// the summands over `F_1, π(F_1), …, π^{m−1}(F_1)` and rewriting
/// everything back into `F_1`.
/// **Exhaustive decomposition** (the reference oracle): search ordered
/// tuples of factor-base points for `R = P_1 + … + P_m`.
///
/// Costs `|F|^{m−1}` group operations per target whether or not a
/// decomposition exists.  [`groebner_decompose`] answers the same
/// question algebraically; the two are cross-checked in the tests.
pub fn enumerate_decompose(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    target: &BinaryPoint,
    m: usize,
) -> Option<Vec<usize>> {
    decompose(kc, fb, index_of, target, m, 0)
}

fn decompose(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    target: &BinaryPoint,
    m: usize,
    start: usize,
) -> Option<Vec<usize>> {
    if m == 0 {
        return if *target == BinaryPoint::Infinity {
            Some(Vec::new())
        } else {
            None
        };
    }
    if m == 1 {
        let idx = *index_of.get(&point_key(target))?;
        return if idx >= start { Some(vec![idx]) } else { None };
    }
    for i in start..fb.points.len() {
        let rest = kc.add(target, &point_neg(&fb.points[i]));
        if let Some(mut tail) = decompose(kc, fb, index_of, &rest, m - 1, i) {
            let mut out = vec![i];
            out.append(&mut tail);
            return Some(out);
        }
    }
    None
}

/// How the decomposition oracle answers "is `R` a sum of `m`
/// factor-base points?".
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DecompositionStrategy {
    /// Weil-restrict the Semaev condition `S₃ = 0` to a low-degree
    /// Boolean system over the factor-base subspace and solve it with a
    /// Gröbner basis — the step that makes the real algorithm
    /// sub-exponential.  See
    /// [`crate::cryptanalysis::koblitz_groebner`].
    Groebner,
    /// Exhaustive ordered-tuple search over the materialised factor
    /// base.  Costs `|F|^{m−1}` group operations per target whether or
    /// not a decomposition exists; kept as the reference oracle the
    /// algebraic one is tested against.
    Enumerate,
    /// The same Semaev system, handed to the CDCL SAT solver instead of
    /// a Gröbner engine
    /// ([`crate::cryptanalysis::semaev_sat::encode_boolean_system`]).
    /// Clause learning rather than degree growth: the two blow up on
    /// different systems, which is the comparison Soos–Nohl–Castelluccia
    /// opened and this module lets you measure.
    Sat,
}

/// Lift a tuple of candidate `x`-coordinates to factor-base points
/// whose sum really is `target`.
///
/// `S₃ = 0` constrains abscissae only, so it pins the summands down to
/// sign: each `x` admits up to two points `(x, y)` and `(x, x+y)`, both
/// of which lie in the (negation-closed) factor base.  This walks the
/// `≤ 2^m` sign choices and returns the first that closes the group
/// identity, so a spurious root of the polynomial system can never
/// become a relation.
fn lift_candidate(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    xs: &[F2mElement],
    target: &BinaryPoint,
) -> Option<Vec<usize>> {
    fn walk(
        kc: &KoblitzCurve,
        fb: &FrobeniusFactorBase,
        index_of: &HashMap<(BigUint, BigUint), usize>,
        xs: &[F2mElement],
        depth: usize,
        acc: &BinaryPoint,
        chosen: &mut Vec<usize>,
        target: &BinaryPoint,
    ) -> bool {
        if depth == xs.len() {
            return acc == target;
        }
        for p in points_with_x(&kc.curve, &xs[depth]) {
            let idx = match index_of.get(&point_key(&p)) {
                Some(i) => *i,
                None => continue,
            };
            chosen.push(idx);
            let next = kc.add(acc, &fb.points[idx]);
            if walk(kc, fb, index_of, xs, depth + 1, &next, chosen, target) {
                return true;
            }
            chosen.pop();
        }
        false
    }

    let mut chosen = Vec::with_capacity(xs.len());
    if walk(
        kc,
        fb,
        index_of,
        xs,
        0,
        &BinaryPoint::Infinity,
        &mut chosen,
        target,
    ) {
        chosen.sort_unstable();
        Some(chosen)
    } else {
        None
    }
}

/// **Algebraic decomposition**: solve the Semaev system for
/// `R = P_1 + … + P_m` over the factor base.
///
/// Returns the factor-base indices of a decomposition, plus the solver
/// statistics, or `None` when the system has no root that lifts.  A
/// `None` from a completed solve is a *proof* that no decomposition
/// exists over this factor base; a `None` with `stats.exhausted` set
/// only means the node budget ran out.
pub fn groebner_decompose(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    engine: SolverEngine,
    node_budget: usize,
) -> (Option<Vec<usize>>, SolveStats) {
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => return (None, SolveStats::default()),
    };
    let sys = match build_decomposition_system(&fb.subspace_basis, &x_r, &kc.curve.b, m, st) {
        Some(sys) => sys,
        None => return (None, SolveStats::default()),
    };
    // A root of S₃ fixes the summands only up to sign, so some roots do
    // not lift.  Lift each as it is found and stop at the first that
    // closes the group identity.
    let opts = SolveOptions {
        engine,
        max_solutions: usize::MAX,
        node_budget,
    };
    let mut found: Option<Vec<usize>> = None;
    let (_, stats) = solve_boolean_system_filtered(&sys.equations, sys.n_vars, &opts, |root| {
        let xs: Vec<F2mElement> = (0..m)
            .map(|i| sys.summand_x(&fb.subspace_basis, root, i, kc.n))
            .collect();
        match lift_candidate(kc, fb, index_of, &xs, target) {
            Some(idxs) => {
                found = Some(idxs);
                true
            }
            None => false,
        }
    });
    (found, stats)
}

/// What a SAT decomposition attempt cost and concluded.
#[derive(Clone, Debug, Default)]
pub struct SatDecompositionStats {
    /// CDCL solve calls (one per model examined, plus the final one).
    pub solver_calls: usize,
    /// Roots the solver produced, including any that failed to lift.
    pub models: usize,
    /// The search ended in UNSAT, so "no decomposition" is *proven* —
    /// the SAT analogue of the Gröbner infeasibility certificate.
    pub refuted: bool,
    /// The model cap or the solver's conflict budget was hit, so a
    /// `None` result is inconclusive rather than a refutation.
    pub exhausted: bool,
    /// Models that did not satisfy the original system — always zero
    /// unless the CNF encoding is wrong.
    pub spurious: usize,
}

/// **SAT decomposition**: the same Semaev system as
/// [`groebner_decompose`], solved by CDCL.
///
/// Models are enumerated by blocking clause: each root that fails to
/// lift to a genuine point decomposition is excluded and the instance
/// re-solved, so a final UNSAT proves no decomposition exists.  Every
/// model is checked against the original equations before use, and the
/// lifted points are re-checked in the group, so neither an encoding
/// bug nor a sign ambiguity can produce a false relation.
pub fn sat_decompose(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    max_models: usize,
) -> (Option<Vec<usize>>, SatDecompositionStats) {
    let mut stats = SatDecompositionStats::default();
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => return (None, stats),
    };
    let sys = match build_decomposition_system(&fb.subspace_basis, &x_r, &kc.curve.b, m, st) {
        Some(sys) => sys,
        None => return (None, stats),
    };

    let mut blocked: Vec<u64> = Vec::new();
    loop {
        let mut enc = encode_boolean_system(sys.n_vars, &sys.equations, &blocked);
        stats.solver_calls += 1;
        if enc.trivially_unsat {
            stats.refuted = true;
            return (None, stats);
        }
        match enc.solver.solve() {
            SolveResult::Unsat => {
                stats.refuted = true;
                return (None, stats);
            }
            SolveResult::Unknown => {
                stats.exhausted = true;
                return (None, stats);
            }
            SolveResult::Sat => {
                let root = enc.model_assignment();
                if !sys.equations.iter().all(|e| e.eval(root) == 0) {
                    // Cannot happen with a correct encoding; block it
                    // rather than loop, and report it.
                    stats.spurious += 1;
                    blocked.push(root);
                    continue;
                }
                stats.models += 1;
                let xs: Vec<F2mElement> = (0..m)
                    .map(|i| sys.summand_x(&fb.subspace_basis, root, i, kc.n))
                    .collect();
                if let Some(idxs) = lift_candidate(kc, fb, index_of, &xs, target) {
                    return (Some(idxs), stats);
                }
                blocked.push(root);
                if blocked.len() >= max_models {
                    stats.exhausted = true;
                    return (None, stats);
                }
            }
        }
    }
}

/// Turn a decomposition into a relation row over the orbit unknowns.
///
/// With `x_o := log_G ([h]·rep_o)` and `P_i = π^{k_i}(rep_{o_i})`,
/// multiplying `R = Σ_i P_i` by the cofactor `h` gives
///
/// ```text
///     h·a + h·b·d  ≡  Σ_i λ^{k_i} · x_{o_i}   (mod r),
/// ```
///
/// because `π` acts as `[λ]` on the order-`r` subgroup and `[h]P_i`
/// lands in that subgroup whatever `P_i` was.
fn relation_from_decomposition(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    idxs: &[usize],
    coef_a: &BigUint,
    coef_b: &BigUint,
) -> KoblitzRelation {
    let r = &kc.subgroup_order;
    let mut row = vec![BigUint::zero(); fb.orbits.len()];
    let mut summands = Vec::with_capacity(idxs.len());
    for &i in idxs {
        let (o, k) = fb.orbit_of[i];
        let coeff = kc.lambda.modpow(&BigUint::from(k), r);
        row[o] = (&row[o] + coeff) % r;
        summands.push((o, k));
    }
    KoblitzRelation {
        coef_a: coef_a.clone(),
        coef_b: coef_b.clone(),
        summands,
        row,
    }
}

// ── Driver ─────────────────────────────────────────────────────────

/// Tuning knobs for [`koblitz_index_calculus_dlp`].
#[derive(Clone, Debug)]
pub struct KoblitzIcOptions {
    /// Number of factor-base points a relation decomposes into.
    pub m: usize,
    /// Index of the irreducible factor of `x^n − 1` defining the
    /// factor base.
    pub factor_index: usize,
    /// Relations to collect beyond the number of orbit unknowns.
    pub extra_relations: usize,
    /// Cap on random `(a, b)` trials during relation collection.
    pub max_trials: usize,
    /// Seed for the `(a, b)` sampler, so runs are reproducible.
    pub seed: u64,
    /// How to answer the decomposition question.
    pub strategy: DecompositionStrategy,
    /// Which algebraic engine reduces the Semaev system.  Ignored by
    /// [`DecompositionStrategy::Enumerate`].
    pub engine: SolverEngine,
    /// Splitting nodes (Gröbner-basis computations) one decomposition
    /// may spend before it gives up.  Ignored by
    /// [`DecompositionStrategy::Enumerate`].
    pub node_budget: usize,
    /// Models one [`DecompositionStrategy::Sat`] decomposition may
    /// examine before giving up.
    pub max_models: usize,
}

impl Default for KoblitzIcOptions {
    fn default() -> Self {
        Self {
            m: 2,
            factor_index: 0,
            extra_relations: 4,
            max_trials: 20_000,
            seed: 0x4b_6f_62_6c_69_74_7a_00, // "Koblitz\0"
            strategy: DecompositionStrategy::Groebner,
            engine: SolverEngine::default(),
            node_budget: 4096,
            max_models: 64,
        }
    }
}

/// What a run of the algorithm did.
#[derive(Clone, Debug)]
pub struct KoblitzIcReport {
    /// `|F|`.
    pub factor_base_size: usize,
    /// Number of `π`-orbits, i.e. the number of unknowns actually
    /// solved for.
    pub orbit_count: usize,
    /// `ℓ = ord_n(2)`, the dimension of the invariant subspace.
    pub ell: u32,
    /// Relations collected.
    pub relations: usize,
    /// Random `(a, b)` pairs tried.
    pub trials: usize,
    /// Recovered discrete logarithm, if the run succeeded.
    pub log: Option<BigUint>,
    /// Algebraic reductions (F4 passes or Gröbner bases) computed
    /// across every decomposition attempt.
    pub reductions: usize,
    /// Decomposition branches closed by a basis reducing to `{1}` —
    /// targets rejected algebraically instead of by search.
    pub infeasible_branches: usize,
    /// CDCL solve calls across every decomposition attempt.
    pub sat_calls: usize,
    /// Targets the SAT oracle refuted outright (UNSAT).
    pub sat_refutations: usize,
}

/// **Solve `Q = [d]·G`** on a Koblitz curve by index calculus over a
/// Frobenius-invariant factor base.
///
/// Returns the report; `report.log` carries the answer on success.  The
/// result is always verified by recomputing `[d]G` before it is
/// returned, so a `Some` is never a false positive.
pub fn koblitz_index_calculus_dlp(
    kc: &KoblitzCurve,
    q: &BinaryPoint,
    opts: &KoblitzIcOptions,
) -> Option<KoblitzIcReport> {
    let fb = build_frobenius_factor_base(kc, opts.factor_index)?;
    let r = &kc.subgroup_order;
    let g = kc.generator().clone();

    let index_of = fb.index_map();

    let mut report = KoblitzIcReport {
        factor_base_size: fb.points.len(),
        orbit_count: fb.orbits.len(),
        ell: fb.ell,
        relations: 0,
        trials: 0,
        log: None,
        reductions: 0,
        infeasible_branches: 0,
        sat_calls: 0,
        sat_refutations: 0,
    };
    if fb.points.is_empty() {
        return Some(report);
    }

    let field = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let wanted = fb.orbits.len() + opts.extra_relations;
    let mut relations: Vec<KoblitzRelation> = Vec::with_capacity(wanted);
    let mut rng = StdRng::seed_from_u64(opts.seed);
    let r_u64 = r.to_u64_digits().first().copied().unwrap_or(1).max(2);

    while relations.len() < wanted && report.trials < opts.max_trials {
        report.trials += 1;
        let a = BigUint::from(rng.gen_range(1..r_u64));
        let b = BigUint::from(rng.gen_range(1..r_u64));
        let target = kc.add(&kc.mul(&g, &a), &kc.mul(q, &b));

        // R = O is the lucky case: d = −a/b directly.
        if target == BinaryPoint::Infinity {
            let d = solve_for_d(&a, &b, r)?;
            if kc.mul(&g, &d) == *q {
                report.log = Some(d);
                report.relations = relations.len();
                return Some(report);
            }
            continue;
        }

        let found = match opts.strategy {
            DecompositionStrategy::Enumerate => decompose(kc, &fb, &index_of, &target, opts.m, 0),
            DecompositionStrategy::Groebner => {
                let (idxs, stats) = groebner_decompose(
                    kc,
                    &fb,
                    &index_of,
                    &field,
                    &target,
                    opts.m,
                    opts.engine,
                    opts.node_budget,
                );
                report.reductions += stats.reductions;
                report.infeasible_branches += stats.infeasible_branches;
                idxs
            }
            DecompositionStrategy::Sat => {
                let (idxs, stats) =
                    sat_decompose(kc, &fb, &index_of, &field, &target, opts.m, opts.max_models);
                report.sat_calls += stats.solver_calls;
                report.sat_refutations += usize::from(stats.refuted);
                idxs
            }
        };
        if let Some(idxs) = found {
            relations.push(relation_from_decomposition(kc, &fb, &idxs, &a, &b));
        }
    }
    report.relations = relations.len();
    if relations.len() < fb.orbits.len() + 1 {
        return Some(report);
    }

    // Unknowns: x_1 … x_s (orbit logs) and d, in the last column.
    //   Σ_o c_o x_o  −  (h·b)·d  ≡  h·a   (mod r)
    let s = fb.orbits.len();
    let h = &kc.cofactor % r;
    let mut matrix: Vec<Vec<BigUint>> = Vec::with_capacity(relations.len());
    let mut rhs: Vec<BigUint> = Vec::with_capacity(relations.len());
    for rel in &relations {
        let mut row = rel.row.clone();
        row.push((r - (&h * &rel.coef_b) % r) % r);
        matrix.push(row);
        rhs.push((&h * &rel.coef_a) % r);
    }
    let sol = gaussian_eliminate_mod_n(&mut matrix, &mut rhs, r)?;
    let d = sol.get(s)?.clone();
    if kc.mul(&g, &d) == *q {
        report.log = Some(d);
    }
    Some(report)
}

/// `d = −a / b (mod r)`, the degenerate relation `[a]G + [b]Q = O`.
fn solve_for_d(a: &BigUint, b: &BigUint, r: &BigUint) -> Option<BigUint> {
    let b_inv = mod_inverse(&(b % r), r)?;
    Some(((r - (a % r)) * b_inv) % r)
}

// ── Cost model ─────────────────────────────────────────────────────

/// Predicted savings from the Frobenius structure, in the units the
/// GGMP paper states them.
#[derive(Clone, Debug)]
pub struct KoblitzSpeedup {
    /// Extension degree.
    pub n: u32,
    /// Factor-base size `|F|`.
    pub factor_base_size: usize,
    /// Unknowns after collapsing `π`-orbits.
    pub unknowns: usize,
    /// `|F| / unknowns` — the measured relation-collection speed-up,
    /// which the theory predicts is `≈ n`.
    pub relation_collection_speedup: f64,
    /// `(|F| / unknowns)²` — the measured linear-algebra speed-up for a
    /// quadratic-cost sparse solve, predicted `≈ n²`.
    pub linear_algebra_speedup: f64,
    /// `m!`, the symmetry-breaking saving in the decomposition search.
    pub symmetry_breaking_factor: f64,
    /// The Gallant–Lambert–Vanstone / Wiener–Zuccherato speed-up that
    /// Pollard rho already gets from the same endomorphism, `√(2n)`.
    pub rho_speedup: f64,
}

/// Cost model for a Frobenius-invariant index calculus on `n`, a factor
/// base of `fb_size` points collapsing to `unknowns` orbits, with
/// `m`-point decompositions.
///
/// The point of the table is the comparison in the paper's conclusion:
/// index calculus gains more from the Koblitz structure than rho does
/// (`n` and `n²` versus `√(2n)`) — and yet stays the worse attack for
/// the curves anyone deploys, because its un-accelerated cost is so
/// much higher to begin with.
pub fn koblitz_speedup_model(n: u32, fb_size: usize, unknowns: usize, m: usize) -> KoblitzSpeedup {
    let ratio = if unknowns == 0 {
        1.0
    } else {
        fb_size as f64 / unknowns as f64
    };
    let fact = (1..=m).map(|i| i as f64).product::<f64>();
    KoblitzSpeedup {
        n,
        factor_base_size: fb_size,
        unknowns,
        relation_collection_speedup: ratio,
        linear_algebra_speedup: ratio * ratio,
        symmetry_breaking_factor: fact,
        rho_speedup: (2.0 * n as f64).sqrt(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn irreducibility_test_matches_known_polynomials() {
        // x² + x + 1, x³ + x + 1, x⁴ + x + 1 are irreducible.
        assert!(is_irreducible_f2(0b111));
        assert!(is_irreducible_f2(0b1011));
        assert!(is_irreducible_f2(0b10011));
        // x² + 1 = (x + 1)², x³ + 1 = (x + 1)(x² + x + 1).
        assert!(!is_irreducible_f2(0b101));
        assert!(!is_irreducible_f2(0b1001));
    }

    #[test]
    fn find_irreducible_returns_degree_n_irreducible() {
        for n in 3..=16u32 {
            let irr = find_irreducible(n).expect("degree-n irreducible exists");
            assert_eq!(irr.degree, n);
            let mut mask = 1u64 << n;
            for t in &irr.low_terms {
                mask |= 1 << t;
            }
            assert!(is_irreducible_f2(mask), "n = {n}");
        }
    }

    #[test]
    fn x_n_minus_1_factors_have_degree_ord_2_mod_n() {
        // x⁷ − 1 = (x − 1)(x³ + x + 1)(x³ + x² + 1); ord_7(2) = 3.
        assert_eq!(order_of_2_mod_n(7), Some(3));
        let f7 = factor_x_n_minus_1(7);
        assert_eq!(f7.len(), 2);
        assert!(f7.contains(&0b1011));
        assert!(f7.contains(&0b1101));

        // ord_9(2) = 6, so a single degree-6 factor x⁶ + x³ + 1.
        assert_eq!(order_of_2_mod_n(9), Some(6));
        let f9 = factor_x_n_minus_1(9);
        assert_eq!(f9.len(), 1);
        assert_eq!(f9[0], 0b1001001);
    }

    #[test]
    fn point_counts_match_direct_enumeration() {
        for a in 0..=1u8 {
            for n in 3..=9u32 {
                let irr = find_irreducible(n).unwrap();
                let curve = BinaryCurve {
                    m: n,
                    irreducible: irr,
                    a: if a == 0 {
                        F2mElement::zero(n)
                    } else {
                        F2mElement::one(n)
                    },
                    b: F2mElement::one(n),
                    generator: BinaryPoint::Infinity,
                    order: BigUint::zero(),
                    cofactor: BigUint::one(),
                };
                let mut count = 1u64; // O
                for raw in 0..(1u64 << n) {
                    let x = F2mElement::from_biguint(&BigUint::from(raw), n);
                    count += points_with_x(&curve, &x).len() as u64;
                }
                assert_eq!(
                    BigUint::from(count),
                    koblitz_point_count(a, n),
                    "a = {a}, n = {n}"
                );
            }
        }
    }

    #[test]
    fn frobenius_acts_as_lambda_on_the_subgroup() {
        let kc = KoblitzCurve::new(0, 9).expect("K_0 over F_512");
        assert_eq!(kc.subgroup_order, BigUint::from(127u32));
        let g = kc.generator().clone();
        // π = [λ] must hold across the whole subgroup, not just on G.
        for k in 1..20u32 {
            let p = kc.mul(&g, &BigUint::from(k));
            assert_eq!(kc.frobenius(&p), kc.mul(&p, &kc.lambda));
        }
        // λ is a root of the characteristic polynomial λ² − tλ + 2.
        let r = &kc.subgroup_order;
        let t = r - BigUint::one(); // t = −1 for a = 0
        let lhs = (&kc.lambda * &kc.lambda + BigUint::from(2u32)) % r;
        assert_eq!(lhs, (&t * &kc.lambda) % r);
    }

    #[test]
    fn factor_base_is_frobenius_invariant() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        assert_eq!(fb.ell, 6);
        assert_eq!(fb.subspace.len(), 64);

        let keys: std::collections::HashSet<_> = fb.points.iter().map(point_key).collect();
        for p in &fb.points {
            // π(F) = F, and the base is closed under negation too.
            assert!(keys.contains(&point_key(&kc.frobenius(p))));
            assert!(keys.contains(&point_key(&point_neg(p))));
        }

        // Orbit bookkeeping: point = π^k(representative).
        for (i, p) in fb.points.iter().enumerate() {
            let (o, k) = fb.orbit_of[i];
            let mut rep = fb.points[fb.orbits[o][0]].clone();
            for _ in 0..k {
                rep = kc.frobenius(&rep);
            }
            assert_eq!(&rep, p);
        }
        // The orbits partition the factor base.
        assert_eq!(
            fb.orbits.iter().map(|o| o.len()).sum::<usize>(),
            fb.points.len()
        );
        // Collapsing them is where the factor-n saving comes from.
        assert!(fb.unknowns() * 2 < fb.points.len());
    }

    #[test]
    fn subspace_is_closed_under_squaring() {
        let kc = KoblitzCurve::new(1, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let set: std::collections::HashSet<BigUint> =
            fb.subspace.iter().map(|v| v.to_biguint()).collect();
        for v in &fb.subspace {
            let sq = v.square(&kc.curve.irreducible);
            assert!(set.contains(&sq.to_biguint()));
        }
    }

    #[test]
    fn relation_rows_are_consistent_with_the_orbit_logs() {
        // A relation must satisfy h·a + h·b·d ≡ Σ λ^k x_o with
        // x_o = log_G([h]·rep_o); check that against brute-forced logs.
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let r = &kc.subgroup_order;
        let g = kc.generator().clone();

        // Brute-force log table for the (tiny) subgroup.
        let mut logs: HashMap<(BigUint, BigUint), BigUint> = HashMap::new();
        let mut acc = BinaryPoint::Infinity;
        let mut i = BigUint::zero();
        while &i < r {
            logs.insert(point_key(&acc), i.clone());
            acc = kc.add(&acc, &g);
            i += BigUint::one();
        }

        let mut index_of = HashMap::new();
        for (idx, p) in fb.points.iter().enumerate() {
            index_of.insert(point_key(p), idx);
        }
        let a = BigUint::from(11u32);
        let b = BigUint::from(23u32);
        let d = BigUint::from(37u32);
        let q = kc.mul(&g, &d);
        let target = kc.add(&kc.mul(&g, &a), &kc.mul(&q, &b));
        let idxs = decompose(&kc, &fb, &index_of, &target, 2, 0)
            .expect("a 2-point decomposition exists for this target");
        let rel = relation_from_decomposition(&kc, &fb, &idxs, &a, &b);

        let h = &kc.cofactor % r;
        let lhs = (&h * &a + &h * &b * &d) % r;
        let mut rhs = BigUint::zero();
        for (o, coeff) in rel.row.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            let rep = fb.points[fb.orbits[o][0]].clone();
            let scaled = kc.mul(&rep, &kc.cofactor);
            let x_o = logs.get(&point_key(&scaled)).expect("in ⟨G⟩ after ×h");
            rhs = (rhs + coeff * x_o) % r;
        }
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn all_three_oracles_answer_every_target_identically() {
        // Exhaustive search, Gröbner and SAT must agree on *both*
        // answers: the same decomposability verdict, and whatever they
        // return must be a genuine decomposition.
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let index_of = fb.index_map();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let g = kc.generator().clone();

        let mut agreed = 0;
        let mut decomposable = 0;
        for k in 1..40u32 {
            let target = kc.mul(&g, &BigUint::from(k));
            let by_search = enumerate_decompose(&kc, &fb, &index_of, &target, 2);
            let (by_algebra, stats) = groebner_decompose(
                &kc,
                &fb,
                &index_of,
                &st,
                &target,
                2,
                SolverEngine::default(),
                20_000,
            );
            assert!(!stats.exhausted, "budget should suffice at n = 9");
            let (by_sat, sat_stats) = sat_decompose(&kc, &fb, &index_of, &st, &target, 2, 64);
            assert!(!sat_stats.exhausted, "model cap should suffice at n = 9");
            assert_eq!(sat_stats.spurious, 0, "the CNF encoding must be exact");
            assert_eq!(
                by_search.is_some(),
                by_algebra.is_some(),
                "search and Gröbner disagree on [{k}]G"
            );
            assert_eq!(
                by_search.is_some(),
                by_sat.is_some(),
                "search and SAT disagree on [{k}]G"
            );
            if let Some(idxs) = by_sat {
                let mut acc = BinaryPoint::Infinity;
                for i in &idxs {
                    acc = kc.add(&acc, &fb.points[*i]);
                }
                assert_eq!(acc, target, "SAT decomposition of [{k}]G is wrong");
            }
            if let Some(idxs) = by_algebra {
                decomposable += 1;
                // Whatever it returned must actually sum to the target.
                let mut acc = BinaryPoint::Infinity;
                for i in &idxs {
                    acc = kc.add(&acc, &fb.points[*i]);
                }
                assert_eq!(acc, target, "algebraic decomposition of [{k}]G is wrong");
            }
            agreed += 1;
        }
        assert_eq!(agreed, 39);
        assert!(decomposable > 0, "some targets must decompose");
    }

    #[test]
    fn chained_s3_handles_three_summands() {
        // m ≥ 3 chains S₃ over intermediate field unknowns instead of
        // resolving S_{m+1}; check the chain against exhaustive search.
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let index_of = fb.index_map();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let g = kc.generator().clone();

        for k in [5u32, 29, 88] {
            let target = kc.mul(&g, &BigUint::from(k));
            let by_search = enumerate_decompose(&kc, &fb, &index_of, &target, 3);
            let (by_algebra, _) = groebner_decompose(
                &kc,
                &fb,
                &index_of,
                &st,
                &target,
                3,
                SolverEngine::default(),
                50_000,
            );
            assert_eq!(by_search.is_some(), by_algebra.is_some(), "3-point, [{k}]G");
            if let Some(idxs) = by_algebra {
                assert_eq!(idxs.len(), 3);
                let mut acc = BinaryPoint::Infinity;
                for i in &idxs {
                    acc = kc.add(&acc, &fb.points[*i]);
                }
                assert_eq!(acc, target);
            }
        }
    }

    #[test]
    fn the_two_algebraic_engines_agree() {
        // Matrix-F4 and textbook Buchberger must reach the same verdict;
        // only the cost differs.
        let kc = KoblitzCurve::new(1, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let index_of = fb.index_map();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let g = kc.generator().clone();

        for k in 1..6u32 {
            let target = kc.mul(&g, &BigUint::from(k));
            let (f4, _) = groebner_decompose(
                &kc,
                &fb,
                &index_of,
                &st,
                &target,
                2,
                SolverEngine::MatrixF4 { max_degree: 3 },
                20_000,
            );
            let (bb, _) = groebner_decompose(
                &kc,
                &fb,
                &index_of,
                &st,
                &target,
                2,
                SolverEngine::Buchberger,
                20_000,
            );
            assert_eq!(f4.is_some(), bb.is_some(), "engines disagree on [{k}]G");
        }
    }

    #[test]
    fn undecomposable_targets_are_rejected_algebraically() {
        // The point of the Gröbner step: a target with no decomposition
        // is *refuted* by the algebra — the reduction produces the
        // constant 1 — rather than by searching the factor base.
        //
        // `K_0 / F_2^7` is the clean case: its invariant subspace has 8
        // elements but only one of them is the abscissa of a curve
        // point, so no target decomposes at all.
        let kc = KoblitzCurve::new(0, 7).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        assert_eq!(fb.points.len(), 1, "the degenerate factor base");
        let index_of = fb.index_map();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let g = kc.generator().clone();

        let mut certified = 0;
        for k in 1..29u32 {
            let target = kc.mul(&g, &BigUint::from(k));
            assert!(enumerate_decompose(&kc, &fb, &index_of, &target, 2).is_none());
            let (out, stats) = groebner_decompose(
                &kc,
                &fb,
                &index_of,
                &st,
                &target,
                2,
                SolverEngine::default(),
                20_000,
            );
            assert!(out.is_none(), "[{k}]G must not decompose");
            assert!(
                stats.infeasible_branches > 0,
                "[{k}]G should be closed by an infeasibility certificate"
            );
            assert_eq!(stats.splits, 0, "no search should have been needed");
            certified += 1;
        }
        assert_eq!(certified, 28);
    }

    #[test]
    fn sat_refutes_undecomposable_targets() {
        // The SAT analogue of the Gröbner certificate: on the
        // degenerate `K_0 / F_2^7` factor base nothing decomposes, and
        // every instance comes back UNSAT rather than merely unsolved.
        let kc = KoblitzCurve::new(0, 7).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let index_of = fb.index_map();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let g = kc.generator().clone();

        for k in 1..29u32 {
            let target = kc.mul(&g, &BigUint::from(k));
            let (out, stats) = sat_decompose(&kc, &fb, &index_of, &st, &target, 2, 64);
            assert!(out.is_none(), "[{k}]G must not decompose");
            assert!(stats.refuted, "[{k}]G should be refuted by UNSAT");
            assert!(!stats.exhausted);
            assert_eq!(stats.spurious, 0);
        }
    }

    #[test]
    fn sat_handles_the_cubic_chained_system() {
        // With m ≥ 3 the chained S₃ is cubic, so the encoding needs
        // Tseitin auxiliaries above degree 2.
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let index_of = fb.index_map();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let g = kc.generator().clone();

        for k in [5u32, 29] {
            let target = kc.mul(&g, &BigUint::from(k));
            let by_search = enumerate_decompose(&kc, &fb, &index_of, &target, 3);
            let (by_sat, stats) = sat_decompose(&kc, &fb, &index_of, &st, &target, 3, 64);
            assert_eq!(stats.spurious, 0);
            assert_eq!(by_search.is_some(), by_sat.is_some(), "3-point, [{k}]G");
            if let Some(idxs) = by_sat {
                let mut acc = BinaryPoint::Infinity;
                for i in &idxs {
                    acc = kc.add(&acc, &fb.points[*i]);
                }
                assert_eq!(acc, target);
            }
        }
    }

    #[test]
    fn solves_the_dlp_with_the_sat_oracle() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let g = kc.generator().clone();
        let d = BigUint::from(53u32);
        let q = kc.mul(&g, &d);
        let opts = KoblitzIcOptions {
            strategy: DecompositionStrategy::Sat,
            ..KoblitzIcOptions::default()
        };
        let report = koblitz_index_calculus_dlp(&kc, &q, &opts).unwrap();
        assert_eq!(report.log, Some(d));
        assert!(report.sat_calls > 0, "the SAT path must have run");
        assert_eq!(report.reductions, 0, "no Gröbner work on this path");
    }

    #[test]
    fn solves_the_dlp_by_enumeration_too() {
        // The reference oracle drives the same solver to the same log.
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let g = kc.generator().clone();
        let d = BigUint::from(53u32);
        let q = kc.mul(&g, &d);
        let opts = KoblitzIcOptions {
            strategy: DecompositionStrategy::Enumerate,
            ..KoblitzIcOptions::default()
        };
        let report = koblitz_index_calculus_dlp(&kc, &q, &opts).unwrap();
        assert_eq!(report.log, Some(d));
        // No algebra ran on this path.
        assert_eq!(report.reductions, 0);
    }

    #[test]
    fn solves_the_dlp_on_k0_over_f512() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let g = kc.generator().clone();
        for d in [
            BigUint::from(5u32),
            BigUint::from(61u32),
            BigUint::from(118u32),
        ] {
            let q = kc.mul(&g, &d);
            let opts = KoblitzIcOptions::default();
            let report = koblitz_index_calculus_dlp(&kc, &q, &opts).unwrap();
            assert_eq!(report.log, Some(d.clone()), "report = {report:?}");
        }
    }

    #[test]
    fn solves_the_dlp_on_k1_over_f512() {
        let kc = KoblitzCurve::new(1, 9).unwrap();
        let g = kc.generator().clone();
        let d = BigUint::from(29u32);
        let q = kc.mul(&g, &d);
        let report = koblitz_index_calculus_dlp(&kc, &q, &KoblitzIcOptions::default()).unwrap();
        assert_eq!(report.log, Some(d));
    }

    #[test]
    fn speedup_model_tracks_the_predicted_factors() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let model = koblitz_speedup_model(kc.n, fb.points.len(), fb.unknowns(), 2);
        // Orbits have size n except for the short ones, so the measured
        // relation-collection saving sits just below n.
        assert!(model.relation_collection_speedup <= kc.n as f64);
        assert!(model.relation_collection_speedup > kc.n as f64 * 0.5);
        assert!(
            (model.linear_algebra_speedup
                - model.relation_collection_speedup * model.relation_collection_speedup)
                .abs()
                < 1e-9
        );
        assert_eq!(model.symmetry_breaking_factor, 2.0);
        // The paper's headline comparison: more than rho's √(2n).
        assert!(model.relation_collection_speedup > model.rho_speedup);
    }
}
