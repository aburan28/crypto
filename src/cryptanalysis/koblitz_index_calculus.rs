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
//! Negation adds `log_P(-Q) = -log_P(Q)`, so the `|F|` unknowns of a
//! classical index calculus collapse to one per signed Frobenius orbit,
//! generically `≈ |F| / (2n)`. Consequences:
//!
//! - **Relation collection:** we need one relation per signed orbit
//!   unknown instead of one per factor-base point.
//! - **Linear algebra:** Frobenius shrinks the relation matrix by a factor `n` in
//!   *both* dimensions, so a quadratic (sparse-linear-algebra) solve
//!   costs `n²` less; negation removes the redundant opposite columns.
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
//! The implementation also supports smaller invariant nonlinear sets:
//! a seed subspace can be closed under Frobenius, then optionally under
//! translation by rational 2-torsion.  These domains are represented by
//! their exact finite coordinate set.  SAT constrains every summand with
//! a prefix trie for that set; it does not replace the nonlinear domain
//! with its full linear span.
//!
//! ## Honest scope
//!
//! - **Three decomposition oracles.**  "Is `R` a sum of `m` factor-base
//!   points?" can be answered by
//!   [`DecompositionStrategy::Groebner`] — Semaev's `S₃` Weil-restricted
//!   to a low-degree Boolean system over the invariant subspace and
//!   solved with matrix-F4 (see
//!   [`crate::cryptanalysis::koblitz_groebner`]) — by
//!   [`DecompositionStrategy::Sat`], which encodes the descended
//!   equations with native XOR rows and exact factor-base membership,
//!   or by
//!   [`DecompositionStrategy::Enumerate`], a table-driven search over
//!   ordered tuples costing `|F|^{m−1}` group operations.  They are
//!   cross-checked against each other on the toy cases in the tests.
//!   At the sizes this module can reach the search is still 10–100×
//!   *faster* in wall-clock terms: `|F|` is a few dozen points, so
//!   `|F|^{m−1}` is nothing, while the Macaulay matrix already has
//!   thousands of columns.  A completed Gröbner or SAT refutation proves
//!   that the encoded decomposition does not exist; an exhausted budget
//!   remains inconclusive.  The exact factor-base coverage results do not
//!   establish a SAT speedup or a solving-complexity bound.
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

use std::collections::{HashMap, HashSet};

use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use crate::binary_ecc::curve::{point_add, point_neg, scalar_mul};
use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crate::cryptanalysis::binary_semaev::solve_artin_schreier;
use crate::cryptanalysis::ec_index_calculus::{gaussian_eliminate_mod_n, sqrt_mod_p};
use crate::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, matrix_f4_f2, solve_boolean_system_filtered, FieldStructure,
    SolveOptions, SolveStats, SolverEngine,
};
use crate::cryptanalysis::koblitz_relation_solver::{IncrementalRelationSolver, RowStatus};
use crate::cryptanalysis::sat::SolveResult;
use crate::cryptanalysis::semaev_sat::{encode_boolean_system_with, XorEncoding};
use crate::utils::mod_inverse;

/// Largest extension degree this module will build a curve for.  The
/// factor base and the point-counting/factoring helpers are all
/// materialised, so this is a deliberate guard rail, not a limit of
/// the mathematics.  Curve construction costs trial division to
/// `√#E ≈ 2^{n/2}` and a sparse irreducible search; both are cheap to
/// `n = 40`.  What actually bounds a run is the `2^dim` factor base
/// and, for the meet-in-the-middle oracle, its `|F|²` pair table.
pub const MAX_N: u32 = 40;

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

/// **Sparse irreducible search**: the smallest-mask irreducible
/// polynomial of degree `n` over `F_2` among those with at most four
/// terms below `z^n` — i.e. trinomials and pentanomials.
///
/// [`find_irreducible`] scans all `2^n` masks, which is fine to `n ≈ 24`
/// and hopeless beyond.  For every degree in range the smallest
/// irreducible polynomial *is* sparse, so this returns the same answer
/// far faster (a test pins the agreement for `n ≤ 20`), and it is what
/// the measurement harness uses to reach `n = 63`.
pub fn find_irreducible_sparse(n: u32) -> Option<IrreduciblePoly> {
    if n == 0 || n >= 64 {
        return None;
    }
    if n == 1 {
        return Some(IrreduciblePoly {
            degree: 1,
            low_terms: vec![0],
        });
    }
    // An irreducible polynomial of degree ≥ 1 has a non-zero constant
    // term, so bit 0 is always set; try 0, 1, 2 further bits below n.
    let mut candidates: Vec<u64> = Vec::new();
    candidates.push(1);
    for i in 1..n {
        candidates.push(1 | (1 << i));
        for j in (i + 1)..n {
            candidates.push(1 | (1 << i) | (1 << j));
            for k in (j + 1)..n {
                candidates.push(1 | (1 << i) | (1 << j) | (1 << k));
            }
        }
    }
    candidates.sort_unstable();
    let hi = 1u64 << n;
    for low in candidates {
        if is_irreducible_f2(hi | low) {
            return Some(IrreduciblePoly {
                degree: n,
                low_terms: (0..n).filter(|i| (low >> i) & 1 == 1).collect(),
            });
        }
    }
    None
}

/// The `index`-th Frobenius-invariant subspace of `F_{2^n}`, as an
/// `F_2`-basis — the factor base's `x`-coordinate support, without
/// building the curve, counting its points, or materialising `F`.
///
/// Point counting and factoring cap [`KoblitzCurve::new`] at
/// [`MAX_N`]; the *system* the decomposition oracles solve needs only
/// the field and this subspace, so it can be measured much further out.
/// Returns `None` when `x^n − 1` has no `index`-th non-trivial factor.
pub fn invariant_subspace_basis(
    n: u32,
    index: usize,
) -> Option<(IrreduciblePoly, Vec<F2mElement>)> {
    let irr = find_irreducible_sparse(n)?;
    let f_j = *factor_x_n_minus_1(n).get(index)?;
    let ell = 63 - f_j.leading_zeros();
    let exps: Vec<u32> = (0..=ell).filter(|k| (f_j >> k) & 1 == 1).collect();
    let basis = linearised_kernel_basis(&exps, n, &irr);
    if basis.len() != ell as usize {
        return None;
    }
    Some((irr, basis))
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

/// **The degree-`ord_n(2)` irreducible factors of `x^n − 1`** over
/// `F_2`, found by scanning the `2^ℓ` monic polynomials of that degree.
///
/// **Not a complete factorisation.**  `x^n − 1` has one irreducible
/// factor per 2-cyclotomic coset mod `n`, of degree equal to that
/// coset's size, and `ord_n(2)` is only the *largest* of those sizes.
/// At `n = 9` the cosets are `{0}`, `{3, 6}` and the six-element rest,
/// so `x^9 − 1` also has a degree-2 factor that this function does not
/// return.  Use [`all_factors_of_x_n_minus_1`] for the full list; this
/// one is kept as-is because
/// [`build_frobenius_factor_base`] indexes into it.
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
        // Identical to the exhaustive `find_irreducible` wherever both
        // are defined (a test pins that for every n ≤ 24) and far
        // cheaper past it: the exhaustive scan would touch 2^n masks.
        let irreducible = find_irreducible_sparse(n)?;
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

/// Algebraic description of the allowed abscissae.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FactorBaseDomain {
    LinearSubspace,
    /// A subset of the signed Frobenius orbits of a linear subspace
    /// base, selected by the factor-base search.  The coordinates still
    /// live in the subspace (so the algebraic system keeps its `ℓ`
    /// variables per summand); SAT additionally constrains each summand
    /// to the retained abscissae with the exact domain trie.
    SubspaceSubset {
        retained_orbits: usize,
    },
    /// Union of Frobenius translates of an ell-dimensional seed space.
    /// The union need not be a linear subspace.
    FrobeniusUnion {
        seed_dimension: usize,
    },
    /// Explicit abscissa representatives, closed under Frobenius.
    ExplicitFrobeniusOrbits {
        representatives: usize,
    },
    /// Reciprocal closure from translation by rational 2-torsion.
    TwoTorsionSaturation,
}

/// A Frobenius-invariant factor base and its orbit structure.
#[derive(Clone, Debug)]
pub struct FrobeniusFactorBase {
    pub domain: FactorBaseDomain,
    /// Dimension of the linear ambient coordinates used for descent.
    /// For a nonlinear union this is n, not log2 of the base size.
    pub ell: u32,
    /// The chosen `f_j`, as an `F_2[x]` bitmask; zero for a nonlinear union.
    pub f_j: u64,
    /// Exponents `k` with `f_{j,k} = 1`, i.e. the linearised polynomial
    /// is `F_j(X) = Σ X^{2^k}` over these `k`.
    pub linearised_exponents: Vec<u32>,
    /// Allowed abscissae, closed under squaring. For LinearSubspace
    /// these are the 2^ell roots of F_j; for FrobeniusUnion they are a union.
    pub subspace: Vec<F2mElement>,
    /// Basis of the ambient linear space used to encode coordinates.
    /// A union uses the full field basis and requires domain constraints.
    pub subspace_basis: Vec<F2mElement>,
    /// The factor base itself: every point whose abscissa is a root.
    pub points: Vec<BinaryPoint>,
    /// `π`-orbits, as lists of indices into `points`.  Orbit `o` is
    /// `[i_0, i_1, …]` with `points[i_{k}] = π^k(points[i_0])`.
    pub orbits: Vec<Vec<usize>>,
    /// For each point index: `(orbit, k)` with `point = π^k(representative)`.
    pub orbit_of: Vec<(usize, u32)>,
    /// Orbits under Frobenius and negation. These are the relation
    /// columns because `log(-P) = -log(P)` adds no independent unknown.
    pub signed_orbits: Vec<Vec<usize>>,
    /// For each point: `(signed orbit, k, negated)` with
    /// `point = (-1)^negated π^k(representative)`.
    pub signed_orbit_of: Vec<(usize, u32, bool)>,
}

impl FrobeniusFactorBase {
    /// Number of unknowns the linear algebra actually carries — one per
    /// signed `π`-orbit, versus `points.len()` for a non-invariant base.
    pub fn unknowns(&self) -> usize {
        self.signed_orbits.len()
    }

    /// Whether the abscissae are encoded in the full polynomial basis
    /// of `F_{2^n}` (nonlinear unions and saturations) rather than in
    /// the `ℓ`-dimensional basis of an invariant subspace.
    pub fn uses_ambient_basis(&self) -> bool {
        !matches!(
            self.domain,
            FactorBaseDomain::LinearSubspace | FactorBaseDomain::SubspaceSubset { .. }
        )
    }

    /// Canonical abscissa of each signed orbit: the smallest `x` (as an
    /// integer) among its points.  Stable across rebuilds, so a search
    /// can name orbits by it.
    pub fn signed_orbit_abscissa_representatives(&self) -> Vec<BigUint> {
        self.signed_orbits
            .iter()
            .map(|orbit| {
                orbit
                    .iter()
                    .filter_map(|&i| match &self.points[i] {
                        BinaryPoint::Affine { x, .. } => Some(x.to_biguint()),
                        BinaryPoint::Infinity => None,
                    })
                    .min()
                    .expect("signed orbits contain affine points")
            })
            .collect()
    }

    /// The class `[r]P` of each factor-base point in the `h`-torsion,
    /// where `h = #E / r` is the cofactor.
    ///
    /// A target `R ∈ ⟨G⟩` has `[r]R = O`, so a decomposition
    /// `R = Σ P_i` forces `Σ [r]P_i = O`.  When the factor base does not
    /// meet `⟨G⟩` — which is the common case, since `x = 0` lies in
    /// every invariant subspace and carries the 2-torsion point — these
    /// classes are non-trivial and constrain which `m` can work at all.
    pub fn cofactor_classes(&self, kc: &KoblitzCurve) -> Vec<BinaryPoint> {
        self.points
            .iter()
            .map(|p| kc.mul(p, &kc.subgroup_order))
            .collect()
    }

    /// Distinct h-torsion classes `[r]P`, derived with one scalar
    /// multiplication per signed Frobenius orbit.
    ///
    /// Multiplication by `r` commutes with Frobenius and negation, so the
    /// remaining classes in an orbit are recovered with cheap public group
    /// operations. Multiplicity is irrelevant to exact-summand reachability
    /// because decomposition permits repeated factor-base points.
    pub fn distinct_cofactor_classes(&self, kc: &KoblitzCurve) -> Vec<BinaryPoint> {
        let mut classes = HashMap::new();
        for orbit in &self.signed_orbits {
            let Some(&representative) = orbit.first() else {
                continue;
            };
            let mut current = kc.mul(&self.points[representative], &kc.subgroup_order);
            for _ in 0..kc.n {
                for candidate in [current.clone(), point_neg(&current)] {
                    classes.entry(point_key(&candidate)).or_insert(candidate);
                }
                current = kc.frobenius(&current);
            }
        }
        classes.into_values().collect()
    }

    /// **Can an `m`-point decomposition of a target in `⟨G⟩` exist?**
    ///
    /// Not a statement about size — about cosets.  Every summand
    /// contributes its `h`-torsion class, and the classes must cancel.
    /// On `K_1 / F_2^7` every factor-base point sits in the non-trivial
    /// class of a cofactor-2 curve, so odd `m` decomposes *nothing* however
    /// large the factor base is, and even `m` decomposes everything.
    /// Checking this before a sweep costs one scalar multiplication per
    /// point and saves searching for decompositions that cannot exist.
    pub fn m_can_decompose(&self, kc: &KoblitzCurve, m: usize) -> bool {
        if self.points.is_empty() || m == 0 {
            return m == 0;
        }
        let classes = self.distinct_cofactor_classes(kc);
        if m == 1 {
            return classes.iter().any(|c| *c == BinaryPoint::Infinity);
        }
        // `Σ c_i = O` with `m` summands ⇔ `−c_m` lies in the (m − 1)-fold
        // sumset, and the classes are closed under negation, so it is
        // enough to build `m − 1` layers and intersect with the classes.
        //
        // Every layer is closed under Frobenius and negation (the
        // classes are, and both commute with addition), so a layer is
        // generated by adding the classes to one representative per
        // signed Frobenius orbit of the previous layer and closing the
        // result under squaring and negation — cheap field squarings in
        // place of the `|layer| · |classes|` point additions the naive
        // walk pays.  Each layer holds at most `h` points.
        let identity = pack_point(&BinaryPoint::Infinity);
        let class_keys: HashSet<u64> = classes.iter().map(pack_point).collect();
        let mut layer: Vec<BinaryPoint> = vec![BinaryPoint::Infinity];
        let mut layer_keys: HashSet<u64> = HashSet::from([identity]);
        for _ in 1..m {
            let reps = signed_frobenius_orbit_representatives(kc, &layer);
            let mut next_keys: HashSet<u64> = HashSet::new();
            let mut next: Vec<BinaryPoint> = Vec::new();
            for q in &reps {
                for c in &classes {
                    let seed = kc.add(q, c);
                    if next_keys.contains(&pack_point(&seed)) {
                        continue;
                    }
                    // Close the orbit of `seed` under π and negation.
                    let mut current = seed;
                    for _ in 0..kc.n {
                        for candidate in [current.clone(), point_neg(&current)] {
                            if next_keys.insert(pack_point(&candidate)) {
                                next.push(candidate);
                            }
                        }
                        current = kc.frobenius(&current);
                    }
                }
            }
            layer = next;
            layer_keys = next_keys;
        }
        layer_keys.iter().any(|k| class_keys.contains(k))
    }

    /// The naive `m`-fold sumset walk over the cofactor classes —
    /// `|layer| · |classes|` point additions per layer.  Kept as the
    /// reference the orbit-generated walk is tested against.
    #[cfg(test)]
    fn m_can_decompose_reference(&self, kc: &KoblitzCurve, m: usize) -> bool {
        if self.points.is_empty() || m == 0 {
            return m == 0;
        }
        let classes = self.distinct_cofactor_classes(kc);
        let identity = pack_point(&BinaryPoint::Infinity);
        let mut layer: Vec<BinaryPoint> = vec![BinaryPoint::Infinity];
        for j in 1..=m {
            let mut seen: HashSet<u64> = HashSet::new();
            let mut next: Vec<BinaryPoint> = Vec::new();
            for acc in &layer {
                for c in &classes {
                    let sum = kc.add(acc, c);
                    if seen.insert(pack_point(&sum)) {
                        next.push(sum);
                    }
                }
            }
            if j == m {
                return seen.contains(&identity);
            }
            layer = next;
        }
        false
    }

    /// The summand counts in `2 ..= max_m` that [`Self::m_can_decompose`]
    /// admits — the ones worth trying.
    pub fn admissible_summand_counts(&self, kc: &KoblitzCurve, max_m: usize) -> Vec<usize> {
        (2..=max_m)
            .filter(|&m| self.m_can_decompose(kc, m))
            .collect()
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

/// One representative per orbit of a point set under Frobenius and
/// negation (the set is assumed closed under both).
fn signed_frobenius_orbit_representatives(kc: &KoblitzCurve, points: &[BinaryPoint]) -> Vec<BinaryPoint> {
    let mut seen: HashSet<u64> = HashSet::with_capacity(points.len());
    let mut reps = Vec::new();
    for p in points {
        if seen.contains(&pack_point(p)) {
            continue;
        }
        reps.push(p.clone());
        let mut current = p.clone();
        for _ in 0..kc.n {
            seen.insert(pack_point(&current));
            seen.insert(pack_point(&point_neg(&current)));
            current = kc.frobenius(&current);
        }
    }
    reps
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

/// **The 2-cyclotomic cosets mod `n`**: orbits of `s ↦ 2s` on `Z/n`.
///
/// These classify every Frobenius-stable `F_2`-subspace of `F_{2^n}`.
/// A subspace closed under squaring is an `F_2[x]/(x^n − 1)`-submodule,
/// i.e. a binary cyclic code of length `n`, i.e. a divisor of
/// `x^n − 1`; and the irreducible factors of `x^n − 1` correspond one
/// for one with these cosets, a coset of size `d` giving a factor of
/// degree `d`.
///
/// So the *available dimensions* of an invariant factor base are
/// exactly the subset sums of the coset sizes — see
/// [`available_subspace_dimensions`].
pub fn cyclotomic_cosets(n: u32) -> Vec<Vec<u32>> {
    let mut seen = vec![false; n as usize];
    let mut out = Vec::new();
    for s in 0..n {
        if seen[s as usize] {
            continue;
        }
        let mut coset = Vec::new();
        let mut x = s;
        while !seen[x as usize] {
            seen[x as usize] = true;
            coset.push(x);
            x = (x * 2) % n;
        }
        coset.sort_unstable();
        out.push(coset);
    }
    out
}

/// **Every dimension a Frobenius-invariant factor base can have** over
/// `F_{2^n}`, in increasing order.
///
/// The subset sums of the cyclotomic coset sizes.  This is a complete
/// classification, not a heuristic: a stable subspace *is* a divisor of
/// `x^n − 1`, so no other dimension is achievable.
///
/// It is also the answer to "which `n` are worth attacking this way".
/// When `2` is primitive mod `n` — `n = 131` and `n = 163`, the sizes
/// that matter — there are only two cosets, `{0}` and everything else,
/// so the only dimensions are `0, 1, n − 1, n`: nothing usable.  When
/// `n` has many small cosets the sizes are dense, and the factor base
/// can be tuned to whatever `|F| = 2^dim` the sizing condition wants.
pub fn available_subspace_dimensions(n: u32) -> Vec<u32> {
    let mut reach = std::collections::BTreeSet::from([0u32]);
    for c in cyclotomic_cosets(n) {
        let d = c.len() as u32;
        for r in reach.clone() {
            reach.insert(r + d);
        }
    }
    reach.into_iter().collect()
}

/// Product of two `F_2[x]` polynomials given as bitmasks, without
/// reduction.  `None` if the product would not fit a `u64`.
fn poly_mul_full(a: u64, b: u64) -> Option<u64> {
    let (da, db) = (poly_deg(a), poly_deg(b));
    if let (Some(da), Some(db)) = (da, db) {
        if da + db >= 64 {
            return None;
        }
    }
    let mut acc = 0u64;
    let mut b = b;
    let mut shift = 0;
    while b != 0 {
        if b & 1 == 1 {
            acc ^= a << shift;
        }
        b >>= 1;
        shift += 1;
    }
    Some(acc)
}

/// **Every irreducible factor of `x^n − 1`** over `F_2`, including
/// `x + 1`, as bitmasks.
///
/// [`factor_x_n_minus_1`] returns only the non-trivial factors; a
/// divisor may use `x + 1` too, so the divisor-based constructions take
/// their indices into this list.
pub fn all_factors_of_x_n_minus_1(n: u32) -> Vec<u64> {
    let mut out = vec![0b11u64]; // x + 1
                                 // The degrees that occur are exactly the cyclotomic coset sizes;
                                 // the size-1 coset is {0}, already covered by x + 1 above.
    let mut degrees: Vec<u32> = cyclotomic_cosets(n)
        .iter()
        .map(|c| c.len() as u32)
        .filter(|d| *d > 1)
        .collect();
    degrees.sort_unstable();
    degrees.dedup();
    for d in degrees {
        if d > 24 {
            break; // the scan below is 2^d; beyond this it is not worth it
        }
        out.extend(irreducible_factors_of_degree(n, d));
    }
    out
}

/// Every irreducible factor of `x^n − 1` of degree exactly `d`.
///
/// Scans the `2^d` monic polynomials of that degree.  Used by
/// [`all_factors_of_x_n_minus_1`]; [`factor_x_n_minus_1`] is the
/// `d = ord_n(2)` case of the same search.
fn irreducible_factors_of_degree(n: u32, d: u32) -> Vec<u64> {
    let mut out = Vec::new();
    let hi = 1u64 << d;
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

/// An `F_2`-basis of the Frobenius-invariant subspace belonging to the
/// divisor `Π factors[i]`, `i ∈ indices` — dimension = total degree.
///
/// The single-factor case ([`invariant_subspace_basis`]) is stuck at
/// `dim = ord_n(2)`, which is the only dimension one irreducible factor
/// can give.  Products give every dimension in
/// [`available_subspace_dimensions`], so the factor base can be sized
/// to the instance instead of the instance to the factor base.
pub fn subspace_basis_for_divisor(
    n: u32,
    indices: &[usize],
    irr: &IrreduciblePoly,
) -> Option<Vec<F2mElement>> {
    let factors = all_factors_of_x_n_minus_1(n);
    let mut f = 1u64;
    for &i in indices {
        f = poly_mul_full(f, *factors.get(i)?)?;
    }
    let deg = poly_deg(f)?;
    if deg == 0 || deg >= n {
        return None;
    }
    let exps: Vec<u32> = (0..=deg).filter(|k| (f >> k) & 1 == 1).collect();
    let basis = linearised_kernel_basis(&exps, n, irr);
    if basis.len() != deg as usize {
        return None;
    }
    Some(basis)
}

/// **Build a Frobenius-invariant factor base from a divisor** of
/// `x^n − 1`, rather than from a single irreducible factor.
///
/// `indices` select factors from [`all_factors_of_x_n_minus_1`]; the
/// subspace has dimension equal to their total degree, so `|F| ≈ 2^dim`
/// is tunable.  That matters because the summand count `m` a
/// decomposition needs falls as `|F|` grows — `m ≈ n/dim` — and `m ≥ 3`
/// is what forces the chained system and its `(m − 2)·n` extra
/// unknowns.  A large enough invariant subspace buys `m = 2` and skips
/// the chaining entirely.
pub fn build_frobenius_factor_base_from_divisor(
    kc: &KoblitzCurve,
    indices: &[usize],
) -> Option<FrobeniusFactorBase> {
    let subspace_basis = subspace_basis_for_divisor(kc.n, indices, &kc.curve.irreducible)?;
    let ell = subspace_basis.len() as u32;
    let factors = all_factors_of_x_n_minus_1(kc.n);
    let mut f_j = 1u64;
    for &i in indices {
        f_j = poly_mul_full(f_j, *factors.get(i)?)?;
    }
    let exps: Vec<u32> = (0..=ell).filter(|k| (f_j >> k) & 1 == 1).collect();
    finish_factor_base(kc, ell, f_j, exps, subspace_basis)
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
    if subspace_basis.len() != ell as usize {
        // Kernel dimension must equal deg f_j; anything else means the
        // factor did not divide x^n − 1 after all.
        return None;
    }
    finish_factor_base(kc, ell, f_j, exps, subspace_basis)
}

/// Shared tail of the factor-base constructors: span the subspace,
/// collect the points over it, and walk the `π`-orbits.
fn finish_factor_base(
    kc: &KoblitzCurve,
    ell: u32,
    f_j: u64,
    exps: Vec<u32>,
    subspace_basis: Vec<F2mElement>,
) -> Option<FrobeniusFactorBase> {
    let subspace = span_f2(&subspace_basis, kc.n);
    if subspace.len() != (1usize << ell) {
        return None;
    }

    finish_factor_base_domain(
        kc,
        ell,
        f_j,
        exps,
        subspace_basis,
        subspace,
        FactorBaseDomain::LinearSubspace,
    )
}

/// Materialise the union of all Frobenius translates of a seed space.
/// This permits small invariant *sets* even when no small invariant
/// linear subspace exists. Construction is charged: at most n*2^seed_dim
/// abscissae before deduplication; SAT uses an explicit domain trie.
/// This is a toy construction, not a claim of scalable relation solving.
pub fn build_frobenius_union_factor_base(
    kc: &KoblitzCurve,
    seed_basis: &[F2mElement],
) -> Option<FrobeniusFactorBase> {
    if seed_basis.is_empty() || seed_basis.len() > 12 {
        return None;
    }
    let seed = span_f2(seed_basis, kc.n);
    let unique: std::collections::HashSet<_> = seed.iter().map(|x| x.to_biguint()).collect();
    if unique.len() != (1usize << seed_basis.len()) {
        return None;
    }
    let mut xs = std::collections::BTreeMap::new();
    for mut x in seed {
        for _ in 0..kc.n {
            xs.entry(x.to_biguint()).or_insert_with(|| x.clone());
            x = x.square(&kc.curve.irreducible);
        }
    }
    let ambient: Vec<_> = (0..kc.n)
        .map(|i| F2mElement::from_bit_positions(&[i], kc.n))
        .collect();
    finish_factor_base_domain(
        kc,
        kc.n,
        0,
        Vec::new(),
        ambient,
        xs.into_values().collect(),
        FactorBaseDomain::FrobeniusUnion {
            seed_dimension: seed_basis.len(),
        },
    )
}

/// Build a factor base from explicit abscissa-orbit representatives.
///
/// Each supplied coordinate is closed under the `2`-power Frobenius;
/// rational points above the resulting coordinates are materialised and
/// constrained exactly by the SAT domain trie. This is the constructor
/// used for factor bases selected by a finite orbit search, where the
/// union need not be the Frobenius closure of a small linear seed space.
pub fn build_explicit_frobenius_orbit_factor_base(
    kc: &KoblitzCurve,
    representatives: &[F2mElement],
) -> Option<FrobeniusFactorBase> {
    if representatives.is_empty() {
        return None;
    }
    let mut xs = std::collections::BTreeMap::new();
    for representative in representatives {
        let mut x = F2mElement::from_biguint(&representative.to_biguint(), kc.n);
        for _ in 0..kc.n {
            xs.entry(x.to_biguint()).or_insert_with(|| x.clone());
            x = x.square(&kc.curve.irreducible);
        }
    }
    let ambient = (0..kc.n)
        .map(|i| F2mElement::from_bit_positions(&[i], kc.n))
        .collect();
    finish_factor_base_domain(
        kc,
        kc.n,
        0,
        Vec::new(),
        ambient,
        xs.into_values().collect(),
        FactorBaseDomain::ExplicitFrobeniusOrbits {
            representatives: representatives.len(),
        },
    )
}

/// **Keep only the listed signed orbits** of a factor base.
///
/// The retained set is still closed under Frobenius and negation, so
/// every relation identity survives; what changes is the relation
/// column count (down to `keep.len()`) and, for the enumeration and
/// pair-table oracles, the per-target cost.  A subspace base stays in
/// its subspace coordinates ([`FactorBaseDomain::SubspaceSubset`]); a
/// nonlinear base becomes an explicit orbit union.  Returns `None` for
/// an empty or out-of-range selection.
pub fn restrict_factor_base_to_orbits(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    keep: &[usize],
) -> Option<FrobeniusFactorBase> {
    if keep.is_empty() || keep.iter().any(|&o| o >= fb.signed_orbits.len()) {
        return None;
    }
    let mut keep_sorted = keep.to_vec();
    keep_sorted.sort_unstable();
    keep_sorted.dedup();
    let mut xs = std::collections::BTreeMap::new();
    for &o in &keep_sorted {
        for &i in &fb.signed_orbits[o] {
            if let BinaryPoint::Affine { x, .. } = &fb.points[i] {
                xs.entry(x.to_biguint()).or_insert_with(|| x.clone());
            }
        }
    }
    if fb.uses_ambient_basis() {
        // One representative per retained orbit; the constructor closes
        // it under Frobenius again, which is a no-op here.
        let representatives: Vec<F2mElement> = keep_sorted
            .iter()
            .map(|&o| {
                let i = fb.signed_orbits[o][0];
                match &fb.points[i] {
                    BinaryPoint::Affine { x, .. } => x.clone(),
                    BinaryPoint::Infinity => unreachable!("factor bases hold affine points"),
                }
            })
            .collect();
        return build_explicit_frobenius_orbit_factor_base(kc, &representatives);
    }
    finish_factor_base_domain(
        kc,
        fb.ell,
        fb.f_j,
        fb.linearised_exponents.clone(),
        fb.subspace_basis.clone(),
        xs.into_values().collect(),
        FactorBaseDomain::SubspaceSubset {
            retained_orbits: keep_sorted.len(),
        },
    )
}

/// B union (B+T), omitting infinity, where T=(0,1) on K_a.
/// Since [h]T=O for even h, this changes the available torsion lifts
/// without adding projected points. Translation commutes with Frobenius.
/// Its x-coordinate action away from T is x -> 1/x.
pub fn saturate_factor_base_two_torsion(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
) -> Option<FrobeniusFactorBase> {
    let t = points_with_x(&kc.curve, &F2mElement::zero(kc.n))
        .into_iter()
        .next()?;
    if kc.mul(&t, &kc.cofactor) != BinaryPoint::Infinity {
        return None;
    }
    let mut xs = std::collections::BTreeMap::new();
    for p in &fb.points {
        for q in [p.clone(), kc.add(p, &t)] {
            if let BinaryPoint::Affine { x, .. } = q {
                xs.entry(x.to_biguint()).or_insert(x);
            }
        }
    }
    let ambient = (0..kc.n)
        .map(|i| F2mElement::from_bit_positions(&[i], kc.n))
        .collect();
    finish_factor_base_domain(
        kc,
        kc.n,
        0,
        Vec::new(),
        ambient,
        xs.into_values().collect(),
        FactorBaseDomain::TwoTorsionSaturation,
    )
}

fn finish_factor_base_domain(
    kc: &KoblitzCurve,
    ell: u32,
    f_j: u64,
    exps: Vec<u32>,
    subspace_basis: Vec<F2mElement>,
    subspace: Vec<F2mElement>,
    domain: FactorBaseDomain,
) -> Option<FrobeniusFactorBase> {
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

    let mut signed_orbit_of = vec![(usize::MAX, 0, false); points.len()];
    let mut signed_orbits = Vec::new();
    for start in 0..points.len() {
        if signed_orbit_of[start].0 != usize::MAX {
            continue;
        }
        let signed_orbit = signed_orbits.len();
        let mut members = Vec::new();
        let mut current = points[start].clone();
        for k in 0..kc.n {
            for (negated, point) in [(false, current.clone()), (true, point_neg(&current))] {
                let index = *index_of.get(&point_key(&point))?;
                if signed_orbit_of[index].0 == usize::MAX {
                    signed_orbit_of[index] = (signed_orbit, k, negated);
                    members.push(index);
                } else if signed_orbit_of[index].0 != signed_orbit {
                    return None;
                }
            }
            current = kc.frobenius(&current);
        }
        if current != points[start] {
            return None;
        }
        signed_orbits.push(members);
    }
    if signed_orbits.iter().map(Vec::len).sum::<usize>() != points.len() {
        return None;
    }

    Some(FrobeniusFactorBase {
        domain,
        ell,
        f_j,
        linearised_exponents: exps,
        subspace,
        subspace_basis,
        points,
        orbits,
        orbit_of,
        signed_orbits,
        signed_orbit_of,
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
    /// The decomposition, as `(signed orbit, k)` pairs. Consult
    /// [`Self::summand_negated`] for the sign of each summand.
    pub summands: Vec<(usize, u32)>,
    /// Whether each summand is `-π^k(rep)` rather than `π^k(rep)`.
    pub summand_negated: Vec<bool>,
    /// Dense row over the signed-orbit unknowns: entry `o` is
    /// `Σ_i (-1)^negated_i λ^{k_i}` over its summands, reduced mod `r`.
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
    /// Gröbner basis. Subexponential cost requires additional, unproved
    /// assumptions about solving these systems. See
    /// [`crate::cryptanalysis::koblitz_groebner`].
    Groebner,
    /// Exhaustive ordered-tuple search over the materialised factor
    /// base.  Costs `|F|^{m−1}` group operations per target whether or
    /// not a decomposition exists; kept as the reference oracle the
    /// algebraic one is tested against.
    Enumerate,
    /// The same Semaev system, handed to the CDCL SAT solver instead of
    /// a Gröbner engine
    /// ([`crate::cryptanalysis::semaev_sat::encode_boolean_system_with`]).
    /// The default keeps descended parity equations as native XOR rows;
    /// the CNF path remains available as a controlled comparison.
    Sat,
    /// **Meet in the middle** over a precomputed table of all pair sums
    /// `P_i + P_j` ([`PairSumTable`]).  `m = 2` is one lookup, `m = 3`
    /// is `|F|` lookups, `m = 4` is `|F|²` lookups — against
    /// `|F|^{m−1}` group operations for [`Self::Enumerate`].  Exact and
    /// complete like enumeration; costs `|F|²` memory once per run.
    PairTable,
}

/// Packed identity of a point for hashing and sorting, in one `u64`:
/// `0` for `O`, otherwise `2(x + 1) + s` where the sign bit `s`
/// separates `P = (x, y)` from `−P = (x, x + y)` by comparing `y` with
/// `x + y` as integers (they coincide only at `x = 0`, where `P = −P`).
/// Exact for every `n ≤ 62`.
pub fn pack_point(p: &BinaryPoint) -> u64 {
    match p {
        BinaryPoint::Infinity => 0,
        BinaryPoint::Affine { x, y } => {
            let xb = x.raw_bits().first().copied().unwrap_or(0);
            let yb = y.raw_bits().first().copied().unwrap_or(0);
            let sign = u64::from(yb > (xb ^ yb));
            ((xb + 1) << 1) | sign
        }
    }
}

/// **Every pair sum of a factor base**, `P_i + P_j` for `i ≤ j`, sorted
/// by packed point identity so a sum can be looked up by binary search.
///
/// This is the birthday-paradox half of a decomposition search: with
/// the table built once, "is `R` a sum of two base points" is one
/// lookup, "of three" is `|F|` lookups of `R − P_k`, and "of four" is
/// `|F|²` lookups of `R − (P_k + P_l)` walked off the table itself.
/// The factor-base search uses the same table to count *every*
/// witness of a target, which is what exact yield needs.
#[derive(Clone, Debug)]
pub struct PairSumTable {
    /// `(packed sum, i, j)` with `i ≤ j`, sorted by the packed sum.
    entries: Vec<(u64, u32, u32)>,
}

impl PairSumTable {
    /// Build the table with `|F|(|F|+1)/2` point additions, in parallel
    /// over the first index; 16 bytes per entry.  Returns `None` when the
    /// field is too wide to pack a point into a `u64`.
    pub fn build(kc: &KoblitzCurve, fb: &FrobeniusFactorBase) -> Option<Self> {
        if kc.n > 62 || fb.points.len() > u32::MAX as usize {
            return None;
        }
        let mut entries: Vec<(u64, u32, u32)> = (0..fb.points.len())
            .into_par_iter()
            .flat_map_iter(|i| {
                let kc = kc;
                let fb = fb;
                (i..fb.points.len()).map(move |j| {
                    let sum = kc.add(&fb.points[i], &fb.points[j]);
                    (pack_point(&sum), i as u32, j as u32)
                })
            })
            .collect();
        entries.par_sort_unstable();
        Some(Self { entries })
    }

    /// Number of stored pair sums (with multiplicity).
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the table is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// All `(i, j)` with `P_i + P_j` equal to the packed point.
    pub fn lookup(&self, key: u64) -> &[(u64, u32, u32)] {
        let start = self.entries.partition_point(|e| e.0 < key);
        let end = start + self.entries[start..].partition_point(|e| e.0 == key);
        &self.entries[start..end]
    }

    /// **Decompose** `target` into exactly `m ∈ {2, 3, 4}` base points,
    /// returning sorted indices, or `None` when no decomposition exists
    /// (a complete search) or `m` is unsupported.  The returned sum is
    /// re-checked in the group before it is handed back.
    pub fn decompose(
        &self,
        kc: &KoblitzCurve,
        fb: &FrobeniusFactorBase,
        target: &BinaryPoint,
        m: usize,
    ) -> Option<Vec<usize>> {
        let mut found: Option<Vec<usize>> = None;
        self.witnesses(kc, fb, target, m, &mut |witness| {
            found = Some(witness.to_vec());
            false
        });
        let idxs = found?;
        let sum = idxs
            .iter()
            .fold(BinaryPoint::Infinity, |s, &i| kc.add(&s, &fb.points[i]));
        (sum == *target).then_some(idxs)
    }

    /// **Enumerate every sorted witness** `i_1 ≤ … ≤ i_m` with
    /// `Σ P_{i_k} = target`, calling `sink` for each; the sink returns
    /// `true` to continue and `false` to stop.  Supports `m ∈ {2, 3, 4}`.
    pub fn witnesses(
        &self,
        kc: &KoblitzCurve,
        fb: &FrobeniusFactorBase,
        target: &BinaryPoint,
        m: usize,
        sink: &mut dyn FnMut(&[usize]) -> bool,
    ) {
        match m {
            2 => {
                for &(_, i, j) in self.lookup(pack_point(target)) {
                    if !sink(&[i as usize, j as usize]) {
                        return;
                    }
                }
            }
            3 => {
                for (k, p) in fb.points.iter().enumerate() {
                    let rest = kc.add(target, &point_neg(p));
                    for &(_, i, j) in self.lookup(pack_point(&rest)) {
                        if j as usize <= k && !sink(&[i as usize, j as usize, k]) {
                            return;
                        }
                    }
                }
            }
            4 => {
                // Walk the table as the second half: R − (P_k + P_l).
                let mut last_key: Option<u64> = None;
                let mut rest = BinaryPoint::Infinity;
                for &(key, k, l) in &self.entries {
                    if last_key != Some(key) {
                        let pair = kc.add(&fb.points[k as usize], &fb.points[l as usize]);
                        rest = kc.add(target, &point_neg(&pair));
                        last_key = Some(key);
                    }
                    for &(_, i, j) in self.lookup(pack_point(&rest)) {
                        if j <= k && !sink(&[i as usize, j as usize, k as usize, l as usize]) {
                            return;
                        }
                    }
                }
            }
            _ => {}
        }
    }
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
    /// Implied Macaulay rows added to the CNF (see
    /// [`sat_decompose`]'s `macaulay_degree`).
    pub implied_rows: usize,
    /// Conflicts across every solve call — a machine-independent
    /// measure of search effort.
    pub conflicts: u64,
}

/// Encode a finite set of coordinates without enumerating its complement.
/// Missing branches of the binary trie become forbidden-prefix clauses.
/// No auxiliary variables are needed; at most O(ell * |codes|) clauses.
fn add_coordinate_domain(
    solver: &mut crate::cryptanalysis::sat::Solver,
    offset: usize,
    ell: usize,
    codes: &[u64],
) {
    fn visit(
        solver: &mut crate::cryptanalysis::sat::Solver,
        offset: usize,
        bit: usize,
        codes: &[u64],
        prefix: &mut Vec<i32>,
    ) {
        if codes.is_empty() {
            solver.add_clause(prefix.clone());
            return;
        }
        if bit == 0 {
            return;
        }
        let b = bit - 1;
        let split = codes.partition_point(|code| (code >> b) & 1 == 0);
        let lit = (offset + b + 1) as i32;
        prefix.push(lit); // forbid prefix with this bit zero
        visit(solver, offset, b, &codes[..split], prefix);
        *prefix.last_mut().unwrap() = -lit;
        visit(solver, offset, b, &codes[split..], prefix);
        prefix.pop();
    }
    visit(solver, offset, ell, codes, &mut Vec::new());
}

/// Absolute trace F_(2^n) -> F_2 in the configured field representation.
fn absolute_trace_bit(x: &F2mElement, n: u32, irr: &IrreduciblePoly) -> bool {
    let mut t = F2mElement::zero(n);
    let mut power = x.clone();
    for _ in 0..n {
        t = t.add(&power);
        power = power.square(irr);
    }
    debug_assert!(t.is_zero() || t == F2mElement::one(n));
    !t.is_zero()
}

/// SAT controls for a paired comparison on the same decomposition system.
#[derive(Clone, Copy, Debug)]
pub struct SatDecompositionOptions {
    pub encoding: XorEncoding,
    /// Prioritize summand coordinates; all other variables remain eligible.
    pub branch_on_summands: bool,
    /// Exclude summand coordinates without a rational factor-base point.
    /// Costs O(m * ell * 2^ell) preprocessing; for materialised bases only.
    pub restrict_to_factor_base: bool,
    /// Add the group homomorphism Tr(x(P) + a) as one linear XOR row.
    pub trace_constraint: bool,
    /// Cumulative conflict cap across every model of one target.
    pub conflict_budget: u64,
    /// Order the summands lexicographically by subspace code,
    /// `code(x_1) ≤ code(x_2) ≤ … ≤ code(x_m)`.  The summands are
    /// interchangeable, so if any decomposition exists a sorted one
    /// does (for a chained system the intermediates simply follow the
    /// new order), and the constraint is exact: no verdict changes.
    ///
    /// **Off by default, because it does not pay here.**  Measured with
    /// `koblitz_sat_symmetry_ablation` (conflict budget 400 000, verdicts
    /// cross-checked against search on every target):
    ///
    /// | instance | targets | conflicts off | conflicts on |
    /// |:---------|--------:|--------------:|-------------:|
    /// | K_0/2^9, m = 3 | 8 | 5 512 | 6 412 |
    /// | K_1/2^9, m = 2 | 16 | 428 | 823 |
    /// | K_0/2^7, m = 2 | 16 | 70 | 79 |
    /// | K_1/2^15, m = 3 | 8 | 0 | 0 |
    ///
    /// The degree-2 Macaulay rows and the trace row already leave these
    /// systems almost propagation-closed, and the `ℓ − 1` "equal so far"
    /// auxiliaries per adjacent pair add decisions without removing
    /// search.  Kept as a control for larger chained instances.
    pub symmetry_breaking: bool,
}

impl Default for SatDecompositionOptions {
    fn default() -> Self {
        Self {
            encoding: XorEncoding::Native,
            branch_on_summands: false,
            restrict_to_factor_base: false,
            trace_constraint: true,
            conflict_budget: u64::MAX,
            symmetry_breaking: false,
        }
    }
}

/// Install `code(a) ≤ code(b)` (unsigned, bit `width − 1` most
/// significant) between two blocks of `width` solver variables that
/// start at 1-indexed `a_first` and `b_first`.
///
/// One "equal so far" auxiliary per bit position below the top; the
/// auxiliaries are only ever *forced true* by equality, never forced
/// false, so the constraint admits exactly the assignments with
/// `a ≤ b` and nothing else is excluded.
fn add_lex_leq(
    solver: &mut crate::cryptanalysis::sat::Solver,
    a_first: u32,
    b_first: u32,
    width: usize,
) {
    if width == 0 {
        return;
    }
    let aux = solver.add_vars(width as u32 - 1);
    let aux: Vec<i32> = aux.map(|v| v as i32).collect();
    let a = |t: usize| (a_first + t as u32) as i32;
    let b = |t: usize| (b_first + t as u32) as i32;
    // Top bit: a_top → b_top.
    let top = width - 1;
    solver.add_clause(vec![-a(top), b(top)]);
    // e_t means "a and b agree on every bit above t".
    for t in (0..top).rev() {
        let e_t = aux[t];
        // e_t is forced by equality at bit t+1 given e_{t+1} (or the
        // top level, which is unconditional).
        let above: Option<i32> = if t + 1 == top { None } else { Some(aux[t + 1]) };
        for (la, lb) in [(a(t + 1), b(t + 1)), (-a(t + 1), -b(t + 1))] {
            let mut clause = vec![la, lb, e_t];
            if let Some(e_above) = above {
                clause.push(-e_above);
            }
            solver.add_clause(clause);
        }
        // Under e_t: a_t → b_t.
        solver.add_clause(vec![-e_t, -a(t), b(t)]);
    }
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
    macaulay_degree: Option<u32>,
) -> (Option<Vec<usize>>, SatDecompositionStats) {
    sat_decompose_with(
        kc,
        fb,
        index_of,
        st,
        target,
        m,
        max_models,
        macaulay_degree,
        SatDecompositionOptions::default(),
    )
}

/// Configurable SAT decomposition, retaining CNF as an explicit control.
/// Pass `macaulay_degree = None` for a SAT-only run without F4 preprocessing.
/// Three-point unions use the wide S4 encoder (no Macaulay preprocessing);
/// its optional group trace row is native even with CNF polynomial equations.
pub fn sat_decompose_with(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    max_models: usize,
    macaulay_degree: Option<u32>,
    options: SatDecompositionOptions,
) -> (Option<Vec<usize>>, SatDecompositionStats) {
    if m == 3 && fb.uses_ambient_basis() {
        return sat_decompose_union_s4(kc, fb, index_of, target, max_models, options);
    }
    let mut stats = SatDecompositionStats::default();
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => return (None, stats),
    };
    let sys = match build_decomposition_system(&fb.subspace_basis, &x_r, &kc.curve.b, m, st) {
        Some(sys) => sys,
        None => return (None, stats),
    };

    // Algebraic preprocessing: hand the solver the degree-`D` Macaulay
    // consequences of the system alongside the system itself.  Each row
    // is an `F_2`-combination of multiples of the equations, so it is
    // implied and adding it cannot change the answer — but it can save
    // the search from rediscovering it.  Measured on the refuting
    // instances this cuts conflicts several-fold where the system is
    // overdetermined, and does nothing where it is not.
    //
    // The rows are **added, never substituted**.  `matrix_f4_f2` skips
    // input polynomials of degree above `D`, so replacing the system
    // with its degree-2 rows silently drops every cubic equation of a
    // chained (`m ≥ 3`) system — which turns UNSAT into SAT.
    let mut equations = sys.equations.clone();
    if options.trace_constraint {
        use crate::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
        // For a rational point, delta(P)=Tr(x(P)+a), delta(O)=0,
        // is a homomorphism to F_2 (Kosters--Yeo, arXiv:1503.08001).
        // Thus sum Tr(x_i) = Tr(x_R) + (m+1)Tr(a). This is a
        // necessary group condition, not an assumption about polynomial roots.
        let rhs = absolute_trace_bit(&x_r, kc.n, &kc.curve.irreducible)
            ^ (m % 2 == 0 && kc.n % 2 == 1 && kc.a == 1);
        let mut terms = Vec::new();
        for (j, basis_element) in fb.subspace_basis.iter().enumerate() {
            if absolute_trace_bit(basis_element, kc.n, &kc.curve.irreducible) {
                for i in 0..m {
                    terms.push(F2BoolMono::var((i * fb.subspace_basis.len() + j) as u32));
                }
            }
        }
        if rhs {
            terms.push(F2BoolMono::one());
        }
        equations.push(F2BoolPoly::from_monos(terms, sys.n_vars));
    }
    if let Some(d) = macaulay_degree {
        if let Some(rows) = matrix_f4_f2(&sys.equations, sys.n_vars, d) {
            stats.implied_rows = rows.len();
            equations.extend(rows);
        }
    }

    // Enumerate models incrementally: encode once, then add a blocking
    // clause per rejected root and re-solve.  Re-encoding from scratch
    // each time costs a full Tseitin/XOR build per model examined, and
    // with `max_models` in the dozens that dominated the loop.
    let mut enc = encode_boolean_system_with(sys.n_vars, &equations, &[], options.encoding);
    enc.solver.conflict_budget = options.conflict_budget;
    if options.symmetry_breaking {
        let ell = fb.subspace_basis.len();
        for i in 1..m {
            add_lex_leq(
                &mut enc.solver,
                ((i - 1) * ell + 1) as u32,
                (i * ell + 1) as u32,
                ell,
            );
        }
    }
    if options.branch_on_summands {
        enc.solver
            .set_branch_priority(&(1..=(m * fb.subspace_basis.len()) as u32).collect::<Vec<_>>());
    }
    if options.restrict_to_factor_base || fb.domain != FactorBaseDomain::LinearSubspace {
        let ell = fb.subspace_basis.len();
        let legal_x: std::collections::HashSet<BigUint> = fb
            .points
            .iter()
            .filter_map(|p| match p {
                BinaryPoint::Affine { x, .. } => Some(x.to_biguint()),
                BinaryPoint::Infinity => None,
            })
            .collect();
        let mut codes: Vec<u64> = if fb.uses_ambient_basis() {
            // Union construction uses the full polynomial field basis.
            legal_x
                .iter()
                .map(|x| x.to_u64_digits().first().copied().unwrap_or(0))
                .collect()
        } else {
            (0..(1u64 << ell))
                .filter(|&code| {
                    legal_x.contains(
                        &sys.summand_x(&fb.subspace_basis, code, 0, kc.n)
                            .to_biguint(),
                    )
                })
                .collect()
        };
        codes.sort_unstable();
        codes.dedup();
        for i in 0..m {
            add_coordinate_domain(&mut enc.solver, i * ell, ell, &codes);
        }
    }
    let mut blocked: Vec<u64> = Vec::new();
    loop {
        stats.solver_calls += 1;
        if enc.trivially_unsat {
            stats.refuted = true;
            return (None, stats);
        }
        let outcome = enc.solver.solve();
        // The solver counter is cumulative across incremental calls.
        stats.conflicts = enc.solver.conflicts();
        match outcome {
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
                    // Fail closed: a bad model invalidates this attempt.
                    // Blocking it and later reporting UNSAT would hide a bug.
                    stats.spurious += 1;
                    stats.exhausted = true;
                    return (None, stats);
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
                // Forbid this assignment and search again.
                enc.solver.reset_search();
                let clause: Vec<i32> = (0..sys.n_vars)
                    .map(|i| {
                        let lit = (i + 1) as i32;
                        if (root >> i) & 1 == 1 {
                            -lit
                        } else {
                            lit
                        }
                    })
                    .collect();
                enc.solver.add_clause(clause);
            }
        }
    }
}

/// The wide symmetrised S4 encoder already supports more than 64
/// problem variables. Reuse it for three-point union decompositions;
/// the finite coordinate domain keeps each X inside the union.
fn sat_decompose_union_s4(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    target: &BinaryPoint,
    max_models: usize,
    options: SatDecompositionOptions,
) -> (Option<Vec<usize>>, SatDecompositionStats) {
    use crate::cryptanalysis::semaev_sat::encode_semaev_s4;
    let mut stats = SatDecompositionStats::default();
    let BinaryPoint::Affine { x: x_r, .. } = target else {
        stats.exhausted = true;
        return (None, stats);
    };
    let mut enc = encode_semaev_s4(
        kc.n,
        kc.n,
        &kc.curve.irreducible,
        &kc.curve.b,
        x_r,
        options.encoding,
    );
    enc.solver.conflict_budget = options.conflict_budget;
    let mut codes: Vec<_> = fb
        .points
        .iter()
        .filter_map(|p| match p {
            BinaryPoint::Affine { x, .. } => Some(x.raw_bits().first().copied().unwrap_or(0)),
            BinaryPoint::Infinity => None,
        })
        .collect();
    codes.sort_unstable();
    codes.dedup();
    for i in 0..3 {
        add_coordinate_domain(&mut enc.solver, i * kc.n as usize, kc.n as usize, &codes);
    }
    if options.trace_constraint {
        let mut vars = Vec::new();
        for j in 0..kc.n {
            if absolute_trace_bit(
                &F2mElement::from_bit_positions(&[j], kc.n),
                kc.n,
                &kc.curve.irreducible,
            ) {
                for i in 0..3 {
                    vars.push(i * kc.n + j + 1);
                }
            }
        }
        let rhs = absolute_trace_bit(x_r, kc.n, &kc.curve.irreducible);
        // This optional group row is native, including when the
        // polynomial equations use the CNF control encoding.
        enc.solver.add_xor(&vars, rhs);
    }
    loop {
        stats.solver_calls += 1;
        let result = enc.solver.solve();
        stats.conflicts = enc.solver.conflicts();
        match result {
            SolveResult::Unsat => {
                stats.refuted = true;
                return (None, stats);
            }
            SolveResult::Unknown => {
                stats.exhausted = true;
                return (None, stats);
            }
            SolveResult::Sat => {
                let xs = enc.decode();
                if !crate::cryptanalysis::binary_semaev::binary_semaev_s4(
                    &xs[0],
                    &xs[1],
                    &xs[2],
                    x_r,
                    &kc.curve.b,
                    &kc.curve.irreducible,
                )
                .is_zero()
                {
                    stats.spurious += 1;
                    stats.exhausted = true;
                    return (None, stats);
                }
                stats.models += 1;
                if let Some(ids) = lift_candidate(kc, fb, index_of, &xs, target) {
                    return (Some(ids), stats);
                }
                if stats.models >= max_models {
                    stats.exhausted = true;
                    return (None, stats);
                }
                let model = enc.solver.model();
                enc.solver.reset_search();
                // The x coordinates determine all S4 auxiliaries, so
                // block the x tuple, including every failed rational lift.
                let clause = (0..3 * kc.n as usize)
                    .map(|i| {
                        let lit = (i + 1) as i32;
                        if model[i] {
                            -lit
                        } else {
                            lit
                        }
                    })
                    .collect();
                enc.solver.add_clause(clause);
            }
        }
    }
}

/// Turn a decomposition into a relation row over the orbit unknowns.
///
/// With `x_o := log_G ([h]·rep_o)` and
/// `P_i = (-1)^{s_i}π^{k_i}(rep_{o_i})`,
/// multiplying `R = Σ_i P_i` by the cofactor `h` gives
///
/// ```text
///     h·a + h·b·d  ≡  Σ_i (-1)^{s_i} λ^{k_i} · x_{o_i}   (mod r),
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
    relation_from_decomposition_with_mode(kc, fb, idxs, coef_a, coef_b, true, None)
}

fn relation_from_decomposition_with_mode(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    idxs: &[usize],
    coef_a: &BigUint,
    coef_b: &BigUint,
    collapse_negation: bool,
    projected_orbits: Option<&ProjectedSignedOrbitMap>,
) -> KoblitzRelation {
    let r = &kc.subgroup_order;
    let unknowns = projected_orbits.map_or_else(
        || {
            if collapse_negation {
                fb.unknowns()
            } else {
                fb.orbits.len()
            }
        },
        |projected| projected.representatives.len(),
    );
    let mut row = vec![BigUint::zero(); unknowns];
    let mut summands = Vec::with_capacity(idxs.len());
    let mut summand_negated = Vec::with_capacity(idxs.len());
    for &i in idxs {
        let location = if let Some(projected) = projected_orbits {
            projected.orbit_of[i]
        } else if collapse_negation {
            Some(fb.signed_orbit_of[i])
        } else {
            let (o, k) = fb.orbit_of[i];
            Some((o, k, false))
        };
        let Some((o, k, negated)) = location else {
            // [h]P = O, so this summand contributes zero after the
            // relation is projected into the prime-order subgroup.
            continue;
        };
        let mut coeff = kc.lambda.modpow(&BigUint::from(k), r);
        if negated && !coeff.is_zero() {
            coeff = r - coeff;
        }
        row[o] = (&row[o] + coeff) % r;
        summands.push((o, k));
        summand_negated.push(negated);
    }
    KoblitzRelation {
        coef_a: coef_a.clone(),
        coef_b: coef_b.clone(),
        summands,
        summand_negated,
        row,
    }
}

// ── Driver ─────────────────────────────────────────────────────────

#[derive(Clone, Debug)]
struct ProjectedSignedOrbitMap {
    /// For each factor-base point P, the signed Frobenius location of
    /// [h]P. `None` means [h]P = O and therefore contributes no column.
    orbit_of: Vec<Option<(usize, u32, bool)>>,
    representatives: Vec<BinaryPoint>,
}

/// Build the quotient factor-base columns after public cofactor projection.
///
/// The ordinary factor-base orbit table can contain several columns whose
/// cofactor projections are the same point, or signed Frobenius translates
/// of one another, in the prime-order subgroup. Those columns are
/// algebraically dependent before any relation is collected. Merging them
/// uses only point multiplication, equality, negation and Frobenius; it does
/// not compute or attach a discrete logarithm.
fn projected_signed_orbit_map(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
) -> ProjectedSignedOrbitMap {
    let projected: Vec<_> = fb
        .points
        .iter()
        .map(|point| kc.mul(point, &kc.cofactor))
        .collect();
    let mut representatives: Vec<BinaryPoint> = Vec::new();
    let mut seen = HashSet::new();

    for point in &projected {
        if *point == BinaryPoint::Infinity || seen.contains(&point_key(point)) {
            continue;
        }
        let mut current = point.clone();
        let mut canonical = point.clone();
        let mut canonical_key = point_key(point);
        for _ in 0..kc.n {
            for candidate in [current.clone(), point_neg(&current)] {
                let key = point_key(&candidate);
                seen.insert(key.clone());
                if key < canonical_key {
                    canonical = candidate;
                    canonical_key = key;
                }
            }
            current = kc.frobenius(&current);
        }
        representatives.push(canonical);
    }
    representatives.sort_by_key(point_key);

    let mut location_by_key = HashMap::new();
    for (orbit, representative) in representatives.iter().enumerate() {
        let mut current = representative.clone();
        for k in 0..kc.n {
            location_by_key
                .entry(point_key(&current))
                .or_insert((orbit, k, false));
            location_by_key
                .entry(point_key(&point_neg(&current)))
                .or_insert((orbit, k, true));
            current = kc.frobenius(&current);
        }
    }
    let orbit_of = projected
        .iter()
        .map(|point| {
            if *point == BinaryPoint::Infinity {
                return None;
            }
            Some(
                *location_by_key
                    .get(&point_key(point))
                    .expect("cofactor projection was not found in its canonical orbit"),
            )
        })
        .collect();

    ProjectedSignedOrbitMap {
        orbit_of,
        representatives,
    }
}

/// Number of nonzero signed-Frobenius columns remaining after every
/// factor-base point is projected into the prime-order subgroup by `[h]`.
/// This is a public algebraic property of the factor base, not a log rank.
pub fn projected_signed_orbit_count(kc: &KoblitzCurve, fb: &FrobeniusFactorBase) -> usize {
    projected_signed_orbit_map(kc, fb).representatives.len()
}

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
    /// Macaulay degree whose implied rows are handed to the SAT solver
    /// alongside the system, or `None` to encode the system alone.
    /// Degree 2 is cheap and helps on overdetermined instances; degree
    /// 3 buys fewer conflicts but far more clauses.
    pub sat_macaulay_degree: Option<u32>,
    /// Native-XOR/CNF, branching, domain, trace, and conflict controls
    /// for [`DecompositionStrategy::Sat`].
    pub sat_options: SatDecompositionOptions,
    /// Identify `P` and `-P` as one signed Frobenius relation unknown.
    /// Disable only for a matched Frobenius-only control.
    pub collapse_negation: bool,
    /// Attempt the relation solve as soon as `unknowns + 1` equations
    /// are available and stop only when the recovered scalar verifies.
    /// Disable to retain the fixed-surplus collection control.
    pub stop_on_verified_rank: bool,
    /// Independent SAT targets launched per deterministic relation batch.
    /// Values above one use Rayon; non-SAT strategies remain serial.
    pub relation_batch_size: usize,
    /// Permit `aG+bQ=O` to recover the target scalar without the factor-base
    /// relation matrix. Disable in index-calculus benchmarks so this generic
    /// direct relation is counted and skipped rather than credited as a solve.
    pub allow_direct_relation: bool,
    /// Merge factor-base columns whose cofactor projections are equal up to
    /// signed Frobenius. This removes public algebraic dependencies without
    /// computing any factor-base logarithm.
    pub collapse_projected_orbits: bool,
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
            sat_macaulay_degree: Some(2),
            sat_options: SatDecompositionOptions::default(),
            collapse_negation: true,
            stop_on_verified_rank: true,
            relation_batch_size: 1,
            allow_direct_relation: true,
            collapse_projected_orbits: false,
        }
    }
}

/// What a run of the algorithm did.
#[derive(Clone, Debug)]
pub struct KoblitzIcReport {
    /// `|F|`.
    pub factor_base_size: usize,
    /// Number of signed `π`-orbits, i.e. the number of unknowns
    /// actually solved for.
    pub orbit_count: usize,
    /// `ℓ = ord_n(2)`, the dimension of the invariant subspace.
    pub ell: u32,
    /// Relations collected.
    pub relations: usize,
    /// Relations that raised the rank of the reduced relation matrix.
    pub independent_relations: usize,
    /// Relations that were linear combinations of earlier ones.
    pub dependent_relations: usize,
    /// Relations contradicting earlier ones.  Any nonzero value means a
    /// wrong relation was produced; the run fails closed.
    pub inconsistent_relations: usize,
    /// Pinned scalars that did not verify as `[d]G = Q`.  Any nonzero
    /// value likewise invalidates the run.
    pub verification_failures: usize,
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
    /// Inconclusive attempts (model/conflict caps or invalid models).
    pub sat_unknowns: usize,
    /// Encoding/model verification failures; any nonzero value invalidates a run.
    pub sat_invalid_models: usize,
    /// SAT models examined across relation collection.
    pub sat_models: usize,
    /// Cumulative CDCL conflicts across relation collection.
    pub sat_conflicts: u64,
    /// Time spent generating candidate targets and decomposing them.
    pub relation_collection_ns: u128,
    /// Time spent solving the final modular relation matrix.
    pub linear_algebra_ns: u128,
    /// Modular relation solves attempted, including rank-deficient ones.
    pub linear_solve_attempts: usize,
    /// Whether negation was folded into signed Frobenius columns.
    pub collapse_negation: bool,
    /// Requested SAT target batch size.
    pub relation_batch_size: usize,
    /// Relation batches actually launched.
    pub relation_batches: usize,
    /// Log recovered from R = O directly, bypassing the relation matrix.
    pub direct_relation: bool,
    /// Direct `aG+bQ=O` trials skipped because the benchmark forbade the
    /// relation-matrix bypass.
    pub direct_relations_skipped: usize,
    /// Whether the factor base's cofactor classes can sum to zero with
    /// exactly the requested number of summands.
    pub m_cofactor_admissible: bool,
    /// Whether public cofactor-projection dependencies were merged.
    pub collapse_projected_orbits: bool,
    /// Time spent constructing the public signed-Frobenius projection map.
    pub projected_orbit_construction_ns: u128,
    /// Time spent proving the requested summand count can reach the subgroup.
    pub cofactor_admission_ns: u128,
    /// Pair sums stored by the meet-in-the-middle oracle (zero for
    /// every other strategy).
    pub pair_table_entries: usize,
    /// Time spent building that table.
    pub pair_table_ns: u128,
}

fn solve_relation_system(
    relations: &[KoblitzRelation],
    relation_unknowns: usize,
    cofactor: &BigUint,
    modulus: &BigUint,
) -> Option<BigUint> {
    let h = cofactor % modulus;
    let mut matrix = Vec::with_capacity(relations.len());
    let mut rhs = Vec::with_capacity(relations.len());
    for relation in relations {
        let mut row = relation.row.clone();
        row.push((modulus - (&h * &relation.coef_b) % modulus) % modulus);
        matrix.push(row);
        rhs.push((&h * &relation.coef_a) % modulus);
    }
    let solution = gaussian_eliminate_mod_n(&mut matrix, &mut rhs, modulus)?;
    solution.get(relation_unknowns).cloned()
}

enum RelationAttemptOutcome {
    Direct,
    Enumerated(Option<Vec<usize>>),
    Groebner(Option<Vec<usize>>, SolveStats),
    Sat(Option<Vec<usize>>, SatDecompositionStats),
}

/// Live milestones from the existing small-curve pipeline.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum KoblitzIcEvent {
    FactorBaseStarted,
    FactorBaseReady {
        points: usize,
        orbits: usize,
    },
    /// The meet-in-the-middle pair table was built.
    PairTableReady {
        entries: usize,
    },
    RelationCollectionStarted {
        wanted: usize,
    },
    RelationProgress {
        collected: usize,
        wanted: usize,
        trials: usize,
    },
    RelationCollectionFinished {
        collected: usize,
        trials: usize,
    },
    LinearAlgebraStarted {
        rows: usize,
        columns: usize,
    },
    LinearAlgebraFinished,
    /// The current relation matrix did not yield a candidate.
    LinearAlgebraIncomplete,
    LinearAlgebraSkipped,
    VerificationStarted,
    VerificationFinished {
        verified: bool,
    },
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
    koblitz_index_calculus_dlp_with_progress(kc, q, opts, &mut |_| {})
}

/// The same small-curve solver with synchronous progress notifications.
/// Events follow real operations and can repeat when rank checks resume collection.
pub fn koblitz_index_calculus_dlp_with_progress(
    kc: &KoblitzCurve,
    q: &BinaryPoint,
    opts: &KoblitzIcOptions,
    progress: &mut dyn FnMut(KoblitzIcEvent),
) -> Option<KoblitzIcReport> {
    progress(KoblitzIcEvent::FactorBaseStarted);
    let fb = build_frobenius_factor_base(kc, opts.factor_index)?;
    koblitz_index_calculus_dlp_observed(kc, q, &fb, opts, progress)
}

/// Run index calculus with a caller-supplied invariant factor base.
///
/// This admits explicit nonlinear orbit unions selected by a finite
/// search while retaining the same enumeration, Gröbner, and SAT
/// decomposition controls as [`koblitz_index_calculus_dlp`].
pub fn koblitz_index_calculus_dlp_with_factor_base(
    kc: &KoblitzCurve,
    q: &BinaryPoint,
    fb: &FrobeniusFactorBase,
    opts: &KoblitzIcOptions,
) -> Option<KoblitzIcReport> {
    koblitz_index_calculus_dlp_observed(kc, q, fb, opts, &mut |_| {})
}

/// Run index calculus with a caller-supplied factor base, reporting
/// progress.  Emits [`KoblitzIcEvent::FactorBaseReady`] but not
/// `FactorBaseStarted`, since the base was built by the caller.
pub fn koblitz_index_calculus_dlp_with_factor_base_and_progress(
    kc: &KoblitzCurve,
    q: &BinaryPoint,
    fb: &FrobeniusFactorBase,
    opts: &KoblitzIcOptions,
    progress: &mut dyn FnMut(KoblitzIcEvent),
) -> Option<KoblitzIcReport> {
    koblitz_index_calculus_dlp_observed(kc, q, fb, opts, progress)
}

fn koblitz_index_calculus_dlp_observed(
    kc: &KoblitzCurve,
    q: &BinaryPoint,
    fb: &FrobeniusFactorBase,
    opts: &KoblitzIcOptions,
    progress: &mut dyn FnMut(KoblitzIcEvent),
) -> Option<KoblitzIcReport> {
    let r = &kc.subgroup_order;
    let g = kc.generator().clone();
    let projection_start = std::time::Instant::now();
    let projected_orbits = opts
        .collapse_projected_orbits
        .then(|| projected_signed_orbit_map(kc, fb));
    let projected_orbit_construction_ns = projection_start.elapsed().as_nanos();
    let relation_unknowns = projected_orbits.as_ref().map_or_else(
        || {
            if opts.collapse_negation {
                fb.unknowns()
            } else {
                fb.orbits.len()
            }
        },
        |projected| projected.representatives.len(),
    );
    let admission_start = std::time::Instant::now();
    let m_cofactor_admissible = fb.m_can_decompose(kc, opts.m);
    let cofactor_admission_ns = admission_start.elapsed().as_nanos();

    let index_of = fb.index_map();
    progress(KoblitzIcEvent::FactorBaseReady {
        points: fb.points.len(),
        orbits: relation_unknowns,
    });

    let mut report = KoblitzIcReport {
        factor_base_size: fb.points.len(),
        orbit_count: relation_unknowns,
        ell: fb.ell,
        relations: 0,
        independent_relations: 0,
        dependent_relations: 0,
        inconsistent_relations: 0,
        verification_failures: 0,
        trials: 0,
        log: None,
        reductions: 0,
        infeasible_branches: 0,
        sat_calls: 0,
        sat_refutations: 0,
        sat_unknowns: 0,
        sat_invalid_models: 0,
        sat_models: 0,
        sat_conflicts: 0,
        relation_collection_ns: 0,
        linear_algebra_ns: 0,
        linear_solve_attempts: 0,
        collapse_negation: opts.collapse_negation,
        relation_batch_size: opts.relation_batch_size.max(1),
        relation_batches: 0,
        direct_relation: false,
        direct_relations_skipped: 0,
        m_cofactor_admissible,
        collapse_projected_orbits: opts.collapse_projected_orbits,
        projected_orbit_construction_ns,
        cofactor_admission_ns,
        pair_table_entries: 0,
        pair_table_ns: 0,
    };
    if fb.points.is_empty() || !m_cofactor_admissible {
        return Some(report);
    }

    // The meet-in-the-middle oracle needs its table once per run.
    let pair_table = if opts.strategy == DecompositionStrategy::PairTable {
        let start = std::time::Instant::now();
        let table = PairSumTable::build(kc, fb)?;
        report.pair_table_ns = start.elapsed().as_nanos();
        report.pair_table_entries = table.len();
        progress(KoblitzIcEvent::PairTableReady {
            entries: table.len(),
        });
        Some(table)
    } else {
        None
    };

    let field = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let wanted = relation_unknowns + opts.extra_relations.max(1);
    let mut relations: Vec<KoblitzRelation> = Vec::with_capacity(wanted);
    // Incremental reduced echelon form over Z/rZ; the dense big-integer
    // solver is the fallback for a modulus wider than 64 bits.
    let mut echelon = IncrementalRelationSolver::new(relation_unknowns, r);
    let mut rng = StdRng::seed_from_u64(opts.seed);
    let r_u64 = r.to_u64_digits().first().copied().unwrap_or(1).max(2);
    let relation_start = std::time::Instant::now();

    // Extract, announce and verify a candidate scalar from the echelon
    // form.  `Some(true)` means solved; `Some(false)` means a pinned
    // scalar failed verification, which only a wrong relation can cause.
    let finish_incremental = |echelon: &IncrementalRelationSolver,
                                  report: &mut KoblitzIcReport,
                                  progress: &mut dyn FnMut(KoblitzIcEvent)|
     -> Option<bool> {
        let d = echelon.target_biguint()?;
        progress(KoblitzIcEvent::LinearAlgebraStarted {
            rows: echelon.rank(),
            columns: relation_unknowns + 1,
        });
        report.linear_solve_attempts += 1;
        progress(KoblitzIcEvent::LinearAlgebraFinished);
        progress(KoblitzIcEvent::VerificationStarted);
        let verified = kc.mul(&g, &d) == *q;
        progress(KoblitzIcEvent::VerificationFinished { verified });
        if verified {
            report.log = Some(d);
        } else {
            report.verification_failures += 1;
        }
        Some(verified)
    };

    progress(KoblitzIcEvent::RelationCollectionStarted { wanted });
    while relations.len() < wanted && report.trials < opts.max_trials {
        let remaining_trials = opts.max_trials - report.trials;
        let batch_size = opts.relation_batch_size.max(1).min(remaining_trials);
        let attempts: Vec<_> = (0..batch_size)
            .map(|_| {
                let a = BigUint::from(rng.gen_range(1..r_u64));
                let b = BigUint::from(rng.gen_range(1..r_u64));
                let target = kc.add(&kc.mul(&g, &a), &kc.mul(q, &b));
                (a, b, target)
            })
            .collect();
        report.trials += attempts.len();
        report.relation_batches += 1;

        let evaluate = |(_, _, target): &(BigUint, BigUint, BinaryPoint)| {
            if *target == BinaryPoint::Infinity {
                return RelationAttemptOutcome::Direct;
            }
            match opts.strategy {
                DecompositionStrategy::Enumerate => RelationAttemptOutcome::Enumerated(decompose(
                    kc, fb, &index_of, target, opts.m, 0,
                )),
                DecompositionStrategy::PairTable => RelationAttemptOutcome::Enumerated(
                    pair_table
                        .as_ref()
                        .and_then(|table| table.decompose(kc, fb, target, opts.m)),
                ),
                DecompositionStrategy::Groebner => {
                    let (idxs, stats) = groebner_decompose(
                        kc,
                        fb,
                        &index_of,
                        &field,
                        target,
                        opts.m,
                        opts.engine,
                        opts.node_budget,
                    );
                    RelationAttemptOutcome::Groebner(idxs, stats)
                }
                DecompositionStrategy::Sat => {
                    let (idxs, stats) = sat_decompose_with(
                        kc,
                        fb,
                        &index_of,
                        &field,
                        target,
                        opts.m,
                        opts.max_models,
                        opts.sat_macaulay_degree,
                        opts.sat_options,
                    );
                    RelationAttemptOutcome::Sat(idxs, stats)
                }
            }
        };
        // Targets are drawn serially above and consumed in order below,
        // so the outcome is independent of thread scheduling.
        let outcomes: Vec<_> = if batch_size > 1 {
            attempts.par_iter().map(evaluate).collect()
        } else {
            attempts.iter().map(evaluate).collect()
        };

        for outcome in &outcomes {
            match outcome {
                RelationAttemptOutcome::Groebner(_, stats) => {
                    report.reductions += stats.reductions;
                    report.infeasible_branches += stats.infeasible_branches;
                }
                RelationAttemptOutcome::Sat(_, stats) => {
                    report.sat_calls += stats.solver_calls;
                    report.sat_refutations += usize::from(stats.refuted);
                    report.sat_unknowns += usize::from(stats.exhausted);
                    report.sat_invalid_models += stats.spurious;
                    report.sat_models += stats.models;
                    report.sat_conflicts += stats.conflicts;
                }
                RelationAttemptOutcome::Direct | RelationAttemptOutcome::Enumerated(_) => {}
            }
        }
        if report.sat_invalid_models != 0 {
            report.relations = relations.len();
            report.relation_collection_ns = relation_start
                .elapsed()
                .as_nanos()
                .saturating_sub(report.linear_algebra_ns);
            return Some(report);
        }

        let mut inconsistent = false;
        for ((a, b, _target), outcome) in attempts.into_iter().zip(outcomes) {
            let found = match outcome {
                RelationAttemptOutcome::Direct => {
                    if !opts.allow_direct_relation {
                        report.direct_relations_skipped += 1;
                        continue;
                    }
                    let d = solve_for_d(&a, &b, r)?;
                    progress(KoblitzIcEvent::RelationCollectionFinished {
                        collected: relations.len(),
                        trials: report.trials,
                    });
                    progress(KoblitzIcEvent::LinearAlgebraSkipped);
                    progress(KoblitzIcEvent::VerificationStarted);
                    if kc.mul(&g, &d) == *q {
                        report.log = Some(d);
                        report.direct_relation = true;
                        report.relations = relations.len();
                        report.independent_relations =
                            echelon.as_ref().map_or(0, IncrementalRelationSolver::rank);
                        report.relation_collection_ns = relation_start
                            .elapsed()
                            .as_nanos()
                            .saturating_sub(report.linear_algebra_ns);
                        progress(KoblitzIcEvent::VerificationFinished { verified: true });
                        return Some(report);
                    }
                    progress(KoblitzIcEvent::VerificationFinished { verified: false });
                    progress(KoblitzIcEvent::RelationCollectionStarted { wanted });
                    None
                }
                RelationAttemptOutcome::Enumerated(idxs)
                | RelationAttemptOutcome::Groebner(idxs, _)
                | RelationAttemptOutcome::Sat(idxs, _) => idxs,
            };
            if let Some(idxs) = found {
                let relation = relation_from_decomposition_with_mode(
                    kc,
                    fb,
                    &idxs,
                    &a,
                    &b,
                    opts.collapse_negation,
                    projected_orbits.as_ref(),
                );
                if let Some(echelon) = echelon.as_mut() {
                    let linear_start = std::time::Instant::now();
                    let status = echelon.add_relation(&relation, &kc.cofactor);
                    report.linear_algebra_ns += linear_start.elapsed().as_nanos();
                    match status {
                        RowStatus::Independent => report.independent_relations += 1,
                        RowStatus::Dependent => report.dependent_relations += 1,
                        RowStatus::Inconsistent => {
                            report.inconsistent_relations += 1;
                            inconsistent = true;
                        }
                    }
                }
                relations.push(relation);
            }
        }

        progress(KoblitzIcEvent::RelationProgress {
            collected: relations.len(),
            wanted,
            trials: report.trials,
        });
        if inconsistent {
            // A relation contradicts the others: only a wrong
            // decomposition-to-row rewrite can do that.  Fail closed.
            break;
        }
        if opts.stop_on_verified_rank {
            if let Some(echelon) = echelon.as_ref() {
                let linear_start = std::time::Instant::now();
                let outcome = finish_incremental(echelon, &mut report, progress);
                report.linear_algebra_ns += linear_start.elapsed().as_nanos();
                match outcome {
                    Some(true) => {
                        report.relations = relations.len();
                        report.relation_collection_ns = relation_start
                            .elapsed()
                            .as_nanos()
                            .saturating_sub(report.linear_algebra_ns);
                        return Some(report);
                    }
                    // Pinned but wrong: a wrong relation slipped in.
                    Some(false) => break,
                    None => {}
                }
            } else if relations.len() >= relation_unknowns + 1 {
                progress(KoblitzIcEvent::LinearAlgebraStarted {
                    rows: relations.len(),
                    columns: relation_unknowns + 1,
                });
                report.linear_solve_attempts += 1;
                let linear_start = std::time::Instant::now();
                let candidate =
                    solve_relation_system(&relations, relation_unknowns, &kc.cofactor, r);
                report.linear_algebra_ns += linear_start.elapsed().as_nanos();
                if candidate.is_none() {
                    progress(KoblitzIcEvent::LinearAlgebraIncomplete);
                }
                let candidate = candidate.filter(|d| {
                    progress(KoblitzIcEvent::LinearAlgebraFinished);
                    progress(KoblitzIcEvent::VerificationStarted);
                    let verified = kc.mul(&g, d) == *q;
                    progress(KoblitzIcEvent::VerificationFinished { verified });
                    verified
                });
                if let Some(d) = candidate {
                    report.log = Some(d);
                    report.relations = relations.len();
                    report.relation_collection_ns = relation_start
                        .elapsed()
                        .as_nanos()
                        .saturating_sub(report.linear_algebra_ns);
                    return Some(report);
                }
                progress(KoblitzIcEvent::RelationCollectionStarted { wanted });
            }
        }
    }
    report.relation_collection_ns = relation_start
        .elapsed()
        .as_nanos()
        .saturating_sub(report.linear_algebra_ns);
    report.relations = relations.len();
    progress(KoblitzIcEvent::RelationCollectionFinished {
        collected: relations.len(),
        trials: report.trials,
    });
    if report.inconsistent_relations > 0 || report.verification_failures > 0 {
        return Some(report);
    }

    // Unknowns: x_1 … x_s (orbit logs) and d, in the last column.
    //   Σ_o c_o x_o  −  (h·b)·d  ≡  h·a   (mod r)
    if let Some(echelon) = echelon.as_ref() {
        if relations.is_empty() {
            return Some(report);
        }
        let linear_start = std::time::Instant::now();
        let outcome = finish_incremental(echelon, &mut report, progress);
        report.linear_algebra_ns += linear_start.elapsed().as_nanos();
        if outcome.is_none() {
            progress(KoblitzIcEvent::LinearAlgebraStarted {
                rows: echelon.rank(),
                columns: relation_unknowns + 1,
            });
            progress(KoblitzIcEvent::LinearAlgebraIncomplete);
        }
        return Some(report);
    }
    if relations.len() < relation_unknowns + 1 {
        return Some(report);
    }
    progress(KoblitzIcEvent::LinearAlgebraStarted {
        rows: relations.len(),
        columns: relation_unknowns + 1,
    });
    report.linear_solve_attempts += 1;
    let linear_start = std::time::Instant::now();
    let Some(d) = solve_relation_system(&relations, relation_unknowns, &kc.cofactor, r) else {
        report.linear_algebra_ns += linear_start.elapsed().as_nanos();
        progress(KoblitzIcEvent::LinearAlgebraIncomplete);
        return Some(report);
    };
    report.linear_algebra_ns += linear_start.elapsed().as_nanos();
    progress(KoblitzIcEvent::LinearAlgebraFinished);
    progress(KoblitzIcEvent::VerificationStarted);
    let verified = kc.mul(&g, &d) == *q;
    if verified {
        report.log = Some(d);
    }
    progress(KoblitzIcEvent::VerificationFinished { verified });
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
    /// Unknowns after collapsing signed `π`-orbits.
    pub unknowns: usize,
    /// `|F| / unknowns` — the measured relation-collection speed-up,
    /// which is generically `≈ 2n` after Frobenius and negation.
    pub relation_collection_speedup: f64,
    /// `(|F| / unknowns)²` — the measured linear-algebra speed-up for a
    /// quadratic-cost sparse solve, generically `≈ 4n²` after both
    /// identifications. The Frobenius contribution alone is `n²`.
    pub linear_algebra_speedup: f64,
    /// `m!`, the symmetry-breaking saving in the decomposition search.
    pub symmetry_breaking_factor: f64,
    /// The Gallant–Lambert–Vanstone / Wiener–Zuccherato speed-up that
    /// Pollard rho already gets from the same endomorphism, `√(2n)`.
    pub rho_speedup: f64,
}

/// Cost model for a Frobenius-invariant index calculus on `n`, a factor
/// base of `fb_size` points collapsing to `unknowns` signed orbits, with
/// `m`-point decompositions.
///
/// The point of the table is the comparison in the paper's conclusion:
/// index calculus gains more from the Koblitz structure than rho does
/// (`n` and `n²` from Frobenius, plus the standard negation quotient,
/// versus `√(2n)`) — and yet stays the worse attack for
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
    fn cyclotomic_cosets_partition_and_classify() {
        for n in [7u32, 9, 15, 21, 31, 63] {
            let cosets = cyclotomic_cosets(n);
            let mut all: Vec<u32> = cosets.iter().flatten().copied().collect();
            all.sort_unstable();
            assert_eq!(
                all,
                (0..n).collect::<Vec<_>>(),
                "cosets must partition Z/{n}"
            );
            // Each coset is closed under doubling.
            for c in &cosets {
                for s in c {
                    assert!(c.contains(&((s * 2) % n)));
                }
            }
        }
        // 2 is primitive mod 131 and mod 163, so there are only two
        // cosets and no usable invariant subspace exists at all — the
        // sizes that matter are exactly the ones this construction
        // cannot serve.
        for n in [131u32, 163] {
            assert_eq!(cyclotomic_cosets(n).len(), 2);
            assert_eq!(available_subspace_dimensions(n), vec![0, 1, n - 1, n]);
        }
    }

    #[test]
    fn the_factorisation_is_complete() {
        // One irreducible factor per cyclotomic coset, of that coset's
        // degree — so the degrees sum to n, and the dimensions the
        // divisors reach are exactly the classification's.
        for n in [7u32, 9, 15, 21, 31, 63] {
            let factors = all_factors_of_x_n_minus_1(n);
            assert_eq!(factors.len(), cyclotomic_cosets(n).len(), "n = {n}");
            let degrees: Vec<u32> = factors.iter().map(|f| 63 - f.leading_zeros()).collect();
            assert_eq!(
                degrees.iter().sum::<u32>(),
                n,
                "degrees must sum to n = {n}"
            );

            let mut reach = std::collections::BTreeSet::from([0u32]);
            for d in &degrees {
                for r in reach.clone() {
                    reach.insert(r + d);
                }
            }
            assert_eq!(
                reach.into_iter().collect::<Vec<_>>(),
                available_subspace_dimensions(n),
                "reachable dimensions must match the coset classification at n = {n}"
            );
        }
    }

    #[test]
    fn divisor_factor_bases_are_invariant_and_sized_to_the_divisor() {
        let kc = KoblitzCurve::new(1, 15).unwrap();
        for indices in [vec![1usize], vec![1, 2], vec![1, 2, 3]] {
            let fb = build_frobenius_factor_base_from_divisor(&kc, &indices).unwrap();
            let expected: u32 = indices
                .iter()
                .map(|&i| 63 - all_factors_of_x_n_minus_1(15)[i].leading_zeros())
                .sum();
            assert_eq!(fb.ell, expected, "dimension is the divisor's degree");
            assert_eq!(fb.subspace.len(), 1usize << expected);

            let keys: std::collections::HashSet<_> = fb.points.iter().map(point_key).collect();
            for p in &fb.points {
                assert!(keys.contains(&point_key(&kc.frobenius(p))), "π(F) ⊆ F");
                assert!(keys.contains(&point_key(&point_neg(p))), "−F ⊆ F");
            }
            assert_eq!(
                fb.orbits.iter().map(|o| o.len()).sum::<usize>(),
                fb.points.len()
            );
        }
    }

    #[test]
    fn a_bigger_invariant_subspace_grows_the_base_the_divisor_asks_for() {
        // The point of divisor bases: |F| is tunable, so the summand
        // count m ≈ n/dim a decomposition needs can be brought down to
        // 2 — which also replaces the cubic chained system with a
        // quadratic one.
        let kc = KoblitzCurve::new(1, 7).unwrap();
        let small = build_frobenius_factor_base_from_divisor(&kc, &[1]).unwrap();
        let large = build_frobenius_factor_base_from_divisor(&kc, &[1, 2]).unwrap();
        assert_eq!((small.ell, large.ell), (3, 6));
        assert_eq!(small.points.len(), 15);
        assert_eq!(large.points.len(), 71);
        assert!(large.orbits.len() > small.orbits.len());
    }

    #[test]
    fn the_cofactor_class_decides_which_m_can_decompose() {
        // Not a size effect: on K_1/F_2^7 (cofactor 2) no factor-base
        // point lies in ⟨G⟩, so a sum of m of them reaches ⟨G⟩ only for
        // even m.  Odd m decomposes nothing however big the base is.
        let kc = KoblitzCurve::new(1, 7).unwrap();
        let g = kc.generator().clone();
        for indices in [vec![1usize], vec![1, 2]] {
            let fb = build_frobenius_factor_base_from_divisor(&kc, &indices).unwrap();
            let index_of = fb.index_map();
            assert!(
                fb.cofactor_classes(&kc)
                    .iter()
                    .all(|c| *c != BinaryPoint::Infinity),
                "no point of this base lies in ⟨G⟩"
            );
            assert_eq!(fb.admissible_summand_counts(&kc, 4), vec![2, 4]);

            // …and the predicate matches what search actually finds.
            for m in 2..=4usize {
                let hits = (1..12u32)
                    .filter(|k| {
                        let t = kc.mul(&g, &BigUint::from(*k));
                        enumerate_decompose(&kc, &fb, &index_of, &t, m).is_some()
                    })
                    .count();
                if fb.m_can_decompose(&kc, m) {
                    assert_eq!(hits, 11, "m = {m} is admissible and should decompose");
                } else {
                    assert_eq!(
                        hits, 0,
                        "m = {m} is inadmissible and must decompose nothing"
                    );
                }
            }
        }
    }

    #[test]
    fn orbit_derived_cofactor_classes_match_pointwise_projection() {
        for (a, n, indices) in [(1, 15, vec![2usize, 4]), (0, 23, vec![0usize, 2])] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let fb = build_frobenius_factor_base_from_divisor(&kc, &indices).unwrap();
            let pointwise: HashSet<_> = fb.cofactor_classes(&kc).iter().map(point_key).collect();
            let orbit_derived: HashSet<_> = fb
                .distinct_cofactor_classes(&kc)
                .iter()
                .map(point_key)
                .collect();
            assert_eq!(orbit_derived, pointwise, "a={a}, n={n}");
        }
    }

    #[test]
    fn projected_orbit_map_is_algebraic_and_exact() {
        let kc = KoblitzCurve::new(1, 15).unwrap();
        let fb = build_frobenius_factor_base_from_divisor(&kc, &[0, 2]).unwrap();
        let projected = projected_signed_orbit_map(&kc, &fb);
        assert_eq!(
            projected.representatives.len(),
            projected_signed_orbit_count(&kc, &fb)
        );
        assert!(projected.representatives.len() < fb.unknowns());
        for (index, point) in fb.points.iter().enumerate() {
            let expected = kc.mul(point, &kc.cofactor);
            let Some((orbit, k, negated)) = projected.orbit_of[index] else {
                assert_eq!(expected, BinaryPoint::Infinity);
                continue;
            };
            let mut reconstructed = projected.representatives[orbit].clone();
            for _ in 0..k {
                reconstructed = kc.frobenius(&reconstructed);
            }
            if negated {
                reconstructed = point_neg(&reconstructed);
            }
            assert_eq!(reconstructed, expected);
        }
    }

    #[test]
    fn inadmissible_summand_count_stops_before_relation_search() {
        let kc = KoblitzCurve::new(1, 7).unwrap();
        let fb = build_frobenius_factor_base_from_divisor(&kc, &[1]).unwrap();
        let report = koblitz_index_calculus_dlp_with_factor_base(
            &kc,
            kc.generator(),
            &fb,
            &KoblitzIcOptions {
                m: 3,
                strategy: DecompositionStrategy::Sat,
                allow_direct_relation: false,
                collapse_projected_orbits: true,
                ..KoblitzIcOptions::default()
            },
        )
        .unwrap();
        assert!(!report.m_cofactor_admissible);
        assert_eq!(
            (report.trials, report.sat_calls, report.relations),
            (0, 0, 0)
        );
        assert_eq!(report.log, None);
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
        // Signed-orbit bookkeeping additionally identifies -P with P.
        for (i, p) in fb.points.iter().enumerate() {
            let (o, k, negated) = fb.signed_orbit_of[i];
            let mut rep = fb.points[fb.signed_orbits[o][0]].clone();
            for _ in 0..k {
                rep = kc.frobenius(&rep);
            }
            if negated {
                rep = point_neg(&rep);
            }
            assert_eq!(&rep, p);
        }
        // The orbits partition the factor base.
        assert_eq!(
            fb.orbits.iter().map(|o| o.len()).sum::<usize>(),
            fb.points.len()
        );
        assert_eq!(
            fb.signed_orbits.iter().map(|o| o.len()).sum::<usize>(),
            fb.points.len()
        );
        assert!(fb.signed_orbits.len() <= fb.orbits.len());
        let negated_index = fb
            .signed_orbit_of
            .iter()
            .position(|location| location.2)
            .expect("the negation-closed base has a negative representative");
        let (signed_orbit, k, _) = fb.signed_orbit_of[negated_index];
        let relation = relation_from_decomposition(
            &kc,
            &fb,
            &[negated_index],
            &BigUint::zero(),
            &BigUint::zero(),
        );
        let lambda_k = kc.lambda.modpow(&BigUint::from(k), &kc.subgroup_order);
        assert_eq!(relation.row[signed_orbit], &kc.subgroup_order - lambda_k);
        assert_eq!(relation.summand_negated, vec![true]);
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
    fn selected_explicit_orbit_base_has_four_relation_columns() {
        let kc = KoblitzCurve::new(1, 19).unwrap();
        let representatives: Vec<_> = [16795u64, 1315, 8461, 6685]
            .into_iter()
            .map(|value| F2mElement::from_biguint(&BigUint::from(value), kc.n))
            .collect();
        let fb = build_explicit_frobenius_orbit_factor_base(&kc, &representatives).unwrap();
        assert_eq!(
            fb.domain,
            FactorBaseDomain::ExplicitFrobeniusOrbits { representatives: 4 }
        );
        assert_eq!(fb.subspace.len(), 76);
        assert_eq!(fb.points.len(), 152);
        assert_eq!(fb.orbits.len(), 8);
        assert_eq!(fb.signed_orbits.len(), 4);
        assert_eq!(fb.unknowns(), 4);
    }

    #[test]
    #[ignore = "seconds; deliberate n=19 selected-base SAT round trip"]
    fn selected_explicit_orbit_base_sat_round_trip_n19() {
        let kc = KoblitzCurve::new(1, 19).unwrap();
        let representatives: Vec<_> = [16795u64, 1315, 8461, 6685]
            .into_iter()
            .map(|value| F2mElement::from_biguint(&BigUint::from(value), kc.n))
            .collect();
        let fb = build_explicit_frobenius_orbit_factor_base(&kc, &representatives).unwrap();
        let index = fb.index_map();
        let target = [203u64, 2143, 2901]
            .into_iter()
            .map(|scalar| kc.mul(kc.generator(), &BigUint::from(scalar)))
            .fold(BinaryPoint::Infinity, |sum, point| kc.add(&sum, &point));
        let field = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let options = SatDecompositionOptions {
            branch_on_summands: true,
            restrict_to_factor_base: true,
            conflict_budget: 2_000_000,
            ..Default::default()
        };
        let (decomposition, stats) =
            sat_decompose_with(&kc, &fb, &index, &field, &target, 3, 1024, None, options);
        eprintln!("selected explicit-orbit SAT stats: {stats:?}");
        assert_eq!(stats.spurious, 0);
        assert!(!stats.exhausted, "selected-base SAT exhausted: {stats:?}");
        let decomposition = decomposition.expect("planted target must decompose");
        let sum = decomposition
            .iter()
            .fold(BinaryPoint::Infinity, |sum, &i| kc.add(&sum, &fb.points[i]));
        assert_eq!(sum, target);
    }

    #[test]
    fn relation_rows_are_consistent_with_the_orbit_logs() {
        // A relation must satisfy h·a + h·b·d ≡ Σ +-λ^k x_o with
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
            let rep = fb.points[fb.signed_orbits[o][0]].clone();
            let scaled = kc.mul(&rep, &kc.cofactor);
            let x_o = logs.get(&point_key(&scaled)).expect("in ⟨G⟩ after ×h");
            rhs = (rhs + coeff * x_o) % r;
        }
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn reciprocal_saturation_preserves_projection_and_contains_the_base() {
        let kc = KoblitzCurve::new(1, 15).unwrap();
        let fb = build_frobenius_factor_base_from_divisor(&kc, &[0, 3]).unwrap();
        let saturated = saturate_factor_base_two_torsion(&kc, &fb).unwrap();
        let keys: std::collections::HashSet<_> = saturated.points.iter().map(point_key).collect();
        for p in &fb.points {
            assert!(keys.contains(&point_key(p)));
        }
        let project = |base: &FrobeniusFactorBase| -> std::collections::HashSet<_> {
            base.points
                .iter()
                .map(|p| point_key(&kc.mul(p, &kc.cofactor)))
                .collect()
        };
        assert_eq!(project(&fb), project(&saturated));
        let t = points_with_x(&kc.curve, &F2mElement::zero(kc.n)).remove(0);
        for p in &saturated.points {
            assert!(keys.contains(&point_key(&kc.frobenius(p))));
            if let BinaryPoint::Affine { x, .. } = p {
                if !x.is_zero() {
                    let BinaryPoint::Affine {
                        x: translated_x, ..
                    } = kc.add(p, &t)
                    else {
                        panic!("non-torsion point");
                    };
                    assert_eq!(
                        x.mul(&translated_x, &kc.curve.irreducible),
                        F2mElement::one(kc.n)
                    );
                }
            }
        }
    }

    #[test]
    fn domain_trie_accepts_exactly_the_allowed_coordinates() {
        use crate::cryptanalysis::sat::{SolveResult, Solver};
        for ell in 1..=5usize {
            for mode in 0..4 {
                let codes: Vec<_> = (0..(1u64 << ell))
                    .filter(|x| match mode {
                        0 => false,
                        1 => true,
                        2 => x % 3 == 0,
                        _ => x & 1 == 1,
                    })
                    .collect();
                for x in 0..(1u64 << ell) {
                    let mut solver = Solver::new(ell as u32 + 2);
                    add_coordinate_domain(&mut solver, 2, ell, &codes);
                    for j in 0..ell {
                        let lit = (j + 3) as i32;
                        solver.add_clause(vec![if (x >> j) & 1 == 1 { lit } else { -lit }]);
                    }
                    assert_eq!(solver.solve() == SolveResult::Sat, codes.contains(&x));
                }
            }
        }
    }

    #[test]
    fn frobenius_union_domain_is_invariant_and_sat_matches_search() {
        let kc = KoblitzCurve::new(1, 7).unwrap();
        let basis = vec![F2mElement::one(7), F2mElement::from_bit_positions(&[1], 7)];
        let fb = build_frobenius_union_factor_base(&kc, &basis).unwrap();
        let keys: std::collections::HashSet<_> = fb.points.iter().map(point_key).collect();
        for p in &fb.points {
            assert!(keys.contains(&point_key(&kc.frobenius(p))));
            assert!(keys.contains(&point_key(&point_neg(p))));
        }
        assert!(matches!(
            fb.domain,
            FactorBaseDomain::FrobeniusUnion { seed_dimension: 2 }
        ));
        let index = fb.index_map();
        let st = FieldStructure::new(7, &kc.curve.irreducible);
        for m in [2, 3] {
            for k in 1u32..=6 {
                let target = kc.mul(kc.generator(), &BigUint::from(k));
                let reference = enumerate_decompose(&kc, &fb, &index, &target, m);
                let (out, stats) = sat_decompose_with(
                    &kc,
                    &fb,
                    &index,
                    &st,
                    &target,
                    m,
                    128,
                    None,
                    SatDecompositionOptions {
                        conflict_budget: 100_000,
                        ..Default::default()
                    },
                );
                assert_eq!(stats.spurious, 0);
                assert!(!stats.exhausted);
                assert_eq!(out.is_some(), reference.is_some());
            }
        }
    }

    #[test]
    fn trace_bit_is_a_group_homomorphism_including_torsion() {
        for a in [0, 1] {
            let kc = KoblitzCurve::new(a, 7).unwrap();
            let delta = |p: &BinaryPoint| match p {
                BinaryPoint::Infinity => false,
                BinaryPoint::Affine { x, .. } => {
                    absolute_trace_bit(x, kc.n, &kc.curve.irreducible) ^ (a == 1)
                }
            };
            let torsion = points_with_x(&kc.curve, &F2mElement::zero(kc.n)).remove(0);
            for x in 0..(1u32 << kc.n) {
                let x = F2mElement::from_biguint(&BigUint::from(x), kc.n);
                for p in points_with_x(&kc.curve, &x) {
                    for q in [&p, &point_neg(&p), kc.generator(), &torsion] {
                        assert_eq!(delta(&kc.add(&p, q)), delta(&p) ^ delta(q));
                    }
                }
            }
        }
    }

    #[test]
    fn native_sat_domain_filter_preserves_group_decompositions() {
        for (a, n, m) in [(1, 7, 2), (0, 9, 3)] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let fb = build_frobenius_factor_base(&kc, 0).unwrap();
            let index = fb.index_map();
            let st = FieldStructure::new(n, &kc.curve.irreducible);
            for k in 1u32..=6 {
                let target = kc.mul(kc.generator(), &BigUint::from(k));
                let reference = enumerate_decompose(&kc, &fb, &index, &target, m);
                for restrict in [false, true] {
                    let (out, stats) = sat_decompose_with(
                        &kc,
                        &fb,
                        &index,
                        &st,
                        &target,
                        m,
                        128,
                        None,
                        SatDecompositionOptions {
                            restrict_to_factor_base: restrict,
                            conflict_budget: 1_000_000,
                            ..Default::default()
                        },
                    );
                    assert!(!stats.exhausted);
                    assert_eq!(stats.spurious, 0);
                    assert_eq!(out.is_some(), reference.is_some());
                    assert_eq!(stats.refuted, out.is_none());
                    if let Some(ids) = out {
                        let sum = ids
                            .iter()
                            .fold(BinaryPoint::Infinity, |p, &i| kc.add(&p, &fb.points[i]));
                        assert_eq!(sum, target);
                    }
                }
            }
        }
    }

    #[test]
    fn sparse_and_exhaustive_irreducible_searches_agree_up_to_the_old_cap() {
        // KoblitzCurve::new now uses the sparse search; every field
        // representation a fixture could have been generated with before
        // is unchanged.
        for n in 3..=24u32 {
            let full = find_irreducible(n).expect("exhaustive search finds one");
            let sparse = find_irreducible_sparse(n).expect("sparse search finds one");
            assert_eq!((full.degree, full.low_terms), (sparse.degree, sparse.low_terms), "n = {n}");
        }
    }

    #[test]
    fn pair_table_agrees_with_enumeration_for_two_three_and_four_summands() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let index = fb.index_map();
        let table = PairSumTable::build(&kc, &fb).unwrap();
        assert_eq!(table.len(), fb.points.len() * (fb.points.len() + 1) / 2);
        let r = kc.subgroup_order.to_u64_digits()[0];
        for m in [2usize, 3, 4] {
            // m = 4 costs |F|³ per enumerated target; keep that one short.
            let scalars: Vec<u64> = if m == 4 {
                (1..=12).collect()
            } else {
                (1..r).collect()
            };
            for k in scalars {
                let target = kc.mul(kc.generator(), &BigUint::from(k));
                let reference = enumerate_decompose(&kc, &fb, &index, &target, m);
                let tabled = table.decompose(&kc, &fb, &target, m);
                assert_eq!(tabled.is_some(), reference.is_some(), "m = {m}, k = {k}");
                if let Some(idxs) = tabled {
                    assert_eq!(idxs.len(), m);
                    assert!(idxs.windows(2).all(|w| w[0] <= w[1]));
                    let sum = idxs
                        .iter()
                        .fold(BinaryPoint::Infinity, |s, &i| kc.add(&s, &fb.points[i]));
                    assert_eq!(sum, target);
                }
            }
        }
        // Every witness the table enumerates is genuine and distinct.
        let target = kc.mul(kc.generator(), &BigUint::from(53u32));
        let mut seen = HashSet::new();
        table.witnesses(&kc, &fb, &target, 3, &mut |idxs| {
            let sum = idxs
                .iter()
                .fold(BinaryPoint::Infinity, |s, &i| kc.add(&s, &fb.points[i]));
            assert_eq!(sum, target);
            assert!(seen.insert(idxs.to_vec()), "duplicate witness {idxs:?}");
            true
        });
    }

    #[test]
    fn packed_point_keys_are_injective() {
        let kc = KoblitzCurve::new(1, 11).unwrap();
        let mut keys = HashSet::new();
        let mut count = 1usize;
        keys.insert(pack_point(&BinaryPoint::Infinity));
        for raw in 0..(1u64 << kc.n) {
            let x = F2mElement::from_biguint(&BigUint::from(raw), kc.n);
            for p in points_with_x(&kc.curve, &x) {
                assert!(keys.insert(pack_point(&p)), "collision at {p:?}");
                count += 1;
            }
        }
        assert_eq!(BigUint::from(count), kc.group_order);
    }

    #[test]
    fn pruned_subspace_bases_keep_every_oracle_honest() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let parent = build_frobenius_factor_base_from_divisor(&kc, &[1, 2]).unwrap();
        assert!(parent.signed_orbits.len() >= 4);
        let keep: Vec<usize> = (0..parent.signed_orbits.len()).step_by(2).collect();
        let fb = restrict_factor_base_to_orbits(&kc, &parent, &keep).unwrap();
        assert_eq!(
            fb.domain,
            FactorBaseDomain::SubspaceSubset {
                retained_orbits: keep.len()
            }
        );
        assert_eq!(fb.signed_orbits.len(), keep.len());
        assert_eq!(fb.subspace_basis, parent.subspace_basis);
        assert!(!fb.uses_ambient_basis());
        let parent_keys: HashSet<_> = parent.points.iter().map(point_key).collect();
        assert!(fb.points.iter().all(|p| parent_keys.contains(&point_key(p))));
        let index = fb.index_map();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let table = PairSumTable::build(&kc, &fb).unwrap();
        let mut found = 0;
        for m in [2usize, 3] {
            for k in 1u32..=40 {
                let target = kc.mul(kc.generator(), &BigUint::from(k));
                let reference = enumerate_decompose(&kc, &fb, &index, &target, m);
                assert_eq!(
                    table.decompose(&kc, &fb, &target, m).is_some(),
                    reference.is_some()
                );
                let (sat, stats) = sat_decompose_with(
                    &kc,
                    &fb,
                    &index,
                    &st,
                    &target,
                    m,
                    256,
                    Some(2),
                    SatDecompositionOptions {
                        conflict_budget: 1_000_000,
                        ..Default::default()
                    },
                );
                assert_eq!(stats.spurious, 0);
                assert!(!stats.exhausted, "m = {m}, k = {k}: {stats:?}");
                assert_eq!(sat.is_some(), reference.is_some(), "m = {m}, k = {k}");
                if let Some(idxs) = sat {
                    // Only retained points may appear.
                    assert!(idxs.iter().all(|&i| i < fb.points.len()));
                    let sum = idxs
                        .iter()
                        .fold(BinaryPoint::Infinity, |s, &i| kc.add(&s, &fb.points[i]));
                    assert_eq!(sum, target);
                    found += 1;
                }
                let (groebner, _) = groebner_decompose(
                    &kc,
                    &fb,
                    &index,
                    &st,
                    &target,
                    m,
                    SolverEngine::default(),
                    20_000,
                );
                assert_eq!(groebner.is_some(), reference.is_some(), "F4 m = {m}, k = {k}");
            }
        }
        assert!(found > 0, "the pruned base must still decompose something");
    }

    #[test]
    fn lex_leader_constraint_admits_exactly_the_sorted_pairs() {
        use crate::cryptanalysis::sat::Solver;
        for width in 1..=4usize {
            let mut solver = Solver::new(2 * width as u32);
            add_lex_leq(&mut solver, 1, width as u32 + 1, width);
            let mut models = HashSet::new();
            loop {
                if solver.solve() != SolveResult::Sat {
                    break;
                }
                let model = solver.model();
                let code = |first: usize| {
                    (0..width).fold(0u32, |acc, t| acc | (u32::from(model[first + t]) << t))
                };
                let (a, b) = (code(0), code(width));
                assert!(a <= b, "width {width}: {a} > {b}");
                assert!(models.insert((a, b)));
                solver.reset_search();
                let clause = (0..2 * width)
                    .map(|i| if model[i] { -((i + 1) as i32) } else { (i + 1) as i32 })
                    .collect();
                solver.add_clause(clause);
            }
            let total = 1u32 << width;
            assert_eq!(models.len() as u32, total * (total + 1) / 2, "width {width}");
        }
    }

    #[test]
    fn symmetry_breaking_changes_no_sat_verdict() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let index = fb.index_map();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        for m in [2usize, 3] {
            for k in 1u32..=24 {
                let target = kc.mul(kc.generator(), &BigUint::from(k));
                let reference = enumerate_decompose(&kc, &fb, &index, &target, m);
                let (out, stats) = sat_decompose_with(
                    &kc,
                    &fb,
                    &index,
                    &st,
                    &target,
                    m,
                    256,
                    Some(2),
                    SatDecompositionOptions {
                        symmetry_breaking: true,
                        conflict_budget: 2_000_000,
                        ..Default::default()
                    },
                );
                assert_eq!(stats.spurious, 0);
                assert!(!stats.exhausted, "m = {m}, k = {k}: {stats:?}");
                assert_eq!(out.is_some(), reference.is_some(), "m = {m}, k = {k}");
                if let Some(idxs) = out {
                    let sum = idxs
                        .iter()
                        .fold(BinaryPoint::Infinity, |s, &i| kc.add(&s, &fb.points[i]));
                    assert_eq!(sum, target);
                }
            }
        }
    }

    #[test]
    fn curves_past_the_old_cap_construct_with_a_consistent_frobenius_eigenvalue() {
        // Only these degrees above 24 have a prime-order subgroup larger
        // than its cofactor; the others are rejected by the constructor.
        for (a, n, r) in [(1u8, 29u32, 42_457u64), (0, 31, 1_439_393)] {
            let kc = KoblitzCurve::new(a, n).expect("usable curve");
            assert_eq!(kc.subgroup_order, BigUint::from(r));
            let g = kc.generator().clone();
            assert_eq!(kc.mul(&g, &kc.subgroup_order), BinaryPoint::Infinity);
            for k in [2u32, 3, 1000] {
                let p = kc.mul(&g, &BigUint::from(k));
                assert_eq!(kc.frobenius(&p), kc.mul(&p, &kc.lambda));
            }
        }
        for (a, n) in [(0u8, 25u32), (1, 27), (1, 33)] {
            assert!(KoblitzCurve::new(a, n).is_none(), "K_{a}/2^{n}");
        }
    }

    #[test]
    fn orbit_generated_admissibility_walk_matches_the_naive_one() {
        for (a, n, indices) in [
            (1u8, 7u32, vec![1usize]),
            (1, 7, vec![0, 1]),
            (0, 9, vec![1]),
            (0, 9, vec![1, 2]),
            (1, 15, vec![2]),
            (1, 15, vec![0, 2]),
            (1, 15, vec![0, 1, 2]),
            (0, 13, vec![0]),
        ] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            let fb = build_frobenius_factor_base_from_divisor(&kc, &indices).unwrap();
            for m in 1..=4 {
                assert_eq!(
                    fb.m_can_decompose(&kc, m),
                    fb.m_can_decompose_reference(&kc, m),
                    "K_{a}/2^{n} divisor {indices:?} m = {m}"
                );
            }
            // And on a nonlinear union, which has more cofactor classes.
            let union = build_frobenius_union_factor_base(
                &kc,
                &[F2mElement::one(n), F2mElement::z(n)],
            )
            .unwrap();
            for m in 1..=4 {
                assert_eq!(
                    union.m_can_decompose(&kc, m),
                    union.m_can_decompose_reference(&kc, m),
                    "K_{a}/2^{n} union m = {m}"
                );
            }
        }
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
            let (by_sat, sat_stats) =
                sat_decompose(&kc, &fb, &index_of, &st, &target, 2, 64, Some(2));
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
            let (out, stats) = sat_decompose(&kc, &fb, &index_of, &st, &target, 2, 64, Some(2));
            assert!(out.is_none(), "[{k}]G must not decompose");
            assert!(stats.refuted, "[{k}]G should be refuted by UNSAT");
            assert!(!stats.exhausted);
            assert_eq!(stats.spurious, 0);
        }
    }

    #[test]
    fn macaulay_rows_may_be_added_but_never_substituted() {
        // The trap this preprocessing has to avoid, pinned as a test.
        //
        // `matrix_f4_f2` skips input polynomials whose degree exceeds
        // the Macaulay degree it is asked for.  So the degree-2 rows of
        // a *cubic* system — every chained m ≥ 3 system — describe only
        // its quadratic part.  Adding them is sound; replacing the
        // system with them drops constraints, and on a refuting
        // instance that turns UNSAT into SAT.
        let kc = KoblitzCurve::new(1, 15).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let g = kc.generator().clone();
        let target = kc.mul(&g, &BigUint::from(6u32));
        let x_r = match &target {
            BinaryPoint::Affine { x, .. } => x.clone(),
            BinaryPoint::Infinity => unreachable!(),
        };
        let sys =
            build_decomposition_system(&fb.subspace_basis, &x_r, &kc.curve.b, 3, &st).unwrap();

        let deg = |p: &crate::cryptanalysis::pq_groebner_f2::F2BoolPoly| {
            p.terms
                .iter()
                .map(|t| t.mask.count_ones())
                .max()
                .unwrap_or(0)
        };
        assert_eq!(
            sys.equations.iter().map(deg).max(),
            Some(3),
            "the chained system is cubic"
        );

        // The mechanism: no degree-2 row can carry a cubic monomial, so
        // the cubic equations are simply absent from that row set.
        let rows = matrix_f4_f2(&sys.equations, sys.n_vars, 2).unwrap();
        assert!(
            rows.iter().all(|r| deg(r) <= 2),
            "a degree-2 Macaulay row cannot contain a cubic term"
        );
        assert!(
            sys.equations.iter().any(|e| deg(e) == 3),
            "and the system does contain cubic equations, which those rows drop"
        );
    }

    #[test]
    fn macaulay_rows_are_implied_by_the_system() {
        // The other half: adding rows is sound.  Every row is an
        // F_2-combination of multiples of the equations, so it vanishes
        // wherever the system does.  Verified exhaustively on a
        // quadratic m = 2 system, which is small enough to scan.
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let g = kc.generator().clone();
        let target = kc.mul(&g, &BigUint::from(11u32));
        let x_r = match &target {
            BinaryPoint::Affine { x, .. } => x.clone(),
            BinaryPoint::Infinity => unreachable!(),
        };
        let sys =
            build_decomposition_system(&fb.subspace_basis, &x_r, &kc.curve.b, 2, &st).unwrap();
        assert_eq!(sys.n_vars, 12);

        let mut solutions = 0;
        for d in [2u32, 3] {
            let rows = matrix_f4_f2(&sys.equations, sys.n_vars, d).unwrap();
            solutions = 0;
            for pt in 0..(1u64 << sys.n_vars) {
                if sys.equations.iter().all(|e| e.eval(pt) == 0) {
                    solutions += 1;
                    assert!(
                        rows.iter().all(|e| e.eval(pt) == 0),
                        "a degree-{d} row rejected a real solution"
                    );
                }
            }
        }
        assert!(solutions > 0, "the instance must have solutions to check");
    }

    #[test]
    fn macaulay_preprocessing_does_not_change_any_verdict() {
        // Whatever it does to the cost, it must not touch the answers.
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let index_of = fb.index_map();
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let g = kc.generator().clone();
        for k in 1..12u32 {
            let target = kc.mul(&g, &BigUint::from(k));
            let (plain, s1) = sat_decompose(&kc, &fb, &index_of, &st, &target, 2, 64, None);
            let (hybrid, s2) = sat_decompose(&kc, &fb, &index_of, &st, &target, 2, 64, Some(2));
            assert_eq!(
                plain.is_some(),
                hybrid.is_some(),
                "verdict changed at [{k}]G"
            );
            assert_eq!(s1.refuted, s2.refuted, "refutation changed at [{k}]G");
            assert_eq!(s1.implied_rows, 0);
            assert_eq!((s1.spurious, s2.spurious), (0, 0));
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
            let (by_sat, stats) = sat_decompose(&kc, &fb, &index_of, &st, &target, 3, 64, Some(2));
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
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let mut signed_fixed_relations = None;
        let mut signed_early_relations = None;
        for (collapse_negation, stop_on_verified_rank) in
            [(false, false), (true, false), (true, true)]
        {
            let opts = KoblitzIcOptions {
                strategy: DecompositionStrategy::Sat,
                collapse_negation,
                stop_on_verified_rank,
                ..KoblitzIcOptions::default()
            };
            let report = koblitz_index_calculus_dlp_with_factor_base(&kc, &q, &fb, &opts).unwrap();
            assert_eq!(report.log, Some(d.clone()));
            assert_eq!(report.collapse_negation, collapse_negation);
            assert_eq!(
                report.orbit_count,
                if collapse_negation {
                    fb.unknowns()
                } else {
                    fb.orbits.len()
                }
            );
            assert!(report.sat_calls > 0, "the SAT path must have run");
            assert_eq!(report.reductions, 0, "no Gröbner work on this path");
            if !stop_on_verified_rank {
                assert_eq!(report.relations, report.orbit_count + opts.extra_relations);
            }
            if collapse_negation && stop_on_verified_rank {
                signed_early_relations = Some(report.relations);
            } else if collapse_negation {
                signed_fixed_relations = Some(report.relations);
            }
        }
        assert!(signed_early_relations.unwrap() < signed_fixed_relations.unwrap());
    }

    #[test]
    fn projected_orbit_columns_recover_unknown_scalar_without_labels() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let d = BigUint::from(53u32);
        let q = kc.mul(kc.generator(), &d);
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let report = koblitz_index_calculus_dlp_with_factor_base(
            &kc,
            &q,
            &fb,
            &KoblitzIcOptions {
                strategy: DecompositionStrategy::Sat,
                allow_direct_relation: false,
                collapse_projected_orbits: true,
                ..KoblitzIcOptions::default()
            },
        )
        .unwrap();
        assert!(report.m_cofactor_admissible);
        assert!(report.collapse_projected_orbits);
        assert_eq!(report.orbit_count, projected_signed_orbit_count(&kc, &fb));
        assert_eq!(report.log, Some(d));
        assert!(!report.direct_relation);
        assert_eq!(report.sat_invalid_models, 0);
    }

    #[test]
    fn parallel_sat_relation_batches_recover_the_same_log() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let d = BigUint::from(53u32);
        let q = kc.mul(kc.generator(), &d);
        let fb = build_frobenius_factor_base(&kc, 0).unwrap();
        let opts = KoblitzIcOptions {
            strategy: DecompositionStrategy::Sat,
            relation_batch_size: 4,
            ..KoblitzIcOptions::default()
        };
        let report = koblitz_index_calculus_dlp_with_factor_base(&kc, &q, &fb, &opts).unwrap();
        assert_eq!(report.log, Some(d));
        assert_eq!(report.relation_batch_size, 4);
        assert!(report.relation_batches < report.trials);
        assert_eq!(report.sat_invalid_models, 0);
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
        // Signed orbits have size up to 2n. Short or self-negative
        // orbits make the measured saving smaller than that ceiling.
        assert!(model.relation_collection_speedup <= 2.0 * kc.n as f64);
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
