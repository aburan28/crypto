//! Subfield structure of a binary curve `E_{a,b} : y² + xy = x³ + ax² + b`
//! over `GF(2^n)`.
//!
//! Two degrees, because they differ:
//!
//! - the **j-field degree** `k_j`: the least `k | n` with
//!   `j = 1/b ∈ GF(2^k)`, the field of moduli;
//! - the **definition degree** `k_d`: the least `k | n` such that `E` is
//!   `GF(2^n)`-isomorphic to the base change of a curve over `GF(2^k)`.
//!
//! `E_{a,b} ≅ E_{a',b'}` over `GF(2^n)` iff `b = b'` and
//! `Tr_n(a) = Tr_n(a')`, and for `a' ∈ GF(2^k)`,
//! `Tr_n(a') = (n/k)·Tr_k(a')`.  So `E` descends to `GF(2^k)` iff
//! `b ∈ GF(2^k)` and either `n/k` is odd or `Tr_n(a) = 0`: when `n/k` is
//! even, the quadratic twist (`Tr_n(a) = 1`) does not descend even though
//! its `j` does.  A Koblitz curve has `k_d = 1`; the quadratic twist of
//! one over `GF(2^{2m})` has `k_j = 1` but `k_d = 2`.
//!
//! When `k_d < n` the descended model's trace `t_k` over `GF(2^k)` is a
//! size-independent invariant of the family: `t_n = V_{n/k}(t_k, 2^k)`
//! (a Lucas sequence).  It is found as the unique Hasse-bounded root of
//! that equation, and by counting points over `GF(2^k)` when the root
//! is not unique.  When `n/k` is even both quadratic twists over
//! `GF(2^k)` descend, so only `|t_k|` is an invariant.

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

use super::arith::lucas_v;
use super::Status;
use crate::binary_ecc::{F2mElement, IrreduciblePoly};

/// Above this subfield degree the Lucas root search is not attempted.
pub const LUCAS_MAX_K: u32 = 40;

/// Above this subfield degree an ambiguous root is not resolved by
/// counting points.
pub const COUNT_MAX_K: u32 = 16;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Subfield {
    pub j_field_degree: u32,
    pub definition_degree: u32,
    /// `t_k` of the descended model; `|t_k|` when `sign_free`.
    pub base_trace: Option<BigInt>,
    pub sign_free: bool,
    pub status: Status,
    pub method: &'static str,
}

fn divisors(n: u32) -> Vec<u32> {
    (1..=n).filter(|k| n % k == 0).collect()
}

fn in_subfield(x: &F2mElement, k: u32, irr: &IrreduciblePoly) -> bool {
    x.square_k_times(k, irr) == *x
}

/// `Σ_{i<count} x^{2^{step·i}}`: the absolute trace for `step = 1`,
/// `count = n`; the relative trace `GF(2^n) → GF(2^step)` for
/// `count = n/step`.
fn trace_sum(x: &F2mElement, step: u32, count: u32, irr: &IrreduciblePoly) -> F2mElement {
    let mut acc = x.clone();
    let mut y = x.clone();
    for _ in 1..count {
        y = y.square_k_times(step, irr);
        acc = acc.add(&y);
    }
    acc
}

fn abs_trace_is_one(x: &F2mElement, n: u32, irr: &IrreduciblePoly) -> bool {
    let t = trace_sum(x, 1, n, irr);
    debug_assert!(t.is_zero() || t == F2mElement::one(n));
    !t.is_zero()
}

/// `(k_j, k_d, Tr_n(a) = 1)` for `E_{a,b}/GF(2^n)`: the j-field degree,
/// the definition degree and the twist class.
pub fn degrees(n: u32, irr: &IrreduciblePoly, a: &F2mElement, b: &F2mElement) -> (u32, u32, bool) {
    let divs = divisors(n);
    let j_deg = *divs
        .iter()
        .find(|&&k| in_subfield(b, k, irr))
        .expect("b ∈ GF(2^n)");
    let tr_a = abs_trace_is_one(a, n, irr);
    let def_deg = *divs
        .iter()
        .find(|&&k| k % j_deg == 0 && ((n / k) % 2 == 1 || !tr_a))
        .expect("k = n qualifies");
    (j_deg, def_deg, tr_a)
}

/// The subfield degrees of `E_{a,b}/GF(2^n)` and, when it descends, the
/// trace of the descended model.  `t_n` is the curve's trace over
/// `GF(2^n)`.  `Err` when no Hasse-bounded `t_k` reproduces `t_n`: the
/// recorded trace is then wrong.
pub fn analyse(
    n: u32,
    irr: &IrreduciblePoly,
    a: &F2mElement,
    b: &F2mElement,
    t_n: &BigInt,
) -> Result<Subfield, String> {
    let (j_deg, def_deg, tr_a) = degrees(n, irr, a, b);
    let mut out = Subfield {
        j_field_degree: j_deg,
        definition_degree: def_deg,
        base_trace: None,
        sign_free: false,
        status: Status::NotApplicable,
        method: "no proper subfield of definition",
    };
    if def_deg == n {
        return Ok(out);
    }
    let (k, m) = (def_deg, n / def_deg);
    out.sign_free = m % 2 == 0;
    if k > LUCAS_MAX_K {
        out.status = Status::NotEvaluated;
        out.method = "subfield degree above the Lucas search limit";
        return Ok(out);
    }
    let qk = BigInt::one() << k;
    // |s| ≤ 2√(2^k) ⟺ s² ≤ 2^{k+2}; an ordinary binary curve has odd trace.
    let bound = (BigInt::one() << (k + 2)).sqrt();
    let mut roots: Vec<BigInt> = Vec::new();
    let mut s = -bound.clone();
    while s <= bound {
        if s.is_odd() && lucas_v(&s, &qk, m) == *t_n {
            roots.push(s.clone());
        }
        s += 1;
    }
    if out.sign_free {
        roots.retain(|r| r.is_positive());
    }
    match roots.len() {
        0 => Err(format!(
            "no |t_k| ≤ 2√(2^{k}) has V_{m}(t_k, 2^{k}) = {t_n}: the recorded trace is inconsistent"
        )),
        1 => {
            out.base_trace = roots.pop();
            out.status = Status::Proved;
            out.method = "unique Hasse-bounded root of V_{n/k}(t_k, 2^k) = t_n";
            Ok(out)
        }
        _ if k <= COUNT_MAX_K => {
            let t = descended_trace_by_count(n, k, irr, b, tr_a);
            out.base_trace = Some(if out.sign_free { t.abs() } else { t });
            out.status = Status::Proved;
            out.method = "point count of the descended model over GF(2^k)";
            Ok(out)
        }
        _ => {
            out.status = Status::Unknown;
            out.method = "several Hasse-bounded roots; subfield too large to count";
            Ok(out)
        }
    }
}

/// An `F₂`-basis of `GF(2^k) ⊂ GF(2^n)`, as the span of the relative
/// traces `Tr_{n/k}(z^i)`, with an element of absolute trace one in it.
fn subfield_basis(n: u32, k: u32, irr: &IrreduciblePoly) -> (Vec<F2mElement>, F2mElement) {
    let mut basis: Vec<F2mElement> = Vec::new();
    // Echelon form keyed by leading bit, to test independence.
    let mut echelon: Vec<BigUint> = Vec::new();
    let mut trace_one = None;
    for i in 0..n {
        let c = trace_sum(&F2mElement::from_bit_positions(&[i], n), k, n / k, irr);
        if trace_one.is_none() && abs_trace_is_one(&c, k, irr) {
            trace_one = Some(c.clone());
        }
        let mut v = c.to_biguint();
        for row in &echelon {
            if v.bit(row.bits() - 1) {
                v ^= row;
            }
        }
        if !v.is_zero() {
            echelon.push(v);
            echelon.sort_by_key(|r| std::cmp::Reverse(r.bits()));
            basis.push(c);
        }
        if basis.len() == k as usize && trace_one.is_some() {
            break;
        }
    }
    assert_eq!(basis.len(), k as usize, "Tr_{{n/k}} maps onto GF(2^k)");
    (basis, trace_one.expect("the absolute trace is onto"))
}

/// The trace over `GF(2^k)` of the descended model `E_{a',b}`, by
/// counting its points: `x = 0` gives one point, each other `x` two or
/// none as `Tr_k(x + a' + b/x²)` is `0` or `1`.
/// `tr_a` is `Tr_n(a)` of the curve over `GF(2^n)`.
pub fn descended_trace_by_count(
    n: u32,
    k: u32,
    irr: &IrreduciblePoly,
    b: &F2mElement,
    tr_a: bool,
) -> BigInt {
    let (basis, trace_one) = subfield_basis(n, k, irr);
    let m = n / k;
    // a' ∈ GF(2^k) with (n/k)·Tr_k(a') = Tr_n(a).
    let a_k = if m % 2 == 1 && tr_a {
        trace_one
    } else {
        F2mElement::zero(n)
    };
    let mut points: u64 = 2; // O and (0, √b)
    let mut x = F2mElement::zero(n);
    for gray in 1u64..(1u64 << k) {
        x = x.add(&basis[gray.trailing_zeros() as usize]);
        let inv = x.flt_inverse(irr).expect("x ≠ 0");
        let w = x.add(&a_k).add(&b.mul(&inv.square(irr), irr));
        if !abs_trace_is_one(&w, k, irr) {
            points += 2;
        }
    }
    BigInt::from((1u64 << k) + 1) - BigInt::from(points)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn irr(modulus: u64) -> IrreduciblePoly {
        let degree = 63 - modulus.leading_zeros();
        IrreduciblePoly {
            degree,
            low_terms: (0..degree).filter(|i| modulus >> i & 1 == 1).collect(),
        }
    }

    fn el(v: u64, n: u32) -> F2mElement {
        F2mElement::from_biguint(&BigUint::from(v), n)
    }

    #[test]
    fn koblitz_curves_descend_to_gf2() {
        // K_0 over GF(2^7) (modulus 0x83), trace 13 in the registry.
        let f = irr(0x83);
        let s = analyse(7, &f, &el(0, 7), &el(1, 7), &BigInt::from(13)).unwrap();
        assert_eq!((s.j_field_degree, s.definition_degree), (1, 1));
        assert_eq!(s.base_trace, Some(BigInt::from(-1)));
        assert!(!s.sign_free);
    }

    #[test]
    fn a_gf4_curve_over_gf2_6_descends_to_gf4() {
        // E_{0,ω}/GF(4) over GF(2^6), modulus 0x43, b = 0x3a, trace −11:
        // V_3(t₂, 4) = t₂³ − 12·t₂ = −11 has the one root t₂ = 1.
        let f = irr(0x43);
        let (a, b) = (el(0, 6), el(0x3a, 6));
        let s = analyse(6, &f, &a, &b, &BigInt::from(-11)).unwrap();
        assert_eq!((s.j_field_degree, s.definition_degree), (2, 2));
        assert_eq!(s.base_trace, Some(BigInt::one()));
        assert!(!s.sign_free);
        // The count agrees with the Lucas root.
        assert_eq!(descended_trace_by_count(6, 2, &f, &b, false), BigInt::one());
    }

    #[test]
    fn the_twist_of_a_koblitz_curve_over_an_even_degree_descends_only_to_gf4() {
        // E_{ω,1} over GF(2^18) (modulus 0x40009, a = 0x9208 = ω): j = 1,
        // but Tr_18(ω) = 1 and 18 is even, so GF(2) is not a field of
        // definition.  GF(4) is: E_{ω,1}/GF(4) is the quadratic twist of
        // K_0/GF(4) (t = −3), so t₂ = 3 and t₁₈ = V_9(3, 4) = 999, the
        // negative of K_0's trace over GF(2^18).
        let f = irr(0x40009);
        let (a, b) = (el(0x9208, 18), el(1, 18));
        assert!(in_subfield(&a, 2, &f) && !in_subfield(&a, 1, &f));
        let s = analyse(18, &f, &a, &b, &BigInt::from(999)).unwrap();
        assert_eq!((s.j_field_degree, s.definition_degree), (1, 2));
        assert_eq!(s.base_trace, Some(BigInt::from(3)));
        assert_eq!(
            descended_trace_by_count(18, 2, &f, &b, true),
            BigInt::from(3)
        );
    }

    #[test]
    fn a_generic_curve_has_no_proper_subfield_of_definition() {
        let f = irr(0x201b);
        let s = analyse(13, &f, &el(1, 13), &el(0x1503, 13), &BigInt::from(11)).unwrap();
        assert_eq!((s.j_field_degree, s.definition_degree), (13, 13));
        assert_eq!(s.status, Status::NotApplicable);
    }

    #[test]
    fn an_inconsistent_trace_is_an_error() {
        let f = irr(0x83);
        assert!(analyse(7, &f, &el(0, 7), &el(1, 7), &BigInt::from(15)).is_err());
    }
}
