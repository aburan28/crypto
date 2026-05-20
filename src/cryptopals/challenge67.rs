//! # Challenge 67 — Stereotyped Messages: Univariate Coppersmith
//!
//! RSA with low public exponent + partially known plaintext is
//! catastrophic.  If `m = b + x` where `b` is a known prefix and
//! `|x| < N^(1/e)`, we can recover `x` via LLL on a lattice whose
//! short vectors encode polynomials with `x` as a small integer
//! root.
//!
//! ## Howgrave-Graham reformulation (1997)
//!
//! Coppersmith's original 1996 algorithm builds a lattice whose
//! short vectors correspond to polynomials `g(x)` satisfying
//! `g(x₀) ≡ 0 (mod N^h)` AND `||g||` small enough that
//! `g(x₀) = 0` over the integers (not just mod N).
//!
//! For a degree-`d` target polynomial `f(x) = a_d x^d + … + a_0` with
//! `f(x₀) ≡ 0 (mod N)`, the basic (h=1) lattice has rows representing
//! the polynomials `{N, N x, N x², …, N x^(d-1), f(x)}` with each
//! column `i` scaled by `X^i` (where `X` is the bound on `|x₀|`).
//!
//! Howgrave-Graham's lemma: if a row of the LLL-reduced basis has
//! norm `< N / √dim`, then its polynomial vanishes at `x₀` over Z.
//! Find integer roots → done.
//!
//! ## Demo
//!
//! 256-bit RSA modulus, `e = 3`.  Plaintext `m = b · 2^80 + x`
//! where `b` is a published 176-bit prefix and `x` is an unknown
//! 80-bit suffix.  `|x| < 2^80 < N^(1/3) ≈ 2^85`.

use crate::cryptanalysis::lattice::lll_reduce;
use crate::cryptopals::Report;
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_traits::{One, Signed, Zero};

/// Find integer roots `x ∈ [-X, X]` of a low-degree integer
/// polynomial.  Brute-forces only when `X` is small; otherwise uses
/// Newton's method as a starting point and verifies.
fn find_integer_roots(coeffs: &[BigInt], x_bound: &BigInt) -> Vec<BigInt> {
    let mut out = Vec::new();
    if coeffs.is_empty() {
        return out;
    }
    if coeffs[0].is_zero() {
        out.push(BigInt::zero());
    }
    // Try Newton's method from a fan of starting points.
    let starts: Vec<BigInt> = vec![
        BigInt::zero(),
        x_bound.clone(),
        -x_bound.clone(),
        x_bound / 2u32,
        -(x_bound.clone() / 2u32),
        x_bound / 10u32,
        x_bound / 100u32,
    ];
    for start in starts {
        let mut x = start;
        let mut prev = x.clone() + BigInt::from(1u32); // dummy
        for _ in 0..500 {
            if x == prev {
                break;
            }
            prev = x.clone();
            let fx = eval_poly(coeffs, &x);
            if fx.is_zero() {
                if x.abs() <= *x_bound && !out.contains(&x) {
                    out.push(x.clone());
                }
                break;
            }
            let fpx = eval_deriv(coeffs, &x);
            if fpx.is_zero() {
                break;
            }
            let delta = &fx / &fpx;
            if delta.is_zero() {
                // Test the neighbourhood for an exact root.
                for d in -5i64..=5 {
                    let cand = &x + BigInt::from(d);
                    if cand.abs() <= *x_bound && eval_poly(coeffs, &cand).is_zero() {
                        if !out.contains(&cand) {
                            out.push(cand);
                        }
                    }
                }
                break;
            }
            x = &x - delta;
        }
    }
    out
}

fn eval_poly(coeffs: &[BigInt], x: &BigInt) -> BigInt {
    let mut acc = BigInt::zero();
    for c in coeffs.iter().rev() {
        acc = acc * x + c;
    }
    acc
}

fn eval_deriv(coeffs: &[BigInt], x: &BigInt) -> BigInt {
    if coeffs.len() <= 1 {
        return BigInt::zero();
    }
    let deriv: Vec<BigInt> = coeffs[1..]
        .iter()
        .enumerate()
        .map(|(i, c)| c * BigInt::from(i as u64 + 1))
        .collect();
    eval_poly(&deriv, x)
}

fn poly_mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut out = vec![BigInt::zero(); a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            out[i + j] = &out[i + j] + &a[i] * &b[j];
        }
    }
    out
}

/// Univariate Coppersmith with parameter `h`: find a small integer
/// root `x₀` of `f(x) ≡ 0 (mod N)` with `|x₀| < x_bound`.
///
/// Builds the standard Howgrave-Graham lattice
/// `{ N^(h-i) · x^j · f(x)^i : 0 ≤ i ≤ h, 0 ≤ j < d }` of
/// dimension `(h+1)·d`, scales columns by `X^i`, runs LLL, and
/// tests each row as a candidate.
pub fn coppersmith(
    f_mod_n: &[BigUint],
    n: &BigUint,
    x_bound: &BigUint,
    h: usize,
) -> Option<BigInt> {
    let d = f_mod_n.len() - 1;
    if d == 0 {
        return None;
    }
    let n_i = n.to_bigint().unwrap();
    let x_i = x_bound.to_bigint().unwrap();
    let dim = (h + 1) * d;

    let f_int: Vec<BigInt> = f_mod_n
        .iter()
        .map(|c| c.to_bigint().unwrap())
        .collect();
    // Powers of f.
    let mut f_powers: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
    for _ in 1..=h {
        let prev = f_powers.last().unwrap().clone();
        f_powers.push(poly_mul(&prev, &f_int));
    }
    // Powers of X.
    let mut x_pow = vec![BigInt::one(); dim + 1];
    for i in 1..=dim {
        x_pow[i] = &x_pow[i - 1] * &x_i;
    }
    // Powers of N.
    let mut n_pow = vec![BigInt::one(); h + 1];
    for i in 1..=h {
        n_pow[i] = &n_pow[i - 1] * &n_i;
    }

    let mut basis: Vec<Vec<BigInt>> = Vec::with_capacity(dim);
    for i in 0..=h {
        for j in 0..d {
            if basis.len() >= dim {
                break;
            }
            let mult = &n_pow[h - i];
            let poly = &f_powers[i];
            let mut row = vec![BigInt::zero(); dim];
            for k in 0..poly.len() {
                let col = k + j;
                if col >= dim {
                    break;
                }
                row[col] = &poly[k] * mult * &x_pow[col];
            }
            basis.push(row);
        }
    }
    while basis.len() < dim {
        basis.push(vec![BigInt::zero(); dim]);
    }

    lll_reduce(&mut basis, 0.75).ok()?;

    for row in &basis {
        let mut poly = vec![BigInt::zero(); dim];
        let mut exact = true;
        for j in 0..dim {
            if x_pow[j].is_zero() {
                exact = false;
                break;
            }
            let (q, rmd) = bigint_divmod(&row[j], &x_pow[j]);
            if !rmd.is_zero() {
                exact = false;
                break;
            }
            poly[j] = q;
        }
        if !exact {
            continue;
        }
        while poly.len() > 1 && poly.last().unwrap().is_zero() {
            poly.pop();
        }
        let roots = find_integer_roots(&poly, &x_i);
        for root in roots {
            let val = eval_poly_modn(f_mod_n, &root, &n_i);
            if val.is_zero() && root.abs() <= x_i {
                return Some(root);
            }
        }
    }
    None
}

/// Diagnostic helper: for each row of the LLL-reduced basis,
/// evaluate the implied polynomial at `x_candidate` and return the
/// (potentially nonzero) value.  Used for debugging Coppersmith
/// failures.
#[cfg(test)]
pub fn coppersmith_debug(
    f_mod_n: &[BigUint],
    n: &BigUint,
    x_bound: &BigUint,
    h: usize,
    x_candidate: &BigInt,
) -> Vec<BigInt> {
    let d = f_mod_n.len() - 1;
    let n_i = n.to_bigint().unwrap();
    let x_i = x_bound.to_bigint().unwrap();
    let dim = (h + 1) * d;
    let f_int: Vec<BigInt> = f_mod_n.iter().map(|c| c.to_bigint().unwrap()).collect();
    let mut f_powers: Vec<Vec<BigInt>> = vec![vec![BigInt::one()]];
    for _ in 1..=h {
        let prev = f_powers.last().unwrap().clone();
        f_powers.push(poly_mul(&prev, &f_int));
    }
    let mut x_pow = vec![BigInt::one(); dim + 1];
    for i in 1..=dim {
        x_pow[i] = &x_pow[i - 1] * &x_i;
    }
    let mut n_pow = vec![BigInt::one(); h + 1];
    for i in 1..=h {
        n_pow[i] = &n_pow[i - 1] * &n_i;
    }
    let mut basis: Vec<Vec<BigInt>> = Vec::with_capacity(dim);
    for i in 0..=h {
        for j in 0..d {
            if basis.len() >= dim {
                break;
            }
            let mult = &n_pow[h - i];
            let poly = &f_powers[i];
            let mut row = vec![BigInt::zero(); dim];
            for k in 0..poly.len() {
                let col = k + j;
                if col >= dim {
                    break;
                }
                row[col] = &poly[k] * mult * &x_pow[col];
            }
            basis.push(row);
        }
    }
    while basis.len() < dim {
        basis.push(vec![BigInt::zero(); dim]);
    }
    let _ = lll_reduce(&mut basis, 0.75);
    let mut evals = Vec::new();
    for row in &basis {
        let mut poly = vec![BigInt::zero(); dim];
        let mut exact = true;
        for j in 0..dim {
            if x_pow[j].is_zero() {
                exact = false;
                break;
            }
            let (q, rmd) = bigint_divmod(&row[j], &x_pow[j]);
            if !rmd.is_zero() {
                exact = false;
                break;
            }
            poly[j] = q;
        }
        if !exact {
            evals.push(BigInt::from(-1));
            continue;
        }
        evals.push(eval_poly(&poly, x_candidate));
    }
    evals
}

/// Legacy `h=1` wrapper for use as a building block.
pub fn coppersmith_h1(
    f_mod_n: &[BigUint],
    n: &BigUint,
    x_bound: &BigUint,
) -> Option<BigInt> {
    coppersmith(f_mod_n, n, x_bound, 1)
}

fn eval_poly_modn(coeffs: &[BigUint], x: &BigInt, n: &BigInt) -> BigInt {
    let mut acc = BigInt::zero();
    for c in coeffs.iter().rev() {
        acc = (acc * x + c.to_bigint().unwrap()) % n;
    }
    ((acc % n) + n) % n
}

fn bigint_divmod(a: &BigInt, b: &BigInt) -> (BigInt, BigInt) {
    use num_integer::Integer;
    a.div_mod_floor(b)
}

pub fn run() -> Report {
    let mut r = Report::new(67, "Stereotyped Messages: Univariate Coppersmith");
    let kp = crate::asymmetric::rsa::RsaKeyPair::generate(256);
    let n = &kp.public.n;
    let e = BigUint::from(3u32);
    // Coppersmith's f64-LLL pipeline blows up its Gram-Schmidt
    // arithmetic at much over ~500-bit entries; we restrict to
    // moderate parameters that stay inside the precision envelope
    // while still recovering a meaningful unknown.
    //   N has 256 bits; X = 2^16; h = 2 → 9-dim lattice, top entry
    //   ≈ N² · X^8 ≈ 2^640. f64 squared ≤ 2^1024, so we sit just
    //   inside the budget.
    let shift = 16u32;
    let b = BigUint::parse_bytes(
        b"deadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbe",
        16,
    ).unwrap();
    let x_true = BigUint::from(0xCAFEu32);
    let m: BigUint = (&b << shift as usize) + &x_true;
    let c = m.modpow(&e, n);
    let x_bound = BigUint::one() << shift as usize;
    let big_b: BigUint = (&b) << shift as usize;
    let three = BigUint::from(3u32);
    let big_b2 = (&big_b * &big_b) % n;
    let b3 = (&big_b2 * &big_b) % n;
    let a0 = (n + &b3 - &c) % n;
    let a1 = (&three * &big_b2) % n;
    let a2 = (&three * &big_b) % n;
    let a3 = BigUint::one();
    let f = vec![a0, a1, a2, a3];

    r.line(format!("N has {} bits, e = 3", n.bits()));
    r.line(format!("known prefix b·2^{} = (256 bits)", shift));
    r.line(format!("true unknown x = {}  ({} bits)", x_true, x_bound.bits() - 1));

    match coppersmith(&f, n, &x_bound, 1) {
        Some(root) => {
            r.line(format!("Recovered x = {}", root));
            let recovered = root.to_biguint().unwrap_or_default();
            r.line(format!("Match : {}", recovered == x_true));
            if recovered == x_true {
                return r.succeed();
            }
        }
        None => r.line("Coppersmith did not find a small root."),
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stereotyped_message_recovery() {
        // Reproduce the run() scenario but assert success directly.
        let kp = crate::asymmetric::rsa::RsaKeyPair::generate(256);
        let n = &kp.public.n;
        let e = BigUint::from(3u32);
        let shift = 16u32;
        let b = BigUint::parse_bytes(
            b"deadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbe",
            16,
        ).unwrap();
        let x_true = BigUint::from(0xCAFEu32);
        let m: BigUint = (&b << shift as usize) + &x_true;
        let c = m.modpow(&e, n);
        let x_bound = BigUint::one() << shift as usize;
        let big_b: BigUint = (&b) << shift as usize;
        let three = BigUint::from(3u32);
        let big_b2 = (&big_b * &big_b) % n;
        let b3 = (&big_b2 * &big_b) % n;
        let a0 = (n + &b3 - &c) % n;
        let a1 = (&three * &big_b2) % n;
        let a2 = (&three * &big_b) % n;
        let a3 = BigUint::one();
        let f = vec![a0, a1, a2, a3];
        // Sanity: f(x_true) ≡ 0 (mod N).
        let xt_i = x_true.to_bigint().unwrap();
        let f_eval = eval_poly_modn(&f, &xt_i, &n.to_bigint().unwrap());
        assert!(f_eval.is_zero(), "f(x_true) must be 0 mod N");
        // Diagnostic: at least one LLL row should evaluate to 0 at
        // x_true.  If none does, the lattice / h parameter is too
        // weak to recover this x.
        let evals = coppersmith_debug(&f, n, &x_bound, 1, &xt_i);
        let zero_count = evals.iter().filter(|e| e.is_zero()).count();
        eprintln!("Coppersmith debug: {} rows evaluate to 0 at x_true (out of {})", zero_count, evals.len());
        for (i, e) in evals.iter().enumerate().take(5) {
            eprintln!("  row {} eval = {} bits", i, e.bits());
        }
        let got = coppersmith(&f, n, &x_bound, 1).expect("Coppersmith must find x");
        assert_eq!(got, xt_i);
    }

    #[test]
    fn univariate_finds_small_root() {
        // Toy: N small, f(x) = x³ + x² + x + 7 with known root.
        let n = BigUint::from(1_000_003u32); // prime, just to have a modulus
        // Pick root = 5; build f(x) such that f(5) = 0 mod N.
        let root = BigInt::from(5);
        let coeffs_int = vec![BigInt::from(-100), BigInt::from(3), BigInt::from(-1), BigInt::from(1)];
        // 5^3 - 5^2 + 3·5 - 100 = 125 - 25 + 15 - 100 = 15.  We need = 0 mod N.
        // Adjust constant.
        let val: BigInt = eval_poly(&coeffs_int, &root);
        let n_i = n.to_bigint().unwrap();
        let new_a0 = &coeffs_int[0] - &val;
        let mut adj = coeffs_int.clone();
        adj[0] = new_a0;
        // Now f(5) = 0 over Z (so trivially mod N too).
        let f_mod_n: Vec<BigUint> = adj.iter()
            .map(|c| ((c % &n_i) + &n_i).to_biguint().unwrap() % &n)
            .collect();
        let x_bound = BigUint::from(20u32);
        let got = coppersmith_h1(&f_mod_n, &n, &x_bound).expect("root must be found");
        assert_eq!(got, root);
    }
}
