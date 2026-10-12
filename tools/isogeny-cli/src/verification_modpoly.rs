//! The classical modular polynomial `Φ_ℓ(X, Y) mod p`, computed natively
//! from the `q`-expansion of `j`.
//!
//! `Φ_ℓ(X, Y) = X^{ℓ+1} + Y^{ℓ+1} − X^ℓY^ℓ + Σ c_{ab} X^a Y^b` with
//! `c_{ab} = c_{ba}` and `a, b ≤ ℓ`, and `Φ_ℓ(j(q), j(q^ℓ)) = 0`.  The
//! symmetric unknown `c_{ab}` (`a ≤ b`) first appears at `q^{−(a + ℓb)}`
//! with coefficient one, and those pole orders are distinct, so the
//! system is triangular: walking the exponents from `q^{−ℓ(ℓ+1)}` upward
//! determines each unknown at its own exponent without a division.  Every
//! other exponent, and a band of positive ones, must then vanish; they are
//! checked, so a wrong assumption about the form of `Φ_ℓ` or a wrong `j`
//! series fails the construction instead of yielding a wrong polynomial.
//! Because no division happens, the result is the reduction of the
//! integral `Φ_ℓ` for every `p`.
//!
//! `j = 1728·E4³ / (E4³ − E6²)` is computed mod `p` from the divisor sums
//! in `E4 = 1 + 240 Σ σ₃(n) qⁿ` and `E6 = 1 − 504 Σ σ₅(n) qⁿ`.

use super::field::{Fe, Field};
use super::poly::{self, Poly};

/// Extra positive exponents checked beyond those that determine `Φ_ℓ`.
const CHECK_BAND: usize = 8;

/// `Φ_ℓ mod p` as a dense `(ℓ+2) × (ℓ+2)` coefficient table.
#[derive(Clone, Debug)]
pub struct ModPoly {
    pub ell: u64,
    /// `c[a][b]`, the coefficient of `X^a Y^b`.
    pub c: Vec<Vec<Fe>>,
    /// Exponents whose coefficient was required to vanish and did.
    pub checks_passed: usize,
}

/// `J_k` with `j(q) = q^{-1} Σ_{k≥0} J_k q^k`, for `k < len`.
pub fn j_series(f: &Field, len: usize) -> Vec<Fe> {
    let n = len + 1;
    let sigma = |k: u32, m: usize| -> u128 {
        (1..=m)
            .filter(|d| m.is_multiple_of(*d))
            .map(|d| (d as u128).pow(k))
            .sum()
    };
    let mut e4 = vec![f.zero(); n];
    let mut e6 = vec![f.zero(); n];
    e4[0] = f.one();
    e6[0] = f.one();
    for m in 1..n {
        e4[m] = f.mul(&f.from_u64(240), &f.from_u128(sigma(3, m)));
        e6[m] = f.neg(&f.mul(&f.from_u64(504), &f.from_u128(sigma(5, m))));
    }
    let tr = |a: Poly| -> Poly {
        let mut a = a;
        a.resize(n, f.zero());
        a.truncate(n);
        a
    };
    let e4_3 = tr(poly::mul(f, &tr(poly::mul(f, &e4, &e4)), &e4));
    let e6_2 = tr(poly::mul(f, &e6, &e6));
    // E4³ − E6² = 1728 q Π(1 − qⁿ)^24 has no constant term.
    let diff: Vec<Fe> = (0..n).map(|i| f.sub(&e4_3[i], &e6_2[i])).collect();
    debug_assert!(f.is_zero(&diff[0]));
    let den: Vec<Fe> = diff[1..].to_vec(); // (E4³ − E6²)/q, constant 1728
                                           // Series inverse of den to `len` terms.
    let inv0 = f.inv(&den[0]).expect("1728 is a unit for p > 3");
    let mut inv = vec![f.zero(); len];
    inv[0] = inv0;
    for k in 1..len {
        let mut acc = f.zero();
        for i in 1..=k.min(den.len() - 1) {
            acc = f.add(&acc, &f.mul(&den[i], &inv[k - i]));
        }
        inv[k] = f.neg(&f.mul(&acc, &inv0));
    }
    let num = poly::mul(f, &e4_3[..len].to_vec(), &inv);
    let k1728 = f.from_u64(1728);
    (0..len)
        .map(|k| f.mul(&k1728, num.get(k).unwrap_or(&f.zero())))
        .collect()
}

/// Truncated product of two series.
fn mul_trunc(f: &Field, a: &[Fe], b: &[Fe], len: usize) -> Vec<Fe> {
    let mut out = vec![f.zero(); len];
    for (i, ai) in a.iter().enumerate().take(len) {
        if f.is_zero(ai) {
            continue;
        }
        for (j, bj) in b.iter().enumerate().take(len - i) {
            out[i + j] = f.add(&out[i + j], &f.mul(ai, bj));
        }
    }
    out
}

/// `Φ_ℓ mod p` for a prime `ℓ < p`.  Returns `Err` if any check fails.
pub fn modular_polynomial(f: &Field, ell: u64) -> Result<ModPoly, String> {
    let l = ell as usize;
    let top = l * (l + 1); // the deepest pole, of Y^{ℓ+1} and X^ℓY^ℓ
    let len = top + CHECK_BAND + 1; // exponents −top ..= CHECK_BAND
    let jser = j_series(f, len);
    // A[a] = (Σ J_k q^k)^a, B[b] = (Σ J_k q^{ℓk})^b, both truncated.
    let mut apow = vec![{
        let mut one = vec![f.zero(); len];
        one[0] = f.one();
        one
    }];
    for a in 1..=l + 1 {
        let next = mul_trunc(f, &apow[a - 1], &jser, len);
        apow.push(next);
    }
    let mut jl = vec![f.zero(); len];
    for k in 0..len {
        if k * l < len {
            jl[k * l] = jser[k];
        }
    }
    let mut bpow = vec![apow[0].clone()];
    for b in 1..=l + 1 {
        let next = mul_trunc(f, &bpow[b - 1], &jl, len);
        bpow.push(next);
    }
    // R[e + top] accumulates the coefficient of q^e, e ∈ [−top, CHECK_BAND].
    let mut r = vec![f.zero(); len];
    // X^aY^b = q^{−pole} · A[a]·B[b]; its q^e lands at index e + top, so
    // only the first `pole + CHECK_BAND + 1` terms of the product matter.
    // B[b] is a series in q^ℓ: iterate its nonzero terms only.
    let add_monomial = |r: &mut Vec<Fe>, a: usize, b: usize, coef: &Fe| {
        let pole = a + l * b;
        let need = (pole + CHECK_BAND + 1).min(len);
        let base = top - pole;
        for i in (0..need).step_by(l) {
            let bi = bpow[b][i];
            if f.is_zero(&bi) {
                continue;
            }
            let cb = f.mul(coef, &bi);
            for (k, ak) in apow[a].iter().enumerate().take(need - i) {
                let idx = base + i + k;
                r[idx] = f.add(&r[idx], &f.mul(&cb, ak));
            }
        }
    };
    let one = f.one();
    let minus_one = f.neg(&one);
    add_monomial(&mut r, l + 1, 0, &one);
    add_monomial(&mut r, 0, l + 1, &one);
    add_monomial(&mut r, l, l, &minus_one);
    // The unknown whose leading pole is P, if any.
    let mut by_pole: Vec<Option<(usize, usize)>> = vec![None; top + 1];
    for b in 0..=l {
        for a in 0..=b {
            if a == l && b == l {
                continue;
            }
            let pole = a + l * b;
            assert!(by_pole[pole].is_none(), "pole orders are distinct");
            by_pole[pole] = Some((a, b));
        }
    }
    let mut c = vec![vec![f.zero(); l + 2]; l + 2];
    c[l + 1][0] = one;
    c[0][l + 1] = one;
    c[l][l] = minus_one;
    let mut checks = 0usize;
    for idx in 0..len {
        let pole = top as isize - idx as isize; // e = −pole
        let slot = if pole >= 0 {
            by_pole[pole as usize]
        } else {
            None
        };
        match slot {
            Some((a, b)) => {
                let coef = f.neg(&r[idx]);
                c[a][b] = coef;
                c[b][a] = coef;
                add_monomial(&mut r, a, b, &coef);
                if a != b {
                    add_monomial(&mut r, b, a, &coef);
                }
                if !f.is_zero(&r[idx]) {
                    return Err(format!(
                        "Φ_{ell}: unknown ({a},{b}) left q^{} nonzero",
                        -pole
                    ));
                }
            }
            None => {
                if !f.is_zero(&r[idx]) {
                    return Err(format!("Φ_{ell}: q^{} does not vanish", -pole));
                }
                checks += 1;
            }
        }
    }
    Ok(ModPoly {
        ell,
        c,
        checks_passed: checks,
    })
}

impl ModPoly {
    /// `Φ_ℓ(x, Y)` as a polynomial in `Y`.
    pub fn at_x(&self, f: &Field, x: &Fe) -> Poly {
        let n = self.c.len();
        let mut xp = vec![f.one(); n];
        for i in 1..n {
            xp[i] = f.mul(&xp[i - 1], x);
        }
        let out = (0..n)
            .map(|b| {
                (0..n).fold(f.zero(), |acc, a| {
                    f.add(&acc, &f.mul(&self.c[a][b], &xp[a]))
                })
            })
            .collect();
        poly::trim(f, out)
    }

    /// `∂^{dx+dy} Φ / ∂X^{dx} ∂Y^{dy}` at `(x, y)`, for `dx + dy ≤ 2`.
    pub fn partial(&self, f: &Field, x: &Fe, y: &Fe, dx: u32, dy: u32) -> Fe {
        let n = self.c.len();
        let falling =
            |i: usize, d: u32| -> u64 { (0..d as usize).map(|t| (i - t) as u64).product::<u64>() };
        let mut acc = f.zero();
        for a in (dx as usize)..n {
            let xa = f.mul(
                &f.from_u64(falling(a, dx)),
                &f.pow(x, &((a - dx as usize) as u64).into()),
            );
            for b in (dy as usize)..n {
                if f.is_zero(&self.c[a][b]) {
                    continue;
                }
                let yb = f.mul(
                    &f.from_u64(falling(b, dy)),
                    &f.pow(y, &((b - dy as usize) as u64).into()),
                );
                acc = f.add(&acc, &f.mul(&self.c[a][b], &f.mul(&xa, &yb)));
            }
        }
        acc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::modular_polynomial::{phi_2, phi_3};
    use num_bigint::{BigInt, BigUint};
    use num_traits::Signed;

    fn reduce(f: &Field, v: &BigInt) -> Fe {
        let m = f.from_big(&v.abs().to_biguint().unwrap());
        if v.is_negative() {
            f.neg(&m)
        } else {
            m
        }
    }

    #[test]
    fn j_series_matches_the_integer_expansion() {
        let f = Field::new(&BigUint::from(1_000_000_007u64)).unwrap();
        let j = j_series(&f, 4);
        // j = q^-1 + 744 + 196884 q + 21493760 q^2 + 864299970 q^3.
        let want = [1u64, 744, 196_884, 21_493_760, 864_299_970];
        for (k, w) in want.iter().take(4).enumerate() {
            assert_eq!(f.to_big(&j[k]), BigUint::from(*w) % 1_000_000_007u64);
        }
    }

    #[test]
    fn phi2_and_phi3_match_the_literature_tables() {
        let bits = super::super::field::FIELD_BITS;
        let p = if bits == 64 {
            (BigUint::from(1u8) << 61) - 1u8
        } else if bits == 128 {
            (BigUint::from(1u8) << 127) - 1u8
        } else if bits == 192 {
            (BigUint::from(1u8) << 192) - (BigUint::from(1u8) << 64) - 1u8
        } else {
            crate::catalog::registry()["curves"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|r| r["family"] == "prime")
                .map(|r| {
                    BigUint::parse_bytes(
                        crate::catalog::integer(&r["params"]["p"])
                            .to_string()
                            .as_bytes(),
                        10,
                    )
                    .unwrap()
                })
                .filter(|p| p.bits() <= bits as u64)
                .max_by_key(BigUint::bits)
                .unwrap()
        };
        let f = Field::new(&p).unwrap();
        for table in [phi_2(), phi_3()] {
            let m = modular_polynomial(&f, table.l).unwrap();
            let l = table.l as u32;
            for a in 0..=l + 1 {
                for b in 0..=l + 1 {
                    assert_eq!(
                        m.c[a as usize][b as usize],
                        reduce(&f, &table.coeff(a, b)),
                        "Φ_{l} coefficient ({a},{b})"
                    );
                }
            }
        }
    }
}
