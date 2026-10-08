//! Koblitz curves `K_a: y² + xy = x³ + a x² + 1` over `GF(2^n)` for
//! `63 < n ≤ 127`: the wide family the m = 83 confidence gate needs
//! (AGENTS.md §8a), past the one-word
//! [`KoblitzCurve`](crate::cryptanalysis::koblitz_index_calculus::KoblitzCurve),
//! whose constructor factors `#E` by trial division and packs elements in a
//! `u64`.
//!
//! Nothing here is searched for that could be chosen differently:
//!
//! - **Modulus**: the caller's (ecbench passes the registry's, so the ICV1
//!   slug is the registered one; for m = 83 that is §8a's
//!   `z^83 + z^45 + z^2 + z + 1`).
//! - **Order**: the Koblitz recurrence.  With `t = (−1)^{1−a}` and
//!   `V_0 = 2, V_1 = t, V_k = t·V_{k−1} − 2·V_{k−2}`,
//!   `#E(F_{2^n}) = 2^n + 1 − V_n`.  `#E(F_2)` is `4` for `a = 0` and `2`
//!   for `a = 1`, and divides `#E`; the rest must be a prime `r`
//!   (Miller–Rabin, 32 fixed bases), or the curve is refused.
//! - **Generator**: the first `x = 1, 2, 3, …` (read as a polynomial) whose
//!   `c = x + a + x^{−2}` has trace zero, `y = x·H(c)` with `H` the
//!   half-trace, then `G = [h]P`, which must be nonzero and killed by `r`.
//! - **λ**: the root of `λ² − tλ + 2 ≡ 0 (mod r)` with `φ(G) = [λ]G`, where
//!   `φ(x, y) = (x², y²)`.

use num_bigint::{BigInt, BigUint};
use num_traits::{One, ToPrimitive};

use crate::binary_ecc::IrreduciblePoly;
use crate::cryptanalysis::curve_id::{self, CurveId};
use crate::cryptanalysis::koblitz_strong_rho::{RawPointG, RhoScalar, WideStrongRho};
use crate::cryptanalysis::wide_gf2m::WideGf2;

/// A wide Koblitz curve with its prime subgroup.
#[derive(Clone, Debug)]
pub struct WideKoblitz {
    pub a: u8,
    pub n: u32,
    pub irreducible: IrreduciblePoly,
    /// Frobenius trace over `F_2`: `−1` for `a = 0`, `+1` for `a = 1`.
    pub trace: i64,
    pub group_order: u128,
    pub r: u128,
    pub cofactor: u128,
    /// `φ(Q) = [λ]Q` on `⟨G⟩`.
    pub lambda: u128,
    pub generator: (u128, u128),
    /// The `x` the generator rule stopped at (its `P` gave `G = [h]P`).
    pub generator_seed_x: u128,
}

/// `a·b mod m` through the strong rho's exact `u128` routine.
fn mulm(a: u128, b: u128, m: u128) -> u128 {
    <u128 as RhoScalar>::mul_mod(a, b, m)
}

fn powm(mut b: u128, mut e: u128, m: u128) -> u128 {
    let mut acc = 1u128 % m;
    b %= m;
    while e > 0 {
        if e & 1 == 1 {
            acc = mulm(acc, b, m);
        }
        b = mulm(b, b, m);
        e >>= 1;
    }
    acc
}

/// Miller–Rabin with the first 32 primes as bases: deterministic far past
/// `2^127` for any practical purpose here (the bound for the first 13
/// primes alone is above `3.3·10^24`, and these orders are certified again
/// by the generator having order exactly `r`).
pub fn is_probable_prime(n: u128) -> bool {
    const BASES: [u128; 32] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113, 127, 131,
    ];
    if n < 2 {
        return false;
    }
    for p in BASES {
        if n == p {
            return true;
        }
        if n.is_multiple_of(p) {
            return false;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    'base: for a in BASES {
        let mut x = powm(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulm(x, x, n);
            if x == n - 1 {
                continue 'base;
            }
        }
        return false;
    }
    true
}

/// A square root of `v` modulo the prime `p` (Tonelli–Shanks), if any.
fn sqrt_mod(v: u128, p: u128) -> Option<u128> {
    let v = v % p;
    if v == 0 {
        return Some(0);
    }
    if powm(v, (p - 1) / 2, p) != 1 {
        return None;
    }
    let s = (p - 1).trailing_zeros();
    let q = (p - 1) >> s;
    let mut z = 2u128;
    while powm(z, (p - 1) / 2, p) != p - 1 {
        z += 1;
    }
    let (mut m, mut c, mut t, mut r) = (s, powm(z, q, p), powm(v, q, p), powm(v, q.div_ceil(2), p));
    while t != 1 {
        let mut i = 0;
        let mut t2 = t;
        while t2 != 1 {
            t2 = mulm(t2, t2, p);
            i += 1;
        }
        let b = powm(c, 1u128 << (m - i - 1), p);
        m = i;
        c = mulm(b, b, p);
        t = mulm(t, c, p);
        r = mulm(r, b, p);
    }
    Some(r)
}

/// `#E(F_{2^n})` for `K_a`.
pub fn koblitz_order(a: u8, n: u32) -> BigUint {
    let t = BigInt::from(if a == 0 { -1 } else { 1 });
    let (mut v0, mut v1) = (BigInt::from(2), t.clone());
    for _ in 1..n {
        let v2 = &t * &v1 - BigInt::from(2) * &v0;
        v0 = v1;
        v1 = v2;
    }
    let order: BigInt = (BigInt::one() << n) + BigInt::one() - v1;
    order.to_biguint().expect("a curve order is positive")
}

impl WideKoblitz {
    /// Build `K_a` over the field of `irr` (`64 ≤ n ≤ 127`).
    pub fn new(a: u8, irr: &IrreduciblePoly) -> Result<Self, String> {
        let n = irr.degree;
        if a > 1 || !(64..=127).contains(&n) {
            return Err(format!(
                "a wide Koblitz curve takes a ∈ {{0,1}} and 64 ≤ n ≤ 127, not a={a} n={n}"
            ));
        }
        let field = WideGf2::new(irr);
        let order = koblitz_order(a, n)
            .to_u128()
            .ok_or("the group order does not fit a u128")?;
        let cofactor: u128 = if a == 0 { 4 } else { 2 };
        if order % cofactor != 0 {
            return Err(format!("#E(F_2) = {cofactor} does not divide #E = {order}"));
        }
        let r = order / cofactor;
        if r >= 1u128 << 127 || !is_probable_prime(r) {
            return Err(format!("#E / {cofactor} = {r} is not a usable prime"));
        }
        let trace: i64 = if a == 0 { -1 } else { 1 };
        // λ: roots of λ² − tλ + 2 mod r.
        let t_mod = if trace < 0 { r - 1 } else { 1 };
        let disc = (mulm(t_mod, t_mod, r) + r - 8 % r) % r;
        let root = sqrt_mod(disc, r).ok_or("t² − 8 is not a square mod r")?;
        let inv2 = r.div_ceil(2);
        let roots = [
            mulm((t_mod + root) % r, inv2, r),
            mulm((t_mod + r - root) % r, inv2, r),
        ];

        let mut out = Self {
            a,
            n,
            irreducible: irr.clone(),
            trace,
            group_order: order,
            r,
            cofactor,
            lambda: 0,
            generator: (0, 0),
            generator_seed_x: 0,
        };
        let ops = WideStrongRho::from_parts_unchecked(
            field.clone(),
            u128::from(a),
            r,
            RawPointG::Infinity,
        );
        // The generator rule.
        let mut x = 0u128;
        let g = loop {
            x += 1;
            if x > 1 << 20 {
                return Err("no generator in 2^20 abscissae".into());
            }
            let xi = field.inv(x);
            let c = x ^ u128::from(a) ^ field.sqr(xi);
            let Some(z) = field.solve_quadratic(c) else {
                continue;
            };
            let p = RawPointG::Affine {
                x,
                y: field.mul(x, z),
            };
            let g = ops.scalar_mul(p, cofactor);
            if g == RawPointG::Infinity {
                continue;
            }
            if ops.scalar_mul(g, r) != RawPointG::Infinity {
                return Err("[r][h]P is not the identity: the order is wrong".into());
            }
            break g;
        };
        let RawPointG::Affine { x: gx, y: gy } = g else {
            unreachable!()
        };
        let frob = RawPointG::Affine {
            x: field.sqr(gx),
            y: field.sqr(gy),
        };
        let lambda = roots
            .into_iter()
            .find(|&l| ops.scalar_mul(g, l) == frob)
            .ok_or("neither root of λ² − tλ + 2 is the Frobenius eigenvalue on ⟨G⟩")?;
        out.lambda = lambda;
        out.generator = (gx, gy);
        out.generator_seed_x = x;
        Ok(out)
    }

    pub fn field(&self) -> WideGf2 {
        WideGf2::new(&self.irreducible)
    }

    /// The strong reference (and the group arithmetic every wide method
    /// uses) on this curve.
    pub fn strong_rho(&self) -> WideStrongRho {
        WideStrongRho::from_parts(
            self.field(),
            u128::from(self.a),
            self.r,
            self.lambda,
            RawPointG::Affine {
                x: self.generator.0,
                y: self.generator.1,
            },
        )
    }

    /// The curve's ICV1 identity, by the same rule as every one-word curve.
    pub fn curve_id(&self) -> CurveId {
        curve_id::binary(
            self.n,
            &curve_id::modulus_integer(&self.irreducible),
            &BigUint::from(self.a),
            &BigUint::one(),
            &BigUint::from(self.group_order),
            Some(-7),
        )
        .expect("a constructed curve is non-singular and inside the Hasse interval")
    }

    /// The points with abscissa `x` (none, one when `x = 0`, or two).
    pub fn points_with_x(&self, x: u128) -> Vec<(u128, u128)> {
        let f = self.field();
        if x == 0 {
            // y² = 1.
            return vec![(0, 1)];
        }
        let c = x ^ u128::from(self.a) ^ f.sqr(f.inv(x));
        match f.solve_quadratic(c) {
            Some(z) => {
                let y = f.mul(x, z);
                vec![(x, y), (x, y ^ x)]
            }
            None => vec![],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m83() -> IrreduciblePoly {
        IrreduciblePoly {
            degree: 83,
            low_terms: vec![45, 2, 1, 0],
        }
    }

    #[test]
    fn the_m83_gate_curve_has_the_agents_md_order() {
        let k = WideKoblitz::new(0, &m83()).unwrap();
        assert_eq!(k.group_order, 4 * 2417851639230796216685689u128);
        assert_eq!(k.r, 2417851639230796216685689);
        assert_eq!(k.cofactor, 4);
        assert_eq!(k.curve_id().slug, "icv1-f2m83-tm6151469093347-debefd74");
        // λ has order dividing n and acts as Frobenius.
        assert_eq!(powm(k.lambda, 83, k.r), 1);
    }

    #[test]
    fn the_wide_order_agrees_with_the_one_word_constructor_where_both_exist() {
        use crate::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
        for (a, n) in [(0u8, 41u32), (1, 47), (0, 53), (1, 59), (0, 61)] {
            let kc = KoblitzCurve::new(a, n).unwrap();
            assert_eq!(koblitz_order(a, n), kc.group_order, "a={a} n={n}");
        }
    }

    #[test]
    fn group_law_holds_at_m83() {
        let k = WideKoblitz::new(0, &m83()).unwrap();
        let rho = k.strong_rho();
        let g = rho.generator();
        let a = rho.scalar_mul(g, 123456789);
        let b = rho.scalar_mul(g, 987654321);
        assert_eq!(rho.add(a, b), rho.scalar_mul(g, 123456789 + 987654321));
        assert_eq!(rho.scalar_mul(g, k.r), RawPointG::Infinity);
        assert_eq!(rho.scalar_mul(g, k.r - 1), {
            let RawPointG::Affine { x, y } = g else {
                unreachable!()
            };
            RawPointG::Affine { x, y: y ^ x }
        });
    }

    #[test]
    fn primality_and_square_roots() {
        assert!(is_probable_prime(2417851639230796216685689));
        assert!(!is_probable_prime(2417851639230796216685689 * 3));
        assert!(is_probable_prime((1u128 << 127) - 1));
        let p = 2417851639230796216685689u128;
        for v in [2u128, 3, 5, 12345678901234567890] {
            if let Some(s) = sqrt_mod(v, p) {
                assert_eq!(mulm(s, s, p), v % p);
            }
        }
    }
}
