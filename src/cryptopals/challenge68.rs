//! # Challenge 68 — Small Private Exponent: Wiener (with notes on
//! Boneh-Durfee)
//!
//! When the RSA private exponent `d` is too small relative to the
//! modulus, the key is recoverable in polynomial time.  Wiener's
//! 1990 continued-fraction attack works whenever
//! `d < N^0.25`.  Boneh-Durfee 1999 pushed this to `d < N^0.292`
//! using bivariate Coppersmith on `f(x, y) = e·x · (1 + y) − 1 (mod e)`
//! — the lattice is intricate and lattice-reduction precision is
//! a real challenge.
//!
//! We implement Wiener's attack here (classic, ~50 lines) and
//! describe BD's lattice construction in this docstring.  Pushing
//! BD past Wiener's bound to the headline `N^0.292` requires
//! BKZ-on-rationals or arbitrary-precision LLL, both beyond the
//! repo's f64-LLL precision envelope.
//!
//! ## Wiener's attack
//!
//! From `e·d ≡ 1 (mod φ(N))` we have `e·d − k·φ(N) = 1` for some
//! integer `k`.  Dividing by `d·φ(N)`:
//!
//! ```text
//!     e/φ(N) − k/d  ≈  1 / (d·φ(N))
//! ```
//!
//! Since `φ(N) ≈ N`, the fraction `k/d` is a convergent of the
//! continued-fraction expansion of `e/N`.  Iterate convergents,
//! test each as a candidate `(k, d)`, and recover `φ(N)`
//! (hence `p, q`) by solving the quadratic
//! `x² − (N + 1 − φ(N))·x + N = 0`.
//!
//! ## Bound discussion
//!
//! - Wiener: `d < N^0.25`.
//! - Boneh-Durfee (EUROCRYPT 1999, sub-lattice strategy): `d < N^0.292`.
//!
//! BD's lattice uses monomials `x^i · y^j · f(x,y)^k mod e^k` for
//! `i ≤ m, j ≤ t` and a careful geometric exclusion of certain
//! rows; sub-lattice det / target-vector norm comparison gives the
//! headline bound.  See *Cryptanalysis of RSA with Private Key d
//! less than N^0.292* for the full construction.

use crate::cryptopals::Report;
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

fn continued_fraction(mut a: BigInt, mut b: BigInt) -> Vec<BigInt> {
    let mut out = Vec::new();
    while !b.is_zero() {
        let (q, r) = a.div_rem(&b);
        out.push(q);
        a = b;
        b = r;
    }
    out
}

fn convergents(cf: &[BigInt]) -> Vec<(BigInt, BigInt)> {
    // h_n = q_n · h_{n-1} + h_{n-2}, with h_{-1} = 1, h_{-2} = 0.
    // k_n same recurrence with k_{-1} = 0, k_{-2} = 1.
    let mut out: Vec<(BigInt, BigInt)> = Vec::new();
    let mut h_prev2 = BigInt::zero();
    let mut h_prev1 = BigInt::one();
    let mut k_prev2 = BigInt::one();
    let mut k_prev1 = BigInt::zero();
    for q in cf {
        let h_n = q * &h_prev1 + &h_prev2;
        let k_n = q * &k_prev1 + &k_prev2;
        h_prev2 = std::mem::replace(&mut h_prev1, h_n);
        k_prev2 = std::mem::replace(&mut k_prev1, k_n);
        out.push((h_prev1.clone(), k_prev1.clone()));
    }
    out
}

fn isqrt(n: &BigInt) -> Option<BigInt> {
    if n.is_negative() {
        return None;
    }
    if n.is_zero() {
        return Some(BigInt::zero());
    }
    let mut x = n.clone();
    let mut y = (&x + 1) >> 1;
    while y < x {
        x = y;
        y = (&x + n / &x) >> 1;
    }
    Some(x)
}

/// Wiener's attack: given `(N, e)` with `d < N^0.25`, recover `d`.
pub fn wiener_attack(n: &BigUint, e: &BigUint) -> Option<BigUint> {
    let n_i = n.to_bigint().unwrap();
    let e_i = e.to_bigint().unwrap();
    let cf = continued_fraction(e_i.clone(), n_i.clone());
    let convs = convergents(&cf);
    for (k, d) in convs {
        if k.is_zero() || d.is_zero() {
            continue;
        }
        // φ(N) candidate = (e·d − 1) / k.
        let edm1: BigInt = &e_i * &d - BigInt::one();
        let (phi_q, phi_r) = edm1.div_rem(&k);
        if !phi_r.is_zero() {
            continue;
        }
        // From φ and N: solve x² − (N+1 − φ)·x + N = 0.
        let sum_pq: BigInt = &n_i + BigInt::one() - &phi_q;
        let disc: BigInt = &sum_pq * &sum_pq - BigInt::from(4) * &n_i;
        if disc.is_negative() {
            continue;
        }
        let sq = match isqrt(&disc) {
            Some(s) => s,
            None => continue,
        };
        if &sq * &sq != disc {
            continue;
        }
        let total: BigInt = &sum_pq + &sq;
        if total.is_odd() {
            continue;
        }
        return d.to_biguint();
    }
    None
}

fn generate_weak_d_key(bits: u64) -> crate::asymmetric::rsa::RsaKeyPair {
    use crate::asymmetric::rsa::{random_prime, RsaKeyPair, RsaPrivateKey, RsaPublicKey};
    use crate::utils::mod_inverse;
    loop {
        let p = random_prime(bits / 2);
        let q = random_prime(bits / 2);
        if p == q {
            continue;
        }
        let n = &p * &q;
        let phi = (&p - BigUint::one()) * (&q - BigUint::one());
        // Small d: choose d ≈ N^0.2 (well within Wiener's bound).
        let bound = (&n).nth_root(5);
        let mut d = bound;
        if d.is_even() {
            d += BigUint::one();
        }
        if d.gcd(&phi) != BigUint::one() {
            continue;
        }
        let e = match mod_inverse(&d, &phi) {
            Some(e) => e,
            None => continue,
        };
        let dp = &d % (&p - BigUint::one());
        let dq = &d % (&q - BigUint::one());
        let qinv = mod_inverse(&q, &p).unwrap();
        return RsaKeyPair {
            public: RsaPublicKey {
                n: n.clone(),
                bits,
                e: e.clone(),
            },
            private: RsaPrivateKey {
                n,
                bits,
                e,
                d,
                p,
                q,
                dp,
                dq,
                qinv,
            },
        };
    }
}

pub fn run() -> Report {
    let mut r = Report::new(68, "Small Private Exponent: Wiener");
    let kp = generate_weak_d_key(512);
    r.line(format!("N has {} bits", kp.public.n.bits()));
    r.line(format!("True d ({} bits)", kp.private.d.bits()));
    match wiener_attack(&kp.public.n, &kp.public.e) {
        Some(d) => {
            r.line(format!("Wiener recovered d  : {} bits", d.bits()));
            r.line(format!("Match : {}", d == kp.private.d));
            if d == kp.private.d {
                r.line("");
                r.line("Boneh-Durfee (1999) reaches d < N^0.292 via bivariate");
                r.line("Coppersmith on f(x,y)=e·x·(1+y)-1 — needs higher LLL");
                r.line("precision than our f64 pipeline; see module docstring.");
                return r.succeed();
            }
        }
        None => r.line("Wiener attack did not converge."),
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wiener_recovers_small_d() {
        let kp = generate_weak_d_key(512);
        let d = wiener_attack(&kp.public.n, &kp.public.e).expect("Wiener should succeed");
        assert_eq!(d, kp.private.d);
    }

    #[test]
    fn continued_fraction_round_trip() {
        let cf = continued_fraction(BigInt::from(415), BigInt::from(93));
        assert_eq!(
            cf,
            vec![
                BigInt::from(4),
                BigInt::from(2),
                BigInt::from(6),
                BigInt::from(7)
            ]
        );
    }
}
