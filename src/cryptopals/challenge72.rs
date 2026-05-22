//! # Challenge 72 — Bellcore Fault Attack on RSA-CRT
//!
//! RSA-CRT signing:
//!
//! ```text
//!     s_p = m^{d_p} mod p
//!     s_q = m^{d_q} mod q
//!     s   = CRT(s_p, s_q)         (mod N)
//! ```
//!
//! If a transient hardware fault corrupts `s_p` (a glitch, a laser,
//! a Rowhammer), the resulting signature `s'` satisfies:
//!
//! - `s' ≡ s_p'  (mod p)`   ← wrong
//! - `s' ≡ s_q   (mod q)`   ← correct
//!
//! Therefore `s'^e − m ≡ 0 (mod q)` but not `mod p`.  Computing
//! `gcd(s'^e − m, N)` reveals `q`, factoring `N`.
//!
//! One faulty signature + one verification of correctness → full
//! key recovery.  Boneh-DeMillo-Lipton 1997, the eternal reason
//! production RSA-CRT impls recompute and verify before releasing
//! a signature.

use crate::cryptopals::Report;
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_integer::Integer;
use num_traits::One;

fn mod_inv(a: &BigUint, n: &BigUint) -> BigUint {
    let a_i = a.to_bigint().unwrap();
    let n_i = n.to_bigint().unwrap();
    let g = a_i.extended_gcd(&n_i);
    (((g.x % &n_i) + &n_i) % &n_i).to_biguint().unwrap()
}

/// CRT-combine two residues mod `p` and mod `q` into one mod `p·q`.
fn crt(s_p: &BigUint, s_q: &BigUint, p: &BigUint, q: &BigUint) -> BigUint {
    let n = p * q;
    let q_inv_mod_p = mod_inv(q, p);
    let h = ((s_p + p - (s_q % p)) * &q_inv_mod_p) % p;
    (s_q + h * q) % n
}

/// RSA-CRT signature (correct version).
pub fn rsa_crt_sign(m: &BigUint, kp: &crate::asymmetric::rsa::RsaPrivateKey) -> BigUint {
    let s_p = m.modpow(&kp.dp, &kp.p);
    let s_q = m.modpow(&kp.dq, &kp.q);
    crt(&s_p, &s_q, &kp.p, &kp.q)
}

/// Faulty RSA-CRT: flips one bit of `s_p` before CRT.
pub fn rsa_crt_sign_faulty(
    m: &BigUint,
    kp: &crate::asymmetric::rsa::RsaPrivateKey,
    fault_bit: u32,
) -> BigUint {
    let s_p_clean = m.modpow(&kp.dp, &kp.p);
    let s_p_dirty = &s_p_clean ^ (BigUint::one() << fault_bit as usize);
    let s_q = m.modpow(&kp.dq, &kp.q);
    crt(&s_p_dirty, &s_q, &kp.p, &kp.q)
}

/// Bellcore attack: factor `N` from one faulty signature.
pub fn bellcore_factor(
    m: &BigUint,
    s_faulty: &BigUint,
    e: &BigUint,
    n: &BigUint,
) -> Option<BigUint> {
    let lhs = s_faulty.modpow(e, n);
    let diff: BigUint = if lhs >= *m {
        lhs - m
    } else {
        m + n - lhs
    };
    let g = diff.gcd(n);
    if g > BigUint::one() && g < *n {
        Some(g)
    } else {
        None
    }
}

pub fn run() -> Report {
    let mut r = Report::new(72, "Bellcore Fault on RSA-CRT");
    let kp = crate::asymmetric::rsa::RsaKeyPair::generate(512);
    let m = BigUint::from_bytes_be(b"the eagle has landed");
    // Correct signature is fine.
    let s_good = rsa_crt_sign(&m, &kp.private);
    let verify_good = s_good.modpow(&kp.public.e, &kp.public.n) == &m % &kp.public.n;
    r.line(format!("Honest CRT signature verifies   : {}", verify_good));
    // Inject a fault at bit 17.
    let s_bad = rsa_crt_sign_faulty(&m, &kp.private, 17);
    let verify_bad = s_bad.modpow(&kp.public.e, &kp.public.n) == &m % &kp.public.n;
    r.line(format!("Faulty signature verifies       : {}", verify_bad));
    assert!(verify_good);
    assert!(!verify_bad);
    // Now factor N.
    let factor = bellcore_factor(&m, &s_bad, &kp.public.e, &kp.public.n);
    match factor {
        Some(f) => {
            let other = &kp.public.n / &f;
            let ok = (f == kp.private.p && other == kp.private.q)
                || (f == kp.private.q && other == kp.private.p);
            r.line(format!("Recovered prime ({} bits)", f.bits()));
            r.line(format!("Matches one of p, q             : {}", ok));
            if ok {
                return r.succeed();
            }
        }
        None => r.line("GCD didn't expose a factor — fault too small or arithmetic mistake."),
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fault_reveals_factor() {
        let kp = crate::asymmetric::rsa::RsaKeyPair::generate(384);
        let m = BigUint::from_bytes_be(b"hi mom");
        let s = rsa_crt_sign_faulty(&m, &kp.private, 11);
        let f = bellcore_factor(&m, &s, &kp.public.e, &kp.public.n).expect("must factor");
        let other = &kp.public.n / &f;
        assert!(
            (f == kp.private.p && other == kp.private.q)
                || (f == kp.private.q && other == kp.private.p)
        );
    }

    #[test]
    fn honest_signature_verifies() {
        let kp = crate::asymmetric::rsa::RsaKeyPair::generate(384);
        let m = BigUint::from_bytes_be(b"verifies");
        let s = rsa_crt_sign(&m, &kp.private);
        assert_eq!(s.modpow(&kp.public.e, &kp.public.n), m % &kp.public.n);
    }
}
