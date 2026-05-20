//! # Challenge 61 — Duplicate-Signature Key Selection (DSKS)
//!
//! Given a published `(message, signature)` pair signed by Alice, an
//! attacker constructs a *different* public-key/parameter set that
//! also "verifies" the signature.  Two flavours:
//!
//! ## ECDSA DSKS
//!
//! ECDSA's verification equation `R = u₁·G + u₂·Q` doesn't pin down
//! `(G, Q)` uniquely — it's a single equation in two unknowns.  Eve
//! picks her own `d'`, computes `t = u₁ + u₂·d'`, sets
//! `G' = t⁻¹·R` and `Q' = d'·G'`, and ships `(G', Q')` as her
//! "domain parameters + public key".  Alice's signature verifies
//! against them.  Eve learns `d'` (she chose it), so she has
//! plausible deniability about who signed the message.
//!
//! ## RSA DSKS
//!
//! RSA verification: `s^e ≡ pad(m) mod N`.  Eve wants to pick
//! `(e', N')` such that this still holds, with `s`, `m` fixed.  Eve
//! does this by treating the equation as a discrete-log problem:
//! find the dlog of `pad(m)` in `Z*_{N'}` to the base `s`.  She
//! picks `N' = p'·q'` with `p' − 1` and `q' − 1` both smooth, so
//! Pohlig-Hellman makes the dlog tractable.

use crate::cryptopals::Report;
use crate::cryptopals::set8_util::{crt_combine, crt_pair, parse_big};
use crate::ecc::curve::CurveParams;
use crate::ecc::field::FieldElement;
use crate::ecc::point::Point;
use crate::hash::sha256::sha256;
use num_bigint::{BigInt, BigUint, RandBigInt, ToBigInt};
use num_integer::Integer;
use num_traits::{One, Zero};
use rand::{Rng, SeedableRng};

fn mod_inv(a: &BigUint, n: &BigUint) -> Option<BigUint> {
    let a_i = a.to_bigint().unwrap();
    let n_i = n.to_bigint().unwrap();
    let g = a_i.extended_gcd(&n_i);
    if g.gcd != BigInt::one() {
        return None;
    }
    Some((((g.x % &n_i) + &n_i) % &n_i).to_biguint().unwrap())
}

/// ── Part 1: ECDSA DSKS ───────────────────────────────────────────
///
/// Forge a `(domain_params, Q')` such that Alice's signature
/// `(r, s)` over `msg` verifies under those parameters.

#[derive(Debug, Clone)]
pub struct ForgedEcdsaParams {
    pub curve: CurveParams,
    pub g_prime: Point,
    pub q_prime: Point,
    pub d_prime: BigUint,
}

pub fn forge_ecdsa(
    curve: &CurveParams,
    r: &BigUint,
    s: &BigUint,
    msg: &[u8],
    seed: u64,
) -> ForgedEcdsaParams {
    use rand::rngs::StdRng;
    let n = &curve.n;
    let z = {
        let h = sha256(msg);
        let mut b = BigUint::from_bytes_be(&h);
        if curve.order_bits() < 256 {
            b >>= 256 - curve.order_bits();
        }
        b % n
    };
    let s_inv = mod_inv(s, n).expect("s must be coprime to n");
    let u1 = (&z * &s_inv) % n;
    let u2 = (r * &s_inv) % n;
    // Reconstruct R: any point with x = r works.  We just need
    // R = u1·G + u2·Q where Q = d·G with Alice's secret d.  Since
    // we don't have d, we pick d' ourselves.
    // Crucially, we know R must satisfy x(R) = r — but we get to
    // *choose* which point has that x-coordinate.  Easiest: pick
    // any d', compute Q', then derive G' that makes things work.
    // The cryptopals recipe:
    //   1. Pick d'.
    //   2. t := u1 + u2·d' (mod n).
    //   3. R := u1·G + u2·Q (compute it from the original public key).
    //   4. G' := t⁻¹·R.
    //   5. Q' := d'·G'.
    let mut rng = StdRng::seed_from_u64(seed);
    let d_prime = rng.gen_biguint_below(n);
    // We don't have Alice's Q here — but we don't need it.  The
    // forgery construction uses R = (u1 + u2·d)·G ≡ t·G.  We don't
    // need the *original* d, just a t we can build with our d'.
    //
    // Concretely: we have x(R) = r.  R is fixed by Alice.  We
    // reconstruct R by lifting r to a point on the curve and using
    // it.  Then everything below works.
    let r_point = lift_x_to_point(&curve, r);
    let t = (&u1 + (&u2 * &d_prime) % n) % n;
    let t_inv = mod_inv(&t, n).expect("t must be coprime to n");
    let g_prime = r_point.scalar_mul(&t_inv, &curve.a_fe());
    let q_prime = g_prime.scalar_mul(&d_prime, &curve.a_fe());
    ForgedEcdsaParams {
        curve: curve.clone(),
        g_prime,
        q_prime,
        d_prime,
    }
}

/// Lift an x-coordinate on the curve to *some* point with that x.
/// We don't care which sign of y — only the x matters for the dsks
/// reconstruction.
fn lift_x_to_point(curve: &CurveParams, x: &BigUint) -> Point {
    use crate::cryptanalysis::ec_index_calculus::sqrt_mod_p;
    let p = &curve.p;
    let xb = curve.fe(x.clone());
    // y² = x³ + a x + b
    let x3 = xb.mul(&xb).mul(&xb);
    let ax = curve.a_fe().mul(&xb);
    let b_fe = curve.fe(curve.b.clone());
    let rhs = x3.add(&ax).add(&b_fe);
    let y_val = sqrt_mod_p(&rhs.value, p).expect("x lies on curve");
    Point::Affine {
        x: curve.fe(x.clone()),
        y: curve.fe(y_val),
    }
}

/// Custom verify with caller-supplied generator and order.  Equivalent
/// in shape to a standard ECDSA verify; needed because our existing
/// `verify` hard-codes `curve.G`.
pub fn verify_with_g(
    curve: &CurveParams,
    g: &Point,
    q: &Point,
    msg: &[u8],
    r: &BigUint,
    s: &BigUint,
) -> bool {
    let n = &curve.n;
    if r.is_zero() || s.is_zero() || r >= n || s >= n {
        return false;
    }
    let z = {
        let h = sha256(msg);
        let mut b = BigUint::from_bytes_be(&h);
        if curve.order_bits() < 256 {
            b >>= 256 - curve.order_bits();
        }
        b % n
    };
    let s_inv = match mod_inv(s, n) {
        Some(v) => v,
        None => return false,
    };
    let u1 = (&z * &s_inv) % n;
    let u2 = (r * &s_inv) % n;
    let a_fe = curve.a_fe();
    let p1 = g.scalar_mul(&u1, &a_fe);
    let p2 = q.scalar_mul(&u2, &a_fe);
    let pt = p1.add(&p2, &a_fe);
    match pt.x_coord() {
        None => false,
        Some(x) => &(x % n) == r,
    }
}

/// ── Part 2: RSA DSKS ─────────────────────────────────────────────
///
/// Given Alice's RSA signature `s` over message `m` with public key
/// `(e, N)`, Eve forges `(e', N')` with `s^e' ≡ pad(m) mod N'`.
///
/// We restrict to a *toy* RSA size so Pohlig-Hellman over `p'-1` and
/// `q'-1` is tractable.  Real-world attack against 2048-bit RSA
/// needs a careful smooth-prime search; for the demo we use small
/// smooth primes whose discrete logs we can solve with BSGS.

/// EMSA-PKCS1-v1_5 padding for SHA-256 — minimal stub used here.
/// Returns the BigUint encoding of `0x00 01 FF .. FF 00 || DIGEST_INFO || H`.
pub fn rsa_pkcs1_pad(msg: &[u8], modulus_bytes: usize) -> BigUint {
    let h = sha256(msg);
    // DigestInfo prefix for SHA-256: 19-byte ASN.1 wrapper.
    let digest_info: [u8; 19] = [
        0x30, 0x31, 0x30, 0x0d, 0x06, 0x09, 0x60, 0x86, 0x48, 0x01, 0x65, 0x03, 0x04,
        0x02, 0x01, 0x05, 0x00, 0x04, 0x20,
    ];
    let t_len = digest_info.len() + 32;
    assert!(modulus_bytes >= t_len + 11);
    let pad_len = modulus_bytes - t_len - 3;
    let mut out = Vec::with_capacity(modulus_bytes);
    out.push(0x00);
    out.push(0x01);
    out.extend(std::iter::repeat(0xff).take(pad_len));
    out.push(0x00);
    out.extend_from_slice(&digest_info);
    out.extend_from_slice(&h);
    BigUint::from_bytes_be(&out)
}

/// Pohlig-Hellman discrete log of `y` to base `g` modulo prime `p`
/// where `g` is a primitive root.  Returns `x` with `g^x = y mod p`.
///
/// Uses BSGS in each prime-power subgroup of `p-1`, then CRT.
pub fn pohlig_hellman_dlog(g: &BigUint, y: &BigUint, p: &BigUint) -> Option<BigUint> {
    use crate::cryptopals::set8_util::small_factors;
    let order = p - BigUint::one();
    // Factor order completely (small factors only — the cryptopals
    // demo uses smooth moduli where this works).
    let mut factors: Vec<(BigUint, u32)> = Vec::new();
    let mut n = order.clone();
    let mut pf: u64 = 2;
    while pf < 1u64 << 24 && &n > &BigUint::one() {
        let pf_b = BigUint::from(pf);
        let mut e: u32 = 0;
        while (&n % &pf_b).is_zero() {
            n /= &pf_b;
            e += 1;
        }
        if e > 0 {
            factors.push((pf_b, e));
        }
        pf = if pf == 2 { 3 } else { pf + 2 };
    }
    if !n.is_one() {
        // Residual large factor — can't BSGS that cheaply.
        return None;
    }
    let _ = small_factors;
    // Per prime-power q^e dividing order, solve dlog in subgroup of order q^e.
    let mut residues: Vec<(BigUint, BigUint)> = Vec::new();
    for (q, e) in &factors {
        let q_e = q.modpow(&BigUint::from(*e), &(p * BigUint::from(2u32))); // safe upper
        let cofactor = &order / &q_e;
        let g_sub = g.modpow(&cofactor, p);
        let y_sub = y.modpow(&cofactor, p);
        let x = bsgs(&g_sub, &y_sub, p, &q_e)?;
        residues.push((x, q_e));
    }
    let (x, _) = crt_combine(&residues);
    Some(x)
}

/// Baby-step giant-step in a cyclic subgroup of order `n` ≤ ~2^24.
pub fn bsgs(g: &BigUint, y: &BigUint, p: &BigUint, n: &BigUint) -> Option<BigUint> {
    use std::collections::HashMap;
    let m = n.sqrt() + BigUint::one();
    let m_u = m.to_u64_digits().get(0).copied().unwrap_or(0) as usize;
    let mut table: HashMap<Vec<u8>, u64> = HashMap::with_capacity(m_u);
    let mut step = BigUint::one();
    for j in 0..m_u as u64 {
        table.entry(step.to_bytes_be()).or_insert(j);
        step = (&step * g) % p;
    }
    // factor = g^(-m) mod p
    let g_inv = mod_inv(g, p)?;
    let factor = g_inv.modpow(&m, p);
    let mut gamma = y.clone();
    for i in 0..m_u as u64 {
        if let Some(j) = table.get(&gamma.to_bytes_be()) {
            return Some(BigUint::from(i) * &m + BigUint::from(*j));
        }
        gamma = (&gamma * &factor) % p;
    }
    None
}

/// Find a small smooth prime `p` near `bits_target` with
/// `p - 1 = 2 · q_1 · q_2 · …` where every `q_i < bound`.  Used to
/// build Eve's forged modulus.  Returns `(p, factors_of_p-1)`.
pub fn smooth_prime(bits_target: u32, bound: u64, seed: u64) -> Option<(BigUint, Vec<BigUint>)> {
    use crate::asymmetric::rsa::is_prime;
    use rand::rngs::StdRng;
    let mut rng = StdRng::seed_from_u64(seed);
    let half = bits_target as usize / 2;
    for _ in 0..5000 {
        // Build q-1/2 = 2 · ∏ small primes near target.
        let mut acc = BigUint::from(2u32);
        let mut factors = vec![BigUint::from(2u32)];
        let small_primes: Vec<u64> = (3..bound).filter(|n| is_small_prime(*n)).collect();
        while (acc.bits() as usize) < (half + bits_target as usize - half - 1) {
            let q = small_primes[rng.gen_range(0..small_primes.len())];
            let qb = BigUint::from(q);
            if factors.contains(&qb) {
                continue;
            }
            acc *= &qb;
            factors.push(qb);
            if acc.bits() as usize >= bits_target as usize - 1 {
                break;
            }
        }
        let candidate = &acc + BigUint::one();
        if is_prime(&candidate) {
            return Some((candidate, factors));
        }
    }
    None
}

fn is_small_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n < 4 {
        return true;
    }
    if n % 2 == 0 {
        return false;
    }
    let mut i = 3u64;
    while i * i <= n {
        if n % i == 0 {
            return false;
        }
        i += 2;
    }
    true
}

/// Build an RSA DSKS forgery for `(m, s)` originally signed under
/// `(e, N)`.  Returns `(e', N')` with `s^e' ≡ pad(m) mod N'`.
pub fn forge_rsa(
    m: &[u8],
    s: &BigUint,
    bits_target: u32,
    seed: u64,
) -> Option<(BigUint, BigUint)> {
    // 1. Find two distinct smooth primes p', q' such that p'·q'
    //    has at least `bits_target` bits.
    let (p_prime, _) = smooth_prime(bits_target / 2 + 1, 1 << 16, seed)?;
    let (q_prime, _) = smooth_prime(bits_target / 2 + 1, 1 << 16, seed.wrapping_add(1000))?;
    if p_prime == q_prime {
        return None;
    }
    let n_prime = &p_prime * &q_prime;
    let modulus_bytes = (n_prime.bits() as usize + 7) / 8;
    let pad = rsa_pkcs1_pad(m, modulus_bytes);
    let s_mod_p = s % &p_prime;
    let pad_mod_p = &pad % &p_prime;
    let s_mod_q = s % &q_prime;
    let pad_mod_q = &pad % &q_prime;
    // 2. Solve discrete log in each subgroup.
    let ep = pohlig_hellman_dlog(&s_mod_p, &pad_mod_p, &p_prime)?;
    let eq = pohlig_hellman_dlog(&s_mod_q, &pad_mod_q, &q_prime)?;
    // 3. CRT-combine mod (p'-1)·(q'-1) — but the simplest stitch
    //    is the per-subgroup order.  We CRT on the full p-1, q-1.
    let (e_prime, _) = crt_pair(
        (&ep, &(&p_prime - BigUint::one())),
        (&eq, &(&q_prime - BigUint::one())),
    );
    Some((e_prime, n_prime))
}

pub fn run() -> Report {
    let mut r = Report::new(61, "Duplicate-Signature Key Selection (DSKS)");

    // ── ECDSA portion ──
    let curve = CurveParams::secp256k1();
    let msg = b"crazy flamboyant for the rap enjoyment";
    let _alice = crate::ecc::keys::EccKeyPair::generate(&curve);
    let sig =
        crate::ecc::ecdsa::sign(msg, &_alice.private, &curve);
    r.line("ECDSA DSKS — Alice signs, Eve forges (G', Q'):");
    let forged = forge_ecdsa(&curve, &sig.r, &sig.s, msg, 5);
    let ok = verify_with_g(
        &forged.curve,
        &forged.g_prime,
        &forged.q_prime,
        msg,
        &sig.r,
        &sig.s,
    );
    r.line(format!("  Alice's signature (r, s) verifies under (G', Q'): {}", ok));
    assert!(ok);

    // ── RSA portion ──
    r.line("");
    r.line("RSA DSKS — toy modulus + Pohlig-Hellman:");
    // Alice signs with a small RSA key so the demo runs in seconds.
    let alice_rsa = crate::asymmetric::rsa::RsaKeyPair::generate(256);
    let msg2 = b"crazy flamboyant for the rap enjoyment";
    let s = crate::asymmetric::rsa::rsa_sign(msg2, &alice_rsa.private);
    r.line(format!("Alice (N, e) = ({}, {})", alice_rsa.public.n, alice_rsa.public.e));
    let _ = parse_big;
    match forge_rsa(msg2, &s, 256, 99) {
        Some((e_prime, n_prime)) => {
            let modulus_bytes = (n_prime.bits() as usize + 7) / 8;
            let pad = rsa_pkcs1_pad(msg2, modulus_bytes);
            let lhs = s.modpow(&e_prime, &n_prime);
            let ok = lhs == pad;
            r.line(format!("  Forged N' bits = {}", n_prime.bits()));
            r.line(format!("  s^e' = pad(m) mod N' : {}", ok));
            if ok {
                return r.succeed();
            }
        }
        None => r.line("  Smooth-prime search did not converge in budget."),
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ecdsa_dsks_verifies() {
        let curve = CurveParams::secp256k1();
        let msg = b"hello";
        let kp = crate::ecc::keys::EccKeyPair::generate(&curve);
        let sig = crate::ecc::ecdsa::sign(msg, &kp.private, &curve);
        let forged = forge_ecdsa(&curve, &sig.r, &sig.s, msg, 11);
        assert!(verify_with_g(
            &forged.curve,
            &forged.g_prime,
            &forged.q_prime,
            msg,
            &sig.r,
            &sig.s,
        ));
    }

    #[test]
    #[ignore]
    fn rsa_dsks_verifies() {
        // Slow: smooth-prime search + Pohlig-Hellman.
        let kp = crate::asymmetric::rsa::RsaKeyPair::generate(256);
        let msg = b"hello";
        let s = crate::asymmetric::rsa::rsa_sign(msg, &kp.private);
        let (e_prime, n_prime) = forge_rsa(msg, &s, 256, 1).expect("forgery must succeed");
        let modulus_bytes = (n_prime.bits() as usize + 7) / 8;
        let pad = rsa_pkcs1_pad(msg, modulus_bytes);
        assert_eq!(s.modpow(&e_prime, &n_prime), pad);
    }
}
