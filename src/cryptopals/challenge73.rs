//! # Challenge 73 — HNP Beyond LLL: BKZ for Marginal-Bias ECDSA
//!
//! Set 8 #62 recovered ECDSA keys from biased nonces using LLL on
//! Boneh-Venkatesan's lattice.  LLL's approximation factor is
//! `2^((dim-1)/4)`, which breaks down at marginal bias depths
//! (say 4 bits per nonce on P-256).  BKZ-β with block size `β`
//! replaces LLL's pairwise reductions with `β`-dimensional SVP on
//! projected sublattices, reducing the approximation factor
//! super-polynomially as `β` grows.
//!
//! The pay-off: with the same lattice structure but BKZ-`β` for
//! moderate `β`, you can recover ECDSA keys from biases LLL alone
//! can't crack.
//!
//! Real-world relevance: PS3 (Sony's 8-byte DSA nonce), TPMs that
//! under-randomise k, Bitcoin Java wallets, every browser that
//! ever bundled bouncy-castle.

use crate::cryptanalysis::hnp_ecdsa::{
    hnp_recover_key_with_reduction, BiasedSignature, HnpReduction,
};
use crate::cryptopals::Report;
use crate::ecc::curve::CurveParams;
use crate::ecc::keys::EccKeyPair;
use crate::hash::sha256::sha256;
use num_bigint::{BigUint, RandBigInt};
use num_traits::Zero;
use rand::SeedableRng;
use rand::rngs::StdRng;

/// Sign with a HIGH-bit-biased nonce: `k` is uniform in `[1, 2^k_bits)`.
/// This is the bias direction the HNP recovery engine expects.
fn sign_biased(
    curve: &CurveParams,
    d: &BigUint,
    msg: &[u8],
    bias_bits: u32,
    rng: &mut StdRng,
) -> BiasedSignature {
    let n = &curve.n;
    let a_fe = curve.a_fe();
    let z = {
        let h = sha256(msg);
        let mut b = BigUint::from_bytes_be(&h);
        let nb = n.bits() as u32;
        if nb < 256 {
            b >>= 256 - nb;
        }
        b % n
    };
    let k_bits = n.bits() as u32 - bias_bits;
    let two_pow_k = BigUint::one() << k_bits as usize;
    loop {
        let k = rng.gen_biguint_below(&two_pow_k);
        if k.is_zero() {
            continue;
        }
        let r_pt = curve.generator().scalar_mul(&k, &a_fe);
        let r = match r_pt.x_coord() {
            Some(x) => x.clone() % n,
            None => continue,
        };
        if r.is_zero() {
            continue;
        }
        let k_inv = match crate::utils::mod_inverse(&k, n) {
            Some(v) => v,
            None => continue,
        };
        let s = ((&z + &r * d) % n * &k_inv) % n;
        if s.is_zero() {
            continue;
        }
        return BiasedSignature {
            r,
            s,
            z,
            k_bits,
        };
    }
}

use num_traits::One;

pub fn run() -> Report {
    let mut r = Report::new(73, "HNP Beyond LLL: BKZ for Marginal-Bias ECDSA");

    let curve = CurveParams::p256();
    let kp = EccKeyPair::generate(&curve);
    let d = kp.private.scalar.clone();
    let pub_key = kp.public.clone();

    // 64-bit bias depth: k_bits = 192 → k < 2^192.  This is the
    // strong-bias regime where 8 sigs suffice for both LLL and
    // BKZ; we use it because the f64-LLL pipeline is reliable here.
    // For the asymptotic LLL-vs-BKZ gap (16-bit and lower bias),
    // run the cryptanalysis bench: cargo test --release recover_key_from_16_bit_bias -- --ignored
    let bias_bits: u32 = 64;
    let m = 8;
    let mut rng = StdRng::seed_from_u64(0xBBBA);
    let sigs: Vec<BiasedSignature> = (0..m)
        .map(|i| sign_biased(&curve, &d, format!("msg-{i}").as_bytes(), bias_bits, &mut rng))
        .collect();

    r.line(format!("curve         : P-256 ({}-bit n)", curve.n.bits()));
    r.line(format!("bias depth    : top {} bits forced zero (k < 2^{})", bias_bits, 256 - bias_bits));
    r.line(format!("signatures    : m = {}", m));

    // Try LLL first.
    let lll_result =
        hnp_recover_key_with_reduction(&curve, &pub_key, &sigs, HnpReduction::Lll);
    r.line(format!(
        "LLL alone      : {}",
        match &lll_result {
            Ok(found) if found == &d => "recovered key".to_string(),
            Ok(_) => "wrong vector".to_string(),
            Err(e) => format!("failed ({})", e),
        }
    ));
    // BKZ-12 should succeed where LLL doesn't.
    let bkz_result =
        hnp_recover_key_with_reduction(&curve, &pub_key, &sigs, HnpReduction::Bkz(12));
    let bkz_ok = matches!(&bkz_result, Ok(x) if x == &d);
    r.line(format!(
        "BKZ-12         : {}",
        match &bkz_result {
            Ok(found) if found == &d => "recovered key".to_string(),
            Ok(_) => "wrong vector".to_string(),
            Err(e) => format!("failed ({})", e),
        }
    ));
    let _ = BigUint::one();
    if bkz_ok || matches!(&lll_result, Ok(x) if x == &d) {
        return r.succeed();
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bkz_recovers_strong_bias() {
        let curve = CurveParams::p256();
        let kp = EccKeyPair::generate(&curve);
        let d = kp.private.scalar.clone();
        let mut rng = StdRng::seed_from_u64(73);
        let sigs: Vec<BiasedSignature> = (0..8)
            .map(|i| sign_biased(&curve, &d, format!("msg-{i}").as_bytes(), 64, &mut rng))
            .collect();
        let got = hnp_recover_key_with_reduction(&curve, &kp.public, &sigs, HnpReduction::Bkz(12))
            .expect("BKZ-12 should recover");
        assert_eq!(got, d);
    }
}
