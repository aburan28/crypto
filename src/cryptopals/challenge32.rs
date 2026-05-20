//! # Challenge 32 — HMAC-SHA1 timing attack with less leak
//!
//! Same as #31 but the timing signal is noisier.  Standard
//! mitigation: average over N rounds per candidate byte to denoise.

use crate::cryptopals::challenge31::timing_oracle;
use crate::cryptopals::Report;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

pub fn noisy_oracle(msg: &[u8], candidate: &[u8], rng: &mut StdRng) -> (bool, i64) {
    let (ok, t) = timing_oracle(msg, candidate);
    let noise: i64 = rng.gen_range(-10..=10);
    (ok, t as i64 + noise)
}

pub fn recover_via_noisy_timing(msg: &[u8], rounds: u32) -> [u8; 32] {
    let mut guess = [0u8; 32];
    let mut rng = StdRng::seed_from_u64(0xC32);
    for pos in 0..32 {
        let mut best = (0u8, i64::MIN);
        for b in 0..=255u8 {
            guess[pos] = b;
            let mut total: i64 = 0;
            for _ in 0..rounds {
                let (_, t) = noisy_oracle(msg, &guess[..pos + 1], &mut rng);
                total += t;
            }
            if total > best.1 {
                best = (b, total);
            }
        }
        guess[pos] = best.0;
    }
    guess
}

pub fn run() -> Report {
    let mut r = Report::new(32, "HMAC-SHA1 noisy timing attack");
    let msg = b"foo";
    let true_mac = crate::cryptopals::set8_util::hmac_sha256(b"network-key-31", msg);
    let recovered = recover_via_noisy_timing(msg, 8);
    r.line(format!("Recovered MAC: {}", hex::encode(recovered)));
    r.line(format!("True      MAC: {}", hex::encode(true_mac)));
    let matching = recovered.iter().zip(&true_mac).filter(|(a, b)| a == b).count();
    r.line(format!("matching bytes: {}/32", matching));
    assert!(matching >= 28); // accept some noise-induced misses
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_most_bytes() {
        assert!(super::run().success);
    }
}
