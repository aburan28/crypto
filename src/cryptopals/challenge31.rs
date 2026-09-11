//! # Challenge 31 — HMAC-SHA1 timing attack
//!
//! A "validate signature" endpoint does byte-by-byte HMAC compare,
//! returning as soon as a mismatch is found.  An attacker measures
//! response time per byte to recover the MAC one byte at a time.
//!
//! We simulate the timing oracle in-process with a deterministic
//! "sleep" that's proportional to the number of correct leading
//! bytes — slow enough to be measurable but not actually sleeping.

use crate::cryptopals::Report;
use crate::cryptopals::set8_util::hmac_sha256;

const KEY: &[u8] = b"network-key-31";

/// HMAC-SHA1 of (key, msg) — actually we use HMAC-SHA256 to reuse
/// existing code.  The attack is shape-identical: a 32-byte MAC.
fn server_hmac(msg: &[u8]) -> [u8; 32] {
    hmac_sha256(KEY, msg)
}

/// Timing oracle.  Returns `(was_valid, simulated_elapsed_us)`.
/// The simulated elapsed time is proportional to how many leading
/// bytes of the candidate match the true MAC.
pub fn timing_oracle(msg: &[u8], candidate: &[u8]) -> (bool, u32) {
    let mac = server_hmac(msg);
    let mut matched = 0u32;
    for (i, &cb) in candidate.iter().enumerate() {
        if i < mac.len() && cb == mac[i] {
            matched += 1;
        } else {
            break;
        }
    }
    let elapsed = matched * 50; // 50 µs per matching byte
    (matched as usize == mac.len() && candidate.len() == mac.len(), elapsed)
}

pub fn recover_via_timing(msg: &[u8]) -> [u8; 32] {
    let mut guess = [0u8; 32];
    for pos in 0..32 {
        let mut best = (0u8, 0u32);
        for b in 0..=255u8 {
            guess[pos] = b;
            let (_, t) = timing_oracle(msg, &guess[..pos + 1]);
            if t > best.1 {
                best = (b, t);
            }
        }
        guess[pos] = best.0;
    }
    guess
}

pub fn run() -> Report {
    let mut r = Report::new(31, "HMAC-SHA1 timing attack");
    let msg = b"foo";
    let true_mac = server_hmac(msg);
    let recovered = recover_via_timing(msg);
    r.line(format!("True MAC      : {}", hex::encode(true_mac)));
    r.line(format!("Recovered MAC : {}", hex::encode(recovered)));
    assert_eq!(recovered, true_mac);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_mac_via_timing() {
        assert!(super::run().success);
    }
}
