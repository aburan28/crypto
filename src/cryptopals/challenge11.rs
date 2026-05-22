//! # Challenge 11 — ECB / CBC detection oracle
//!
//! An encryption oracle picks ECB or CBC at random and pre/post-pads
//! the plaintext with 5-10 random bytes.  Detect the mode the oracle
//! used.
//!
//! Detection: feed the oracle a long run of identical bytes (e.g.
//! 16·4 = 64 `A`s).  In ECB, repeated input blocks yield repeated
//! ciphertext blocks; CBC, with its IV chaining, randomises every
//! block.

use crate::cryptopals::challenge10::cbc_encrypt_no_iv_prefix;
use crate::cryptopals::Report;
use crate::symmetric::aes::{encrypt_block, AesKey};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Ecb,
    Cbc,
}

/// The encryption oracle: random key, random IV, random 5-10 byte
/// pre/post pads, and a 50/50 ECB vs CBC pick.  Returns
/// `(ciphertext, actual_mode_for_grading)`.
pub fn oracle(pt: &[u8], rng: &mut StdRng) -> (Vec<u8>, Mode) {
    let mut key = [0u8; 16];
    rng.fill(&mut key);
    let key = AesKey::new(&key).unwrap();
    let pre = rng.gen_range(5..=10);
    let post = rng.gen_range(5..=10);
    let mut msg = Vec::with_capacity(pre + pt.len() + post);
    msg.extend((0..pre).map(|_| rng.gen::<u8>()));
    msg.extend_from_slice(pt);
    msg.extend((0..post).map(|_| rng.gen::<u8>()));
    let use_ecb: bool = rng.gen();
    if use_ecb {
        let padded = crate::cryptopals::low_util::pkcs7_pad(&msg, 16);
        let mut out = Vec::new();
        for chunk in padded.chunks_exact(16) {
            let mut b = [0u8; 16];
            b.copy_from_slice(chunk);
            out.extend_from_slice(&encrypt_block(&b, &key));
        }
        (out, Mode::Ecb)
    } else {
        let mut iv = [0u8; 16];
        rng.fill(&mut iv);
        (cbc_encrypt_no_iv_prefix(&msg, &key, &iv), Mode::Cbc)
    }
}

/// Detect mode by feeding the oracle 64 A's and looking for
/// duplicate ciphertext blocks.
pub fn detect_mode<F: FnMut(&[u8]) -> Vec<u8>>(mut oracle: F) -> Mode {
    let pt = vec![b'A'; 64];
    let ct = oracle(&pt);
    let mut blocks: Vec<&[u8]> = ct.chunks_exact(16).collect();
    blocks.sort();
    blocks.dedup();
    if blocks.len() < ct.len() / 16 {
        Mode::Ecb
    } else {
        Mode::Cbc
    }
}

pub fn run() -> Report {
    let mut r = Report::new(11, "ECB / CBC detection oracle");
    let mut rng = StdRng::seed_from_u64(0xDEADBEEF);
    let trials = 20;
    let mut correct = 0;
    for _ in 0..trials {
        let mut local_rng = StdRng::seed_from_u64(rng.gen());
        let mut actual: Option<Mode> = None;
        let guess = detect_mode(|pt| {
            let (ct, mode) = oracle(pt, &mut local_rng);
            actual = Some(mode);
            ct
        });
        if actual == Some(guess) {
            correct += 1;
        }
    }
    r.line(format!("Correct guesses: {}/{}", correct, trials));
    assert_eq!(correct, trials);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn always_detects_correctly() {
        assert!(super::run().success);
    }
}
