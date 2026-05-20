//! # Challenge 23 — Clone MT19937 from 624 outputs
//!
//! MT19937's `temper` step is fully invertible.  Tap 624 consecutive
//! outputs, invert each one to recover the underlying state word,
//! and you have a clone of the RNG.

use crate::cryptopals::challenge21::Mt19937;
use crate::cryptopals::Report;

fn untemper(mut y: u32) -> u32 {
    // Invert y ^= y >> 18.
    y ^= y >> 18;
    // Invert y ^= (y << 15) & 0xefc60000.
    y ^= (y << 15) & 0xefc6_0000;
    // Invert y ^= (y << 7) & 0x9d2c5680.  This one needs an iterative
    // approach: each 7-bit chunk depends on the previous.
    let mut a = y;
    let mask = 0x9d2c_5680;
    let mut result: u32 = 0;
    for i in 0..32 {
        let bit_idx = i;
        let known = if bit_idx < 7 {
            (a >> bit_idx) & 1
        } else {
            let xor = (result >> (bit_idx - 7)) & 1 & ((mask >> bit_idx) & 1);
            ((a >> bit_idx) & 1) ^ xor
        };
        result |= known << bit_idx;
    }
    a = result;
    // Invert y ^= y >> 11.
    let mut r2: u32 = 0;
    for i in (0..32).rev() {
        let bit = (a >> i) & 1;
        if i + 11 < 32 {
            let known = bit ^ ((r2 >> (i + 11)) & 1);
            r2 |= known << i;
        } else {
            r2 |= bit << i;
        }
    }
    r2
}

pub fn clone_mt(outputs: &[u32; 624]) -> Mt19937 {
    let mut state = [0u32; 624];
    for i in 0..624 {
        state[i] = untemper(outputs[i]);
    }
    Mt19937 { state, idx: 624 }
}

pub fn run() -> Report {
    let mut r = Report::new(23, "Clone MT19937 from 624 outputs");
    let mut rng = Mt19937::new(0xC0FFEE);
    let mut outs = [0u32; 624];
    for i in 0..624 {
        outs[i] = rng.next_u32();
    }
    let mut clone = clone_mt(&outs);
    for _ in 0..10 {
        let a = rng.next_u32();
        let b = clone.next_u32();
        assert_eq!(a, b, "clone diverged");
    }
    r.line("Cloned RNG matches source for 10 subsequent outputs.");
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn clone_matches() {
        assert!(super::run().success);
    }
}
