//! # Challenge 21 — Implement MT19937
//!
//! Standard 32-bit Mersenne Twister.  Round-trip a known seed and
//! emit a few outputs against the reference vectors.

use crate::cryptopals::Report;

pub struct Mt19937 {
    pub state: [u32; 624],
    pub idx: usize,
}

impl Mt19937 {
    pub fn new(seed: u32) -> Self {
        let mut state = [0u32; 624];
        state[0] = seed;
        for i in 1..624 {
            state[i] = (1_812_433_253u32)
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Mt19937 { state, idx: 624 }
    }

    fn twist(&mut self) {
        for i in 0..624 {
            let y = (self.state[i] & 0x8000_0000)
                | (self.state[(i + 1) % 624] & 0x7fff_ffff);
            let mut next = self.state[(i + 397) % 624] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= 0x9908_b0df;
            }
            self.state[i] = next;
        }
        self.idx = 0;
    }

    pub fn next_u32(&mut self) -> u32 {
        if self.idx >= 624 {
            self.twist();
        }
        let mut y = self.state[self.idx];
        self.idx += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }
}

pub fn run() -> Report {
    let mut r = Report::new(21, "Implement MT19937");
    let mut rng = Mt19937::new(5489); // reference vector seed
    let outputs: Vec<u32> = (0..5).map(|_| rng.next_u32()).collect();
    r.line(format!("MT19937(5489)[0..5] = {:?}", outputs));
    // Known reference: first output for seed 5489 is 3499211612.
    assert_eq!(outputs[0], 3499211612);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn reference_vector_matches() {
        assert!(super::run().success);
    }
}
