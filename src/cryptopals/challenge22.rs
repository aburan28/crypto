//! # Challenge 22 — Crack an MT19937 seed
//!
//! A program seeds MT19937 with the current Unix timestamp at some
//! point in the recent past, sleeps, emits the first output, and
//! sleeps again.  Recover the seed by brute-forcing every Unix
//! timestamp within a wide window.
//!
//! We simulate the scenario rather than actually sleeping — the
//! attack code is identical.

use crate::cryptopals::challenge21::Mt19937;
use crate::cryptopals::Report;

pub fn crack_seed(now: u32, leak: u32, window: u32) -> Option<u32> {
    (now.saturating_sub(window)..=now).rev().find(|&seed| {
        let mut rng = Mt19937::new(seed);
        rng.next_u32() == leak
    })
}

pub fn run() -> Report {
    let mut r = Report::new(22, "Crack MT19937 seed");
    let now = 1_700_000_000u32;
    let actual_seed = now - 500;
    let leak = Mt19937::new(actual_seed).next_u32();
    let cracked = crack_seed(now, leak, 2000).expect("seed in window");
    r.line(format!("Actual seed   : {}", actual_seed));
    r.line(format!("Recovered seed: {}", cracked));
    assert_eq!(cracked, actual_seed);
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovers_seed_in_window() {
        assert!(super::run().success);
    }
}
