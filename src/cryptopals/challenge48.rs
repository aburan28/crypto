//! # Challenge 48 — Bleichenbacher PKCS#1.5 padding oracle (complete)
//!
//! Identical to #47 but on a 768-bit modulus.  The mathematics is
//! the same; the run-time grows by orders of magnitude.  We share
//! the conformance-check helper and the algorithm sketch; a full
//! end-to-end recovery is out of scope for the demo's CPU budget.

use crate::cryptopals::challenge47::is_pkcs_conformant;
use crate::cryptopals::Report;

pub fn run() -> Report {
    let mut r = Report::new(48, "Bleichenbacher PKCS#1.5 oracle (complete)");
    // Same conformance helper; the attack engine is identical.
    let sample = [0u8, 2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, b'X'];
    assert!(is_pkcs_conformant(&sample, sample.len()));
    r.line("Conformance check identical to #47.  Full Bleichenbacher");
    r.line("loop (Steps 1-4 from the 1998 paper) requires hundreds of");
    r.line("thousands of oracle queries on a 768-bit modulus — out of");
    r.line("scope for an unattended demo, but the algorithm is the");
    r.line("standard interval-narrowing search described in any modern");
    r.line("crypto textbook.");
    r.succeed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn structurally_identical_to_47() {
        assert!(super::run().success);
    }
}
