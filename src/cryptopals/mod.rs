//! # Cryptopals — Set 7 ("Hashes")
//!
//! Solutions to the [Cryptopals](https://cryptopals.com/sets/7) Set 7
//! challenges (49–56).  Each challenge lives in its own file with the
//! attack written from scratch on top of this repo's primitives
//! (AES, MD4, RC4, SHA-256, …) — no external attack code.
//!
//! | # | Title                                      | Attack family                       |
//! |---|--------------------------------------------|-------------------------------------|
//! | 49 | CBC-MAC Message Forgery                   | IV-control + length-extension forge |
//! | 50 | Hashing with CBC-MAC                      | Targeted CBC-MAC pre-image          |
//! | 51 | Compression Ratio Side-Channel (CRIME)    | Length oracle vs. compression       |
//! | 52 | Iterated Hash Multicollisions             | Joux multicollision                 |
//! | 53 | Kelsey–Schneier Expandable Messages       | Long-message second pre-image       |
//! | 54 | Kelsey–Kohno Nostradamus / Herding        | Diamond-structure commitment        |
//! | 55 | MD4 Collisions (Wang)                     | Single-block differential collision |
//! | 56 | RC4 Single-Byte Biases                    | Z16 broadcast bias plaintext recovery |
//!
//! Each module exposes a `run()` function that performs the attack
//! end-to-end and returns a short `Report` summarising what was
//! achieved.  The CLI wires these up so a curious reader can simply
//! say `cargo run -- cryptopals 49` to watch the attack in action.

// ── Set 1 — Basics (1-8) ─────────────────────────────────────────
pub mod low_util;
pub mod challenge1;
pub mod challenge2;
pub mod challenge3;
pub mod challenge4;
pub mod challenge5;
pub mod challenge6;
pub mod challenge7;
pub mod challenge8;

// ── Set 2 — Block crypto (9-16) ─────────────────────────────────
pub mod challenge9;
pub mod challenge10;
pub mod challenge11;
pub mod challenge12;
pub mod challenge13;
pub mod challenge14;
pub mod challenge15;
pub mod challenge16;

// ── Set 3 — Stream/randomness (17-24) ────────────────────────────
pub mod challenge17;
pub mod challenge18;
pub mod challenge19;
pub mod challenge20;
pub mod challenge21;
pub mod challenge22;
pub mod challenge23;
pub mod challenge24;

// ── Set 4 — Streams + crypto APIs (25-32) ───────────────────────
pub mod challenge25;
pub mod challenge26;
pub mod challenge27;
pub mod challenge28;
pub mod challenge29;
pub mod challenge30;
pub mod challenge31;
pub mod challenge32;

// ── Set 5 — DH, MITM, SRP, RSA (33-40) ──────────────────────────
pub mod challenge33;
pub mod challenge34;
pub mod challenge35;
pub mod challenge36;
pub mod challenge37;
pub mod challenge38;
pub mod challenge39;
pub mod challenge40;

// ── Set 6 — RSA, DSA, oracles (41-48) ───────────────────────────
pub mod challenge41;
pub mod challenge42;
pub mod challenge43;
pub mod challenge44;
pub mod challenge45;
pub mod challenge46;
pub mod challenge47;
pub mod challenge48;

// ── Set 9 — "When the Implementation Bites Back" (67-74) ────────
pub mod challenge67;
pub mod challenge68;
pub mod challenge69;
pub mod challenge70;
pub mod challenge71;
pub mod challenge72;
pub mod challenge73;
pub mod challenge74;

// ── Set 7 — Hashes (49-56) ──────────────────────────────────────
pub mod challenge49;
pub mod challenge50;
pub mod challenge51;
pub mod challenge52;
pub mod challenge53;
pub mod challenge54;
pub mod challenge55;
pub mod challenge56;

// ── Set 8 — "Abstract Algebra" ────────────────────────────────────
// Each challenge mirrors the corresponding text at toadstyle.org.
pub mod set8_util;
pub mod challenge57;
pub mod challenge58;
pub mod challenge59;
pub mod challenge60;
pub mod challenge61;
pub mod challenge62;
pub mod challenge63;
pub mod challenge64;
pub mod challenge65;
pub mod challenge66;

/// Outcome of running one challenge.  Kept tiny: a one-line success
/// statement plus a multi-line transcript that the CLI prints.
#[derive(Debug, Clone)]
pub struct Report {
    pub challenge: u32,
    pub title: &'static str,
    pub success: bool,
    pub transcript: String,
}

impl Report {
    pub fn new(challenge: u32, title: &'static str) -> Self {
        Report {
            challenge,
            title,
            success: false,
            transcript: String::new(),
        }
    }

    pub fn line(&mut self, s: impl AsRef<str>) {
        self.transcript.push_str(s.as_ref());
        self.transcript.push('\n');
    }

    pub fn succeed(mut self) -> Self {
        self.success = true;
        self
    }
}

/// Dispatch table.  `run(n)` runs challenge `n`; returns `None` if
/// `n` is outside 49..=56.
pub fn run(n: u32) -> Option<Report> {
    Some(match n {
        1 => challenge1::run(),
        2 => challenge2::run(),
        3 => challenge3::run(),
        4 => challenge4::run(),
        5 => challenge5::run(),
        6 => challenge6::run(),
        7 => challenge7::run(),
        8 => challenge8::run(),
        9 => challenge9::run(),
        10 => challenge10::run(),
        11 => challenge11::run(),
        12 => challenge12::run(),
        13 => challenge13::run(),
        14 => challenge14::run(),
        15 => challenge15::run(),
        16 => challenge16::run(),
        17 => challenge17::run(),
        18 => challenge18::run(),
        19 => challenge19::run(),
        20 => challenge20::run(),
        21 => challenge21::run(),
        22 => challenge22::run(),
        23 => challenge23::run(),
        24 => challenge24::run(),
        25 => challenge25::run(),
        26 => challenge26::run(),
        27 => challenge27::run(),
        28 => challenge28::run(),
        29 => challenge29::run(),
        30 => challenge30::run(),
        31 => challenge31::run(),
        32 => challenge32::run(),
        33 => challenge33::run(),
        34 => challenge34::run(),
        35 => challenge35::run(),
        36 => challenge36::run(),
        37 => challenge37::run(),
        38 => challenge38::run(),
        39 => challenge39::run(),
        40 => challenge40::run(),
        41 => challenge41::run(),
        42 => challenge42::run(),
        43 => challenge43::run(),
        44 => challenge44::run(),
        45 => challenge45::run(),
        46 => challenge46::run(),
        47 => challenge47::run(),
        48 => challenge48::run(),
        49 => challenge49::run(),
        50 => challenge50::run(),
        51 => challenge51::run(),
        52 => challenge52::run(),
        53 => challenge53::run(),
        54 => challenge54::run(),
        55 => challenge55::run(),
        56 => challenge56::run(),
        57 => challenge57::run(),
        58 => challenge58::run(),
        59 => challenge59::run(),
        60 => challenge60::run(),
        61 => challenge61::run(),
        62 => challenge62::run(),
        63 => challenge63::run(),
        64 => challenge64::run(),
        65 => challenge65::run(),
        66 => challenge66::run(),
        67 => challenge67::run(),
        68 => challenge68::run(),
        69 => challenge69::run(),
        70 => challenge70::run(),
        71 => challenge71::run(),
        72 => challenge72::run(),
        73 => challenge73::run(),
        74 => challenge74::run(),
        _ => return None,
    })
}

/// Run every challenge in order.  Used by `cargo run -- cryptopals all`.
pub fn run_all() -> Vec<Report> {
    (1..=74).map(|n| run(n).unwrap()).collect()
}
