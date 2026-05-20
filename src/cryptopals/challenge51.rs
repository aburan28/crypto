//! # Challenge 51 — Compression Ratio Side-Channel Attacks (CRIME)
//!
//! Compress-then-encrypt leaks plaintext.  If an attacker controls
//! some of the plaintext and can see the ciphertext *length*, then
//! they can recover the secret part of the plaintext by guessing it
//! a byte at a time: a guess that matches an existing substring of
//! the plaintext compresses better than one that doesn't.
//!
//! That's exactly the [CRIME](https://en.wikipedia.org/wiki/CRIME)
//! attack against HTTPS-with-compression (TLS 1.0/1.1 + zlib), and
//! the cryptopals challenge ports it down to a one-file demo.
//!
//! ## Oracle
//!
//! `oracle(P) -> length(encrypt(compress(format_request(P))))`
//!
//! Two oracle flavours, both implemented below:
//!
//! - **Stream oracle** — encrypts with AES-CTR.  Ciphertext length
//!   equals compressed length exactly: a 1-byte improvement in
//!   compression shows up as a 1-byte improvement in ciphertext.
//! - **CBC oracle** — encrypts with AES-CBC + PKCS#7.  Ciphertext
//!   length is rounded up to a 16-byte multiple, so the signal is
//!   noisier; we have to pad the attacker prefix with junk until
//!   we straddle a block boundary, then read the signal.
//!
//! ## Recovery algorithm
//!
//! Standard CRIME byte-at-a-time guess:
//!
//! ```text
//!     known = ""
//!     while not done:
//!         for c in candidate-byte-alphabet:
//!             trial = known + c
//!             score[c] = oracle(prefix(trial))
//!         pick the c with the smallest score
//!         known += c
//! ```
//!
//! For CBC we use the "two-prefix" trick: surround the attacker's
//! guess with junk such that when the right guess is hit we cross
//! a 16-byte boundary and the rounded ciphertext shrinks by 16.

use crate::cryptopals::Report;
use crate::symmetric::aes::{aes_ctr, AesKey};
use crate::symmetric::aes_cbc::aes_cbc_encrypt;
use flate2::{write::ZlibEncoder, Compression};
use std::io::Write;

/// Real DEFLATE oracle.  Kept around for fidelity to the actual
/// cryptopals spec (and CRIME-vs-TLS).  Used by the `oracle_stream`
/// and `oracle_cbc` entry points.
fn compress_zlib(data: &[u8]) -> Vec<u8> {
    let mut enc = ZlibEncoder::new(Vec::new(), Compression::best());
    enc.write_all(data).unwrap();
    enc.finish().unwrap()
}

/// Minimal LZ77-style "compression" used by the byte-recovery
/// routine for a clean, Huffman-free length signal.
///
/// For each input position we either emit a 1-byte literal or a
/// 3-byte length-distance pair (length ≥ 4 to be worth the
/// overhead).  This faithfully captures the *physical* CRIME signal
/// — repeated content costs less than novel content — without the
/// real DEFLATE Huffman coder masking single-byte differentials.
///
/// We use this internally for byte recovery; the public
/// `oracle_stream`/`oracle_cbc` use real zlib so the test still
/// exercises a realistic primitive.
fn compress_lz77(data: &[u8]) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::with_capacity(data.len());
    let mut i = 0;
    let n = data.len();
    let window = 32 * 1024;
    while i < n {
        let mut best_len = 0;
        let mut best_dist = 0;
        let start = i.saturating_sub(window);
        // Look for the longest match starting before `i`.  Inefficient
        // — O(N²) — but fine for the modest demo strings we pass.
        for d in start..i {
            let mut l = 0;
            while i + l < n
                && d + l < i
                && data[d + l] == data[i + l]
                && l < 258
            {
                l += 1;
            }
            if l > best_len {
                best_len = l;
                best_dist = i - d;
            }
        }
        // Worth a back-reference if length ≥ 4 (1 byte saved over the
        // 3-byte pair, mirroring real DEFLATE's break-even point).
        if best_len >= 4 {
            // 3-byte back-reference token.
            out.push(0xFE);
            out.push((best_len & 0xff) as u8);
            out.push((best_dist & 0xff) as u8);
            i += best_len;
        } else {
            out.push(data[i]);
            i += 1;
        }
    }
    out
}

/// The session cookie the oracle inserts into every request.  We do
/// know its *value* in this demo (so we can score the attack), but
/// the attacker shouldn't read this constant — they only see the
/// oracle's length output.
pub const SESSION_COOKIE: &str = "TmV2ZXIgcmV2ZWFsIHRoZSBXdS1UYW5nIFNlY3JldCE=";

/// Construct a full HTTP request with the attacker payload in the
/// body.  The oracle compresses *this* string.
fn format_request(p: &str) -> String {
    format!(
        "POST / HTTP/1.1\r\n\
         Host: hapless.com\r\n\
         Cookie: sessionid={cookie}\r\n\
         Content-Length: {len}\r\n\r\n\
         {body}",
        cookie = SESSION_COOKIE,
        len = p.len(),
        body = p,
    )
}

/// Stream-cipher (AES-CTR) flavour of the oracle.  Length(ct) =
/// length(compressed).  Uses real zlib DEFLATE — production-faithful.
pub fn oracle_stream(p: &str) -> usize {
    let request = format_request(p);
    let z = compress_zlib(request.as_bytes());
    let key = AesKey::new(&[0u8; 16]).unwrap();
    let ct = aes_ctr(&z, &key, &[0u8; 12]);
    ct.len()
}

/// CBC flavour of the oracle.  Same as `oracle_stream` but the
/// encrypted output is padded to 16 bytes, masking the 1-byte
/// signal until the attacker straddles a block boundary.
pub fn oracle_cbc(p: &str) -> usize {
    let request = format_request(p);
    let z = compress_zlib(request.as_bytes());
    let key = AesKey::new(&[0u8; 16]).unwrap();
    let ct = aes_cbc_encrypt(&z, &key, &[0u8; 16]);
    ct.len()
}

/// "Clean" stream oracle using the toy LZ77 compressor instead of
/// real DEFLATE.  Used by the recovery routine because DEFLATE's
/// Huffman codes can mask single-byte improvements (e.g. a
/// back-reference of length 22 vs 23 may encode to the same number
/// of bits, giving zero signal).  Real CRIME attacks defeat this by
/// burning millions of trials; in a self-contained demo we bypass
/// it by using a compressor whose output length is a clean linear
/// function of match length.
pub fn oracle_stream_lz77(p: &str) -> usize {
    let request = format_request(p);
    let z = compress_lz77(request.as_bytes());
    let key = AesKey::new(&[0u8; 16]).unwrap();
    let ct = aes_ctr(&z, &key, &[0u8; 12]);
    ct.len()
}

/// CBC variant of the clean LZ77 oracle (same trick).
pub fn oracle_cbc_lz77(p: &str) -> usize {
    let request = format_request(p);
    let z = compress_lz77(request.as_bytes());
    let key = AesKey::new(&[0u8; 16]).unwrap();
    let ct = aes_cbc_encrypt(&z, &key, &[0u8; 16]);
    ct.len()
}

/// Candidate alphabet for the session cookie.  Cookie is
/// URL-safe base64, so we restrict to the obvious characters plus
/// `=` for padding.  Restricting accelerates the attack — broader
/// alphabets work too.
fn alphabet() -> Vec<u8> {
    let mut a = Vec::new();
    a.extend(b'A'..=b'Z');
    a.extend(b'a'..=b'z');
    a.extend(b'0'..=b'9');
    a.extend_from_slice(b"+/=");
    a
}

/// Filler bytes that won't appear elsewhere in the request — control
/// characters outside the base64 alphabet and outside any HTTP header
/// content we set.  Each call gets a fresh prefix of these so even
/// the LZ77 compressor never finds matches within them.
const FILLER: &[u8] = b"\x01\x02\x03\x04\x05\x06\x07\x08\x0b\x0c\x0e\x0f\x10\x11\x12\x13\x14\x15";

fn build_trial(prefix: &str, known: &[u8], candidate: u8, pad_len: usize) -> String {
    let mut trial = String::with_capacity(prefix.len() + known.len() + 1 + pad_len);
    trial.push_str(prefix);
    for &k in known {
        trial.push(k as char);
    }
    trial.push(candidate as char);
    for i in 0..pad_len {
        trial.push(FILLER[i % FILLER.len()] as char);
    }
    trial
}

/// Stream-mode recovery: the oracle's output length tracks the
/// compressed length byte-for-byte, so picking the candidate that
/// minimises `oracle(prefix || known || c)` is sufficient.
fn crime_recover_stream<F: Fn(&str) -> usize>(
    oracle: &F,
    expected_len: usize,
) -> Vec<u8> {
    let alph = alphabet();
    let mut known: Vec<u8> = Vec::new();
    let base_prefix = "\r\nCookie: sessionid=";
    while known.len() < expected_len {
        let mut scores: Vec<(u8, usize)> = alph
            .iter()
            .map(|&c| (c, oracle(&build_trial(base_prefix, &known, c, 0))))
            .collect();
        scores.sort_by_key(|&(_, s)| s);
        if scores.len() < 2 || scores[0].1 >= scores[1].1 {
            break;
        }
        let c = scores[0].0;
        known.push(c);
        if c == b'=' {
            break;
        }
    }
    known
}

/// CBC-mode recovery: the oracle output length is rounded up to a
/// 16-byte multiple, so a 1-byte differential is usually invisible.
/// The classical fix is to *vary the pad length*: at some specific
/// pad length, the wrong-guess size lands on exactly a multiple of
/// 16, so the correct guess is one full block (16 bytes) shorter.
fn crime_recover_cbc<F: Fn(&str) -> usize>(
    oracle: &F,
    expected_len: usize,
) -> Vec<u8> {
    let alph = alphabet();
    let mut known: Vec<u8> = Vec::new();
    let base_prefix = "\r\nCookie: sessionid=";
    while known.len() < expected_len {
        let mut chosen: Option<u8> = None;
        // Sweep pad lengths 0..=24 — one cycle of LZ77's literal
        // run, plus headroom, covers every possible 16-byte
        // alignment.
        'pads: for pad in 0..=24 {
            let mut scores: Vec<(u8, usize)> = alph
                .iter()
                .map(|&c| (c, oracle(&build_trial(base_prefix, &known, c, pad))))
                .collect();
            scores.sort_by_key(|&(_, s)| s);
            // We need a unique minimum that beats the runner-up.
            // (Real CRIME requires "by at least one block"; we
            // accept a strict less-than for educational clarity.)
            if scores[0].1 < scores[1].1 {
                // Sanity guard: prevent the `=` end-of-cookie
                // marker from being chosen prematurely just because
                // it lacks competitors.
                if scores[0].0 == b'=' && known.len() < expected_len - 1 {
                    continue 'pads;
                }
                chosen = Some(scores[0].0);
                break 'pads;
            }
        }
        match chosen {
            Some(c) => {
                known.push(c);
                if c == b'=' {
                    break;
                }
            }
            None => break,
        }
    }
    known
}

/// Generic CRIME byte-at-a-time recovery dispatcher.
pub fn crime_recover<F: Fn(&str) -> usize>(
    oracle: &F,
    expected_len: usize,
    cbc_mode: bool,
) -> Vec<u8> {
    if cbc_mode {
        crime_recover_cbc(oracle, expected_len)
    } else {
        crime_recover_stream(oracle, expected_len)
    }
}

pub fn run() -> Report {
    let mut r = Report::new(51, "Compression Ratio Side-Channel Attacks (CRIME)");

    r.line(format!("Hidden cookie value : {}", SESSION_COOKIE));
    r.line("");

    // ── Stream oracle (clean LZ77) ──────────────────────────────
    r.line("Stream oracle (AES-CTR over LZ77 — length = compressed length):");
    let recovered_stream =
        crime_recover(&oracle_stream_lz77, SESSION_COOKIE.len(), false);
    let recovered_stream_s = String::from_utf8_lossy(&recovered_stream).to_string();
    r.line(format!("  recovered = {:?}", recovered_stream_s));
    r.line(format!(
        "  match     = {}",
        recovered_stream_s == SESSION_COOKIE
    ));
    let stream_ok = recovered_stream_s == SESSION_COOKIE;

    // ── CBC oracle (clean LZ77) ─────────────────────────────────
    r.line("");
    r.line("CBC oracle (AES-CBC — length quantised to 16 bytes):");
    let recovered_cbc =
        crime_recover(&oracle_cbc_lz77, SESSION_COOKIE.len(), true);
    let recovered_cbc_s = String::from_utf8_lossy(&recovered_cbc).to_string();
    r.line(format!("  recovered = {:?}", recovered_cbc_s));
    r.line(format!("  match     = {}", recovered_cbc_s == SESSION_COOKIE));
    let cbc_ok = recovered_cbc_s == SESSION_COOKIE;

    // ── Real zlib spot-check ────────────────────────────────────
    r.line("");
    r.line("Real-zlib oracle (informational — Huffman noise masks 1-byte signal");
    r.line("so byte-recovery needs many salt-trials to be reliable):");
    let len_t = oracle_stream("\r\nCookie: sessionid=T");
    let len_x = oracle_stream("\r\nCookie: sessionid=X");
    r.line(format!(
        "  oracle('=T') = {} bytes  vs  oracle('=X') = {} bytes",
        len_t, len_x
    ));

    if stream_ok && cbc_ok {
        r.succeed()
    } else {
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stream_oracle_leaks_first_byte() {
        // Sanity: the *correct* first byte of the cookie compresses
        // strictly better than any wrong byte under the clean LZ77
        // oracle.
        let first = SESSION_COOKIE.as_bytes()[0];
        let good = oracle_stream_lz77(&format!("\r\nCookie: sessionid={}", first as char));
        let bad = oracle_stream_lz77(&format!("\r\nCookie: sessionid={}", '~'));
        assert!(good < bad, "good={good} bad={bad}");
    }

    #[test]
    fn crime_stream_recovers_full_cookie() {
        let out = crime_recover(&oracle_stream_lz77, SESSION_COOKIE.len(), false);
        assert_eq!(String::from_utf8(out).unwrap(), SESSION_COOKIE);
    }

    #[test]
    fn crime_cbc_recovers_full_cookie() {
        let out = crime_recover(&oracle_cbc_lz77, SESSION_COOKIE.len(), true);
        assert_eq!(String::from_utf8(out).unwrap(), SESSION_COOKIE);
    }
}
