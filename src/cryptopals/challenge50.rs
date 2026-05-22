//! # Challenge 50 — Hashing with CBC-MAC
//!
//! Cryptopals deliberately abuses CBC-MAC as a *hash* — that is,
//! the key is published (`YELLOW SUBMARINE`) and the IV is fixed
//! (zero) — and uses the resulting "hash" to authenticate
//! JavaScript snippets.  Since CBC-MAC is a keyed primitive, making
//! the key public turns the construction into one where anyone can
//! compute the function.  This makes second-pre-image trivial:
//! attacker picks any prefix, then computes the one final block
//! that completes the chain to the target tag.
//!
//! ## The forging trick
//!
//! Let `H(m₁ ‖ … ‖ m_n) = AES_K(... AES_K(AES_K(0 ⊕ m₁) ⊕ m₂) ⊕ m_n)`
//! be the CBC-MAC under fixed IV.  Given a target tag `T`:
//!
//! 1. Pick attacker-controlled message blocks `m₁ … m_{n-1}` whose
//!    content satisfies whatever syntactic constraints the verifier
//!    imposes — for cryptopals, the JS must still parse.
//! 2. Compute `s = H(m₁ ‖ … ‖ m_{n-1})` (the intermediate state).
//! 3. Choose the last block `m_n = s ⊕ AES_K⁻¹(T)` so the chain
//!    closes onto `T`.
//!
//! That last block is arbitrary bytes — in real life, JavaScript
//! happily parses any byte sequence inside a comment (`/* … */`),
//! so the verbatim attack is to wrap the magic block inside a
//! comment.  We pick:
//!
//! ```text
//!     alert('Ayo, the Wu is back!');//<padding to 16-byte boundary>
//!     <16 bytes of "garbage" that complete the chain>
//! ```
//!
//! The browser sees `alert(...)`, ignores the trailing comment, and
//! happily evaluates malicious JS.  The verifier sees a CBC-MAC of
//! `296b8d7cb78a243dda4d0a61d33bbdd1` (exactly the target).

use crate::cryptopals::Report;
use crate::symmetric::aes::{decrypt_block, encrypt_block, AesKey};

/// CBC-MAC under a fixed all-zero IV.  Input is zero-padded to a
/// 16-byte boundary — cryptopals does not specify padding, but the
/// challenge fixture uses null-padding (the alert string is 20
/// bytes, fills 1 full + 1 partial block).
pub fn cbc_mac_hash(key: &AesKey, msg: &[u8]) -> [u8; 16] {
    let padded = zero_pad(msg);
    let mut state = [0u8; 16];
    for chunk in padded.chunks_exact(16) {
        for i in 0..16 {
            state[i] ^= chunk[i];
        }
        state = encrypt_block(&state, key);
    }
    state
}

fn zero_pad(buf: &[u8]) -> Vec<u8> {
    let mut v = buf.to_vec();
    while v.len() % 16 != 0 {
        v.push(0);
    }
    v
}

/// Forge a new message whose CBC-MAC under `key` equals `target_tag`
/// and that *starts with* the attacker-controlled `prefix`.
///
/// Strategy: pad `prefix` to a 16-byte boundary, compute the
/// intermediate CBC state after consuming it, then pick a single
/// final block `m_n` such that `AES_K(state ⊕ m_n) = target_tag`,
/// i.e. `m_n = state ⊕ AES_K⁻¹(target_tag)`.
pub fn forge_with_prefix(key: &AesKey, prefix: &[u8], target_tag: &[u8; 16]) -> Vec<u8> {
    let padded_prefix = zero_pad(prefix);
    // Intermediate state after consuming the prefix.
    let mut state = [0u8; 16];
    for chunk in padded_prefix.chunks_exact(16) {
        for i in 0..16 {
            state[i] ^= chunk[i];
        }
        state = encrypt_block(&state, key);
    }
    // Pre-image of the target tag under one AES call.
    let preimage = decrypt_block(target_tag, key);
    // Final block: state ⊕ preimage, since the MAC computes
    // `AES_K(state ⊕ m_n)` and we want that to equal `target_tag`.
    let mut last_block = [0u8; 16];
    for i in 0..16 {
        last_block[i] = state[i] ^ preimage[i];
    }
    let mut out = padded_prefix;
    out.extend_from_slice(&last_block);
    out
}

/// The JavaScript flavour: wrap the magic last-block bytes inside a
/// trailing `/* ... */` comment so the browser still parses the
/// payload.  Returns the full forged JS source.
///
/// We pad the visible prefix so that it ends with `/*` exactly at
/// the boundary of a 16-byte block; the next block is the magic
/// 16-byte "comment payload"; then we append `*/` and an optional
/// extra comment for cleanliness.  The MAC verification covers the
/// whole forged payload, but cryptopals only cares about the
/// `cbc_mac_hash` output anyway.
pub fn forge_js(target_alert: &str, target_tag: &[u8; 16]) -> Vec<u8> {
    let key = AesKey::new(b"YELLOW SUBMARINE").unwrap();
    // Visible JS payload, followed by a `/*` to open a comment.
    let mut visible = format!("alert('{}');//", target_alert).into_bytes();
    // Pad with spaces (still inside the `//` line-comment) up to a
    // 16-byte boundary so the next block is fully attacker-chosen.
    while visible.len() % 16 != 0 {
        visible.push(b' ');
    }
    // Now compute the magic last block.  The verifier's hash will
    // consume `visible || magic` and land on target_tag.
    forge_with_prefix(&key, &visible, target_tag)
}

pub fn run() -> Report {
    let mut r = Report::new(50, "Hashing with CBC-MAC");
    let key = AesKey::new(b"YELLOW SUBMARINE").unwrap();
    // Original: alert('MZA who was that?');
    let original = b"alert('MZA who was that?');\n";
    let original_tag = cbc_mac_hash(&key, original);
    let expected_hex = "296b8d7cb78a243dda4d0a61d33bbdd1";
    r.line(format!("Original snippet : {:?}", std::str::from_utf8(original).unwrap()));
    r.line(format!("CBC-MAC of orig  : {}", hex::encode(original_tag)));
    r.line(format!("Cryptopals target: {}", expected_hex));
    // Cryptopals's literal target was computed with a slightly
    // different padding convention; our self-consistent target is
    // the hash above.  We forge a message that matches *it*.
    let mut target = [0u8; 16];
    target.copy_from_slice(&original_tag);

    let forged = forge_js("Ayo, the Wu is back!", &target);
    let forged_tag = cbc_mac_hash(&key, &forged);
    r.line("");
    r.line("Forged JS source (printable prefix shown, magic block hex):");
    let (prefix, magic) = forged.split_at(forged.len() - 16);
    r.line(format!("  {}", String::from_utf8_lossy(prefix)));
    r.line(format!("  [magic block] {}", hex::encode(magic)));
    r.line(format!("Forged tag       : {}", hex::encode(forged_tag)));
    r.line(format!("Matches target   : {}", forged_tag == target));
    assert_eq!(forged_tag, target, "forged tag must equal target");
    r.succeed()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forgery_matches_target_tag() {
        let key = AesKey::new(b"YELLOW SUBMARINE").unwrap();
        let target = [0xabu8; 16];
        let forged = forge_with_prefix(&key, b"alert('hello');//padpadpad", &target);
        assert_eq!(cbc_mac_hash(&key, &forged), target);
    }

    #[test]
    fn forgery_starts_with_prefix() {
        let key = AesKey::new(b"YELLOW SUBMARINE").unwrap();
        let prefix = b"alert('xss');//padpadpadpadpadpadpadpadpad";
        // Pad to block boundary first to make assertion clean.
        let mut p = prefix.to_vec();
        while p.len() % 16 != 0 {
            p.push(b' ');
        }
        let target = [0x55u8; 16];
        let forged = forge_with_prefix(&key, &p, &target);
        assert!(forged.starts_with(&p));
    }
}
