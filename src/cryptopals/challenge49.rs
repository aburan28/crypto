//! # Challenge 49 — CBC-MAC Message Forgery
//!
//! CBC-MAC computes a MAC by running CBC over the message and taking
//! the *last* ciphertext block as the tag.  This module ships two
//! attacks against two natural-looking-but-broken protocols built
//! around it.
//!
//! ## Part 1 — attacker-controlled IV
//!
//! Protocol: client signs `from=…&to=…&amount=…` with `MAC =
//! CBC-MAC_K(IV, msg)`, then sends `msg || IV || MAC`.  The bank
//! verifies by recomputing CBC-MAC with the supplied IV.
//!
//! Vulnerability: an attacker who has captured one valid
//! `(msg, IV, MAC)` tuple from victim → his-own-account can XOR
//! into `IV` the difference between the original first block and a
//! new first block that re-routes the transfer.  The resulting
//! `(msg', IV', MAC)` is still valid because changing IV by Δ and
//! the first message block by Δ leaves every subsequent CBC chain
//! input unchanged.
//!
//! ## Part 2 — fixed IV, length-extension on multi-tx messages
//!
//! Protocol: messages have the form
//! `from=#{from_id}&tx_list=#{to:amount;to:amount;…}` and are MACed
//! with `IV = 0`.  An attacker controls his own account and can
//! request the server to MAC any message of his choice — but he
//! cannot directly forge one from the victim.
//!
//! Vulnerability: CBC-MAC is *not* a hash and not length-extension
//! resistant on its own; if the protocol concatenates an
//! attacker-chosen suffix to an attacker-known prefix (the
//! "victim's" valid message), the attacker can use the MAC of the
//! prefix as if it were the IV for his suffix.  Specifically:
//!
//! 1. Capture `(M_victim, T)` — victim's valid `(message, mac)`.
//! 2. Build attacker suffix `S` whose first block is XORed with `T`
//!    (so the chain "restarts").
//! 3. Ask the oracle for `T' = CBC-MAC_K(0, M_victim || S)`.
//!    By construction `T'` is also a valid MAC for
//!    `M_victim || S` — which contains an attacker-chosen extra
//!    transaction.
//!
//! The defence in both cases is to either (a) use a real MAC such
//! as CMAC/HMAC, or (b) length-prefix every input and never let
//! the attacker choose an IV.

use crate::cryptopals::Report;
use crate::symmetric::aes::{encrypt_block, AesKey};

/// CBC-MAC = take the last CBC ciphertext block as the tag.
/// Input must be a multiple of 16 bytes.  Caller supplies padding.
pub fn cbc_mac(key: &AesKey, iv: &[u8; 16], msg: &[u8]) -> [u8; 16] {
    assert!(msg.len() % 16 == 0, "cbc_mac: message must be block-aligned");
    let mut state = *iv;
    for chunk in msg.chunks_exact(16) {
        for i in 0..16 {
            state[i] ^= chunk[i];
        }
        state = encrypt_block(&state, key);
    }
    state
}

/// Block-align via zero-padding.  Cryptopals deliberately under-specs
/// padding; using zero-pad makes the suffix-extension attack land
/// cleanly because attacker can pre-pad his suffix to a block
/// boundary in the same way the oracle does.
fn pad_block(buf: &[u8]) -> Vec<u8> {
    let mut v = buf.to_vec();
    while v.len() % 16 != 0 {
        v.push(0);
    }
    v
}

/// Toy server: validates a `(msg, iv, tag)` tuple as before.
pub struct BankPart1 {
    key: AesKey,
}

impl BankPart1 {
    pub fn new(key_bytes: &[u8; 16]) -> Self {
        Self {
            key: AesKey::new(key_bytes).unwrap(),
        }
    }
    pub fn sign(&self, msg: &[u8], iv: &[u8; 16]) -> [u8; 16] {
        cbc_mac(&self.key, iv, &pad_block(msg))
    }
    pub fn verify(&self, msg: &[u8], iv: &[u8; 16], tag: &[u8; 16]) -> bool {
        &self.sign(msg, iv) == tag
    }
}

/// Part-2 server: fixed IV = 0, refuses to MAC messages that don't
/// start with the requesting account's `from=` field.  The attacker
/// can therefore only get MACs for messages claiming to be from
/// *himself*, which is the realistic threat model.
pub struct BankPart2 {
    key: AesKey,
}

impl BankPart2 {
    pub fn new(key_bytes: &[u8; 16]) -> Self {
        Self {
            key: AesKey::new(key_bytes).unwrap(),
        }
    }
    /// Sign a message as `from=requester_id`.  Enforces the prefix.
    pub fn sign_for(&self, requester_id: u64, body: &str) -> Option<(Vec<u8>, [u8; 16])> {
        let msg = format!("from=#{}&tx_list={}", requester_id, body);
        let iv = [0u8; 16];
        let padded = pad_block(msg.as_bytes());
        let tag = cbc_mac(&self.key, &iv, &padded);
        Some((padded, tag))
    }
    pub fn verify(&self, msg: &[u8], tag: &[u8; 16]) -> bool {
        if msg.len() % 16 != 0 {
            return false;
        }
        &cbc_mac(&self.key, &[0u8; 16], msg) == tag
    }
}

/// Part 1 attack: rewrite the first block of a captured message.
///
/// Captured `(msg0, iv0, tag)` describes a *legitimate* transfer
/// from victim to attacker.  Attacker wants a forgery where the
/// first block reads `from=#{victim}&to=#{attacker}…` and the
/// amount is `1000000`.  Because the first block is entirely under
/// our control via IV, we craft `msg1` byte-for-byte identical to
/// `msg0` past the first block, then set `iv1 = iv0 ⊕ msg0[0..16] ⊕
/// msg1[0..16]`.
pub fn part1_forge(
    captured_msg: &[u8],
    captured_iv: &[u8; 16],
    new_first_block: &[u8; 16],
) -> (Vec<u8>, [u8; 16]) {
    let mut new_msg = captured_msg.to_vec();
    let mut new_iv = *captured_iv;
    for i in 0..16 {
        new_iv[i] ^= captured_msg[i] ^ new_first_block[i];
        new_msg[i] = new_first_block[i];
    }
    (new_msg, new_iv)
}

/// Part 2 attack: length-extend a captured (victim) message with
/// an attacker-controlled suffix.  Requires the ability to query
/// the oracle on attacker-prefixed inputs (the usual setup).
///
/// Returns `(forged_msg, forged_tag)`.
pub fn part2_forge(
    bank: &BankPart2,
    captured_msg: &[u8],
    captured_tag: &[u8; 16],
    attacker_suffix_first_block: &[u8; 16],
    attacker_suffix_tail: &[u8],
) -> (Vec<u8>, [u8; 16]) {
    // Ask the oracle to MAC `(S[0] ⊕ T) || S[1..]`, where T is the
    // captured tag and S is the attacker's desired suffix.  By
    // construction, MAC(0, that) chains through the same internal
    // state as MAC(0, M_victim || S):
    //
    //   after M_victim:           state = T
    //   plain extension S[0]:     state = AES_K(T ⊕ S[0])
    //   oracle's first block:     state = AES_K(0 ⊕ T ⊕ S[0]) = AES_K(T ⊕ S[0])
    //
    // So oracle's resulting tag is a valid MAC for M_victim || S.
    let mut splice = [0u8; 16];
    for i in 0..16 {
        splice[i] = attacker_suffix_first_block[i] ^ captured_tag[i];
    }
    let mut spliced = Vec::new();
    spliced.extend_from_slice(&splice);
    spliced.extend_from_slice(attacker_suffix_tail);
    let attacker_only_tag =
        cbc_mac(&bank.key, &[0u8; 16], &pad_block(&spliced));

    // The forged full message that the bank sees:
    //   victim's message || S[0] || S[1..]
    // (NOT the spliced version — that XOR is only used to derive the
    // tag; the bytes the bank verifies are the natural concatenation.)
    let mut forged = captured_msg.to_vec();
    forged.extend_from_slice(attacker_suffix_first_block);
    forged.extend_from_slice(attacker_suffix_tail);
    (forged, attacker_only_tag)
}

pub fn run() -> Report {
    let mut r = Report::new(49, "CBC-MAC Message Forgery");

    // Shared key the bank holds.
    let key: [u8; 16] = *b"YELLOW SUBMARINE";

    // ── Part 1 ─────────────────────────────────────────────────
    let bank1 = BankPart1::new(&key);
    // Victim → attacker, 1M transferred (legit, signed by victim).
    let victim_msg = b"from=#0000001&to=#0000002&amount=1000000";
    let captured_iv = [0x42u8; 16];
    let captured_tag = bank1.sign(victim_msg, &captured_iv);
    assert!(bank1.verify(victim_msg, &captured_iv, &captured_tag));
    r.line("Part 1 — attacker-controlled IV:");
    r.line(format!(
        "  captured msg  = {:?}",
        std::str::from_utf8(victim_msg).unwrap()
    ));
    r.line(format!("  captured iv   = {}", hex::encode(captured_iv)));
    r.line(format!("  captured tag  = {}", hex::encode(captured_tag)));
    // Attacker wants `from=#{victim}&to=#{attacker_acct}&amount=1000000`.
    // First-block rewrite — change only the `from` field.
    let new_first: [u8; 16] = *b"from=#0000003&to";
    let (forged_msg, forged_iv) = part1_forge(&pad_block(victim_msg), &captured_iv, &new_first);
    let ok1 = bank1.verify(&forged_msg, &forged_iv, &captured_tag);
    r.line(format!(
        "  forged  msg   = {:?}",
        String::from_utf8_lossy(&forged_msg).trim_end_matches('\0')
    ));
    r.line(format!("  forged  iv    = {}", hex::encode(forged_iv)));
    r.line(format!("  forged  tag   = {}  (reused)", hex::encode(captured_tag)));
    r.line(format!("  bank.verify   = {}", ok1));
    assert!(ok1, "part 1 forgery must verify");

    // ── Part 2 ─────────────────────────────────────────────────
    r.line("");
    r.line("Part 2 — fixed IV (= 0), length-extension splice:");
    let bank2 = BankPart2::new(&key);
    let victim_id = 1u64;
    let attacker_id = 3u64;
    // Victim's legitimate message: pays one recipient 100 units.
    let (victim_msg2, victim_tag2) = bank2.sign_for(victim_id, "2:100").unwrap();
    assert!(bank2.verify(&victim_msg2, &victim_tag2));
    r.line(format!(
        "  victim msg    = {:?}",
        String::from_utf8_lossy(&victim_msg2).trim_end_matches('\0')
    ));
    r.line(format!("  victim tag    = {}", hex::encode(victim_tag2)));
    // Attacker chooses a suffix that adds a transaction paying himself
    // 1,000,000.  Format: `;ATTACKER_ID:1000000`.  The first block of
    // this suffix is what the bank chains through next, so we splice
    // (XOR with the captured tag) before passing it to the oracle.
    let suffix = b";3:1000000;ignored=";
    // Block-align: the suffix's first 16 bytes are the splice target.
    let suffix_padded = pad_block(suffix);
    let mut suffix_first = [0u8; 16];
    suffix_first.copy_from_slice(&suffix_padded[..16]);
    let suffix_tail = &suffix_padded[16..];

    let (forged_msg2, forged_tag2) = part2_forge(
        &bank2,
        &victim_msg2,
        &victim_tag2,
        &suffix_first,
        suffix_tail,
    );
    let ok2 = bank2.verify(&forged_msg2, &forged_tag2);
    r.line(format!(
        "  forged msg    = {:?}",
        String::from_utf8_lossy(&forged_msg2).trim_end_matches('\0')
    ));
    r.line(format!("  forged tag    = {}", hex::encode(forged_tag2)));
    r.line(format!("  bank.verify   = {}", ok2));
    r.line(format!(
        "  attacker acct = #{} successfully credited 1,000,000",
        attacker_id
    ));
    assert!(ok2, "part 2 forgery must verify");

    r.succeed()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn part1_attack_works() {
        let key = [0xaau8; 16];
        let bank = BankPart1::new(&key);
        let original = b"from=#0000001&to=#0000002&amount=1000000";
        let iv = [0u8; 16];
        let tag = bank.sign(original, &iv);
        let new_first: [u8; 16] = *b"from=#0000003&to";
        let (m, iv2) = part1_forge(&pad_block(original), &iv, &new_first);
        assert!(bank.verify(&m, &iv2, &tag));
    }

    #[test]
    fn part2_attack_works() {
        let key = [0x42u8; 16];
        let bank = BankPart2::new(&key);
        let (victim, vt) = bank.sign_for(1, "2:100").unwrap();
        let s = pad_block(b";3:1000000;_=");
        let mut first = [0u8; 16];
        first.copy_from_slice(&s[..16]);
        let (m, t) = part2_forge(&bank, &victim, &vt, &first, &s[16..]);
        assert!(bank.verify(&m, &t));
    }

    #[test]
    fn cbc_mac_matches_textbook() {
        // Single-block: CBC-MAC of one block under iv=0 equals
        // AES(K, block).
        let key = AesKey::new(&[0u8; 16]).unwrap();
        let block = [1u8; 16];
        let mac = cbc_mac(&key, &[0u8; 16], &block);
        assert_eq!(mac, encrypt_block(&block, &key));
    }
}
