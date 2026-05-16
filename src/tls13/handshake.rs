//! **TLS 1.3 handshake-message encoding/decoding** (RFC 8446 §4).
//!
//! Every handshake message has a 4-byte header:
//!
//! ```text
//!     struct {
//!         HandshakeType msg_type;
//!         uint24 length;
//!         opaque body[length];
//!     } Handshake;
//! ```
//!
//! We support a minimal subset sufficient for an ECDHE-X25519 +
//! AES-128-GCM-only handshake:
//!
//! - `ClientHello` / `ServerHello`
//! - `EncryptedExtensions` (mostly empty in this minimal version)
//! - `Finished`
//!
//! Certificate / CertificateVerify are out of scope for this v1 —
//! we model an anonymous ECDHE handshake (= no server-side
//! authentication beyond proving knowledge of the key share).  For a
//! real implementation, see RustCrypto's `rustls`.

use crate::utils::random::random_bytes_vec;

/// Handshake message types we recognise.
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HandshakeType {
    ClientHello = 1,
    ServerHello = 2,
    EncryptedExtensions = 8,
    Certificate = 11,
    CertificateVerify = 15,
    Finished = 20,
}

impl HandshakeType {
    pub fn from_u8(b: u8) -> Option<Self> {
        match b {
            1 => Some(HandshakeType::ClientHello),
            2 => Some(HandshakeType::ServerHello),
            8 => Some(HandshakeType::EncryptedExtensions),
            11 => Some(HandshakeType::Certificate),
            15 => Some(HandshakeType::CertificateVerify),
            20 => Some(HandshakeType::Finished),
            _ => None,
        }
    }
}

/// Wrap a handshake-message body with the 4-byte header.
pub fn wrap_handshake(msg_type: HandshakeType, body: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(4 + body.len());
    out.push(msg_type as u8);
    out.extend_from_slice(&[
        ((body.len() >> 16) & 0xFF) as u8,
        ((body.len() >> 8) & 0xFF) as u8,
        (body.len() & 0xFF) as u8,
    ]);
    out.extend_from_slice(body);
    out
}

/// Parsed handshake header.
#[derive(Clone, Debug)]
pub struct ParsedHandshake {
    pub msg_type: HandshakeType,
    pub body: Vec<u8>,
    pub consumed: usize,
}

/// Parse a single handshake message from a buffer.
pub fn parse_handshake(buf: &[u8]) -> Option<ParsedHandshake> {
    if buf.len() < 4 {
        return None;
    }
    let msg_type = HandshakeType::from_u8(buf[0])?;
    let len =
        ((buf[1] as usize) << 16) | ((buf[2] as usize) << 8) | (buf[3] as usize);
    if buf.len() < 4 + len {
        return None;
    }
    Some(ParsedHandshake {
        msg_type,
        body: buf[4..4 + len].to_vec(),
        consumed: 4 + len,
    })
}

// ── ClientHello ──────────────────────────────────────────────────────

/// Minimal TLS 1.3 ClientHello.
///
/// ```text
///     struct {
///         ProtocolVersion legacy_version = 0x0303;
///         Random random;
///         opaque legacy_session_id<0..32>;
///         CipherSuite cipher_suites<2..2^16-2>;
///         opaque legacy_compression_methods<1..2^8-1>;
///         Extension extensions<8..2^16-1>;
///     } ClientHello;
/// ```
///
/// We hard-code:
/// - `cipher_suites = [TLS_AES_128_GCM_SHA256 (0x1301)]`
/// - `legacy_compression_methods = [0]` (null)
/// - extensions: `supported_versions = [TLS 1.3]`, `supported_groups = [x25519]`,
///   `key_share` with one X25519 entry, `signature_algorithms = [ed25519]`.
#[derive(Clone, Debug)]
pub struct ClientHello {
    pub random: [u8; 32],
    pub legacy_session_id: Vec<u8>,
    pub client_x25519_public: [u8; 32],
}

impl ClientHello {
    /// Generate a fresh ClientHello with random nonce + session id.
    pub fn new(client_x25519_public: [u8; 32]) -> Self {
        let mut random = [0u8; 32];
        random.copy_from_slice(&random_bytes_vec(32));
        let legacy_session_id = random_bytes_vec(32);
        Self {
            random,
            legacy_session_id,
            client_x25519_public,
        }
    }

    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::new();
        // legacy_version = 0x0303
        out.extend_from_slice(&[0x03, 0x03]);
        // random
        out.extend_from_slice(&self.random);
        // legacy_session_id
        out.push(self.legacy_session_id.len() as u8);
        out.extend_from_slice(&self.legacy_session_id);
        // cipher_suites: just TLS_AES_128_GCM_SHA256
        out.extend_from_slice(&[0x00, 0x02, 0x13, 0x01]);
        // legacy_compression_methods: null
        out.extend_from_slice(&[0x01, 0x00]);
        // Extensions
        let exts = encode_client_extensions(&self.client_x25519_public);
        out.extend_from_slice(&(exts.len() as u16).to_be_bytes());
        out.extend_from_slice(&exts);
        out
    }

    pub fn decode(body: &[u8]) -> Option<Self> {
        let mut p = 0;
        if body.len() < 2 + 32 + 1 {
            return None;
        }
        // skip legacy_version
        p += 2;
        let mut random = [0u8; 32];
        random.copy_from_slice(&body[p..p + 32]);
        p += 32;
        // session_id
        let sid_len = *body.get(p)? as usize;
        p += 1;
        if body.len() < p + sid_len {
            return None;
        }
        let legacy_session_id = body[p..p + sid_len].to_vec();
        p += sid_len;
        // cipher_suites
        if body.len() < p + 2 {
            return None;
        }
        let cs_len = u16::from_be_bytes([body[p], body[p + 1]]) as usize;
        p += 2 + cs_len;
        // compression_methods
        if body.len() < p + 1 {
            return None;
        }
        let cm_len = body[p] as usize;
        p += 1 + cm_len;
        // Extensions
        if body.len() < p + 2 {
            return None;
        }
        let ext_len = u16::from_be_bytes([body[p], body[p + 1]]) as usize;
        p += 2;
        if body.len() < p + ext_len {
            return None;
        }
        let ext_bytes = &body[p..p + ext_len];
        let client_x25519_public = parse_key_share_x25519(ext_bytes)?;
        Some(Self {
            random,
            legacy_session_id,
            client_x25519_public,
        })
    }
}

// ── ServerHello ──────────────────────────────────────────────────────

#[derive(Clone, Debug)]
pub struct ServerHello {
    pub random: [u8; 32],
    pub legacy_session_id_echo: Vec<u8>,
    pub server_x25519_public: [u8; 32],
}

impl ServerHello {
    pub fn new(client_session_id: &[u8], server_x25519_public: [u8; 32]) -> Self {
        let mut random = [0u8; 32];
        random.copy_from_slice(&random_bytes_vec(32));
        Self {
            random,
            legacy_session_id_echo: client_session_id.to_vec(),
            server_x25519_public,
        }
    }

    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::new();
        // legacy_version
        out.extend_from_slice(&[0x03, 0x03]);
        // random
        out.extend_from_slice(&self.random);
        // session_id echo
        out.push(self.legacy_session_id_echo.len() as u8);
        out.extend_from_slice(&self.legacy_session_id_echo);
        // chosen cipher_suite
        out.extend_from_slice(&[0x13, 0x01]);
        // legacy_compression_method
        out.push(0x00);
        // Extensions: supported_versions + key_share
        let exts = encode_server_extensions(&self.server_x25519_public);
        out.extend_from_slice(&(exts.len() as u16).to_be_bytes());
        out.extend_from_slice(&exts);
        out
    }

    pub fn decode(body: &[u8]) -> Option<Self> {
        let mut p = 0;
        p += 2; // legacy_version
        if body.len() < p + 32 {
            return None;
        }
        let mut random = [0u8; 32];
        random.copy_from_slice(&body[p..p + 32]);
        p += 32;
        let sid_len = *body.get(p)? as usize;
        p += 1;
        let legacy_session_id_echo = body[p..p + sid_len].to_vec();
        p += sid_len;
        p += 2; // cipher_suite
        p += 1; // legacy_compression_method
        // Extensions
        if body.len() < p + 2 {
            return None;
        }
        let ext_len = u16::from_be_bytes([body[p], body[p + 1]]) as usize;
        p += 2;
        let ext_bytes = &body[p..p + ext_len];
        let server_x25519_public = parse_key_share_x25519(ext_bytes)?;
        Some(Self {
            random,
            legacy_session_id_echo,
            server_x25519_public,
        })
    }
}

// ── Extension encoding ───────────────────────────────────────────────

const EXT_SUPPORTED_VERSIONS: u16 = 43;
const EXT_SUPPORTED_GROUPS: u16 = 10;
const EXT_KEY_SHARE: u16 = 51;
const EXT_SIGNATURE_ALGORITHMS: u16 = 13;
const NAMED_GROUP_X25519: u16 = 29;
const SIG_SCHEME_ED25519: u16 = 0x0807;

fn encode_extension(ext_type: u16, body: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(4 + body.len());
    out.extend_from_slice(&ext_type.to_be_bytes());
    out.extend_from_slice(&(body.len() as u16).to_be_bytes());
    out.extend_from_slice(body);
    out
}

fn encode_client_extensions(x25519_pub: &[u8; 32]) -> Vec<u8> {
    let mut out = Vec::new();
    // supported_versions = [TLS 1.3 = 0x0304]
    out.extend_from_slice(&encode_extension(
        EXT_SUPPORTED_VERSIONS,
        &[0x02, 0x03, 0x04],
    ));
    // supported_groups = [x25519]
    let mut sg = Vec::new();
    sg.extend_from_slice(&(2u16).to_be_bytes());
    sg.extend_from_slice(&NAMED_GROUP_X25519.to_be_bytes());
    out.extend_from_slice(&encode_extension(EXT_SUPPORTED_GROUPS, &sg));
    // signature_algorithms = [ed25519]
    let mut sa = Vec::new();
    sa.extend_from_slice(&(2u16).to_be_bytes());
    sa.extend_from_slice(&SIG_SCHEME_ED25519.to_be_bytes());
    out.extend_from_slice(&encode_extension(EXT_SIGNATURE_ALGORITHMS, &sa));
    // key_share: one X25519 entry
    let mut ks = Vec::new();
    // ClientShares vector
    let mut shares = Vec::new();
    shares.extend_from_slice(&NAMED_GROUP_X25519.to_be_bytes());
    shares.extend_from_slice(&(32u16).to_be_bytes());
    shares.extend_from_slice(x25519_pub);
    ks.extend_from_slice(&(shares.len() as u16).to_be_bytes());
    ks.extend_from_slice(&shares);
    out.extend_from_slice(&encode_extension(EXT_KEY_SHARE, &ks));
    out
}

fn encode_server_extensions(x25519_pub: &[u8; 32]) -> Vec<u8> {
    let mut out = Vec::new();
    // supported_versions = TLS 1.3 (selected)
    out.extend_from_slice(&encode_extension(EXT_SUPPORTED_VERSIONS, &[0x03, 0x04]));
    // key_share: server's single chosen share
    let mut ks = Vec::new();
    ks.extend_from_slice(&NAMED_GROUP_X25519.to_be_bytes());
    ks.extend_from_slice(&(32u16).to_be_bytes());
    ks.extend_from_slice(x25519_pub);
    out.extend_from_slice(&encode_extension(EXT_KEY_SHARE, &ks));
    out
}

/// Scan an extensions blob for `key_share` and return the X25519
/// public point.
fn parse_key_share_x25519(ext_bytes: &[u8]) -> Option<[u8; 32]> {
    let mut p = 0;
    while p + 4 <= ext_bytes.len() {
        let ext_type = u16::from_be_bytes([ext_bytes[p], ext_bytes[p + 1]]);
        let ext_len =
            u16::from_be_bytes([ext_bytes[p + 2], ext_bytes[p + 3]]) as usize;
        if p + 4 + ext_len > ext_bytes.len() {
            return None;
        }
        let body = &ext_bytes[p + 4..p + 4 + ext_len];
        if ext_type == EXT_KEY_SHARE {
            return parse_key_share_body(body);
        }
        p += 4 + ext_len;
    }
    None
}

fn parse_key_share_body(body: &[u8]) -> Option<[u8; 32]> {
    // Body format may be ClientHello (length-prefixed list of shares)
    // or ServerHello (single share, no outer length).  We try both.
    if body.len() >= 36 && body[0] == 0 && body[1] == NAMED_GROUP_X25519 as u8 {
        // Looks like a server-style single share.
        let kx_len = u16::from_be_bytes([body[2], body[3]]) as usize;
        if kx_len == 32 && body.len() >= 4 + 32 {
            let mut k = [0u8; 32];
            k.copy_from_slice(&body[4..36]);
            return Some(k);
        }
    }
    if body.len() < 2 {
        return None;
    }
    // Try client-style: length-prefixed shares list.
    let total = u16::from_be_bytes([body[0], body[1]]) as usize;
    if body.len() < 2 + total {
        return None;
    }
    let mut q = 2;
    while q + 4 <= 2 + total {
        let group = u16::from_be_bytes([body[q], body[q + 1]]);
        let kx_len =
            u16::from_be_bytes([body[q + 2], body[q + 3]]) as usize;
        if q + 4 + kx_len > 2 + total {
            return None;
        }
        if group == NAMED_GROUP_X25519 && kx_len == 32 {
            let mut k = [0u8; 32];
            k.copy_from_slice(&body[q + 4..q + 4 + 32]);
            return Some(k);
        }
        q += 4 + kx_len;
    }
    None
}

// ── Finished ─────────────────────────────────────────────────────────

/// `Finished` body = `verify_data` (HMAC of transcript hash with
/// the appropriate `_finished_key`).  See RFC 8446 §4.4.4.
#[derive(Clone, Debug)]
pub struct Finished {
    pub verify_data: Vec<u8>,
}

impl Finished {
    pub fn encode(&self) -> Vec<u8> {
        self.verify_data.clone()
    }
    pub fn decode(body: &[u8]) -> Self {
        Self {
            verify_data: body.to_vec(),
        }
    }
}

// ── Tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn handshake_header_round_trips() {
        let body = b"abcd".to_vec();
        let wire = wrap_handshake(HandshakeType::Finished, &body);
        assert_eq!(wire[0], HandshakeType::Finished as u8);
        let parsed = parse_handshake(&wire).expect("parse");
        assert_eq!(parsed.msg_type, HandshakeType::Finished);
        assert_eq!(parsed.body, body);
        assert_eq!(parsed.consumed, wire.len());
    }

    #[test]
    fn client_hello_round_trip() {
        let ch = ClientHello::new([0x42u8; 32]);
        let body = ch.encode();
        let ch2 = ClientHello::decode(&body).expect("decode");
        assert_eq!(ch2.random, ch.random);
        assert_eq!(ch2.legacy_session_id, ch.legacy_session_id);
        assert_eq!(ch2.client_x25519_public, ch.client_x25519_public);
    }

    #[test]
    fn server_hello_round_trip() {
        let sh = ServerHello::new(b"session-id", [0x88u8; 32]);
        let body = sh.encode();
        let sh2 = ServerHello::decode(&body).expect("decode");
        assert_eq!(sh2.random, sh.random);
        assert_eq!(sh2.legacy_session_id_echo, sh.legacy_session_id_echo);
        assert_eq!(sh2.server_x25519_public, sh.server_x25519_public);
    }

    #[test]
    fn finished_round_trip() {
        let f = Finished { verify_data: vec![1, 2, 3, 4, 5, 6, 7, 8] };
        let body = f.encode();
        let f2 = Finished::decode(&body);
        assert_eq!(f2.verify_data, f.verify_data);
    }
}
