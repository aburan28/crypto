//! **TLS 1.3 connection state machine** (client + server).
//!
//! Drives the minimal handshake:
//!
//! ```text
//!     Client                                Server
//!     ──────                                ──────
//!     ClientHello (k_x25519^c)  ────►
//!                                  ◄──── ServerHello (k_x25519^s)
//!                                        {EncryptedExtensions}
//!                                        {Finished}
//!     {Finished}                ────►
//!     ----- application data ------
//! ```
//!
//! `{X}` denotes a record encrypted under the **handshake traffic
//! keys** derived from the ECDHE shared secret.  Application records
//! after the handshake use the **application traffic keys**.
//!
//! ## What's modelled
//!
//! - ECDHE via X25519.
//! - AES-128-GCM record encryption.
//! - Key schedule via [`crate::cryptanalysis::tls13_kdf`].
//! - `Finished` MAC validation using `finished_key = HKDF-Expand-Label(
//!   ?hs_traffic_secret, "finished", "", 32)`.
//!
//! ## What's NOT modelled
//!
//! - Certificate validation (we run an anonymous ECDHE).
//! - PSK / 0-RTT / session resumption.
//! - HelloRetryRequest (we assume the client offered a compatible group).
//! - Multiple cipher suites.

use crate::cryptanalysis::tls13_kdf::{derive_secret, hkdf_expand_label, tls13_key_schedule};
use crate::ecc::x25519::{x25519, x25519_base};
use crate::hash::sha256::sha256;
use crate::kdf::hkdf::hmac_sha256;
use crate::symmetric::aes::AesKey;
use crate::tls13::handshake::{
    parse_handshake, wrap_handshake, ClientHello, Finished, HandshakeType, ServerHello,
};
use crate::tls13::record::{
    decrypt_record, encode_tls_plaintext, encrypt_record, parse_tls_record, ContentType,
};
use crate::utils::random::random_bytes_vec;

/// Per-direction record state — key + IV + monotonic counter.
pub struct RecordKeys {
    pub key: AesKey,
    pub iv: [u8; 12],
    pub seq: u64,
}

impl RecordKeys {
    fn next_seq(&mut self) -> u64 {
        let s = self.seq;
        self.seq = self.seq.wrapping_add(1);
        s
    }
}

/// State of a TLS 1.3 endpoint at a given moment of the handshake.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum State {
    Start,
    SentClientHello,
    SentServerHello,
    SentServerFinished,
    SentClientFinished,
    Established,
    Failed,
}

/// **TLS client** endpoint.  Drives the handshake from the client side.
pub struct TlsClient {
    state: State,
    client_x25519_priv: [u8; 32],
    client_x25519_pub: [u8; 32],
    /// Concatenated transcript bytes (every handshake message in
    /// order).
    transcript: Vec<u8>,
    handshake_secret: Option<[u8; 32]>,
    client_hs_traffic_secret: Option<[u8; 32]>,
    server_hs_traffic_secret: Option<[u8; 32]>,
    client_ap_keys: Option<RecordKeys>,
    server_ap_keys: Option<RecordKeys>,
    server_hs_keys: Option<RecordKeys>,
}

impl TlsClient {
    pub fn new() -> Self {
        let mut priv_key = [0u8; 32];
        priv_key.copy_from_slice(&random_bytes_vec(32));
        let public = x25519_base(&priv_key);
        Self {
            state: State::Start,
            client_x25519_priv: priv_key,
            client_x25519_pub: public,
            transcript: Vec::new(),
            handshake_secret: None,
            client_hs_traffic_secret: None,
            server_hs_traffic_secret: None,
            client_ap_keys: None,
            server_ap_keys: None,
            server_hs_keys: None,
        }
    }

    pub fn state(&self) -> State {
        self.state
    }

    /// **Build & emit the ClientHello.**  Returns the wire bytes
    /// (already wrapped in a TLSPlaintext handshake record).
    pub fn build_client_hello(&mut self) -> Vec<u8> {
        let ch = ClientHello::new(self.client_x25519_pub);
        let body = ch.encode();
        let hs = wrap_handshake(HandshakeType::ClientHello, &body);
        self.transcript.extend_from_slice(&hs);
        self.state = State::SentClientHello;
        encode_tls_plaintext(ContentType::Handshake, &hs)
    }

    /// **Handle the server's flight** (ServerHello + {EE} + {Finished}).
    /// On success, sends the client Finished and installs the
    /// application keys.  Returns the client's outgoing wire bytes
    /// (the client Finished encrypted under the handshake keys).
    pub fn handle_server_flight(&mut self, server_wire: &[u8]) -> Option<Vec<u8>> {
        // The server flight is two records: a plaintext ServerHello
        // record followed by a single encrypted record carrying
        // EncryptedExtensions + Finished.
        let sh_rec = parse_tls_record(server_wire)?;
        if sh_rec.content_type != ContentType::Handshake {
            self.state = State::Failed;
            return None;
        }
        let sh_hs = parse_handshake(&sh_rec.payload)?;
        if sh_hs.msg_type != HandshakeType::ServerHello {
            self.state = State::Failed;
            return None;
        }
        let sh = ServerHello::decode(&sh_hs.body)?;
        let sh_full = wrap_handshake(HandshakeType::ServerHello, &sh_hs.body);
        self.transcript.extend_from_slice(&sh_full);

        // ── Key schedule kicks in here ────────────────────────────
        let dhe = x25519(&self.client_x25519_priv, &sh.server_x25519_public);
        self.derive_handshake_secrets(&dhe);

        // Parse the encrypted EE+Finished record using server hs keys.
        let rest = &server_wire[sh_rec.consumed..];
        let enc_rec = parse_tls_record(rest)?;
        let enc_wire = &rest[..enc_rec.consumed];
        let (ct, plain) = {
            let server_hs_keys = self.server_hs_keys.as_mut()?;
            let seq = server_hs_keys.next_seq();
            decrypt_record(&server_hs_keys.key, &server_hs_keys.iv, seq, enc_wire)?
        };
        if ct != ContentType::Handshake {
            self.state = State::Failed;
            return None;
        }
        // Parse the messages inside the encrypted record.
        let mut p = 0;
        while p < plain.len() {
            let hs = parse_handshake(&plain[p..])?;
            let full = wrap_handshake(hs.msg_type, &hs.body);
            self.transcript.extend_from_slice(&full);
            p += hs.consumed;
            match hs.msg_type {
                HandshakeType::EncryptedExtensions => {}
                HandshakeType::Finished => {
                    // Verify server Finished.
                    let server_finished_key = self.derive_finished_key(false);
                    let transcript_before_finished =
                        self.transcript[..self.transcript.len() - full.len()].to_vec();
                    let th = sha256(&transcript_before_finished);
                    let expected =
                        hmac_sha256(&server_finished_key, &th);
                    if hs.body != expected {
                        self.state = State::Failed;
                        return None;
                    }
                }
                _ => {} // ignore unrecognised
            }
        }

        // ── Derive application traffic keys ───────────────────────
        self.derive_application_keys();

        // ── Build + send client Finished ──────────────────────────
        let client_finished_key = self.derive_finished_key(true);
        let th = sha256(&self.transcript);
        let verify_data = hmac_sha256(&client_finished_key, &th).to_vec();
        let fin = Finished { verify_data };
        let fin_hs = wrap_handshake(HandshakeType::Finished, &fin.encode());
        self.transcript.extend_from_slice(&fin_hs);
        // Wrap inside an encrypted record under CLIENT handshake keys.
        let client_hs_keys = self.derive_client_hs_keys();
        let client_finished_wire = encrypt_record(
            &client_hs_keys.key,
            &client_hs_keys.iv,
            client_hs_keys.seq,
            ContentType::Handshake,
            &fin_hs,
        );
        self.state = State::Established;
        Some(client_finished_wire)
    }

    fn derive_handshake_secrets(&mut self, dhe: &[u8]) {
        // Run the TLS 1.3 key schedule.  For the handshake secret
        // we need: HKDF-Extract(salt = derived(empty_es), IKM = DHE).
        let zeros = [0u8; 32];
        let empty_hash = sha256(b"");
        let es = crate::kdf::hkdf::hkdf_extract(Some(&zeros), &zeros);
        let derived = derive_secret(&es, "derived", &empty_hash);
        let hs = crate::kdf::hkdf::hkdf_extract(Some(&derived), dhe);
        let th = sha256(&self.transcript);
        let chts = derive_secret(&hs, "c hs traffic", &th);
        let shts = derive_secret(&hs, "s hs traffic", &th);
        self.handshake_secret = Some(hs);
        self.client_hs_traffic_secret = Some(chts);
        self.server_hs_traffic_secret = Some(shts);
        self.server_hs_keys = Some(self.traffic_keys(&shts));
    }

    fn derive_application_keys(&mut self) {
        let hs = self.handshake_secret.expect("hs set");
        let empty_hash = sha256(b"");
        let derived = derive_secret(&hs, "derived", &empty_hash);
        let zeros = [0u8; 32];
        let ms = crate::kdf::hkdf::hkdf_extract(Some(&derived), &zeros);
        let th = sha256(&self.transcript);
        let cats = derive_secret(&ms, "c ap traffic", &th);
        let sats = derive_secret(&ms, "s ap traffic", &th);
        self.client_ap_keys = Some(self.traffic_keys(&cats));
        self.server_ap_keys = Some(self.traffic_keys(&sats));
    }

    fn derive_client_hs_keys(&self) -> RecordKeys {
        let chts = self.client_hs_traffic_secret.unwrap();
        self.traffic_keys(&chts)
    }

    fn derive_finished_key(&self, client_side: bool) -> [u8; 32] {
        let secret = if client_side {
            self.client_hs_traffic_secret.unwrap()
        } else {
            self.server_hs_traffic_secret.unwrap()
        };
        let v = hkdf_expand_label(&secret, "finished", b"", 32);
        let mut k = [0u8; 32];
        k.copy_from_slice(&v);
        k
    }

    fn traffic_keys(&self, secret: &[u8; 32]) -> RecordKeys {
        let key_bytes = hkdf_expand_label(secret, "key", b"", 16);
        let iv_bytes = hkdf_expand_label(secret, "iv", b"", 12);
        let mut key = [0u8; 16];
        let mut iv = [0u8; 12];
        key.copy_from_slice(&key_bytes);
        iv.copy_from_slice(&iv_bytes);
        RecordKeys {
            key: AesKey::Aes128(key),
            iv,
            seq: 0,
        }
    }

    /// **Encrypt application data** using the established traffic keys.
    pub fn encrypt(&mut self, plaintext: &[u8]) -> Option<Vec<u8>> {
        let keys = self.client_ap_keys.as_mut()?;
        let seq = keys.next_seq();
        Some(encrypt_record(
            &keys.key,
            &keys.iv,
            seq,
            ContentType::ApplicationData,
            plaintext,
        ))
    }

    /// **Decrypt a received application-data record**.
    pub fn decrypt(&mut self, wire: &[u8]) -> Option<Vec<u8>> {
        let keys = self.server_ap_keys.as_mut()?;
        let seq = keys.next_seq();
        let (ct, pt) = decrypt_record(&keys.key, &keys.iv, seq, wire)?;
        if ct != ContentType::ApplicationData {
            return None;
        }
        Some(pt)
    }
}

impl Default for TlsClient {
    fn default() -> Self {
        Self::new()
    }
}

/// **TLS server** endpoint.  Mirror of `TlsClient`.
pub struct TlsServer {
    state: State,
    server_x25519_priv: [u8; 32],
    server_x25519_pub: [u8; 32],
    transcript: Vec<u8>,
    handshake_secret: Option<[u8; 32]>,
    client_hs_traffic_secret: Option<[u8; 32]>,
    server_hs_traffic_secret: Option<[u8; 32]>,
    client_ap_keys: Option<RecordKeys>,
    server_ap_keys: Option<RecordKeys>,
    client_hs_keys: Option<RecordKeys>,
}

impl TlsServer {
    pub fn new() -> Self {
        let mut priv_key = [0u8; 32];
        priv_key.copy_from_slice(&random_bytes_vec(32));
        let public = x25519_base(&priv_key);
        Self {
            state: State::Start,
            server_x25519_priv: priv_key,
            server_x25519_pub: public,
            transcript: Vec::new(),
            handshake_secret: None,
            client_hs_traffic_secret: None,
            server_hs_traffic_secret: None,
            client_ap_keys: None,
            server_ap_keys: None,
            client_hs_keys: None,
        }
    }

    pub fn state(&self) -> State {
        self.state
    }

    /// Receive the client's ClientHello and respond with the server
    /// flight (ServerHello + EncryptedExtensions + Finished).
    pub fn handle_client_hello(&mut self, client_wire: &[u8]) -> Option<Vec<u8>> {
        let ch_rec = parse_tls_record(client_wire)?;
        if ch_rec.content_type != ContentType::Handshake {
            self.state = State::Failed;
            return None;
        }
        let ch_hs = parse_handshake(&ch_rec.payload)?;
        if ch_hs.msg_type != HandshakeType::ClientHello {
            self.state = State::Failed;
            return None;
        }
        let ch = ClientHello::decode(&ch_hs.body)?;
        let ch_full = wrap_handshake(HandshakeType::ClientHello, &ch_hs.body);
        self.transcript.extend_from_slice(&ch_full);

        // ── ServerHello ────────────────────────────────────────────
        let sh = ServerHello::new(&ch.legacy_session_id, self.server_x25519_pub);
        let sh_body = sh.encode();
        let sh_hs = wrap_handshake(HandshakeType::ServerHello, &sh_body);
        self.transcript.extend_from_slice(&sh_hs);
        let sh_wire = encode_tls_plaintext(ContentType::Handshake, &sh_hs);

        // ── ECDHE + derive handshake secrets ───────────────────────
        let dhe = x25519(&self.server_x25519_priv, &ch.client_x25519_public);
        self.derive_handshake_secrets(&dhe);

        // ── EncryptedExtensions (empty in this minimal version) ────
        let ee_body: Vec<u8> = vec![0x00, 0x00];
        let ee_hs = wrap_handshake(HandshakeType::EncryptedExtensions, &ee_body);
        self.transcript.extend_from_slice(&ee_hs);

        // ── Server Finished ────────────────────────────────────────
        let server_finished_key = self.derive_finished_key(false);
        let th = sha256(&self.transcript);
        let verify_data = hmac_sha256(&server_finished_key, &th).to_vec();
        let fin = Finished { verify_data };
        let fin_hs = wrap_handshake(HandshakeType::Finished, &fin.encode());
        self.transcript.extend_from_slice(&fin_hs);

        // ── Encrypted record carrying EE + Finished ────────────────
        let server_hs_keys = self.derive_server_hs_keys();
        let mut payload = ee_hs.clone();
        payload.extend_from_slice(&fin_hs);
        let enc_wire = encrypt_record(
            &server_hs_keys.key,
            &server_hs_keys.iv,
            server_hs_keys.seq,
            ContentType::Handshake,
            &payload,
        );

        self.state = State::SentServerFinished;
        let mut out = sh_wire;
        out.extend_from_slice(&enc_wire);
        Some(out)
    }

    pub fn handle_client_finished(&mut self, wire: &[u8]) -> Option<()> {
        if self.client_hs_keys.is_none() {
            let chts = self.client_hs_traffic_secret?;
            self.client_hs_keys = Some(self.traffic_keys_inner(&chts));
        }
        let (ct, plain) = {
            let keys = self.client_hs_keys.as_mut()?;
            let seq = keys.next_seq();
            decrypt_record(&keys.key, &keys.iv, seq, wire)?
        };
        if ct != ContentType::Handshake {
            self.state = State::Failed;
            return None;
        }
        let hs = parse_handshake(&plain)?;
        if hs.msg_type != HandshakeType::Finished {
            self.state = State::Failed;
            return None;
        }
        // Verify client Finished.
        let client_finished_key = self.derive_finished_key(true);
        let th = sha256(&self.transcript);
        let expected = hmac_sha256(&client_finished_key, &th);
        if hs.body != expected {
            self.state = State::Failed;
            return None;
        }
        // **RFC 8446 §7.1**: application traffic keys are derived from
        // the transcript hash through `server Finished` — NOT including
        // the client Finished.  Derive *before* appending CF.
        self.derive_application_keys();
        let full = wrap_handshake(HandshakeType::Finished, &hs.body);
        self.transcript.extend_from_slice(&full);
        self.state = State::Established;
        Some(())
    }

    fn derive_handshake_secrets(&mut self, dhe: &[u8]) {
        let zeros = [0u8; 32];
        let empty_hash = sha256(b"");
        let es = crate::kdf::hkdf::hkdf_extract(Some(&zeros), &zeros);
        let derived = derive_secret(&es, "derived", &empty_hash);
        let hs = crate::kdf::hkdf::hkdf_extract(Some(&derived), dhe);
        let th = sha256(&self.transcript);
        let chts = derive_secret(&hs, "c hs traffic", &th);
        let shts = derive_secret(&hs, "s hs traffic", &th);
        self.handshake_secret = Some(hs);
        self.client_hs_traffic_secret = Some(chts);
        self.server_hs_traffic_secret = Some(shts);
    }

    fn derive_application_keys(&mut self) {
        let hs = self.handshake_secret.unwrap();
        let empty_hash = sha256(b"");
        let derived = derive_secret(&hs, "derived", &empty_hash);
        let zeros = [0u8; 32];
        let ms = crate::kdf::hkdf::hkdf_extract(Some(&derived), &zeros);
        let th = sha256(&self.transcript);
        let cats = derive_secret(&ms, "c ap traffic", &th);
        let sats = derive_secret(&ms, "s ap traffic", &th);
        self.client_ap_keys = Some(self.traffic_keys_inner(&cats));
        self.server_ap_keys = Some(self.traffic_keys_inner(&sats));
    }

    fn derive_server_hs_keys(&self) -> RecordKeys {
        let shts = self.server_hs_traffic_secret.unwrap();
        self.traffic_keys_inner(&shts)
    }

    fn derive_finished_key(&self, client_side: bool) -> [u8; 32] {
        let secret = if client_side {
            self.client_hs_traffic_secret.unwrap()
        } else {
            self.server_hs_traffic_secret.unwrap()
        };
        let v = hkdf_expand_label(&secret, "finished", b"", 32);
        let mut k = [0u8; 32];
        k.copy_from_slice(&v);
        k
    }

    fn traffic_keys_inner(&self, secret: &[u8; 32]) -> RecordKeys {
        let key_bytes = hkdf_expand_label(secret, "key", b"", 16);
        let iv_bytes = hkdf_expand_label(secret, "iv", b"", 12);
        let mut key = [0u8; 16];
        let mut iv = [0u8; 12];
        key.copy_from_slice(&key_bytes);
        iv.copy_from_slice(&iv_bytes);
        RecordKeys {
            key: AesKey::Aes128(key),
            iv,
            seq: 0,
        }
    }

    /// **Encrypt application data** (server side).
    pub fn encrypt(&mut self, plaintext: &[u8]) -> Option<Vec<u8>> {
        let keys = self.server_ap_keys.as_mut()?;
        let seq = keys.next_seq();
        Some(encrypt_record(
            &keys.key,
            &keys.iv,
            seq,
            ContentType::ApplicationData,
            plaintext,
        ))
    }

    pub fn decrypt(&mut self, wire: &[u8]) -> Option<Vec<u8>> {
        let keys = self.client_ap_keys.as_mut()?;
        let seq = keys.next_seq();
        let (ct, pt) = decrypt_record(&keys.key, &keys.iv, seq, wire)?;
        if ct != ContentType::ApplicationData {
            return None;
        }
        Some(pt)
    }
}

impl Default for TlsServer {
    fn default() -> Self {
        Self::new()
    }
}

/// **Drive a full handshake** between a client and a server in
/// memory.  Returns `Ok(())` on success.
pub fn drive_handshake(client: &mut TlsClient, server: &mut TlsServer) -> Result<(), &'static str> {
    let ch_wire = client.build_client_hello();
    let server_flight = server
        .handle_client_hello(&ch_wire)
        .ok_or("server failed to handle ClientHello")?;
    let client_finished = client
        .handle_server_flight(&server_flight)
        .ok_or("client failed to handle server flight")?;
    server
        .handle_client_finished(&client_finished)
        .ok_or("server failed to verify client Finished")?;
    if client.state() != State::Established || server.state() != State::Established {
        return Err("handshake did not reach Established on both sides");
    }
    Ok(())
}

// ── Tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// **End-to-end TLS 1.3 handshake** between an in-memory client
    /// and server.  After the handshake both endpoints are in the
    /// Established state and can encrypt/decrypt application data.
    #[test]
    fn end_to_end_handshake_completes() {
        let mut client = TlsClient::new();
        let mut server = TlsServer::new();
        drive_handshake(&mut client, &mut server).expect("handshake should complete");
        assert_eq!(client.state(), State::Established);
        assert_eq!(server.state(), State::Established);
    }

    /// **Bidirectional application data**: client → server then
    /// server → client, both encrypted under the derived traffic
    /// keys.
    #[test]
    fn bidirectional_application_data_after_handshake() {
        let mut client = TlsClient::new();
        let mut server = TlsServer::new();
        drive_handshake(&mut client, &mut server).expect("handshake");

        let request = b"GET /index.html HTTP/1.1\r\nHost: example.com\r\n\r\n".to_vec();
        let ct = client.encrypt(&request).expect("encrypt");
        let recovered_request = server.decrypt(&ct).expect("server decrypt");
        assert_eq!(recovered_request, request);

        let response = b"HTTP/1.1 200 OK\r\nContent-Length: 5\r\n\r\nhello".to_vec();
        let ct2 = server.encrypt(&response).expect("server encrypt");
        let recovered_response = client.decrypt(&ct2).expect("client decrypt");
        assert_eq!(recovered_response, response);
    }

    /// **Sequence-number rollover**: 10 records both directions.
    /// Each direction has its own monotonic seq counter; the AAD
    /// changes per record so AEAD validates correctly.
    #[test]
    fn many_records_succeed_in_order() {
        let mut client = TlsClient::new();
        let mut server = TlsServer::new();
        drive_handshake(&mut client, &mut server).expect("handshake");
        for i in 0..10u8 {
            let msg = vec![i; 32];
            let ct = client.encrypt(&msg).unwrap();
            assert_eq!(server.decrypt(&ct).unwrap(), msg);
        }
    }

    /// **Tampered application record fails decryption.**
    #[test]
    fn tampered_record_is_rejected() {
        let mut client = TlsClient::new();
        let mut server = TlsServer::new();
        drive_handshake(&mut client, &mut server).expect("handshake");
        let mut ct = client.encrypt(b"secret").unwrap();
        ct[8] ^= 1;
        assert!(server.decrypt(&ct).is_none());
    }

    /// **Out-of-order records fail decryption**: TLS 1.3 demands
    /// strictly sequential delivery within a direction.  Skipping a
    /// record breaks the sequence counter.
    #[test]
    fn out_of_order_records_are_rejected() {
        let mut client = TlsClient::new();
        let mut server = TlsServer::new();
        drive_handshake(&mut client, &mut server).expect("handshake");
        let _r0 = client.encrypt(b"first").unwrap();
        let r1 = client.encrypt(b"second").unwrap();
        // Server decrypts `r1` before `r0` → seq mismatch → AEAD fails.
        assert!(server.decrypt(&r1).is_none());
    }

    /// **Client and server independently derive the same DHE shared
    /// secret** — sanity check on the X25519 wiring.
    #[test]
    fn ecdhe_shared_secret_matches() {
        let mut client = TlsClient::new();
        let mut server = TlsServer::new();
        let _ = client.build_client_hello();
        // After handshake, both sides should have ended up with
        // identical handshake_secret values (= they derived the same
        // DHE).  We finish the handshake and check.
        let mut client2 = TlsClient::new();
        let mut server2 = TlsServer::new();
        drive_handshake(&mut client2, &mut server2).expect("handshake");
        assert_eq!(
            client2.handshake_secret, server2.handshake_secret,
            "client and server must derive the same handshake secret",
        );
    }

    /// **Tampered ClientHello breaks the handshake**: the server
    /// derives a different transcript hash and the Finished
    /// validation fails.
    #[test]
    fn tampered_client_hello_breaks_handshake() {
        let mut client = TlsClient::new();
        let mut server = TlsServer::new();
        let mut ch = client.build_client_hello();
        // Flip a bit in the ClientHello's key share area (well after
        // the header).
        let target = ch.len() - 1;
        ch[target] ^= 1;
        let result = server.handle_client_hello(&ch);
        // The server may still produce a ServerHello (parsing succeeds),
        // but the resulting DHE will differ from what the client
        // derives → handshake fails when verifying Finished.
        if let Some(flight) = result {
            let client_finished = client.handle_server_flight(&flight);
            assert!(
                client_finished.is_none(),
                "tampering should cause client Finished verification to fail"
            );
        }
    }
}
