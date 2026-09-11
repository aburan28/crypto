//! **TLS 1.3 protocol handler** (RFC 8446).
//!
//! Minimal-but-working TLS 1.3 implementation:
//!
//! - **ECDHE** via X25519 (`named_group = 0x001D`).
//! - **AEAD** via AES-128-GCM (`cipher_suite = TLS_AES_128_GCM_SHA256`).
//! - **Key schedule** via [`crate::cryptanalysis::tls13_kdf`].
//! - **Record layer** with TLS 1.3 per-record AAD + nonce derivation
//!   (RFC 8446 §5.2-5.3).
//! - **Handshake state machine** for both client and server, ending
//!   in a verified `Finished` exchange.
//!
//! ## What's modelled
//!
//! - 1-RTT handshake with `ClientHello` → `ServerHello` →
//!   `{EncryptedExtensions}` → `{Finished}` → `{Finished}`.
//! - Application-data records encrypted under the application
//!   traffic keys.
//! - Sequence-number-bound AEAD (out-of-order delivery rejected).
//!
//! ## What's NOT modelled
//!
//! - Certificate / CertificateVerify (anonymous ECDHE only — see
//!   `connection.rs` docstring for the rationale).
//! - PSK, 0-RTT early data, session resumption.
//! - HelloRetryRequest (assumes client and server agree on x25519).
//! - Multiple cipher suites or named groups.
//! - Alerts (failures just put the state machine in `Failed`).
//!
//! ## Example
//!
//! ```
//! use crypto::tls13::connection::{drive_handshake, TlsClient, TlsServer};
//! let mut client = TlsClient::new();
//! let mut server = TlsServer::new();
//! drive_handshake(&mut client, &mut server).unwrap();
//! let ciphertext = client.encrypt(b"hello, server").unwrap();
//! let plaintext = server.decrypt(&ciphertext).unwrap();
//! assert_eq!(plaintext, b"hello, server");
//! ```

pub mod connection;
pub mod handshake;
pub mod record;

pub use connection::{drive_handshake, RecordKeys, State, TlsClient, TlsServer};
pub use handshake::{ClientHello, Finished, HandshakeType, ServerHello};
pub use record::{ContentType, decrypt_record, encrypt_record, parse_tls_record};
