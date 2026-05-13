//! Encoding / serialisation primitives — DER, PEM, ASN.1.
//!
//! Bitcoin / TLS / X.509 / PKCS-anything all use DER (Distinguished
//! Encoding Rules).  Without DER parsing the library cannot ingest
//! real-world signatures or certificates; the
//! [`crate::cryptanalysis::signature_corpus`] analyzer needs this
//! to consume Bitcoin transaction data without an external parser.

pub mod der;

pub use der::{decode_ecdsa_signature, encode_ecdsa_signature, DerError};
