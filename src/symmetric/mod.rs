//! Symmetric-key cryptography.
//!
//! ## Block ciphers
//!
//! - **AES-128 / AES-256** — `aes::AesKey` + `encrypt_block` /
//!   `decrypt_block` (FIPS 197).
//! - **Serpent** (`serpent`), **SM4** (`sm4`), **Magma** (64-bit
//!   GOST `gost_magma`), **Kuznyechik** (128-bit GOST `kuznyechik`),
//!   **Threefish** (`threefish`).
//!
//! ## Stream / authenticated modes
//!
//! - **AES-CTR / GCM** — see [`aes`].
//! - **AES-CBC** — see [`aes_cbc`] (includes a deliberate padding-oracle
//!   demo).
//! - **AES-CMAC** — see [`cmac`].
//! - **ChaCha20-Poly1305 AEAD** — see [`chacha20`].
//!
//! ## Generic modes of operation
//!
//! [`modes`] hosts generic block-cipher modes that work over the
//! [`modes::BlockCipher<N>`] trait — so the same code wraps AES, SM4,
//! Serpent, Kuznyechik, etc.  Includes:
//!
//! - Confidentiality-only: [`modes::ecb`], [`modes::cfb`], [`modes::ofb`].
//! - AEAD: [`modes::ccm`] (RFC 3610), [`modes::siv`] (RFC 5297,
//!   misuse-resistant).
//! - Specialised: [`modes::kw`] (RFC 3394 key wrap),
//!   [`modes::xts`] (IEEE P1619 disk encryption).

pub mod aes;
pub mod aes_cbc;
pub mod chacha20;
pub mod cmac;
pub mod des;
pub mod des3;
pub mod gost_magma;
pub mod kuznyechik;
pub mod modes;
pub mod rc4;
pub mod rc5;
pub mod serpent;
pub mod sm4;
pub mod threefish;
pub mod visualize;

pub use aes::{aes_ctr, aes_gcm_decrypt, aes_gcm_encrypt, decrypt_block, encrypt_block, AesKey};
pub use aes_cbc::{
    aes_cbc_decrypt, aes_cbc_encrypt, padding_oracle_attack, pkcs7_pad as cbc_pkcs7_pad,
    pkcs7_unpad as cbc_pkcs7_unpad,
};
pub use chacha20::{
    chacha20_block, chacha20_poly1305_decrypt, chacha20_poly1305_encrypt, chacha20_xor, poly1305,
};
pub use cmac::{aes_cmac, CmacTag};
pub use des::Des;
pub use des3::TripleDes;
pub use gost_magma::Magma;
pub use kuznyechik::Kuznyechik;
pub use modes::{
    ccm_decrypt, ccm_encrypt, cfb_decrypt, cfb_encrypt, ecb_decrypt, ecb_encrypt, kw_unwrap,
    kw_wrap, ofb_apply, siv_decrypt, siv_encrypt, xts_decrypt, xts_encrypt, BlockCipher,
    BlockCipher128, CcmTagLen,
};
pub use rc4::{rc4, Rc4};
pub use rc5::Rc5;
pub use sm4::Sm4;
