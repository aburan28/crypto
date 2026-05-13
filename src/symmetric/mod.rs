//! Symmetric-key cryptography: AES (block + CTR + GCM) and ChaCha20-Poly1305.

pub mod aes;
pub mod aes_cbc;
pub mod chacha20;
pub mod cmac;
pub mod gost_magma;
pub mod kuznyechik;
pub mod serpent;
pub mod sm4;
pub mod threefish;

pub use aes::{aes_ctr, aes_gcm_decrypt, aes_gcm_encrypt, decrypt_block, encrypt_block, AesKey};
pub use aes_cbc::{
    aes_cbc_decrypt, aes_cbc_encrypt, padding_oracle_attack, pkcs7_pad, pkcs7_unpad,
};
pub use chacha20::{
    chacha20_block, chacha20_poly1305_decrypt, chacha20_poly1305_encrypt, chacha20_xor, poly1305,
};
pub use cmac::{aes_cmac, CmacTag};
pub use gost_magma::Magma;
pub use kuznyechik::Kuznyechik;
pub use sm4::Sm4;
