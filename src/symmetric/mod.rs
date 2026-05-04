//! Symmetric-key cryptography: AES (block + CTR + GCM) and ChaCha20-Poly1305.

pub mod aes;
pub mod aes_cbc;
pub mod chacha20;
pub mod serpent;
pub mod threefish;

pub use aes::{
    AesKey, aes_ctr, aes_gcm_decrypt, aes_gcm_encrypt,
    decrypt_block, encrypt_block,
};
pub use aes_cbc::{
    aes_cbc_decrypt, aes_cbc_encrypt, padding_oracle_attack,
    pkcs7_pad, pkcs7_unpad,
};
pub use chacha20::{
    chacha20_block, chacha20_poly1305_decrypt, chacha20_poly1305_encrypt, chacha20_xor,
    poly1305,
};
