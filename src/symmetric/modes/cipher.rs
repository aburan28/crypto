//! **Block-cipher abstraction** for generic modes of operation.
//!
//! Each cipher implements [`BlockCipher`] (or [`BlockCipher128`] for
//! the AES-class 16-byte ciphers).  Modes (ECB, CFB, OFB, CCM, …) are
//! generic over this trait so the same code wraps AES, SM4, Serpent,
//! Kuznyechik, etc.
//!
//! ## Why a trait instead of free functions
//!
//! Modes share a lot of code (CTR chaining, GHASH, CMAC, …).  Pushing
//! the cipher behind a trait means each mode is written once, against
//! `impl BlockCipher<N>`, and the user picks the underlying cipher at
//! call site.  At toy-library scale we accept the indirection cost;
//! production crates monomorphise via generics, which the compiler
//! handles for us.

/// Generic block cipher with a compile-time-known block size of `N`
/// bytes.  Implementors expose `encrypt_block` and `decrypt_block` on
/// fixed-size byte arrays; everything else (ECB, CFB, CTR, …) is
/// derived in [`super`].
pub trait BlockCipher<const N: usize> {
    /// Encrypt one block in place.
    fn encrypt_block(&self, block: &mut [u8; N]);
    /// Decrypt one block in place.
    fn decrypt_block(&self, block: &mut [u8; N]);
}

/// Convenience marker for 16-byte (128-bit) block ciphers — the
/// AES-family ciphers (AES, SM4, Serpent, Kuznyechik).
///
/// CCM, SIV, GCM, GCM-SIV, OCB3, XTS, KW, and EAX all require a
/// 16-byte block, so we collect them under this alias.
pub trait BlockCipher128: BlockCipher<16> {}
impl<T: BlockCipher<16>> BlockCipher128 for T {}

// ── BlockCipher impls for our shipping primitives ────────────────────
//
// NOTE: many cipher modules referenced by the original trait-impl
// file (Aria, Camellia, Twofish, Speck, Simon, …) are placeholders in
// this worktree.  Only the ciphers that actually have real
// implementations get a `BlockCipher` impl here.  Restore the rest
// once those modules land.

use crate::symmetric::aes::{decrypt_block, encrypt_block, AesKey};
use crate::symmetric::des::Des;
use crate::symmetric::des3::TripleDes;
use crate::symmetric::gost_magma::Magma;
use crate::symmetric::kuznyechik::Kuznyechik;
use crate::symmetric::rc5::Rc5;
use crate::symmetric::serpent::{serpent_decrypt, serpent_encrypt, SerpentKey};
use crate::symmetric::sm4::Sm4;

impl BlockCipher<16> for AesKey {
    fn encrypt_block(&self, block: &mut [u8; 16]) {
        *block = encrypt_block(block, self);
    }
    fn decrypt_block(&self, block: &mut [u8; 16]) {
        *block = decrypt_block(block, self);
    }
}

impl BlockCipher<16> for Sm4 {
    fn encrypt_block(&self, block: &mut [u8; 16]) {
        Sm4::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 16]) {
        Sm4::decrypt_block(self, block);
    }
}

impl BlockCipher<16> for Kuznyechik {
    fn encrypt_block(&self, block: &mut [u8; 16]) {
        Kuznyechik::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 16]) {
        Kuznyechik::decrypt_block(self, block);
    }
}

impl BlockCipher<16> for SerpentKey {
    fn encrypt_block(&self, block: &mut [u8; 16]) {
        *block = serpent_encrypt(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 16]) {
        *block = serpent_decrypt(self, block);
    }
}

impl BlockCipher<8> for Magma {
    fn encrypt_block(&self, block: &mut [u8; 8]) {
        Magma::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 8]) {
        Magma::decrypt_block(self, block);
    }
}

impl BlockCipher<8> for Des {
    fn encrypt_block(&self, block: &mut [u8; 8]) {
        Des::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 8]) {
        Des::decrypt_block(self, block);
    }
}

impl BlockCipher<8> for TripleDes {
    fn encrypt_block(&self, block: &mut [u8; 8]) {
        TripleDes::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 8]) {
        TripleDes::decrypt_block(self, block);
    }
}

impl BlockCipher<8> for Rc5 {
    fn encrypt_block(&self, block: &mut [u8; 8]) {
        Rc5::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 8]) {
        Rc5::decrypt_block(self, block);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aes_blockcipher_round_trip() {
        let key = AesKey::new(&[0u8; 16]).unwrap();
        let mut block = *b"YELLOW SUBMARINE";
        let original = block;
        BlockCipher::<16>::encrypt_block(&key, &mut block);
        assert_ne!(block, original);
        BlockCipher::<16>::decrypt_block(&key, &mut block);
        assert_eq!(block, original);
    }
}
