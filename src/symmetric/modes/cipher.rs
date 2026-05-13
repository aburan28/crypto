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

// ── Tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symmetric::aes::AesKey;

    /// AES round-trip via the BlockCipher trait.
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

    /// SM4 round-trip via the BlockCipher trait.
    #[test]
    fn sm4_blockcipher_round_trip() {
        let cipher = Sm4::new(&[0u8; 16]);
        let mut block = [0x01u8; 16];
        let original = block;
        BlockCipher::<16>::encrypt_block(&cipher, &mut block);
        assert_ne!(block, original);
        BlockCipher::<16>::decrypt_block(&cipher, &mut block);
        assert_eq!(block, original);
    }

    /// Generic helper accepting `impl BlockCipher<16>` works across
    /// cipher types — proves the trait actually unifies them.
    fn round_trip<C: BlockCipher<16>>(cipher: &C, mut block: [u8; 16]) -> bool {
        let original = block;
        cipher.encrypt_block(&mut block);
        let mid_differs = block != original;
        cipher.decrypt_block(&mut block);
        mid_differs && block == original
    }

    #[test]
    fn generic_round_trip_works_for_aes_sm4_kuznyechik() {
        let aes = AesKey::new(&[1u8; 32]).unwrap();
        let sm4 = Sm4::new(&[1u8; 16]);
        let kuz = Kuznyechik::new(&[1u8; 32]);
        assert!(round_trip(&aes, [2u8; 16]));
        assert!(round_trip(&sm4, [2u8; 16]));
        assert!(round_trip(&kuz, [2u8; 16]));
    }

    /// Magma is the only 8-byte-block cipher we ship; verify it
    /// implements `BlockCipher<8>`.
    #[test]
    fn magma_blockcipher_round_trip() {
        let magma = Magma::new(&[0x42u8; 32]);
        let mut block = [0xAAu8; 8];
        let original = block;
        BlockCipher::<8>::encrypt_block(&magma, &mut block);
        assert_ne!(block, original);
        BlockCipher::<8>::decrypt_block(&magma, &mut block);
        assert_eq!(block, original);
    }
}
