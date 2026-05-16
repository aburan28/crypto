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
use crate::symmetric::aria::{Aria128, Aria192, Aria256};
use crate::symmetric::blowfish::Blowfish;
use crate::symmetric::camellia::{Camellia128, Camellia192, Camellia256};
use crate::symmetric::cast5::Cast5;
use crate::symmetric::des::Des;
use crate::symmetric::des3::TripleDes;
use crate::symmetric::feal::Feal8;
use crate::symmetric::gift::{Gift128, Gift64};
use crate::symmetric::gost_magma::Magma;
use crate::symmetric::hight::Hight;
use crate::symmetric::idea::Idea;
use crate::symmetric::kasumi::Kasumi;
use crate::symmetric::kuznyechik::Kuznyechik;
use crate::symmetric::lucifer::Lucifer;
use crate::symmetric::mars::Mars;
use crate::symmetric::misty1::Misty1;
use crate::symmetric::present::{Present128, Present80};
use crate::symmetric::rc5::Rc5;
use crate::symmetric::rc6::Rc6;
use crate::symmetric::safer::SaferK64;
use crate::symmetric::seed::Seed;
use crate::symmetric::serpent::{serpent_decrypt, serpent_encrypt, SerpentKey};
use crate::symmetric::simon::{Simon128_128, Simon128_256, Simon64_128};
use crate::symmetric::skinny::{Skinny128_128, Skinny64_128};
use crate::symmetric::skipjack::Skipjack;
use crate::symmetric::sm4::Sm4;
use crate::symmetric::speck::{Speck128_128, Speck128_256, Speck64_128};
use crate::symmetric::square::Square;
use crate::symmetric::tea::Tea;
use crate::symmetric::twofish::Twofish;
use crate::symmetric::xtea::Xtea;

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

impl BlockCipher<16> for Speck128_128 {
    fn encrypt_block(&self, block: &mut [u8; 16]) {
        Speck128_128::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 16]) {
        Speck128_128::decrypt_block(self, block);
    }
}

impl BlockCipher<16> for Speck128_256 {
    fn encrypt_block(&self, block: &mut [u8; 16]) {
        Speck128_256::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 16]) {
        Speck128_256::decrypt_block(self, block);
    }
}

impl BlockCipher<8> for Speck64_128 {
    fn encrypt_block(&self, block: &mut [u8; 8]) {
        Speck64_128::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 8]) {
        Speck64_128::decrypt_block(self, block);
    }
}

impl BlockCipher<16> for Simon128_128 {
    fn encrypt_block(&self, block: &mut [u8; 16]) {
        Simon128_128::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 16]) {
        Simon128_128::decrypt_block(self, block);
    }
}

impl BlockCipher<16> for Simon128_256 {
    fn encrypt_block(&self, block: &mut [u8; 16]) {
        Simon128_256::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 16]) {
        Simon128_256::decrypt_block(self, block);
    }
}

impl BlockCipher<8> for Simon64_128 {
    fn encrypt_block(&self, block: &mut [u8; 8]) {
        Simon64_128::encrypt_block(self, block);
    }
    fn decrypt_block(&self, block: &mut [u8; 8]) {
        Simon64_128::decrypt_block(self, block);
    }
}

// ── 128-bit-block ciphers (BlockCipher<16>) ─────────────────────────

macro_rules! bc16 {
    ($t:ty) => {
        impl BlockCipher<16> for $t {
            fn encrypt_block(&self, block: &mut [u8; 16]) {
                <$t>::encrypt_block(self, block);
            }
            fn decrypt_block(&self, block: &mut [u8; 16]) {
                <$t>::decrypt_block(self, block);
            }
        }
    };
}

bc16!(Aria128);
bc16!(Aria192);
bc16!(Aria256);
bc16!(Camellia128);
bc16!(Camellia192);
bc16!(Camellia256);
bc16!(Gift128);
bc16!(Lucifer);
bc16!(Mars);
bc16!(Rc6);
bc16!(Seed);
bc16!(Skinny128_128);
bc16!(Square);
bc16!(Twofish);

// ── 64-bit-block ciphers (BlockCipher<8>) ───────────────────────────

macro_rules! bc8 {
    ($t:ty) => {
        impl BlockCipher<8> for $t {
            fn encrypt_block(&self, block: &mut [u8; 8]) {
                <$t>::encrypt_block(self, block);
            }
            fn decrypt_block(&self, block: &mut [u8; 8]) {
                <$t>::decrypt_block(self, block);
            }
        }
    };
}

bc8!(Blowfish);
bc8!(Cast5);
bc8!(Feal8);
bc8!(Gift64);
bc8!(Hight);
bc8!(Idea);
bc8!(Kasumi);
bc8!(Misty1);
bc8!(Present80);
bc8!(Present128);
bc8!(SaferK64);
bc8!(Skinny64_128);
bc8!(Skipjack);
bc8!(Tea);
bc8!(Xtea);

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

    /// Simon and Speck round-trip through the generic `BlockCipher`
    /// trait, covering both 8- and 16-byte block sizes.
    #[test]
    fn simon_speck_blockcipher_round_trip() {
        let s128_128 = Speck128_128::new(&[0x11u8; 16]);
        let s128_256 = Speck128_256::new(&[0x22u8; 32]);
        let s64_128 = Speck64_128::new(&[0x33u8; 16]);
        assert!(round_trip(&s128_128, [3u8; 16]));
        assert!(round_trip(&s128_256, [4u8; 16]));

        let m128_128 = Simon128_128::new(&[0x55u8; 16]);
        let m128_256 = Simon128_256::new(&[0x66u8; 32]);
        let m64_128 = Simon64_128::new(&[0x77u8; 16]);
        assert!(round_trip(&m128_128, [7u8; 16]));
        assert!(round_trip(&m128_256, [8u8; 16]));

        // 8-byte-block variants
        for cipher in [
            (&s64_128 as &dyn BlockCipher<8>),
            (&m64_128 as &dyn BlockCipher<8>),
        ] {
            let mut block = [0x5Au8; 8];
            let orig = block;
            cipher.encrypt_block(&mut block);
            assert_ne!(block, orig);
            cipher.decrypt_block(&mut block);
            assert_eq!(block, orig);
        }
    }

    fn rt8<C: BlockCipher<8>>(c: &C) {
        let mut block = [0x5Au8; 8];
        let orig = block;
        c.encrypt_block(&mut block);
        assert_ne!(block, orig, "encrypt was identity");
        c.decrypt_block(&mut block);
        assert_eq!(block, orig, "decrypt failed to invert");
    }

    fn rt16<C: BlockCipher<16>>(c: &C) {
        let mut block = [0x5Au8; 16];
        let orig = block;
        c.encrypt_block(&mut block);
        assert_ne!(block, orig, "encrypt was identity");
        c.decrypt_block(&mut block);
        assert_eq!(block, orig, "decrypt failed to invert");
    }

    /// Every 128-bit-block cipher we ship round-trips through the trait.
    #[test]
    fn all_128bit_block_ciphers_via_trait() {
        rt16(&Aria128::new(&[0u8; 16]));
        rt16(&Aria192::new(&[0u8; 24]));
        rt16(&Aria256::new(&[0u8; 32]));
        rt16(&Camellia128::new(&[0u8; 16]));
        rt16(&Camellia192::new(&[0u8; 24]));
        rt16(&Camellia256::new(&[0u8; 32]));
        rt16(&Gift128::new(&[0u8; 16]));
        rt16(&Lucifer::new(&[0u8; 16]));
        rt16(&Mars::new(&[0u8; 16]).unwrap());
        rt16(&Rc6::new(&[0u8; 16]).unwrap());
        rt16(&Seed::new(&[0u8; 16]));
        rt16(&Skinny128_128::new(&[0u8; 16]));
        rt16(&Square::new(&[0u8; 16]));
        rt16(&Twofish::new(&[0u8; 16]).unwrap());
    }

    /// Every 64-bit-block cipher we ship round-trips through the trait.
    #[test]
    fn all_64bit_block_ciphers_via_trait() {
        rt8(&Blowfish::new(&[0u8; 16]).unwrap());
        rt8(&Cast5::new(&[0u8; 16]).unwrap());
        rt8(&Feal8::new(&[0u8; 8]));
        rt8(&Gift64::new(&[0u8; 16]));
        rt8(&Hight::new(&[0u8; 16]));
        rt8(&Idea::new(&[0u8; 16]));
        rt8(&Kasumi::new(&[0u8; 16]));
        rt8(&Misty1::new(&[0u8; 16]));
        rt8(&Present80::new(&[0u8; 10]));
        rt8(&Present128::new(&[0u8; 16]));
        rt8(&SaferK64::new(&[0u8; 8]));
        rt8(&Skinny64_128::new(&[0u8; 16]));
        rt8(&Skipjack::new(&[0u8; 10]));
        rt8(&Tea::new(&[0u8; 16]));
        rt8(&Xtea::new(&[0u8; 16]));
    }
}
