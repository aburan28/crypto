//! Message encryption built from the Eigen-Séta KEM and ChaCha20-Poly1305.
//!
//! Wire format: fixed-length KEM ciphertext || 12-byte random nonce ||
//! authenticated ciphertext (including its 16-byte tag). Associated data is
//! caller-supplied and must match during decryption. This is research code;
//! the KEM's conjectural security limits the PKE's security.

use crate::kem::{Ciphertext as KemCiphertext, EigenSeta, PublicKey, SecretKey};
use crate::params::ParamSet;
use chacha20poly1305::aead::{Aead, KeyInit, Payload};
use chacha20poly1305::{ChaCha20Poly1305, Nonce};
use rand_core::{CryptoRng, RngCore};
use zeroize::Zeroize;

pub const NONCE_BYTES: usize = 12;
pub const TAG_BYTES: usize = 16;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PkeCiphertext {
    pub bytes: Vec<u8>,
}

pub struct EigenSetaPke<P: ParamSet<N>, const N: usize>(core::marker::PhantomData<P>);

impl<P: ParamSet<N>, const N: usize> EigenSetaPke<P, N> {
    /// Encrypt a message under a validated public key. Returns `None` if the
    /// public key is invalid or AEAD encryption fails.
    pub fn encrypt<R: RngCore + CryptoRng>(
        pk: &PublicKey,
        plaintext: &[u8],
        associated_data: &[u8],
        rng: &mut R,
    ) -> Option<PkeCiphertext> {
        let (kem_ct, mut key) = EigenSeta::<P, N>::encaps(pk, rng)?;
        let cipher = ChaCha20Poly1305::new((&key).into());
        let mut nonce = [0u8; NONCE_BYTES];
        rng.fill_bytes(&mut nonce);
        let sealed = cipher.encrypt(
            Nonce::from_slice(&nonce),
            Payload {
                msg: plaintext,
                aad: associated_data,
            },
        );
        key.zeroize();
        let sealed = sealed.ok()?;
        let mut bytes = Vec::with_capacity(kem_ct.bytes.len() + NONCE_BYTES + sealed.len());
        bytes.extend_from_slice(&kem_ct.bytes);
        bytes.extend_from_slice(&nonce);
        bytes.extend_from_slice(&sealed);
        Some(PkeCiphertext { bytes })
    }

    /// Authenticate and decrypt. Any malformed KEM ciphertext, nonce, tag,
    /// body, or mismatched associated data returns `None`.
    pub fn decrypt(
        sk: &SecretKey,
        ciphertext: &PkeCiphertext,
        associated_data: &[u8],
    ) -> Option<Vec<u8>> {
        let kem_len = EigenSeta::<P, N>::ct_bytes();
        let min_len = kem_len + NONCE_BYTES + TAG_BYTES;
        if ciphertext.bytes.len() < min_len {
            return None;
        }
        let kem_ct = KemCiphertext {
            bytes: ciphertext.bytes[..kem_len].to_vec(),
        };
        let nonce = Nonce::from_slice(&ciphertext.bytes[kem_len..kem_len + NONCE_BYTES]);
        let mut key = EigenSeta::<P, N>::decaps(sk, &kem_ct);
        let cipher = ChaCha20Poly1305::new((&key).into());
        let opened = cipher.decrypt(
            nonce,
            Payload {
                msg: &ciphertext.bytes[kem_len + NONCE_BYTES..],
                aad: associated_data,
            },
        );
        key.zeroize();
        opened.ok()
    }
}
