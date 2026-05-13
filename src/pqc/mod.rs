//! Post-quantum cryptography: lattice-based (ML-KEM/Kyber, NTRU,
//! FrodoKEM), code-based (McEliece, HQC), and hybrid (X-Wing)
//! key-encapsulation mechanisms.
//!
//! See individual modules for algorithm documentation and limitations.

pub mod bike;
pub mod classic_mceliece;
pub mod csidh;
pub mod frodo;
pub mod hqc;
pub mod kyber;
pub mod mceliece;
pub mod ntru;
pub mod ntru_prime;
pub mod x_wing;

pub use kyber::{
    kyber_decapsulate, kyber_encapsulate, kyber_keygen, KyberCiphertext, KyberPrivateKey,
    KyberPublicKey,
};

pub use mceliece::{
    mceliece_decrypt, mceliece_encrypt, McElieceKeyPair, McEliecePrivateKey, McElieceePublicKey,
};

pub use ntru::{
    ntru_decapsulate, ntru_encapsulate, ntru_keygen, NtruKeyPair, NtruPrivateKey, NtruPublicKey,
};

pub use frodo::{
    frodo_decapsulate, frodo_encapsulate, frodo_keygen, FrodoCiphertext, FrodoKeyPair,
    FrodoPrivateKey, FrodoPublicKey,
};

pub use hqc::{
    hqc_decapsulate, hqc_encapsulate, hqc_keygen, HqcCiphertext, HqcKeyPair, HqcPrivateKey,
    HqcPublicKey,
};

pub use x_wing::{
    x_wing_decapsulate, x_wing_encapsulate, x_wing_keygen, XWingCiphertext, XWingKeyPair,
    XWingPrivateKey, XWingPublicKey,
};
