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
    KyberCiphertext, KyberPrivateKey, KyberPublicKey,
    kyber_decapsulate, kyber_encapsulate, kyber_keygen,
};

pub use mceliece::{
    McElieceKeyPair, McEliecePrivateKey, McElieceePublicKey,
    mceliece_decrypt, mceliece_encrypt,
};

pub use ntru::{
    NtruKeyPair, NtruPrivateKey, NtruPublicKey,
    ntru_decapsulate, ntru_encapsulate, ntru_keygen,
};

pub use frodo::{
    FrodoCiphertext, FrodoKeyPair, FrodoPrivateKey, FrodoPublicKey,
    frodo_decapsulate, frodo_encapsulate, frodo_keygen,
};

pub use hqc::{
    HqcCiphertext, HqcKeyPair, HqcPrivateKey, HqcPublicKey,
    hqc_decapsulate, hqc_encapsulate, hqc_keygen,
};

pub use x_wing::{
    XWingCiphertext, XWingKeyPair, XWingPrivateKey, XWingPublicKey,
    x_wing_decapsulate, x_wing_encapsulate, x_wing_keygen,
};
