//! Post-quantum cryptography.
//!
//! NIST standards: ML-KEM (FIPS 203, `ml_kem`), ML-DSA (FIPS 204,
//! `ml_dsa`), SLH-DSA (FIPS 205, `slh_dsa`) — all validated against
//! official test vectors where the module docs say so.
//!
//! Other families: lattice KEMs (toy Kyber, NTRU, NTRU Prime,
//! FrodoKEM), code-based KEMs (McEliece, Classic McEliece, BIKE, HQC),
//! hybrid (X-Wing), and isogeny-based schemes (CSIDH key exchange,
//! SQIsign signatures — both with toy parameters).
//!
//! See individual modules for algorithm documentation and limitations.

pub mod bike;
pub mod classic_mceliece;
pub mod csidh;
pub mod frodo;
pub mod hqc;
pub mod kyber;
pub mod mceliece;
pub mod ml_dsa;
pub mod ml_kem;
pub mod ntru;
pub mod ntru_prime;
pub mod slh_dsa;
pub mod sqisign;
pub mod x_wing;

pub use kyber::{
    kyber_decapsulate, kyber_encapsulate, kyber_keygen, KyberCiphertext, KyberPrivateKey,
    KyberPublicKey,
};

pub use ml_kem::{
    ml_kem_decaps, ml_kem_encaps, ml_kem_encaps_internal, ml_kem_keygen, ml_kem_keygen_internal,
    MlKemDecapsKey, MlKemEncapsKey, MlKemParams, ML_KEM_1024, ML_KEM_512, ML_KEM_768,
};

pub use ml_dsa::{ml_dsa_65_keygen, ml_dsa_65_sign, ml_dsa_65_verify, MlDsaPublicKey, MlDsaSecretKey};

pub use sqisign::{
    sqisign_keygen, sqisign_sign, sqisign_verify, SqiSignPublicKey, SqiSignSecretKey, SqiSignature,
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
