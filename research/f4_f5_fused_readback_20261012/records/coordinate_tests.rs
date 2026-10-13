//! Compile the actual coordinate module's tests natively before rebuilding
//! the large root test target. All other modules come from the locked rlib.
pub use crypto_lib::{binary_ecc, ecc, hash, isogeny};
pub mod cryptanalysis {
    pub use crypto_lib::cryptanalysis::*;
    #[path = "/workspace/crypto/src/cryptanalysis/curve_traits/mod.rs"]
    pub mod curve_traits;
}
