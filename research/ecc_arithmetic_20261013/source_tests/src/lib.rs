#[path = "../../../../src/ecc/f25519.rs"]
mod f25519;
#[path = "../../../../src/ecc/curve25519.rs"]
pub mod curve25519;
#[path = "../../../../src/ecc/x25519.rs"]
pub mod x25519;
#[path = "../../../../src/ecc/p521_field.rs"]
pub mod p521_field;

// Match the original reference path used by x25519.rs tests.
pub mod ecc {
    pub use crate::curve25519;
}
