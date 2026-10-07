use num_bigint::{BigInt, BigUint};

use crate::Result;

pub const CURVE_UID: &str =
    "urn:ec-record:1:sha256:5531c4a08bdb64b6e86a6e30e9a08aa57edef7af15ac5f6d4d2a83a53bf2f646";
pub const P_DEC: &str = "6277101735386680763835789423207666416083908700390324961279";
pub const N_DEC: &str = "6277101735386680763835789423176059013767194773182842284081";
pub const T_DEC: &str = "31607402316713927207482677199";
pub const D_DEC: &str = "-24109379060336110122544161233113975664949272517896865359515";
pub const ALGEBRAIC_BOUND: u64 = 65_521;
pub const MAPPABLE_BOUND: u64 = 113;
pub const LARGE_PRIME_BOUND: u64 = 2_147_483_647;
pub const REF0_SHELL_ID: &str = "REF-0";
pub const REF0_COUNT: u64 = 1_025;
pub const LOGICAL_SHARD_RECORDS: u64 = 1_048_576;

pub fn bigint(text: &str, name: &str) -> Result<BigInt> {
    BigInt::parse_bytes(text.as_bytes(), 10).ok_or_else(|| format!("invalid decimal {name}"))
}

pub fn biguint(text: &str, name: &str) -> Result<BigUint> {
    BigUint::parse_bytes(text.as_bytes(), 10)
        .ok_or_else(|| format!("invalid unsigned decimal {name}"))
}

pub fn p() -> BigInt {
    bigint(P_DEC, "p").expect("frozen p")
}

pub fn t() -> BigInt {
    bigint(T_DEC, "t").expect("frozen t")
}

pub fn d() -> BigInt {
    bigint(D_DEC, "D").expect("frozen D")
}

pub fn subgroup_order() -> BigInt {
    bigint(N_DEC, "n").expect("frozen n")
}
