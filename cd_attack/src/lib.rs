//! Castryck-Decru attack scaffolding.
//!
//! Session 1: genus-2 Jacobian arithmetic.
//!   - F_p / F_{p^2} field arithmetic
//!   - F_{p^2}[x] polynomial ring (add/sub/mul/div/gcd/ext_gcd)
//!   - Hyperelliptic curve C: y^2 = f(x), deg f = 5, over F_{p^2}
//!   - Mumford-rep'd divisor classes on Jac(C)
//!   - Cantor's algorithm for the group law

pub mod field;
pub mod poly;
pub mod jacobian;
pub mod richelot;
pub mod glue;
pub mod glue_bridge;
pub mod ec;
pub mod uvtable;
pub mod attack;
