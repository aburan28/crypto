//! # `bitop/v1`: one unit of work for every ECDLP method.
//!
//! Every cost in this repository used to be charged in a unit that depended
//! on the method or the host: group-addition equivalents converted with
//! nanosecond ratios from one machine, solver work priced by wall time or not
//! at all, Callgrind instructions, reference picoseconds.  Two numbers in two
//! of those units cannot be divided, and an exponent fitted in one is not an
//! exponent in another.
//!
//! This module is the replacement: a **bit operation** is one two-input
//! Boolean gate (AND, XOR, OR) in the textbook circuit for the primitive.
//! Each primitive the methods count — a field multiplication, an affine
//! addition, a multiply-subtract modulo `r`, a 64-bit word operation, a hash
//! probe — has a frozen gate count that is a function of the operand width
//! alone, so
//!
//! ```text
//!     work = Σ  count(primitive) · gates(primitive, width)
//! ```
//!
//! is an exact integer that depends on the counts and this file and on
//! nothing else: not the host, not the compiler, not the run.  The
//! conventions and their derivations are in `docs/units/BITOPS.md`; the
//! test vectors in `docs/units/bitop-v1.vectors.json` are what any other
//! implementation (the cryptanalysis repository's Python harness, cairn's
//! evaluator) must reproduce.
//!
//! What the unit deliberately leaves out — memory traffic, table size, branch
//! and loop overhead — is reported beside it, never folded in.  Work a method
//! does that no counter sees is caught by the instruction audit
//! (`ecbench audit-ir`): the Callgrind instructions per counted bit operation
//! of a phase are compared with the reference walk's, and a phase far above
//! it is doing work the counters miss.
//!
//! A new convention is a new version (`bitop/v2`), never an edit of this one:
//! a committed figure must reprice to the same integer forever.

use serde::{Deserialize, Serialize};

/// The unit's name, as records carry it.
pub const UNIT: &str = "bitop/v1";

/// The field the curve is defined over.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Field {
    /// `F_p`, `p` of `bits` bits.
    Prime { bits: u32 },
    /// `F_{2^n}` in a polynomial basis.
    Binary { n: u32 },
    /// `F_{p^4}`, `p` of `bits` bits.
    Fp4 { bits: u32 },
}

/// Gates of a `w`-bit ripple-carry addition: `w` full adders of five gates.
pub fn int_add(w: u32) -> u128 {
    5 * w as u128
}

/// Gates of a `w × w`-bit schoolbook product: `w²` partial-product ANDs and
/// `w − 1` additions of a `w`-bit row.
pub fn int_mul(w: u32) -> u128 {
    let w = w as u128;
    w * w + 5 * w * w.saturating_sub(1)
}

/// The primitives every counter is expressed in.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Prim {
    /// Addition or subtraction in the base field.
    FieldAdd,
    FieldMul,
    FieldSqr,
    FieldInv,
    /// Square root in `F_p` (Tonelli–Shanks).
    FieldSqrt,
    /// Legendre (Jacobi) symbol in `F_p`.
    Legendre,
    /// Solve `u² + u = c` in `F_{2^n}` (half-trace as a matrix–vector product).
    AsSolve,
    /// Frobenius-orbit canonicalisation of an abscissa: normal-basis
    /// conversion and the least of `n` rotations.
    Canon,
    /// Affine group addition on the curve.
    GroupAdd,
    /// Affine doubling.
    GroupDouble,
    /// Multiplication modulo the subgroup order `r`.
    ModrMul,
    /// Addition modulo `r`.
    ModrAdd,
    /// One 64-bit word operation (XOR, AND, OR, compare) of a bit-packed
    /// algorithm: 64 gates.
    WordOp,
    /// A hash-table probe or insert of a key as wide as a field element:
    /// a fold-XOR hash, a key comparison and a select, `4κ` gates.  The
    /// keys are field elements, already uniform, so no multiplier is needed;
    /// what a table really costs is memory, which the unit reports beside
    /// the gates and never folds in.
    HashProbe,
}

impl Prim {
    pub const ALL: [Prim; 14] = [
        Prim::FieldAdd,
        Prim::FieldMul,
        Prim::FieldSqr,
        Prim::FieldInv,
        Prim::FieldSqrt,
        Prim::Legendre,
        Prim::AsSolve,
        Prim::Canon,
        Prim::GroupAdd,
        Prim::GroupDouble,
        Prim::ModrMul,
        Prim::ModrAdd,
        Prim::WordOp,
        Prim::HashProbe,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Prim::FieldAdd => "field_add",
            Prim::FieldMul => "field_mul",
            Prim::FieldSqr => "field_sqr",
            Prim::FieldInv => "field_inv",
            Prim::FieldSqrt => "field_sqrt",
            Prim::Legendre => "legendre",
            Prim::AsSolve => "as_solve",
            Prim::Canon => "canon",
            Prim::GroupAdd => "group_add",
            Prim::GroupDouble => "group_double",
            Prim::ModrMul => "modr_mul",
            Prim::ModrAdd => "modr_add",
            Prim::WordOp => "word_op",
            Prim::HashProbe => "hash_probe",
        }
    }

    pub fn from_name(s: &str) -> Option<Prim> {
        Prim::ALL.into_iter().find(|p| p.name() == s)
    }
}

/// What a price depends on: the field, and the width of `r`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Widths {
    pub field: Field,
    /// Bits of the subgroup order `r`.
    pub r_bits: u32,
}

/// Why a primitive has no price in a field (an Artin–Schreier solve over
/// `F_p`, say): the counter that asked for it is misclassified.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NoPrice {
    pub prim: Prim,
    pub field: Field,
}

impl std::fmt::Display for NoPrice {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} has no bitop/v1 price over {:?}",
            self.prim.name(),
            self.field
        )
    }
}

/// Prime-field primitives, `p` of `b` bits, Montgomery representation.
mod fp {
    use super::{int_add, int_mul};
    /// Addition, a conditional subtraction of `p`, and a `b`-bit select
    /// (three gates a bit).
    pub fn add(b: u32) -> u128 {
        2 * int_add(b) + 3 * b as u128
    }
    /// Montgomery product: the product, the reduction's product, and the
    /// final conditional subtraction.
    pub fn mul(b: u32) -> u128 {
        2 * int_mul(b) + add(b)
    }
    /// A squaring is charged as a product.
    pub fn sqr(b: u32) -> u128 {
        mul(b)
    }
    /// Binary extended Euclid: `2b` steps of a subtraction, a comparison and
    /// a select on `b`-bit words.
    pub fn inv(b: u32) -> u128 {
        2 * b as u128 * add(b)
    }
    /// Tonelli–Shanks: an exponentiation (`b` squarings and `b/2`
    /// products) and an allowance of `b/2` products for the 2-adic loop.
    pub fn sqrt(b: u32) -> u128 {
        2 * b as u128 * mul(b)
    }
    /// The Jacobi algorithm, a gcd: charged as an inversion.
    pub fn legendre(b: u32) -> u128 {
        inv(b)
    }
}

/// Binary-field primitives, `F_{2^n}` in a polynomial basis reduced by a
/// pentanomial (four low terms; a trinomial would be cheaper, and the
/// convention takes the dearer one so no curve is flattered by its modulus).
mod f2 {
    const LOW_TERMS: u128 = 4;
    pub fn add(n: u32) -> u128 {
        n as u128
    }
    /// Schoolbook carry-less product (`n²` ANDs, `(n−1)²` XORs) and the
    /// reduction of its `n − 1` high bits.
    pub fn mul(n: u32) -> u128 {
        let n = n as u128;
        n * n + (n - 1) * (n - 1) + LOW_TERMS * (n - 1)
    }
    /// Squaring spreads bits (wiring) and reduces.
    pub fn sqr(n: u32) -> u128 {
        LOW_TERMS * (n as u128 - 1)
    }
    /// Polynomial extended Euclid: `2n` steps, each a shift-XOR of two
    /// `n`-bit pairs and an `n`-bit select.
    pub fn inv(n: u32) -> u128 {
        let n = n as u128;
        2 * n * (2 * n + 3 * n)
    }
    /// Half-trace as a precomputed `n × n` matrix over `F_2` times a vector.
    pub fn as_solve(n: u32) -> u128 {
        let n = n as u128;
        n * n + n * (n - 1)
    }
    /// Normal-basis conversion (a matrix–vector product) and the least of
    /// `n` rotations, two gates a compared bit.
    pub fn canon(n: u32) -> u128 {
        as_solve(n) + 2 * (n as u128) * (n as u128)
    }
}

/// Field operations behind one group operation.
struct GroupFormula {
    inv: u128,
    mul: u128,
    sqr: u128,
    add: u128,
}

/// `y² = x³ + ax + b`, affine: `λ = (y₂−y₁)/(x₂−x₁)`, `x₃ = λ² − x₁ − x₂`,
/// `y₃ = λ(x₁ − x₃) − y₁`.
const PRIME_ADD: GroupFormula = GroupFormula {
    inv: 1,
    mul: 2,
    sqr: 1,
    add: 6,
};
/// `λ = (3x₁² + a)/2y₁`, `x₃ = λ² − 2x₁`, `y₃ = λ(x₁ − x₃) − y₁`.
const PRIME_DOUBLE: GroupFormula = GroupFormula {
    inv: 1,
    mul: 2,
    sqr: 2,
    add: 8,
};
/// `y² + xy = x³ + ax² + b`, affine: `λ = (y₁+y₂)/(x₁+x₂)`,
/// `x₃ = λ² + λ + x₁ + x₂ + a`, `y₃ = λ(x₁ + x₃) + x₃ + y₁`.
const BINARY_ADD: GroupFormula = GroupFormula {
    inv: 1,
    mul: 2,
    sqr: 1,
    add: 9,
};
/// `λ = x₁ + y₁/x₁`, `x₃ = λ² + λ + a`, `y₃ = x₁² + (λ + 1)x₃`.
const BINARY_DOUBLE: GroupFormula = GroupFormula {
    inv: 1,
    mul: 2,
    sqr: 2,
    add: 5,
};
/// `E(F_{p⁴})`: one group operation is 97 `F_p` products, the convention of
/// `ecbench::fp4` (`C_ADD`), which charges the extension inversion through
/// the norm.
const FP4_GROUP_OP_FP_MULS: u128 = 97;
/// An `F_{p⁴}` product is 16 `F_p` products (schoolbook), as `ecbench::fp4`
/// charges it.
const FP4_MUL_FP_MULS: u128 = 16;

fn group(f: &GroupFormula, inv: u128, mul: u128, sqr: u128, add: u128) -> u128 {
    f.inv * inv + f.mul * mul + f.sqr * sqr + f.add * add
}

impl Widths {
    /// Gates of one `prim` at these widths.
    pub fn gates(&self, prim: Prim) -> Result<u128, NoPrice> {
        let none = || NoPrice {
            prim,
            field: self.field,
        };
        let w = self.r_bits.max(1);
        // Arithmetic modulo r is prime-field arithmetic at r's width.
        match prim {
            Prim::ModrMul => return Ok(fp::mul(w)),
            Prim::ModrAdd => return Ok(fp::add(w)),
            Prim::WordOp => return Ok(64),
            _ => {}
        }
        Ok(match self.field {
            Field::Prime { bits: b } => match prim {
                Prim::FieldAdd => fp::add(b),
                Prim::FieldMul => fp::mul(b),
                Prim::FieldSqr => fp::sqr(b),
                Prim::FieldInv => fp::inv(b),
                Prim::FieldSqrt => fp::sqrt(b),
                Prim::Legendre => fp::legendre(b),
                Prim::GroupAdd => group(&PRIME_ADD, fp::inv(b), fp::mul(b), fp::sqr(b), fp::add(b)),
                Prim::GroupDouble => group(
                    &PRIME_DOUBLE,
                    fp::inv(b),
                    fp::mul(b),
                    fp::sqr(b),
                    fp::add(b),
                ),
                Prim::HashProbe => 4 * b as u128,
                Prim::AsSolve | Prim::Canon => return Err(none()),
                Prim::ModrMul | Prim::ModrAdd | Prim::WordOp => unreachable!(),
            },
            Field::Binary { n } => match prim {
                Prim::FieldAdd => f2::add(n),
                Prim::FieldMul => f2::mul(n),
                Prim::FieldSqr => f2::sqr(n),
                Prim::FieldInv => f2::inv(n),
                Prim::AsSolve => f2::as_solve(n),
                Prim::Canon => f2::canon(n),
                Prim::GroupAdd => {
                    group(&BINARY_ADD, f2::inv(n), f2::mul(n), f2::sqr(n), f2::add(n))
                }
                Prim::GroupDouble => group(
                    &BINARY_DOUBLE,
                    f2::inv(n),
                    f2::mul(n),
                    f2::sqr(n),
                    f2::add(n),
                ),
                Prim::HashProbe => 4 * n as u128,
                Prim::FieldSqrt | Prim::Legendre => return Err(none()),
                Prim::ModrMul | Prim::ModrAdd | Prim::WordOp => unreachable!(),
            },
            Field::Fp4 { bits: b } => match prim {
                Prim::FieldAdd => 4 * fp::add(b),
                Prim::FieldMul | Prim::FieldSqr => FP4_MUL_FP_MULS * fp::mul(b),
                Prim::GroupAdd | Prim::GroupDouble => FP4_GROUP_OP_FP_MULS * fp::mul(b),
                Prim::HashProbe => 16 * b as u128,
                Prim::FieldInv | Prim::FieldSqrt | Prim::Legendre | Prim::AsSolve | Prim::Canon => {
                    return Err(none())
                }
                Prim::ModrMul | Prim::ModrAdd | Prim::WordOp => unreachable!(),
            },
        })
    }

    /// Bits of a field element: the hash key width and the field size the
    /// tables report.
    pub fn element_bits(&self) -> u32 {
        match self.field {
            Field::Prime { bits } => bits,
            Field::Binary { n } => n,
            Field::Fp4 { bits } => 4 * bits,
        }
    }
}

/// `log₂` of an exact gate count.
pub fn log2(x: u128) -> f64 {
    if x == 0 {
        return f64::NEG_INFINITY;
    }
    let lz = x.leading_zeros();
    let shift = (128 - lz).saturating_sub(53);
    ((x >> shift) as f64).log2() + shift as f64
}

/// The test vectors every implementation of `bitop/v1` must reproduce:
/// each primitive's gates at a spread of widths in each field family.
pub fn vectors() -> serde_json::Value {
    let mut out = Vec::new();
    let fields = [
        Field::Prime { bits: 17 },
        Field::Prime { bits: 32 },
        Field::Prime { bits: 256 },
        Field::Binary { n: 23 },
        Field::Binary { n: 41 },
        Field::Binary { n: 131 },
        Field::Fp4 { bits: 13 },
    ];
    for field in fields {
        for r_bits in [16u32, 40, 130] {
            let w = Widths { field, r_bits };
            for p in Prim::ALL {
                if let Ok(g) = w.gates(p) {
                    out.push(serde_json::json!({
                        "field": field,
                        "r_bits": r_bits,
                        "prim": p.name(),
                        "gates": g.to_string(),
                    }));
                }
            }
        }
    }
    serde_json::json!({
        "schema": "bitop.vectors/v1",
        "unit": UNIT,
        "what": "Gates of each bitop/v1 primitive at the listed widths. Generated by `ecbench bitop-vectors` from src/cryptanalysis/bitops.rs; any other implementation must reproduce every row exactly.",
        "vectors": out,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_widths_by_hand() {
        // 8-bit integers: 64 ANDs and 7 eight-bit additions of 40 gates.
        assert_eq!(int_mul(8), 64 + 7 * 40);
        // F_{2^3}: 9 ANDs, 4 XORs, reduction 4·2.
        assert_eq!(f2::mul(3), 9 + 4 + 8);
        let w = Widths {
            field: Field::Binary { n: 3 },
            r_bits: 3,
        };
        // I + 2M + S + 9A = 2·3·15 + 2·21 + 8 + 27.
        assert_eq!(w.gates(Prim::GroupAdd).unwrap(), 90 + 42 + 8 + 27);
    }

    #[test]
    fn every_primitive_prices_where_it_should() {
        let p = Widths {
            field: Field::Prime { bits: 32 },
            r_bits: 32,
        };
        let b = Widths {
            field: Field::Binary { n: 41 },
            r_bits: 39,
        };
        for prim in Prim::ALL {
            let in_prime = p.gates(prim).is_ok();
            let in_binary = b.gates(prim).is_ok();
            match prim {
                Prim::AsSolve | Prim::Canon => assert!(!in_prime && in_binary),
                Prim::FieldSqrt | Prim::Legendre => assert!(in_prime && !in_binary),
                _ => assert!(in_prime && in_binary, "{}", prim.name()),
            }
            assert_eq!(Prim::from_name(prim.name()), Some(prim));
        }
    }

    #[test]
    fn group_operations_dominate_field_products() {
        for n in [13u32, 41, 83, 131] {
            let w = Widths {
                field: Field::Binary { n },
                r_bits: n,
            };
            let add = w.gates(Prim::GroupAdd).unwrap();
            assert!(add > 4 * w.gates(Prim::FieldMul).unwrap());
            assert!(add > w.gates(Prim::ModrMul).unwrap());
        }
    }

    #[test]
    fn log2_is_exact_on_powers_and_close_elsewhere() {
        assert_eq!(log2(1 << 100), 100.0);
        assert!((log2(3) - 3f64.log2()).abs() < 1e-15);
        assert!((log2(u128::MAX) - 128.0).abs() < 1e-12);
    }
}
