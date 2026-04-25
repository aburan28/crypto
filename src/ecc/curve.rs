//! Elliptic curve parameter sets.
//!
//! A `CurveParams` bundles the Weierstrass curve coefficients (a, b),
//! the field prime p, the generator point G, the group order n, and the
//! cofactor h.  All arithmetic is performed using `FieldElement` (for
//! field operations) or plain `BigUint` (for scalar operations mod n).

use super::field::FieldElement;
use super::point::Point;
use num_bigint::BigUint;

/// Parameters for a short Weierstrass curve y² ≡ x³ + ax + b (mod p).
#[derive(Clone, Debug)]
pub struct CurveParams {
    /// Human-readable name
    pub name: &'static str,
    /// Field prime
    pub p: BigUint,
    /// Curve coefficient a
    pub a: BigUint,
    /// Curve coefficient b
    pub b: BigUint,
    /// Generator x-coordinate
    pub gx: BigUint,
    /// Generator y-coordinate
    pub gy: BigUint,
    /// Group order (number of points on the curve)
    pub n: BigUint,
    /// Cofactor (usually 1 for NIST/secp curves)
    pub h: u32,
}

impl CurveParams {
    // ── Convenience constructors ──────────────────────────────────────────────

    /// secp256k1 — the curve used by Bitcoin and Ethereum.
    ///
    /// y² = x³ + 7 (mod p), a=0, b=7.
    /// Defined in SEC 2, section 2.4.1.
    pub fn secp256k1() -> Self {
        CurveParams {
            name: "secp256k1",
            p: BigUint::parse_bytes(
                b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16,
            ).unwrap(),
            a: BigUint::from(0u32),
            b: BigUint::from(7u32),
            gx: BigUint::parse_bytes(
                b"79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16,
            ).unwrap(),
            gy: BigUint::parse_bytes(
                b"483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16,
            ).unwrap(),
            n: BigUint::parse_bytes(
                b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16,
            ).unwrap(),
            h: 1,
        }
    }

    /// NIST P-256 (secp256r1) — widely used in TLS and FIDO2.
    ///
    /// y² = x³ - 3x + b (mod p), a = p - 3.
    pub fn p256() -> Self {
        CurveParams {
            name: "P-256",
            p: BigUint::parse_bytes(
                b"ffffffff00000001000000000000000000000000ffffffffffffffffffffffff", 16,
            ).unwrap(),
            a: BigUint::parse_bytes(
                b"ffffffff00000001000000000000000000000000fffffffffffffffffffffffc", 16,
            ).unwrap(),
            b: BigUint::parse_bytes(
                b"5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b", 16,
            ).unwrap(),
            gx: BigUint::parse_bytes(
                b"6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296", 16,
            ).unwrap(),
            gy: BigUint::parse_bytes(
                b"4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5", 16,
            ).unwrap(),
            n: BigUint::parse_bytes(
                b"ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551", 16,
            ).unwrap(),
            h: 1,
        }
    }

    // ── Helpers ───────────────────────────────────────────────────────────────

    /// Wrap a raw value as a `FieldElement` over this curve's prime field.
    pub fn fe(&self, v: BigUint) -> FieldElement {
        FieldElement::new(v, self.p.clone())
    }

    /// Return the curve's `a` coefficient as a `FieldElement`.
    pub fn a_fe(&self) -> FieldElement {
        self.fe(self.a.clone())
    }

    /// Return the generator point G.
    pub fn generator(&self) -> Point {
        Point::Affine {
            x: self.fe(self.gx.clone()),
            y: self.fe(self.gy.clone()),
        }
    }

    /// Check whether a point satisfies the curve equation y² ≡ x³ + ax + b (mod p).
    pub fn is_on_curve(&self, point: &Point) -> bool {
        match point {
            Point::Infinity => true,
            Point::Affine { x, y } => {
                let lhs = y.mul(y); // y²
                let a = self.a_fe();
                let rhs = x.mul(&x.mul(x))   // x³
                    .add(&a.mul(x))           // + ax
                    .add(&self.fe(self.b.clone())); // + b
                lhs == rhs
            }
        }
    }
}
