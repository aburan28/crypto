//! Fixed-width arithmetic used by the P-192 singular-companion experiment.
//!
//! This module is intentionally independent of the repository's textbook
//! `BigUint` point formulas.  Field elements are three-limb Montgomery values,
//! the legitimate P-192 group uses Jacobian coordinates, and the singular
//! companion is represented by its nodal torus.  The experiment can therefore
//! compare the result of the real low-level victim call with a separately
//! implemented arithmetic path.

use crate::ct_bignum::{MontgomeryContext, Uint};
use crate::ecc::curve::CurveParams;
use crate::ecc::point::Point;
use num_bigint::BigUint;
use num_traits::One;
use serde::{Deserialize, Serialize};

pub type U192 = Uint<3>;

pub const P_HEX: &str = "fffffffffffffffffffffffffffffffeffffffffffffffff";
pub const N_HEX: &str = "ffffffffffffffffffffffff99def836146bc9b1b4d22831";
pub const B_HEX: &str = "64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1";
pub const GX_HEX: &str = "188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012";
pub const GY_HEX: &str = "07192b95ffc8da78631011ed6b24cdd573f977a11e794811";
pub const SINGULAR_ORDER_DEC: &str = "6277101735386680763835789423207666416083908700390324961280";
pub const M_DEC: &str = "93297598515282145091939986791922449178951680";
pub const QMAX: u64 = 67_280_421_310_721;

pub const SINGULAR_DISTINCT_FACTORS: [u64; 10] =
    [2, 3, 5, 17, 257, 641, 65_537, 274_177, 6_700_417, QMAX];

fn parse_hex(s: &str) -> BigUint {
    BigUint::parse_bytes(s.as_bytes(), 16).expect("hard-coded hexadecimal integer")
}

fn parse_dec(s: &str) -> BigUint {
    BigUint::parse_bytes(s.as_bytes(), 10).expect("hard-coded decimal integer")
}

pub fn p192_prime() -> BigUint {
    parse_hex(P_HEX)
}

pub fn p192_order() -> BigUint {
    parse_hex(N_HEX)
}

pub fn singular_order() -> BigUint {
    parse_dec(SINGULAR_ORDER_DEC)
}

pub fn residue_modulus() -> BigUint {
    parse_dec(M_DEC)
}

/// Field/group-operation accounting.  These are deterministic native
/// operation counts, not wall-clock proxies.
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct NativeOpCounts {
    pub field_additions: u64,
    pub field_subtractions: u64,
    pub field_multiplications: u64,
    pub field_squarings: u64,
    pub field_inversions: u64,
    pub point_additions: u64,
    pub point_doublings: u64,
    pub scalar_multiplications: u64,
}

impl NativeOpCounts {
    pub fn merge(&mut self, other: Self) {
        self.field_additions = self.field_additions.saturating_add(other.field_additions);
        self.field_subtractions = self
            .field_subtractions
            .saturating_add(other.field_subtractions);
        self.field_multiplications = self
            .field_multiplications
            .saturating_add(other.field_multiplications);
        self.field_squarings = self.field_squarings.saturating_add(other.field_squarings);
        self.field_inversions = self.field_inversions.saturating_add(other.field_inversions);
        self.point_additions = self.point_additions.saturating_add(other.point_additions);
        self.point_doublings = self.point_doublings.saturating_add(other.point_doublings);
        self.scalar_multiplications = self
            .scalar_multiplications
            .saturating_add(other.scalar_multiplications);
    }
}

/// A canonical field element stored in Montgomery form.
#[derive(Clone, Copy, Debug)]
pub struct Fp(pub U192);

impl PartialEq for Fp {
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.0.ct_eq_full(&other.0))
    }
}

impl Eq for Fp {}

/// Runtime Montgomery context and its local operation counters.
#[derive(Clone, Debug)]
pub struct P192Arithmetic {
    ctx: MontgomeryContext<3>,
    pub counts: NativeOpCounts,
}

impl Default for P192Arithmetic {
    fn default() -> Self {
        Self::new()
    }
}

impl P192Arithmetic {
    pub fn new() -> Self {
        let p = U192::from_biguint(&p192_prime());
        Self {
            ctx: MontgomeryContext::new(p).expect("P-192 modulus is odd"),
            counts: NativeOpCounts::default(),
        }
    }

    pub fn modulus(&self) -> U192 {
        self.ctx.n
    }

    pub fn zero(&self) -> Fp {
        Fp(U192::ZERO)
    }

    pub fn one(&self) -> Fp {
        Fp(self.ctx.r_mod_n)
    }

    pub fn from_u64(&self, value: u64) -> Fp {
        Fp(self.ctx.to_montgomery(&Uint([value, 0, 0])))
    }

    pub fn from_biguint(&self, value: &BigUint) -> Result<Fp, String> {
        if value >= &p192_prime() {
            return Err("P-192 field element is not canonical".into());
        }
        Ok(Fp(self.ctx.to_montgomery(&U192::from_biguint(value))))
    }

    pub fn from_hex(&self, value: &str) -> Fp {
        self.from_biguint(&parse_hex(value))
            .expect("hard-coded P-192 field element")
    }

    pub fn canonical(&self, value: Fp) -> U192 {
        self.ctx.from_montgomery(&value.0)
    }

    pub fn to_biguint(&self, value: Fp) -> BigUint {
        self.canonical(value).to_biguint()
    }

    pub fn to_bytes(&self, value: Fp) -> [u8; 24] {
        let limbs = self.canonical(value).0;
        let mut out = [0u8; 24];
        for (i, limb) in limbs.iter().enumerate() {
            let offset = 24 - 8 * (i + 1);
            out[offset..offset + 8].copy_from_slice(&limb.to_be_bytes());
        }
        out
    }

    pub fn add(&mut self, left: Fp, right: Fp) -> Fp {
        self.counts.field_additions += 1;
        Fp(left.0.add_mod(&right.0, &self.ctx.n))
    }

    pub fn sub(&mut self, left: Fp, right: Fp) -> Fp {
        self.counts.field_subtractions += 1;
        Fp(left.0.sub_mod(&right.0, &self.ctx.n))
    }

    pub fn neg(&mut self, value: Fp) -> Fp {
        self.sub(self.zero(), value)
    }

    pub fn mul(&mut self, left: Fp, right: Fp) -> Fp {
        self.counts.field_multiplications += 1;
        Fp(self.ctx.mont_mul(&left.0, &right.0))
    }

    pub fn square(&mut self, value: Fp) -> Fp {
        self.counts.field_squarings += 1;
        Fp(self.ctx.mont_sqr(&value.0))
    }

    pub fn mul_small(&mut self, value: Fp, scalar: u64) -> Fp {
        let scalar = self.from_u64(scalar);
        self.mul(value, scalar)
    }

    pub fn pow(&mut self, base: Fp, exponent: U192) -> Fp {
        let mut acc = self.one();
        for bit in (0..192).rev() {
            acc = self.square(acc);
            if (exponent.0[bit / 64] >> (bit % 64)) & 1 == 1 {
                acc = self.mul(acc, base);
            }
        }
        acc
    }

    pub fn inverse(&mut self, value: Fp) -> Option<Fp> {
        if value == self.zero() {
            return None;
        }
        self.counts.field_inversions += 1;
        let exponent = &p192_prime() - BigUint::from(2u8);
        Some(self.pow(value, U192::from_biguint(&exponent)))
    }

    pub fn reset_counts(&mut self) {
        self.counts = NativeOpCounts::default();
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AffinePoint {
    pub x: Fp,
    pub y: Fp,
}

impl AffinePoint {
    pub fn generator(arithmetic: &P192Arithmetic) -> Self {
        Self {
            x: arithmetic.from_hex(GX_HEX),
            y: arithmetic.from_hex(GY_HEX),
        }
    }

    pub fn singular_generator(arithmetic: &P192Arithmetic) -> Self {
        Self {
            x: arithmetic.from_u64(66),
            y: arithmetic.from_u64(536),
        }
    }

    pub fn neg(self, arithmetic: &mut P192Arithmetic) -> Self {
        Self {
            x: self.x,
            y: arithmetic.neg(self.y),
        }
    }

    pub fn to_repo_point(self, arithmetic: &P192Arithmetic, curve: &CurveParams) -> Point {
        Point::Affine {
            x: curve.fe(arithmetic.to_biguint(self.x)),
            y: curve.fe(arithmetic.to_biguint(self.y)),
        }
    }

    pub fn from_repo_point(
        point: &Point,
        arithmetic: &P192Arithmetic,
    ) -> Result<Option<Self>, String> {
        match point {
            Point::Infinity => Ok(None),
            Point::Affine { x, y } => Ok(Some(Self {
                x: arithmetic.from_biguint(&x.value)?,
                y: arithmetic.from_biguint(&y.value)?,
            })),
        }
    }

    pub fn sec1_24(self, arithmetic: &P192Arithmetic) -> [u8; 49] {
        let mut out = [0u8; 49];
        out[0] = 4;
        out[1..25].copy_from_slice(&arithmetic.to_bytes(self.x));
        out[25..49].copy_from_slice(&arithmetic.to_bytes(self.y));
        out
    }
}

/// Jacobian `(X:Y:Z)` with affine interpretation `(X/Z²,Y/Z³)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JacobianPoint {
    pub x: Fp,
    pub y: Fp,
    pub z: Fp,
}

impl JacobianPoint {
    pub fn infinity(arithmetic: &P192Arithmetic) -> Self {
        Self {
            x: arithmetic.zero(),
            y: arithmetic.one(),
            z: arithmetic.zero(),
        }
    }

    pub fn from_affine(point: AffinePoint, arithmetic: &P192Arithmetic) -> Self {
        Self {
            x: point.x,
            y: point.y,
            z: arithmetic.one(),
        }
    }

    pub fn is_infinity(&self, arithmetic: &P192Arithmetic) -> bool {
        self.z == arithmetic.zero()
    }

    pub fn double(self, arithmetic: &mut P192Arithmetic) -> Self {
        arithmetic.counts.point_doublings += 1;
        if self.is_infinity(arithmetic) || self.y == arithmetic.zero() {
            return Self::infinity(arithmetic);
        }

        // E = 3(X-Z²)(X+Z²), the a=-3 specialization.
        let z2 = arithmetic.square(self.z);
        let x_minus_z2 = arithmetic.sub(self.x, z2);
        let x_plus_z2 = arithmetic.add(self.x, z2);
        let e0 = arithmetic.mul(x_minus_z2, x_plus_z2);
        let e = arithmetic.mul_small(e0, 3);
        let a = arithmetic.square(self.x);
        let b = arithmetic.square(self.y);
        let c = arithmetic.square(b);
        let x_plus_b = arithmetic.add(self.x, b);
        let x_plus_b2 = arithmetic.square(x_plus_b);
        let x_plus_b2_minus_a = arithmetic.sub(x_plus_b2, a);
        let d0 = arithmetic.sub(x_plus_b2_minus_a, c);
        let d = arithmetic.mul_small(d0, 2);
        let f = arithmetic.square(e);
        let two_d = arithmetic.mul_small(d, 2);
        let x3 = arithmetic.sub(f, two_d);
        let d_minus_x3 = arithmetic.sub(d, x3);
        let e_times = arithmetic.mul(e, d_minus_x3);
        let eight_c = arithmetic.mul_small(c, 8);
        let y3 = arithmetic.sub(e_times, eight_c);
        let yz = arithmetic.mul(self.y, self.z);
        let z3 = arithmetic.mul_small(yz, 2);
        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub fn add_affine(self, other: AffinePoint, arithmetic: &mut P192Arithmetic) -> Self {
        arithmetic.counts.point_additions += 1;
        if self.is_infinity(arithmetic) {
            return Self::from_affine(other, arithmetic);
        }

        let z1z1 = arithmetic.square(self.z);
        let u2 = arithmetic.mul(other.x, z1z1);
        let z1_cubed = arithmetic.mul(self.z, z1z1);
        let s2 = arithmetic.mul(other.y, z1_cubed);
        let h = arithmetic.sub(u2, self.x);
        let s2_minus_y1 = arithmetic.sub(s2, self.y);
        if h == arithmetic.zero() {
            if s2_minus_y1 == arithmetic.zero() {
                // Avoid double-counting this exceptional addition as a
                // separate point addition in the operation certificate.
                arithmetic.counts.point_additions -= 1;
                return self.double(arithmetic);
            }
            return Self::infinity(arithmetic);
        }

        let hh = arithmetic.square(h);
        let i = arithmetic.mul_small(hh, 4);
        let j = arithmetic.mul(h, i);
        let r = arithmetic.mul_small(s2_minus_y1, 2);
        let v = arithmetic.mul(self.x, i);
        let r2 = arithmetic.square(r);
        let r2_minus_j = arithmetic.sub(r2, j);
        let two_v = arithmetic.mul_small(v, 2);
        let x3 = arithmetic.sub(r2_minus_j, two_v);
        let v_minus_x3 = arithmetic.sub(v, x3);
        let r_times = arithmetic.mul(r, v_minus_x3);
        let y1j = arithmetic.mul(self.y, j);
        let y1j2 = arithmetic.mul_small(y1j, 2);
        let y3 = arithmetic.sub(r_times, y1j2);
        let z_plus_h = arithmetic.add(self.z, h);
        let z_plus_h2 = arithmetic.square(z_plus_h);
        let z_intermediate = arithmetic.sub(z_plus_h2, z1z1);
        let z3 = arithmetic.sub(z_intermediate, hh);
        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub fn to_affine(self, arithmetic: &mut P192Arithmetic) -> Option<AffinePoint> {
        if self.is_infinity(arithmetic) {
            return None;
        }
        let z_inv = arithmetic.inverse(self.z)?;
        let z2 = arithmetic.square(z_inv);
        let z3 = arithmetic.mul(z2, z_inv);
        Some(AffinePoint {
            x: arithmetic.mul(self.x, z2),
            y: arithmetic.mul(self.y, z3),
        })
    }
}

pub fn scalar_mul(
    base: AffinePoint,
    scalar: &BigUint,
    arithmetic: &mut P192Arithmetic,
) -> JacobianPoint {
    arithmetic.counts.scalar_multiplications += 1;
    let mut acc = JacobianPoint::infinity(arithmetic);
    for bit in (0..scalar.bits()).rev() {
        acc = acc.double(arithmetic);
        if scalar.bit(bit) {
            acc = acc.add_affine(base, arithmetic);
        }
    }
    acc
}

pub fn scalar_mul_u64(
    base: AffinePoint,
    scalar: u64,
    arithmetic: &mut P192Arithmetic,
) -> JacobianPoint {
    scalar_mul(base, &BigUint::from(scalar), arithmetic)
}

/// Normalize a batch with one inversion (Montgomery's trick).  Infinity
/// entries stay `None` and do not enter the product.
pub fn batch_to_affine(
    points: &[JacobianPoint],
    arithmetic: &mut P192Arithmetic,
) -> Vec<Option<AffinePoint>> {
    let mut prefixes = Vec::with_capacity(points.len());
    let mut product = arithmetic.one();
    for point in points {
        prefixes.push(product);
        if !point.is_infinity(arithmetic) {
            product = arithmetic.mul(product, point.z);
        }
    }
    if points.iter().all(|point| point.is_infinity(arithmetic)) {
        return vec![None; points.len()];
    }
    let mut inverse = arithmetic
        .inverse(product)
        .expect("product of nonzero P-192 coordinates is invertible");
    let mut out = vec![None; points.len()];
    for index in (0..points.len()).rev() {
        let point = points[index];
        if point.is_infinity(arithmetic) {
            continue;
        }
        let z_inv = arithmetic.mul(inverse, prefixes[index]);
        inverse = arithmetic.mul(inverse, point.z);
        let z2 = arithmetic.square(z_inv);
        let z3 = arithmetic.mul(z2, z_inv);
        out[index] = Some(AffinePoint {
            x: arithmetic.mul(point.x, z2),
            y: arithmetic.mul(point.y, z3),
        });
    }
    out
}

pub fn affine_equal(
    left: JacobianPoint,
    right: JacobianPoint,
    arithmetic: &mut P192Arithmetic,
) -> bool {
    if left.is_infinity(arithmetic) || right.is_infinity(arithmetic) {
        return left.is_infinity(arithmetic) && right.is_infinity(arithmetic);
    }
    let z1_2 = arithmetic.square(left.z);
    let z2_2 = arithmetic.square(right.z);
    let lx = arithmetic.mul(left.x, z2_2);
    let rx = arithmetic.mul(right.x, z1_2);
    if lx != rx {
        return false;
    }
    let z1_3 = arithmetic.mul(z1_2, left.z);
    let z2_3 = arithmetic.mul(z2_2, right.z);
    arithmetic.mul(left.y, z2_3) == arithmetic.mul(right.y, z1_3)
}

/// Homogeneous nodal parameter `t=u/v`.  `(1:0)` is the identity and
/// multiplication is `(u1u2-3v1v2 : u1v2+u2v1)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TorusPoint {
    pub u: Fp,
    pub v: Fp,
}

impl TorusPoint {
    pub fn identity(arithmetic: &P192Arithmetic) -> Self {
        Self {
            u: arithmetic.one(),
            v: arithmetic.zero(),
        }
    }

    pub fn from_t(t: u64, arithmetic: &P192Arithmetic) -> Self {
        Self {
            u: arithmetic.from_u64(t),
            v: arithmetic.one(),
        }
    }

    pub fn is_identity(&self, arithmetic: &P192Arithmetic) -> bool {
        self.v == arithmetic.zero()
    }

    pub fn compose(self, other: Self, arithmetic: &mut P192Arithmetic) -> Self {
        arithmetic.counts.point_additions += 1;
        let uu = arithmetic.mul(self.u, other.u);
        let vv = arithmetic.mul(self.v, other.v);
        let three_vv = arithmetic.mul_small(vv, 3);
        let u = arithmetic.sub(uu, three_vv);
        let uv = arithmetic.mul(self.u, other.v);
        let vu = arithmetic.mul(other.u, self.v);
        let v = arithmetic.add(uv, vu);
        Self { u, v }
    }

    pub fn square(self, arithmetic: &mut P192Arithmetic) -> Self {
        arithmetic.counts.point_doublings += 1;
        let u2 = arithmetic.square(self.u);
        let v2 = arithmetic.square(self.v);
        let three_v2 = arithmetic.mul_small(v2, 3);
        let u = arithmetic.sub(u2, three_v2);
        let uv = arithmetic.mul(self.u, self.v);
        let v = arithmetic.mul_small(uv, 2);
        Self { u, v }
    }

    pub fn negate(self, arithmetic: &mut P192Arithmetic) -> Self {
        Self {
            u: self.u,
            v: arithmetic.neg(self.v),
        }
    }

    pub fn scalar_mul(self, scalar: &BigUint, arithmetic: &mut P192Arithmetic) -> Self {
        arithmetic.counts.scalar_multiplications += 1;
        let mut acc = Self::identity(arithmetic);
        for bit in (0..scalar.bits()).rev() {
            acc = acc.square(arithmetic);
            if scalar.bit(bit) {
                acc = acc.compose(self, arithmetic);
            }
        }
        acc
    }

    pub fn equivalent(self, other: Self, arithmetic: &mut P192Arithmetic) -> bool {
        arithmetic.mul(self.u, other.v) == arithmetic.mul(other.u, self.v)
    }

    pub fn to_curve(self, arithmetic: &mut P192Arithmetic) -> Option<AffinePoint> {
        if self.is_identity(arithmetic) {
            return None;
        }
        let v_inv = arithmetic.inverse(self.v)?;
        let t = arithmetic.mul(self.u, v_inv);
        let t2 = arithmetic.square(t);
        let x = arithmetic.add(t2, arithmetic.from_u64(2));
        let x_plus_one = arithmetic.add(x, arithmetic.one());
        let y = arithmetic.mul(t, x_plus_one);
        Some(AffinePoint { x, y })
    }
}

/// Convert homogeneous torus parameters to singular-curve points with one
/// inversion for the whole batch.  Identity entries remain `None`.
pub fn batch_torus_to_curve(
    points: &[TorusPoint],
    arithmetic: &mut P192Arithmetic,
) -> Vec<Option<AffinePoint>> {
    let mut prefixes = Vec::with_capacity(points.len());
    let mut product = arithmetic.one();
    for point in points {
        prefixes.push(product);
        if !point.is_identity(arithmetic) {
            product = arithmetic.mul(product, point.v);
        }
    }
    if points.iter().all(|point| point.is_identity(arithmetic)) {
        return vec![None; points.len()];
    }
    let mut inverse = arithmetic
        .inverse(product)
        .expect("product of nonzero torus denominators is invertible");
    let mut out = vec![None; points.len()];
    for index in (0..points.len()).rev() {
        let point = points[index];
        if point.is_identity(arithmetic) {
            continue;
        }
        let v_inv = arithmetic.mul(inverse, prefixes[index]);
        inverse = arithmetic.mul(inverse, point.v);
        let t = arithmetic.mul(point.u, v_inv);
        let t2 = arithmetic.square(t);
        let x = arithmetic.add(t2, arithmetic.from_u64(2));
        let x_plus_one = arithmetic.add(x, arithmetic.one());
        let y = arithmetic.mul(t, x_plus_one);
        out[index] = Some(AffinePoint { x, y });
    }
    out
}

pub fn is_on_p192(point: AffinePoint, arithmetic: &mut P192Arithmetic) -> bool {
    let y2 = arithmetic.square(point.y);
    let x2 = arithmetic.square(point.x);
    let x3 = arithmetic.mul(x2, point.x);
    let three_x = arithmetic.mul_small(point.x, 3);
    let rhs0 = arithmetic.sub(x3, three_x);
    let b = arithmetic.from_hex(B_HEX);
    y2 == arithmetic.add(rhs0, b)
}

pub fn is_on_singular_companion(point: AffinePoint, arithmetic: &mut P192Arithmetic) -> bool {
    let y2 = arithmetic.square(point.y);
    let x2 = arithmetic.square(point.x);
    let x3 = arithmetic.mul(x2, point.x);
    let three_x = arithmetic.mul_small(point.x, 3);
    let x3_minus_three_x = arithmetic.sub(x3, three_x);
    let rhs = arithmetic.sub(x3_minus_three_x, arithmetic.from_u64(2));
    y2 == rhs
}

/// A stable 64-bit fingerprint of a canonical x-coordinate.  `BabyTable`
/// hashes this value again and every hit is replayed exactly.
pub fn x_fingerprint(point: AffinePoint, arithmetic: &P192Arithmetic) -> u64 {
    let limbs = arithmetic.canonical(point.x).0;
    let mut h = 0x6a09_e667_f3bc_c909u64;
    for limb in limbs {
        h ^= limb.wrapping_add(0x9e37_79b9_7f4a_7c15);
        h ^= h >> 30;
        h = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
        h ^= h >> 27;
        h = h.wrapping_mul(0x94d0_49bb_1331_11eb);
        h ^= h >> 31;
    }
    h
}

pub fn fixed_width_biguint(value: &BigUint) -> Result<[u8; 24], String> {
    let bytes = value.to_bytes_be();
    if bytes.len() > 24 {
        return Err("integer does not fit a P-192 coordinate".into());
    }
    let mut out = [0u8; 24];
    out[24 - bytes.len()..].copy_from_slice(&bytes);
    Ok(out)
}

pub fn singular_factor_product() -> BigUint {
    (BigUint::one() << 64usize)
        * BigUint::from(3u8)
        * BigUint::from(5u8)
        * BigUint::from(17u8)
        * BigUint::from(257u16)
        * BigUint::from(641u16)
        * BigUint::from(65_537u32)
        * BigUint::from(274_177u32)
        * BigUint::from(6_700_417u32)
        * BigUint::from(QMAX)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;

    #[test]
    fn montgomery_field_matches_biguint() {
        let mut arithmetic = P192Arithmetic::new();
        let p = p192_prime();
        let a = parse_hex("123456789abcdef00112233445566778899aabbccddeeff");
        let b = parse_hex("fedcba98765432100123456789abcdef0123456789abcd");
        let fa = arithmetic.from_biguint(&a).unwrap();
        let fb = arithmetic.from_biguint(&b).unwrap();
        let sum = arithmetic.add(fa, fb);
        let product = arithmetic.mul(fa, fb);
        let inv = arithmetic.inverse(fa).unwrap();
        assert_eq!(arithmetic.to_biguint(sum), (&a + &b) % &p);
        assert_eq!(arithmetic.to_biguint(product), (&a * &b) % &p);
        assert_eq!(arithmetic.mul(fa, inv), arithmetic.one());
    }

    #[test]
    fn jacobian_generator_order_and_repo_agree() {
        let mut arithmetic = P192Arithmetic::new();
        let curve = CurveParams::p192();
        let g = AffinePoint::generator(&arithmetic);
        assert!(is_on_p192(g, &mut arithmetic));
        let two = scalar_mul_u64(g, 2, &mut arithmetic)
            .to_affine(&mut arithmetic)
            .unwrap();
        let repo_two = curve
            .generator()
            .scalar_mul(&BigUint::from(2u8), &curve.a_fe());
        assert_eq!(
            Some(two),
            AffinePoint::from_repo_point(&repo_two, &arithmetic).unwrap()
        );
        assert!(scalar_mul(g, &p192_order(), &mut arithmetic).is_infinity(&arithmetic));
    }

    #[test]
    fn torus_parameterization_is_independent_and_exact() {
        let mut arithmetic = P192Arithmetic::new();
        let r = TorusPoint::from_t(8, &arithmetic);
        let curve_r = r.to_curve(&mut arithmetic).unwrap();
        assert_eq!(curve_r, AffinePoint::singular_generator(&arithmetic));
        assert!(is_on_singular_companion(curve_r, &mut arithmetic));
        assert!(!is_on_p192(curve_r, &mut arithmetic));
        let doubled = r.square(&mut arithmetic);
        let composed = r.compose(r, &mut arithmetic);
        assert!(doubled.equivalent(composed, &mut arithmetic));
        assert!(r
            .scalar_mul(&singular_order(), &mut arithmetic)
            .is_identity(&arithmetic));
        for factor in SINGULAR_DISTINCT_FACTORS {
            let exponent = singular_order() / BigUint::from(factor);
            assert!(!r
                .scalar_mul(&exponent, &mut arithmetic)
                .is_identity(&arithmetic));
        }
    }

    #[test]
    fn batch_normalization_matches_individual_normalization() {
        let mut arithmetic = P192Arithmetic::new();
        let g = AffinePoint::generator(&arithmetic);
        let points: Vec<_> = (0u64..32)
            .map(|i| scalar_mul_u64(g, i, &mut arithmetic))
            .collect();
        let individual: Vec<_> = points
            .iter()
            .copied()
            .map(|point| point.to_affine(&mut arithmetic))
            .collect();
        let batch = batch_to_affine(&points, &mut arithmetic);
        assert_eq!(batch, individual);
    }

    #[test]
    fn frozen_singular_factorization_matches() {
        assert_eq!(singular_factor_product(), singular_order());
        assert_eq!(residue_modulus() * BigUint::from(QMAX), singular_order());
        assert_eq!(p192_order(), parse_hex(N_HEX));
        assert!(!p192_prime().is_zero());
    }
}
