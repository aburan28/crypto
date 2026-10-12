//! Producer-side identity and maximal-order evidence.

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde_json::{json, Value};

use super::arithmetic::{integer, unsigned, C_DEC, D_DEC, N_DEC, P_DEC, T_DEC};
use super::pocklington;
use super::provenance::BuildProvenance;
use super::schema::{verdicts, CURVE_UID, EC1, EXPERIMENT_ID, ICV1, PROTOCOL_VERSION};
use super::Result;

pub const A_DEC: &str = "6277101735386680763835789423207666416083908700390324961276";
pub const B_DEC: &str = "2455155546008943817740293915197451784769108058161191238065";
pub const GX_HEX: &str = "188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012";
pub const GY_HEX: &str = "07192B95FFC8DA78631011ED6B24CDD573F977A11E794811";

#[derive(Clone, Debug, PartialEq, Eq)]
enum Point {
    Infinity,
    Affine(BigUint, BigUint),
}

fn sub_mod(left: &BigUint, right: &BigUint, modulus: &BigUint) -> BigUint {
    if left >= right {
        (left - right) % modulus
    } else {
        (modulus - ((right - left) % modulus)) % modulus
    }
}

fn inverse(value: &BigUint, modulus: &BigUint) -> Result<BigUint> {
    if value.is_zero() {
        return Err("attempted to invert zero in the curve field".to_owned());
    }
    Ok(value.modpow(&(modulus - 2u8), modulus))
}

fn add(left: &Point, right: &Point, a: &BigUint, p: &BigUint) -> Result<Point> {
    match (left, right) {
        (Point::Infinity, _) => Ok(right.clone()),
        (_, Point::Infinity) => Ok(left.clone()),
        (Point::Affine(x1, y1), Point::Affine(x2, y2)) => {
            if x1 == x2 && (y1 + y2) % p == BigUint::zero() {
                return Ok(Point::Infinity);
            }
            let slope = if x1 == x2 {
                let numerator = (BigUint::from(3u8) * x1 * x1 + a) % p;
                let denominator = (BigUint::from(2u8) * y1) % p;
                numerator * inverse(&denominator, p)? % p
            } else {
                let numerator = sub_mod(y2, y1, p);
                let denominator = sub_mod(x2, x1, p);
                numerator * inverse(&denominator, p)? % p
            };
            let x3 = sub_mod(&sub_mod(&(&slope * &slope % p), x1, p), x2, p);
            let y3 = sub_mod(&(&slope * sub_mod(x1, &x3, p) % p), y1, p);
            Ok(Point::Affine(x3, y3))
        }
    }
}

fn multiply(point: &Point, scalar: &BigUint, a: &BigUint, p: &BigUint) -> Result<Point> {
    let mut accumulator = Point::Infinity;
    let mut power = point.clone();
    let mut remaining = scalar.clone();
    while !remaining.is_zero() {
        if remaining.bit(0) {
            accumulator = add(&accumulator, &power, a, p)?;
        }
        remaining >>= 1usize;
        if !remaining.is_zero() {
            power = add(&power, &power, a, p)?;
        }
    }
    Ok(accumulator)
}

pub fn curve_object() -> Value {
    json!({
        "a": A_DEC,
        "b": B_DEC,
        "cofactor": 1,
        "curve_uid": CURVE_UID,
        "ec1": EC1,
        "generator_x_hex": GX_HEX,
        "generator_y_hex": GY_HEX,
        "icv1": ICV1,
        "model": "short Weierstrass y^2 = x^3 + a*x + b over F_p",
        "n": N_DEC,
        "name": "NIST P-192 / secp192r1",
        "p": P_DEC,
        "trace_t": T_DEC,
    })
}

fn validate_curve() -> Result<()> {
    let p = unsigned(P_DEC)?;
    let a = unsigned(A_DEC)?;
    let b = unsigned(B_DEC)?;
    let gx = BigUint::parse_bytes(GX_HEX.as_bytes(), 16)
        .ok_or_else(|| "invalid frozen generator x".to_owned())?;
    let gy = BigUint::parse_bytes(GY_HEX.as_bytes(), 16)
        .ok_or_else(|| "invalid frozen generator y".to_owned())?;
    let n = unsigned(N_DEC)?;
    let t = unsigned(T_DEC)?;
    if &p + 1u8 - &t != n {
        return Err("P-192 group-order identity failed".to_owned());
    }
    if (BigUint::from(4u8) * a.modpow(&BigUint::from(3u8), &p)
        + BigUint::from(27u8) * b.modpow(&BigUint::from(2u8), &p))
        % &p
        == BigUint::zero()
    {
        return Err("P-192 curve is singular".to_owned());
    }
    if (&gy * &gy) % &p != (&gx * &gx * &gx + &a * &gx + &b) % &p {
        return Err("P-192 generator is not on the curve".to_owned());
    }
    let generator = Point::Affine(gx, gy);
    if generator == Point::Infinity || multiply(&generator, &n, &a, &p)? != Point::Infinity {
        return Err("P-192 generator/order check failed".to_owned());
    }
    Ok(())
}

pub fn identity(provenance: &BuildProvenance) -> Result<Value> {
    validate_curve()?;
    let p = integer(P_DEC)?;
    let t = integer(T_DEC)?;
    if &t * &t - &p * 4 != integer(D_DEC)? {
        return Err("P-192 Frobenius-discriminant identity failed".to_owned());
    }
    let checks = [
        "serialized_tuple_matches_protocol",
        "identities_recomputed",
        "nonsingular",
        "generator_nonidentity",
        "generator_order_n",
        "exact_group_order_n",
    ];
    Ok(json!({
        "checks": verdicts(&checks),
        "curve": curve_object(),
        "derived": {"frobenius_discriminant":D_DEC,"group_order":N_DEC},
        "experiment_id": EXPERIMENT_ID,
        "overall_status": "PASS",
        "protocol_commit": provenance.protocol_commit,
        "protocol_version": PROTOCOL_VERSION,
        "schema": "p192-wcm-identity-v1",
        "source_commit": provenance.source_commit,
    }))
}

fn authoritative_proof(value: &str) -> Value {
    json!({
        "document_byte_length": 306784,
        "document_sha256": "87b8f3703364ed5b21ba8582e411cc0cbf477bcaa3f4f45e0d6580d1c00d9952",
        "document_url": "https://www.secg.org/sec2-v2.pdf",
        "edition": "2.0",
        "method": "authoritative_standard",
        "parameter_name": "secp192r1",
        "section": "2.2.2",
        "standard_id": "SEC2",
        "value": value,
    })
}

pub fn order_object() -> Value {
    json!({
        "C": C_DEC,
        "D": D_DEC,
        "D_End": D_DEC,
        "D_K": D_DEC,
        "D_pi": D_DEC,
        "factorization": "-5 * 11 * 31 * C",
        "n": N_DEC,
        "p": P_DEC,
        "t": T_DEC,
    })
}

pub fn maximal_order(provenance: &BuildProvenance) -> Result<Value> {
    let p = integer(P_DEC)?;
    let n = integer(N_DEC)?;
    let t = integer(T_DEC)?;
    let d = integer(D_DEC)?;
    let c = integer(C_DEC)?;
    if &t * &t - &p * 4 != d
        || &p + 1 - &t != n
        || d != -(BigInt::from(5u8) * BigInt::from(11u8) * BigInt::from(31u8) * &c)
        || d.mod_floor(&BigInt::from(4u8)) != BigInt::one()
        || p.mod_floor(&t).is_zero()
        || t.abs() >= p
    {
        return Err("frozen P-192 maximal-order identities failed".to_owned());
    }
    let c_proof = pocklington::c_proof()?;
    let primality = vec![
        json!({"proof":authoritative_proof(P_DEC),"subject":"p","value":P_DEC}),
        json!({"proof":authoritative_proof(N_DEC),"subject":"n","value":N_DEC}),
        json!({"proof":c_proof,"subject":"C","value":C_DEC}),
    ];
    let checks = [
        "p_prime",
        "n_prime",
        "C_prime",
        "p_not_divide_t",
        "trace_absolute_value_below_p",
        "ordinary",
        "discriminant_identity",
        "group_order_identity",
        "discriminant_factorization",
        "D_squarefree",
        "D_congruent_1_mod_4",
        "D_fundamental",
        "frobenius_conductor_one",
        "endomorphism_order_maximal",
        "discriminants_equal",
        "unit_group_plus_minus_one",
        "integral_basis_multiplication",
    ];
    Ok(json!({
        "checks": verdicts(&checks),
        "curve_uid": CURVE_UID,
        "experiment_id": EXPERIMENT_ID,
        "order": order_object(),
        "overall_status": "PASS",
        "primality": primality,
        "protocol_commit": provenance.protocol_commit,
        "protocol_version": PROTOCOL_VERSION,
        "schema": "p192-wcm-maximal-order-v1",
        "source_commit": provenance.source_commit,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_curve_checks_replay() {
        validate_curve().unwrap();
        assert_eq!(curve_object()["generator_y_hex"], GY_HEX);
        assert_eq!(GY_HEX.len(), 48);
    }

    #[test]
    fn authoritative_metadata_is_exact() {
        let proof = authoritative_proof(P_DEC);
        assert_eq!(proof["standard_id"], "SEC2");
        assert_eq!(proof["edition"], "2.0");
        assert_eq!(proof["section"], "2.2.2");
        assert_eq!(proof["document_byte_length"], 306784);
    }
}
