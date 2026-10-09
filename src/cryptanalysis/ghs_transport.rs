//! Checked characteristic-two trace transport over composite-degree fields.
//!
//! This is the executable `m = 1` branch of the binary GHS route.  It
//! verifies the field, points, supplied subgroup annihilator, and the
//! Frobenius-fixed curve model before applying the trace.  A nonzero image
//! preserves a DLP modulo a *certified prime* subgroup order; certification
//! of that supplied order and the relation `Q = [d]P` are external inputs.
//! Higher-magic GHS covers need a smooth model and a norm-conorm map; this
//! module does not claim that the trace computes either one.

use crate::binary_ecc::F2mElement;
use crate::cryptanalysis::ec_trapdoor::FieldTower;
use crate::cryptanalysis::ghs_descent::{trace_map, ECurve, Pt};
use crate::cryptanalysis::ghs_screen::{validate_input, GhsCurveInput};
use num_bigint::BigUint;
use num_traits::One;

/// Affine coordinates are canonical polynomial-basis bitsets.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BinaryPointInput {
    pub x: BigUint,
    pub y: BigUint,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TraceStatus {
    /// The specified Frobenius does not fix both curve coefficients.
    CurveNotDefinedOverSubfield,
    /// The trace killed the source generator, so this map loses its DLP.
    GeneratorKilledByTrace,
    /// Point transport is executable and its image has the supplied
    /// annihilator. If that annihilator is certified prime, the nonzero
    /// image has the same prime order as the source generator.
    NonzeroSubgroupImage,
}

#[derive(Clone, Debug)]
pub struct TraceTransport {
    pub relative_degree: u32,
    pub base_degree: u32,
    pub status: TraceStatus,
    pub generator_image: Option<Pt>,
    pub target_image: Option<Pt>,
}

fn decode_point(input: &BinaryPointInput, m: u32) -> Result<Pt, String> {
    if input.x.bits() > u64::from(m) || input.y.bits() > u64::from(m) {
        return Err("point coordinates are outside the binary field".into());
    }
    Ok(Pt::Aff {
        x: F2mElement::from_biguint(&input.x, m),
        y: F2mElement::from_biguint(&input.y, m),
    })
}

/// Execute the trace branch for `F_(2^N)/F_(2^l)` with `N = n*l`.
///
/// `order` is a claimed subgroup order, not a primality certificate. The
/// function checks `[order]P = [order]Q = O` and, for a nonzero trace image,
/// `[order]Tr(P) = [order]Tr(Q) = O`. For a certified prime `order` and a
/// valid DLP input `Q = [d]P`, a nonzero image preserves `d mod order`.
pub fn trace_transport(
    curve_input: &GhsCurveInput,
    relative_degree: u32,
    generator: &BinaryPointInput,
    target: &BinaryPointInput,
    order: &BigUint,
) -> Result<TraceTransport, String> {
    let irr = validate_input(curve_input).map_err(|error| error.to_string())?;
    let m = curve_input.absolute_degree;
    if relative_degree < 2 || !m.is_multiple_of(relative_degree) {
        return Err("relative degree must be a nontrivial divisor of the field degree".into());
    }
    if order <= &BigUint::one() || order.bits() > u64::from(m) + 2 {
        return Err("subgroup annihilator is outside the expected curve-order range".into());
    }
    let base_degree = m / relative_degree;
    let tower = FieldTower::new(m, relative_degree, base_degree, irr.clone());
    let a = F2mElement::from_biguint(&curve_input.a, m);
    let b = F2mElement::from_biguint(&curve_input.b, m);
    let curve = ECurve::new(m, irr, a.clone(), b.clone());
    let p = decode_point(generator, m)?;
    let q = decode_point(target, m)?;
    if !curve.is_on_curve(&p) || !curve.is_on_curve(&q) {
        return Err("a supplied point is not on the curve".into());
    }
    if curve.scalar_mul(&p, order) != Pt::Inf || curve.scalar_mul(&q, order) != Pt::Inf {
        return Err("the supplied annihilator does not kill both points".into());
    }
    if !tower.is_in_subfield(&a, base_degree) || !tower.is_in_subfield(&b, base_degree) {
        return Ok(TraceTransport {
            relative_degree,
            base_degree,
            status: TraceStatus::CurveNotDefinedOverSubfield,
            generator_image: None,
            target_image: None,
        });
    }
    let image_p = trace_map(&curve, &tower, &p);
    let image_q = trace_map(&curve, &tower, &q);
    if !curve.is_on_curve(&image_p) || !curve.is_on_curve(&image_q) {
        return Err("trace produced a point off the original curve".into());
    }
    if curve.scalar_mul(&image_p, order) != Pt::Inf || curve.scalar_mul(&image_q, order) != Pt::Inf
    {
        return Err("trace image failed the subgroup annihilator check".into());
    }
    let status = if image_p == Pt::Inf {
        TraceStatus::GeneratorKilledByTrace
    } else {
        TraceStatus::NonzeroSubgroupImage
    };
    Ok(TraceTransport {
        relative_degree,
        base_degree,
        status,
        generator_image: Some(image_p),
        target_image: Some(image_q),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn curve(a: u32) -> GhsCurveInput {
        GhsCurveInput {
            absolute_degree: 6,
            modulus: BigUint::from(0x43u32), // z^6 + z + 1
            a: BigUint::from(a),
            b: BigUint::one(),
        }
    }

    fn two_torsion() -> BinaryPointInput {
        BinaryPointInput {
            x: BigUint::from(0u8),
            y: BigUint::one(),
        }
    }

    #[test]
    fn transports_a_binary_subgroup_through_a_composite_tower() {
        let p = two_torsion();
        let result = trace_transport(&curve(0), 3, &p, &p, &BigUint::from(2u8)).unwrap();
        assert_eq!(result.status, TraceStatus::NonzeroSubgroupImage);
        assert_eq!(result.base_degree, 2);
        assert_eq!(result.generator_image, result.target_image);
        assert!(matches!(result.generator_image, Some(Pt::Aff { .. })));
    }

    #[test]
    fn reports_frobenius_and_trace_limits() {
        let p = two_torsion();
        let not_fixed = trace_transport(&curve(2), 3, &p, &p, &BigUint::from(2u8)).unwrap();
        assert_eq!(not_fixed.status, TraceStatus::CurveNotDefinedOverSubfield);
        let killed = trace_transport(&curve(0), 2, &p, &p, &BigUint::from(2u8)).unwrap();
        assert_eq!(killed.status, TraceStatus::GeneratorKilledByTrace);
        assert_eq!(killed.generator_image, Some(Pt::Inf));
    }

    #[test]
    fn rejects_invalid_points_and_annihilators() {
        let p = two_torsion();
        let off_curve = BinaryPointInput {
            x: BigUint::from(0u8),
            y: BigUint::from(0u8),
        };
        assert!(trace_transport(&curve(0), 3, &p, &off_curve, &BigUint::from(2u8)).is_err());
        assert!(trace_transport(&curve(0), 3, &p, &p, &BigUint::from(3u8)).is_err());
        assert!(trace_transport(&curve(0), 4, &p, &p, &BigUint::from(2u8)).is_err());
    }
}
