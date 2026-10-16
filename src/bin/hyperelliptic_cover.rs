//! Standalone, structured interface to the repository's explicit
//! hyperelliptic-cover certificates and checked Jacobian arithmetic.
//!
//! This is deliberately not a discrete-log solver.  The prime-field cover
//! emitted by the constructor has two points at infinity, so the current
//! one-infinity odd-characteristic Jacobian implementation is not silently
//! applied to it.  Binary cover arithmetic is supported because that family
//! has one rational point at infinity and a checked characteristic-two
//! Mumford implementation.

#[path = "curve_cover_check/checker.rs"]
mod checker;
#[path = "curve_cover_check/models.rs"]
mod models;

use clap::{Parser, Subcommand};
use crypto_lib::binary_ecc::{
    F2mElement, F2mPoly, HyperellipticCurve, MumfordDivisor, OrdinaryBinaryCover,
};
use crypto_lib::prime_hyperelliptic::{FpPoly, HyperellipticCurveP, MumfordDivisorP};
use num_bigint::BigUint;
use num_traits::Zero;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
};

const RESPONSE_SCHEMA: &str = "hyperelliptic-infrastructure-command/v1";
const ARITHMETIC_SCHEMA: &str = "hyperelliptic-jacobian-operation/v1";
const MAX_INPUT_BYTES: u64 = 32 * 1024 * 1024;

#[derive(Parser)]
#[command(
    name = "hyperelliptic-cover",
    about = "Construct and verify explicit H -> E covers and run checked Jacobian arithmetic",
    long_about = "A bounded, certificate-oriented interface. It constructs the repository's fixed same-field cover families, replays their certificates, runs checked arithmetic where the infinity model is supported, and exports catalog evidence. It is not an index-calculus or key-recovery solver."
)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Construct and replay a cover certificate for one exact elliptic model object.
    Construct {
        /// Exact supported elliptic model shape, or `-` for stdin; no registry lookup is implied.
        #[arg(long)]
        input: PathBuf,
        /// Response path; omit for stdout.
        #[arg(long)]
        output: Option<PathBuf>,
    },
    /// Replay a supplied `{model, certificate}` request.
    Verify {
        /// Verification request JSON, or `-` for stdin.
        #[arg(long)]
        input: PathBuf,
        /// Response path; omit for stdout.
        #[arg(long)]
        output: Option<PathBuf>,
    },
    /// Run checked arithmetic on supported binary-cover or odd-prime one-infinity Jacobians.
    Arithmetic {
        /// Arithmetic request JSON, or `-` for stdin.
        #[arg(long)]
        input: PathBuf,
        /// Response path; omit for stdout.
        #[arg(long)]
        output: Option<PathBuf>,
    },
    /// Rebuild the deterministic cover catalog from the exact model registry.
    CatalogExport {
        /// Exact EC1/ICV1 registry JSON (required; release archives do not bundle it).
        #[arg(long)]
        registry: PathBuf,
        /// Catalog path; omit to write the complete catalog to stdout.
        #[arg(long)]
        output: Option<PathBuf>,
    },
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct VerifyRequest {
    model: Value,
    certificate: checker::Certificate,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct DivisorWire {
    /// Coefficients in ascending powers, as canonical unsigned integers.
    u: Vec<String>,
    /// Coefficients in ascending powers, as canonical unsigned integers.
    v: Vec<String>,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(rename_all = "snake_case")]
enum ArithmeticOperation {
    Validate,
    Identity,
    Equal,
    Negate,
    Add,
    ScalarMultiply,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ArithmeticRequest {
    schema: String,
    /// Exact target elliptic model.  The source H is reconstructed from the
    /// verified fixed-family certificate, rather than accepted by assertion.
    #[serde(default)]
    model: Option<Value>,
    /// A separate odd-characteristic, odd-degree one-infinity curve.  This is
    /// intentionally not conflated with the even-sextic prime cover.
    #[serde(default)]
    odd_curve: Option<OddCurveWire>,
    operation: ArithmeticOperation,
    #[serde(default)]
    left: Option<DivisorWire>,
    #[serde(default)]
    right: Option<DivisorWire>,
    #[serde(default)]
    scalar: Option<String>,
}

#[derive(Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct OddCurveWire {
    family: String,
    p: String,
    genus: u32,
    /// `f(x)` coefficients in ascending powers for `y^2=f(x)`.
    f: Vec<String>,
    infinity: String,
}

fn read_bytes(path: &Path) -> Result<Vec<u8>, String> {
    let (reader, source): (Box<dyn Read>, String) = if path.as_os_str() == "-" {
        (Box::new(std::io::stdin()), "stdin".into())
    } else {
        let file = fs::File::open(path).map_err(|e| format!("read {}: {e}", path.display()))?;
        (Box::new(file), path.display().to_string())
    };
    let mut bytes = Vec::new();
    reader
        .take(MAX_INPUT_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| format!("read {source}: {e}"))?;
    if bytes.len() as u64 > MAX_INPUT_BYTES {
        return Err(format!("input exceeds the {MAX_INPUT_BYTES}-byte limit"));
    }
    Ok(bytes)
}

fn parse_json(path: &Path) -> Result<Value, String> {
    serde_json::from_slice(&read_bytes(path)?).map_err(|e| format!("parse JSON: {e}"))
}

fn response(command: &str, status: &str, body: Value) -> Value {
    json!({
        "schema": RESPONSE_SCHEMA,
        "command": command,
        "status": status,
        "result": body,
    })
}

fn model_error(command: &str, error: checker::ModelError) -> (Value, i32) {
    match error {
        checker::ModelError::Unsupported(reason) => (
            response(
                command,
                "unsupported",
                json!({"reason": reason, "existence": "not_established_by_this_command"}),
            ),
            2,
        ),
        checker::ModelError::Invalid(reason) => (
            response(command, "invalid_input", json!({"reason": reason})),
            3,
        ),
    }
}

fn strict_model(value: &Value) -> Result<checker::Model, checker::ModelError> {
    let object = value.as_object().ok_or_else(|| {
        checker::ModelError::Invalid("elliptic model must be a JSON object".into())
    })?;
    let allowed: &[&str] = match value.get("form").and_then(Value::as_str) {
        Some("y^2=x^3+a*x+b") => &["v", "form", "p", "field", "a", "b"],
        Some("y^2+xy=x^3+a*x^2+b") => &["v", "form", "modulus", "field", "a", "b"],
        _ => &["v", "form", "p", "modulus", "field", "a", "b"],
    };
    if let Some(unexpected) = object.keys().find(|name| !allowed.contains(&name.as_str())) {
        return Err(checker::ModelError::Invalid(format!(
            "unexpected elliptic model field `{unexpected}`"
        )));
    }
    checker::model(value)
}

fn proof_receipt(model: &checker::Model) -> Result<Value, String> {
    match &model.field {
        checker::Field::Prime(_) => Ok(json!({
            "construction": "explicit",
            "field_of_definition": "same_as_target",
            "field_validation": "fixed-base_probable_prime_screen",
            "target_smoothness": "verified_by_discriminant",
            "source_smoothness": "verified_by_squarefree_sextic",
            "genus": "verified_genus_2",
            "defining_equation_identity": "verified_exact_polynomial_identity",
            "degree": "verified_degree_2",
            "separability": "verified_characteristic_not_2",
            "infinity": "two_base_field_rational_points",
            "jacobian_arithmetic": "unsupported_even_degree_two_infinity_model",
            "pullback": "not_implemented_for_this_cover_family",
            "pushforward": "not_implemented_for_this_cover_family",
            "subgroup_transfer": "unknown",
        })),
        // The same sextic cover over `GF(p^k)`; the checker also proves the
        // modulus irreducible over `GF(p)` by Rabin's test.
        checker::Field::Extension { .. } => Ok(json!({
            "construction": "explicit",
            "field_of_definition": "same_as_target",
            "field_validation": "fixed-base_probable_prime_screen_of_p_and_rabin_irreducible_modulus",
            "target_smoothness": "verified_by_discriminant",
            "source_smoothness": "verified_by_squarefree_sextic",
            "genus": "verified_genus_2",
            "defining_equation_identity": "verified_exact_polynomial_identity",
            "degree": "verified_degree_2",
            "separability": "verified_characteristic_not_2",
            "infinity": "two_base_field_rational_points",
            "jacobian_arithmetic": "unsupported_even_degree_two_infinity_model",
            "pullback": "not_implemented_for_this_cover_family",
            "pushforward": "not_implemented_for_this_cover_family",
            "subgroup_transfer": "unknown",
        })),
        checker::Field::Binary(irreducible) => {
            let cover = OrdinaryBinaryCover::new(
                irreducible.degree,
                irreducible.clone(),
                F2mElement::from_biguint(&model.a, irreducible.degree),
                F2mElement::from_biguint(&model.b, irreducible.degree),
            )
            .map_err(|error| format!("binary cover runtime validation failed: {error}"))?;
            cover
                .verify_certificate()
                .map_err(|error| format!("binary cover certificate replay failed: {error}"))?;
            let certificate = cover.certificate();
            Ok(json!({
                "construction": "explicit",
                "field_of_definition": "same_as_target",
                "field_validation": "irreducibility_proved_by_rabin",
                "target_smoothness": "verified_b_nonzero",
                "source_smoothness": "verified_affine_and_infinity_fixed_family_criteria",
                "genus": "verified_genus_3_from_two_order_3_artin_schreier_poles",
                "defining_equation_identity": "verified_exact_polynomial_identity",
                "degree": "verified_degree_3",
                "separability": "verified_characteristic_2_does_not_divide_3",
                "infinity": "one_base_field_rational_point",
                "ramification": {
                    "infinity_index": certificate.infinity_ramification_index,
                    "exceptional_index": certificate.exceptional_ramification_index,
                    "riemann_hurwitz_total": certificate.riemann_hurwitz_ramification_total,
                },
                "jacobian_arithmetic": "implemented_checked_one_infinity_model",
                "pullback": "implemented_for_rational_target_point_classes",
                "pushforward": "implemented_for_rational_source_point_divisors_and_typed_pullback_witnesses",
                "composition": "verified_phi_push_phi_pull_equals_multiplication_by_3_on_typed_rational_point_class_domain",
                "arbitrary_jacobian_pushforward": "unsupported",
                "subgroup_transfer": "conditional_on_an_independently_verified_exact_subgroup_order_coprime_to_3; this command does not accept or verify an order",
            }))
        }
    }
}

fn construct(model_value: Value) -> (Value, i32) {
    let model = match strict_model(&model_value) {
        Ok(model) => model,
        Err(error) => return model_error("construct", error),
    };
    let certificate = checker::construct(&model);
    if let Err(reason) = checker::verify(&model, &certificate) {
        return (
            response(
                "construct",
                "error",
                json!({"reason": format!("internal certificate replay failed: {reason}")}),
            ),
            1,
        );
    }
    let verification = match proof_receipt(&model) {
        Ok(verification) => verification,
        Err(reason) => return (response("construct", "error", json!({"reason": reason})), 1),
    };
    (
        response(
            "construct",
            "ok",
            json!({
                "model": model_value,
                "certificate": certificate,
                "verification": verification,
                "scope": "explicit same-field cover; no minimality, unconditional subgroup-transfer, or DLP-advantage claim",
            }),
        ),
        0,
    )
}

fn verify(request: VerifyRequest) -> (Value, i32) {
    let model = match strict_model(&request.model) {
        Ok(model) => model,
        Err(error) => return model_error("verify", error),
    };
    match checker::verify(&model, &request.certificate) {
        Ok(()) => match proof_receipt(&model) {
            Ok(verification) => (
                response(
                    "verify",
                    "ok",
                    json!({
                        "certificate": request.certificate,
                        "verification": verification,
                    }),
                ),
                0,
            ),
            Err(reason) => (response("verify", "error", json!({"reason": reason})), 1),
        },
        Err(reason) => (
            response("verify", "invalid_input", json!({"reason": reason})),
            3,
        ),
    }
}

fn parse_number(text: &str) -> Result<BigUint, String> {
    if text.len() > 1300 {
        return Err("integer exceeds field-size limit".into());
    }
    let (digits, radix) = text
        .strip_prefix("0x")
        .map_or((text, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("invalid unsigned integer: {text}"))
}

fn binary_poly(encoded: &[String], field: &checker::Field, m: u32) -> Result<F2mPoly, String> {
    let coefficients = checker::decoded(encoded, field)?
        .into_iter()
        .map(|coefficient| F2mElement::from_biguint(&coefficient, m))
        .collect();
    Ok(F2mPoly::from_coeffs(coefficients, m))
}

fn divisor_from_wire(
    wire: &DivisorWire,
    field: &checker::Field,
    curve: &HyperellipticCurve,
) -> Result<MumfordDivisor, String> {
    let u = binary_poly(&wire.u, field, curve.m)?;
    let v = binary_poly(&wire.v, field, curve.m)?;
    MumfordDivisor::try_new(curve, u, v).map_err(|e| format!("{e:?}"))
}

fn encoded_poly(poly: &F2mPoly) -> Vec<String> {
    poly.coeffs
        .iter()
        .map(|coefficient| format!("0x{:x}", coefficient.to_biguint()))
        .collect()
}

fn divisor_to_wire(divisor: &MumfordDivisor) -> DivisorWire {
    DivisorWire {
        u: encoded_poly(&divisor.u),
        v: encoded_poly(&divisor.v),
    }
}

fn required_divisor<'a>(
    value: &'a Option<DivisorWire>,
    name: &str,
) -> Result<&'a DivisorWire, String> {
    value
        .as_ref()
        .ok_or_else(|| format!("operation requires `{name}` divisor"))
}

fn operation_name(operation: ArithmeticOperation) -> &'static str {
    match operation {
        ArithmeticOperation::Validate => "validate",
        ArithmeticOperation::Identity => "identity",
        ArithmeticOperation::Equal => "equal",
        ArithmeticOperation::Negate => "negate",
        ArithmeticOperation::Add => "add",
        ArithmeticOperation::ScalarMultiply => "scalar_multiply",
    }
}

fn arithmetic_result(
    operation: ArithmeticOperation,
    context: Value,
    computed: Result<Value, String>,
) -> (Value, i32) {
    let operation = operation_name(operation);
    match computed {
        Ok(result) => (
            response(
                "arithmetic",
                "ok",
                json!({
                    "operation": operation,
                    "context": context,
                    "result": result,
                }),
            ),
            0,
        ),
        Err(reason) => (
            response(
                "arithmetic",
                "invalid_input",
                json!({"operation": operation, "reason": reason}),
            ),
            3,
        ),
    }
}

fn binary_arithmetic(request: &ArithmeticRequest, model_value: &Value) -> (Value, i32) {
    let model = match strict_model(model_value) {
        Ok(model) => model,
        Err(error) => return model_error("arithmetic", error),
    };
    let certificate = checker::construct(&model);
    if let Err(reason) = checker::verify(&model, &certificate) {
        return (
            response("arithmetic", "error", json!({"reason": reason})),
            1,
        );
    }
    let irr = match &model.field {
        checker::Field::Extension { .. } => {
            return (
                response(
                    "arithmetic",
                    "unsupported",
                    json!({
                        "reason": "the explicit GF(p^k) cover is an even-degree sextic with two rational points at infinity; the checked odd-characteristic implementation currently supports only odd-degree one-infinity models over prime fields",
                        "construction_status": "explicit_verified",
                        "arithmetic_status": "unsupported_infinity_configuration",
                    }),
                ),
                2,
            )
        }
        checker::Field::Prime(_) => {
            return (
                response(
                    "arithmetic",
                    "unsupported",
                    json!({
                        "reason": "the explicit prime cover is an even-degree sextic with two rational points at infinity; the checked odd-characteristic implementation currently supports only odd-degree one-infinity models",
                        "construction_status": "explicit_verified",
                        "arithmetic_status": "unsupported_infinity_configuration",
                    }),
                ),
                2,
            )
        }
        checker::Field::Binary(irr) => irr,
    };
    let h = match binary_poly(&certificate.h, &model.field, irr.degree) {
        Ok(poly) => poly,
        Err(reason) => {
            return (
                response("arithmetic", "error", json!({"reason": reason})),
                1,
            )
        }
    };
    let f = match binary_poly(&certificate.f, &model.field, irr.degree) {
        Ok(poly) => poly,
        Err(reason) => {
            return (
                response("arithmetic", "error", json!({"reason": reason})),
                1,
            )
        }
    };
    let curve = match HyperellipticCurve::try_new(irr.degree, irr.clone(), h, f, certificate.genus)
    {
        Ok(curve) => curve,
        Err(error) => {
            return (
                response(
                    "arithmetic",
                    "error",
                    json!({"reason": format!("constructed curve failed checked initialization: {error:?}")}),
                ),
                1,
            )
        }
    };

    let parse = |wire: &Option<DivisorWire>, name: &str| -> Result<MumfordDivisor, String> {
        divisor_from_wire(required_divisor(wire, name)?, &model.field, &curve)
    };
    let computed: Result<Value, String> = match request.operation {
        ArithmeticOperation::Identity => Ok(json!({
            "divisor": divisor_to_wire(&MumfordDivisor::identity(curve.m)),
        })),
        ArithmeticOperation::Validate => parse(&request.left, "left")
            .map(|divisor| json!({"valid": true, "divisor": divisor_to_wire(&divisor)})),
        ArithmeticOperation::Equal => (|| {
            let left = parse(&request.left, "left")?;
            let right = parse(&request.right, "right")?;
            let equal = left.try_eq(&right, &curve).map_err(|e| format!("{e:?}"))?;
            Ok(json!({"equal": equal}))
        })(),
        ArithmeticOperation::Negate => (|| {
            let left = parse(&request.left, "left")?;
            let result = left.try_neg(&curve).map_err(|e| format!("{e:?}"))?;
            Ok(json!({"divisor": divisor_to_wire(&result)}))
        })(),
        ArithmeticOperation::Add => (|| {
            let left = parse(&request.left, "left")?;
            let right = parse(&request.right, "right")?;
            let result = left.try_add(&right, &curve).map_err(|e| format!("{e:?}"))?;
            Ok(json!({"divisor": divisor_to_wire(&result)}))
        })(),
        ArithmeticOperation::ScalarMultiply => (|| {
            let left = parse(&request.left, "left")?;
            let scalar = request
                .scalar
                .as_deref()
                .ok_or_else(|| "operation requires `scalar`".to_string())
                .and_then(parse_number)?;
            let result = left
                .try_scalar_mul(&scalar, &curve)
                .map_err(|e| format!("{e:?}"))?;
            Ok(json!({"divisor": divisor_to_wire(&result), "scalar": format!("0x{scalar:x}")}))
        })(),
    };
    arithmetic_result(
        request.operation,
        json!({
            "family": "binary_cubic_pullback_v1",
            "curve_construction": certificate,
            "infinity_configuration": "one_rational_point",
        }),
        computed,
    )
}

fn prime_poly(encoded: &[String], p: &BigUint, name: &str) -> Result<FpPoly, String> {
    if p < &BigUint::from(2u32) {
        return Err("prime-field modulus must be at least 2".into());
    }
    if encoded.len() > 8194 {
        return Err(format!(
            "{name} exceeds the supported polynomial-size limit"
        ));
    }
    let coefficients = encoded
        .iter()
        .map(|text| parse_number(text))
        .collect::<Result<Vec<_>, _>>()?;
    if coefficients.iter().any(|coefficient| coefficient >= p) {
        return Err(format!("{name} has a noncanonical coefficient"));
    }
    if coefficients.last().is_some_and(Zero::is_zero) {
        return Err(format!("{name} has a noncanonical trailing zero"));
    }
    Ok(FpPoly::from_coeffs(coefficients, p.clone()))
}

fn prime_divisor_from_wire(
    wire: &DivisorWire,
    curve: &HyperellipticCurveP,
) -> Result<MumfordDivisorP, String> {
    let u = prime_poly(&wire.u, &curve.p, "u")?;
    let v = prime_poly(&wire.v, &curve.p, "v")?;
    MumfordDivisorP::try_from_polynomials(curve, u, v).map_err(|e| format!("{e:?}"))
}

fn encoded_prime_poly(poly: &FpPoly) -> Vec<String> {
    poly.coeffs
        .iter()
        .map(|coefficient| format!("0x{coefficient:x}"))
        .collect()
}

fn prime_divisor_to_wire(divisor: &MumfordDivisorP) -> DivisorWire {
    DivisorWire {
        u: encoded_prime_poly(&divisor.u),
        v: encoded_prime_poly(&divisor.v),
    }
}

fn odd_arithmetic(request: &ArithmeticRequest, supplied: &OddCurveWire) -> (Value, i32) {
    if supplied.family != "odd_prime_one_infinity_v1" || supplied.infinity != "one_rational" {
        return (
            response(
                "arithmetic",
                "unsupported",
                json!({
                    "reason": "odd arithmetic requires family odd_prime_one_infinity_v1 and infinity one_rational",
                    "arithmetic_status": "unsupported_model",
                }),
            ),
            2,
        );
    }
    let p = match parse_number(&supplied.p) {
        Ok(p) => p,
        Err(reason) => {
            return (
                response("arithmetic", "invalid_input", json!({"reason": reason})),
                3,
            )
        }
    };
    if p < BigUint::from(3u32) || !p.bit(0) {
        return (
            response(
                "arithmetic",
                "invalid_input",
                json!({"reason": "odd one-infinity arithmetic requires an odd modulus at least 3"}),
            ),
            3,
        );
    }
    let f = match prime_poly(&supplied.f, &p, "f") {
        Ok(f) => f,
        Err(reason) => {
            return (
                response("arithmetic", "invalid_input", json!({"reason": reason})),
                3,
            )
        }
    };
    let curve = match HyperellipticCurveP::try_new(p, f, supplied.genus) {
        Ok(curve) => curve,
        Err(error) => {
            return (
                response(
                    "arithmetic",
                    "invalid_input",
                    json!({"reason": format!("{error:?}")}),
                ),
                3,
            )
        }
    };
    let parse = |wire: &Option<DivisorWire>, name: &str| -> Result<MumfordDivisorP, String> {
        prime_divisor_from_wire(required_divisor(wire, name)?, &curve)
    };
    let computed: Result<Value, String> = match request.operation {
        ArithmeticOperation::Identity => MumfordDivisorP::try_identity(&curve)
            .map(|divisor| json!({"divisor": prime_divisor_to_wire(&divisor)}))
            .map_err(|e| format!("{e:?}")),
        ArithmeticOperation::Validate => parse(&request.left, "left")
            .map(|divisor| json!({"valid": true, "divisor": prime_divisor_to_wire(&divisor)})),
        ArithmeticOperation::Equal => (|| {
            let left = parse(&request.left, "left")?;
            let right = parse(&request.right, "right")?;
            let equal = left
                .checked_eq(&right, &curve)
                .map_err(|e| format!("{e:?}"))?;
            Ok(json!({"equal": equal}))
        })(),
        ArithmeticOperation::Negate => (|| {
            let left = parse(&request.left, "left")?;
            let result = left.checked_neg(&curve).map_err(|e| format!("{e:?}"))?;
            Ok(json!({"divisor": prime_divisor_to_wire(&result)}))
        })(),
        ArithmeticOperation::Add => (|| {
            let left = parse(&request.left, "left")?;
            let right = parse(&request.right, "right")?;
            let result = left
                .checked_add(&right, &curve)
                .map_err(|e| format!("{e:?}"))?;
            Ok(json!({"divisor": prime_divisor_to_wire(&result)}))
        })(),
        ArithmeticOperation::ScalarMultiply => (|| {
            let left = parse(&request.left, "left")?;
            let scalar = request
                .scalar
                .as_deref()
                .ok_or_else(|| "operation requires `scalar`".to_string())
                .and_then(parse_number)?;
            let result = left
                .checked_scalar_mul(&scalar, &curve)
                .map_err(|e| format!("{e:?}"))?;
            Ok(
                json!({"divisor": prime_divisor_to_wire(&result), "scalar": format!("0x{scalar:x}")}),
            )
        })(),
    };
    arithmetic_result(
        request.operation,
        json!({
            "family": supplied.family,
            "field": {"characteristic": supplied.p, "representation": "prime_residue"},
            "curve": supplied,
            "infinity_configuration": "one_rational_point",
        }),
        computed,
    )
}

fn arithmetic(request: ArithmeticRequest) -> (Value, i32) {
    if request.schema != ARITHMETIC_SCHEMA {
        return (
            response(
                "arithmetic",
                "unsupported",
                json!({"reason": "unsupported arithmetic request schema"}),
            ),
            2,
        );
    }
    match (&request.model, &request.odd_curve) {
        (Some(model), None) => binary_arithmetic(&request, model),
        (None, Some(curve)) => odd_arithmetic(&request, curve),
        _ => (
            response(
                "arithmetic",
                "invalid_input",
                json!({"reason": "supply exactly one of `model` or `odd_curve`"}),
            ),
            3,
        ),
    }
}

fn execute(command: Command) -> (Value, i32, Option<PathBuf>) {
    match command {
        Command::Construct { input, output } => match parse_json(&input) {
            Ok(model) => {
                let (value, code) = construct(model);
                (value, code, output)
            }
            Err(reason) => (
                response("construct", "invalid_input", json!({"reason": reason})),
                3,
                output,
            ),
        },
        Command::Verify { input, output } => {
            let parsed = read_bytes(&input).and_then(|bytes| {
                serde_json::from_slice::<VerifyRequest>(&bytes)
                    .map_err(|e| format!("parse verification request: {e}"))
            });
            match parsed {
                Ok(request) => {
                    let (value, code) = verify(request);
                    (value, code, output)
                }
                Err(reason) => (
                    response("verify", "invalid_input", json!({"reason": reason})),
                    3,
                    output,
                ),
            }
        }
        Command::Arithmetic { input, output } => {
            let parsed = read_bytes(&input).and_then(|bytes| {
                serde_json::from_slice::<ArithmeticRequest>(&bytes)
                    .map_err(|e| format!("parse arithmetic request: {e}"))
            });
            match parsed {
                Ok(request) => {
                    let (value, code) = arithmetic(request);
                    (value, code, output)
                }
                Err(reason) => (
                    response("arithmetic", "invalid_input", json!({"reason": reason})),
                    3,
                    output,
                ),
            }
        }
        Command::CatalogExport { registry, output } => match read_bytes(&registry) {
            Ok(bytes) => match checker::catalog(&bytes) {
                Ok(mut catalog) => {
                    catalog["generated_by"] = json!("hyperelliptic-cover catalog-export");
                    let producer_sources = [
                        include_str!("curve_cover_check/checker.rs"),
                        include_str!("hyperelliptic_cover.rs"),
                    ]
                    .concat();
                    catalog["producer_source_sha256"] =
                        json!(checker::digest(producer_sources.as_bytes()));
                    (catalog, 0, output)
                }
                Err(reason) => (
                    response("catalog_export", "invalid_input", json!({"reason": reason})),
                    3,
                    output,
                ),
            },
            Err(reason) => (
                response("catalog_export", "invalid_input", json!({"reason": reason})),
                3,
                output,
            ),
        },
    }
}

fn write_value(value: &Value, output: Option<PathBuf>) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value).map_err(|e| e.to_string())? + "\n";
    if let Some(path) = output {
        fs::write(&path, text).map_err(|e| format!("write {}: {e}", path.display()))
    } else {
        std::io::stdout()
            .lock()
            .write_all(text.as_bytes())
            .map_err(|e| format!("write stdout: {e}"))
    }
}

fn main() {
    let args = Args::parse();
    let (value, mut code, output) = execute(args.command);
    if let Err(reason) = write_value(&value, output) {
        eprintln!(
            "{}",
            serde_json::to_string(&response("write", "error", json!({"reason": reason})))
                .expect("static response serializes")
        );
        code = 1;
    }
    if code != 0 {
        std::process::exit(code);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn binary_model() -> Value {
        let modulus = BigUint::from(0xbu32);
        let hash = checker::digest(format!("f2m-modulus:0x{modulus:x}").as_bytes());
        json!({
            "v": "1",
            "form": "y^2+xy=x^3+a*x^2+b",
            "modulus": "0xb",
            "field": format!("f2m-3-{}", &hash[..8]),
            "a": "1",
            "b": "1",
        })
    }

    #[test]
    fn construction_has_separate_verified_claims() {
        let (answer, code) = construct(binary_model());
        assert_eq!(code, 0);
        assert_eq!(answer["status"], "ok");
        assert_eq!(
            answer["result"]["verification"]["defining_equation_identity"],
            "verified_exact_polynomial_identity"
        );
        assert_eq!(
            answer["result"]["verification"]["jacobian_arithmetic"],
            "implemented_checked_one_infinity_model"
        );
    }

    #[test]
    fn construction_rejects_unknown_model_fields() {
        let mut model = binary_model();
        model
            .as_object_mut()
            .unwrap()
            .insert("unexpected".into(), json!(true));
        let (answer, code) = construct(model);
        assert_eq!(code, 3);
        assert_eq!(answer["status"], "invalid_input");
    }

    #[test]
    fn prime_cover_arithmetic_is_structured_unsupported() {
        let request = ArithmeticRequest {
            schema: ARITHMETIC_SCHEMA.into(),
            model: Some(json!({
                "v":"1", "form":"y^2=x^3+a*x+b", "p":"11",
                "field":"fp-11", "a":"1", "b":"1"
            })),
            odd_curve: None,
            operation: ArithmeticOperation::Identity,
            left: None,
            right: None,
            scalar: None,
        };
        let (answer, code) = arithmetic(request);
        assert_eq!(code, 2);
        assert_eq!(answer["status"], "unsupported");
        assert_eq!(
            answer["result"]["arithmetic_status"],
            "unsupported_infinity_configuration"
        );
    }

    #[test]
    fn malformed_divisor_is_rejected_before_arithmetic() {
        let request = ArithmeticRequest {
            schema: ARITHMETIC_SCHEMA.into(),
            model: Some(binary_model()),
            odd_curve: None,
            operation: ArithmeticOperation::Validate,
            left: Some(DivisorWire {
                u: Vec::new(),
                v: Vec::new(),
            }),
            right: None,
            scalar: None,
        };
        let (answer, code) = arithmetic(request);
        assert_eq!(code, 3);
        assert_eq!(answer["status"], "invalid_input");
    }

    #[test]
    fn odd_one_infinity_branch_divisor_is_supported_separately() {
        let request = ArithmeticRequest {
            schema: ARITHMETIC_SCHEMA.into(),
            model: None,
            odd_curve: Some(OddCurveWire {
                family: "odd_prime_one_infinity_v1".into(),
                p: "5".into(),
                genus: 2,
                // y^2 = x^5 + x; derivative is 1, so this is squarefree.
                f: ["0", "1", "0", "0", "0", "1"]
                    .into_iter()
                    .map(str::to_string)
                    .collect(),
                infinity: "one_rational".into(),
            }),
            operation: ArithmeticOperation::Validate,
            // The Weierstrass point (0,0) gives the order-two class (x,0).
            left: Some(DivisorWire {
                u: vec!["0".into(), "1".into()],
                v: Vec::new(),
            }),
            right: None,
            scalar: None,
        };
        let (answer, code) = arithmetic(request);
        assert_eq!(code, 0, "{answer}");
        assert_eq!(answer["status"], "ok");
        assert!(answer["result"]["result"]["valid"]
            .as_bool()
            .expect("valid is a boolean"));
    }

    #[test]
    fn zero_prime_modulus_is_a_structured_error_not_a_panic() {
        let request = ArithmeticRequest {
            schema: ARITHMETIC_SCHEMA.into(),
            model: None,
            odd_curve: Some(OddCurveWire {
                family: "odd_prime_one_infinity_v1".into(),
                p: "0".into(),
                genus: 2,
                f: vec!["1".into()],
                infinity: "one_rational".into(),
            }),
            operation: ArithmeticOperation::Identity,
            left: None,
            right: None,
            scalar: None,
        };
        let (answer, code) = arithmetic(request);
        assert_eq!(code, 3);
        assert_eq!(answer["status"], "invalid_input");
    }
}
