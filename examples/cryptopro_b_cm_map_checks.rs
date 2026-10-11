//! Replay public forward-map correctness checks for CryptoPro-B.
//! All point inputs have known scalars; this example performs no log search.

use crypto_lib::hash::sha256::sha256;
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};
use serde_json::{json, Value};

const MAP_BYTES: &[u8] = include_bytes!(
    "../research/cryptopro_b_cm_map_checks_20261005/endomorphism155.json"
);
const FIXTURE_BYTES: &[u8] = include_bytes!(
    "../research/cryptopro_b_cm_map_checks_20261005/evidence/legacy/results.json"
);
const MAP_SHA256: &str = "51ad7edaf9c931364be2513a1163a95474858381152a715ad9140648fad6a042";
const FIXTURE_SHA256: &str = "e69d3bc1bcd7b8cfccaeb8f0e332dfb4d23c348ea4f944b550387155166d99b5";
const VERIFIER_BYTES: &[u8] = include_bytes!("cryptopro_b_cm_map_checks.rs");

type Check<T> = Result<T, String>;
type Poly = Vec<BigUint>;
type Point = Option<(BigUint, BigUint)>;

fn require(condition: bool, message: &str) -> Check<()> {
    if condition {
        Ok(())
    } else {
        Err(message.to_owned())
    }
}

fn decimal(text: &str) -> Check<BigUint> {
    BigUint::parse_bytes(text.as_bytes(), 10).ok_or_else(|| "Invalid integer".to_owned())
}

fn hexadecimal(text: &str) -> BigUint {
    BigUint::parse_bytes(text.as_bytes(), 16).expect("valid constant")
}

// The archived JSON stores 255-bit coefficients as number tokens. Quote
// those tokens before serde_json sees them, so none passes through f64.
// Strings (including escaped quotes) are copied without interpretation.
fn preserve_number_tokens(input: &str) -> String {
    let mut output = String::with_capacity(input.len());
    let mut chars = input.chars().peekable();
    let mut in_string = false;
    let mut escaped = false;
    while let Some(ch) = chars.next() {
        if in_string {
            output.push(ch);
            if escaped {
                escaped = false;
            } else if ch == '\\' {
                escaped = true;
            } else if ch == '"' {
                in_string = false;
            }
        } else if ch == '"' {
            in_string = true;
            output.push(ch);
        } else if ch == '-' || ch.is_ascii_digit() {
            output.push('"');
            output.push(ch);
            while chars.peek().is_some_and(|c| {
                c.is_ascii_digit() || matches!(c, '.' | 'e' | 'E' | '+' | '-')
            }) {
                output.push(chars.next().expect("peeked character"));
            }
            output.push('"');
        } else {
            output.push(ch);
        }
    }
    output
}

fn integer(value: &Value) -> Check<BigUint> {
    decimal(value.as_str().ok_or_else(|| "Missing integer token".to_owned())?)
}

fn polynomial(value: &Value) -> Check<Poly> {
    value
        .as_array()
        .ok_or_else(|| "Missing polynomial".to_owned())?
        .iter()
        .map(integer)
        .collect()
}

#[derive(Clone)]
struct Component {
    degree: usize,
    source_a: BigUint,
    source_b: BigUint,
    target_a: BigUint,
    target_b: BigUint,
    kernel: Poly,
    numerator: Poly,
    denominator: Poly,
}

struct MapData {
    p: BigUint,
    n: BigUint,
    scalar: BigUint,
    u: BigUint,
    maps: Vec<Component>,
}

fn read_map() -> Check<MapData> {
    require(hex::encode(sha256(MAP_BYTES)) == MAP_SHA256, "Map hash mismatch")?;
    let text = std::str::from_utf8(MAP_BYTES).map_err(|e| e.to_string())?;
    let data: Value = serde_json::from_str(&preserve_number_tokens(text))
        .map_err(|e| e.to_string())?;
    require(data["coefficient_order"] == "ascending powers", "Coefficient order mismatch")?;
    require(data["endomorphism_label"] == "omega", "Map label mismatch")?;
    let mut maps = Vec::new();
    for map in data["maps"].as_array().ok_or("Missing component maps")? {
        maps.push(Component {
            degree: map["degree"].as_str().ok_or("Missing degree")?
                .parse().map_err(|e: std::num::ParseIntError| e.to_string())?,
            source_a: integer(&map["source_A"])? ,
            source_b: integer(&map["source_B"])? ,
            target_a: integer(&map["target_A"])? ,
            target_b: integer(&map["target_B"])? ,
            kernel: polynomial(&map["kernel"])? ,
            numerator: polynomial(&map["x_numerator"])? ,
            denominator: polynomial(&map["x_denominator"])? ,
        });
    }
    Ok(MapData {
        p: integer(&data["p"])? ,
        n: integer(&data["n"])? ,
        scalar: integer(&data["endomorphism_scalar"])? ,
        u: integer(&data["final_isomorphism_u"])? ,
        maps,
    })
}

fn subtract(x: &BigUint, y: &BigUint, modulus: &BigUint) -> BigUint {
    (x + modulus - y % modulus) % modulus
}

fn inverse(x: &BigUint, modulus: &BigUint) -> Check<BigUint> {
    let m = BigInt::from(modulus.clone());
    let gcd = BigInt::from(x.clone()).extended_gcd(&m);
    require(gcd.gcd == BigInt::one(), "Noninvertible denominator")?;
    gcd.x.mod_floor(&m).to_biguint().ok_or_else(|| "Invalid inverse".to_owned())
}

struct Field(BigUint);

impl Field {
    fn trim(&self, mut poly: Poly) -> Poly {
        for coefficient in &mut poly {
            *coefficient %= &self.0;
        }
        while poly.len() > 1 && poly.last().is_some_and(Zero::is_zero) {
            poly.pop();
        }
        if poly.is_empty() { poly.push(BigUint::zero()); }
        poly
    }

    fn plus(&self, lhs: &[BigUint], rhs: &[BigUint]) -> Poly {
        let mut result = vec![BigUint::zero(); lhs.len().max(rhs.len())];
        for (i, coefficient) in lhs.iter().enumerate() { result[i] += coefficient; }
        for (i, coefficient) in rhs.iter().enumerate() { result[i] += coefficient; }
        self.trim(result)
    }

    fn scale(&self, poly: &[BigUint], coefficient: &BigUint) -> Poly {
        self.trim(poly.iter().map(|v| v * coefficient).collect())
    }

    fn times(&self, lhs: &[BigUint], rhs: &[BigUint]) -> Poly {
        let mut result = vec![BigUint::zero(); lhs.len() + rhs.len() - 1];
        for (i, a) in lhs.iter().enumerate() {
            for (j, b) in rhs.iter().enumerate() { result[i + j] += a * b; }
        }
        self.trim(result)
    }

    fn derivative(&self, poly: &[BigUint]) -> Poly {
        self.trim(poly.iter().enumerate().skip(1).map(|(i, v)| v * i).collect())
    }

    fn evaluate(&self, poly: &[BigUint], x: &BigUint) -> BigUint {
        poly.iter().rev().fold(BigUint::zero(), |r, v| (r * x + v) % &self.0)
    }
}

struct Curve {
    field: Field,
    a: BigUint,
    b: BigUint,
}

impl Curve {
    fn on_curve(&self, point: &Point, a: &BigUint, b: &BigUint) -> bool {
        point.as_ref().is_none_or(|(x, y)| {
            (y * y) % &self.field.0 == (x * x * x + a * x + b) % &self.field.0
        })
    }

    fn negate(&self, point: &Point) -> Point {
        point.as_ref().map(|(x, y)| (x.clone(), subtract(&BigUint::zero(), y, &self.field.0)))
    }

    fn add(&self, lhs: &Point, rhs: &Point) -> Check<Point> {
        let (Some((x, y)), Some((u, v))) = (lhs, rhs) else {
            return Ok(if lhs.is_none() { rhs.clone() } else { lhs.clone() });
        };
        let p = &self.field.0;
        let slope = if x == u {
            if (y + v) % p == BigUint::zero() { return Ok(None); }
            (x * x * 3u32 + &self.a) * inverse(&(y * 2u32), p)? % p
        } else {
            subtract(v, y, p) * inverse(&subtract(u, x, p), p)? % p
        };
        let z = subtract(&subtract(&(&slope * &slope), x, p), u, p);
        let w = subtract(&(&slope * subtract(x, &z, p)), y, p);
        Ok(Some((z, w)))
    }

    fn multiply(&self, scalar: &BigUint, point: &Point) -> Check<Point> {
        let mut k = scalar.clone();
        let mut current = point.clone();
        let mut result = None;
        while !k.is_zero() {
            if k.is_odd() { result = self.add(&result, &current)?; }
            k >>= 1;
            if !k.is_zero() { current = self.add(&current, &current)?; }
        }
        Ok(result)
    }

    fn cm_map(&self, data: &MapData, point: &Point) -> Check<Point> {
        let f = &self.field;
        let mut image = point.clone();
        for map in &data.maps {
            require(self.on_curve(&image, &map.source_a, &map.source_b), "Invalid map source")?;
            let Some((x, y)) = image.as_ref() else { continue; };
            let nv = f.evaluate(&map.numerator, x);
            let dv = f.evaluate(&map.denominator, x);
            if dv.is_zero() { image = None; continue; }
            let di = inverse(&dv, &f.0)?;
            let tangent = subtract(
                &(f.evaluate(&f.derivative(&map.numerator), x) * &dv),
                &(&nv * f.evaluate(&f.derivative(&map.denominator), x)),
                &f.0,
            );
            image = Some((&nv * &di % &f.0, y * tangent * &di * &di % &f.0));
            require(self.on_curve(&image, &map.target_a, &map.target_b), "Invalid map target")?;
        }
        if let Some((x, y)) = image {
            image = Some((
                x * inverse(&data.u.modpow(&BigUint::from(2u32), &f.0), &f.0)? % &f.0,
                y * inverse(&data.u.modpow(&BigUint::from(3u32), &f.0), &f.0)? % &f.0,
            ));
        }
        require(self.on_curve(&image, &self.a, &self.b), "Invalid final map point")?;
        Ok(image)
    }
}

fn verify_chain(curve: &Curve, data: &MapData) -> Check<()> {
    let f = &curve.field;
    require(data.maps.iter().map(|m| m.degree).collect::<Vec<_>>() == [5, 31], "Map degrees mismatch")?;
    let mut previous = (&curve.a, &curve.b);
    for map in &data.maps {
        require((&map.source_a, &map.source_b) == previous, "Map chain mismatch")?;
        require(f.times(&map.kernel, &map.kernel) == f.trim(map.denominator.clone()), "Kernel square mismatch")?;
        require(f.trim(map.numerator.clone()).len() == map.degree + 1
            && f.trim(map.denominator.clone()).len() == map.degree, "Polynomial degree mismatch")?;
        let tangent = f.plus(
            &f.times(&f.derivative(&map.numerator), &map.denominator),
            &f.scale(&f.times(&map.numerator, &f.derivative(&map.denominator)), &(&f.0 - 1u32)),
        );
        let source = vec![map.source_b.clone(), map.source_a.clone(), BigUint::zero(), BigUint::one()];
        let lhs = f.times(&source, &f.times(&tangent, &tangent));
        let den_square = f.times(&map.denominator, &map.denominator);
        let inner = f.plus(
            &f.plus(
                &f.times(&f.times(&map.numerator, &map.numerator), &map.numerator),
                &f.scale(&f.times(&map.numerator, &den_square), &map.target_a),
            ),
            &f.scale(&f.times(&den_square, &map.denominator), &map.target_b),
        );
        require(lhs == f.times(&map.denominator, &inner), "Exact polynomial identity failed")?;
        previous = (&map.target_a, &map.target_b);
    }
    require(*previous.0 == &curve.a * data.u.modpow(&BigUint::from(4u32), &f.0) % &f.0
        && *previous.1 == &curve.b * data.u.modpow(&BigUint::from(6u32), &f.0) % &f.0,
        "Final isomorphism failed")
}

fn curve() -> Curve {
    let p = (BigUint::one() << 255usize) + 3225u32;
    Curve {
        a: &p - 3u32,
        b: hexadecimal("3e1af419a269a5f866a7d3c25c3df80ae979259373ff2b182f49d4ce7e1bbc8b"),
        field: Field(p),
    }
}

fn verify() -> Check<Value> {
    require(
        hex::encode(sha256(FIXTURE_BYTES)) == FIXTURE_SHA256,
        "Known-input hash mismatch",
    )?;
    let data = read_map()?;
    let curve = curve();
    let p = &curve.field.0;
    let n = hexadecimal("800000000000000000000000000000015f700cfff1a624e5e497161bcc8a198f");
    require(data.p == *p && data.n == n, "Curve parameter mismatch")?;
    let generator = Some((BigUint::one(), hexadecimal(
        "3fa8124359f96680b83d1c3eb2c070e5c545c9858d03ecfb744bf8d717717efc",
    )));
    require(curve.on_curve(&generator, &curve.a, &curve.b)
        && curve.multiply(&n, &generator)?.is_none(), "Generator/order check failed")?;
    let t = BigInt::from(p.clone()) + 1u32 - BigInt::from(n.clone());
    let v = decimal("4646402506017662432554672533504826433")?;
    let vi = BigInt::from(v.clone());
    require(&t * &t - BigInt::from(p.clone()) * 4 == -619 * &vi * &vi, "CM discriminant failed")?;
    let factors = [103u64, 6217, 17233477, 3989640269, 105533924655391];
    require(factors.iter().fold(BigUint::one(), |product, x| product * x) == v, "Conductor product failed")?;
    let ni = BigInt::from(n.clone());
    require((&t - &vi).is_even(), "Frobenius parity failed")?;
    let lambda_integer = (BigInt::one() - (&t - &vi) / 2u32)
        * BigInt::from(inverse(&v, &n)?);
    let lambda = lambda_integer.mod_floor(&ni).to_biguint().ok_or("Invalid CM scalar")?;
    require(lambda == data.scalar, "Subgroup scalar mismatch")?;
    require(subtract(&(&lambda * &lambda + 155u32), &lambda, &n).is_zero(), "Scalar CM polynomial failed")?;
    let determinant = subtract(&BigUint::from(156u32),
        &(&lambda * subtract(&BigUint::one(), &lambda, &n)), &n);
    require(determinant == BigUint::one(), "Subgroup determinant failed")?;
    verify_chain(&curve, &data)?;

    let fixture: Value = serde_json::from_slice(FIXTURE_BYTES).map_err(|e| e.to_string())?;
    require(fixture["source_map_sha256"] == MAP_SHA256, "Fixture map hash mismatch")?;
    let fixture_records = fixture["records"].as_array().ok_or("Missing known inputs")?;
    require(fixture_records.len() == 39, "Known-input count mismatch")?;
    let mut points = Vec::new();
    let mut records = Vec::new();
    for record in fixture_records {
        let text = record["known_scalar"].as_str().ok_or("Missing known scalar")?;
        let scalar = BigInt::parse_bytes(text.as_bytes(), 10).ok_or("Invalid known scalar")?;
        let mut point = curve.multiply(scalar.magnitude(), &generator)?;
        if scalar.sign() == Sign::Minus { point = curve.negate(&point); }
        let image = curve.cm_map(&data, &point)?;
        require(image == curve.multiply(&lambda, &point)?, "Map/scalar action differs")?;
        let polynomial = curve.add(
            &curve.add(&curve.cm_map(&data, &image)?, &curve.negate(&image))?,
            &curve.multiply(&BigUint::from(155u32), &point)?,
        )?;
        require(polynomial.is_none(), "CM point polynomial failed")?;
        require(curve.cm_map(&data, &curve.negate(&point))? == curve.negate(&image), "Negation failed")?;
        records.push(json!({
            "known_scalar": text, "map_equals_scalar_action": true,
            "cm_map_polynomial": true, "negation_compatibility": true
        }));
        points.push(point);
    }
    let mut pair_checks = 0;
    for (lhs, rhs) in points[..16].iter().zip(&points[points.len() - 16..]) {
        require(curve.cm_map(&data, &curve.add(lhs, rhs)?)?
            == curve.add(&curve.cm_map(&data, lhs)?, &curve.cm_map(&data, rhs)?)?, "Additivity failed")?;
        pair_checks += 1;
    }
    require(Value::Array(records.clone()) == fixture["records"], "Frozen replay differs")?;
    Ok(json!({
        "curve": "id-GostR3410-2001-CryptoPro-B-ParamSet",
        "map": "degree-155 CM endomorphism omega", "implementation": "Rust",
        "source_map_sha256": MAP_SHA256,
        "known_input_fixture_sha256": FIXTURE_SHA256,
        "verifier_source_sha256": hex::encode(sha256(VERIFIER_BYTES)),
        "scope": "forward-map correctness on known public inputs",
        "genus_two_transfer_tested": false, "relation_solver_benchmarked": false,
        "exact_component_map_identities": [
            {"degree": 5, "exact_polynomial_identity": "passed"},
            {"degree": 31, "exact_polynomial_identity": "passed"}
        ],
        "generator_order": "passed", "cm_discriminant_identity": "passed",
        "frobenius_conductor_factor_product": "passed", "scalar_cm_polynomial": "passed",
        "polarization_subgroup_determinant": "passed",
        "known_input_cases": records.len(), "additivity_pairs": pair_checks,
        "failed_checks": 0, "records": records,
        "costs": null, "speedup": null,
        "limitation": "No genus-two transfer or relation-solving cost is tested."
    }))
}

fn main() {
    match verify() {
        Ok(report) => println!("{}", serde_json::to_string_pretty(&report).expect("JSON report")),
        Err(error) => {
            eprintln!("CM correctness check failed: {error}");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forward_checks_replay_all_known_inputs() {
        let report = verify().expect("public forward checks pass");
        assert_eq!(report["known_input_cases"], 39);
        assert_eq!(report["additivity_pairs"], 16);
        assert_eq!(report["failed_checks"], 0);
    }

    #[test]
    fn altered_coefficient_fails_exact_identity() {
        let mut data = read_map().unwrap();
        data.maps[0].numerator[0] += 1u32;
        assert_eq!(verify_chain(&curve(), &data).unwrap_err(), "Exact polynomial identity failed");
    }

    #[test]
    fn altered_final_isomorphism_is_rejected() {
        let mut data = read_map().unwrap();
        data.u += 1u32;
        assert_eq!(verify_chain(&curve(), &data).unwrap_err(), "Final isomorphism failed");
    }

    #[test]
    fn large_json_numbers_and_escaped_strings_are_preserved() {
        let input = r#"{"p":57896044618658097711785492504343953926634992332820282019728792003956564823297,"s":"12 \\\"34","t":-2.5e+4}"#;
        let data: Value = serde_json::from_str(&preserve_number_tokens(input)).unwrap();
        assert_eq!(data["p"], "57896044618658097711785492504343953926634992332820282019728792003956564823297");
        assert_eq!(data["s"], "12 \\\"34");
        assert_eq!(data["t"], "-2.5e+4");
    }
}
