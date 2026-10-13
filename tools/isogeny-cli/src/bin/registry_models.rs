//! Native normalization of explicitly declared prime models for the catalogue.
//! The source records are historical evidence and are never rewritten.
use num_bigint::{BigInt, BigUint};
use num_traits::{One, ToPrimitive, Zero};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{fs, path::Path, process::Command};

fn hash(bytes: &[u8]) -> String {
    hex::encode(Sha256::digest(bytes))
}

fn uint(value: &Value) -> Result<BigUint, String> {
    let text = match value {
        Value::String(text) => text.clone(),
        Value::Number(number) => number.to_string(),
        _ => return Err("expected a nonnegative integer".into()),
    };
    let (digits, radix) = text
        .strip_prefix("0x")
        .map_or((text.as_str(), 10), |s| (s, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix).ok_or_else(|| format!("bad integer {text}"))
}

fn form(record: &Value) -> Option<String> {
    record
        .get("model")
        .or_else(|| record.get("form"))?
        .as_str()
        .map(|text| text.split_whitespace().collect())
}

fn sub(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
    (a + p - b % p) % p
}

/// With u=x+a/3, y²=x³+a x²+b x becomes y²=u³+A u+B.
fn short(p: &BigUint, a: &BigUint, b: &BigUint) -> Result<(BigUint, BigUint, BigUint), String> {
    let prime = p
        .to_u64()
        .ok_or("explicit quadratic model fields above u64 are not supported")?;
    if prime <= 3 || !isogeny_algos::field::is_prime(prime) {
        return Err("quadratic Weierstrass normalization requires a prime above 3".into());
    }
    let inv3 = BigUint::from(3u8).modpow(&(p - 2u8), p);
    if (&inv3 * 3u8) % p != BigUint::one() {
        return Err("3 is not invertible in the declared field".into());
    }
    let shift = a * &inv3 % p;
    let big_a = sub(b, &(a * &shift % p), p);
    let cube = shift.modpow(&BigUint::from(3u8), p);
    let big_b = sub(&(&cube * 2u8 % p), &(b * &shift % p), p);
    Ok((big_a, big_b, shift))
}

fn identity(p: &BigUint, a: &BigUint, b: &BigUint, order: &BigUint) -> Result<Value, String> {
    let trace: BigInt = BigInt::from(p.clone()) + BigInt::one() - BigInt::from(order.clone());
    if &trace * &trace > BigInt::from(p * 4u8) {
        return Err("declared order violates the Hasse bound".into());
    }
    let four_a3 = a.modpow(&BigUint::from(3u8), p) * 4u8 % p;
    let discriminant = (&four_a3 + b * b * 27u8) % p;
    if discriminant.is_zero() {
        return Err("singular normalized model".into());
    }
    let j = BigUint::from(1728u16) * &four_a3 * discriminant.modpow(&(p - 2u8), p) % p;
    let model_json = format!(
        r#"{{"a":"{a}","b":"{b}","field":"fp-{p}","form":"y^2=x^3+a*x+b","p":"{p}","v":"1"}}"#
    );
    let digest = hash(model_json.as_bytes());
    let trace_tag = if trace < BigInt::zero() {
        format!("tm{}", -&trace)
    } else {
        format!("t{trace}")
    };
    let slug = format!("icv1-fp{}-{trace_tag}-{}", p.bits(), &digest[..8]);
    let icv1 = format!(
        "ICV1:fp-{p}:{trace}:{order}:{j}:unk:unk:r:{}",
        &digest[..12]
    );
    Ok(json!({"slug":slug,"icv1":icv1,"model_json":model_json,
        "model_sha256":digest,"trace":serde_json::from_str::<Value>(&trace.to_string()).map_err(|e| e.to_string())?,
        "order":order.to_string(),"j":j.to_string(),"end":"unk"}))
}

fn normalize(record: &Value) -> Result<Value, String> {
    let p = uint(&record["p"])?;
    let a = uint(&record["a"])?;
    let b = uint(&record["b"])?;
    let order = uint(record.get("group_order").unwrap_or(&record["order"]))?;
    let (big_a, big_b, shift) = short(&p, &a, &b)?;
    let id = identity(&p, &big_a, &big_b, &order)?;
    for (source_key, id_key) in [("icv1", "icv1"), ("slug", "slug"), ("icv1_slug", "slug")] {
        if let Some(declared) = record.get(source_key) {
            if declared != &id[id_key] {
                return Err(format!(
                    "declared {source_key} differs from native normalized identity"
                ));
            }
        }
    }
    Ok(
        json!({"identity":id,"params":{"p":p.to_string(),"a":big_a.to_string(),"b":big_b.to_string()},
        "coordinate_map":{"u":"x+a/3","v":"y","x_shift":shift.to_string()},
        "original":record}),
    )
}

fn walk(
    value: &Value,
    source: &str,
    source_hash: &str,
    pointer: &str,
    rows: &mut Vec<Value>,
) -> Result<(), String> {
    match value {
        Value::Object(object) => {
            if form(value).as_deref() == Some("y^2=x^3+ax^2+bx")
                && ["p", "a", "b"].iter().all(|key| object.contains_key(*key))
                && (object.contains_key("order") || object.contains_key("group_order"))
            {
                let mut normalized =
                    normalize(value).map_err(|error| format!("{source}{pointer}: {error}"))?;
                normalized["source"] = json!(source);
                normalized["source_sha256"] = json!(source_hash);
                normalized["pointer"] = json!(pointer);
                rows.push(normalized);
            }
            for (key, child) in object {
                let escaped = key.replace('~', "~0").replace('/', "~1");
                walk(
                    child,
                    source,
                    source_hash,
                    &format!("{pointer}/{escaped}"),
                    rows,
                )?;
            }
        }
        Value::Array(array) => {
            for (index, child) in array.iter().enumerate() {
                walk(
                    child,
                    source,
                    source_hash,
                    &format!("{pointer}/{index}"),
                    rows,
                )?;
            }
        }
        _ => {}
    }
    Ok(())
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let check = args.iter().any(|arg| arg == "--check");
    let out = args
        .windows(2)
        .find(|pair| pair[0] == "--out")
        .map_or("docs/curves/sources/prime-models.json", |pair| {
            pair[1].as_str()
        });
    // Sparse checkouts must not silently drop a previously normalized input.
    if let Ok(previous) = fs::read_to_string(out) {
        let previous: Value = serde_json::from_str(&previous).map_err(|e| e.to_string())?;
        if let Some(records) = previous["records"].as_array() {
            for record in records {
                let source = record["source"]
                    .as_str()
                    .ok_or("missing normalization source")?;
                if !Path::new(source).is_file() {
                    return Err(format!("required normalization input {source} is off disk; materialize it before rebuilding"));
                }
            }
        }
    }
    let names = Command::new("git")
        .args(["ls-files", "-z", "--", "*.json"])
        .output()
        .map_err(|e| e.to_string())?;
    if !names.status.success() {
        return Err("git ls-files failed".into());
    }
    let mut records = Vec::new();
    for bytes in names
        .stdout
        .split(|byte| *byte == 0)
        .filter(|bytes| !bytes.is_empty())
    {
        let source = std::str::from_utf8(bytes).map_err(|e| e.to_string())?;
        if source == "docs/curves/sources/prime-models.json" || source == out {
            continue;
        }
        let path = Path::new(source);
        let Ok(metadata) = path.metadata() else {
            continue;
        };
        if metadata.len() > 8 << 20 {
            continue;
        }
        let bytes = fs::read(path).map_err(|e| e.to_string())?;
        let Ok(value) = serde_json::from_slice(&bytes) else {
            continue;
        };
        walk(&value, source, &hash(&bytes), "", &mut records)?;
    }
    let result = json!({"schema":"prime-models-normalization/v1",
        "generated_by":"registry_models (native Rust)",
        "what_this_is":"Exact short-Weierstrass metadata derived from explicitly declared source models; historical source bytes are unchanged.",
        "records":records});
    let text = serde_json::to_string_pretty(&result).map_err(|e| e.to_string())? + "\n";
    if check {
        if fs::read_to_string(out).map_err(|e| e.to_string())? != text {
            return Err(format!("{out} is stale; run registry_models"));
        }
        println!("{out} is current");
    } else {
        fs::write(out, text).map_err(|e| e.to_string())?;
        println!(
            "wrote {out}: {} normalized records",
            result["records"].as_array().unwrap().len()
        );
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn count(p: u64, a: u64, b: u64, quadratic: bool) -> u64 {
        let mut squares = vec![0u64; p as usize];
        for y in 0..p {
            squares[(y * y % p) as usize] += 1;
        }
        1 + (0..p)
            .map(|x| {
                let rhs = if quadratic {
                    (x * x % p * x % p + a * x % p * x % p + b * x % p) % p
                } else {
                    (x * x % p * x % p + a * x % p + b) % p
                };
                squares[rhs as usize]
            })
            .sum::<u64>()
    }
    #[test]
    fn translated_coordinates_are_a_bijection_on_small_fields() {
        for p in [5u64, 7, 11, 17, 101] {
            for a in 0..p {
                for b in 0..p {
                    let (big_a, big_b, shift) = short(&p.into(), &a.into(), &b.into()).unwrap();
                    let a_short: u64 = big_a.to_string().parse().unwrap();
                    let b_short: u64 = big_b.to_string().parse().unwrap();
                    let shift: u64 = shift.to_string().parse().unwrap();
                    let mut seen = vec![false; p as usize];
                    for x in 0..p {
                        let u = (x + shift) % p;
                        assert!(!seen[u as usize]);
                        seen[u as usize] = true;
                        assert_eq!(
                            (x * x % p * x % p + a * x % p * x % p + b * x % p) % p,
                            (u * u % p * u % p + a_short * u % p + b_short) % p
                        );
                    }
                    assert!(seen.into_iter().all(|x| x));
                    assert_eq!(count(p, a, b, true), count(p, a_short, b_short, false));
                }
            }
        }
    }
    #[test]
    fn both_frozen_agm_identities_and_orders_match_the_native_adapter() {
        let p = 1048583u64;
        for (a, b, big_a, big_b, slug, j) in [
            (
                2u64,
                2u64,
                699056u64,
                737891u64,
                "icv1-fp21-tm700-d62ba3a6",
                "128",
            ),
            (
                1048579,
                1048579,
                699046,
                388354,
                "icv1-fp21-tm700-07d957b2",
                "10976",
            ),
        ] {
            let record = json!({"p":p,"a":a,"b":b,"order":1049284,"model":"y^2 = x^3 + a x^2 + b x","slug":slug});
            let result = normalize(&record).unwrap();
            assert_eq!(result["params"]["a"], big_a.to_string());
            assert_eq!(result["params"]["b"], big_b.to_string());
            assert_eq!(result["identity"]["slug"], slug);
            assert_eq!(result["identity"]["j"], j);
            assert_eq!(count(p, a, b, true), 1049284);
            assert_eq!(count(p, big_a, big_b, false), 1049284);
            assert_ne!(count(p, a, b, false), 1049284);
        }
    }
    #[test]
    fn a_conflicting_historical_identity_is_rejected() {
        let record =
            json!({"p":1048583,"a":2,"b":2,"order":1049284,"slug":"icv1-fp21-tm700-56b6783a"});
        assert!(normalize(&record).unwrap_err().contains("declared slug"));
        let record = json!({"p":1048583,"a":2,"b":2,"order":1049284,"icv1":"wrong"});
        assert!(normalize(&record).unwrap_err().contains("declared icv1"));
    }

    #[test]
    fn composite_fields_and_singular_models_are_rejected() {
        for p in [9u64, 25, 341, 561] {
            let record = json!({"p":p,"a":2,"b":2,"order":p+1});
            assert!(normalize(&record).unwrap_err().contains("requires a prime"));
        }
        let record = json!({"p":101,"a":0,"b":0,"order":102});
        assert!(normalize(&record).unwrap_err().contains("singular"));
    }
}
