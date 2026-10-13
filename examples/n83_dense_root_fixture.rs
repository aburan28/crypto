//! Freeze public N83 target points for paired packed-root table comparisons.
//! Usage: n83_dense_root_fixture <count> <points.jsonl> <answers.jsonl>
//! The producer receives only points.jsonl; answers.jsonl stays with the
//! independent verifier. Generation is outside every measured online interval.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::json;

const DOMAIN: &[u8] = b"n83-dense-root-held-out-20261013-v1";
const SUBGROUP_ORDER: u64 = 8_569_786_107_849_059;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: n83_dense_root_fixture <count> <points.jsonl> <answers.jsonl>".into());
    }
    let count: u32 = args[1].parse()?;
    if count == 0 || count > 4096 {
        return Err("count must be in 1..=4096".into());
    }
    let curve = KoblitzCurve::new(1, 83).ok_or("N83 Koblitz construction failed")?;
    if curve.subgroup_order.to_u64() != Some(SUBGROUP_ORDER)
        || curve.curve.irreducible.low_terms != [0, 2, 4, 7]
    {
        return Err("unexpected curve identity".into());
    }
    let generator = curve.generator();
    let mut points = BufWriter::new(File::create(&args[2])?);
    let mut answers = BufWriter::new(File::create(&args[3])?);
    let mut distinct = HashSet::new();
    for index in 0..count {
        let mut hasher = blake3::Hasher::new();
        hasher.update(DOMAIN);
        hasher.update(&index.to_le_bytes());
        let digest = hasher.finalize();
        let mut scalar_bytes = [0u8; 8];
        scalar_bytes.copy_from_slice(&digest.as_bytes()[..8]);
        let scalar = 1 + u64::from_le_bytes(scalar_bytes) % (SUBGROUP_ORDER - 1);
        let target = curve.mul(generator, &BigUint::from(scalar));
        let BinaryPoint::Affine { x, y } = target else {
            return Err(format!("target {index} is infinity").into());
        };
        let coordinates = [x.to_biguint().to_string(), y.to_biguint().to_string()];
        if !distinct.insert(coordinates.clone()) {
            return Err(format!("duplicate target {index}").into());
        }
        writeln!(points, "{}", json!(coordinates))?;
        writeln!(
            answers,
            "{}",
            json!({
                "index": index,
                "scalar": scalar.to_string(),
                "point": coordinates,
                "digest_blake3": digest.to_hex().to_string(),
            })
        )?;
    }
    points.flush()?;
    answers.flush()?;
    eprintln!(
        "frozen {count} distinct N83 public targets; domain={}",
        String::from_utf8_lossy(DOMAIN)
    );
    Ok(())
}
