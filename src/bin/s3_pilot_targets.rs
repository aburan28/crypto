//! Deterministic public one-target fixtures for the paired S3 study.
//!
//! Usage: s3_pilot_targets N COUNT SEED_LABEL NEW_OUTPUT_DIRECTORY
//! Scalar witnesses stay in `fixtures.json`; a solver receives only its
//! individual `Tnnn/public_target.json` file.

use crypto_lib::binary_ecc::curve::{scalar_mul, BinaryPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

fn fixture(curve: &KoblitzCurve, label: &str, index: u32) -> Result<Value, String> {
    let r = &curve.subgroup_order;
    let mut seed = Vec::with_capacity(label.len() + 10);
    seed.extend_from_slice(b"s3-pair-target-v1:");
    seed.extend_from_slice(label.as_bytes());
    seed.push(b':');
    seed.extend_from_slice(&curve.n.to_be_bytes());
    seed.extend_from_slice(&index.to_be_bytes());
    let digest = sha256(&seed);
    let scalar = BigUint::from(1u8) + BigUint::from_bytes_be(&digest) % (r - BigUint::from(1u8));
    let point = scalar_mul(&curve.curve, curve.generator(), &scalar);
    if !curve.curve.is_on_curve(&point)
        || scalar_mul(&curve.curve, &point, r) != BinaryPoint::Infinity
    {
        return Err("derived point is not in the declared subgroup".into());
    }
    let [x, y] = match point {
        BinaryPoint::Affine { x, y } => [x, y],
        BinaryPoint::Infinity => return Err("derived identity target".into()),
    };
    let words = [
        x.to_biguint().to_u64().ok_or("x does not fit u64")?,
        y.to_biguint().to_u64().ok_or("y does not fit u64")?,
    ];
    Ok(json!({
        "label": format!("T{index:03}"),
        "public_point": words,
        "fixture_scalar": scalar.to_string(),
        "scalar_derivation_sha256": hex::encode(digest),
        "workload_id": null,
    }))
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    fs::write(path, bytes).map_err(|e| e.to_string())
}

fn run() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        return Err("usage: s3_pilot_targets N COUNT SEED_LABEL NEW_OUTPUT_DIRECTORY".into());
    }
    let n: u32 = args[1].parse::<u32>().map_err(|e| e.to_string())?;
    let count: u32 = args[2].parse::<u32>().map_err(|e| e.to_string())?;
    if count == 0 || count > 10_000 {
        return Err("COUNT must be between 1 and 10000".into());
    }
    let label = &args[3];
    if label.is_empty() {
        return Err("SEED_LABEL must be nonempty".into());
    }
    let out = PathBuf::from(&args[4]);
    if out.exists() {
        return Err(format!(
            "refusing to replace existing directory {}",
            out.display()
        ));
    }
    let curve = KoblitzCurve::new(0, n).ok_or("unsupported Koblitz instance")?;
    let mut fixtures = Vec::with_capacity(count as usize);
    let mut scalars = HashSet::new();
    let mut points = HashSet::new();
    for i in 1..=count {
        let record = fixture(&curve, label, i)?;
        let scalar = record["fixture_scalar"].as_str().ok_or("missing scalar")?;
        let point = record["public_point"].as_array().ok_or("missing point")?;
        if !scalars.insert(scalar.to_string()) || !points.insert(point.to_vec()) {
            return Err("duplicate scalar or point; use another seed label".into());
        }
        fixtures.push(record);
    }
    let generator = match curve.generator() {
        BinaryPoint::Affine { x, y } => [
            x.to_biguint()
                .to_u64()
                .ok_or("generator x does not fit u64")?,
            y.to_biguint()
                .to_u64()
                .ok_or("generator y does not fit u64")?,
        ],
        BinaryPoint::Infinity => return Err("generator is identity".into()),
    };
    fs::create_dir(&out).map_err(|e| e.to_string())?;
    for record in &fixtures {
        let name = record["label"].as_str().ok_or("missing label")?;
        let folder = out.join(name);
        fs::create_dir(&folder).map_err(|e| e.to_string())?;
        let bytes = serde_json::to_vec(&record["public_point"]).map_err(|e| e.to_string())?;
        write_new(
            &folder.join("public_target.json"),
            &[bytes, vec![b'\n']].concat(),
        )?;
    }
    let manifest = json!({
        "kind": "s3_pair_public_target_fixtures_v1",
        "field_degree": n,
        "koblitz_a": 0,
        "subgroup_order": curve.subgroup_order.to_string(),
        "generator": generator,
        "scalar_law": "1 + int_be(SHA256('s3-pair-target-v1:' || seed_label || ':' || u32be(n) || u32be(index))) mod (r-1)",
        "seed_label": label,
        "count": count,
        "fixtures": fixtures,
        "solver_input_policy": "only each public_target.json point is passed to the solver; fixture_scalar is for independent checking after measurement"
    });
    let mut bytes = serde_json::to_vec_pretty(&manifest).map_err(|e| e.to_string())?;
    bytes.push(b'\n');
    write_new(&out.join("fixtures.json"), &bytes)?;
    println!("{}", out.display());
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("s3_pilot_targets: {error}");
        std::process::exit(2);
    }
}
