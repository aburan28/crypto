//! N83 dimension-12 source subspace, projected signed-Frobenius closure.
//! Protocol: research/koblitz_n83_factor_base_sweep_20261008/verification/structured-d12-orbit-20261010/PROTOCOL.md.
use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_fast_arith::FastBinaryCurve128;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::error::Error;
use std::fs::{self, OpenOptions};
use std::io::{Read, Write};
use std::path::Path;
use std::process::Command;
use std::time::Instant;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
type Words = (u128, u128);
const STUDY: &str = "koblitz_n83_factor_base_sweep_20261008";
const CURVE: &str = "icv1-f2m83-tm6151469093347-debefd74";
const OLD_HASH: &str = "db3e7f25877c75dfd6f53db46447972ce6a44434cb7ea1e53e450d21af0fe79c";
const S3: &str = "s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v3-structured-orbit/a0";
const DESIGN: &[u8] =
    include_bytes!("../research/koblitz_n83_factor_base_sweep_20261008/structured-orbit-v3.json");
const SOURCE: &[u8] = include_bytes!("koblitz_n83_structured_orbit_export.rs");

fn tagged((x, y): Words) -> [u8; 33] {
    let mut row = [0u8; 33];
    row[0] = 1;
    row[1..17].copy_from_slice(&x.to_le_bytes());
    row[17..33].copy_from_slice(&y.to_le_bytes());
    row
}

fn point_json(p: Words) -> Value {
    json!([p.0.to_string(), p.1.to_string()])
}

fn parse_point(v: &Value) -> Result<Words> {
    let row = v.as_array().ok_or("point array")?;
    if row.len() != 2 {
        return Err("point arity".into());
    }
    let x: u128 = row[0].as_str().ok_or("x string")?.parse()?;
    let y: u128 = row[1].as_str().ok_or("y string")?.parse()?;
    if x >= 1u128 << 83 || y >= 1u128 << 83 {
        return Err("out of field".into());
    }
    Ok((x, y))
}

fn generic((x, y): Words) -> BinaryPoint {
    BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(x), 83),
        y: F2mElement::from_biguint(&BigUint::from(y), 83),
    }
}

fn set_hashes(points: &HashSet<Words>) -> (String, String) {
    let mut old_rows: Vec<_> = points.iter().copied().map(tagged).collect();
    old_rows.sort_unstable();
    let mut old = blake3::Hasher::new();
    for row in old_rows {
        old.update(&row);
    }
    let mut new_rows: Vec<_> = points.iter().copied().collect();
    new_rows.sort_unstable();
    let mut exporter = blake3::Hasher::new();
    exporter.update(b"n83.sorted-polynomial-points/v1");
    for (x, y) in new_rows {
        exporter.update(&x.to_le_bytes());
        exporter.update(&y.to_le_bytes());
    }
    (
        old.finalize().to_hex().to_string(),
        exporter.finalize().to_hex().to_string(),
    )
}

fn select(kc: &KoblitzCurve, fast: &FastBinaryCurve128) -> Result<Vec<(u16, Words, Words)>> {
    let mut signed_seeds = HashSet::new();
    let mut orbit_keys = HashSet::new();
    let mut chosen = Vec::new();
    for x in 0..4096u128 {
        let Some(source) = fast.points_with_x(1, x).first().copied().flatten() else {
            continue;
        };
        let Some(rep) = fast.scalar_mul(Some(source), &kc.cofactor) else {
            continue;
        };
        signed_seeds.insert(tagged(rep).min(tagged((rep.0, rep.0 ^ rep.1))));
        let mut p = rep;
        let mut key = [255u8; 33];
        for _ in 0..83 {
            key = key.min(tagged(p)).min(tagged((p.0, p.0 ^ p.1)));
            p = (fast.gf.sqr(p.0), fast.gf.sqr(p.1));
        }
        if p != rep {
            return Err("Frobenius did not close".into());
        }
        if orbit_keys.insert(key) {
            chosen.push((x as u16, source, rep));
        }
    }
    if signed_seeds.len() != 2027 || chosen.len() != 2001 {
        return Err("historical orbit count".into());
    }
    Ok(chosen)
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<()> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.write_all(bytes)?;
    file.sync_all()?;
    Ok(())
}

fn construct(root: &Path) -> Result<()> {
    let start = Instant::now();
    let kc = KoblitzCurve::known_n83_k0().ok_or("pinned K0")?;
    if kc.label() != CURVE || kc.cofactor != BigUint::from(4u8) {
        return Err("curve binding".into());
    }
    let fast = FastBinaryCurve128::new(&kc.curve.irreducible, 0).ok_or("fast curve")?;
    let chosen = select(&kc, &fast)?;
    let header = json!({
        "schema":"n83.structured-orbit-base/v1", "study":STUDY,
        "curve_slug":CURVE, "n":83, "curve_a":0, "curve_b":"1",
        "modulus_low_terms":[0,1,2,45], "subgroup_order":kc.subgroup_order.to_string(),
        "cofactor":"4", "frobenius_lambda":kc.lambda.to_string(),
        "source_dimension":12, "closure":"signed_frobenius",
        "source_x_codes":chosen.iter().map(|r| r.0).collect::<Vec<_>>(),
        "source_points_before_cofactor_clearing":chosen.iter().map(|r| point_json(r.1)).collect::<Vec<_>>(),
        "representatives":chosen.iter().map(|r| point_json(r.2)).collect::<Vec<_>>(),
        "orbit_columns":2001, "point_count":332166,
        "historical_set_blake3":OLD_HASH
    });
    let mut plain = Vec::new();
    serde_json::to_writer(&mut plain, &header)?;
    plain.write_all(b"\n")?;
    let mut points = HashSet::new();
    for (column, (_, _, rep)) in chosen.iter().enumerate() {
        let mut p = *rep;
        let mut coeff = BigUint::from(1u8);
        for phase in 0..83 {
            for negative in [false, true] {
                let q = if negative { (p.0, p.0 ^ p.1) } else { p };
                if !points.insert(q) {
                    return Err("duplicate point".into());
                }
                let label = if negative {
                    &kc.subgroup_order - &coeff
                } else {
                    coeff.clone()
                };
                serde_json::to_writer(
                    &mut plain,
                    &json!({
                        "point":point_json(q), "column":column, "phase":phase,
                        "negative":negative, "coefficient":label.to_string()
                    }),
                )?;
                plain.write_all(b"\n")?;
            }
            p = (fast.gf.sqr(p.0), fast.gf.sqr(p.1));
            coeff = (&coeff * &kc.lambda) % &kc.subgroup_order;
        }
        if p != *rep {
            return Err("orbit did not close".into());
        }
    }
    if points.len() != 332166 {
        return Err("point count".into());
    }
    let (old_hash, point_set_hash) = set_hashes(&points);
    if old_hash != OLD_HASH {
        return Err("historical set hash".into());
    }
    let plain_hash = blake3::hash(&plain).to_hex().to_string();
    let mut gzip = GzEncoder::new(Vec::new(), Compression::default());
    gzip.write_all(&plain)?;
    let compressed = gzip.finish()?;
    let hash = blake3::hash(&compressed).to_hex().to_string();
    let object = format!("objects/{hash}.jsonl.gz");
    fs::create_dir(root)?;
    fs::create_dir(root.join("objects"))?;
    write_new(&root.join(&object), &compressed)?;
    let git = Command::new("git").args(["rev-parse", "HEAD"]).output()?;
    if !git.status.success() {
        return Err("source commit".into());
    }
    let commit = String::from_utf8(git.stdout)?.trim().to_owned();
    let manifest = json!({
        "schema":"n83.structured-orbit-panel/v1", "study":STUDY,
        "status":"completed_factor_base_object", "source_commit":commit,
        "source_blake3":blake3::hash(SOURCE).to_hex().to_string(),
        "design_blake3":blake3::hash(DESIGN).to_hex().to_string(),
        "object":object, "compressed_blake3":hash, "plain_blake3":plain_hash,
        "compressed_bytes":compressed.len(), "point_set_blake3":point_set_hash,
        "historical_set_blake3":old_hash, "orbit_columns":chosen.len(),
        "point_records":points.len(), "s3_uri":format!("{S3}/{object}"),
        "construction_elapsed_ms_L0":start.elapsed().as_secs_f64()*1000.0,
        "relation_yield":Value::Null, "matrix_rank":Value::Null,
        "total_index_calculus_runtime_ms":Value::Null,
        "selected_best_total_runtime":Value::Null
    });
    let mut manifest_bytes = serde_json::to_vec_pretty(&manifest)?;
    manifest_bytes.push(b'\n');
    write_new(&root.join("manifest.json"), &manifest_bytes)?;
    println!("{manifest}");
    Ok(())
}

fn replay(root: &Path) -> Result<()> {
    let start = Instant::now();
    let manifest_bytes = fs::read(root.join("manifest.json"))?;
    let manifest: Value = serde_json::from_slice(&manifest_bytes)?;
    if manifest["schema"] != "n83.structured-orbit-panel/v1"
        || manifest["status"] != "completed_factor_base_object"
        || manifest["source_blake3"] != blake3::hash(SOURCE).to_hex().as_str()
        || manifest["design_blake3"] != blake3::hash(DESIGN).to_hex().as_str()
    {
        return Err("manifest/source/design binding".into());
    }
    let object = manifest["object"].as_str().ok_or("object")?;
    let hash = manifest["compressed_blake3"].as_str().ok_or("hash")?;
    if object != format!("objects/{hash}.jsonl.gz") {
        return Err("content path".into());
    }
    let compressed = fs::read(root.join(object))?;
    if blake3::hash(&compressed).to_hex().as_str() != hash
        || manifest["compressed_bytes"].as_u64() != Some(compressed.len() as u64)
    {
        return Err("compressed object".into());
    }
    let mut plain = Vec::new();
    GzDecoder::new(&compressed[..]).read_to_end(&mut plain)?;
    if manifest["plain_blake3"] != blake3::hash(&plain).to_hex().as_str() {
        return Err("plain hash".into());
    }
    let text = std::str::from_utf8(&plain)?;
    let mut lines = text.lines();
    let header: Value = serde_json::from_str(lines.next().ok_or("header")?)?;
    let kc = KoblitzCurve::known_n83_k0().ok_or("pinned K0")?;
    if header["schema"] != "n83.structured-orbit-base/v1"
        || header["study"] != STUDY
        || header["curve_slug"] != CURVE
        || header["modulus_low_terms"] != json!([0, 1, 2, 45])
        || header["subgroup_order"] != kc.subgroup_order.to_string()
        || header["cofactor"] != "4"
        || header["frobenius_lambda"] != kc.lambda.to_string()
        || header["source_dimension"] != 12
        || header["orbit_columns"] != 2001
        || header["point_count"] != 332166
        || header["historical_set_blake3"] != OLD_HASH
    {
        return Err("header/curve binding".into());
    }
    let reps = header["representatives"]
        .as_array()
        .ok_or("representatives")?;
    let sources = header["source_points_before_cofactor_clearing"]
        .as_array()
        .ok_or("sources")?;
    let codes = header["source_x_codes"].as_array().ok_or("codes")?;
    if reps.len() != 2001 || sources.len() != reps.len() || codes.len() != reps.len() {
        return Err("representative lengths".into());
    }
    let mut seen = HashSet::new();
    let mut previous = None;
    for (column, ((rep, source), code)) in reps.iter().zip(sources).zip(codes).enumerate() {
        let r = parse_point(rep)?;
        let s = parse_point(source)?;
        let x = code.as_u64().ok_or("source code")?;
        if x >= 4096 || previous.is_some_and(|p| p >= x) || s.0 != u128::from(x) {
            return Err("source x order/domain".into());
        }
        previous = Some(x);
        let rg = generic(r);
        let sg = generic(s);
        if !kc.curve.is_on_curve(&sg)
            || !kc.curve.is_on_curve(&rg)
            || kc.mul(&sg, &kc.cofactor) != rg
            || kc.mul(&rg, &kc.subgroup_order) != BinaryPoint::Infinity
            || kc.mul(&rg, &kc.lambda) != kc.frobenius(&rg)
        {
            return Err("generic source/projection/subgroup/eigenvalue".into());
        }
        let mut expected = rg;
        let mut coeff = BigUint::from(1u8);
        for phase in 0..83 {
            for negative in [false, true] {
                let entry: Value = serde_json::from_str(lines.next().ok_or("point row")?)?;
                let p = parse_point(&entry["point"])?;
                let want = if negative {
                    point_neg(&expected)
                } else {
                    expected.clone()
                };
                let label = if negative {
                    &kc.subgroup_order - &coeff
                } else {
                    coeff.clone()
                };
                if generic(p) != want
                    || !seen.insert(p)
                    || entry["column"].as_u64() != Some(column as u64)
                    || entry["phase"].as_u64() != Some(phase)
                    || entry["negative"].as_bool() != Some(negative)
                    || entry["coefficient"] != label.to_string()
                {
                    return Err("generic point/orbit/label/uniqueness".into());
                }
            }
            expected = kc.frobenius(&expected);
            coeff = (&coeff * &kc.lambda) % &kc.subgroup_order;
        }
        if expected != generic(r) {
            return Err("generic orbit closure".into());
        }
    }
    if lines.next().is_some() || seen.len() != 332166 {
        return Err("trailing rows/point count".into());
    }
    let (old_hash, point_hash) = set_hashes(&seen);
    if old_hash != OLD_HASH
        || manifest["historical_set_blake3"] != old_hash
        || manifest["point_set_blake3"] != point_hash
    {
        return Err("independent set hash".into());
    }
    let receipt = json!({
        "schema":"n83.structured-orbit-replay/v1", "status":"PASS",
        "object":object, "manifest_blake3":blake3::hash(&manifest_bytes).to_hex().to_string(),
        "compressed_blake3":hash, "points_checked":seen.len(),
        "representatives_checked":reps.len(), "historical_set_blake3":old_hash,
        "backend":"generic multi-limb curve operations versus producer Gf2_128",
        "independence":"same host and repository; independent-host replay pending",
        "elapsed_ms_L0":start.elapsed().as_secs_f64()*1000.0,
        "relation_stage_executed":false, "rank_stage_executed":false,
        "total_index_calculus_runtime_ms":Value::Null
    });
    let mut bytes = serde_json::to_vec_pretty(&receipt)?;
    bytes.push(b'\n');
    write_new(&root.join("replay.json"), &bytes)?;
    println!("{receipt}");
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        return Err("usage: example construct|replay ROOT".into());
    }
    match args[1].as_str() {
        "construct" => construct(Path::new(&args[2])),
        "replay" => replay(Path::new(&args[2])),
        _ => Err("expected construct or replay".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn selected_source_orbits_match_historical_counts() {
        let kc = KoblitzCurve::known_n83_k0().unwrap();
        let fast = FastBinaryCurve128::new(&kc.curve.irreducible, 0).unwrap();
        let rows = select(&kc, &fast).unwrap();
        assert_eq!(rows.len(), 2001);
        assert!(rows.windows(2).all(|w| w[0].0 < w[1].0));
    }
}
