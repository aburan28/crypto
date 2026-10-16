//! Export and replay a typed, public N83 polynomial-subspace factor base.
//! Usage: construct NEW_DIRECTORY DIMENSION | verify DIRECTORY
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, FrobeniusFactorBase,
    KoblitzCurve,
};
use flate2::{read::GzDecoder, Compression, GzBuilder};
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::error::Error;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::process::Command;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const CURVE: &str = "icv1-f2m83-tm6151469093347-debefd74";
const PREFIX: &str = "s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/standard-subspace/a0";
const SOURCE: &[u8] = include_bytes!("koblitz_n83_standard_subspace_export.rs");

fn point_json(point: &BinaryPoint) -> Result<Value> {
    match point {
        BinaryPoint::Infinity => Err("factor-base row is identity".into()),
        BinaryPoint::Affine { x, y } => Ok(json!({
            "x": x.to_biguint().to_string(), "y": y.to_biguint().to_string(),
        })),
    }
}

fn point_digest(base: &FrobeniusFactorBase) -> Result<String> {
    let mut pairs = Vec::with_capacity(base.points.len());
    for point in &base.points {
        let BinaryPoint::Affine { x, y } = point else {
            return Err("factor-base point is identity".into());
        };
        pairs.push((x.to_biguint(), y.to_biguint()));
    }
    pairs.sort_unstable();
    if pairs.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err("duplicate factor-base point".into());
    }
    let mut hash = blake3::Hasher::new();
    for (x, y) in pairs {
        hash.update(x.to_string().as_bytes());
        hash.update(b":");
        hash.update(y.to_string().as_bytes());
        hash.update(b"\n");
    }
    Ok(hash.finalize().to_hex().to_string())
}

fn build(dimension: u32) -> Result<(KoblitzCurve, FrobeniusFactorBase, FrobeniusFactorBase)> {
    if !(1..=18).contains(&dimension) {
        return Err("dimension must be 1..=18".into());
    }
    let curve = KoblitzCurve::known_n83_k0().ok_or("pinned K0 curve unavailable")?;
    if curve.label() != CURVE || curve.cofactor != BigUint::from(4u8) {
        return Err("curve identity or cofactor changed".into());
    }
    let source = build_standard_subspace_factor_base(&curve, dimension)?;
    let projected = cofactor_project_factor_base(&curve, &source)?;
    if dimension == 18
        && (
            source.points.len(),
            projected.points.len(),
            projected.signed_orbits.len(),
        ) != (261_447, 261_444, 130_722)
    {
        return Err("dimension-18 prior-art counts changed".into());
    }
    if projected.signed_orbit_of.len() != projected.points.len()
        || projected.signed_orbits.len() * 2 != projected.points.len()
        || projected
            .signed_orbit_of
            .iter()
            .any(|&(_, phase, _)| phase != 0)
    {
        return Err("projected base is not an exact sign-folded point set".into());
    }
    Ok((curve, source, projected))
}

fn header(
    curve: &KoblitzCurve,
    source: &FrobeniusFactorBase,
    projected: &FrobeniusFactorBase,
    dimension: u32,
) -> Result<Value> {
    Ok(json!({
        "schema": "n83.standard-subspace-factor-base/v1",
        "curve_slug": CURVE,
        "field_degree": 83,
        "modulus_low_terms": [0, 1, 2, 45],
        "curve_a": 0, "curve_b": "1",
        "subgroup_order": curve.subgroup_order.to_string(),
        "cofactor": curve.cofactor.to_string(),
        "generator": point_json(curve.generator())?,
        "source_dimension": dimension,
        "source_x_rule": "polynomial_basis_low_bits_0_through_dimension_minus_1",
        "source_point_count": source.points.len(),
        "source_point_set_blake3": point_digest(source)?,
        "projected_point_count": projected.points.len(),
        "projected_point_set_blake3": point_digest(projected)?,
        "signed_columns": projected.signed_orbits.len(),
        "closure": "sign_only_after_public_cofactor_projection",
        "subspace_membership_after_projection": false,
        "source_blake3": blake3::hash(SOURCE).to_hex().to_string(),
        "selected_best_total_runtime": Value::Null,
    }))
}

fn rows(projected: &FrobeniusFactorBase) -> Result<Vec<Value>> {
    let mut out = Vec::with_capacity(projected.points.len());
    for (index, point) in projected.points.iter().enumerate() {
        let (column, phase, negated) = projected.signed_orbit_of[index];
        if phase != 0 || column >= projected.signed_orbits.len() {
            return Err("invalid sign-only label".into());
        }
        out.push(json!({
            "index": index, "point": point_json(point)?,
            "signed_column": column, "negated": negated,
        }));
    }
    Ok(out)
}

fn write_jsonl(path: &Path, header: &Value, rows: &[Value]) -> Result<()> {
    let mut stream = BufWriter::new(OpenOptions::new().write(true).create_new(true).open(path)?);
    serde_json::to_writer(&mut stream, header)?;
    stream.write_all(b"\n")?;
    for row in rows {
        serde_json::to_writer(&mut stream, row)?;
        stream.write_all(b"\n")?;
    }
    stream.flush()?;
    Ok(())
}

fn git_revision(spec: &str) -> Result<String> {
    let output = Command::new("git").args(["rev-parse", spec]).output()?;
    if !output.status.success() {
        return Err(format!("cannot pin Git revision {spec}").into());
    }
    let revision = String::from_utf8(output.stdout)?.trim().to_owned();
    if revision.len() != 40 || !revision.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(format!("invalid Git revision for {spec}").into());
    }
    Ok(revision)
}

fn construct(directory: &Path, dimension: u32) -> Result<Value> {
    fs::create_dir(directory)?;
    fs::create_dir(directory.join("objects"))?;
    let (curve, source, projected) = build(dimension)?;
    let header = header(&curve, &source, &projected, dimension)?;
    let rows = rows(&projected)?;
    let plain = directory.join("factor-base.jsonl");
    write_jsonl(&plain, &header, &rows)?;
    let plain_bytes = fs::read(&plain)?;
    let compressed_temp = directory.join("objects/.partial.jsonl.gz");
    let output = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&compressed_temp)?;
    let mut gzip = GzBuilder::new()
        .mtime(0)
        .write(output, Compression::default());
    std::io::copy(&mut File::open(&plain)?, &mut gzip)?;
    gzip.finish()?;
    let compressed_bytes = fs::read(&compressed_temp)?;
    let compressed_hash = blake3::hash(&compressed_bytes).to_hex().to_string();
    let object = format!("objects/{compressed_hash}.jsonl.gz");
    fs::rename(compressed_temp, directory.join(&object))?;
    let manifest = json!({
        "schema": "n83.standard-subspace-factor-base-manifest/v1",
        "curve_slug": CURVE,
        "source_dimension": dimension,
        "source_blake3": blake3::hash(SOURCE).to_hex().to_string(),
        "source_commit": git_revision("HEAD")?,
        "src_git_tree": git_revision("HEAD:src")?,
        "cargo_manifest_blake3": blake3::hash(include_bytes!("../Cargo.toml")).to_hex().to_string(),
        "cargo_lock_blake3": blake3::hash(include_bytes!("../Cargo.lock")).to_hex().to_string(),
        "object": object,
        "s3_uri": format!("{PREFIX}/{object}"),
        "compressed_blake3": compressed_hash,
        "compressed_bytes": compressed_bytes.len(),
        "plain_blake3": blake3::hash(&plain_bytes).to_hex().to_string(),
        "plain_bytes": plain_bytes.len(),
        "source_point_count": source.points.len(),
        "projected_point_count": projected.points.len(),
        "signed_columns": projected.signed_orbits.len(),
        "projected_point_set_blake3": header["projected_point_set_blake3"],
        "selected_best_total_runtime": Value::Null,
    });
    fs::write(
        directory.join("manifest.json"),
        format!("{}\n", serde_json::to_string_pretty(&manifest)?),
    )?;
    verify(directory)?;
    Ok(manifest)
}

fn verify(directory: &Path) -> Result<Value> {
    let manifest: Value = serde_json::from_slice(&fs::read(directory.join("manifest.json"))?)?;
    let dimension = manifest["source_dimension"]
        .as_u64()
        .ok_or("missing source dimension")? as u32;
    let (curve, source, projected) = build(dimension)?;
    let expected_header = header(&curve, &source, &projected, dimension)?;
    let expected_rows = rows(&projected)?;
    let object = manifest["object"].as_str().ok_or("missing object path")?;
    let object_hash = manifest["compressed_blake3"]
        .as_str()
        .ok_or("missing compressed hash")?;
    let source_commit = manifest["source_commit"]
        .as_str()
        .ok_or("missing source commit")?;
    if object != format!("objects/{object_hash}.jsonl.gz")
        || manifest["s3_uri"] != format!("{PREFIX}/{object}")
        || manifest["schema"] != "n83.standard-subspace-factor-base-manifest/v1"
        || manifest["curve_slug"] != CURVE
        || manifest["source_blake3"] != blake3::hash(SOURCE).to_hex().to_string()
        || manifest["src_git_tree"] != git_revision(&format!("{source_commit}:src"))?
        || manifest["cargo_manifest_blake3"]
            != blake3::hash(include_bytes!("../Cargo.toml"))
                .to_hex()
                .to_string()
        || manifest["cargo_lock_blake3"]
            != blake3::hash(include_bytes!("../Cargo.lock"))
                .to_hex()
                .to_string()
    {
        return Err("manifest identity mismatch".into());
    }
    let compressed = fs::read(directory.join(object))?;
    if blake3::hash(&compressed).to_hex().to_string() != object_hash
        || manifest["compressed_bytes"] != compressed.len()
    {
        return Err("compressed object mismatch".into());
    }
    let mut plain = Vec::new();
    GzDecoder::new(&compressed[..]).read_to_end(&mut plain)?;
    if manifest["plain_blake3"] != blake3::hash(&plain).to_hex().to_string()
        || manifest["plain_bytes"] != plain.len()
    {
        return Err("plain object mismatch".into());
    }
    let optional_plain = directory.join("factor-base.jsonl");
    if optional_plain.exists() && fs::read(optional_plain)? != plain {
        return Err("local plain copy mismatch".into());
    }
    let mut lines = BufReader::new(&plain[..]).lines();
    let actual_header: Value = serde_json::from_str(&lines.next().ok_or("missing header")??)?;
    if actual_header != expected_header
        || manifest["source_point_count"] != source.points.len()
        || manifest["projected_point_count"] != projected.points.len()
        || manifest["signed_columns"] != projected.signed_orbits.len()
        || manifest["projected_point_set_blake3"] != expected_header["projected_point_set_blake3"]
    {
        return Err("header or construction replay mismatch".into());
    }
    for (index, expected) in expected_rows.iter().enumerate() {
        let actual: Value = serde_json::from_str(&lines.next().ok_or("missing point row")??)?;
        if actual != *expected || actual["index"] != index {
            return Err("point or sign-column replay mismatch".into());
        }
    }
    if lines.next().is_some() {
        return Err("extra point row".into());
    }
    Ok(json!({"status": "PASS_EXACT_STANDARD_SUBSPACE_REPLAY",
              "source_dimension": dimension,
              "projected_point_count": projected.points.len(),
              "signed_columns": projected.signed_orbits.len(),
              "object": object,
              "compressed_blake3": object_hash}))
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let output = match args.as_slice() {
        [_, mode, directory, dimension] if mode == "construct" => {
            construct(Path::new(directory), dimension.parse()?)?
        }
        [_, mode, directory] if mode == "verify" => verify(Path::new(directory))?,
        _ => return Err("usage: construct NEW_DIRECTORY DIMENSION | verify DIRECTORY".into()),
    };
    println!("{}", serde_json::to_string(&output)?);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sign_only_archive_round_trip_and_tamper_rejection() {
        let directory = std::env::temp_dir().join(format!(
            "n83-standard-subspace-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("archive")
        ));
        if directory.exists() {
            fs::remove_dir_all(&directory).unwrap();
        }
        let manifest = construct(&directory, 4).unwrap();
        assert_eq!(
            verify(&directory).unwrap()["status"],
            "PASS_EXACT_STANDARD_SUBSPACE_REPLAY"
        );
        fs::remove_file(directory.join("factor-base.jsonl")).unwrap();
        assert_eq!(
            verify(&directory).unwrap()["status"],
            "PASS_EXACT_STANDARD_SUBSPACE_REPLAY"
        );
        let object = directory.join(manifest["object"].as_str().unwrap());
        let mut bytes = fs::read(&object).unwrap();
        let last = bytes.len() - 1;
        bytes[last] ^= 1;
        fs::write(&object, bytes).unwrap();
        assert!(verify(&directory).is_err());
        fs::remove_dir_all(directory).unwrap();
    }
}
