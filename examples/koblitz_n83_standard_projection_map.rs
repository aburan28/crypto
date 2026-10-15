//! Exact source-point to projected signed-Frobenius column map for N83.
//! Usage: construct NEW_JSONL DIMENSION | verify JSONL DIMENSION
use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, point_key, KoblitzCurve,
    ProjectedFactorBase,
};
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fs::OpenOptions;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
type Key = (BigUint, BigUint);
const CURVE: &str = "icv1-f2m83-tm6151469093347-debefd74";
const SOURCE: &[u8] = include_bytes!("koblitz_n83_standard_projection_map.rs");

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Location {
    column: usize,
    phase: u32,
    negated: bool,
}

struct ProjectionMap {
    curve: KoblitzCurve,
    source: crypto_lib::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    projected: crypto_lib::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    representatives: Vec<BinaryPoint>,
    locations: Vec<Location>,
    projected_indices: Vec<Option<usize>>,
}

fn point_json(point: &BinaryPoint) -> Result<Value> {
    match point {
        BinaryPoint::Infinity => Err("factor-base row is identity".into()),
        BinaryPoint::Affine { x, y } => Ok(json!({
            "x": x.to_biguint().to_string(), "y": y.to_biguint().to_string(),
        })),
    }
}

fn point_digest(points: &[BinaryPoint]) -> Result<String> {
    let mut keys: Vec<Key> = points.iter().map(point_key).collect();
    keys.sort_unstable();
    if keys.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err("duplicate point in factor base".into());
    }
    let mut hasher = blake3::Hasher::new();
    for (x_plus_one, y) in keys {
        hasher.update((x_plus_one - BigUint::from(1u8)).to_string().as_bytes());
        hasher.update(b":");
        hasher.update(y.to_string().as_bytes());
        hasher.update(b"\n");
    }
    Ok(hasher.finalize().to_hex().to_string())
}

/// If R = (-1)^s pi^k(Q) is canonical, Q = (-1)^s pi^(n-k)(R).
fn canonical(curve: &KoblitzCurve, point: &BinaryPoint) -> (Key, BinaryPoint, u32, bool) {
    let mut current = point.clone();
    let mut best_key = point_key(point);
    let mut best_point = point.clone();
    let mut best_forward = 0;
    let mut best_negated = false;
    for k in 0..curve.n {
        for (negated, candidate) in [(false, current.clone()), (true, point_neg(&current))] {
            let key = point_key(&candidate);
            if key < best_key {
                best_key = key;
                best_point = candidate;
                best_forward = k;
                best_negated = negated;
            }
        }
        current = curve.frobenius(&current);
    }
    (
        best_key,
        best_point,
        (curve.n - best_forward) % curve.n,
        best_negated,
    )
}

fn build(dimension: u32) -> Result<ProjectionMap> {
    if !(2..=18).contains(&dimension) {
        return Err("dimension must be 2..=18; dimension 1 has empty cofactor image".into());
    }
    let curve = KoblitzCurve::known_n83_k0().ok_or("pinned K0 curve unavailable")?;
    if curve.label() != CURVE || curve.n != 83 || curve.cofactor != BigUint::from(4u8) {
        return Err("curve identity or cofactor changed".into());
    }
    let source = build_standard_subspace_factor_base(&curve, dimension)?;
    let projected = cofactor_project_factor_base(&curve, &source)?;
    if projected.points.len() != 2 * projected.signed_orbits.len() {
        return Err("projected sign orbit count changed".into());
    }
    let mut prelim: Vec<Option<(Key, u32, bool)>> = vec![None; projected.points.len()];
    let mut representatives = BTreeMap::<Key, BinaryPoint>::new();
    for orbit in &projected.signed_orbits {
        if orbit.len() != 2 {
            return Err("projected sign orbit is not a pair".into());
        }
        let first = orbit[0];
        let second = orbit[1];
        if point_neg(&projected.points[first]) != projected.points[second] {
            return Err("projected sign pair is malformed".into());
        }
        let (key, canonical_point, phase, negated) = canonical(&curve, &projected.points[first]);
        representatives
            .entry(key.clone())
            .or_insert(canonical_point);
        prelim[first] = Some((key.clone(), phase, negated));
        prelim[second] = Some((key, phase, !negated));
    }
    let column_of: BTreeMap<Key, usize> = representatives
        .keys()
        .cloned()
        .enumerate()
        .map(|(index, key)| (key, index))
        .collect();
    let locations: Vec<Location> = prelim
        .into_iter()
        .map(|item| {
            let (key, phase, negated) = item.ok_or("unmapped projected point")?;
            Ok(Location {
                column: *column_of.get(&key).ok_or("missing canonical column")?,
                phase,
                negated,
            })
        })
        .collect::<Result<_>>()?;
    let index_of: BTreeMap<Key, usize> = projected
        .points
        .iter()
        .enumerate()
        .map(|(index, point)| (point_key(point), index))
        .collect();
    if index_of.len() != projected.points.len() {
        return Err("duplicate projected point".into());
    }
    let projected_indices: Vec<Option<usize>> = source
        .points
        .iter()
        .map(|point| {
            let image = curve.mul(point, &curve.cofactor);
            if image == BinaryPoint::Infinity {
                Ok(None)
            } else {
                Ok(Some(
                    *index_of
                        .get(&point_key(&image))
                        .ok_or("source projection missing")?,
                ))
            }
        })
        .collect::<Result<_>>()?;
    if dimension == 12 {
        let reference = ProjectedFactorBase::new(&curve, &source);
        if reference.columns() != column_of.len() || reference.columns() != 2_001 {
            return Err("dimension-12 existing projection-map control changed".into());
        }
        for (source_index, projected_index) in projected_indices.iter().enumerate() {
            if reference
                .location_of(source_index)
                .map(|location| location.0)
                != projected_index.map(|index| locations[index].column)
            {
                return Err("dimension-12 source column differs from existing map".into());
            }
        }
    }
    Ok(ProjectionMap {
        curve,
        source,
        projected,
        representatives: representatives.into_values().collect(),
        locations,
        projected_indices,
    })
}

fn header(map: &ProjectionMap, dimension: u32) -> Result<Value> {
    Ok(json!({
        "schema": "n83.standard-subspace-projection-map/v1",
        "curve_slug": CURVE,
        "source_dimension": dimension,
        "source_x_rule": "polynomial_basis_low_bits_0_through_dimension_minus_1",
        "source_point_count": map.source.points.len(),
        "source_point_set_blake3": point_digest(&map.source.points)?,
        "projected_point_count": map.projected.points.len(),
        "projected_point_set_blake3": point_digest(&map.projected.points)?,
        "sign_only_columns": map.projected.signed_orbits.len(),
        "projected_signed_frobenius_columns": map.representatives.len(),
        "excluded_source_points": map.projected_indices.iter().filter(|index| index.is_none()).count(),
        "subgroup_order": map.curve.subgroup_order.to_string(),
        "cofactor": map.curve.cofactor.to_string(),
        "frobenius_lambda": map.curve.lambda.to_string(),
        "source_blake3": blake3::hash(SOURCE).to_hex().to_string(),
        "selected_best_total_runtime": Value::Null,
    }))
}

fn column_row(index: usize, point: &BinaryPoint) -> Result<Value> {
    Ok(json!({"type": "column", "index": index, "representative": point_json(point)?}))
}

fn projected_row(map: &ProjectionMap, index: usize) -> Result<Value> {
    let location = map.locations[index];
    Ok(json!({
        "type": "projected", "index": index,
        "point": point_json(&map.projected.points[index])?,
        "sign_only_column": map.projected.signed_orbit_of[index].0,
        "sign_only_negated": map.projected.signed_orbit_of[index].2,
        "column": location.column, "phase": location.phase, "negated": location.negated,
    }))
}

fn source_row(map: &ProjectionMap, index: usize) -> Result<Value> {
    let point = &map.source.points[index];
    let BinaryPoint::Affine { x, .. } = point else {
        return Err("source point is identity".into());
    };
    let projected_index = map.projected_indices[index];
    let location = projected_index.map(|i| map.locations[i]);
    Ok(json!({
        "type": "source", "index": index,
        "point": point_json(point)?, "source_x_code": x.to_biguint().to_string(),
        "projected_index": projected_index,
        "column": location.map(|l| l.column),
        "phase": location.map(|l| l.phase),
        "negated": location.map(|l| l.negated),
    }))
}

fn rows(map: &ProjectionMap, dimension: u32) -> Result<impl Iterator<Item = Result<Value>> + '_> {
    let header = std::iter::once(Ok(header(map, dimension)?));
    let columns = map
        .representatives
        .iter()
        .enumerate()
        .map(|(i, p)| column_row(i, p));
    let projected = (0..map.projected.points.len()).map(|i| projected_row(map, i));
    let source = (0..map.source.points.len()).map(|i| source_row(map, i));
    Ok(header.chain(columns).chain(projected).chain(source))
}

fn construct(path: &Path, dimension: u32) -> Result<Value> {
    let map = build(dimension)?;
    let file = OpenOptions::new().write(true).create_new(true).open(path)?;
    let mut stream = BufWriter::new(file);
    for row in rows(&map, dimension)? {
        serde_json::to_writer(&mut stream, &row?)?;
        stream.write_all(b"\n")?;
    }
    stream.flush()?;
    verify(path, dimension)
}

fn verify(path: &Path, dimension: u32) -> Result<Value> {
    let map = build(dimension)?;
    let mut actual = BufReader::new(std::fs::File::open(path)?).lines();
    for expected in rows(&map, dimension)? {
        let line = actual.next().ok_or("missing map row")??;
        let parsed: Value = serde_json::from_str(&line)?;
        if parsed != expected? {
            return Err("map row differs from exact construction".into());
        }
    }
    if actual.next().is_some() {
        return Err("extra map row".into());
    }
    Ok(json!({
        "status": "PASS_EXACT_SOURCE_TO_PROJECTED_COLUMN_MAP",
        "source_dimension": dimension,
        "source_point_count": map.source.points.len(),
        "projected_point_count": map.projected.points.len(),
        "projected_signed_frobenius_columns": map.representatives.len(),
        "excluded_source_points": map.projected_indices.iter().filter(|index| index.is_none()).count(),
    }))
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let result = match args.as_slice() {
        [_, mode, path, dimension] if mode == "construct" => {
            construct(Path::new(path), dimension.parse()?)?
        }
        [_, mode, path, dimension] if mode == "verify" => {
            verify(Path::new(path), dimension.parse()?)?
        }
        _ => return Err("usage: construct NEW_JSONL DIMENSION | verify JSONL DIMENSION".into()),
    };
    println!("{}", serde_json::to_string(&result)?);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_labels_reconstruct_projected_points() {
        let map = build(4).unwrap();
        let mut columns = BTreeSet::new();
        for (point, location) in map.projected.points.iter().zip(&map.locations) {
            let mut reconstructed = map.representatives[location.column].clone();
            for _ in 0..location.phase {
                reconstructed = map.curve.frobenius(&reconstructed);
            }
            if location.negated {
                reconstructed = point_neg(&reconstructed);
            }
            assert_eq!(&reconstructed, point);
            columns.insert(location.column);
        }
        assert_eq!(columns.len(), map.representatives.len());
    }
}
