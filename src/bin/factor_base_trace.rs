//! Additive, source-backed export of the historical #1898 factor bases.
//! No point generation, relation search, logarithm recovery or timing run.
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::{canonical, sha256_hex, short_id};
use crypto_lib::cryptanalysis::ecbench::workload::CurveSpec;
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, GroupOps, PrimeCurve, PrimePoint};
use serde_json::{json, Value};

const SOURCE_SHA: &str = "07ee0c22dcb48969581739bdc91e4ff6c78a4ad86d52ba6792d773fa3778e248";
const SEED: &str = "factor-base-evaluation-2026-10-15-v1";
const MAGIC: &[u8; 8] = b"FBTRACE1";
const SCHEMA: &str = "ecbench.historical_factor_base_trace/v1";

type Result<T> = std::result::Result<T, String>;
fn number(v: &Value, k: &str) -> Result<u64> {
    v[k].as_u64().ok_or_else(|| format!("invalid integer {k}"))
}
fn text<'a>(v: &'a Value, k: &str) -> Result<&'a str> {
    v[k].as_str().ok_or_else(|| format!("invalid string {k}"))
}
fn read_json(p: &Path) -> Result<Value> {
    serde_json::from_slice(&fs::read(p).map_err(|e| e.to_string())?).map_err(|e| e.to_string())
}
fn write_new(p: &Path, bytes: &[u8]) -> Result<()> {
    use std::io::Write;
    let mut f = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(p)
        .map_err(|e| format!("{}: {e}", p.display()))?;
    f.write_all(bytes).map_err(|e| e.to_string())
}

// Domain-separated leaves and parents; duplicate the last node at an odd level.
// Leaves bind the curve, full FB1 digest, column and both coordinates.
fn order_root(curve: &str, fb: &str, points: &[Value]) -> Result<String> {
    if points.is_empty() {
        return Err("empty factor base".into());
    }
    let mut level = Vec::new();
    for (col, point) in points.iter().enumerate() {
        let leaf = json!({"schema": SCHEMA, "curve": curve, "fb_sha256": fb,
                          "column": col, "point": point});
        let mut bytes = vec![0];
        bytes.extend(canonical(&leaf)?.as_bytes());
        level.push(crypto_lib::hash::sha256::sha256(&bytes));
    }
    while level.len() > 1 {
        let mut next = Vec::new();
        for pair in level.chunks(2) {
            let mut bytes = vec![1];
            bytes.extend(pair[0]);
            bytes.extend(*pair.get(1).unwrap_or(&pair[0]));
            next.push(crypto_lib::hash::sha256::sha256(&bytes));
        }
        level = next;
    }
    Ok(hex::encode(level[0]))
}

fn recipe(kind: &str) -> Result<&'static str> {
    match kind {
        "coset" => Ok("Enumerate all d elements x=c*h^j mod p, j=0..d-1; lift valid x; canonicalize y; shuffle with the historical base RNG."),
        "low_x" => Ok("Scan x=0,1,...,p-1; lift valid x until B canonical classes; shuffle with the historical base RNG."),
        "random" => Ok("Historical base RNG draws randrange(p); discard previously drawn x before lifting; retain until B classes; shuffle with that same RNG."),
        "low_digit_weight" => Ok("Enumerate 8 radix positions by increasing support; coefficients -2,-1,1,2; sort each support by (L1 height, seeded digit digest, integer value); reduce modulo p, deduplicate x, lift until B; shuffle with the historical base RNG."),
        _ => Err(format!("unknown historical family {kind}")),
    }
}

// Same preimage and native keys as ecbench::methods::factor_base_facts;
// checked against the existing dump constructor in the compatibility test.
fn identify(
    slug: &str,
    family: &str,
    params: &BTreeMap<&str, String>,
    columns: usize,
    keys: &mut [u64],
) -> Result<(String, String, String)> {
    keys.sort_unstable();
    let key_bytes: Vec<u8> = keys.iter().flat_map(|k| k.to_be_bytes()).collect();
    let points_sha = sha256_hex(&key_bytes);
    let identity = json!({"schema": "ecbench.factor_base/v1", "curve": slug,
        "family": family, "params": params, "columns": columns,
        "signed_points": keys.len(), "points_sha256": points_sha});
    let (id, sha) = short_id("FB1", &identity)?;
    Ok((id, sha, points_sha))
}

// The source is an ordered list of canonical +/- representatives, not known
// multiples of a generator. Sign changes retain the original matrix column.
fn make_trace(case: &Value, source_index: usize) -> Result<(Value, Vec<u8>)> {
    let (p, b, r) = (number(case, "p")?, number(case, "b")?, number(case, "r")?);
    if !(5..(1 << 32)).contains(&p) || b >= p || !(3..(1 << 32)).contains(&r) {
        return Err("historical word-sized curve bounds violated".into());
    }
    let raw = case["base"].as_array().ok_or("base is not an array")?;
    if raw.is_empty()
        || raw.len() != number(case, "m")? as usize
        || number(case, "F")? != 2 * raw.len() as u64
    {
        return Err("historical base cardinality mismatch".into());
    }
    let curve = PrimeCurve { p, a: p - 3, b };
    let mut points = Vec::new();
    let mut cache_rows = Vec::new();
    let mut seen = BTreeSet::new();
    let mut keys = Vec::new();
    for (col, packed) in raw.iter().enumerate() {
        let packed = packed.as_u64().ok_or("invalid packed point")?;
        let (x, y) = (packed >> 32, packed & 0xffff_ffff);
        let point = PrimePoint::affine(x, y);
        if x >= p
            || y == 0
            || y > p / 2
            || !curve.is_on_curve(point)
            || !seen.insert(x)
            || !curve.mul(&mut GroupOps::default(), point, r).infinity
        {
            return Err(format!(
                "invalid, duplicate or non-subgroup point at column {col}"
            ));
        }
        points.push(json!({"x": x, "y": y, "packed": packed}));
        for (q, sign) in [(point, 1i8), (curve.neg(point), -1i8)] {
            let key = curve.key(&q);
            keys.push(key);
            cache_rows.push((key, col as u32, sign));
        }
    }
    let spec = CurveSpec::PrimeExplicit {
        p,
        a: p - 3,
        b,
        group_order: r,
        r,
        gx: number(&points[0], "x")?,
        gy: number(&points[0], "y")?,
    };
    let facts = spec.build()?.facts(&spec);
    if !facts.registered {
        return Err("register historical curve identities before exporting traces".into());
    }
    let kind = text(case, "kind")?;
    let family = format!("historical-{kind}");
    let mut params = BTreeMap::new();
    params.insert("source_sha256", SOURCE_SHA.to_string());
    params.insert("source_index", source_index.to_string());
    // FB1 binds the signed set and parameters; order is bound separately below.
    let (fb_id, fb_sha, points_sha) =
        identify(&facts.slug, &family, &params, points.len(), &mut keys)?;
    let root = order_root(&facts.slug, &fb_sha, &points)?;
    cache_rows.sort_unstable();
    let mut cache = MAGIC.to_vec();
    cache.extend(hex::decode(&fb_sha).map_err(|e| e.to_string())?);
    cache.extend(hex::decode(&root).map_err(|e| e.to_string())?);
    cache.extend((cache_rows.len() as u32).to_be_bytes());
    for (key, col, sign) in cache_rows {
        cache.extend(key.to_be_bytes());
        cache.extend(col.to_be_bytes());
        cache.push(sign as u8);
    }
    let trace = json!({"schema": SCHEMA,
        "source": {"pr": 1898, "commit": "1f7e567ab4adad212253d2b56c906c5342dbec3e",
            "archive_member": "results/generic.json", "sha256": SOURCE_SHA, "array_index": source_index,
            "producer_sha256": "3596d59c0b70b410ea6dde5b328f7491d4b45795acece61e632274054f8ecd16",
            "protocol_sha256": "07b1b92a8008d5c71b93d01155622c474372bfdebc79e5889232fa47b301dc12"},
        "case": {"bits": case["bits"], "field": case["field"], "curve": case["curve"],
                 "coset_seed": case["coset_seed"], "kind": kind},
        "curve_spec": spec, "curve": facts,
        "reference_generator": "first archived representative, chosen only for native curve validation; no generator scalar labels are known",
        "factor_base": {"fb_id": fb_id, "fb_sha256": fb_sha, "family": family,
            "params": params, "columns": points.len(), "signed_points": keys.len(),
            "points_sha256": points_sha},
        "derivation": {"seed": SEED, "recipe": recipe(kind)?, "recorded_extra": case["base_extra"],
            "rng": "CPython random.Random(int(SHA256('|'.join(seed,'base',bits,field,curve,coset_seed,kind)),16)); source evaluate.py is authoritative",
            "digit_tie_digest": "SHA256('|'.join(seed,'digit',bits,field,curve,coset_seed,positions_tuple,coefficients_tuple)) using historical Python str formatting",
            "coordinate_lift": "y^2=x^3-3*x+b mod p; retain valid x and y=min(sqrt(rhs),p-sqrt(rhs))",
            "recorded": ["final shuffled point order", "generation recipe", "seed inputs", "base_extra"],
            "not_recorded": ["complete rejected-candidate event stream", "pre-shuffle acceptance order", "per-point scalar logarithms"],
            "validation_scope": "source hash, exact ordered representatives, curve membership, [r]P=O, duplicate rejection, native ICV1/FB1 identities and cache integrity; historical RNG is not rerun"},
        "ordered_representatives": points,
        "order_merkle": {"algorithm": "SHA256; leaf=H(0x00||canonical JSON {schema,curve,fb_sha256,column,point}); parent=H(0x01||left32||right32); duplicate odd last; no leaf sorting", "root": root},
        "cache": {"format": "FBTRACE1: magic8, fb_sha32, root32, count_u32be, sorted (native_key_u64be,column_u32be,sign_i8) records", "bytes": cache.len(), "sha256": sha256_hex(&cache)},
        "claims": {"new_research_run": false, "new_performance_measurement": false,
                   "complete_historical_derivation_replay": false, "known_point_logs": false}});
    Ok((trace, cache))
}

fn source_cases(p: &Path) -> Result<Vec<Value>> {
    let bytes = fs::read(p).map_err(|e| e.to_string())?;
    if sha256_hex(&bytes) != SOURCE_SHA {
        return Err("frozen source SHA-256 mismatch".into());
    }
    let cases: Vec<Value> = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
    if cases.len() != 384 {
        return Err("expected all 384 historical cases".into());
    }
    Ok(cases)
}
fn registry_records(source: &Path, output: &Path) -> Result<()> {
    use crypto_lib::cryptanalysis::curve_id;
    use num_bigint::BigUint;
    let mut rows = BTreeMap::new();
    for case in source_cases(source)? {
        let (p, b, r) = (
            number(&case, "p")?,
            number(&case, "b")?,
            number(&case, "r")?,
        );
        let id = curve_id::prime(
            &BigUint::from(p),
            &BigUint::from(p - 3),
            &BigUint::from(b),
            &BigUint::from(r),
        )
        .ok_or("invalid historical model")?;
        let parts: Vec<_> = id.icv1.split(':').collect();
        let row = json!({"slug":id.slug,"icv1":id.icv1,"family":"prime",
            "params":{"p":p.to_string(),"a":(p-3).to_string(),"b":b.to_string()},
            "model_json":id.model_json,"trace":parts[2].parse::<i64>().map_err(|e|e.to_string())?,
            "order":parts[3],"j":parts[4],"end":parts[5],
            "aliases":[],"standard_names":[],
            "sources":["docs/curves/sources/factor-base-evaluation-20261015.json"],
            "representations":[],
            "ec1_unresolved":"no tracked record carries this model's subgroup and generator"});
        if let Some(previous) = rows.insert(id.slug, row.clone()) {
            if previous != row {
                return Err("historical model identity collision".into());
            }
        }
    }
    if rows.len() != 48 {
        return Err("expected 48 historical curve models".into());
    }
    let document = json!({"schema":"ecbench.historical_curve_registry/v1",
        "source_sha256":SOURCE_SHA,
        "provenance":"Exact p,a,b,order imported from frozen PR 1898; ICV1 arithmetic is the existing native curve_id::prime implementation; no new point count or generator label is claimed",
        "curves":rows.into_values().collect::<Vec<_>>()});
    write_new(
        output,
        &serde_json::to_vec_pretty(&document).map_err(|e| e.to_string())?,
    )?;
    println!("exported 48 native historical curve identities");
    Ok(())
}
fn catalog_entry(name: &str, trace: &Value, bytes: &[u8], cache: &[u8]) -> Value {
    json!({"name": name, "case": trace["case"], "curve": trace["curve"]["slug"],
        "fb_id": trace["factor_base"]["fb_id"], "fb_sha256": trace["factor_base"]["fb_sha256"],
        "columns": trace["factor_base"]["columns"], "trace_sha256": sha256_hex(bytes),
        "cache_sha256": sha256_hex(cache), "order_root": trace["order_merkle"]["root"]})
}
fn export(source: &Path, output: &Path) -> Result<()> {
    let cases = source_cases(source)?;
    fs::create_dir(output).map_err(|e| format!("output must be new: {e}"))?;
    let mut entries = Vec::new();
    for (i, case) in cases.iter().enumerate() {
        let (trace, cache) = make_trace(case, i)?;
        let id = text(&trace["factor_base"], "fb_id")?;
        let name = format!("{i:03}-{id}");
        let bytes = serde_json::to_vec_pretty(&trace).map_err(|e| e.to_string())?;
        write_new(&output.join(format!("{name}.json")), &bytes)?;
        write_new(&output.join(format!("{name}.cache")), &cache)?;
        entries.push(catalog_entry(&name, &trace, &bytes, &cache));
    }
    write_new(
        &output.join("catalog.json"),
        &serde_json::to_vec_pretty(&json!({
        "schema": "ecbench.historical_factor_base_catalog/v1", "source_sha256": SOURCE_SHA,
        "cases": entries.len(), "entries": entries}))
        .map_err(|e| e.to_string())?,
    )?;
    println!(
        "exported {} exact ordered bases; no research run",
        cases.len()
    );
    Ok(())
}
fn verify(source: &Path, output: &Path) -> Result<()> {
    let cases = source_cases(source)?;
    let catalog = read_json(&output.join("catalog.json"))?;
    let entries = catalog["entries"].as_array().ok_or("invalid catalog")?;
    if entries.len() != cases.len()
        || catalog["source_sha256"] != SOURCE_SHA
        || catalog["cases"] != cases.len()
        || catalog["schema"] != "ecbench.historical_factor_base_catalog/v1"
    {
        return Err("catalog source or cardinality mismatch".into());
    }
    for (i, case) in cases.iter().enumerate() {
        let (trace, cache) = make_trace(case, i)?;
        let name = format!("{i:03}-{}", text(&trace["factor_base"], "fb_id")?);
        let bytes = fs::read(output.join(format!("{name}.json"))).map_err(|e| e.to_string())?;
        if serde_json::from_slice::<Value>(&bytes).map_err(|e| e.to_string())? != trace
            || entries[i] != catalog_entry(&name, &trace, &bytes, &cache)
            || fs::read(output.join(format!("{name}.cache"))).map_err(|e| e.to_string())? != cache
        {
            return Err(format!("trace/cache mismatch at case {i}"));
        }
        let loaded = load_cache(&trace, &cache)?;
        for (col, point) in trace["ordered_representatives"]
            .as_array()
            .unwrap()
            .iter()
            .enumerate()
        {
            let (x, y) = (number(point, "x")?, number(point, "y")?);
            if loaded.query(x, y)? != Some((col as u32, 1))
                || loaded.query(x, loaded.curve.p - y)? != Some((col as u32, -1))
            {
                return Err(format!("signed lookup mismatch at case {i}, column {col}"));
            }
        }
    }
    println!(
        "verified all {} source-backed traces and caches",
        cases.len()
    );
    Ok(())
}

struct Lookup<'a> {
    curve: PrimeCurve,
    cache: &'a [u8],
    count: usize,
}
fn load_cache<'a>(trace: &Value, cache: &'a [u8]) -> Result<Lookup<'a>> {
    let spec: CurveSpec =
        serde_json::from_value(trace["curve_spec"].clone()).map_err(|e| e.to_string())?;
    let CurveSpec::PrimeExplicit { p, a, b, .. } = spec else {
        return Err("not prime explicit".into());
    };
    if trace["schema"] != SCHEMA
        || cache.len() < 76
        || &cache[..8] != MAGIC
        || trace["cache"]["sha256"] != sha256_hex(cache)
        || hex::encode(&cache[8..40]) != text(&trace["factor_base"], "fb_sha256")?
        || hex::encode(&cache[40..72]) != text(&trace["order_merkle"], "root")?
    {
        return Err("cache identity/integrity mismatch".into());
    }
    if order_root(
        text(&trace["curve"], "slug")?,
        text(&trace["factor_base"], "fb_sha256")?,
        trace["ordered_representatives"]
            .as_array()
            .ok_or("missing points")?,
    )? != text(&trace["order_merkle"], "root")?
    {
        return Err("ordered point commitment mismatch".into());
    }
    let count = u32::from_be_bytes(cache[72..76].try_into().unwrap()) as usize;
    if count
        != 2 * trace["ordered_representatives"]
            .as_array()
            .ok_or("missing points")?
            .len()
        || cache.len() != 76 + count * 13
    {
        return Err("cache length mismatch".into());
    }
    if !(5..(1 << 32)).contains(&p) || a >= p || b >= p {
        return Err("invalid word-sized curve".into());
    }
    Ok(Lookup {
        curve: PrimeCurve { p, a, b },
        cache,
        count,
    })
}
impl Lookup<'_> {
    fn query(&self, x: u64, y: u64) -> Result<Option<(u32, i8)>> {
        let Self {
            curve,
            cache,
            count,
        } = self;
        let point = PrimePoint::affine(x, y);
        if x >= curve.p || y >= curve.p || !curve.is_on_curve(point) {
            return Err("query is not a curve point".into());
        }
        let key = curve.key(&point);
        let (mut lo, mut hi) = (0, *count);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let row = &cache[76 + mid * 13..76 + (mid + 1) * 13];
            let found = u64::from_be_bytes(row[..8].try_into().unwrap());
            match found.cmp(&key) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    return Ok(Some((
                        u32::from_be_bytes(row[8..12].try_into().unwrap()),
                        row[12] as i8,
                    )))
                }
            }
        }
        Ok(None)
    }
}
#[cfg(test)]
fn lookup(trace: &Value, cache: &[u8], x: u64, y: u64) -> Result<Option<(u32, i8)>> {
    load_cache(trace, cache)?.query(x, y)
}
fn run() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("registry-records") if args.len() == 4 => registry_records(Path::new(&args[2]), Path::new(&args[3])),
        Some("export") if args.len() == 4 => export(Path::new(&args[2]), Path::new(&args[3])),
        Some("verify") if args.len() == 4 => verify(Path::new(&args[2]), Path::new(&args[3])),
        Some("lookup") if args.len() >= 6 && args.len().is_multiple_of(2) => {
            let trace = read_json(Path::new(&args[2]))?;
            let cache = fs::read(&args[3]).map_err(|e| e.to_string())?;
            let loaded = load_cache(&trace, &cache)?;
            for xy in args[4..].chunks_exact(2) {
                let x = xy[0].parse().map_err(|_| "invalid x")?;
                let y = xy[1].parse().map_err(|_| "invalid y")?;
                println!("{}", serde_json::to_string(&loaded.query(x, y)?).unwrap());
            }
            Ok(())
        }
        _ => Err("usage: factor_base_trace registry-records FROZEN_GENERIC_JSON NEW_JSON; export|verify FROZEN_GENERIC_JSON DIRECTORY; lookup TRACE CACHE X Y [X Y ...]".into()),
    }
}
fn main() {
    if let Err(e) = run() {
        eprintln!("factor_base_trace: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> Value {
        json!({"p":40129,"b":12823,"r":40127,"m":3,"F":6,"bits":16,
            "field":0,"curve":0,"coset_seed":0,"kind":"low_digit_weight",
            "base_extra":{"radix":4,"digit_candidates_drawn":41},
            "base":[170153719376651u64,171802986808769u64,3086u64]})
    }
    #[test]
    fn fb1_matches_the_existing_native_dump_identity() {
        use crypto_lib::cryptanalysis::ecbench::methods::dump_factor_base;
        let spec = CurveSpec::PrimeExplicit {
            p: 40129,
            a: 40126,
            b: 12823,
            group_order: 40127,
            r: 40127,
            gx: 0,
            gy: 3086,
        };
        let dump = dump_factor_base(&spec, "prime-abscissa:size=3").unwrap();
        let curve = PrimeCurve {
            p: 40129,
            a: 40126,
            b: 12823,
        };
        let mut keys: Vec<_> = dump
            .points
            .iter()
            .map(|pt| {
                curve.key(&PrimePoint::affine(
                    u64::from_str_radix(pt.x.trim_start_matches("0x"), 16).unwrap(),
                    u64::from_str_radix(pt.y.trim_start_matches("0x"), 16).unwrap(),
                ))
            })
            .collect();
        let params = dump
            .factor_base
            .params
            .iter()
            .map(|(k, v)| (k.as_str(), v.clone()))
            .collect();
        let (id, sha, points_sha) = identify(
            &dump.curve.slug,
            &dump.factor_base.family,
            &params,
            dump.factor_base.columns as usize,
            &mut keys,
        )
        .unwrap();
        assert_eq!(id, dump.factor_base.fb_id);
        assert_eq!(sha, dump.factor_base.fb_sha256);
        assert_eq!(points_sha, dump.factor_base.points_sha256);
    }
    #[test]
    fn exact_order_signs_and_corruption() {
        let (trace, mut cache) = make_trace(&fixture(), 0).unwrap();
        for (col, point) in trace["ordered_representatives"]
            .as_array()
            .unwrap()
            .iter()
            .enumerate()
        {
            let (x, y) = (point["x"].as_u64().unwrap(), point["y"].as_u64().unwrap());
            assert_eq!(lookup(&trace, &cache, x, y).unwrap(), Some((col as u32, 1)));
            assert_eq!(
                lookup(&trace, &cache, x, 40129 - y).unwrap(),
                Some((col as u32, -1))
            );
        }
        assert!(lookup(&trace, &cache, 0, 0).is_err());
        let outside = 35184372093309u64; // Another archived point, outside this 3-point fixture.
        assert_eq!(
            lookup(&trace, &cache, outside >> 32, outside & 0xffff_ffff).unwrap(),
            None
        );
        let mut reordered = trace.clone();
        reordered["ordered_representatives"]
            .as_array_mut()
            .unwrap()
            .swap(0, 1);
        assert!(load_cache(&reordered, &cache).is_err());
        assert!(load_cache(&trace, &cache[..75]).is_err());
        cache[76] ^= 1;
        assert!(lookup(&trace, &cache, 0, 3086).is_err());
    }
    #[test]
    fn order_commitment_is_separate_from_set_identity() {
        let a = fixture();
        let mut b = a.clone();
        b["base"].as_array_mut().unwrap().swap(0, 1);
        let (a, _) = make_trace(&a, 0).unwrap();
        let (b, _) = make_trace(&b, 0).unwrap();
        assert_eq!(a["factor_base"], b["factor_base"]);
        assert_ne!(a["order_merkle"], b["order_merkle"]);
    }
    #[test]
    fn rejects_duplicate_or_invalid_points() {
        let mut a = fixture();
        a["base"][1] = a["base"][0].clone();
        assert!(make_trace(&a, 0).is_err());
        a["base"][1] = json!(1);
        assert!(make_trace(&a, 0).is_err());
    }
}
