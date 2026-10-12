//! Independent source-curve replay for the preregistered folded-table gate.
//! It uses the general binary-curve group law and the least Frobenius
//! conjugate of x, not FrobeniusPairTable's normal-basis canonicaliser.

use crypto_lib::cryptanalysis::ic_boundary::{koblitz_instance, BinaryInstance};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use flate2::read::GzDecoder;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashMap};
use std::fs::{self, OpenOptions};
use std::io::{Read, Write};

const SUPPORT: &str = "research/notes/ecc2k130/n37_four_policy_support_20261003/RESULT.json.gz";
const SUPPORT_SHA: &str = "8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5";
const REFERENCE: &str =
    "research/notes/ecc2k130/n37_four_policy_pdp_20261003/RESULT_CI_FIXED.json.gz";
const REFERENCE_SHA: &str = "c84066a521137cfa61a2f490ed18a3a957b79da6e27f587ad1579b350cb9764a";
const B03: &str =
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b03.points.jsonl";
const B03_SHA: &str = "84400a2914f06e4a001d0f113f0195952f692634d4ff8285e23ade2599e2bde2";
const B04: &str =
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b04.points.jsonl";
const B04_SHA: &str = "72c5b5361a0b846c9626269ce3928820cb7d420aae479d31dc26f5cbf7e78aa6";
const SLUG: &str = "icv1-f2m37-tm534059-32aad96b";
const N: usize = 37;
const K: usize = 42;
const P: usize = 2 * N * K;
const R: u64 = 230_603_167;
const EXPECTED_TABLE_ADDS: u64 = 66_822;

struct Base {
    points: Vec<FastPoint>,
    cols: Vec<usize>,
    coefficients: Vec<u64>,
    index: HashMap<u64, usize>,
}

struct Table {
    by_orbit: HashMap<u64, (usize, usize)>,
    additions: u64,
    canonicalisations: u64,
}

#[derive(Default)]
struct QueryCount {
    additions: u64,
    lookups: u64,
    canonicalisations: u64,
}

fn pinned(path: &str, digest: &str) -> Result<Vec<u8>, String> {
    let bytes = fs::read(path).map_err(|e| format!("read {path}: {e}"))?;
    let actual = hex::encode(sha256(&bytes));
    if actual != digest {
        return Err(format!("{path}: pinned SHA-256 mismatch {actual}"));
    }
    Ok(bytes)
}

fn gzip_json(path: &str, digest: &str) -> Result<Value, String> {
    let bytes = pinned(path, digest)?;
    let mut raw = Vec::new();
    GzDecoder::new(&bytes[..])
        .read_to_end(&mut raw)
        .map_err(|e| format!("inflate {path}: {e}"))?;
    serde_json::from_slice(&raw).map_err(|e| format!("parse {path}: {e}"))
}

fn point(value: &Value) -> Result<FastPoint, String> {
    let xy: [u64; 2] = serde_json::from_value(value.clone()).map_err(|e| format!("point: {e}"))?;
    Ok(FastPoint::affine(xy[0], xy[1]))
}

fn general_add(inst: &BinaryInstance, p: FastPoint, q: FastPoint) -> FastPoint {
    let kc = inst.koblitz.as_ref().expect("Koblitz curve");
    inst.fast
        .lift(&kc.add(&inst.fast.lower(p), &inst.fast.lower(q)))
}

fn min_conjugate(inst: &BinaryInstance, x: u64) -> u64 {
    let mut v = x;
    let mut least = x;
    for _ in 1..N {
        v = inst.gf.sqr(v);
        least = least.min(v);
    }
    least
}

fn parse_base(inst: &BinaryInstance, policy: &Value, name: &str) -> Result<Base, String> {
    if policy["name"] != name
        || policy["curve"] != "source"
        || policy["columns"].as_u64() != Some(K as u64)
        || policy["signed_classes"].as_u64() != Some((K * N) as u64)
        || policy["physical_points"].as_u64() != Some(P as u64)
    {
        return Err(format!("{name}: bad support metadata"));
    }
    let entries = policy["entries"]
        .as_array()
        .ok_or("missing support entries")?;
    if entries.len() != P {
        return Err(format!("{name}: support length changed"));
    }
    let lambda = inst
        .koblitz
        .as_ref()
        .ok_or("no Koblitz data")?
        .lambda
        .to_u64()
        .ok_or("lambda overflow")?;
    let mut powers = [1u64; N];
    for i in 1..N {
        powers[i] = (u128::from(powers[i - 1]) * u128::from(lambda) % u128::from(R)) as u64;
    }
    let mut base = Base {
        points: Vec::with_capacity(P),
        cols: Vec::with_capacity(P),
        coefficients: Vec::with_capacity(P),
        index: HashMap::with_capacity(P),
    };
    for (i, e) in entries.iter().enumerate() {
        let col = i / (2 * N);
        let exp = (i / 2) % N;
        let sign = if i % 2 == 0 { 1 } else { -1 };
        let coef = if sign == 1 {
            powers[exp]
        } else {
            R - powers[exp]
        };
        let p = point(&e["point"])?;
        if e["column"].as_u64() != Some(col as u64)
            || e["exponent"].as_u64() != Some(exp as u64)
            || e["sign"].as_i64() != Some(sign)
            || e["coefficient"].as_u64() != Some(coef)
            || p.infinity
            || !inst.fast.is_on_curve(p)
            || !inst.fast.mul_u64(p, R).infinity
            || base.index.insert(p.pack(), i).is_some()
        {
            return Err(format!("{name}: invalid factor {i}"));
        }
        base.points.push(p);
        base.cols.push(col);
        base.coefficients.push(coef);
    }
    for i in 0..P {
        let c = i / (2 * N);
        let e = ((i / 2) % N + 1) % N;
        let next = 2 * N * c + 2 * e + i % 2;
        if base.points[i ^ 1] != inst.fast.neg(base.points[i])
            || base.points[next] != inst.fast.frobenius_k(base.points[i], 1)
        {
            return Err(format!("{name}: orbit label mismatch {i}"));
        }
    }
    Ok(base)
}

fn build_table(inst: &BinaryInstance, base: &Base) -> Table {
    let mut table = Table {
        by_orbit: HashMap::with_capacity(EXPECTED_TABLE_ADDS as usize),
        additions: 0,
        canonicalisations: 0,
    };
    for col in 0..K {
        let rep = col * 2 * N;
        for j in col * 2 * N..P {
            table.additions += 1;
            let sum = general_add(inst, base.points[rep], base.points[j]);
            if sum.infinity {
                continue;
            }
            table.canonicalisations += 1;
            table
                .by_orbit
                .entry(min_conjugate(inst, sum.x))
                .or_insert((rep, j));
        }
    }
    table
}

fn probe(
    inst: &BinaryInstance,
    base: &Base,
    table: &Table,
    target: FastPoint,
    count: &mut QueryCount,
) -> Result<Option<Vec<usize>>, String> {
    count.canonicalisations += 1;
    count.lookups += 1;
    let Some(&(i, j)) = table.by_orbit.get(&min_conjugate(inst, target.x)) else {
        return Ok(None);
    };
    count.additions += 1; // table's pair-image check
    let sum = general_add(inst, base.points[i], base.points[j]);
    let mut shifted = sum;
    for shift in 0..N {
        let sign = if shifted == target {
            1
        } else if inst.fast.neg(shifted) == target {
            -1
        } else {
            0
        };
        if sign != 0 {
            // Reconstruct images by applying the same shift to both factors.
            let mut pi = inst.fast.frobenius_k(base.points[i], shift as u32);
            let mut pj = inst.fast.frobenius_k(base.points[j], shift as u32);
            if sign == -1 {
                pi = inst.fast.neg(pi);
                pj = inst.fast.neg(pj);
            }
            let ai = *base.index.get(&pi.pack()).ok_or("missing folded image a")?;
            let bj = *base.index.get(&pj.pack()).ok_or("missing folded image b")?;
            if general_add(inst, base.points[ai], base.points[bj]) != target {
                return Err("independent pair image is wrong".into());
            }
            return Ok(Some(vec![ai, bj]));
        }
        shifted = inst.fast.frobenius_k(shifted, 1);
    }
    Err("independent orbit key matched but no pair image did".into())
}

fn query(
    inst: &BinaryInstance,
    base: &Base,
    table: &Table,
    target: FastPoint,
    m: u32,
) -> Result<(Option<Vec<usize>>, QueryCount), String> {
    let mut count = QueryCount::default();
    if target.infinity {
        return Ok((Some(Vec::new()), count));
    }
    count.lookups += 1; // direct factor-base lookup
    if let Some(&i) = base.index.get(&target.pack()) {
        return Ok((Some(vec![i]), count));
    }
    if let Some(pair) = probe(inst, base, table, target, &mut count)? {
        return Ok((Some(pair), count));
    }
    if m == 2 {
        return Ok((None, count));
    }
    for i in 0..P {
        count.additions += 1; // residual subtraction
        let residual = general_add(inst, target, inst.fast.neg(base.points[i]));
        if residual.infinity {
            continue;
        }
        if let Some(mut pair) = probe(inst, base, table, residual, &mut count)? {
            pair.insert(0, i);
            return Ok((Some(pair), count));
        }
    }
    Ok((None, count))
}

fn targets(inst: &BinaryInstance) -> Result<Vec<FastPoint>, String> {
    let mut out = Vec::new();
    for (path, digest) in [(B03, B03_SHA), (B04, B04_SHA)] {
        let bytes = pinned(path, digest)?;
        let mut count = 0usize;
        for line in bytes.split(|&b| b == b'\n').filter(|line| !line.is_empty()) {
            let q = point(&serde_json::from_slice::<Value>(line).map_err(|e| e.to_string())?)?;
            if q.infinity || !inst.fast.is_on_curve(q) || !inst.fast.mul_u64(q, R).infinity {
                return Err(format!("{path}: invalid target {count}"));
            }
            out.push(q);
            count += 1;
        }
        if count != 1024 {
            return Err(format!("{path}: expected 1024 targets"));
        }
    }
    Ok(out)
}

fn row(base: &Base, indices: &[usize]) -> Value {
    let mut values = BTreeMap::<usize, u64>::new();
    for &i in indices {
        let c = base.cols[i];
        values.insert(
            c,
            (values.get(&c).copied().unwrap_or(0) + base.coefficients[i]) % R,
        );
    }
    json!(values
        .into_iter()
        .filter(|(_, v)| *v != 0)
        .collect::<Vec<_>>())
}

fn check_record(
    inst: &BinaryInstance,
    base: &Base,
    table: &Table,
    target: FastPoint,
    recorded: &Value,
    historical: &Value,
    m: u32,
    arm: &str,
    target_index: usize,
) -> Result<(bool, QueryCount), String> {
    let label = format!("{arm} m{m} target {target_index}");
    let (independent, count) = query(inst, base, table, target, m)?;
    let hit = independent.is_some();
    let status = if hit { "hit" } else { "proved_miss" };
    if recorded["status"] != status
        || historical["status"] != status
        || recorded["query_adds"].as_u64() != Some(count.additions)
        || recorded["query_doubles"].as_u64() != Some(0)
        || recorded["query_scalar_mults"].as_u64() != Some(0)
        || recorded["lookups_uncharged"].as_u64() != Some(count.lookups)
        || recorded["canonicalisations_uncharged"].as_u64() != Some(count.canonicalisations)
        || recorded["frobfold_mismatches"].as_u64() != Some(0)
    {
        return Err(format!("{label}: decision or query counter mismatch"));
    }
    match &independent {
        None => {
            if !recorded["indices"].is_null() || !recorded["row"].is_null() {
                return Err(format!("{label}: miss carries a witness"));
            }
        }
        Some(_) => {
            let indices: Vec<usize> = serde_json::from_value(recorded["indices"].clone())
                .map_err(|e| format!("{label}: witness parse {e}"))?;
            if indices.len() > m as usize || indices.iter().any(|&i| i >= P) {
                return Err(format!("{label}: witness size or index invalid"));
            }
            let mut sum = FastPoint::INFINITY;
            for &i in &indices {
                sum = general_add(inst, sum, base.points[i]);
            }
            if sum != target || recorded["row"] != row(base, &indices) {
                return Err(format!("{label}: witness sum or relation row invalid"));
            }
        }
    }
    Ok((hit, count))
}

fn replay(result_path: &str, support_path: &str) -> Result<Value, String> {
    let raw = fs::read(result_path).map_err(|e| format!("read result: {e}"))?;
    let result: Value = serde_json::from_slice(&raw).map_err(|e| format!("result JSON: {e}"))?;
    if result["schema"] != "n37-frobenius-fold-gate/v1"
        || result["status"] != "PASS"
        || result["curve_slug"] != SLUG
        || result["support_gzip_sha256"] != SUPPORT_SHA
        || result["complete_reference_gzip_sha256"] != REFERENCE_SHA
        || result["point_only_sha256"]["b03"] != B03_SHA
        || result["point_only_sha256"]["b04"] != B04_SHA
        || result["physical_points_each"].as_u64() != Some(P as u64)
        || result["columns_each"].as_u64() != Some(K as u64)
        || result["target_count"].as_u64() != Some(2048)
        || result["complete_table_pair_adds"].as_u64() != Some((P * (P + 1) / 2) as u64)
        || result["expected_folded_pair_adds"].as_u64() != Some(EXPECTED_TABLE_ADDS)
        || !result["s"].is_null()
        || !result["rho_ratio"].is_null()
        || !result["speedup"].is_null()
    {
        return Err("result identity or cost scope changed".into());
    }
    let inst = koblitz_instance(0, 37).ok_or("construct source curve")?;
    if inst.curve_id().slug != SLUG || inst.r != R {
        return Err("source curve identity changed".into());
    }
    let support = gzip_json(support_path, SUPPORT_SHA)?;
    if support["schema"] != "n37-four-policy-support-v1"
        || support["status"] != "PASS"
        || support["r"].as_u64() != Some(R)
    {
        return Err("support identity changed".into());
    }
    let reference = gzip_json(REFERENCE, REFERENCE_SHA)?;
    if reference["schema"] != "n37-four-policy-pdp-v1"
        || reference["status"] != "PASS"
        || reference["target_count"].as_u64() != Some(2048)
    {
        return Err("reference identity changed".into());
    }
    let support_policies = support["policies"].as_array().ok_or("support arms")?;
    let old_policies = reference["policies"].as_array().ok_or("reference arms")?;
    let new_policies = result["policies"].as_array().ok_or("result arms")?;
    if support_policies.len() != 4 || old_policies.len() != 4 || new_policies.len() != 2 {
        return Err("arm count changed".into());
    }
    let targets = targets(&inst)?;
    let mut arm_receipts = Vec::new();
    for (out_index, (in_index, name)) in [(0usize, "original"), (3usize, "pullback")]
        .into_iter()
        .enumerate()
    {
        let arm = &new_policies[out_index];
        if arm["name"] != name || old_policies[in_index]["name"] != name {
            return Err(format!("{name}: arm identity mismatch"));
        }
        let base = parse_base(&inst, &support_policies[in_index], name)?;
        let table = build_table(&inst, &base);
        let reported_table = &arm["table"];
        let minimum_retained_bytes = (P * N * std::mem::size_of::<u32>()) as u64
            + table.by_orbit.len() as u64 * std::mem::size_of::<(u64, (u32, u32))>() as u64;
        if table.additions != EXPECTED_TABLE_ADDS
            || reported_table["representatives"].as_u64() != Some(K as u64)
            || reported_table["entries"].as_u64() != Some(table.by_orbit.len() as u64)
            || reported_table["build_adds"].as_u64() != Some(table.additions)
            || reported_table["build_doubles"].as_u64() != Some(0)
            || reported_table["build_scalar_mults"].as_u64() != Some(0)
            || reported_table["canonicalisations_uncharged"].as_u64()
                != Some(table.canonicalisations)
            || reported_table["frobenius_maps_uncharged"].as_u64() != Some((P * (N - 1)) as u64)
            || reported_table["lookups_uncharged"].as_u64() != Some((P * (N - 1)) as u64)
            || reported_table["minimum_retained_bytes"].as_u64() != Some(minimum_retained_bytes)
        {
            return Err(format!("{name}: independent table count mismatch"));
        }
        let records = arm["targets"].as_array().ok_or("result target records")?;
        let historical = old_policies[in_index]["targets"]
            .as_array()
            .ok_or("reference target records")?;
        if records.len() != targets.len() || historical.len() != targets.len() {
            return Err(format!("{name}: target count mismatch"));
        }
        let mut hits = [0u64; 2];
        let mut adds = [0u64; 2];
        let mut lookups = [0u64; 2];
        let mut canons = [0u64; 2];
        for (i, &target) in targets.iter().enumerate() {
            if records[i]["index"].as_u64() != Some(i as u64)
                || historical[i]["index"].as_u64() != Some(i as u64)
            {
                return Err(format!("{name}: target order mismatch at {i}"));
            }
            for (j, m) in [2u32, 3].into_iter().enumerate() {
                let field = format!("m{m}");
                let (hit, count) = check_record(
                    &inst,
                    &base,
                    &table,
                    target,
                    &records[i][&field],
                    &historical[i][&field],
                    m,
                    name,
                    i,
                )?;
                hits[j] += u64::from(hit);
                adds[j] += count.additions;
                lookups[j] += count.lookups;
                canons[j] += count.canonicalisations;
            }
        }
        let summary = &arm["summary"];
        if summary["m2_hits"].as_u64() != Some(hits[0])
            || summary["m3_hits"].as_u64() != Some(hits[1])
            || summary["m2_query_adds"].as_u64() != Some(adds[0])
            || summary["m3_query_adds"].as_u64() != Some(adds[1])
            || summary["m2_lookups_uncharged"].as_u64() != Some(lookups[0])
            || summary["m3_lookups_uncharged"].as_u64() != Some(lookups[1])
            || summary["m2_canonicalisations_uncharged"].as_u64() != Some(canons[0])
            || summary["m3_canonicalisations_uncharged"].as_u64() != Some(canons[1])
        {
            return Err(format!("{name}: query totals mismatch"));
        }
        arm_receipts.push(json!({
            "name":name,
            "table_entries":table.by_orbit.len(),
            "table_additions":table.additions,
            "m2_hits":hits[0],
            "m3_hits":hits[1],
            "m2_query_adds":adds[0],
            "m3_query_adds":adds[1],
            "decisions_checked":2 * targets.len(),
        }));
    }
    Ok(json!({
        "schema":"n37-frobenius-fold-replay/v1",
        "status":"PASS",
        "result_sha256":hex::encode(sha256(&raw)),
        "support_gzip_sha256":SUPPORT_SHA,
        "complete_reference_gzip_sha256":REFERENCE_SHA,
        "decisions_checked":4 * targets.len(),
        "arms":arm_receipts,
    }))
}

fn main() {
    let mut args = std::env::args().skip(1);
    let result = args
        .next()
        .expect("usage: n37_frobenius_fold_replay RESULT.json OUTPUT.json [SUPPORT.json.gz]");
    let output = args.next().expect("missing OUTPUT.json");
    let support = args.next().unwrap_or_else(|| SUPPORT.to_string());
    assert!(args.next().is_none(), "unexpected argument");
    let receipt = replay(&result, &support).expect("independent folded replay failed");
    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&output)
        .expect("refuse to overwrite replay receipt");
    serde_json::to_writer_pretty(&mut file, &receipt).expect("write replay receipt");
    file.write_all(b"\n").expect("newline");
    println!("{output}");
}
