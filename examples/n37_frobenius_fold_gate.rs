//! Preregistered source-curve folded pair-table exactness and cost gate.
//! Protocol: research/notes/ecc2k130/n37_frobenius_fold_gate_20261003/PROTOCOL.md
//! Historical target files contain points only; this producer reads no scalar fixtures.

use crypto_lib::cryptanalysis::ic_boundary::{
    decompose_mitm_frobenius, koblitz_instance, BinaryGroup, CountedGroup, FactorBase,
    FrobeniusPairTable, GroupOps, OracleCounters,
};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use flate2::read::GzDecoder;
use num_traits::ToPrimitive;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs::{self, OpenOptions};
use std::io::{Read, Write};
use std::time::Instant;

const SUPPORT: &str = "research/notes/ecc2k130/n37_four_policy_support_20261003/RESULT.json.gz";
const SUPPORT_SHA: &str = "8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5";
const SUPPORT_RAW_SHA: &str = "05664103a6dab090a8b4298f34c0964e6f253f625cd73d661562d85a2fce87af";
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

fn pinned_bytes(path: &str, digest: &str) -> Result<Vec<u8>, String> {
    let bytes = fs::read(path).map_err(|e| format!("read {path}: {e}"))?;
    let got = hex::encode(sha256(&bytes));
    if got != digest {
        return Err(format!("{path}: SHA-256 {got} differs from {digest}"));
    }
    Ok(bytes)
}

fn pinned_gzip(path: &str, digest: &str) -> Result<Value, String> {
    let bytes = pinned_bytes(path, digest)?;
    let mut raw = Vec::new();
    GzDecoder::new(&bytes[..])
        .read_to_end(&mut raw)
        .map_err(|e| format!("inflate {path}: {e}"))?;
    if path == SUPPORT && hex::encode(sha256(&raw)) != SUPPORT_RAW_SHA {
        return Err("expanded support SHA-256 changed".into());
    }
    serde_json::from_slice(&raw).map_err(|e| format!("parse {path}: {e}"))
}

fn point(value: &Value) -> Result<FastPoint, String> {
    let xy: [u64; 2] =
        serde_json::from_value(value.clone()).map_err(|e| format!("point coordinates: {e}"))?;
    Ok(FastPoint::affine(xy[0], xy[1]))
}

fn targets(
    inst: &crypto_lib::cryptanalysis::ic_boundary::BinaryInstance,
) -> Result<Vec<FastPoint>, String> {
    let mut out = Vec::with_capacity(2048);
    for (path, digest) in [(B03, B03_SHA), (B04, B04_SHA)] {
        let bytes = pinned_bytes(path, digest)?;
        let lines = bytes.split(|&b| b == b'\n').filter(|line| !line.is_empty());
        let mut count = 0;
        for line in lines {
            let value: Value =
                serde_json::from_slice(line).map_err(|e| format!("target in {path}: {e}"))?;
            let q = point(&value)?;
            if q.infinity || !inst.fast.is_on_curve(q) || !inst.fast.mul_u64(q, R).infinity {
                return Err(format!("{path}: invalid target {count}"));
            }
            out.push(q);
            count += 1;
        }
        if count != 1024 {
            return Err(format!("{path}: expected 1024 targets, found {count}"));
        }
    }
    Ok(out)
}

fn base(
    inst: &crypto_lib::cryptanalysis::ic_boundary::BinaryInstance,
    policy: &Value,
    name: &str,
) -> Result<FactorBase<FastPoint>, String> {
    if policy["name"] != name
        || policy["curve"] != "source"
        || policy["columns"].as_u64() != Some(K as u64)
        || policy["signed_classes"].as_u64() != Some((K * N) as u64)
        || policy["physical_points"].as_u64() != Some(P as u64)
    {
        return Err(format!("{name}: frozen policy identity changed"));
    }
    let entries = policy["entries"].as_array().ok_or("no factor entries")?;
    if entries.len() != P {
        return Err(format!("{name}: factor count {} != {P}", entries.len()));
    }
    let kc = inst.koblitz.as_ref().ok_or("source is not Koblitz")?;
    let lambda = kc.lambda.to_u64().ok_or("lambda outside u64")?;
    let mut powers = [0u64; N];
    powers[0] = 1;
    for i in 1..N {
        powers[i] = (u128::from(powers[i - 1]) * u128::from(lambda) % u128::from(R)) as u64;
    }
    if (u128::from(powers[N - 1]) * u128::from(lambda) % u128::from(R)) != 1 {
        return Err("lambda does not close after 37 steps".into());
    }
    let mut points = Vec::with_capacity(P);
    let mut cols = Vec::with_capacity(P);
    let mut coefs = Vec::with_capacity(P);
    for (i, e) in entries.iter().enumerate() {
        let col = i / (2 * N);
        let exponent = (i / 2) % N;
        let sign = if i % 2 == 0 { 1 } else { -1 };
        let coefficient = if sign == 1 {
            powers[exponent]
        } else {
            R - powers[exponent]
        };
        let p = point(&e["point"])?;
        if e["column"].as_u64() != Some(col as u64)
            || e["exponent"].as_u64() != Some(exponent as u64)
            || e["sign"].as_i64() != Some(sign)
            || e["coefficient"].as_u64() != Some(coefficient)
            || p.infinity
            || p.x == 0
            || !inst.fast.is_on_curve(p)
            || !inst.fast.mul_u64(p, R).infinity
        {
            return Err(format!("{name}: invalid point or label {i}"));
        }
        points.push(p);
        cols.push(col);
        coefs.push(coefficient);
    }
    let g = BinaryGroup(&inst.fast);
    let fb = FactorBase::from_column_map(
        format!("frozen n37 {name} source support"),
        points,
        cols,
        coefs,
        K,
        |p| g.key(p),
        |p| g.key(&g.neg(*p)),
        |p| p.x,
    )?;
    for i in 0..P {
        let col = i / (2 * N);
        let next_exp = ((i / 2) % N + 1) % N;
        let expected = col * 2 * N + 2 * next_exp + i % 2;
        if inst.fast.frobenius_k(fb.points[i], 1) != fb.points[expected]
            || fb.neg_index[i] != (i ^ 1)
        {
            return Err(format!(
                "{name}: Frobenius/sign label closure failed at {i}"
            ));
        }
    }
    Ok(fb)
}

fn row(fb: &FactorBase<FastPoint>, indices: &[usize]) -> Value {
    let mut values = BTreeMap::<usize, u64>::new();
    for &i in indices {
        let col = fb.col_of[i];
        let old = values.get(&col).copied().unwrap_or(0);
        values.insert(col, (old + fb.coef_of[i]) % R);
    }
    json!(values
        .into_iter()
        .filter(|(_, v)| *v != 0)
        .collect::<Vec<_>>())
}

fn query(
    g: &BinaryGroup<'_>,
    fb: &FactorBase<FastPoint>,
    table: &FrobeniusPairTable,
    target: FastPoint,
    m: u32,
) -> Result<Value, String> {
    let mut ops = GroupOps::default();
    let mut ctr = OracleCounters::default();
    let witness = if target.infinity {
        Some(Vec::new())
    } else {
        ctr.lookups += 1; // direct factor-base hash lookup
        if let Some(i) = fb.index_of_key(g.key(&target)) {
            Some(vec![i])
        } else {
            let pair = decompose_mitm_frobenius(g, fb, table, 2, &mut ops, &mut ctr, target);
            if pair.is_some() || m == 2 {
                pair
            } else {
                decompose_mitm_frobenius(g, fb, table, 3, &mut ops, &mut ctr, target)
            }
        }
    };
    if ctr.frobfold_mismatches != 0 {
        return Err("folded pair-image mismatch".into());
    }
    if let Some(indices) = &witness {
        if indices.len() > m as usize || indices.iter().any(|&i| i >= fb.points.len()) {
            return Err("bad folded witness length or index".into());
        }
        let mut sum = g.identity();
        for &i in indices {
            sum = g.0.add(sum, fb.points[i]);
        }
        if sum != target {
            return Err("folded witness does not sum to target".into());
        }
    }
    Ok(json!({
        "status": if witness.is_some() { "hit" } else { "proved_miss" },
        "indices": witness,
        "row": witness.as_ref().map(|indices| row(fb, indices)),
        "query_adds": ops.adds,
        "query_doubles": ops.doubles,
        "query_scalar_mults": ops.scalar_mults,
        "lookups_uncharged": ctr.lookups,
        "canonicalisations_uncharged": ctr.canonicalisations,
        "frobfold_mismatches": ctr.frobfold_mismatches,
    }))
}

fn run() -> Result<Value, String> {
    let inst = koblitz_instance(0, 37).ok_or("construct registered source")?;
    if inst.curve_id().slug != SLUG || inst.r != R || inst.n != 37 || inst.a != 0 {
        return Err("source curve identity changed".into());
    }
    let frozen = pinned_gzip(SUPPORT, SUPPORT_SHA)?;
    if frozen["schema"] != "n37-four-policy-support-v1"
        || frozen["status"] != "PASS"
        || frozen["r"].as_u64() != Some(R)
    {
        return Err("support manifest header changed".into());
    }
    let policies = frozen["policies"]
        .as_array()
        .ok_or("support policies missing")?;
    if policies.len() != 4 {
        return Err("support policy count changed".into());
    }
    let reference = pinned_gzip(REFERENCE, REFERENCE_SHA)?;
    if reference["schema"] != "n37-four-policy-pdp-v1"
        || reference["status"] != "PASS"
        || reference["target_count"].as_u64() != Some(2048)
    {
        return Err("complete-table reference identity changed".into());
    }
    let old = reference["policies"]
        .as_array()
        .ok_or("reference policies missing")?;
    if old.len() != 4 {
        return Err("reference arm count changed".into());
    }
    let targets = targets(&inst)?;
    let g = BinaryGroup(&inst.fast);
    let mut arms = Vec::new();
    for (index, name) in [(0usize, "original"), (3usize, "pullback")] {
        let fb = base(&inst, &policies[index], name)?;
        if old[index]["name"] != name {
            return Err(format!("reference arm {index} is not {name}"));
        }
        let old_targets = old[index]["targets"]
            .as_array()
            .ok_or("reference targets missing")?;
        if old_targets.len() != targets.len() {
            return Err(format!("{name}: reference target count changed"));
        }
        let table = FrobeniusPairTable::build(&inst, &fb)
            .ok_or_else(|| format!("{name}: Frobenius table refused closed base"))?;
        if table.representatives != K as u64
            || table.build_ops.adds != EXPECTED_TABLE_ADDS
            || table.build_ops.doubles != 0
            || table.build_ops.scalar_mults != 0
            || table.build_frobenius_maps != (P * (N - 1)) as u64
            || table.build_lookups != table.build_frobenius_maps
        {
            return Err(format!(
                "{name}: table setup ledger differs from preregistration"
            ));
        }
        let mut records = Vec::with_capacity(targets.len());
        let mut hit2 = 0u64;
        let mut hit3 = 0u64;
        let mut adds2 = 0u64;
        let mut adds3 = 0u64;
        let mut lookups2 = 0u64;
        let mut lookups3 = 0u64;
        let mut canons2 = 0u64;
        let mut canons3 = 0u64;
        let start = Instant::now();
        for (i, &target) in targets.iter().enumerate() {
            let q2 = query(&g, &fb, &table, target, 2)?;
            let q3 = query(&g, &fb, &table, target, 3)?;
            for (m, q) in [("m2", &q2), ("m3", &q3)] {
                if q["status"] != old_targets[i][m]["status"] {
                    return Err(format!(
                        "{name}: {m} decision differs from exact reference at {i}"
                    ));
                }
            }
            hit2 += u64::from(q2["status"] == "hit");
            hit3 += u64::from(q3["status"] == "hit");
            adds2 += q2["query_adds"].as_u64().ok_or("m2 additions")?;
            adds3 += q3["query_adds"].as_u64().ok_or("m3 additions")?;
            lookups2 += q2["lookups_uncharged"].as_u64().ok_or("m2 lookups")?;
            lookups3 += q3["lookups_uncharged"].as_u64().ok_or("m3 lookups")?;
            canons2 += q2["canonicalisations_uncharged"]
                .as_u64()
                .ok_or("m2 canons")?;
            canons3 += q3["canonicalisations_uncharged"]
                .as_u64()
                .ok_or("m3 canons")?;
            records.push(json!({"index":i,"m2":q2,"m3":q3}));
        }
        let query_wall_ns = start.elapsed().as_nanos() as u64;
        let minimum_retained_bytes = (P * N * std::mem::size_of::<u32>()) as u64
            + table.entries * std::mem::size_of::<(u64, (u32, u32))>() as u64;
        arms.push(json!({
            "name":name,
            "table":{
                "representatives":table.representatives,
                "entries":table.entries,
                "build_adds":table.build_ops.adds,
                "build_doubles":table.build_ops.doubles,
                "build_scalar_mults":table.build_ops.scalar_mults,
                "canonicalisations_uncharged":table.build_canonicalisations,
                "frobenius_maps_uncharged":table.build_frobenius_maps,
                "lookups_uncharged":table.build_lookups,
                "minimum_retained_bytes":minimum_retained_bytes,
                "build_wall_ns_descriptive":table.build_wall_ns,
            },
            "summary":{
                "m2_hits":hit2,"m3_hits":hit3,
                "m2_query_adds":adds2,"m3_query_adds":adds3,
                "m2_lookups_uncharged":lookups2,"m3_lookups_uncharged":lookups3,
                "m2_canonicalisations_uncharged":canons2,
                "m3_canonicalisations_uncharged":canons3,
                "query_wall_ns_descriptive":query_wall_ns,
            },
            "targets":records,
        }));
    }
    Ok(json!({
        "schema":"n37-frobenius-fold-gate/v1",
        "status":"PASS",
        "curve_slug":SLUG,
        "support_gzip_sha256":SUPPORT_SHA,
        "complete_reference_gzip_sha256":REFERENCE_SHA,
        "point_only_sha256":{"b03":B03_SHA,"b04":B04_SHA},
        "physical_points_each":P,
        "columns_each":K,
        "target_count":targets.len(),
        "complete_table_pair_adds":P * (P + 1) / 2,
        "expected_folded_pair_adds":EXPECTED_TABLE_ADDS,
        "s":Value::Null,
        "rho_ratio":Value::Null,
        "speedup":Value::Null,
        "cost_scope":"frozen-support table and query group additions; native work and base selection unpriced",
        "policies":arms,
    }))
}

fn main() {
    let mut args = std::env::args().skip(1);
    let output = args
        .next()
        .expect("usage: n37_frobenius_fold_gate OUTPUT.json");
    assert!(args.next().is_none(), "unexpected argument");
    let result = run().expect("folded table gate failed");
    let mut file = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&output)
        .expect("refuse to overwrite output");
    serde_json::to_writer_pretty(&mut file, &result).expect("write result");
    file.write_all(b"\n").expect("newline");
    println!("{output}");
}
