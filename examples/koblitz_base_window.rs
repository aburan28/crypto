//! Deterministic accepted-orbit windows for point-defined Koblitz bases.
//!
//! This is a separate input generator, not a change to the measured compact
//! DLP executable. It has no target-point input and never builds an index.
//! Usage:
//!   koblitz_base_window <n> <a> <columns> <window> <raw-x-cap> \
//!       <header.jsonl> <scan.jsonl> <receipt.json>

use crypto_lib::cryptanalysis::koblitz_fast_arith::{FastBinaryCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::json;
use std::collections::{BTreeMap, HashSet};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

struct Spec {
    n: u32,
    a: u8,
    columns: usize,
    window: usize,
    raw_x_cap: u64,
}

fn mulmod(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 * b as u128) % r as u128) as u64
}

fn write_new(path: &str, bytes: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)?;
    file.write_all(bytes)?;
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 9 {
        return Err("usage: koblitz_base_window <n> <a> <columns> <window> <raw-x-cap> <header.jsonl> <scan.jsonl> <receipt.json>".into());
    }
    let spec = Spec {
        n: args[1].parse()?,
        a: args[2].parse()?,
        columns: args[3].parse()?,
        window: args[4].parse()?,
        raw_x_cap: args[5].parse()?,
    };
    if spec.columns == 0 || spec.raw_x_cap == 0 {
        return Err("columns and raw-x-cap must be positive".into());
    }
    let start = spec
        .columns
        .checked_mul(spec.window)
        .ok_or("window start overflows usize")?;
    let end = start
        .checked_add(spec.columns)
        .ok_or("window end overflows usize")?;
    if args[6] == args[7] || args[6] == args[8] || args[7] == args[8] {
        return Err("header, scan and receipt paths must be distinct".into());
    }
    if args[6..=8].iter().any(|path| Path::new(path).exists()) {
        return Err("output paths must not already exist".into());
    }
    let started = Instant::now();
    let curve = KoblitzCurve::new(spec.a, spec.n).ok_or("unsupported Koblitz rung")?;
    let r = curve
        .subgroup_order
        .to_u64()
        .ok_or("subgroup order does not fit u64")?;
    let lambda = curve.lambda.to_u64().ok_or("lambda does not fit u64")?;
    let fast =
        FastBinaryCurve::new(&curve.curve.irreducible, spec.a as u64).ok_or("field too large")?;
    if spec.raw_x_cap > ((1u64 << spec.n) - 1) {
        return Err("raw-x cap exceeds nonzero field coordinates".into());
    }
    let gf = &fast.gf;
    let b = gf.from_element(&curve.curve.b);
    let n = spec.n as usize;

    let mut trace = BufWriter::new(
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&args[7])?,
    );
    let mut trace_hash = blake3::Hasher::new();
    let mut decision_count = 0u64;
    let mut status_counts = BTreeMap::<&'static str, u64>::new();
    let mut emit = |value: serde_json::Value,
                    status: &'static str|
     -> Result<(), Box<dyn std::error::Error>> {
        let mut line = serde_json::to_vec(&value)?;
        line.push(b'\n');
        trace.write_all(&line)?;
        trace_hash.update(&line);
        decision_count += 1;
        *status_counts.entry(status).or_default() += 1;
        Ok(())
    };

    let mut seen = HashSet::<u64>::new();
    let mut accepted_keys = Vec::<u64>::with_capacity(end);
    let mut selected_keys = Vec::<u64>::with_capacity(spec.columns);
    let mut points = Vec::<FastPoint>::with_capacity(spec.columns * 2 * n);
    let mut labels = Vec::<(usize, u64)>::with_capacity(spec.columns * 2 * n);
    let mut reps = Vec::<Option<[u64; 2]>>::with_capacity(spec.columns);
    let mut scanned = 0u64;
    'scan: for raw_x in 1..=spec.raw_x_cap {
        scanned += 1;
        let lifted = fast.points_with_x(b, raw_x);
        if lifted.is_empty() {
            emit(json!({"raw_x":raw_x,"status":"no_lift"}), "no_lift")?;
        }
        for (lift_index, point) in lifted.into_iter().enumerate() {
            let lifted_xy = point.map(|(x, y)| [x, y]);
            let Some((px, py)) = fast.scalar_mul(point, &curve.cofactor) else {
                emit(
                    json!({"raw_x":raw_x,"lift_index":lift_index,"lifted":lifted_xy,"status":"projected_infinity"}),
                    "projected_infinity",
                )?;
                continue;
            };
            if px <= 1 {
                emit(
                    json!({"raw_x":raw_x,"lift_index":lift_index,"lifted":lifted_xy,"projected":[px,py],"status":"small_projected_x"}),
                    "small_projected_x",
                )?;
                continue;
            }
            let mut key = px;
            let mut x = px;
            for _ in 1..n {
                x = gf.sqr(x);
                key = key.min(x);
            }
            if !seen.insert(key) {
                emit(
                    json!({"raw_x":raw_x,"lift_index":lift_index,"lifted":lifted_xy,"projected":[px,py],"orbit_key":key,"status":"duplicate_orbit"}),
                    "duplicate_orbit",
                )?;
                continue;
            }
            let ordinal = accepted_keys.len();
            accepted_keys.push(key);
            let selected = (start..end).contains(&ordinal);
            emit(
                json!({"raw_x":raw_x,"lift_index":lift_index,"lifted":lifted_xy,"projected":[px,py],"orbit_key":key,"accepted_ordinal":ordinal,"selected":selected,"status":"accepted"}),
                "accepted",
            )?;
            if selected {
                let column = reps.len();
                reps.push(Some([px, py]));
                selected_keys.push(key);
                let (mut cx, mut cy) = (px, py);
                let mut coefficient = 1u64;
                for _ in 0..n {
                    points.push(Some((cx, cy)));
                    labels.push((column, coefficient));
                    points.push(Some((cx, cx ^ cy)));
                    labels.push((column, (r - coefficient) % r));
                    cx = gf.sqr(cx);
                    cy = gf.sqr(cy);
                    coefficient = mulmod(coefficient, lambda, r);
                }
            }
            if accepted_keys.len() == end {
                break 'scan;
            }
        }
    }
    trace.flush()?;
    let trace_blake3 = trace_hash.finalize().to_hex().to_string();
    let success = accepted_keys.len() == end;
    let receipt = json!({
        "schema":"koblitz-base-window-generator-receipt-v1",
        "status":if success { "complete" } else { "raw_x_cap" },
        "n":spec.n,
        "a":spec.a,
        "subgroup_order":r,
        "cofactor":curve.cofactor.to_u64().ok_or("cofactor does not fit u64")?,
        "field_modulus_low_terms":curve.curve.irreducible.low_terms,
        "columns":spec.columns,
        "window":spec.window,
        "accepted_start":start,
        "accepted_end":end,
        "raw_x_cap":spec.raw_x_cap,
        "raw_x_scanned":scanned,
        "accepted_orbits":accepted_keys.len(),
        "selected_orbits":reps.len(),
        "accepted_orbit_keys":accepted_keys,
        "selected_orbit_keys":selected_keys,
        "trace_decisions":decision_count,
        "trace_status_counts":status_counts,
        "trace_blake3":trace_blake3,
        "elapsed_wall_seconds":started.elapsed().as_secs_f64(),
    });
    write_new(&args[8], format!("{receipt}\n").as_bytes())?;
    if !success {
        return Err(format!(
            "raw-x cap {} reached after {} accepted orbits; need {}",
            spec.raw_x_cap, receipt["accepted_orbits"], end
        )
        .into());
    }
    if reps.len() != spec.columns || points.len() != 2 * n * spec.columns {
        return Err("selected base has wrong cardinality".into());
    }
    let rep = reps[0].map(|[x, y]| (x, y));
    assert_eq!(
        fast.scalar_mul(rep, &BigUint::from(lambda)),
        rep.map(|(x, y)| (gf.sqr(x), gf.sqr(y))),
        "Frobenius must act as lambda on the subgroup"
    );
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"compact-orbit-constructed-base-v1");
    for (point, &(column, coefficient)) in points.iter().zip(&labels) {
        let (x, y) = point.ok_or("orbit member at infinity")?;
        for word in [x, y, column as u64, coefficient] {
            hasher.update(&word.to_le_bytes());
        }
    }
    let base_hash = hasher.finalize().to_hex().to_string();
    let header = json!({
        "kind":"point_defined_factor_base",
        "n":spec.n,
        "a":spec.a,
        "subgroup_order":r,
        "field_modulus_low_terms":curve.curve.irreducible.low_terms,
        "orbit_columns":spec.columns,
        "base_hash":base_hash,
        "factor_base_points":points.len(),
        "factor_base_point_coordinates":points.iter().map(|p| p.map(|(x,y)| [x,y])).collect::<Vec<_>>(),
        "factor_base_point_labels":labels,
        "factor_base_representatives":reps,
    });
    write_new(&args[6], format!("{header}\n").as_bytes())?;
    Ok(())
}
