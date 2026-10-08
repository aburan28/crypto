//! Produce the preregistered n37 four-policy equal-support census.
//! No target point, target scalar, PDP solver or rho measurement is read.

#![recursion_limit = "256"]

#[path = "support/n37_policy_common.rs"]
mod common;

use common::{point_from_value, point_words, Context, K, N, R};
use crypto_lib::cryptanalysis::binary_velu::{velu_point_map, velu_preimage_point, Curve};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs::OpenOptions;
use std::io::Write;
use std::time::Instant;

#[cfg(unix)]
fn peak_rss_before_output_bytes() -> Option<u64> {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
        return None;
    }
    let rss = unsafe { usage.assume_init() }.ru_maxrss;
    if rss < 0 {
        return None;
    }
    #[cfg(target_os = "macos")]
    {
        Some(rss as u64)
    }
    #[cfg(not(target_os = "macos"))]
    {
        Some((rss as u64).saturating_mul(1024))
    }
}

#[cfg(not(unix))]
fn peak_rss_before_output_bytes() -> Option<u64> {
    None
}

fn fail(message: impl Into<String>) -> Result<Value, String> {
    Err(message.into())
}

fn map_to_leaf(ctx: &Context, point: FastPoint) -> Result<FastPoint, String> {
    velu_point_map(
        ctx.bridge.map_point(point).ok_or("field bridge point")?,
        &ctx.isogeny.kernel,
        37,
        &ctx.archived_source.irr,
    )
    .ok_or("undefined complete degree-73 map".into())
}

fn validate_policy(
    name: &str,
    curve: &Curve,
    positives: &[FastPoint],
    powers: &[u64],
) -> Result<Value, String> {
    if positives.len() != K * N || powers.len() != N {
        return Err(format!("{name}: wrong point or coefficient count"));
    }
    let mut signed_classes = HashSet::new();
    let mut entries = Vec::with_capacity(2 * K * N);
    for (index, &positive) in positives.iter().enumerate() {
        let column = index / N;
        let exponent = index % N;
        let negative = curve.fast.neg(positive);
        for (sign, point, coefficient) in [
            (1i8, positive, powers[exponent]),
            (-1i8, negative, (R - powers[exponent]) % R),
        ] {
            if point.infinity
                || point.x == 0
                || !curve.fast.is_on_curve(point)
                || !curve.fast.mul_u64(point, R).infinity
            {
                return Err(format!(
                    "{name}: invalid subgroup point at {column}/{exponent}/{sign}"
                ));
            }
            entries.push(json!({
                "point": point_words(point)?,
                "column": column,
                "exponent": exponent,
                "sign": sign,
                "coefficient": coefficient,
            }));
        }
        if !signed_classes.insert(positive.x) {
            return Err(format!(
                "{name}: signed orbit collision at column {column}, exponent {exponent}"
            ));
        }
    }
    if signed_classes.len() != K * N {
        return Err(format!("{name}: useful support is not 1554 classes"));
    }
    Ok(json!({
        "name": name,
        "curve": if name == "original" || name == "pullback" { "source" } else { "leaf" },
        "columns": K,
        "signed_classes": signed_classes.len(),
        "physical_points": entries.len(),
        "entries": entries,
    }))
}

fn source_base_matches_archive(
    ctx: &Context,
    positives: &[FastPoint],
    powers: &[u64],
) -> Result<(), String> {
    let coords = ctx.source["factor_base_point_coordinates"]
        .as_array()
        .ok_or("source base point coordinates")?;
    let labels = ctx.source["factor_base_point_labels"]
        .as_array()
        .ok_or("source base point labels")?;
    if coords.len() != 2 * K * N || labels.len() != coords.len() {
        return Err("source base archive point/label count changed".into());
    }
    for (index, &positive) in positives.iter().enumerate() {
        let column = index / N;
        let exponent = index % N;
        let negative = ctx.control.fast.neg(positive);
        for (sign_index, point, coefficient) in [
            (0, positive, powers[exponent]),
            (1, negative, (R - powers[exponent]) % R),
        ] {
            let record = 2 * index + sign_index;
            if coords[record] != json!(point_words(point)?)
                || labels[record] != json!([column, coefficient])
            {
                return Err(format!("source archive member mismatch at point {record}"));
            }
        }
    }
    Ok(())
}

fn run() -> Result<Value, String> {
    let total = Instant::now();
    let input_started = Instant::now();
    let ctx = Context::load()?;
    let input_ms = input_started.elapsed().as_secs_f64() * 1_000.0;
    let powers = ctx.lambda_powers();
    let generator = ctx.control.fast.lift(ctx.kc.generator());
    if !ctx.control.fast.is_on_curve(generator)
        || ctx.control.fast.mul_u64(generator, R) != FastPoint::INFINITY
        || ctx.control.fast.mul_u64(generator, ctx.lambda)
            != ctx.control.fast.frobenius_k(generator, 1)
        || ((powers[N - 1] as u128 * ctx.lambda as u128) % R as u128) as u64 != 1
    {
        return fail("source generator, Frobenius eigenvalue or order-37 action failed");
    }

    let original_started = Instant::now();
    let source_seeds = ctx.source_seeds()?;
    let mut original = Vec::with_capacity(K * N);
    for seed in &source_seeds {
        let mut current = *seed;
        for _ in 0..N {
            original.push(current);
            current = ctx.control.fast.frobenius_k(current, 1);
        }
        if current != *seed {
            return fail("source orbit did not close after 37 steps");
        }
    }
    source_base_matches_archive(&ctx, &original, &powers)?;
    let original_ms = original_started.elapsed().as_secs_f64() * 1_000.0;

    let transported_started = Instant::now();
    let transported: Vec<FastPoint> = original
        .iter()
        .map(|&point| map_to_leaf(&ctx, point))
        .collect::<Result<_, _>>()?;
    for (index, &point) in original.iter().enumerate() {
        if map_to_leaf(&ctx, ctx.control.fast.neg(point))? != ctx.leaf.fast.neg(transported[index])
        {
            return Err(format!("transported sign differs at {index}"));
        }
    }
    let transported_ms = transported_started.elapsed().as_secs_f64() * 1_000.0;

    let native_started = Instant::now();
    let native_seeds = ctx.native_seeds()?;
    let mut native = Vec::with_capacity(K * N);
    for seed in &native_seeds {
        for (exponent, &power) in powers.iter().enumerate() {
            native.push(if exponent == 0 {
                *seed
            } else {
                ctx.leaf.fast.mul_u64(*seed, power)
            });
        }
    }
    let native_ms = native_started.elapsed().as_secs_f64() * 1_000.0;

    let pullback_started = Instant::now();
    let native_rows = ctx.native["rows"].as_array().ok_or("native rows")?;
    let mut pullback_seeds = Vec::with_capacity(K);
    for (index, row) in native_rows.iter().enumerate() {
        let seed = native_seeds[index];
        let archive_preimage = velu_preimage_point(&ctx.archived_source, &ctx.isogeny, seed)?;
        let control_preimage = ctx
            .bridge
            .inverse_point(archive_preimage)
            .ok_or("inverse field bridge point")?;
        if archive_preimage != point_from_value(&row["archive_pullback"])?
            || control_preimage != point_from_value(&row["control_pullback"])?
            || map_to_leaf(&ctx, control_preimage)? != seed
        {
            return Err(format!("native seed {index} has incorrect exact pullback"));
        }
        pullback_seeds.push(control_preimage);
    }
    let mut pullback = Vec::with_capacity(K * N);
    for seed in &pullback_seeds {
        let mut current = *seed;
        for _ in 0..N {
            pullback.push(current);
            current = ctx.control.fast.frobenius_k(current, 1);
        }
        if current != *seed {
            return fail("pullback source orbit did not close after 37 steps");
        }
    }
    for (index, &point) in pullback.iter().enumerate() {
        if map_to_leaf(&ctx, point)? != native[index] {
            return Err(format!("native and pullback full points differ at {index}"));
        }
        if map_to_leaf(&ctx, ctx.control.fast.neg(point))? != ctx.leaf.fast.neg(native[index]) {
            return Err(format!("native and pullback signs differ at {index}"));
        }
    }
    let pullback_ms = pullback_started.elapsed().as_secs_f64() * 1_000.0;

    let exceptions_started = Instant::now();
    let infinity_image = map_to_leaf(&ctx, FastPoint::INFINITY)?;
    let infinity_preimage =
        velu_preimage_point(&ctx.archived_source, &ctx.isogeny, FastPoint::INFINITY)?;
    let torsion = ctx.control.points_with_x(0)[0];
    let torsion_image = map_to_leaf(&ctx, torsion)?;
    let torsion_preimage_archive =
        velu_preimage_point(&ctx.archived_source, &ctx.isogeny, torsion_image)?;
    let torsion_preimage = ctx
        .bridge
        .inverse_point(torsion_preimage_archive)
        .ok_or("inverse bridge on 2-torsion")?;
    if infinity_image != FastPoint::INFINITY
        || infinity_preimage != FastPoint::INFINITY
        || torsion.infinity
        || torsion_image.infinity
        || !ctx.control.fast.double(torsion).infinity
        || !ctx.leaf.fast.double(torsion_image).infinity
        || ctx.control.fast.mul_u64(torsion, R).infinity
        || ctx.leaf.fast.mul_u64(torsion_image, R).infinity
        || torsion_preimage != torsion
    {
        return fail("infinity or 2-torsion exceptional control failed");
    }
    let exceptions_ms = exceptions_started.elapsed().as_secs_f64() * 1_000.0;

    let validation_started = Instant::now();
    let policies = vec![
        validate_policy("original", &ctx.control, &original, &powers)?,
        validate_policy("transported", &ctx.leaf, &transported, &powers)?,
        validate_policy("descendant_native", &ctx.leaf, &native, &powers)?,
        validate_policy("pullback", &ctx.control, &pullback, &powers)?,
    ];
    let validation_ms = validation_started.elapsed().as_secs_f64() * 1_000.0;
    Ok(json!({
        "schema": "n37-four-policy-support-v1",
        "status": "PASS",
        "protocol_commit": "bc19274f",
        "curve_slug": "icv1-f2m37-tm534059-32aad96b",
        "source_sha256": common::SOURCE_SHA,
        "source_base_hash": common::SOURCE_BASE_HASH,
        "native_sha256": common::NATIVE_SHA,
        "archive_sha256": common::ARCHIVE_SHA,
        "r": R,
        "cofactor": common::ORDER / R,
        "lambda": ctx.lambda,
        "lambda_powers": powers,
        "policy_counts": {
            "columns_each": K,
            "signed_classes_each": K * N,
            "physical_points_each": 2 * K * N,
        },
        "operation_counts": {
            "original_frobenius_steps": K * N,
            "transported_map_calls": K * N,
            "transported_negation_audit_map_calls": K * N,
            "native_scalar_actions": K * (N - 1),
            "native_pullback_inverse_calls": K,
            "native_pullback_audit_map_calls": 2 * K * N + K,
            "pullback_frobenius_steps": K * N,
            "subgroup_scalar_checks": 4 * 2 * K * N,
            "exception_map_calls": 2,
            "exception_inverse_calls": 2,
        },
        "cost_status": "stage diagnostic; field and scalar internals unpriced, cold source selection excluded; S and rho ratio unset",
        "memory": {
            "retained_point_array_bytes_lower_bound": 4 * K * N * std::mem::size_of::<FastPoint>(),
            "peak_rss_before_output_bytes": peak_rss_before_output_bytes(),
        },
        "phases_ms_descriptive": {
            "input": input_ms,
            "original": original_ms,
            "transported": transported_ms,
            "descendant_native": native_ms,
            "pullback": pullback_ms,
            "exceptions": exceptions_ms,
            "validation": validation_ms,
            "total": total.elapsed().as_secs_f64() * 1_000.0,
        },
        "host": {"os": std::env::consts::OS, "arch": std::env::consts::ARCH},
        "exceptions": {"infinity_image": "infinity", "infinity_preimage": "infinity", "source_two_torsion": point_words(torsion)?, "leaf_two_torsion": point_words(torsion_image)?},
        "policies": policies,
    }))
}

fn main() {
    let mut args = std::env::args().skip(1);
    let output = args
        .next()
        .expect("usage: n37_four_policy_support OUTPUT.json");
    assert!(args.next().is_none(), "unexpected extra argument");
    let value = match run() {
        Ok(value) => value,
        Err(error) => {
            json!({"schema":"n37-four-policy-support-v1", "status":"FAIL", "error":error})
        }
    };
    let bytes = serde_json::to_vec_pretty(&value).expect("serialize support receipt");
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&output)
        .expect("refuse to overwrite support receipt");
    file.write_all(&bytes).expect("write support receipt");
    file.write_all(b"\n").expect("terminate support receipt");
    let digest = hex::encode(sha256(
        &std::fs::read(&output).expect("read support receipt"),
    ));
    println!("{} {} {}", value["status"], digest, output);
    if value["status"] != "PASS" {
        std::process::exit(1);
    }
}
