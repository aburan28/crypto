//! Independent full-point replay of the frozen n37 four-policy support census.
//! Source and pullback orbits use the general Koblitz group law; the native
//! leaf orbit is derived by map-after-pullback, then compared to scalar action.

#![recursion_limit = "256"]

#[path = "support/n37_policy_common.rs"]
mod common;

use common::{point_from_value, point_words, Context, K, N, R};
use crypto_lib::cryptanalysis::binary_velu::{velu_point_map, velu_preimage_point, Curve};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::hash::sha256::sha256;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::collections::HashSet;
use std::fs::OpenOptions;
use std::io::Write;

fn map_to_leaf(ctx: &Context, point: FastPoint) -> Result<FastPoint, String> {
    let archive_point = ctx.bridge.map_point(point).ok_or("bridge point")?;
    velu_point_map(
        archive_point,
        &ctx.isogeny.kernel,
        37,
        &ctx.archived_source.irr,
    )
    .ok_or("undefined degree-73 point map".into())
}

fn general_source_orbit(ctx: &Context, seed: FastPoint) -> Result<Vec<FastPoint>, String> {
    let mut point = ctx.control.fast.lower(seed);
    let mut orbit = Vec::with_capacity(N);
    for _ in 0..N {
        orbit.push(ctx.control.fast.lift(&point));
        point = ctx.kc.frobenius(&point);
    }
    if point != ctx.control.fast.lower(seed) {
        return Err("general source Frobenius did not close at 37".into());
    }
    Ok(orbit)
}

fn verify_entries(
    result: &Value,
    name: &str,
    curve: &Curve,
    expected: &[FastPoint],
    powers: &[u64],
) -> Result<(), String> {
    if result["name"] != name
        || result["columns"].as_u64() != Some(K as u64)
        || result["signed_classes"].as_u64() != Some((K * N) as u64)
        || result["physical_points"].as_u64() != Some((2 * K * N) as u64)
        || expected.len() != K * N
    {
        return Err(format!("{name}: policy metadata mismatch"));
    }
    let entries = result["entries"].as_array().ok_or("policy entries")?;
    if entries.len() != 2 * K * N {
        return Err(format!("{name}: entry count mismatch"));
    }
    let mut classes = HashSet::new();
    for (index, &positive) in expected.iter().enumerate() {
        let column = index / N;
        let exponent = index % N;
        let negative = curve.fast.neg(positive);
        for (offset, point, sign, coefficient) in [
            (0usize, positive, 1i64, powers[exponent]),
            (1usize, negative, -1i64, (R - powers[exponent]) % R),
        ] {
            let entry = &entries[2 * index + offset];
            if entry["point"] != json!(point_words(point)?)
                || entry["column"].as_u64() != Some(column as u64)
                || entry["exponent"].as_u64() != Some(exponent as u64)
                || entry["sign"].as_i64() != Some(sign)
                || entry["coefficient"].as_u64() != Some(coefficient)
                || point.infinity
                || point.x == 0
                || !curve.fast.is_on_curve(point)
                || !curve.fast.mul_u64(point, R).infinity
            {
                return Err(format!(
                    "{name}: bad full point or label at {column}/{exponent}/{sign}"
                ));
            }
        }
        if !classes.insert(positive.x) {
            return Err(format!(
                "{name}: duplicate signed class at {column}/{exponent}"
            ));
        }
    }
    if classes.len() != K * N {
        return Err(format!("{name}: useful support shortfall"));
    }
    Ok(())
}

fn audit(manifest: &Value) -> Result<Value, String> {
    let ctx = Context::load()?;
    if manifest["schema"] != "n37-four-policy-support-v1"
        || manifest["status"] != "PASS"
        || manifest["protocol_commit"] != "bc19274f"
        || manifest["source_sha256"] != common::SOURCE_SHA
        || manifest["native_sha256"] != common::NATIVE_SHA
        || manifest["archive_sha256"] != common::ARCHIVE_SHA
        || manifest["source_base_hash"] != common::SOURCE_BASE_HASH
        || manifest["r"].as_u64() != Some(R)
        || manifest["lambda"].as_u64() != Some(ctx.lambda)
    {
        return Err("manifest identity does not match frozen inputs".into());
    }
    let expected_counts = json!({
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
    });
    let retained_lower_bound = 4 * K * N * std::mem::size_of::<FastPoint>();
    if manifest["operation_counts"] != expected_counts
        || manifest["policy_counts"]
            != json!({"columns_each": K, "signed_classes_each": K * N, "physical_points_each": 2 * K * N})
        || manifest["memory"]["retained_point_array_bytes_lower_bound"].as_u64()
            != Some(retained_lower_bound as u64)
        || manifest["memory"]["peak_rss_before_output_bytes"]
            .as_u64()
            .is_none_or(|rss| rss < retained_lower_bound as u64)
    {
        return Err("manifest counters or memory accounting do not match the policy".into());
    }
    let mut powers = Vec::with_capacity(N);
    let mut value = 1u64;
    for _ in 0..N {
        powers.push(value);
        value = ((value as u128 * ctx.lambda as u128) % R as u128) as u64;
    }
    if value != 1 || ctx.lambda == 1 || manifest["lambda_powers"] != json!(powers) {
        return Err("lambda does not have the claimed subgroup orbit".into());
    }
    let generator = ctx.kc.generator();
    if ctx.kc.mul(generator, &BigUint::from(ctx.lambda)) != ctx.kc.frobenius(generator) {
        return Err("general group law rejects Frobenius eigenvalue".into());
    }

    let source_seeds = ctx.source_seeds()?;
    let archived_points = ctx.source["factor_base_point_coordinates"]
        .as_array()
        .ok_or("source archive points")?;
    let archived_labels = ctx.source["factor_base_point_labels"]
        .as_array()
        .ok_or("source archive labels")?;
    if archived_points.len() != 2 * K * N || archived_labels.len() != archived_points.len() {
        return Err("source archived point or label count changed".into());
    }
    let mut original = Vec::with_capacity(K * N);
    for &seed in &source_seeds {
        original.extend(general_source_orbit(&ctx, seed)?);
    }
    for (index, &point) in original.iter().enumerate() {
        let column = index / N;
        let exponent = index % N;
        let negative = ctx.control.fast.neg(point);
        if archived_points[2 * index] != json!(point_words(point)?)
            || archived_points[2 * index + 1] != json!(point_words(negative)?)
            || archived_labels[2 * index] != json!([column, powers[exponent]])
            || archived_labels[2 * index + 1] != json!([column, (R - powers[exponent]) % R])
        {
            return Err(format!("source archived header differs at {index}"));
        }
    }

    let transported = original
        .iter()
        .map(|&point| map_to_leaf(&ctx, point))
        .collect::<Result<Vec<_>, _>>()?;
    for (index, &point) in original.iter().enumerate() {
        if map_to_leaf(&ctx, ctx.control.fast.neg(point))? != ctx.leaf.fast.neg(transported[index])
        {
            return Err(format!("source-to-leaf sign mismatch at {index}"));
        }
    }
    let native_seeds = ctx.native_seeds()?;
    let native_rows = ctx.native["rows"].as_array().ok_or("native rows")?;
    let mut pullback = Vec::with_capacity(K * N);
    let mut native = Vec::with_capacity(K * N);
    for (column, &leaf_seed) in native_seeds.iter().enumerate() {
        let archive_preimage = velu_preimage_point(&ctx.archived_source, &ctx.isogeny, leaf_seed)?;
        let source_preimage = ctx
            .bridge
            .inverse_point(archive_preimage)
            .ok_or("inverse basis map")?;
        if archive_preimage != point_from_value(&native_rows[column]["archive_pullback"])?
            || source_preimage != point_from_value(&native_rows[column]["control_pullback"])?
            || map_to_leaf(&ctx, source_preimage)? != leaf_seed
        {
            return Err(format!("independent pullback mismatch at seed {column}"));
        }
        let orbit = general_source_orbit(&ctx, source_preimage)?;
        for (exponent, point) in orbit.into_iter().enumerate() {
            let image = map_to_leaf(&ctx, point)?;
            let scalar_image = ctx.leaf.fast.mul_u64(leaf_seed, powers[exponent]);
            if image != scalar_image {
                return Err(format!(
                    "native action differs from mapped source at {column}/{exponent}"
                ));
            }
            if map_to_leaf(&ctx, ctx.control.fast.neg(point))? != ctx.leaf.fast.neg(image) {
                return Err(format!(
                    "native pullback sign differs at {column}/{exponent}"
                ));
            }
            pullback.push(point);
            native.push(image);
        }
    }

    let policies = manifest["policies"].as_array().ok_or("policy table")?;
    if policies.len() != 4 {
        return Err("manifest lacks four policies".into());
    }
    verify_entries(&policies[0], "original", &ctx.control, &original, &powers)?;
    verify_entries(
        &policies[1],
        "transported",
        &ctx.leaf,
        &transported,
        &powers,
    )?;
    verify_entries(
        &policies[2],
        "descendant_native",
        &ctx.leaf,
        &native,
        &powers,
    )?;
    verify_entries(&policies[3], "pullback", &ctx.control, &pullback, &powers)?;

    let infinity = FastPoint::INFINITY;
    let torsion = ctx.control.points_with_x(0)[0];
    let torsion_image = map_to_leaf(&ctx, torsion)?;
    let torsion_preimage_archive =
        velu_preimage_point(&ctx.archived_source, &ctx.isogeny, torsion_image)?;
    let torsion_preimage = ctx
        .bridge
        .inverse_point(torsion_preimage_archive)
        .ok_or("inverse basis map on 2-torsion")?;
    if map_to_leaf(&ctx, infinity)? != infinity
        || velu_preimage_point(&ctx.archived_source, &ctx.isogeny, infinity)? != infinity
        || torsion_image.infinity
        || !ctx.control.fast.double(torsion).infinity
        || !ctx.leaf.fast.double(torsion_image).infinity
        || torsion_preimage != torsion
        || manifest["exceptions"]["source_two_torsion"] != json!(point_words(torsion)?)
        || manifest["exceptions"]["leaf_two_torsion"] != json!(point_words(torsion_image)?)
    {
        return Err("exceptional controls differ from manifest".into());
    }
    Ok(json!({
        "schema": "n37-four-policy-support-replay-v1",
        "status": "PASS",
        "recomputed_policies": 4,
        "checked_full_points": 4 * 2 * K * N,
        "checked_signed_classes": 4 * K * N,
        "checked_map_pairs": 4 * K * N,
        "checked_pullback_seeds": K,
        "checked_exception_pullbacks": 2,
        "checked_native_scalar_action_pairs": K * N,
        "source_orbit_arithmetic": "general KoblitzCurve Frobenius/group law",
        "native_orbit_arithmetic": "map of general-law pullback, checked against leaf scalar multiplication",
    }))
}

fn main() {
    let mut args = std::env::args().skip(1);
    let input = args
        .next()
        .expect("usage: n37_four_policy_support_replay RESULT.json RECEIPT.json");
    let output = args.next().expect("missing replay receipt path");
    assert!(args.next().is_none(), "unexpected extra argument");
    let bytes = std::fs::read(&input).expect("read producer manifest");
    let digest = hex::encode(sha256(&bytes));
    let manifest: Value = serde_json::from_slice(&bytes).expect("parse producer manifest");
    let receipt = match audit(&manifest) {
        Ok(mut value) => {
            value["manifest_sha256"] = json!(digest);
            value
        }
        Err(error) => json!({
            "schema":"n37-four-policy-support-replay-v1",
            "status":"FAIL",
            "manifest_sha256":digest,
            "error":error,
        }),
    };
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&output)
        .expect("refuse to overwrite replay receipt");
    file.write_all(&serde_json::to_vec_pretty(&receipt).expect("serialize replay receipt"))
        .expect("write replay receipt");
    file.write_all(b"\n").expect("terminate replay receipt");
    println!("{} {}", receipt["status"], output);
    if receipt["status"] != "PASS" {
        std::process::exit(1);
    }
}
