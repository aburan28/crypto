//! Replay-bound bridge from the stored primary N83 points to the generic IC base.
//!
//! The library's explicit-orbit constructor currently requires a u64 field.
//! This bridge uses the verified JSONL orbit order directly and checks every
//! point and signed-Frobenius label before building the public structure.

use super::{
    curve, general, load_object, parse_point, point_json, point_set_hash, words, Result, SEEDS,
    STUDY,
};
use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    koblitz_point_count, point_key, FactorBaseDomain, FrobeniusFactorBase, KoblitzCurve,
};
use num_bigint::BigUint;
use num_traits::One;
use serde_json::{json, Value};
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

pub(super) struct PrimaryBase {
    pub(super) curve: KoblitzCurve,
    pub(super) factor_base: FrobeniusFactorBase,
}

fn from_record(row: &Value, bytes: &[u8]) -> Result<PrimaryBase> {
    let mut lines = std::str::from_utf8(bytes)?.lines();
    let header: Value = serde_json::from_str(lines.next().ok_or("missing base header")?)?;
    let columns = row["columns"].as_u64().ok_or("column count")? as usize;
    if columns == 0 || header["orbit_columns"].as_u64() != Some(columns as u64) {
        return Err("primary orbit-column count mismatch".into());
    }
    let pinned = curve(0)?;
    if row["a"].as_u64() != Some(0)
        || header["curve_a"].as_u64() != Some(0)
        || header["schema"] != "n83.factor-base/v1"
        || header["study"] != STUDY
        || header["curve_slug"] != "icv1-f2m83-tm6151469093347-debefd74"
        || header["n"] != 83
        || header["policy"] != row["policy"]
        || header["seed"] != row["seed"]
        || header["closure"] != "signed_frobenius"
        || header["modulus_low_terms"] != json!([0, 1, 2, 45])
        || header["subgroup_order"] != pinned.subgroup_order.to_string()
        || header["cofactor"] != pinned.cofactor.to_string()
        || header["frobenius_lambda"] != pinned.lambda.to_string()
        || header["generator"] != point_json(words(pinned.generator()))
        || header["point_count"].as_u64() != Some((166 * columns) as u64)
        || row["points"].as_u64() != Some((166 * columns) as u64)
    {
        return Err("primary curve or manifest binding mismatch".into());
    }
    let reps = header["representatives"]
        .as_array()
        .ok_or("representatives")?;
    if reps.len() != columns {
        return Err("representative count mismatch".into());
    }

    let mut points = Vec::with_capacity(166 * columns);
    let mut point_words = Vec::with_capacity(166 * columns);
    let mut orbits = Vec::with_capacity(2 * columns);
    let mut orbit_of = Vec::with_capacity(166 * columns);
    let mut signed_orbits = Vec::with_capacity(columns);
    let mut signed_orbit_of = Vec::with_capacity(166 * columns);
    let mut x_values = BTreeMap::new();
    let mut seen = HashSet::with_capacity(166 * columns);

    for (column, value) in reps.iter().enumerate() {
        let rep_word = parse_point(value)?.ok_or("identity representative")?;
        let rep = general(Some(rep_word));
        if !pinned.curve.is_on_curve(&rep)
            || pinned.mul(&rep, &pinned.subgroup_order) != BinaryPoint::Infinity
            || pinned.mul(&rep, &pinned.lambda) != pinned.frobenius(&rep)
        {
            return Err("primary representative subgroup or eigenvalue mismatch".into());
        }
        let mut current = rep.clone();
        let mut coefficient = BigUint::one();
        let mut positive_orbit = Vec::with_capacity(83);
        let mut negative_orbit = Vec::with_capacity(83);
        let mut signed_members = Vec::with_capacity(166);
        for phase in 0..83u32 {
            if let BinaryPoint::Affine { x, .. } = &current {
                x_values.entry(x.to_biguint()).or_insert_with(|| x.clone());
            } else {
                return Err("identity inside primary orbit".into());
            }
            for negative in [false, true] {
                let point = if negative {
                    point_neg(&current)
                } else {
                    current.clone()
                };
                let entry: Value = serde_json::from_str(lines.next().ok_or("missing point row")?)?;
                let expected_coefficient = if negative {
                    (&pinned.subgroup_order - &coefficient).to_string()
                } else {
                    coefficient.to_string()
                };
                let index = points.len();
                if entry["point"] != point_json(words(&point))
                    || entry["column"].as_u64() != Some(column as u64)
                    || entry["phase"].as_u64() != Some(phase as u64)
                    || entry["negative"].as_bool() != Some(negative)
                    || entry["coefficient"] != expected_coefficient
                    || !seen.insert(point_key(&point))
                {
                    return Err("primary point, label or distinctness mismatch".into());
                }
                point_words.push(words(&point));
                points.push(point);
                orbit_of.push((2 * column + usize::from(negative), phase));
                signed_orbit_of.push((column, phase, negative));
                signed_members.push(index);
                if negative {
                    negative_orbit.push(index);
                } else {
                    positive_orbit.push(index);
                }
            }
            current = pinned.frobenius(&current);
            coefficient = (&coefficient * &pinned.lambda) % &pinned.subgroup_order;
        }
        if current != rep {
            return Err("primary orbit did not close".into());
        }
        orbits.push(positive_orbit);
        orbits.push(negative_orbit);
        signed_orbits.push(signed_members);
    }
    if lines.next().is_some()
        || points.len() != 166 * columns
        || x_values.len() != 83 * columns
        || row["point_set_blake3"] != point_set_hash(point_words)?
    {
        return Err("primary point count or set hash mismatch".into());
    }

    let curve = KoblitzCurve {
        a: 0,
        n: 83,
        trace: -1,
        group_order: koblitz_point_count(0, 83),
        subgroup_order: pinned.subgroup_order,
        cofactor: pinned.cofactor,
        lambda: pinned.lambda,
        frobenius_is_endomorphism: true,
        curve: pinned.curve,
    };
    let factor_base = FrobeniusFactorBase {
        domain: FactorBaseDomain::ExplicitFrobeniusOrbits {
            representatives: columns,
        },
        ell: 83,
        f_j: 0,
        linearised_exponents: Vec::new(),
        subspace: x_values.into_values().collect(),
        subspace_basis: (0..83)
            .map(|bit| F2mElement::from_bit_positions(&[bit], 83))
            .collect(),
        points,
        orbits,
        orbit_of,
        signed_orbits,
        signed_orbit_of,
    };
    Ok(PrimaryBase { curve, factor_base })
}

pub(super) fn check_panel(root: &Path, columns: usize) -> Result<Value> {
    let manifest_bytes = fs::read(root.join("manifest.json"))?;
    let manifest: Value = serde_json::from_slice(&manifest_bytes)?;
    let replay: Value = serde_json::from_slice(&fs::read(root.join("replay.json"))?)?;
    if manifest["schema"] != "n83.factor-base-panel/v1"
        || manifest["study"] != STUDY
        || manifest["status"] != "completed_factor_base_panel"
        || manifest["completed_base_count"] != 54
        || replay["schema"] != "n83.factor-base-replay/v1"
        || replay["status"] != "PASS"
        || replay["panel_manifest_blake3"] != blake3::hash(&manifest_bytes).to_hex().to_string()
    {
        return Err("panel or replay receipt is not complete and bound".into());
    }
    let rows = manifest["bases"].as_array().ok_or("panel base rows")?;
    if rows.len() != 54 {
        return Err("panel base count mismatch".into());
    }
    let mut rows = rows.iter().filter(|row| {
        row["a"].as_u64() == Some(0)
            && row["policy"] == "public_x_hash"
            && row["seed"].as_u64() == Some(SEEDS[0])
            && row["columns"].as_u64() == Some(columns as u64)
    });
    let row = rows.next().ok_or("requested primary base absent")?;
    if rows.next().is_some() {
        return Err("ambiguous primary base selection".into());
    }
    let checks = replay["checks"].as_array().ok_or("replay checks")?;
    if !checks
        .iter()
        .any(|check| check["object"] == row["object"] && check["status"] == "PASS")
    {
        return Err("selected object lacks a successful replay".into());
    }
    let base = from_record(row, &load_object(root, row)?)?;
    Ok(json!({
        "schema": "n83.primary-factor-base-adapter/v1",
        "study": STUDY,
        "status": "PASS",
        "curve_a": 0,
        "subgroup_order": base.curve.subgroup_order.to_string(),
        "object": row["object"],
        "point_set_blake3": row["point_set_blake3"],
        "orbit_columns": columns,
        "points": base.factor_base.points.len(),
        "frobenius_orbits": base.factor_base.orbits.len(),
        "signed_frobenius_orbits": base.factor_base.signed_orbits.len(),
        "solver_stage_executed": false,
        "total_index_calculus_runtime_ms": Value::Null,
        "selected_best_total_runtime": Value::Null
    }))
}

#[cfg(test)]
mod tests {
    use super::super::construct;
    use super::*;
    use crypto_lib::cryptanalysis::koblitz_index_calculus::{
        koblitz_index_calculus_dlp_with_factor_base, KoblitzIcOptions,
    };
    use std::time::{Duration, Instant};

    #[test]
    fn small_primary_export_maps_to_generic_orbits_and_zero_trial_preflight() {
        let (bytes, row) = construct(
            0,
            "public_x_hash",
            2,
            17,
            Instant::now() + Duration::from_secs(60),
        )
        .unwrap();
        let base = from_record(&row, &bytes).unwrap();
        assert_eq!(base.factor_base.points.len(), 332);
        assert_eq!(base.factor_base.orbits.len(), 4);
        assert_eq!(base.factor_base.signed_orbits.len(), 2);
        assert_eq!(base.factor_base.subspace.len(), 166);
        for (index, point) in base.factor_base.points.iter().enumerate() {
            let (column, phase, negative) = base.factor_base.signed_orbit_of[index];
            assert_eq!(column, index / 166);
            assert_eq!(phase as usize, (index % 166) / 2);
            assert_eq!(negative, index % 2 == 1);
            assert_eq!(
                base.factor_base.orbit_of[index].0,
                2 * column + usize::from(negative)
            );
            assert_eq!(
                point,
                &base.factor_base.points[base.factor_base.signed_orbits[column][index % 166]]
            );
        }
        let mut opts = KoblitzIcOptions::default();
        opts.max_trials = 0;
        opts.allow_direct_relation = false;
        let report = koblitz_index_calculus_dlp_with_factor_base(
            &base.curve,
            base.curve.generator(),
            &base.factor_base,
            &opts,
        )
        .unwrap();
        assert_eq!(report.factor_base_size, 332);
        assert_eq!(report.orbit_count, 2);
        assert_eq!(report.trials, 0);
        assert!(report.log.is_none());

        let mut entries: Vec<Value> = std::str::from_utf8(&bytes)
            .unwrap()
            .lines()
            .map(|line| serde_json::from_str(line).unwrap())
            .collect();
        let mut wrong_header = entries.clone();
        wrong_header[0]["seed"] = json!(18);
        let bad_header = wrong_header
            .iter()
            .map(Value::to_string)
            .collect::<Vec<_>>()
            .join("\n");
        assert!(from_record(&row, bad_header.as_bytes()).is_err());
        entries[1]["coefficient"] = json!("0");
        let bad = entries
            .iter()
            .map(Value::to_string)
            .collect::<Vec<_>>()
            .join("\n");
        assert!(from_record(&row, bad.as_bytes()).is_err());
    }
}
