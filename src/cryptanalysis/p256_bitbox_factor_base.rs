//! Exact affine-bitbox factor bases for the registered P-256 model.
//!
//! Membership has the low-degree presentation
//! `x = offset + sum(2^i b_i)`, `b_i^2 = b_i`, followed by the curve
//! equation.  This module uses the same wide FB1 dump and point-key rules as
//! the accepted Dickson-torus inventory.

use std::collections::BTreeMap;

use num_bigint::BigUint;
use num_traits::{One, Zero};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::json;

use crate::cryptanalysis::curve_id;
use crate::cryptanalysis::ecbench::canonical::{sha256_hex, short_id};
use crate::ecc::curve::CurveParams;

use super::p256_dickson_factor_base::{
    curve_facts, curve_y, lower_hex, point_key_bytes, RelationMetrics, WideBuildResult,
    WideFactorBaseDump, WideFactorBaseFacts, WideFbPoint, CURVE_ICV1, CURVE_SLUG, DUMP_SCHEMA,
    FACTOR_BASE_SCHEMA, POINT_KEY_BYTES, POINT_KEY_ENCODING,
};

pub const FAMILY: &str = "affine-bitbox";
pub const BIT_WIDTH: u32 = 18;
pub const CANDIDATES: u32 = 8;
pub const SELECTION: &str = "max-lifts/k-times-width/k=0..7/v1";

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CandidateCount {
    pub offset: String,
    pub columns: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct BitboxBuildResult {
    pub build: WideBuildResult,
    pub selected_offset: String,
    pub candidate_counts: Vec<CandidateCount>,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct CompactBitboxResult {
    pub schema: String,
    pub curve: String,
    pub family: String,
    pub bit_width: u32,
    pub interval_width: u64,
    pub selection: String,
    pub selected_offset: String,
    pub candidate_counts: Vec<CandidateCount>,
    pub fb_id: String,
    pub fb_sha256: String,
    pub fb1_preimage: serde_json::Value,
    pub columns: u64,
    pub signed_points: u64,
    pub points_sha256: String,
    pub relation_m17: RelationMetrics,
    pub verified: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Column {
    x: BigUint,
    low_y: BigUint,
}

fn interval(offset: u64, width: u64, curve: &CurveParams) -> Result<Vec<Column>, String> {
    if BigUint::from(offset + width) > curve.p {
        return Err("affine-bitbox interval wraps the prime field".into());
    }
    let columns: Vec<Column> = (0..width)
        .into_par_iter()
        .filter_map(|delta| {
            let x = BigUint::from(offset + delta);
            curve_y(&x, curve).and_then(|y| {
                if y.is_zero() {
                    None
                } else {
                    let neg_y = &curve.p - &y;
                    let low_y = if y < neg_y { y } else { neg_y };
                    Some(Column { x, low_y })
                }
            })
        })
        .collect();
    Ok(columns)
}

fn build_once() -> Result<BitboxBuildResult, String> {
    let curve = CurveParams::p256();
    let id = curve_id::prime(&curve.p, &curve.a, &curve.b, &curve.n)
        .ok_or("P-256 identity construction failed".to_string())?;
    if id.slug != CURVE_SLUG || id.icv1 != CURVE_ICV1 {
        return Err(format!(
            "registered identity mismatch: got {} / {}",
            id.slug, id.icv1
        ));
    }

    let width = 1u64 << BIT_WIDTH;
    let mut counts = Vec::with_capacity(CANDIDATES as usize);
    let mut selected: Option<(u64, Vec<Column>)> = None;
    for k in 0..CANDIDATES {
        let offset = u64::from(k) * width;
        let columns = interval(offset, width, &curve)?;
        counts.push(CandidateCount {
            offset: lower_hex(&BigUint::from(offset)),
            columns: columns.len() as u64,
        });
        if selected.as_ref().is_none_or(|(best_offset, best)| {
            columns.len() > best.len() || (columns.len() == best.len() && offset < *best_offset)
        }) {
            selected = Some((offset, columns));
        }
    }
    let (offset, mut columns) = selected.ok_or("no affine-bitbox candidates".to_string())?;
    columns.sort_by(|left, right| left.x.cmp(&right.x));
    if columns.windows(2).any(|pair| pair[0].x >= pair[1].x) {
        return Err("affine-bitbox abscissae are not distinct".into());
    }

    let mut point_key_blob = Vec::with_capacity(columns.len() * 2 * POINT_KEY_BYTES);
    let mut points = Vec::with_capacity(columns.len() * 2);
    for (col, column) in columns.iter().enumerate() {
        let high_y = &curve.p - &column.low_y;
        point_key_blob.extend(point_key_bytes(&column.x, false)?);
        point_key_blob.extend(point_key_bytes(&column.x, true)?);
        points.push(WideFbPoint {
            x: lower_hex(&column.x),
            y: lower_hex(&column.low_y),
            col: col as u64,
            coef: "1".into(),
        });
        points.push(WideFbPoint {
            x: lower_hex(&column.x),
            y: lower_hex(&high_y),
            col: col as u64,
            coef: (&curve.n - BigUint::one()).to_string(),
        });
    }
    let points_sha256 = sha256_hex(&point_key_blob);
    let mut params = BTreeMap::new();
    params.insert("bit_width".into(), BIT_WIDTH.to_string());
    params.insert("offset".into(), lower_hex(&BigUint::from(offset)));
    params.insert("point_key_encoding".into(), POINT_KEY_ENCODING.into());
    params.insert("selection".into(), SELECTION.into());
    let signed_points = points.len() as u64;
    let count = columns.len() as u64;
    let identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": CURVE_SLUG,
        "family": FAMILY,
        "params": params,
        "columns": count,
        "signed_points": signed_points,
        "points_sha256": points_sha256,
    });
    let (fb_id, fb_sha256) = short_id("FB1", &identity)?;
    let factor_base = WideFactorBaseFacts {
        fb_id,
        fb_sha256,
        family: FAMILY.into(),
        params,
        description: format!(
            "{BIT_WIDTH}-bit affine interval on {CURVE_SLUG}; negation folded into one column and materialised as both signs"
        ),
        signed_points,
        abscissae: count,
        columns: count,
        dimension: Some(BIT_WIDTH),
        points_sha256,
    };
    Ok(BitboxBuildResult {
        selected_offset: lower_hex(&BigUint::from(offset)),
        candidate_counts: counts,
        build: WideBuildResult {
            dump: WideFactorBaseDump {
                schema: DUMP_SCHEMA.into(),
                curve: curve_facts(&curve, &id),
                factor_base,
                build_adds: 0,
                build_doubles: 0,
                points,
            },
            legacy_enumeration_sha256: String::new(),
        },
    })
}

/// Build the frozen eight-candidate P-256 affine-bitbox factor base.
pub fn build() -> Result<BitboxBuildResult, String> {
    build_once()
}

/// Rebuild and compare every emitted field and point row.
pub fn verify(result: &BitboxBuildResult) -> Result<(), String> {
    let rebuilt = build_once()?;
    if rebuilt != *result {
        return Err("affine-bitbox result differs from native rebuild".into());
    }
    Ok(())
}

pub fn compact(result: &BitboxBuildResult, verified: bool) -> CompactBitboxResult {
    let fb = &result.build.dump.factor_base;
    let identity = json!({
        "schema": FACTOR_BASE_SCHEMA,
        "curve": CURVE_SLUG,
        "family": FAMILY,
        "params": fb.params,
        "columns": fb.columns,
        "signed_points": fb.signed_points,
        "points_sha256": fb.points_sha256,
    });
    CompactBitboxResult {
        schema: "p256.affine_bitbox_result/v1".into(),
        curve: CURVE_SLUG.into(),
        family: FAMILY.into(),
        bit_width: BIT_WIDTH,
        interval_width: 1u64 << BIT_WIDTH,
        selection: SELECTION.into(),
        selected_offset: result.selected_offset.clone(),
        candidate_counts: result.candidate_counts.clone(),
        fb_id: fb.fb_id.clone(),
        fb_sha256: fb.fb_sha256.clone(),
        fb1_preimage: identity,
        columns: fb.columns,
        signed_points: fb.signed_points,
        points_sha256: fb.points_sha256.clone(),
        relation_m17: super::p256_dickson_factor_base::relation_metrics(fb.columns, 17)
            .expect("the frozen bitbox has more than 17 columns"),
        verified,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_width_and_offsets_are_frozen() {
        assert_eq!(BIT_WIDTH, 18);
        assert_eq!(CANDIDATES, 8);
        assert_eq!(1u64 << BIT_WIDTH, 262_144);
        assert_eq!(SELECTION, "max-lifts/k-times-width/k=0..7/v1");
    }

    #[test]
    fn small_interval_points_are_sorted_and_on_curve() {
        let curve = CurveParams::p256();
        let columns = interval(0, 32, &curve).unwrap();
        assert!(columns.windows(2).all(|pair| pair[0].x < pair[1].x));
        for column in columns {
            let rhs =
                ((&column.x * &column.x % &curve.p) * &column.x + &curve.a * &column.x + &curve.b)
                    % &curve.p;
            assert_eq!((&column.low_y * &column.low_y) % &curve.p, rhs);
        }
    }
}
