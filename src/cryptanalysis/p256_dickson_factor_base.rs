//! Native Dickson-torus factor-base dumps for the registered P-256 model.
//!
//! `ecbench.factor_base_dump/v1` is deliberately word-sized: its curve,
//! point-key and coefficient types are `u64`.  This module keeps its object
//! layout and exact `ecbench.factor_base/v1` FB1 preimage while widening the
//! unavoidable integer fields.  The labelled dump schema is therefore
//! `ecbench.factor_base_dump/v1-wide`, never ordinary v1.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};
use serde_json::json;

use crate::cryptanalysis::curve_id;
use crate::cryptanalysis::ecbench::canonical::{sha256_hex, short_id};
use crate::ecc::curve::CurveParams;

pub const CURVE_SLUG: &str = "icv1-fp256-t89188191154553853111372247798585809583-f188c491";
pub const CURVE_ICV1: &str =
    "ICV1:fp-115792089210356248762697446949407573530086143415290314195533631308867097853951:89188191154553853111372247798585809583:115792089210356248762697446949407573529996955224135760342422259061068512044369:7958909377132088453074743217357398615041065282494610304372115906626967530147:unk:unk:r:f188c4914bb7";
pub const CURVE_EC1: &str = "EC1P256Cp256h0523b774e066";
pub const CURVE_UID: &str =
    "urn:ec-record:1:sha256:0523b774e0666cf37ff4dc9cf452eae024192a7f89bdad10ca226d6e64d4dd70";
pub const DUMP_SCHEMA: &str = "ecbench.factor_base_dump/v1-wide";
pub const FACTOR_BASE_SCHEMA: &str = "ecbench.factor_base/v1";
pub const FAMILY: &str = "dickson-torus";
pub const POINT_KEY_ENCODING: &str = "prime-affine-x-plus-one-shift-sign/be33/v1";
pub const POINT_KEY_BYTES: usize = 33;
pub const MAX_ENUMERATION_DEPTH: u32 = 24;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WideCurveFacts {
    pub slug: String,
    pub icv1: String,
    pub family: String,
    pub construction: String,
    pub field_degree: Option<u32>,
    pub field_bits: u32,
    pub group_order: String,
    pub r: String,
    pub cofactor: String,
    pub generator: [String; 2],
    pub automorphisms_available: u32,
    pub registered: bool,
    pub ec1: String,
    pub curve_uid: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WideFactorBaseFacts {
    pub fb_id: String,
    pub fb_sha256: String,
    pub family: String,
    pub params: BTreeMap<String, String>,
    pub description: String,
    pub signed_points: u64,
    pub abscissae: u64,
    pub columns: u64,
    pub dimension: Option<u32>,
    pub points_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WideFbPoint {
    pub x: String,
    pub y: String,
    pub col: u64,
    /// Decimal modulo the subgroup order; a string because `r` is 256-bit.
    pub coef: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WideFactorBaseDump {
    pub schema: String,
    pub curve: WideCurveFacts,
    pub factor_base: WideFactorBaseFacts,
    pub build_adds: u64,
    pub build_doubles: u64,
    pub points: Vec<WideFbPoint>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct WideBuildResult {
    pub dump: WideFactorBaseDump,
    pub legacy_enumeration_sha256: String,
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct RelationMetrics {
    pub columns: u64,
    pub relation_length: u32,
    pub signed_domain: String,
    pub poisson_mean: f64,
    pub poisson_success: f64,
    /// Bound for the saturated selector ideal; not a measured F4 degree.
    pub selector_ideal_regularity_bound: u32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Fp2 {
    re: BigUint,
    im: BigUint,
}

impl Fp2 {
    fn one() -> Self {
        Self {
            re: BigUint::one(),
            im: BigUint::zero(),
        }
    }

    fn mul(&self, other: &Self, p: &BigUint) -> Self {
        let ac = (&self.re * &other.re) % p;
        let bd = (&self.im * &other.im) % p;
        let ad = (&self.re * &other.im) % p;
        let bc = (&self.im * &other.re) % p;
        let re = if ac >= bd { ac - bd } else { ac + p - bd };
        Self {
            re,
            im: (ad + bc) % p,
        }
    }

    fn pow(&self, exponent: &BigUint, p: &BigUint) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();
        let mut e = exponent.clone();
        while !e.is_zero() {
            if e.bit(0) {
                result = result.mul(&base, p);
            }
            base = base.mul(&base, p);
            e >>= 1usize;
        }
        result
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Column {
    x: BigUint,
    low_y: BigUint,
}

pub(crate) fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn fixed_hex(value: &BigUint, digits: usize) -> String {
    format!("0x{:0>width$}", value.to_str_radix(16), width = digits)
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("not a non-negative integer: `{value}`"))
}

fn parse_spec(spec: &str) -> Result<(u32, Option<BigUint>), String> {
    let (family, tail) = spec.split_once(':').unwrap_or((spec, ""));
    if family != FAMILY {
        return Err(format!("wide P-256 builder is `{FAMILY}`, not `{family}`"));
    }
    let mut depth = None;
    let mut root = None;
    for item in tail.split(',').filter(|s| !s.is_empty()) {
        let (name, value) = item
            .split_once('=')
            .ok_or_else(|| format!("factor-base parameter needs `name=value`: `{item}`"))?;
        match name {
            "depth" => {
                depth = Some(
                    value
                        .parse::<u32>()
                        .map_err(|_| format!("depth is not an integer: `{value}`"))?,
                )
            }
            "root_exponent" => root = Some(parse_big(value)?),
            other => return Err(format!("unknown `{FAMILY}` parameter `{other}`")),
        }
    }
    let depth = depth.ok_or("dickson-torus needs `depth`".to_string())?;
    if !(1..=MAX_ENUMERATION_DEPTH).contains(&depth) {
        return Err(format!(
            "materialised depth {depth} outside 1..={MAX_ENUMERATION_DEPTH}"
        ));
    }
    Ok((depth, root))
}

pub(crate) fn curve_facts(curve: &CurveParams, id: &curve_id::CurveId) -> WideCurveFacts {
    WideCurveFacts {
        slug: id.slug.clone(),
        icv1: id.icv1.clone(),
        family: "prime".into(),
        construction: "CurveParams::p256()".into(),
        field_degree: None,
        field_bits: curve.p.bits() as u32,
        group_order: curve.n.to_string(),
        r: curve.n.to_string(),
        cofactor: curve.h.to_string(),
        generator: [lower_hex(&curve.gx), lower_hex(&curve.gy)],
        automorphisms_available: 2,
        registered: true,
        ec1: CURVE_EC1.into(),
        curve_uid: CURVE_UID.into(),
    }
}

fn norm_one_2_primary_generator(p: &BigUint) -> Result<Fp2, String> {
    let seed = Fp2 {
        re: BigUint::one(),
        im: BigUint::from(5u8),
    };
    let torus_generator = seed.pow(&(p - BigUint::one()), p);
    let p_plus_one = p + BigUint::one();
    let odd_part = &p_plus_one >> 96usize;
    let g96 = torus_generator.pow(&odd_part, p);
    let order = BigUint::one() << 96usize;
    if g96.pow(&order, p) != Fp2::one() || g96.pow(&(BigUint::one() << 95usize), p) == Fp2::one() {
        return Err("the norm-one element does not have exact order 2^96".into());
    }
    Ok(g96)
}

pub(crate) fn curve_y(x: &BigUint, curve: &CurveParams) -> Option<BigUint> {
    let p = &curve.p;
    let rhs = (((x * x) % p) * x + &curve.a * x + &curve.b) % p;
    if rhs.is_zero() {
        return Some(BigUint::zero());
    }
    if rhs.modpow(&((p - BigUint::one()) >> 1usize), p) != BigUint::one() {
        return None;
    }
    let y = rhs.modpow(&((p + BigUint::one()) >> 2usize), p);
    if (&y * &y) % p != rhs {
        return None;
    }
    Some(y)
}

pub(crate) fn point_key_bytes(
    x: &BigUint,
    high_sign: bool,
) -> Result<[u8; POINT_KEY_BYTES], String> {
    let key = ((x + BigUint::one()) << 1usize) + BigUint::from(high_sign);
    let raw = key.to_bytes_be();
    if raw.len() > POINT_KEY_BYTES {
        return Err("P-256 point key exceeded 33 bytes".into());
    }
    let mut out = [0u8; POINT_KEY_BYTES];
    out[POINT_KEY_BYTES - raw.len()..].copy_from_slice(&raw);
    Ok(out)
}

fn index_bytes(index: usize, width: usize) -> Vec<u8> {
    let raw = index.to_be_bytes();
    raw[raw.len() - width..].to_vec()
}

/// Build a factor base from a plug-in-style specification such as
/// `dickson-torus:depth=18`.
pub fn build(spec: &str) -> Result<WideBuildResult, String> {
    let (depth, requested_root) = parse_spec(spec)?;
    let curve = CurveParams::p256();
    let id = curve_id::prime(&curve.p, &curve.a, &curve.b, &curve.n)
        .ok_or("P-256 identity construction failed".to_string())?;
    if id.slug != CURVE_SLUG || id.icv1 != CURVE_ICV1 {
        return Err(format!(
            "registered identity mismatch: got {} / {}",
            id.slug, id.icv1
        ));
    }

    let g96 = norm_one_2_primary_generator(&curve.p)?;
    let root_exponent =
        requested_root.unwrap_or_else(|| BigUint::from(3u8) << (94usize - depth as usize));
    let exponent_modulus = BigUint::one() << (96usize - depth as usize);
    let root_exponent = root_exponent % exponent_modulus;
    let root = g96.pow(&root_exponent, &curve.p);
    let step = g96.pow(&(BigUint::one() << (96usize - depth as usize)), &curve.p);
    let terminal_z = root.pow(&(BigUint::one() << depth as usize), &curve.p);
    let terminal = (BigUint::from(2u8) * &terminal_z.re) % &curve.p;
    if terminal == BigUint::from(2u8) || terminal == &curve.p - BigUint::from(2u8) {
        return Err("terminal trace +/-2 gives a degenerate fibre".into());
    }

    let fibre_size = 1usize << depth;
    let index_width = usize::max(1, depth.div_ceil(8) as usize);
    let mut z = root;
    let mut columns = Vec::with_capacity(fibre_size / 2);
    let mut legacy_bytes = Vec::with_capacity(fibre_size / 2 * (index_width + 64));
    for index in 0..fibre_size {
        let x = (BigUint::from(2u8) * &z.re) % &curve.p;
        if let Some(y) = curve_y(&x, &curve) {
            if y.is_zero() {
                return Err("unexpected affine P-256 2-torsion".into());
            }
            let neg_y = &curve.p - &y;
            let low_y = if y < neg_y { y } else { neg_y };
            legacy_bytes.extend(index_bytes(index, index_width));
            let xb = x.to_bytes_be();
            legacy_bytes.extend(std::iter::repeat_n(0, 32 - xb.len()));
            legacy_bytes.extend(xb);
            let yb = low_y.to_bytes_be();
            legacy_bytes.extend(std::iter::repeat_n(0, 32 - yb.len()));
            legacy_bytes.extend(yb);
            columns.push(Column { x, low_y });
        }
        z = z.mul(&step, &curve.p);
    }
    let legacy_enumeration_sha256 = sha256_hex(&legacy_bytes);

    columns.sort_by(|a, b| a.x.cmp(&b.x));
    for pair in columns.windows(2) {
        if pair[0].x >= pair[1].x {
            return Err("factor-base abscissae are not distinct".into());
        }
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
    params.insert("depth".into(), depth.to_string());
    params.insert("point_key_encoding".into(), POINT_KEY_ENCODING.into());
    params.insert("root_exponent".into(), lower_hex(&root_exponent));
    params.insert("terminal".into(), fixed_hex(&terminal, 64));
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
            "depth-{depth} Dickson trace fibre on {CURVE_SLUG}; negation folded into one column and materialised as both signs"
        ),
        signed_points,
        abscissae: count,
        columns: count,
        dimension: None,
        points_sha256,
    };
    Ok(WideBuildResult {
        dump: WideFactorBaseDump {
            schema: DUMP_SCHEMA.into(),
            curve: curve_facts(&curve, &id),
            factor_base,
            build_adds: 0,
            build_doubles: 0,
            points,
        },
        legacy_enumeration_sha256,
    })
}

/// Rebuild and compare every field and point row of a wide dump.
pub fn verify(dump: &WideFactorBaseDump) -> Result<WideBuildResult, String> {
    if dump.schema != DUMP_SCHEMA {
        return Err(format!(
            "expected schema {DUMP_SCHEMA}, got {}",
            dump.schema
        ));
    }
    let depth = dump
        .factor_base
        .params
        .get("depth")
        .ok_or("factor base has no depth".to_string())?;
    let root = dump
        .factor_base
        .params
        .get("root_exponent")
        .ok_or("factor base has no root_exponent".to_string())?;
    let rebuilt = build(&format!("{FAMILY}:depth={depth},root_exponent={root}"))?;
    if rebuilt.dump != *dump {
        return Err("wide factor-base dump differs from native rebuild".into());
    }
    Ok(rebuilt)
}

/// Exact signed distinct-column domain, followed by its Poisson heuristic.
pub fn relation_metrics(columns: u64, relation_length: u32) -> Result<RelationMetrics, String> {
    if relation_length < 2 || u64::from(relation_length) > columns {
        return Err(format!(
            "relation length {relation_length} outside 2..={columns}"
        ));
    }
    let mut combinations = BigUint::one();
    for i in 0..u64::from(relation_length) {
        combinations *= BigUint::from(columns - i);
        combinations /= BigUint::from(i + 1);
    }
    let signed_domain = combinations << relation_length as usize;
    let order = CurveParams::p256().n;
    let poisson_mean = signed_domain
        .to_f64()
        .zip(order.to_f64())
        .map(|(d, n)| d / n)
        .ok_or("relation domain did not convert to f64".to_string())?;
    Ok(RelationMetrics {
        columns,
        relation_length,
        signed_domain: signed_domain.to_string(),
        poisson_mean,
        poisson_success: -(-poisson_mean).exp_m1(),
        selector_ideal_regularity_bound: relation_length + 1,
    })
}

fn sql_text(value: &str) -> String {
    format!("'{}'", value.replace('\'', "''"))
}

fn sql_optional_int<T: ToString>(value: Option<T>) -> String {
    value
        .map(|v| v.to_string())
        .unwrap_or_else(|| "NULL".into())
}

fn sql_real(value: f64) -> String {
    if value.is_finite() {
        format!("{value:e}")
    } else {
        "NULL".into()
    }
}

fn sql_bool(value: bool) -> &'static str {
    if value {
        "1"
    } else {
        "0"
    }
}

/// Emit an idempotent SQLite loader using the repository's canonical tables.
///
/// The wide builder remains separate from the frozen `ecbench` measurement
/// binary, but its inventory joins the same curve, factor-base, and point
/// tables.  A native rebuild is required before any SQL is returned.
pub fn sql_for_dump(dump: &WideFactorBaseDump) -> Result<String, String> {
    verify(dump)?;
    let fb = &dump.factor_base;
    let curve = &dump.curve;
    let r = BigUint::parse_bytes(curve.r.as_bytes(), 10)
        .and_then(|n| n.to_f64())
        .ok_or_else(|| format!("wide subgroup order is not decimal: {}", curve.r))?;
    let mut out = String::from(crate::cryptanalysis::ecbench::db::SCHEMA_SQL);
    if !out.ends_with('\n') {
        out.push('\n');
    }
    let _ = writeln!(
        out,
        "INSERT INTO curves VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (slug) DO NOTHING;",
        sql_text(&curve.slug),
        sql_text(&curve.icv1),
        sql_text(&curve.family),
        sql_optional_int(curve.field_degree),
        curve.field_bits,
        sql_text(&curve.group_order),
        sql_text(&curve.r),
        sql_text(&curve.cofactor),
        sql_real(r.log2()),
        curve.automorphisms_available,
        sql_real(crate::cryptanalysis::ecbench::record::floor_s(
            curve.automorphisms_available,
        )),
        sql_bool(curve.registered),
    );
    let _ = writeln!(
        out,
        "INSERT INTO curve_representations VALUES ({}, {}, {}, {}, {}) ON CONFLICT (slug, generator_x, generator_y) DO NOTHING;",
        sql_text(&curve.slug),
        sql_text(&curve.generator[0]),
        sql_text(&curve.generator[1]),
        sql_text(&curve.ec1),
        sql_text(&curve.curve_uid),
    );
    let _ = writeln!(
        out,
        "INSERT INTO curve_constructions VALUES ({}, {}) ON CONFLICT (slug, construction) DO NOTHING;",
        sql_text(&curve.slug),
        sql_text(&curve.construction),
    );
    let params = serde_json::to_string(&fb.params).map_err(|e| e.to_string())?;
    let _ = writeln!(
        out,
        "INSERT INTO factor_bases VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (fb_id) DO NOTHING;",
        sql_text(&fb.fb_id),
        sql_text(&fb.fb_sha256),
        sql_text(&curve.slug),
        sql_text(&fb.family),
        sql_text(&params),
        sql_text(&fb.description),
        fb.signed_points,
        fb.abscissae,
        fb.columns,
        sql_optional_int(fb.dimension),
        sql_text(&fb.points_sha256),
    );
    for (index, point) in dump.points.iter().enumerate() {
        let _ = writeln!(
            out,
            "INSERT INTO factor_base_points VALUES ({}, {}, {}, {}, {}, {}) ON CONFLICT (fb_id, idx) DO NOTHING;",
            sql_text(&fb.fb_id),
            index,
            sql_text(&point.x),
            sql_text(&point.y),
            point.col,
            sql_text(&point.coef),
        );
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registered_identity_is_recomputed() {
        let curve = CurveParams::p256();
        let id = curve_id::prime(&curve.p, &curve.a, &curve.b, &curve.n).unwrap();
        assert_eq!(id.slug, CURVE_SLUG);
        assert_eq!(id.icv1, CURVE_ICV1);
    }

    #[test]
    fn sampled_depth_nine_vector_matches() {
        let result = build("dickson-torus:depth=9,root_exponent=0x2d48091227bc911483d168").unwrap();
        assert_eq!(result.dump.factor_base.columns, 301);
        assert_eq!(result.dump.factor_base.signed_points, 602);
        assert_eq!(
            result.legacy_enumeration_sha256,
            "15b51c3754ab8ab57c732fc2302485959fceab2079307b66dc5bfe1863c3782d"
        );
        assert_eq!(
            result.dump.factor_base.points_sha256,
            "0052c8d8c7e10846baa3e7618539a23dbe34a9c9f4bede4713ccb5539f0a8f2d"
        );
        assert_eq!(result.dump.factor_base.fb_id, "FB1h6bcde28318a3");
        assert_eq!(verify(&result.dump).unwrap(), result);
    }

    #[test]
    fn depth_eighteen_relation_metrics_match() {
        let metrics = relation_metrics(131_239, 17).unwrap();
        assert!(
            (metrics.poisson_mean - 3.231_334_861_301_661_5).abs() < 1e-14,
            "{metrics:?}"
        );
        assert!(
            (metrics.poisson_success - 0.960_495_269_758_743_2).abs() < 1e-14,
            "{metrics:?}"
        );
        assert_eq!(metrics.selector_ideal_regularity_bound, 18);
    }

    #[test]
    fn wide_dump_writes_repository_sql() {
        let dump = build("dickson-torus:depth=3").unwrap().dump;
        let sql = sql_for_dump(&dump).unwrap();
        assert!(sql.starts_with("-- ecbench database schema"));
        assert!(sql.contains(&dump.factor_base.fb_id));
        assert!(sql.contains(CURVE_EC1));
        assert_eq!(
            sql.matches("INSERT INTO factor_base_points VALUES").count(),
            dump.points.len()
        );
    }
}
