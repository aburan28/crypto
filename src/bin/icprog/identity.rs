//! Canonical identities, natively: the tournament's IC1 adapter
//! (`research/ic_candidate_tournament_20260915/identity.py`) and the EC1
//! curve identity (`tools/curve_identity.py`), with the same records, the
//! same digests and the same refusals.
//!
//! A label is derived from its record and never taken as evidence of it.
//! The factor-base inventory is counted again here, with [`oracle`]'s own
//! arithmetic, from the orbit representatives a report names.

use std::collections::BTreeSet;
use std::path::Path;

use super::json::{self, J};
use super::oracle::Curve;

fn require(cond: bool, message: &str) -> Result<(), String> {
    if cond {
        Ok(())
    } else {
        Err(message.to_string())
    }
}

fn no_floats(v: &J) -> Result<(), String> {
    match v {
        J::Float(_) => Err("canonical records cannot contain floats or non-JSON values".into()),
        J::Arr(items) => items.iter().try_for_each(no_floats),
        J::Obj(kv) => kv.iter().try_for_each(|(_, v)| no_floats(v)),
        _ => Ok(()),
    }
}

/// Stable, lossless JSON: sorted keys, no spaces, raw UTF-8, no floats.
pub fn canonical(v: &J) -> Result<String, String> {
    no_floats(v)?;
    Ok(json::dumps_compact(v, true, false))
}

pub fn sha256_hex(bytes: &[u8]) -> String {
    hex::encode(crypto_lib::hash::sha256::sha256(bytes))
}

pub fn sha256(v: &J) -> Result<String, String> {
    Ok(sha256_hex(canonical(v)?.as_bytes()))
}

fn natural(v: &J, name: &str, positive: bool) -> Result<i128, String> {
    match v {
        J::Int(i) if *i >= i128::from(positive) => Ok(*i),
        _ => Err(format!("invalid {name}")),
    }
}

fn is_digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

fn digest(v: &J, name: &str) -> Result<(), String> {
    require(
        v.as_str().is_some_and(is_digest),
        &format!("invalid {name}"),
    )
}

/// The object's keys are exactly `names`.
fn fields(v: &J, names: &str, name: &str) -> Result<(), String> {
    let want: BTreeSet<&str> = names.split_whitespace().collect();
    let ok = v.as_obj().is_some_and(|kv| {
        let have: BTreeSet<&str> = kv.iter().map(|(k, _)| k.as_str()).collect();
        have == want && have.len() == kv.len()
    });
    require(ok, &format!("missing or extra fields in {name}"))
}

/// `[a-z][a-z0-9]*`.
fn compact_code(s: &str) -> bool {
    let mut cs = s.chars();
    cs.next().is_some_and(|c| c.is_ascii_lowercase())
        && cs.all(|c| c.is_ascii_lowercase() || c.is_ascii_digit())
}

/// Run and provenance fields, and paths, refused in a candidate identity.
fn identity_only(v: &J) -> Result<(), String> {
    const FORBIDDEN: [&str; 14] = [
        "candidate_id",
        "workload_id",
        "run_id",
        "algorithm_seed",
        "run_seed",
        "target_seeds",
        "timestamp",
        "created_at",
        "measurements",
        "wall_ns",
        "total_operations",
        "source_directory",
        "source_root",
        "path",
    ];
    match v {
        J::Obj(kv) => {
            for (k, v) in kv {
                require(
                    !FORBIDDEN.contains(&k.as_str())
                        && k != "paths"
                        && !k.ends_with("_path")
                        && !k.ends_with("_paths"),
                    &format!("run data or path in candidate identity: {k}"),
                )?;
                identity_only(v)?;
            }
            Ok(())
        }
        J::Arr(items) => items.iter().try_for_each(identity_only),
        J::Str(s) => require(
            !s.starts_with('/') && !s.starts_with("file://"),
            "absolute path in candidate identity",
        ),
        _ => Ok(()),
    }
}

fn point_json(p: (u64, u64)) -> J {
    J::Arr(vec![J::Int(p.0.into()), J::Int(p.1.into())])
}

fn points_json(ps: &[(u64, u64)]) -> J {
    J::Arr(ps.iter().map(|&p| point_json(p)).collect())
}

/// The curve's EC1 record and alias, as the IC1 adapter writes them.
pub fn curve_record(c: &Curve) -> Result<J, String> {
    let order = c.h as i128 * c.r as i128;
    let mut curve = vec![
        ("tag".to_string(), J::Str(format!("kb{}", c.a))),
        ("model".into(), J::Str("y^2+x*y=x^3+a2*x^2+a6".into())),
        (
            "coefficients".into(),
            json::obj([
                ("a1", J::Int(1)),
                ("a2", J::Int(c.a.into())),
                ("a3", J::Int(0)),
                ("a4", J::Int(0)),
                ("a6", J::Int(1)),
            ]),
        ),
        ("order".into(), J::Int(order)),
        ("trace".into(), J::Int((1i128 << c.n) + 1 - order)),
        ("r".into(), J::Int(c.r.into())),
        ("cofactor".into(), J::Int(c.h.into())),
        ("generator".into(), point_json(c.g)),
        (
            "target_group".into(),
            J::Str("prime-order-subgroup-generated-by-G".into()),
        ),
    ];
    let field = json::obj([
        ("p", J::Int(2)),
        ("n", J::Int(c.n.into())),
        ("basis", J::Str("polynomial".into())),
        ("modulus", J::Int(c.modulus.into())),
        (
            "element_encoding",
            J::Str("unsigned-integer-polynomial-bits".into()),
        ),
    ]);
    let draft = J::Obj(vec![
        ("field".into(), field.clone()),
        ("curve".into(), J::Obj(curve.clone())),
    ]);
    let id = format!("EC1N{}Ckb{}h{}", c.n, c.a, &sha256(&draft)?[..12]);
    curve.push(("curve_id".into(), J::Str(id)));
    Ok(J::Obj(vec![
        ("field".into(), field),
        ("curve".into(), J::Obj(curve)),
    ]))
}

fn base_points(report: &J, c: &Curve) -> Result<Vec<(u64, u64)>, String> {
    let (flat, orbits) = (report.get("factor_base"), report.get("factor_base_orbits"));
    require(
        flat.is_some() != orbits.is_some(),
        "supply exactly one factor-base representation",
    )?;
    let mut points: Vec<Option<(u64, u64)>> = Vec::new();
    if let Some(flat) = flat {
        for p in flat.as_arr().ok_or("malformed factor base")? {
            points.push(c.decode(p)?);
        }
    } else {
        for encoded in orbits
            .and_then(J::as_arr)
            .ok_or("malformed factor-base orbits")?
        {
            let start = c.decode(encoded)?;
            require(start.is_some(), "identity orbit representative")?;
            let mut p = start;
            for _ in 0..c.n {
                points.push(p);
                points.push(c.neg(p));
                p = c.frob(p);
            }
            require(p == start, "Frobenius orbit does not close")?;
        }
    }
    let distinct: BTreeSet<Option<(u64, u64)>> = points.iter().copied().collect();
    require(
        !points.is_empty() && !points.contains(&None) && distinct.len() == points.len(),
        "empty, duplicate or identity geometric base",
    )?;
    Ok(points.into_iter().map(|p| p.expect("checked")).collect())
}

/// One representative, the least point, of each sign-and-Frobenius orbit.
fn subgroup_orbits(
    points: &BTreeSet<(u64, u64)>,
    c: &Curve,
) -> Result<BTreeSet<(u64, u64)>, String> {
    let mut remaining = points.clone();
    let mut representatives = BTreeSet::new();
    while let Some(&first) = remaining.iter().next() {
        require(
            c.mul(Some(first), c.r as u128).is_none(),
            "non-subgroup factor-base point",
        )?;
        let mut orbit = BTreeSet::new();
        let mut p = Some(first);
        for _ in 0..c.n {
            orbit.insert(p.expect("on the curve"));
            orbit.insert(c.neg(p).expect("on the curve"));
            p = c.frob(p);
        }
        representatives.insert(*orbit.iter().next().expect("nonempty"));
        for q in &orbit {
            remaining.remove(q);
        }
    }
    Ok(representatives)
}

/// Distinct usable points before folding, with the raw geometry kept apart.
pub fn factor_base_inventory(report: &J, c: &Curve) -> Result<J, String> {
    let points = base_points(report, c)?;
    let convention = report
        .get("column_convention")
        .and_then(J::as_str)
        .unwrap_or("cofactor");
    require(
        convention == "representative" || convention == "cofactor",
        "unknown column convention",
    )?;
    let projected: Vec<Option<(u64, u64)>> = points
        .iter()
        .map(|&p| {
            if convention == "cofactor" {
                c.mul(Some(p), c.h as u128)
            } else {
                Some(p)
            }
        })
        .collect();
    let usable: BTreeSet<(u64, u64)> = projected.iter().flatten().copied().collect();
    require(!usable.is_empty(), "base has no usable points")?;
    let representatives = subgroup_orbits(&usable, c)?;
    let columns = natural(report.at("columns")?, "column count", true)?;
    require(
        columns == representatives.len() as i128,
        "claimed folded column count differs from actual orbits",
    )?;
    if let Some(logs) = report.get("column_logs") {
        let logs = logs.as_arr().ok_or("malformed column logs")?;
        let mut column_points = BTreeSet::new();
        for entry in logs {
            let p = c.decode(entry.at("point")?)?;
            require(p.is_some(), "identity column representative")?;
            column_points.insert(p.expect("checked"));
        }
        let column_reps = subgroup_orbits(&column_points, c)?;
        require(
            logs.len() as i128 == columns && column_reps == representatives,
            "column representatives do not cover the usable base exactly",
        )?;
    }
    let mut sorted_points = points.clone();
    sorted_points.sort_unstable();
    let usable_sorted: Vec<(u64, u64)> = usable.iter().copied().collect();
    let identities = projected.iter().filter(|p| p.is_none()).count();
    Ok(J::Obj(vec![
        ("geometric_point_count".into(), J::Int(points.len() as i128)),
        (
            "geometric_set_sha256".into(),
            J::Str(sha256(&points_json(&sorted_points))?),
        ),
        (
            "geometric_order_sha256".into(),
            J::Str(sha256(&points_json(&points))?),
        ),
        ("usable_point_count".into(), J::Int(usable.len() as i128)),
        (
            "usable_set_sha256".into(),
            J::Str(sha256(&points_json(&usable_sorted))?),
        ),
        (
            "subgroup_map".into(),
            J::Str(
                if convention == "cofactor" {
                    "cofactor-multiplication"
                } else {
                    "identity"
                }
                .into(),
            ),
        ),
        ("identity_images".into(), J::Int(identities as i128)),
        (
            "duplicate_nonidentity_images".into(),
            J::Int((projected.len() - identities - usable.len()) as i128),
        ),
        ("quotient".into(), J::Str("sign-and-Frobenius".into())),
        ("effective_columns".into(), J::Int(columns)),
    ]))
}

const STAGE_FIELDS: [(&str, &str); 4] = [
    (
        "point_decomposition",
        "summands solver summation_polynomial encoding equation_order monomial_order \
         internal_matrix_kernel limits cache_policy source_sha256",
    ),
    (
        "relation_collection",
        "collector query_distribution query_rule filtering verification duplicates \
         dependencies stop_rule source_sha256",
    ),
    (
        "relation_linear_algebra",
        "solver modulus matrix_construction orbit_quotient rank_criterion block_parameters \
         preconditioner source_sha256",
    ),
    (
        "target_descent",
        "method policy recursive_solvers success_rule stop_rule source_sha256",
    ),
];

fn resolved(v: &J) -> bool {
    match v {
        J::Null => false,
        J::Str(s) if s.is_empty() => false,
        J::Obj(kv) => kv.iter().all(|(_, v)| resolved(v)),
        J::Arr(items) => items.iter().all(resolved),
        _ => true,
    }
}

fn stage<'a>(method: &'a J, name: &str) -> Result<&'a J, String> {
    method.at(name)
}

/// An explicit method resolved against an independently counted inventory.
pub fn candidate_record(c: &Curve, fixture: &J, inventory: &J, method: &J) -> Result<J, String> {
    fields(
        method,
        "isogeny endomorphism factor_base point_decomposition relation_collection \
         relation_linear_algebra target_descent implementation",
        "method",
    )?;
    require(
        method.at("isogeny")?.as_str() == Some("none"),
        "ISO1 needs a verified route adapter",
    )?;
    let endo = method.at("endomorphism")?;
    fields(
        endo,
        "order_conductor frobenius_order_conductor volcano_levels",
        "endomorphism",
    )?;
    for key in ["order_conductor", "frobenius_order_conductor"] {
        let value = endo.at(key)?;
        if !matches!(value, J::Null) {
            fields(value, "value proof_sha256", key)?;
            natural(value.at("value")?, key, true)?;
            digest(value.at("proof_sha256")?, "conductor proof")?;
        }
    }
    let levels = endo
        .at("volcano_levels")?
        .as_arr()
        .ok_or("invalid volcano levels")?;
    for level in levels {
        fields(level, "ell level proof_sha256", "volcano level")?;
        natural(level.at("ell")?, "ell", true)?;
        natural(level.at("level")?, "volcano level", false)?;
        digest(level.at("proof_sha256")?, "volcano proof")?;
    }
    let fb = method.at("factor_base")?;
    fields(fb, "construction nominal_bound", "factor-base recipe")?;
    require(
        fb.at("construction")?
            .as_obj()
            .is_some_and(|kv| !kv.is_empty()),
        "missing exact factor-base construction",
    )?;
    for (name, names) in STAGE_FIELDS {
        let s = stage(method, name)?;
        fields(s, names, name)?;
        require(resolved(s), &format!("unresolved {name}"))?;
        digest(s.at("source_sha256")?, &format!("{name} source"))?;
    }
    let imp = method.at("implementation")?;
    fields(
        imp,
        "source_manifest_sha256 components flags",
        "implementation",
    )?;
    digest(
        imp.at("source_manifest_sha256")?,
        "implementation source manifest",
    )?;
    require(
        imp.at("flags")?.as_obj().is_some(),
        "invalid implementation flags",
    )?;
    let components = imp
        .at("components")?
        .as_arr()
        .filter(|v| !v.is_empty())
        .ok_or("missing executed components")?;
    let mut roles = BTreeSet::new();
    let mut sources = BTreeSet::new();
    for component in components {
        fields(component, "role sha256", "component")?;
        let role = component.at("role")?.as_str().unwrap_or("");
        require(
            !role.is_empty() && !roles.contains(role),
            "missing or duplicate component role",
        )?;
        digest(component.at("sha256")?, "component digest")?;
        roles.insert(role.to_string());
        sources.insert(component.at("sha256")?.as_str().unwrap_or("").to_string());
    }
    sources.insert(
        imp.at("source_manifest_sha256")?
            .as_str()
            .unwrap_or("")
            .to_string(),
    );
    let all_sourced = STAGE_FIELDS.iter().all(|(name, _)| {
        stage(method, name)
            .ok()
            .and_then(|s| s.get("source_sha256"))
            .and_then(J::as_str)
            .is_some_and(|h| sources.contains(h))
    });
    require(
        all_sourced,
        "stage source is absent from the implementation manifest/components",
    )?;
    let pdp = method.at("point_decomposition")?;
    natural(pdp.at("summands")?, "summand count", true)?;
    let modulus = method.at("relation_linear_algebra")?.at("modulus")?;
    let order = fixture.at("subgroup_order")?;
    let order = match order {
        J::Int(i) => *i,
        J::Str(s) => s.parse().map_err(|_| "subgroup_order".to_string())?,
        _ => return Err("subgroup_order".into()),
    };
    require(
        matches!(modulus, J::Int(m) if *m == order),
        "relation LA must use the subgroup scalar field",
    )?;
    let codes = [
        pdp.at("solver")?,
        method.at("relation_collection")?.at("collector")?,
        method.at("relation_linear_algebra")?.at("solver")?,
        method.at("target_descent")?.at("method")?,
    ];
    require(
        codes.iter().all(|v| v.as_str().is_some_and(compact_code)),
        "invalid compact stage code",
    )?;
    let mut method = method.clone();
    if let J::Obj(kv) = &mut method {
        if let Some((_, J::Obj(fb))) = kv.iter_mut().find(|(k, _)| k == "factor_base") {
            fb.push(("inventory".into(), inventory.clone()));
        }
    }
    identity_only(&method)?;
    let J::Obj(curve_kv) = curve_record(c)? else {
        unreachable!("a record is an object")
    };
    let J::Obj(method_kv) = method else {
        unreachable!("checked above")
    };
    let mut record = vec![("schema_version".to_string(), J::Int(1))];
    record.extend(curve_kv);
    record.extend(method_kv);
    let record = J::Obj(record);
    canonical(&record)?;
    Ok(record)
}

/// The IC1 candidate: its record, the record's digest and its label, with
/// the inventory counted from `report`'s factor base.
pub fn candidate_manifest(c: &Curve, fixture: &J, report: &J, method: &J) -> Result<J, String> {
    let inventory = &factor_base_inventory(report, c)?;
    let record = candidate_record(c, fixture, inventory, method)?;
    let pdp = record.at("point_decomposition")?;
    let s = |v: &J| v.as_str().unwrap_or("").to_string();
    let label = format!(
        "IC1N{}C{}fb{}PDP{}{}RC{}LA{}TD{}ISO0h{}",
        record.at("field")?.at("n")?.as_i128().unwrap_or(0),
        s(record.at("curve")?.at("tag")?),
        inventory.at("usable_point_count")?.as_i128().unwrap_or(0),
        pdp.at("summands")?.as_i128().unwrap_or(0),
        s(pdp.at("solver")?),
        s(record.at("relation_collection")?.at("collector")?),
        s(record.at("relation_linear_algebra")?.at("solver")?),
        s(record.at("target_descent")?.at("method")?),
        &sha256(&record)?[..12]
    );
    let h = sha256(&record)?;
    Ok(J::Obj(vec![
        ("candidate_id".into(), J::Str(label)),
        ("record_sha256".into(), J::Str(h)),
        ("record".into(), record),
    ]))
}

/// The workload: the curve, the targets and how they were drawn, the
/// seeds and the resources.
pub fn workload_manifest(
    c: &Curve,
    fixture: &J,
    input_law: &str,
    algorithm_seed: i128,
    envelope: &J,
    cache_policy: &str,
) -> Result<J, String> {
    require(
        cache_policy == "cold" || cache_policy == "warm",
        "unknown workload cache policy",
    )?;
    require(!input_law.is_empty(), "missing target input law")?;
    require(algorithm_seed >= 0, "invalid algorithm seed")?;
    let targets = fixture
        .at("targets")?
        .as_arr()
        .filter(|v| !v.is_empty())
        .ok_or("empty target workload")?;
    let seeds = fixture.at("target_seeds")?;
    require(
        seeds.as_arr().map(<[J]>::len) == Some(targets.len()),
        "missing target seeds",
    )?;
    require(
        fixture.at("target_scalar_constructed")? == &J::Bool(false),
        "planted target workload",
    )?;
    require(
        envelope.as_obj().is_some_and(|kv| !kv.is_empty()),
        "missing resources",
    )?;
    let mut encoded = Vec::new();
    for t in targets {
        let p = c.decode(t)?;
        require(
            p.is_some() && c.mul(p, c.r as u128).is_none(),
            "target outside subgroup",
        )?;
        encoded.push(point_json(p.expect("checked")));
    }
    let record = J::Obj(vec![
        (
            "curve_id".into(),
            curve_record(c)?.at("curve")?.at("curve_id")?.clone(),
        ),
        ("targets".into(), J::Arr(encoded)),
        ("input_law".into(), J::Str(input_law.into())),
        ("target_seeds".into(), seeds.clone()),
        ("algorithm_seed".into(), J::Int(algorithm_seed)),
        ("cache_policy".into(), J::Str(cache_policy.into())),
        ("target_count".into(), J::Int(targets.len() as i128)),
        ("resource_envelope".into(), envelope.clone()),
    ]);
    let h = sha256(&record)?;
    Ok(J::Obj(vec![
        ("workload_id".into(), J::Str(h[..12].to_string())),
        ("record_sha256".into(), J::Str(h)),
        ("record".into(), record),
    ]))
}

/// Whether `id` is an IC1 candidate label.
fn ic1_label(id: &str) -> bool {
    // IC1N<d>C<code>fb<d>PDP<d><code>RC<code>LA<code>TD<code>ISO[01]h<hex{12,64}>
    fn number(s: &str) -> Option<&str> {
        let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
        (end > 0 && !s.starts_with('0')).then(|| &s[end..])
    }
    fn code<'a>(s: &'a str, stop: &str) -> Option<&'a str> {
        let at = s.find(stop)?;
        compact_code(&s[..at]).then(|| &s[at..])
    }
    let go = || -> Option<()> {
        let s = number(id.strip_prefix("IC1N")?)?;
        let s = code(s.strip_prefix('C')?, "fb")?;
        let s = number(s.strip_prefix("fb")?)?;
        let s = number(s.strip_prefix("PDP")?)?;
        let s = code(s, "RC")?;
        let s = code(s.strip_prefix("RC")?, "LA")?;
        let s = code(s.strip_prefix("LA")?, "TD")?;
        let s = code(s.strip_prefix("TD")?, "ISO")?;
        let s = s.strip_prefix("ISO")?;
        let s = s.strip_prefix('0').or_else(|| s.strip_prefix('1'))?;
        let hex = s.strip_prefix('h')?;
        ((12..=64).contains(&hex.len())
            && hex
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b)))
        .then_some(())
    };
    go().is_some()
}

pub fn run_id(candidate_id: &str, workload_id: &str, number: i128) -> Result<String, String> {
    require(ic1_label(candidate_id), "invalid candidate ID")?;
    require(
        (12..=64).contains(&workload_id.len())
            && workload_id
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b)),
        "invalid workload ID",
    )?;
    require(number >= 0, "invalid run number")?;
    Ok(format!("{candidate_id}W{workload_id}R{number}"))
}

/// A record written once: a reused (truncated) identity must name the same
/// record byte for byte, and is never overwritten.
pub fn write_immutable(path: &Path, value: &J) -> Result<(), String> {
    let data = canonical(value)? + "\n";
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    if path.exists() {
        let old = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        return require(
            old == data.as_bytes(),
            "immutable record differs; extend a colliding digest",
        );
    }
    std::fs::write(path, data).map_err(|e| format!("{}: {e}", path.display()))
}

// ── EC1 (tools/curve_identity.py) ───────────────────────────────────

/// The EC1 digest: the same canonical JSON, with `curve_identity.py`'s
/// refusal.
fn ec1_digest(v: &J) -> Result<String, String> {
    no_floats(v).map_err(|_| "Identity records require JSON without floats.".to_string())?;
    sha256(v)
}

/// The EC1 representation identity of `field` and `curve` under `tag`.
pub fn ec1_curve_identity(field: &J, curve: &J, tag: &str) -> Result<J, String> {
    let (Some(f), Some(c)) = (field.as_obj(), curve.as_obj()) else {
        return Err("Field and curve records must be objects.".into());
    };
    fn get<'a>(kv: &'a [(String, J)], k: &str) -> Option<&'a J> {
        kv.iter().find(|(x, _)| x == k).map(|(_, v)| v)
    }
    let present = |kv: &[(String, J)], k: &str| get(kv, k).is_some_and(|v| !matches!(v, J::Null));
    let ok = [
        "characteristic",
        "degree",
        "representation",
        "element_encoding",
    ]
    .iter()
    .all(|k| present(f, k))
        && [
            "model",
            "subgroup_order",
            "cofactor",
            "generator",
            "target_group",
        ]
        .iter()
        .all(|k| present(c, k));
    require(
        ok,
        "Exact field, curve, subgroup and generator metadata are required.",
    )?;
    let (Some(J::Int(p)), Some(J::Int(n))) = (get(f, "characteristic"), get(f, "degree")) else {
        return Err("Characteristic and degree must be positive integers.".into());
    };
    let (p, n) = (*p, *n);
    require(
        p >= 2 && n >= 1,
        "Characteristic and degree must be positive integers.",
    )?;
    require(compact_code(tag), "Invalid readable curve tag.")?;
    require(
        n == 1
            || ["modulus_exponents", "modulus", "defining_polynomial"]
                .iter()
                .any(|k| present(f, k)),
        "Extension-field defining polynomial is required.",
    )?;
    let prefix = if p == 2 {
        format!("N{n}")
    } else if n == 1 {
        format!("P{}", 128 - p.leading_zeros())
    } else {
        format!("Q{p}D{n}")
    };
    let clean = J::Obj(c.iter().filter(|(k, _)| k != "curve_id").cloned().collect());
    let sha = ec1_digest(&J::Obj(vec![
        ("field".into(), field.clone()),
        ("curve".into(), clean),
    ]))?;
    Ok(J::Obj(vec![
        (
            "curve_id".into(),
            J::Str(format!("EC1{prefix}C{tag}h{}", &sha[..12])),
        ),
        ("curve_sha256".into(), J::Str(sha.clone())),
        (
            "curve_uid".into(),
            J::Str(format!("urn:ec-record:1:sha256:{sha}")),
        ),
        ("field_sha256".into(), J::Str(ec1_digest(field)?)),
    ]))
}

/// A recorded EC1 alias checked against its stored preimage.
pub fn ec1_inspect_manifest(manifest: &J) -> Result<J, String> {
    let (field, curve) = (manifest.at("field")?, manifest.at("curve")?);
    let alias = curve.get("curve_id").and_then(J::as_str).unwrap_or("");
    // EC1(N\d+|P\d+|Q\d+D\d+)C([a-z][a-z0-9]*)h([a-f0-9]{12,64})
    let parse = || -> Option<(String, String)> {
        let rest = alias.strip_prefix("EC1")?;
        let digits = |s: &str| s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
        let rest = match rest.chars().next()? {
            'N' | 'P' => {
                let d = digits(&rest[1..]);
                (d > 0).then(|| &rest[1 + d..])?
            }
            'Q' => {
                let d = digits(&rest[1..]);
                let r = rest[1 + d..].strip_prefix('D').filter(|_| d > 0)?;
                let e = digits(r);
                (e > 0).then(|| &r[e..])?
            }
            _ => return None,
        };
        let rest = rest.strip_prefix('C')?;
        let at = rest.rfind('h')?;
        let (tag, hex) = (&rest[..at], &rest[at + 1..]);
        (compact_code(tag)
            && (12..=64).contains(&hex.len())
            && hex
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b)))
        .then(|| (tag.to_string(), hex.to_string()))
    };
    let (tag, hex) = parse().ok_or("Missing or unsupported EC1 alias.")?;
    let identity = ec1_curve_identity(field, curve, &tag)?;
    let computed = identity.at("curve_id")?.as_str().unwrap_or("");
    let prefix_of = |s: &str| {
        s.rsplit_once('h')
            .map(|(p, _)| p.to_string())
            .unwrap_or_default()
    };
    require(
        prefix_of(alias) == prefix_of(computed)
            && identity
                .at("curve_sha256")?
                .as_str()
                .unwrap_or("")
                .starts_with(&hex),
        "EC1 alias does not match the canonical curve record.",
    )?;
    let J::Obj(mut kv) = identity else {
        unreachable!("an object")
    };
    if let Some((_, v)) = kv.iter_mut().find(|(k, _)| k == "curve_id") {
        *v = J::Str(alias.to_string());
    }
    Ok(J::Obj(kv))
}

/// The global candidate identity of a recorded method configuration.
pub fn ec1_candidate_identity(manifest: &J) -> Result<J, String> {
    ec1_inspect_manifest(manifest)?;
    let sha = ec1_digest(manifest)?;
    Ok(J::Obj(vec![
        (
            "candidate_uid".into(),
            J::Str(format!("urn:ec-candidate:1:sha256:{sha}")),
        ),
        ("candidate_sha256".into(), J::Str(sha)),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_parse_as_the_adapters_pattern() {
        assert!(ic1_label(
            "IC1N41Ckb0fb6560PDP3pairtableRCaimedLAsparseTDwalkISO0he13eaebe8e18"
        ));
        assert!(!ic1_label(
            "IC1N41Ckb0fb6560PDP3pairtableRCaimedLAsparseTDwalkISO2he13eaebe8e18"
        ));
        assert!(!ic1_label("IC1N041Ckb0fb1PDP3xRCaLAbTDcISO0h0123456789ab"));
        assert_eq!(
            run_id(
                "IC1N41Ckb0fb6560PDP3pairtableRCaimedLAsparseTDwalkISO0he13eaebe8e18",
                "0123456789ab",
                2
            )
            .unwrap(),
            "IC1N41Ckb0fb6560PDP3pairtableRCaimedLAsparseTDwalkISO0he13eaebe8e18W0123456789abR2"
        );
    }

    #[test]
    fn canonical_json_sorts_keys_and_refuses_floats() {
        let v = json::parse(r#"{"b": [1, "é"], "a": {"d": null, "c": true}}"#).unwrap();
        assert_eq!(
            canonical(&v).unwrap(),
            r#"{"a":{"c":true,"d":null},"b":[1,"é"]}"#
        );
        assert!(canonical(&json::parse(r#"{"x": 1.5}"#).unwrap()).is_err());
    }
}
