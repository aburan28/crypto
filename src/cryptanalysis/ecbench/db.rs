//! The database: SQL that loads committed evidence into SQLite.
//!
//! `ecbench db sql <paths…> | sqlite3 ecbench.db` builds or extends an
//! index over session directories, comparison files and factor-base
//! dumps.  The schema is `docs/ecbench/schema.sql`, compiled in, so the
//! file reviewers read is the one that runs.  No database driver is
//! linked: the output is plain SQL with every literal escaped here.
//!
//! Every statement is idempotent, and only in the one way intended: an
//! insert tolerates a conflict on its table's *primary key* when the row
//! is the same row again, and nothing else.  A second row claiming an
//! existing run id is a UNIQUE error; a key arriving with a different
//! identity (a curve slug on another ICV1 string, a workload, method,
//! factor-base or host id on another hash) is stopped by the schema's
//! identity triggers; a failed CHECK is an error.  None is a silently
//! dropped row.  A session row is upserted, so a session loaded while
//! running and again once complete ends up complete.  Load with
//! `sqlite3 -bail` to stop at the first error.

use std::fmt::Write as _;
use std::path::Path;

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use serde_json::Value;

use crate::cryptanalysis::ecbench::compare::Comparison;
use crate::cryptanalysis::ecbench::host::{is_virtual, HostCapsule};
use crate::cryptanalysis::ecbench::methods::FactorBaseDump;
use crate::cryptanalysis::ecbench::record::{floor_s, Record};
use crate::cryptanalysis::ecbench::runner::{read_records, read_session, PlanDoc, Session};
use crate::cryptanalysis::p256_dickson_factor_base::WideFactorBaseDump;

pub const SCHEMA_SQL: &str = include_str!("../../../docs/ecbench/schema.sql");

/// A SQL text literal, or NULL.
fn t(s: Option<&str>) -> String {
    match s {
        Some(s) => format!("'{}'", s.replace('\'', "''")),
        None => "NULL".into(),
    }
}

fn ts(s: &str) -> String {
    t(Some(s))
}

/// An integer literal, or NULL.
fn i<T: Into<i128>>(v: Option<T>) -> String {
    v.map(|x| x.into().to_string())
        .unwrap_or_else(|| "NULL".into())
}

/// A real literal; NULL for non-finite.
fn f(v: Option<f64>) -> String {
    match v {
        Some(x) if x.is_finite() => format!("{x:e}"),
        _ => "NULL".into(),
    }
}

fn b(v: bool) -> &'static str {
    if v {
        "1"
    } else {
        "0"
    }
}

fn json_text(v: &impl serde::Serialize) -> String {
    ts(&serde_json::to_string(v).unwrap_or_else(|_| "null".into()))
}

fn curve_rows(out: &mut String, w: &crate::cryptanalysis::ecbench::workload::Workload) {
    let c = &w.curve;
    let _ = writeln!(
        out,
        "INSERT INTO curves VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (slug) DO NOTHING;",
        ts(&c.slug),
        ts(&c.icv1),
        ts(&c.family),
        i(c.field_degree),
        i(c.field_bits),
        ts(&c.group_order.to_string()),
        ts(&c.r.to_string()),
        ts(&c.cofactor.to_string()),
        f(Some((c.r as f64).log2())),
        c.automorphisms_available,
        f(Some(floor_s(c.automorphisms_available))),
        b(c.registered),
    );
    let _ = writeln!(
        out,
        "INSERT INTO curve_representations VALUES ({}, {}, {}, {}, {}) ON CONFLICT (slug, generator_x, generator_y) DO NOTHING;",
        ts(&c.slug),
        ts(&c.generator[0]),
        ts(&c.generator[1]),
        t(c.ec1.as_deref()),
        t(c.curve_uid.as_deref()),
    );
    let _ = writeln!(
        out,
        "INSERT INTO curve_constructions VALUES ({}, {}) ON CONFLICT (slug, construction) DO NOTHING;",
        ts(&c.slug),
        ts(&c.construction),
    );
    let _ = writeln!(
        out,
        "INSERT INTO workloads VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (workload_id) DO NOTHING;",
        ts(&w.workload_id),
        ts(&w.workload_sha256),
        ts(&c.slug),
        ts(&c.generator[0]),
        ts(&c.generator[1]),
        ts(&w.target[0]),
        ts(&w.target[1]),
        ts(&w.target_law),
        ts(&w.target_seed.to_string()),
        w.target_index,
        t(w.planted.map(|k| k.to_string()).as_deref()),
    );
}

fn host_row(out: &mut String, h: &HostCapsule) {
    let s = &h.stable;
    let _ = writeln!(
        out,
        "INSERT INTO hosts VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (env_class_id) DO NOTHING;",
        ts(&h.env_class_id),
        ts(&h.env_class_sha256),
        ts(&s.os),
        ts(&s.arch),
        t(s.kernel_release.as_deref()),
        t(s.cpu_model.as_deref()),
        ts(&s.features.join(" ")),
        s.logical_cpus,
        i(s.physical_cores),
        s.numa_nodes,
        i(is_virtual(s).map(i32::from)),
        t(s.rustc_vv.as_deref().and_then(|v| v.lines().next())),
        json_text(h),
    );
}

fn session_rows(out: &mut String, s: &Session, dir: &Path) {
    let plan = s.cpu_plan.as_ref();
    let _ = writeln!(
        out,
        "INSERT INTO sessions VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (session_id) DO UPDATE SET finished_unix_ms = excluded.finished_unix_ms, status = excluded.status, records_sha256 = excluded.records_sha256, path = excluded.path;",
        ts(&s.session_id),
        ts(&s.label),
        ts(&s.spec_id),
        ts(&s.spec_sha256),
        ts(&s.env_class_id),
        t(s.git_commit.as_deref()),
        i(s.git_dirty.map(i32::from)),
        t(s.binary_sha256.as_deref()),
        json_text(&s.options.cpus),
        i(plan.map(|p| p.run_cpu)),
        t(plan
            .map(|p| p.reserved.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(","))
            .as_deref()),
        i(plan.and_then(|p| p.node)),
        i(s.preflight.as_ref().map(|p| i32::from(p.quiet))),
        i(Some(s.started_unix_ms as i128)),
        i(s.finished_unix_ms.map(|v| v as i128)),
        ts(&s.status),
        t(s.records_sha256.as_deref()),
        ts(&dir.display().to_string()),
    );
}

fn run_rows(out: &mut String, r: &Record) {
    let m = &r.method;
    let _ = writeln!(
        out,
        "INSERT INTO algorithms VALUES ({}, {}, {}, {}, {}, {}) ON CONFLICT (method_id) DO NOTHING;",
        ts(&m.method_id),
        ts(&m.method_sha256),
        ts(&m.id),
        ts(&m.family),
        json_text(&m.params),
        ts(&m.entry),
    );
    curve_rows(out, &r.workload);
    if let Some(fb) = &r.factor_base {
        let _ = writeln!(
            out,
            "INSERT INTO factor_bases VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (fb_id) DO NOTHING;",
            ts(&fb.fb_id),
            ts(&fb.fb_sha256),
            ts(&r.workload.curve.slug),
            ts(&fb.family),
            json_text(&fb.params),
            ts(&fb.description),
            fb.signed_points,
            fb.abscissae,
            fb.columns,
            i(fb.dimension),
            ts(&fb.points_sha256),
        );
    }
    let _ = writeln!(
        out,
        "INSERT INTO arms VALUES ({}, {}, {}, {}) ON CONFLICT (session_id, arm) DO NOTHING;",
        ts(&r.session_id),
        ts(&r.arm),
        ts(&r.role),
        ts(&m.method_id),
    );
    let tm = &r.time;
    let _ = writeln!(
        out,
        "INSERT INTO runs VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (record_id) DO NOTHING;",
        ts(&r.record_id),
        ts(&r.session_id),
        ts(&r.run_id),
        r.seq,
        r.round,
        b(r.warmup),
        ts(&r.arm),
        ts(&m.method_id),
        ts(&r.workload.workload_id),
        t(r.factor_base.as_ref().map(|fb| fb.fb_id.as_str())),
        ts(&r.algorithm_seed),
        ts(&r.outcome.status),
        t(r.outcome.recovered.as_deref()),
        f(r.cost.total_gae),
        f(r.cost.s),
        f(r.boundaries.ratio_to_floor),
        i(r.boundaries.automorphisms_used),
        b(r.cost.lower_bound),
        ts(&r.cost.unpriced.join(" ")),
        b(r.cost.deterministic),
        i(tm.solve_wall_ns.map(|v| v as i128)),
        tm.process_wall_ns,
        tm.user_ns,
        tm.sys_ns,
        tm.max_rss_kib,
        tm.involuntary_switches,
        i(tm.schedstat_solve.map(|s| s.run_delay_ns as i128)),
        ts(r.isolation.level.name()),
        i(r.isolation.foreign_busy_ticks),
        i(r.isolation.steal_ticks.map(|v| v as i128)),
        t(r.outcome.error.as_deref()),
    );
    for (ord, p) in r.phases.iter().enumerate() {
        let _ = writeln!(
            out,
            "INSERT INTO phases VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (record_id, ord) DO NOTHING;",
            ts(&r.record_id),
            ord,
            ts(&p.name),
            p.adds,
            p.doubles,
            p.scalar_mults,
            f(Some(p.gae)),
            i(p.wall_ns.map(|v| v as i128)),
            json_text(&p.native),
        );
    }
    for (k, v) in &r.counters {
        let _ = writeln!(
            out,
            "INSERT INTO run_counters VALUES ({}, {}, {}) ON CONFLICT (record_id, name) DO NOTHING;",
            ts(&r.record_id),
            ts(k),
            v
        );
    }
    for (ord, bl) in r.isolation.blockers.iter().enumerate() {
        let _ = writeln!(
            out,
            "INSERT INTO isolation_blockers VALUES ({}, {}, {}) ON CONFLICT (record_id, ord) DO NOTHING;",
            ts(&r.record_id),
            ord,
            ts(bl)
        );
    }
}

fn comparison_row(out: &mut String, c: &Comparison) {
    let _ = writeln!(
        out,
        "INSERT OR REPLACE INTO comparisons VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});",
        ts(&c.comparison_id),
        ts(&c.a.session_id),
        ts(&c.a.arm),
        ts(&c.b.session_id),
        ts(&c.b.arm),
        ts(&c.ops.status),
        f(c.ops.ratio_b_over_a),
        f(c.ops.ci95.map(|x| x.0)),
        f(c.ops.ci95.map(|x| x.1)),
        b(c.ops.bounded),
        ts(&c.wall.status),
        f(c.wall.median_ratio),
        f(c.wall.ci95.map(|x| x.0)),
        f(c.wall.ci95.map(|x| x.1)),
        ts(&c.verdict),
        json_text(c),
    );
}

/// A `vs_rho` claim (`ecbench claim build`), checked again at load time.
fn claim_row(out: &mut String, c: &Value) -> Result<(), String> {
    let s = |k: &str| c.get(k).and_then(Value::as_str);
    let need = |k: &str| s(k).ok_or_else(|| format!("claim without {k}"));
    let num = |k: &str| {
        c.get(k)
            .and_then(Value::as_f64)
            .ok_or_else(|| format!("claim without {k}"))
    };
    let check = crate::cryptanalysis::ecbench::claim::check_vs_rho(c)?;
    let _ = writeln!(
        out,
        "INSERT OR REPLACE INTO claims VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {});",
        ts(need("run_id")?),
        ts(need("candidate_id")?),
        ts(need("candidate_manifest_sha256")?),
        ts(need("workload_id")?),
        ts(need("workload_manifest_sha256")?),
        ts(c["ecbench_records"]["session"]
            .as_str()
            .ok_or("claim without its session")?),
        ts(c["ecbench_records"]["ic"]
            .as_str()
            .ok_or("claim without its IC record")?),
        ts(c["ecbench_records"]["rho"]
            .as_str()
            .ok_or("claim without its rho record")?),
        i(c.get("n_or_bits").and_then(Value::as_u64)),
        f(Some(num("ic_online_wall_ms")?)),
        f(Some(num("rho_online_wall_ms")?)),
        f(Some(num("online_speedup")?)),
        b(c.get("independent_validation") == Some(&Value::Bool(true))),
        t(s("ic_replay_certificate_sha256")),
        ts(check["status"].as_str().unwrap_or("FAIL")),
        ts(need("verdict")?),
        json_text(c),
    );
    Ok(())
}

fn dump_rows(out: &mut String, d: &FactorBaseDump) {
    let fb = &d.factor_base;
    let c = &d.curve;
    let _ = writeln!(
        out,
        "INSERT INTO curves VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (slug) DO NOTHING;",
        ts(&c.slug),
        ts(&c.icv1),
        ts(&c.family),
        i(c.field_degree),
        i(c.field_bits),
        ts(&c.group_order.to_string()),
        ts(&c.r.to_string()),
        ts(&c.cofactor.to_string()),
        f(Some((c.r as f64).log2())),
        c.automorphisms_available,
        f(Some(floor_s(c.automorphisms_available))),
        b(c.registered),
    );
    let _ = writeln!(
        out,
        "INSERT INTO factor_bases VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (fb_id) DO NOTHING;",
        ts(&fb.fb_id),
        ts(&fb.fb_sha256),
        ts(&c.slug),
        ts(&fb.family),
        json_text(&fb.params),
        ts(&fb.description),
        fb.signed_points,
        fb.abscissae,
        fb.columns,
        i(fb.dimension),
        ts(&fb.points_sha256),
    );
    for (idx, p) in d.points.iter().enumerate() {
        let _ = writeln!(
            out,
            "INSERT INTO factor_base_points VALUES ({}, {}, {}, {}, {}, {}) ON CONFLICT (fb_id, idx) DO NOTHING;",
            ts(&fb.fb_id),
            idx,
            ts(&p.x),
            ts(&p.y),
            p.col,
            ts(&p.coef.to_string()),
        );
    }
}

fn wide_dump_rows(out: &mut String, d: &WideFactorBaseDump) -> Result<(), String> {
    let fb = &d.factor_base;
    let c = &d.curve;
    let r = BigUint::parse_bytes(c.r.as_bytes(), 10)
        .and_then(|n| n.to_f64())
        .ok_or_else(|| format!("wide subgroup order is not decimal: {}", c.r))?;
    let _ = writeln!(
        out,
        "INSERT INTO curves VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (slug) DO NOTHING;",
        ts(&c.slug),
        ts(&c.icv1),
        ts(&c.family),
        i(c.field_degree),
        i(Some(c.field_bits)),
        ts(&c.group_order),
        ts(&c.r),
        ts(&c.cofactor),
        f(Some(r.log2())),
        c.automorphisms_available,
        f(Some(floor_s(c.automorphisms_available))),
        b(c.registered),
    );
    let _ = writeln!(
        out,
        "INSERT INTO factor_bases VALUES ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}) ON CONFLICT (fb_id) DO NOTHING;",
        ts(&fb.fb_id),
        ts(&fb.fb_sha256),
        ts(&c.slug),
        ts(&fb.family),
        json_text(&fb.params),
        ts(&fb.description),
        fb.signed_points,
        fb.abscissae,
        fb.columns,
        i(fb.dimension),
        ts(&fb.points_sha256),
    );
    for (idx, p) in d.points.iter().enumerate() {
        let _ = writeln!(
            out,
            "INSERT INTO factor_base_points VALUES ({}, {}, {}, {}, {}, {}) ON CONFLICT (fb_id, idx) DO NOTHING;",
            ts(&fb.fb_id),
            idx,
            ts(&p.x),
            ts(&p.y),
            p.col,
            ts(&p.coef),
        );
    }
    Ok(())
}

/// SQL for one path: a session directory, a comparison file or a
/// factor-base dump, recognised by content.
pub fn sql_for_path(path: &Path) -> Result<String, String> {
    let mut out = String::new();
    if path.is_dir() {
        let session = read_session(path)?;
        let host: HostCapsule = serde_json::from_str(
            &std::fs::read_to_string(path.join("host.json")).map_err(|e| e.to_string())?,
        )
        .map_err(|e| format!("host.json: {e}"))?;
        host_row(&mut out, &host);
        session_rows(&mut out, &session, path);
        // Arms that never produced a record still belong to the session.
        if let Ok(text) = std::fs::read_to_string(path.join("plan.json")) {
            if let Ok(plan) = serde_json::from_str::<PlanDoc>(&text) {
                for w in &plan.workloads {
                    curve_rows(&mut out, w);
                }
            }
        }
        for r in read_records(path)? {
            run_rows(&mut out, &r);
        }
        if let Ok(rd) = std::fs::read_dir(path.join("comparisons")) {
            let mut files: Vec<_> = rd.filter_map(|e| e.ok()).map(|e| e.path()).collect();
            files.sort();
            for f in files {
                out.push_str(&sql_for_path(&f)?);
            }
        }
        return Ok(out);
    }
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let v: Value = serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    match v.get("schema").and_then(Value::as_str) {
        Some(crate::cryptanalysis::ecbench::compare::COMPARISON_SCHEMA) => {
            let c: Comparison = serde_json::from_value(v).map_err(|e| e.to_string())?;
            comparison_row(&mut out, &c);
        }
        Some(crate::cryptanalysis::ecbench::methods::FB_DUMP_SCHEMA) => {
            let d: FactorBaseDump = serde_json::from_value(v).map_err(|e| e.to_string())?;
            dump_rows(&mut out, &d);
        }
        Some(crate::cryptanalysis::p256_dickson_factor_base::DUMP_SCHEMA) => {
            let d: WideFactorBaseDump = serde_json::from_value(v).map_err(|e| e.to_string())?;
            wide_dump_rows(&mut out, &d)?;
        }
        Some(crate::cryptanalysis::ecbench::claim::CLAIM_SCHEMA) => {
            claim_row(&mut out, &v).map_err(|e| format!("{}: {e}", path.display()))?;
        }
        other => {
            return Err(format!(
                "{}: not an ecbench file (schema {other:?})",
                path.display()
            ))
        }
    }
    Ok(out)
}

/// The schema, then every path's rows, in one transaction.
pub fn sql(paths: &[&Path]) -> Result<String, String> {
    let mut out = String::from(SCHEMA_SQL);
    out.push_str("\nBEGIN;\n");
    for p in paths {
        out.push_str(&sql_for_path(p)?);
    }
    out.push_str("COMMIT;\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn literals_escape() {
        assert_eq!(ts("it's"), "'it''s'");
        assert_eq!(f(Some(f64::NAN)), "NULL");
        assert_eq!(i::<i64>(None), "NULL");
    }

    #[test]
    fn wide_factor_base_dump_loads() {
        let d = crate::cryptanalysis::p256_dickson_factor_base::build("dickson-torus:depth=3")
            .unwrap()
            .dump;
        let mut sql = String::new();
        wide_dump_rows(&mut sql, &d).unwrap();
        assert!(sql.contains(&d.factor_base.fb_id));
        assert!(sql.contains(&d.curve.slug));
        assert_eq!(
            sql.matches("INSERT INTO factor_base_points").count(),
            d.points.len()
        );
    }
}
