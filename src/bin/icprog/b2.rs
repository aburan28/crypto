//! B2's own measurements, natively (its protocol's measurements 6 and 7,
//! which had no harness yet):
//! - `estimate` (measurement 6): B2's F2 estimate of the index calculus's
//!   cold cost at each suite size, on `M1-T01`'s v2 translation by B2's arm
//!   with `method.fidelity` set to `F2`. F2 runs nothing, so each process
//!   is untimed.
//! - `stepcost` (measurement 7): `rho-bignum`'s one-thread step cost at
//!   three field widths (binary `n = 100` and `127`, prime `p ≈ 2^127`),
//!   each one isolated run of `examples/rho_bignum_rate.rs` built from B2's
//!   arm.
//! - `analyse`: the estimate's error against v0's measured cold times
//!   (R01), as a table of ratios with their geometric mean and range, and
//!   the step costs beside the provisional ones they replace.  Both are
//!   reported, not gated.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use super::bench::{self, Bench};
use super::bround;
use super::json::{self, obj, J};
use super::rounds::speed;
use super::runs;
use super::stats;
use super::suite::{self, Row};

/// The widths measurement 7 declares, and the steps of each run.
const WIDTHS: [(&str, u32); 3] = [("binary", 100), ("binary", 127), ("prime", 127)];
pub const STEPS: u64 = 20_000;

/// `M1-T01` at every suite size.
fn first_rows(r: &bround::Run) -> Result<Vec<Row>, String> {
    Ok(speed::m1_rows(&r.ctx)?
        .into_iter()
        .filter(|row| row.id.ends_with("/M1-T01"))
        .collect())
}

/// The row's v2 translation at `fidelity: F2`, written once.
fn f2_document(r: &bround::Run, cand: &Path, row: &Row) -> Result<PathBuf, String> {
    let path = r
        .ctx
        .runs
        .join("estimate")
        .join("docs")
        .join(format!("{}.json", row.id));
    if !path.exists() {
        let J::Obj(mut kv) = bround::translate_row(r, cand, row)? else {
            return Err(format!("{}: the translation is not one document", row.id));
        };
        for (k, v) in kv.iter_mut() {
            if let (true, J::Obj(m)) = (k == "method", v) {
                match m.iter_mut().find(|(mk, _)| mk == "fidelity") {
                    Some((_, f)) => *f = J::Str("F2".into()),
                    None => m.push(("fidelity".into(), J::Str("F2".into()))),
                }
            }
        }
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
        }
        std::fs::write(&path, json::dumps(&J::Obj(kv), 1) + "\n")
            .map_err(|e| format!("{}: {e}", path.display()))?;
    }
    Ok(path)
}

fn report_path(r: &bround::Run, row: &Row) -> PathBuf {
    r.ctx
        .runs
        .join("estimate")
        .join(format!("{}.report.json", row.id))
}

pub fn estimate(r: &bround::Run, cand: &Path) -> Result<J, String> {
    for row in first_rows(r)? {
        let doc = f2_document(r, cand, &row)?;
        let out = report_path(r, &row);
        if !out.exists() {
            Command::new("taskset")
                .args(["-c", bench::CPUS])
                .arg(cand)
                .arg("price")
                .arg("--params")
                .arg(&doc)
                .args(["--json", "--out"])
                .arg(&out)
                .env("RAYON_NUM_THREADS", "1")
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .status()
                .map_err(|e| format!("taskset: {e}"))?;
        }
        let rep = runs::load(&out);
        println!(
            "{}: {} {}",
            row.id,
            rep.get("status").and_then(J::as_str).unwrap_or("None"),
            rep.get("estimate")
                .and_then(|e| e.get("ic"))
                .and_then(|e| e.get("seconds"))
                .map_or("None".into(), |s| json::dumps_line(s, false))
        );
    }
    Ok(J::Null)
}

fn stepcost_path(r: &bround::Run, kind: &str, width: u32) -> PathBuf {
    r.ctx
        .runs
        .join("stepcost")
        .join(format!("{kind}-{width}.json"))
}

/// Each width once, isolated, with the runner's retries: the first clean
/// run is the width's figure.
pub fn stepcost(r: &bround::Run, b: &Bench, example: &Path) -> Result<J, String> {
    let d = r.ctx.runs.join("stepcost");
    for (kind, width) in WIDTHS {
        let out = stepcost_path(r, kind, width);
        for k in 0..=runs::RETRIES {
            let attempt = runs::attempt(&out, k);
            if !attempt.exists() {
                if let Some(parent) = attempt.parent() {
                    std::fs::create_dir_all(parent)
                        .map_err(|e| format!("{}: {e}", parent.display()))?;
                }
                let cmd = vec![
                    example.to_string_lossy().into_owned(),
                    kind.to_string(),
                    width.to_string(),
                    STEPS.to_string(),
                ];
                b.launch_to(&cmd, &attempt, &d, true)?;
            }
            if runs::clean(&attempt)? {
                break;
            }
        }
        println!(
            "{kind} {width}: {}",
            std::fs::read_to_string(&out).unwrap_or_default().trim()
        );
    }
    Ok(J::Null)
}

/// v0's measured cold milliseconds by size, from R01's profile.
fn v0_cold_ms(programme: &Path) -> Result<Vec<((u8, u32), f64)>, String> {
    let doc = json::read(&programme.join("rounds/R01-baseline-v0/analysis.json"))?;
    let mut out = Vec::new();
    for p in doc
        .at("profile")?
        .as_arr()
        .ok_or("`profile` is not a list")?
    {
        let size = p.at("size")?.as_str().ok_or("`size`")?;
        let a: u8 = size[1..2].parse().map_err(|_| format!("size `{size}`"))?;
        let n: u32 = size[3..].parse().map_err(|_| format!("size `{size}`"))?;
        let ms = p
            .at("cold_ms_median")?
            .as_f64()
            .ok_or("`cold_ms_median` is not a number")?;
        out.push(((a, n), ms));
    }
    Ok(out)
}

pub fn analyse(r: &bround::Run) -> Result<J, String> {
    let v0 = v0_cold_ms(&r.ctx.programme)?;
    let mut sizes = Vec::new();
    let mut ratios = Vec::new();
    let mut rows = first_rows(r)?;
    rows.sort_by(|x, y| x.r.partial_cmp(&y.r).expect("no NaN"));
    for row in rows {
        let rep = runs::load(&report_path(r, &row));
        let seconds = rep
            .get("estimate")
            .and_then(|e| e.get("ic"))
            .and_then(|e| e.get("seconds"))
            .and_then(J::as_f64);
        let cold = v0
            .iter()
            .find(|(s, _)| *s == (row.a, row.n))
            .map(|(_, ms)| *ms);
        let ratio = match (seconds, cold) {
            (Some(s), Some(ms)) if ms > 0.0 => Some(s * 1e3 / ms),
            _ => None,
        };
        if let Some(q) = ratio {
            ratios.push(q);
        }
        sizes.push(obj([
            ("slug", J::Str(suite::curve_slug(row.a, row.n)?)),
            ("n", J::Int(row.n.into())),
            ("status", rep.get("status").cloned().unwrap_or(J::Null)),
            (
                "estimate_ms",
                seconds.map_or(J::Null, |s| J::Float(s * 1e3)),
            ),
            ("v0_cold_ms", cold.map_or(J::Null, J::Float)),
            ("ratio", ratio.map_or(J::Null, J::Float)),
        ]));
    }
    let ci = stats::geo_ci(&ratios, 0.95);
    let error = obj([
        (
            "what",
            J::Str(
                "B2's F2 estimate of the index calculus's cold cost at each suite size over v0's \
                 measured cold time there (R01): a model's error, recorded, not gated"
                    .into(),
            ),
        ),
        ("sizes", J::Arr(sizes)),
        ("geomean", ci.geomean.map_or(J::Null, J::Float)),
        ("min", ci.min.map_or(J::Null, J::Float)),
        ("max", ci.max.map_or(J::Null, J::Float)),
    ]);
    let provisional = json::read(
        &r.ctx
            .programme
            .join("rounds/B2-fields-forms-estimates/estimates.json"),
    )?;
    let mut widths = Vec::new();
    for (kind, width) in WIDTHS {
        let out = stepcost_path(r, kind, width);
        let line = std::fs::read_to_string(&out)
            .ok()
            .and_then(|t| json::parse(t.trim()).ok())
            .unwrap_or(J::Null);
        let record = runs::run_record(&out)?.unwrap_or(J::Null);
        let before = provisional
            .get("rho-bignum")
            .and_then(|p| p.get(kind))
            .and_then(|k| k.get("points"))
            .and_then(J::as_arr)
            .and_then(|pts| {
                pts.iter()
                    .find(|p| p.get("width").and_then(J::as_i128) == Some(width.into()))
            })
            .and_then(|p| p.get("ns_per_step").cloned())
            .unwrap_or(J::Null);
        widths.push(obj([
            ("kind", J::Str(kind.into())),
            ("width", J::Int(width.into())),
            (
                "ns_per_step",
                line.get("ns_per_step").cloned().unwrap_or(J::Null),
            ),
            ("steps", line.get("steps").cloned().unwrap_or(J::Null)),
            (
                "contended",
                record.get("contended").cloned().unwrap_or(J::Null),
            ),
            ("provisional_ns_per_step", before),
        ]));
    }
    Ok(obj([
        ("estimate_error", error),
        (
            "rho_bignum_step_cost",
            obj([
                (
                    "what",
                    J::Str(
                        "one-thread walk steps, isolated, one run a width: the step costs \
                         rho-bignum's estimates use from B2's acceptance on"
                            .into(),
                    ),
                ),
                ("widths", J::Arr(widths)),
            ]),
        ),
    ]))
}
