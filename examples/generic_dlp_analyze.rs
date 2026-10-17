//! Analysis for `research/generic_dlp_memory_parallelism_20261017`
//! (PROTOCOL.md "Analysis"): reads a run's `raw/` directory and writes
//! `summary.json` beside it, plus a Markdown table on stdout.
//!
//! ```sh
//! cargo run --release --example generic_dlp_analyze -- <run_dir>
//! ```
//!
//! For every configuration (part, bits, kernel) it reports, per arm (A, A2,
//! B), the median and minimum of each timer; the ratio of medians A/B and
//! A/A2; the minimum and maximum of the per-round paired ratios; and the
//! verdict of the protocol's rule: **improves** when every per-round A/B
//! ratio lies above the largest per-round A/A2 ratio, **regresses** when
//! every one lies below the smallest, otherwise **no measurable change**.
//! It checks the pinned outputs (every counted field identical across the
//! three arms at one thread; the answer at three) and counts contended runs
//! from the isolation records.  PARI/GP rows are joined by curve and size.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde_json::{json, Map, Value};

fn read_jsonl(path: &Path) -> Vec<Value> {
    let Ok(text) = std::fs::read_to_string(path) else {
        return Vec::new();
    };
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("{}: {e}", path.display())))
        .collect()
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n == 0 {
        f64::NAN
    } else if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

fn min(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::INFINITY, f64::min)
}

fn max(v: &[f64]) -> f64 {
    v.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
}

/// `part-bits[-kernel]-arm-rROUND` split into its fields.
struct Label {
    part: String,
    bits: u32,
    kernel: String,
    arm: String,
    round: u32,
}

fn parse_label(s: &str) -> Label {
    let f: Vec<&str> = s.split('-').collect();
    let round = f.last().unwrap().trim_start_matches('r').parse().unwrap();
    let arm = f[f.len() - 2].to_string();
    let kernel = if f.len() == 5 { f[2].to_string() } else { "-".into() };
    Label {
        part: f[0].into(),
        bits: f[1].parse().unwrap(),
        kernel,
        arm,
        round,
    }
}

const BSGS_PINNED: [&str; 6] = [
    "answers_digest",
    "table_entries",
    "baby_steps",
    "giant_steps",
    "candidates",
    "false_candidates",
];
const RHO_PINNED: [&str; 4] = ["log", "verified", "total_steps", "distinguished_points"];

/// The secret behind a one-target `answers_digest` of
/// `lit2026_w11_bsgs_time`: the digest is one FNV-1a step
/// `(0xcbf29ce484222325 ^ k) · 0x100000001b3 mod 2^64`, and the prime is
/// odd, so it inverts.
fn secret_of_digest(digest: &str) -> u64 {
    let d = u64::from_str_radix(digest, 16).expect("hex digest");
    let prime: u64 = 0x100_0000_01b3;
    // Inverse of an odd number mod 2^64 by Newton iteration.
    let mut inv = prime;
    for _ in 0..6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(prime.wrapping_mul(inv)));
    }
    assert_eq!(prime.wrapping_mul(inv), 1);
    d.wrapping_mul(inv) ^ 0xcbf2_9ce4_8422_2325
}

/// `targets <run_dir>`: one `role bits k` line per curve, for the PARI/GP
/// arm, from the recorded one-thread samples.
fn print_targets(raw: &Path) {
    let mut seen = BTreeMap::new();
    for v in read_jsonl(&raw.join("bsgs1.jsonl")) {
        assert_eq!(v["targets"], json!(1), "one-target samples only");
        let k = secret_of_digest(v["answers_digest"].as_str().unwrap());
        seen.entry(("bsgs", v["bits"].as_u64().unwrap())).or_insert(k);
    }
    for v in read_jsonl(&raw.join("rho1.jsonl")) {
        seen.entry(("rho", v["bits"].as_u64().unwrap()))
            .or_insert(v["k"].as_u64().unwrap());
    }
    for ((role, bits), k) in seen {
        println!("{role} {bits} {k}");
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("targets") {
        print_targets(&PathBuf::from(&args[1]).join("raw"));
        return;
    }
    let run: PathBuf = args.first().expect("run directory").into();
    let raw = run.join("raw");
    let iso: BTreeMap<String, Value> = read_jsonl(&raw.join("isolation.jsonl"))
        .into_iter()
        .map(|r| (r["label"].as_str().unwrap_or("").to_string(), r))
        .collect();

    let mut configs: BTreeMap<(String, u32, String), Vec<(Label, Value)>> = BTreeMap::new();
    for part in ["bsgs1", "bsgs3", "rho1", "rho3"] {
        for v in read_jsonl(&raw.join(format!("{part}.jsonl"))) {
            let l = parse_label(v["label"].as_str().unwrap());
            configs
                .entry((l.part.clone(), l.bits, l.kernel.clone()))
                .or_default()
                .push((l, v));
        }
    }

    let mut rows = Vec::new();
    println!("| part | bits | kernel | timer | A median | A2 median | B median | A/B (ratio of medians) | A/B per round [min, max] | A/A2 per round [min, max] | verdict | pinned | contended (A, A2, B) |");
    println!("|---|---:|---|---|---:|---:|---:|---:|---|---|---|---|---|");
    for ((part, bits, kernel), samples) in &configs {
        let is_bsgs = part.starts_with("bsgs");
        let one_thread = part.ends_with('1');
        let timers: &[&str] = if is_bsgs {
            &["total_s", "build_s", "search_s"]
        } else {
            &["solve_s", "setup_s"]
        };
        // Pinned outputs: every counted field identical across all samples
        // at one thread; the answer at three.
        let fields: Vec<&str> = match (is_bsgs, one_thread) {
            (true, true) => BSGS_PINNED.to_vec(),
            (true, false) => vec!["answers_digest"],
            (false, true) => RHO_PINNED.to_vec(),
            (false, false) => vec!["log", "verified"],
        };
        let first = &samples[0].1;
        let mut mismatches = Vec::new();
        for (l, v) in samples {
            for f in &fields {
                if v[*f] != first[*f] {
                    mismatches.push(format!("{} {}={} (want {})", v["label"], f, v[*f], first[*f]));
                }
            }
            if !is_bsgs && v["verified"] != json!(true) {
                mismatches.push(format!("{} not verified", l.arm));
            }
            if is_bsgs && v["solved"] != v["targets"] {
                mismatches.push(format!("{} unsolved", v["label"]));
            }
        }
        let contended = |arm: &str| {
            samples
                .iter()
                .filter(|(l, v)| {
                    l.arm == arm
                        && iso
                            .get(v["label"].as_str().unwrap())
                            .map(|r| r["run"]["contended"].as_bool().unwrap_or(true))
                            .unwrap_or(true)
                })
                .count()
        };
        let missing_iso = samples
            .iter()
            .filter(|(_, v)| !iso.contains_key(v["label"].as_str().unwrap()))
            .count();
        // Contended runs (or runs without an isolation record) are reported
        // and left out: they are never pooled with uncontended ones, and a
        // round's ratio needs both of its arms uncontended.
        let clean = |v: &Value| {
            iso.get(v["label"].as_str().unwrap())
                .map(|r| r["run"]["contended"] == json!(false) && r["run"]["exit_status"] == json!(0))
                .unwrap_or(false)
        };
        let mut timer_rows = Map::new();
        for t in timers {
            let by_arm = |arm: &str| -> BTreeMap<u32, f64> {
                samples
                    .iter()
                    .filter(|(l, v)| l.arm == arm && clean(v))
                    .map(|(l, v)| (l.round, v[*t].as_f64().unwrap()))
                    .collect()
            };
            let (a, a2, b) = (by_arm("A"), by_arm("A2"), by_arm("B"));
            let vals = |m: &BTreeMap<u32, f64>| m.values().cloned().collect::<Vec<_>>();
            let paired = |x: &BTreeMap<u32, f64>, y: &BTreeMap<u32, f64>| -> Vec<f64> {
                x.iter()
                    .filter_map(|(r, xv)| y.get(r).map(|yv| xv / yv))
                    .collect()
            };
            let ab = paired(&a, &b);
            let aa = paired(&a, &a2);
            let verdict = if ab.is_empty() || aa.is_empty() {
                "incomplete"
            } else if min(&ab) > max(&aa) && min(&ab) > 1.0 {
                "improves"
            } else if max(&ab) < min(&aa) && max(&ab) < 1.0 {
                "regresses"
            } else {
                "no measurable change"
            };
            let (ma, ma2, mb) = (median(&vals(&a)), median(&vals(&a2)), median(&vals(&b)));
            println!(
                "| {part} | {bits} | {kernel} | {t} | {ma:.4} | {ma2:.4} | {mb:.4} | {:.3} | [{:.3}, {:.3}] | [{:.3}, {:.3}] | {verdict} | {} | ({}, {}, {}) |",
                ma / mb,
                min(&ab),
                max(&ab),
                min(&aa),
                max(&aa),
                if mismatches.is_empty() { "identical" } else { "MISMATCH" },
                contended("A"),
                contended("A2"),
                contended("B"),
            );
            timer_rows.insert(
                t.to_string(),
                json!({
                    "rounds": ab.len(),
                    "A": {"median": ma, "min": min(&vals(&a)), "values": vals(&a)},
                    "A2": {"median": ma2, "min": min(&vals(&a2)), "values": vals(&a2)},
                    "B": {"median": mb, "min": min(&vals(&b)), "values": vals(&b)},
                    "ratio_of_medians_A_over_B": ma / mb,
                    "ratio_of_medians_A_over_A2": ma / ma2,
                    "per_round_A_over_B": {"min": min(&ab), "max": max(&ab), "values": ab},
                    "per_round_A_over_A2": {"min": min(&aa), "max": max(&aa), "values": aa},
                    "verdict": verdict,
                }),
            );
        }
        let counted: Map<String, Value> = fields
            .iter()
            .map(|f| (f.to_string(), first[*f].clone()))
            .collect();
        rows.push(json!({
            "part": part,
            "bits": bits,
            "kernel": kernel,
            "threads": if one_thread { 1 } else { 3 },
            "samples": samples.len(),
            "pinned_fields": fields,
            "pinned_values": counted,
            "pinned_identical": mismatches.is_empty(),
            "mismatches": mismatches,
            "contended": {"A": contended("A"), "A2": contended("A2"), "B": contended("B")},
            "missing_isolation_records": missing_iso,
            "timers": timer_rows,
        }));
    }

    // External reference (AGENTS.md §15): PARI/GP `elllog` against the
    // native arms on the same curve and target, both single-threaded.
    let pari = read_jsonl(&raw.join("pari.jsonl"));
    let mut by_curve: BTreeMap<(String, u64), Vec<&Value>> = BTreeMap::new();
    for v in &pari {
        by_curve
            .entry((v["role"].as_str().unwrap().to_string(), v["bits"].as_u64().unwrap()))
            .or_default()
            .push(v);
    }
    let native = |part: &str, bits: u64, kernel: &str, timer: &str, arm: &str| -> Option<f64> {
        rows.iter()
            .find(|r| r["part"] == json!(part) && r["bits"] == json!(bits) && r["kernel"] == json!(kernel))
            .and_then(|r| r["timers"][timer][arm]["median"].as_f64())
    };
    let mut external = Vec::new();
    if !pari.is_empty() {
        println!();
        println!("| curve role | bits | PARI elllog median s (n) | native A median s | native B median s | B / PARI | A / PARI | PARI answers verified |");
        println!("|---|---:|---:|---:|---:|---:|---:|---|");
    }
    for ((role, bits), vs) in &by_curve {
        let secs: Vec<f64> = vs.iter().map(|v| v["elllog_ms"].as_f64().unwrap() / 1e3).collect();
        let verified = vs.iter().all(|v| v["verified"] == json!(true));
        let (part, kernel, timer) = if role == "bsgs" {
            ("bsgs1", "paired", "total_s")
        } else {
            ("rho1", "-", "solve_s")
        };
        let (a, b) = (native(part, *bits, kernel, timer, "A"), native(part, *bits, kernel, timer, "B"));
        let mp = median(&secs);
        println!(
            "| {role} | {bits} | {mp:.3} ({}) | {} | {} | {} | {} | {verified} |",
            secs.len(),
            a.map_or("-".into(), |x| format!("{x:.4}")),
            b.map_or("-".into(), |x| format!("{x:.4}")),
            b.map_or("-".into(), |x| format!("{:.4}", x / mp)),
            a.map_or("-".into(), |x| format!("{:.4}", x / mp)),
        );
        external.push(json!({
            "role": role,
            "bits": bits,
            "tool": vs[0]["tool"],
            "pari_elllog_s": {"median": mp, "min": min(&secs), "values": secs},
            "pari_verified": verified,
            "native_arm": {"part": part, "kernel": kernel, "timer": timer},
            "native_A_median_s": a,
            "native_B_median_s": b,
            "B_over_pari": b.map(|x| x / mp),
            "A_over_pari": a.map(|x| x / mp),
        }));
    }
    let opcount = read_jsonl(&raw.join("opcount.jsonl"));
    let summary = json!({
        "schema": "generic-dlp-summary/v1",
        "run": run.file_name().map(|s| s.to_string_lossy().to_string()),
        "rule": "improves: every per-round A/B ratio above the largest per-round A/A2 ratio (and above 1); regresses: every one below the smallest (and below 1); otherwise no measurable change",
        "configs": rows,
        "external": external,
        "pari": pari,
        "opcount": opcount,
    });
    std::fs::write(
        run.join("summary.json"),
        serde_json::to_string_pretty(&summary).unwrap() + "\n",
    )
    .expect("write summary.json");
    eprintln!("wrote {}", run.join("summary.json").display());
}
