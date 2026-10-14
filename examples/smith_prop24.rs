//! Smith's Propositions 2 and 4, natively (`cryptanalysis::smith_prop24`):
//! the Example 1 positive control, random curves over a large prime
//! (Question 1) and the exhaustive per-isogeny-class census over small primes
//! (Question 2) of `research/smith_prop24_native_20261014/PROTOCOL.md`.
//!
//! ```text
//! cargo run --release --example smith_prop24 -- control --out DIR
//! cargo run --release --example smith_prop24 -- random --p 2147483647 --n 200000 --seed 20261014 --out DIR
//! cargo run --release --example smith_prop24 -- census --p 7 --perms 199 --seed 20261014 --out DIR [--limit CODES]
//! cargo run --release --example smith_prop24 -- resummarize --dir DIR
//! ```
//!
//! `random` writes `random_p<p>_rows.csv.gz` (one row per curve) and
//! `random_p<p>_summary.json`; `census` writes `census_p<p>_classes.csv`,
//! `census_p<p>_types.csv`, `census_p<p>_raw.json` and
//! `census_p<p>_summary.json`; each run appends to `run_manifest.json`.
//! `resummarize` recomputes every summary from the committed rows and fails
//! if any differs.

use crypto_lib::cryptanalysis::smith_prop24 as sp;
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use rayon::prelude::*;
use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::{fs, process::Command, time::Instant};

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter()
        .position(|a| a == k)
        .and_then(|i| args.get(i + 1).cloned())
}

fn git(args: &[&str]) -> Option<String> {
    Command::new("git")
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
}

fn write_json(path: &Path, v: &impl serde::Serialize) {
    fs::write(path, serde_json::to_string_pretty(v).unwrap() + "\n").unwrap();
}

/// The stored file is byte-identical to the freshly serialised summary.
fn same_text(path: &Path, v: &impl serde::Serialize) -> bool {
    fs::read_to_string(path).ok() == Some(serde_json::to_string_pretty(v).unwrap() + "\n")
}

fn append_manifest(out: &Path, entry: Value) {
    let path = out.join("run_manifest.json");
    let mut runs: Vec<Value> = fs::read_to_string(&path)
        .ok()
        .and_then(|s| serde_json::from_str(&s).ok())
        .unwrap_or_default();
    runs.push(entry);
    write_json(&path, &runs);
}

fn manifest_entry(args: &[String], run: &str, secs: f64) -> Value {
    json!({
        "run": run,
        "command": args,
        "git_head": git(&["rev-parse", "HEAD"]),
        "git_dirty_paths": git(&["status", "--porcelain", "--", "src", "examples"]),
        "protocol_sha256": Command::new("sha256sum")
            .arg("research/smith_prop24_native_20261014/PROTOCOL.md")
            .output().ok().filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).split_whitespace().next().unwrap_or("").to_string()),
        "threads": rayon::current_num_threads(),
        "seconds": secs,
        "rustc": Command::new("rustc").arg("--version").output().ok()
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string()),
    })
}

fn boot_seed(seed: u64) -> u64 {
    seed ^ 0xB007_5EED_0000_0001
}

const BOOT: u64 = 999;

fn read_random_rows(path: &Path) -> Vec<sp::RandomRow> {
    let r = BufReader::new(GzDecoder::new(fs::File::open(path).unwrap()));
    r.lines()
        .skip(1)
        .map(|l| sp::RandomRow::parse(&l.unwrap()).expect("row"))
        .collect()
}

fn census_files(out: &Path, p: u64) -> [PathBuf; 4] {
    [
        out.join(format!("census_p{p}_classes.csv")),
        out.join(format!("census_p{p}_types.csv")),
        out.join(format!("census_p{p}_raw.json")),
        out.join(format!("census_p{p}_summary.json")),
    ]
}

fn census_summary_from_files(out: &Path, p: u64) -> (sp::CensusSummary, Value) {
    let [cl, ty, raw, sum] = census_files(out, p);
    let classes: Vec<sp::ClassRow> = fs::read_to_string(cl)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| sp::ClassRow::parse(l).expect("class row"))
        .collect();
    let types: Vec<sp::TypeRow> = fs::read_to_string(ty)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| sp::TypeRow::parse(l).expect("type row"))
        .collect();
    let raw: Value = serde_json::from_str(&fs::read_to_string(raw).unwrap()).unwrap();
    let stored: Value = serde_json::from_str(&fs::read_to_string(sum).unwrap()).unwrap();
    let strat: Vec<u64> = serde_json::from_value(raw["stratum_successes"].clone()).unwrap();
    let diag: sp::Diag = serde_json_diag(&raw["diag"]);
    let s = sp::summarize_census(
        p,
        raw["octics_enumerated"].as_u64().unwrap(),
        &classes,
        &types,
        &strat.try_into().unwrap(),
        diag,
        stored["null"]["permutations"].as_u64().unwrap(),
        stored["null"]["seed"].as_u64().unwrap(),
    );
    (s, stored)
}

fn serde_json_diag(v: &Value) -> sp::Diag {
    let g = |k: &str| v[k].as_u64().unwrap();
    sp::Diag {
        maps: g("maps"),
        m_rank_deficient: g("m_rank_deficient"),
        disc_zero: g("disc_zero"),
        pencil_in_grassmannian: g("pencil_in_grassmannian"),
        non_normal_form: g("non_normal_form"),
        degree_drop: g("degree_drop"),
        pair_condition_fail: g("pair_condition_fail"),
        eval_mismatch: g("eval_mismatch"),
        lemma4_fail: g("lemma4_fail"),
        lambda_disagree: g("lambda_disagree"),
        construction_failed: g("construction_failed"),
        x_agree: g("x_agree"),
        x_disagree: g("x_disagree"),
        x_not_applicable: g("x_not_applicable"),
        x_count_error: g("x_count_error"),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).cloned().unwrap_or_default();
    match cmd.as_str() {
        "control" => {
            let out = PathBuf::from(arg(&args, "--out").expect("--out DIR"));
            fs::create_dir_all(&out).unwrap();
            let t0 = Instant::now();
            let rep = sp::example1();
            let secs = t0.elapsed().as_secs_f64();
            write_json(&out.join("control_example1.json"), &rep);
            append_manifest(&out, manifest_entry(&args, "control_example1", secs));
            eprintln!("Example 1 control passed: {}", rep.passed);
            if !rep.passed {
                std::process::exit(2);
            }
        }
        "random" => {
            let out = PathBuf::from(arg(&args, "--out").expect("--out DIR"));
            let p: u64 = arg(&args, "--p").expect("--p").parse().unwrap();
            let n: u64 = arg(&args, "--n").expect("--n").parse().unwrap();
            let seed: u64 = arg(&args, "--seed").expect("--seed").parse().unwrap();
            fs::create_dir_all(&out).unwrap();
            let t0 = Instant::now();
            let rows: Vec<sp::RandomRow> = (0..n)
                .into_par_iter()
                .map(|i| sp::random_row(p, seed, i))
                .collect();
            let secs = t0.elapsed().as_secs_f64();
            let mut gz = GzEncoder::new(
                fs::File::create(out.join(format!("random_p{p}_rows.csv.gz"))).unwrap(),
                Compression::best(),
            );
            writeln!(gz, "{}", sp::RANDOM_HEADER).unwrap();
            for r in &rows {
                writeln!(gz, "{}", r.csv()).unwrap();
            }
            gz.finish().unwrap();
            let summary = sp::summarize_random(&rows, p, seed, BOOT, boot_seed(seed));
            write_json(&out.join(format!("random_p{p}_summary.json")), &summary);
            append_manifest(&out, manifest_entry(&args, &format!("random_p{p}"), secs));
            eprintln!("random p={p} n={n}: {secs:.1}s");
        }
        "census" => {
            let out = PathBuf::from(arg(&args, "--out").expect("--out DIR"));
            let p: u64 = arg(&args, "--p").expect("--p").parse().unwrap();
            let perms: u64 = arg(&args, "--perms").expect("--perms").parse().unwrap();
            let seed: u64 = arg(&args, "--seed").expect("--seed").parse().unwrap();
            let limit: Option<u64> = arg(&args, "--limit").map(|s| s.parse().unwrap());
            let xcheck_every: u64 = arg(&args, "--xcheck-every")
                .map(|s| s.parse().unwrap())
                .unwrap_or(1);
            fs::create_dir_all(&out).unwrap();
            let t0 = Instant::now();
            let data = sp::census(p, xcheck_every, limit);
            let secs = t0.elapsed().as_secs_f64();
            let [cl, ty, raw, sum] = census_files(&out, p);
            let mut s = String::from(sp::CLASS_HEADER);
            s.push('\n');
            for c in &data.classes {
                s.push_str(&c.csv());
                s.push('\n');
            }
            fs::write(cl, s).unwrap();
            let mut s = String::from(sp::TYPE_HEADER);
            s.push('\n');
            for t in &data.types {
                s.push_str(&t.csv());
                s.push('\n');
            }
            fs::write(ty, s).unwrap();
            write_json(
                &raw,
                &json!({
                    "p": p,
                    "octics_enumerated": data.octics_enumerated,
                    "xcheck_every": xcheck_every,
                    "twist_constant": data.twist_constant,
                    "stratum_values": sp::STRATA,
                    "stratum_successes": data.stratum_successes,
                    "diag": data.diag,
                }),
            );
            let summary = sp::summarize_census(
                p,
                data.octics_enumerated,
                &data.classes,
                &data.types,
                &data.stratum_successes,
                data.diag,
                perms,
                seed,
            );
            write_json(&sum, &summary);
            append_manifest(&out, manifest_entry(&args, &format!("census_p{p}"), secs));
            eprintln!("census p={p}: {secs:.1}s");
        }
        "resummarize" => {
            let dir = PathBuf::from(arg(&args, "--dir").expect("--dir DIR"));
            let mut ok = true;
            for e in fs::read_dir(&dir).unwrap() {
                let name = e.unwrap().file_name().to_string_lossy().to_string();
                if let Some(p) = name
                    .strip_prefix("random_p")
                    .and_then(|s| s.strip_suffix("_rows.csv.gz"))
                {
                    let p: u64 = p.parse().unwrap();
                    let stored: Value = serde_json::from_str(
                        &fs::read_to_string(dir.join(format!("random_p{p}_summary.json"))).unwrap(),
                    )
                    .unwrap();
                    let rows = read_random_rows(&dir.join(&name));
                    let seed = stored["seed"].as_u64().unwrap();
                    let s = sp::summarize_random(&rows, p, seed, BOOT, boot_seed(seed));
                    let same = same_text(&dir.join(format!("random_p{p}_summary.json")), &s);
                    eprintln!(
                        "random p={p}: {}",
                        if same { "identical" } else { "DIFFERS" }
                    );
                    if !same {
                        write_json(
                            &dir.join(format!("random_p{p}_summary.recomputed.json")),
                            &s,
                        );
                    }
                    ok &= same;
                }
                if let Some(p) = name
                    .strip_prefix("census_p")
                    .and_then(|s| s.strip_suffix("_classes.csv"))
                {
                    let p: u64 = p.parse().unwrap();
                    let (s, _) = census_summary_from_files(&dir, p);
                    let same = same_text(&census_files(&dir, p)[3], &s);
                    eprintln!(
                        "census p={p}: {}",
                        if same { "identical" } else { "DIFFERS" }
                    );
                    if !same {
                        write_json(
                            &dir.join(format!("census_p{p}_summary.recomputed.json")),
                            &s,
                        );
                    }
                    ok &= same;
                }
            }
            if !ok {
                std::process::exit(1);
            }
        }
        "verdict" => {
            // Protocol Section 5 verdicts, from the recorded files only.
            let dir = PathBuf::from(arg(&args, "--dir").expect("--dir DIR"));
            let ctl: Value = serde_json::from_str(
                &fs::read_to_string(dir.join("control_example1.json")).unwrap(),
            )
            .unwrap();
            let q1: Value = serde_json::from_str(
                &fs::read_to_string(dir.join("random_p2147483647_summary.json")).unwrap(),
            )
            .unwrap();
            let census: Vec<sp::CensusSummary> = [7u64, 11]
                .iter()
                .filter(|p| census_files(&dir, **p)[0].exists())
                .map(|&p| census_summary_from_files(&dir, p).0)
                .collect();
            let v = json!({
                "control_example1_passed": ctl["passed"],
                "q1_overall": q1["overall_verdict"],
                "q1_fractions": q1["fractions"].as_array().unwrap().iter().map(|f| json!({
                    "name": f["name"], "value": f["value"], "cp_lo": f["cp_lo"], "cp_hi": f["cp_hi"],
                    "prediction": f["prediction"], "verdict": f["verdict"]})).collect::<Vec<_>>(),
                "q2_per_prime": census.iter().map(|s| json!({
                    "p": s.p, "D": s.direct_full, "C": s.class_full, "U": s.class_tractable,
                    "r": s.r, "r_band": s.r_band, "verdict": s.verdict,
                    "null_r_median": s.null.r_median})).collect::<Vec<_>>(),
                "q2_overall": sp::overall_q2(&census),
            });
            write_json(&dir.join("verdict.json"), &v);
            eprintln!("{}", serde_json::to_string_pretty(&v).unwrap());
        }
        "describe" => {
            // Post hoc, descriptive (protocol amendment 1): how the
            // full_success curves spread over the classes that hold a curve
            // with a rational tractable subgroup.  No verdict rests on it.
            let dir = PathBuf::from(arg(&args, "--dir").expect("--dir DIR"));
            let mut out = Vec::new();
            for p in [7u64, 11] {
                let [cl, ..] = census_files(&dir, p);
                if !cl.exists() {
                    continue;
                }
                let classes: Vec<sp::ClassRow> = fs::read_to_string(cl)
                    .unwrap()
                    .lines()
                    .skip(1)
                    .map(|l| sp::ClassRow::parse(l).expect("class row"))
                    .collect();
                // Bins of share = full_success / tractable curves in the class:
                // exactly 0, (0, 0.1], (0.1, 0.2], ..., (0.9, 1].
                let mut classes_by_bin = [0u64; 11];
                let mut curves_by_bin = [0u64; 11];
                let mut zero_sizes: Vec<u64> = Vec::new();
                let mut pos_sizes: Vec<u64> = Vec::new();
                for c in classes.iter().filter(|c| c.tractable_curves > 0) {
                    let share = c.full_success_curves as f64 / c.tractable_curves as f64;
                    let b = if c.full_success_curves == 0 {
                        0
                    } else {
                        1 + (((share * 10.0).ceil() as usize).clamp(1, 10) - 1)
                    };
                    classes_by_bin[b] += 1;
                    curves_by_bin[b] += c.curves;
                    if c.full_success_curves == 0 {
                        zero_sizes.push(c.curves);
                    } else {
                        pos_sizes.push(c.curves);
                    }
                }
                let median = |v: &mut Vec<u64>| -> Option<u64> {
                    v.sort_unstable();
                    v.get(v.len() / 2).copied()
                };
                let full_total: u64 = classes.iter().map(|c| c.full_success_curves).sum();
                let tract_total: u64 = classes.iter().map(|c| c.tractable_curves).sum();
                out.push(json!({
                    "p": p,
                    "bins": ["0", "(0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", "(0.4,0.5]",
                             "(0.5,0.6]", "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]", "(0.9,1]"],
                    "classes_by_bin": classes_by_bin,
                    "curves_by_bin": curves_by_bin,
                    "tractable_classes": classes_by_bin.iter().sum::<u64>(),
                    "tractable_classes_without_full_success": classes_by_bin[0],
                    "curves_in_them": curves_by_bin[0],
                    "median_class_size_without_full_success": median(&mut zero_sizes),
                    "median_class_size_with_full_success": median(&mut pos_sizes),
                    "pooled_share_full_over_tractable": full_total as f64 / tract_total as f64,
                }));
            }
            write_json(&dir.join("describe_classes.json"), &out);
            eprintln!("{}", serde_json::to_string_pretty(&out).unwrap());
        }
        _ => {
            eprintln!("usage: smith_prop24 control|random|census|resummarize|verdict|describe ...");
            std::process::exit(64);
        }
    }
}
