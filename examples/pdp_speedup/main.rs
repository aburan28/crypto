//! Harness for the PDP-oracle speedup program
//! (`research/pdp_oracle_speedup_20261011/`).
//!
//! ```bash
//! cargo run --release --example pdp_speedup -- corpus  --l 8 --out corpus_l8.txt [--threads 3]
//! cargo run --release --example pdp_speedup -- verify  --l 8 --corpus corpus_l8.txt
//! cargo run --release --example pdp_speedup -- brute   --l 8 --corpus corpus_l8.txt [--limit N]
//! cargo run --release --example pdp_speedup -- time    --l 11 --corpus corpus_l11.txt --arm baseline
//! cargo run --release --example pdp_speedup -- one     --l 8 --corpus corpus_l8.txt   # for callgrind
//! cargo run --release --example pdp_speedup --features pdp-opcount -- count --l 8 --corpus corpus_l8.txt
//! ```
//!
//! This is an oracle-stage measurement (`AGENTS.md` §8): it times the
//! pairs-and-solve decision of one abscissa, never an end-to-end ECDLP.

mod algebra;
mod arms;
mod bitslice;
mod bspairs;
mod common;
mod corpus;
mod gatecount;
mod group;
mod kernel;

use arms::Arm;
use common::*;
use corpus::*;
use crypto_lib::cryptanalysis::semaev_decomp::{decompose, eval_f3};
use std::hint::black_box;
use std::time::Instant;

fn arg<T: std::str::FromStr>(args: &[String], name: &str) -> Option<T> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
}

fn arg_s(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

fn load(args: &[String]) -> (Rung, Loaded) {
    let l: u32 = arg(args, "--l").expect("--l");
    let path = arg_s(args, "--corpus").expect("--corpus");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let loaded = read_corpus(&text);
    let rung = Rung::new(l);
    assert_eq!(loaded.l, l, "corpus is for a different rung");
    assert_eq!(loaded.modulus, rung.modulus(), "corpus modulus differs");
    (rung, loaded)
}

fn cmd_corpus(args: &[String]) {
    let l: u32 = arg(args, "--l").expect("--l");
    let out = arg_s(args, "--out").expect("--out");
    let threads: usize = arg(args, "--threads").unwrap_or(3);
    let rung = Rung::new(l);
    let t0 = Instant::now();
    let recs = build_records(&rung, threads);
    let text = write_corpus(&rung, &recs);
    std::fs::write(&out, &text).unwrap();
    let found = recs.iter().filter(|r| !r.roots.is_empty()).count();
    let roots: usize = recs.iter().map(|r| r.roots.len()).sum();
    let rational: usize = recs
        .iter()
        .map(|r| r.rational.iter().filter(|&&b| b).count())
        .sum();
    println!(
        "l={l} n={} targets={} found={found} refuted={} roots={roots} rational_roots={rational} enumerate_secs={:.1} threads={threads}",
        rung.n,
        recs.len(),
        recs.len() - found,
        t0.elapsed().as_secs_f64()
    );
    println!("digest_verdicts {}", digest_verdicts(&recs));
    println!("digest_roots {}", digest_roots(&recs));
}

fn cmd_verify(args: &[String]) {
    let (rung, loaded) = load(args);
    let recs = &loaded.recs;
    let mut bad = 0usize;
    // Integrity of the corpus file itself.
    let dv_file = digest_verdicts(recs);
    let dr_file = digest_roots(recs);
    if dv_file != loaded.file_digest_verdicts || dr_file != loaded.file_digest_roots {
        println!("FAIL corpus digests do not match its own contents");
        bad += 1;
    }
    // Replay the baseline oracle and rebuild the verdict stream from what it said.
    let mut replay = recs.to_vec();
    let t0 = Instant::now();
    for rec in replay.iter_mut() {
        let got = decompose(rec.t.xr, rung.l, &rung.gf);
        let truth = !rec.roots.is_empty();
        if got.is_some() != truth {
            println!(
                "FAIL target {} verdict {:?} vs corpus found={truth}",
                rec.t.idx, got
            );
            bad += 1;
        }
        if let Some(w) = got {
            let s = sorted3(w[0], w[1], w[2]);
            if rec.roots.binary_search(&s).is_err() {
                println!(
                    "FAIL target {} witness {:?} not in the root set",
                    rec.t.idx, w
                );
                bad += 1;
            }
            if eval_f3(w[0], w[1], w[2], rec.t.xr, &rung.gf) != 0 {
                println!("FAIL target {} witness does not vanish", rec.t.idx);
                bad += 1;
            }
        }
        // The verdict stream is built from the baseline's answer.
        if got.is_none() {
            rec.roots.clear();
        } else if rec.roots.is_empty() {
            rec.roots.push([0, 0, 0]);
        }
    }
    let replay_dv = digest_verdicts(&replay);
    if replay_dv != loaded.file_digest_verdicts {
        println!("FAIL replayed verdict digest {replay_dv} != corpus digest");
        bad += 1;
    }
    println!(
        "l={} targets={} replay_secs={:.3} digest_verdicts {} digest_roots {} {}",
        rung.l,
        recs.len(),
        t0.elapsed().as_secs_f64(),
        replay_dv,
        dr_file,
        if bad == 0 { "OK" } else { "FAILED" }
    );
    if bad != 0 {
        std::process::exit(1);
    }
}

fn cmd_brute(args: &[String]) {
    let (rung, loaded) = load(args);
    let limit: usize = arg(args, "--limit").unwrap_or(usize::MAX);
    let threads: usize = arg(args, "--threads").unwrap_or(3);
    let mut bad = 0;
    for rec in loaded.recs.iter().take(limit) {
        let got = brute_roots(&rung, rec.t.xr, threads);
        let same = got == rec.roots;
        if !same {
            bad += 1;
        }
        println!(
            "target {:>2} {} brute_roots={} corpus_roots={} {}",
            rec.t.idx,
            if rec.t.planted { "planted" } else { "uniform" },
            got.len(),
            rec.roots.len(),
            if same { "same" } else { "DIFFERENT" }
        );
    }
    println!("brute-force cross-check: {} mismatching targets", bad);
    if bad != 0 {
        std::process::exit(1);
    }
}

fn cmd_time(args: &[String]) {
    let (rung, loaded) = load(args);
    let arm_name = arg_s(args, "--arm").unwrap_or_else(|| "baseline".into());
    let k: usize = arg(args, "--targets").unwrap_or(16);
    let passes_min: usize = arg(args, "--passes").unwrap_or(1);
    let min_secs: f64 = arg(args, "--min-secs").unwrap_or(0.0);
    let refuted: Vec<(u64, u64)> = loaded
        .recs
        .iter()
        .filter(|r| r.roots.is_empty())
        .map(|r| (r.t.xr, r.t.yr))
        .take(k)
        .collect();
    assert!(!refuted.is_empty(), "no refuted targets in the corpus");
    // Setup (tables of the arm) is built once, outside the timed region, and reported.
    let ts = Instant::now();
    let arm = Arm::build(&arm_name, &rung);
    let setup_secs = ts.elapsed().as_secs_f64();
    // One untimed call on the first target warms the instruction cache and the
    // page tables; it is reported, not hidden.
    black_box(arm.decide(refuted[0].0, refuted[0].1));
    let t0 = Instant::now();
    let mut passes = 0usize;
    loop {
        for &(xr, yr) in &refuted {
            assert!(!black_box(arm.decide(black_box(xr), black_box(yr))));
        }
        passes += 1;
        if passes >= passes_min && t0.elapsed().as_secs_f64() >= min_secs {
            break;
        }
    }
    let secs = t0.elapsed().as_secs_f64();
    let calls = (refuted.len() * passes) as u64;
    let pairs = arm.pairs_total() * calls;
    let label = arg_s(args, "--label").unwrap_or_default();
    let line = format!(
        "{{\"label\":\"{label}\",\"arm\":\"{arm_name}\",\"l\":{},\"n\":{},\"targets\":{},\"passes\":{passes},\"pairs_per_target\":{},\"pairs\":{pairs},\"secs\":{secs:.6},\"setup_secs\":{setup_secs:.6},\"ns_per_pair\":{:.3},\"ns_per_target\":{:.1},\"kernel\":\"{}\"}}",
        rung.l,
        rung.n,
        refuted.len(),
        arm.pairs_total(),
        secs * 1e9 / pairs as f64,
        secs * 1e9 / calls as f64,
        rung.gf.kernel_name()
    );
    println!("{line}");
    if let Some(out) = arg_s(args, "--out") {
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(out)
            .unwrap();
        writeln!(f, "{line}").unwrap();
    }
}

/// `verify-arm --arm A --l L --corpus F`: the arm's verdict on every target
/// against the corpus, and, for arms that can list roots, the full root set.
fn cmd_verify_arm(args: &[String]) {
    let (rung, loaded) = load(args);
    let arm_name = arg_s(args, "--arm").expect("--arm");
    let arm = Arm::build(&arm_name, &rung);
    let mut bad = 0usize;
    let mut replay = loaded.recs.clone();
    let mut roots_checked = 0usize;
    let t0 = Instant::now();
    for rec in replay.iter_mut() {
        let truth = !rec.roots.is_empty();
        let got = arm.decide(rec.t.xr, rec.t.yr);
        if got != truth {
            println!("FAIL target {} verdict {got} vs corpus {truth}", rec.t.idx);
            bad += 1;
        }
        if let Some(all) = arm.all_roots(rec.t.xr, rec.t.yr) {
            roots_checked += 1;
            if all != rec.roots {
                println!("FAIL target {} root set differs", rec.t.idx);
                bad += 1;
            }
        }
        // The verdict stream is built from the arm's answer.
        if !got {
            rec.roots.clear();
        } else if rec.roots.is_empty() {
            rec.roots.push([0, 0, 0]);
        }
    }
    let dv = digest_verdicts(&replay);
    if dv != loaded.file_digest_verdicts {
        println!("FAIL verdict digest {dv} != corpus digest");
        bad += 1;
    }
    println!(
        "{{\"arm\":\"{arm_name}\",\"l\":{},\"targets\":{},\"root_sets_checked\":{roots_checked},\"secs\":{:.3},\"digest_verdicts\":\"{dv}\",\"digest_roots_corpus\":\"{}\",\"mismatches\":{bad}}}",
        rung.l,
        loaded.recs.len(),
        t0.elapsed().as_secs_f64(),
        loaded.file_digest_roots
    );
    if bad != 0 {
        std::process::exit(1);
    }
}

#[inline(never)]
fn run_one(rung: &Rung, xr: u64) -> bool {
    decompose(xr, rung.l, &rung.gf).is_some()
}

/// One refuted target, once, for `valgrind --tool=callgrind`.
fn cmd_one(args: &[String]) {
    let (rung, loaded) = load(args);
    let rec = loaded
        .recs
        .iter()
        .find(|r| r.roots.is_empty())
        .expect("a refuted target");
    let found = run_one(&rung, rec.t.xr);
    println!("l={} pairs={} found={found}", rung.l, rung.pairs_total());
}

/// Median of a slice (sorted copy).
fn median(v: &[f64]) -> f64 {
    let mut w = v.to_vec();
    w.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = w.len();
    if n % 2 == 1 {
        w[n / 2]
    } else {
        0.5 * (w[n / 2 - 1] + w[n / 2])
    }
}

/// Percentile bootstrap (2.5%, 97.5%) of the median of `v`, resampling
/// rounds with replacement; deterministic (SplitMix64 seeded 20261011).
fn bootstrap_median_ci(v: &[f64]) -> (f64, f64) {
    let mut rng = SplitMix64(20261011);
    let mut meds = Vec::new();
    for _ in 0..10_000 {
        let r: Vec<f64> = (0..v.len())
            .map(|_| v[(rng.next() % v.len() as u64) as usize])
            .collect();
        meds.push(median(&r));
    }
    meds.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (meds[250], meds[9750])
}

/// Join the harness's own timings with `isolated_bench.py` records by label
/// and print the A/A or A/B table for every rung.  Labels are
/// `<arm tag>-l<L>-r<round>`; `--a TAG --b TAG` pick the two sides.
fn cmd_summarize(args: &[String]) {
    let times = std::fs::read_to_string(arg_s(args, "--times").expect("--times")).unwrap();
    let iso = std::fs::read_to_string(arg_s(args, "--iso").expect("--iso")).unwrap();
    let metric = arg_s(args, "--metric").unwrap_or_else(|| "ns_per_pair".into());
    let (a_tag, b_tag) = (
        arg_s(args, "--a").expect("--a"),
        arg_s(args, "--b").expect("--b"),
    );
    use std::collections::BTreeMap;
    let mut contended: BTreeMap<String, bool> = BTreeMap::new();
    for line in iso.lines() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        let label = v["label"].as_str().unwrap_or("").to_string();
        let c = v["run"]["contended"].as_bool().unwrap_or(true);
        contended.insert(label, c);
    }
    // (l, round) -> (a ns, b ns)
    let mut rows: BTreeMap<(u32, u32), (Option<f64>, Option<f64>)> = BTreeMap::new();
    let mut n_cont: BTreeMap<u32, usize> = BTreeMap::new();
    for line in times.lines() {
        let v: serde_json::Value = serde_json::from_str(line).unwrap();
        let label = v["label"].as_str().unwrap().to_string();
        let l = v["l"].as_u64().unwrap() as u32;
        let ns = v[metric.as_str()].as_f64().unwrap();
        // label = TAG-lL-rR
        let parts: Vec<&str> = label.rsplitn(3, '-').collect();
        let round: u32 = parts[0].trim_start_matches('r').parse().unwrap();
        let tag = parts[2];
        if contended.get(&label).copied().unwrap_or(true) {
            *n_cont.entry(l).or_default() += 1;
        }
        let e = rows.entry((l, round)).or_default();
        if tag == a_tag {
            e.0 = Some(ns);
        } else if tag == b_tag {
            e.1 = Some(ns);
        }
    }
    println!("| l | rounds | contended runs | A median {metric} | B median {metric} | A run-to-run spread | B/A median | B/A min..max | B/A 95% bootstrap CI |");
    println!("|--:|--:|--:|--:|--:|--:|--:|:--|:--|");
    let ls: Vec<u32> = rows
        .keys()
        .map(|k| k.0)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    for l in ls {
        let pairs: Vec<(f64, f64)> = rows
            .iter()
            .filter(|(k, v)| k.0 == l && v.0.is_some() && v.1.is_some())
            .map(|(_, v)| (v.0.unwrap(), v.1.unwrap()))
            .collect();
        let a: Vec<f64> = pairs.iter().map(|p| p.0).collect();
        let b: Vec<f64> = pairs.iter().map(|p| p.1).collect();
        let ratio: Vec<f64> = pairs.iter().map(|p| p.1 / p.0).collect();
        let (lo, hi) = bootstrap_median_ci(&ratio);
        let amin = a.iter().cloned().fold(f64::MAX, f64::min);
        let amax = a.iter().cloned().fold(0.0, f64::max);
        println!(
            "| {l} | {} | {} | {:.1} | {:.1} | {:.2}% | {:.4} | {:.4}..{:.4} | {:.4}..{:.4} |",
            pairs.len(),
            n_cont.get(&l).copied().unwrap_or(0),
            median(&a),
            median(&b),
            100.0 * (amax - amin) / median(&a),
            median(&ratio),
            ratio.iter().cloned().fold(f64::MAX, f64::min),
            ratio.iter().cloned().fold(0.0, f64::max),
            lo,
            hi
        );
    }
}

/// `count --l L --corpus F [--arm A] [--targets K] [--what pairs|formation]`
/// (build with `--features pdp-opcount`): field operations per pair of an arm
/// on `K` refuted targets, or per call of the quartic formation
/// `quartic_in_x3` (the part an incremental coefficient update would replace).
#[cfg(feature = "pdp-opcount")]
fn cmd_count(args: &[String]) {
    use crypto_lib::cryptanalysis::semaev_decomp::{opcount, quartic_in_x3};
    let (rung, loaded) = load(args);
    let k: usize = arg(args, "--targets").unwrap_or(4);
    let arm_name = arg_s(args, "--arm").unwrap_or_else(|| "baseline".into());
    let what = arg_s(args, "--what").unwrap_or_else(|| "pairs".into());
    let refuted: Vec<(u64, u64)> = loaded
        .recs
        .iter()
        .filter(|r| r.roots.is_empty())
        .map(|r| (r.t.xr, r.t.yr))
        .take(k)
        .collect();
    assert!(rung.gf.folds(), "counts are defined for the folding path");
    if what == "formation" {
        let calls = 100_000u64;
        let mut rng = SplitMix64(20261011);
        let (xr, _) = refuted[0];
        opcount::reset();
        let mut sink = 0u64;
        for _ in 0..calls {
            let x1 = rng.next() & (rung.span() - 1);
            let x2 = rng.next() & (rung.span() - 1);
            let q = quartic_in_x3(x1, x2, xr, &rung.gf);
            sink ^= q.c[0] ^ q.c[4];
        }
        black_box(sink);
        let c = opcount::snapshot();
        println!(
            "{{\"l\":{},\"what\":\"formation\",\"calls\":{calls},\"counts\":{},\"per_call\":{}}}",
            rung.l,
            c.json(),
            c.per_pair_json(calls as f64)
        );
        return;
    }
    let arm = Arm::build(&arm_name, &rung);
    opcount::reset();
    for &(xr, yr) in &refuted {
        assert!(!arm.decide(xr, yr));
    }
    let c = opcount::snapshot();
    let pairs = (arm.pairs_total() * refuted.len() as u64) as f64;
    println!(
        "{{\"l\":{},\"arm\":\"{arm_name}\",\"targets\":{},\"pairs\":{},\"counts\":{},\"per_pair\":{}}}",
        rung.l,
        refuted.len(),
        pairs as u64,
        c.json(),
        c.per_pair_json(pairs)
    );
}

#[cfg(not(feature = "pdp-opcount"))]
fn cmd_count(_args: &[String]) {
    eprintln!("build with --features pdp-opcount");
    std::process::exit(2);
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("corpus") => cmd_corpus(&args),
        Some("verify") => cmd_verify(&args),
        Some("brute") => cmd_brute(&args),
        Some("time") => cmd_time(&args),
        Some("one") => cmd_one(&args),
        Some("verify-arm") => cmd_verify_arm(&args),
        Some("count") => cmd_count(&args),
        Some("summarize") => cmd_summarize(&args),
        Some("degrees") => algebra::cmd_degrees(&args),
        Some("agree") => algebra::cmd_agree(&args),
        Some("sample") => algebra::cmd_sample(&args),
        Some("gates") => gatecount::cmd_gates(&args),
        Some("export-program") => gatecount::cmd_export(&args),
        Some("torsion") => algebra::cmd_torsion(&args),
        Some("corpus-check") => algebra::cmd_corpus_check(&args),
        Some("bs-check") => gatecount::cmd_bscheck(&args),
        _ => {
            eprintln!("usage: pdp_speedup corpus|verify|brute|time|one|count|summarize ...");
            std::process::exit(2);
        }
    }
}
