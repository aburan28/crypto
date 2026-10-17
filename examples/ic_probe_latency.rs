//! Null objects and claim arithmetic for the index-calculus probe-latency
//! study (`research/ic_probe_latency_20261017/`).
//!
//! The question is whether the folded pair-state index calculus' online time
//! per probe is the host's random-access latency over the root table's bytes,
//! and whether that latency is what a cross-method ratio carries. Three
//! subcommands:
//!
//! - `chase`: dependent pointer chase and independent random gather over an
//!   array of a given byte size, on 4 KiB or 2 MiB (transparent huge) pages.
//!   No arithmetic: the null object for a table lookup.
//! - `geomfit`: the per-target probe-count distribution of an IC run against
//!   a geometric law with the same mean.
//! - `claims`: the preregistered claims from the run directories.
//!
//! Usage:
//!   ic_probe_latency chase --bytes N --pages 4k|2m [--accesses A] [--reps R] [--seed S] --out FILE
//!   ic_probe_latency geomfit --label L --targets targets_result.jsonl [--label L --targets F ...] --out FILE
//!   ic_probe_latency claims --stage0 FILE --archived FILE --answers FILE --a1 DIR --b1 DIR
//!       --a2 DIR --b2 DIR --rho DIR,DIR --out FILE

use serde_json::{json, Value};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

const LINE: usize = 64;
const WORDS_PER_LINE: usize = LINE / 8;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (cmd, rest) = args.split_first().unwrap_or_else(|| usage());
    let opts = Opts::parse(rest);
    let out = match cmd.as_str() {
        "chase" => chase(&opts),
        "geomfit" => geomfit(&opts),
        "claims" => claims(&opts),
        _ => usage(),
    };
    let path = opts.one("out");
    fs::write(&path, serde_json::to_string_pretty(&out).unwrap() + "\n")
        .unwrap_or_else(|e| panic!("write {path}: {e}"));
    eprintln!("wrote {path}");
}

fn usage() -> ! {
    eprintln!("usage: ic_probe_latency chase|geomfit|claims ... --out FILE (see the module docs)");
    std::process::exit(2);
}

/// `--key value` pairs; a key may repeat.
struct Opts(Vec<(String, String)>);

impl Opts {
    fn parse(args: &[String]) -> Self {
        let mut pairs = Vec::new();
        let mut it = args.iter();
        while let Some(key) = it.next() {
            let key = key
                .strip_prefix("--")
                .unwrap_or_else(|| panic!("expected --key, got {key}"));
            let value = it.next().unwrap_or_else(|| panic!("--{key} needs a value"));
            pairs.push((key.to_string(), value.clone()));
        }
        Self(pairs)
    }

    fn all(&self, key: &str) -> Vec<String> {
        self.0
            .iter()
            .filter(|(k, _)| k == key)
            .map(|(_, v)| v.clone())
            .collect()
    }

    fn get(&self, key: &str) -> Option<String> {
        let mut all = self.all(key);
        assert!(all.len() <= 1, "--{key} given more than once");
        all.pop()
    }

    fn one(&self, key: &str) -> String {
        self.get(key)
            .unwrap_or_else(|| panic!("--{key} is required"))
    }

    fn num(&self, key: &str, default: u64) -> u64 {
        self.get(key).map_or(default, |v| parse_size(&v))
    }
}

/// Accepts plain integers and the suffixes Ki, Mi, Gi.
fn parse_size(s: &str) -> u64 {
    let (digits, scale) = match s {
        _ if s.ends_with("Gi") => (&s[..s.len() - 2], 1 << 30),
        _ if s.ends_with("Mi") => (&s[..s.len() - 2], 1 << 20),
        _ if s.ends_with("Ki") => (&s[..s.len() - 2], 1 << 10),
        _ => (s, 1),
    };
    digits
        .parse::<u64>()
        .unwrap_or_else(|_| panic!("bad size {s}"))
        * scale
}

// ---------------------------------------------------------------- chase

struct SplitMix(u64);

impl SplitMix {
    #[inline(always)]
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

#[cfg(target_os = "linux")]
fn advise(ptr: *mut u8, bytes: usize, huge: bool) {
    const PAGE: usize = 4096;
    let start = (ptr as usize).div_ceil(PAGE) * PAGE;
    let end = (ptr as usize + bytes) / PAGE * PAGE;
    let advice = if huge {
        libc::MADV_HUGEPAGE
    } else {
        libc::MADV_NOHUGEPAGE
    };
    // SAFETY: [start, end) lies inside the live allocation; madvise only
    // changes the paging policy for that range.
    let rc = unsafe { libc::madvise(start as *mut libc::c_void, end - start, advice) };
    assert_eq!(rc, 0, "madvise failed");
}

#[cfg(not(target_os = "linux"))]
fn advise(_ptr: *mut u8, _bytes: usize, _huge: bool) {
    panic!("chase needs Linux madvise");
}

fn anon_huge_pages_bytes() -> Option<u64> {
    let text = fs::read_to_string("/proc/self/smaps_rollup").ok()?;
    let line = text.lines().find(|l| l.starts_with("AnonHugePages:"))?;
    Some(line.split_whitespace().nth(1)?.parse::<u64>().ok()? * 1024)
}

fn host() -> Value {
    let cpu = fs::read_to_string("/proc/cpuinfo").ok().and_then(|t| {
        t.lines()
            .find(|l| l.starts_with("model name"))
            .map(|l| l.split(':').nth(1).unwrap_or("").trim().to_string())
    });
    let thp = fs::read_to_string("/sys/kernel/mm/transparent_hugepage/enabled")
        .ok()
        .map(|s| s.trim().to_string());
    json!({"cpu_model": cpu, "thp_enabled": thp,
           "available_parallelism": std::thread::available_parallelism().map(|n| n.get()).ok()})
}

/// One 64-byte node per cache line; word 0 of node i holds the next node.
/// Sattolo's shuffle makes the successor map a single cycle over all nodes,
/// so a chase from node 0 visits every line before repeating.
fn chase(o: &Opts) -> Value {
    let bytes = o.num("bytes", 0) as usize;
    let pages = o.one("pages");
    let huge = match pages.as_str() {
        "4k" => false,
        "2m" => true,
        _ => panic!("--pages is 4k or 2m"),
    };
    let accesses = o.num("accesses", 1 << 24) as usize;
    let reps = o.num("reps", 7) as usize;
    let seed = o.num("seed", 20261017);
    let nodes = bytes / LINE;
    assert!(
        nodes.is_power_of_two() && nodes >= 2,
        "--bytes must be a power of two >= 128"
    );

    let words = nodes * WORDS_PER_LINE;
    let mut array: Vec<u64> = Vec::with_capacity(words);
    advise(array.as_mut_ptr() as *mut u8, words * 8, huge);
    let t_fill = Instant::now();
    array.resize(words, 0);
    for i in 0..nodes {
        array[i * WORDS_PER_LINE] = i as u64;
    }
    let mut rng = SplitMix(seed);
    for i in (1..nodes).rev() {
        let j = (rng.next() % i as u64) as usize;
        array.swap(i * WORDS_PER_LINE, j * WORDS_PER_LINE);
    }
    let fill_ms = t_fill.elapsed().as_secs_f64() * 1e3;
    let huge_bytes = anon_huge_pages_bytes();

    let mask = (nodes - 1) as u64;
    let run_chase = |n: usize, start: u64| {
        let mut p = start;
        for _ in 0..n {
            p = array[p as usize * WORDS_PER_LINE];
        }
        p
    };
    let run_gather = |n: usize, s: u64| {
        let mut r = SplitMix(s);
        let mut acc = 0u64;
        for _ in 0..n {
            acc = acc.wrapping_add(array[(r.next() & mask) as usize * WORDS_PER_LINE]);
        }
        acc
    };
    // Warm up: page faults are already taken by the fill; this warms the TLB
    // and the branch predictor without timing.
    let mut sink = run_chase(accesses / 4, 0) ^ run_gather(accesses / 4, seed ^ 1);
    let (mut chase_ns, mut gather_ns) = (Vec::new(), Vec::new());
    let mut cursor = sink & mask;
    for rep in 0..reps {
        let t = Instant::now();
        cursor = run_chase(accesses, cursor);
        chase_ns.push(t.elapsed().as_nanos() as f64 / accesses as f64);
        let t = Instant::now();
        sink ^= run_gather(accesses, seed.wrapping_add(rep as u64 + 2));
        gather_ns.push(t.elapsed().as_nanos() as f64 / accesses as f64);
    }
    sink ^= cursor;
    json!({
        "kind": "ic_probe_latency_chase",
        "schema": 1,
        "bytes": bytes, "pages": pages, "nodes": nodes, "line_bytes": LINE,
        "accesses_per_rep": accesses, "reps": reps, "seed": seed,
        "fill_ms": fill_ms,
        "anon_huge_pages_bytes": huge_bytes,
        "huge_fraction": huge_bytes.map(|h| h as f64 / bytes as f64),
        "chase_ns_per_access": chase_ns,
        "gather_ns_per_access": gather_ns,
        "chase_ns_median": median(&chase_ns),
        "gather_ns_median": median(&gather_ns),
        "sink": sink,
        "host": host(),
    })
}

// ---------------------------------------------------------------- geomfit

fn median(xs: &[f64]) -> f64 {
    assert!(!xs.is_empty());
    let mut v = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn read_jsonl(path: &Path) -> Vec<Value> {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("{}: {e}", path.display())))
        .collect()
}

/// Per-target rows of an IC `targets_result.jsonl`, keyed by fixture index.
struct Target {
    index: u64,
    probes: f64,
    target_ms: f64,
    recovered: Option<u64>,
}

fn targets(path: &Path) -> Vec<Target> {
    read_jsonl(path)
        .into_iter()
        .filter(|r| r["exit_code"] == 0 && r["probes"].is_u64())
        .map(|r| Target {
            index: r["fixture_index"].as_u64().unwrap(),
            probes: r["probes"].as_u64().unwrap() as f64,
            target_ms: r["target_ms"].as_f64().unwrap(),
            recovered: r["recovered_scalar"].as_u64(),
        })
        .collect()
}

/// Distribution summary of probe counts against Geometric(1/mean) on
/// {1, 2, ...}: CV, median/mean, max/mean, and the Kolmogorov distance.
fn probe_stats(probes: &[f64]) -> Value {
    let n = probes.len() as f64;
    let mean = probes.iter().sum::<f64>() / n;
    let var = probes.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
    let mut sorted = probes.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let p = 1.0 / mean;
    let cdf = |k: f64| 1.0 - (1.0 - p).powf(k);
    let mut ks = 0.0f64;
    for (i, &x) in sorted.iter().enumerate() {
        let f = cdf(x);
        ks = ks
            .max((f - i as f64 / n).abs())
            .max(((i + 1) as f64 / n - f).abs());
    }
    let max = *sorted.last().unwrap();
    json!({
        "n": probes.len(),
        "mean": mean,
        "cv": var.sqrt() / mean,
        "median_over_mean": median(probes) / mean,
        "max_over_mean": max / mean,
        "max_over_mean_bound_ln_n_plus_1": n.ln() + 1.0,
        "ks_distance_vs_geometric": ks,
        "ks_5pct_reference_1p36_over_sqrt_n": 1.36 / n.sqrt(),
    })
}

fn geomfit(o: &Opts) -> Value {
    let labels = o.all("label");
    let files = o.all("targets");
    assert_eq!(labels.len(), files.len(), "one --label per --targets");
    let rows: Vec<Value> = labels
        .iter()
        .zip(&files)
        .map(|(label, file)| {
            let t = targets(Path::new(file));
            let probes: Vec<f64> = t.iter().map(|t| t.probes).collect();
            let mut row = probe_stats(&probes);
            row["label"] = json!(label);
            row["targets"] = json!(file);
            row
        })
        .collect();
    json!({"kind": "ic_probe_latency_geomfit", "schema": 1, "bases": rows})
}

// ---------------------------------------------------------------- claims

fn read_json(path: &Path) -> Value {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The summary record a producer run printed (last JSON line of stdout).
fn summary(dir: &Path) -> Value {
    read_jsonl(&dir.join("stdout.jsonl"))
        .into_iter()
        .rev()
        .find(|r| r["kind"] == "compact_orbit_dlp_summary")
        .unwrap_or_else(|| panic!("{}: no summary record", dir.display()))
}

struct IcRun {
    dir: PathBuf,
    ns_per_probe: f64,
    targets: Vec<Target>,
    chase: Value,
    summary: Value,
}

fn ic_run(dir: &str) -> IcRun {
    let dir = PathBuf::from(dir);
    let targets = targets(&dir.join("targets_result.jsonl"));
    assert!(!targets.is_empty(), "{}: no solved targets", dir.display());
    let per: Vec<f64> = targets
        .iter()
        .map(|t| t.target_ms * 1e6 / t.probes)
        .collect();
    IcRun {
        ns_per_probe: median(&per),
        chase: read_json(&dir.join("chase.json")),
        summary: summary(&dir),
        targets,
        dir,
    }
}

fn verdict(pass: bool) -> &'static str {
    if pass {
        "PASS"
    } else {
        "FAIL"
    }
}

fn claims(o: &Opts) -> Value {
    let stage0 = read_json(Path::new(&o.one("stage0")));
    let archived = targets(Path::new(&o.one("archived")));
    let answers: HashMap<u64, u64> = read_jsonl(Path::new(&o.one("answers")))
        .into_iter()
        .map(|r| {
            let s = r["scalar"].as_str().unwrap().parse::<u64>().unwrap();
            (r["index"].as_u64().unwrap(), s)
        })
        .collect();
    let runs: Vec<(&str, IcRun)> = ["a1", "b1", "a2", "b2"]
        .iter()
        .map(|k| (*k, ic_run(&o.one(k))))
        .collect();
    let run = |k: &str| &runs.iter().find(|(n, _)| *n == k).unwrap().1;
    let mut out = Vec::new();

    // G: the same workload. Every on-host target recovers the published
    // scalar, and its probe count equals the archived K1798 run's.
    let arch: HashMap<u64, f64> = archived.iter().map(|t| (t.index, t.probes)).collect();
    let mut scalar_ok = true;
    let mut probes_equal = true;
    for (_, r) in &runs {
        scalar_ok &= r.targets.len() == answers.len();
        for t in &r.targets {
            scalar_ok &= t.recovered == answers.get(&t.index).copied();
            probes_equal &= arch.get(&t.index) == Some(&t.probes);
        }
    }
    out.push(json!({"id": "G.scalars", "statement": "every target in all four IC runs recovers its published scalar",
        "verdict": verdict(scalar_ok)}));
    out.push(json!({"id": "G.probes", "statement": "per-target probe counts equal the archived K1798 dense panel",
        "verdict": verdict(probes_equal)}));

    // H1: geometric probe counts, per archived base and on host.
    let mut bases: Vec<Value> = stage0["bases"].as_array().unwrap().clone();
    let host_probes: Vec<f64> = run("a1").targets.iter().map(|t| t.probes).collect();
    let mut host = probe_stats(&host_probes);
    host["label"] = json!("k1798_dense_on_host_a1");
    bases.push(host);
    for b in &bases {
        let cv = b["cv"].as_f64().unwrap();
        let mm = b["median_over_mean"].as_f64().unwrap();
        let support = (0.9..=1.1).contains(&cv) && (0.6..=0.8).contains(&mm);
        let falsified = !(0.6..=1.5).contains(&cv);
        out.push(json!({"id": format!("H1.{}", b["label"].as_str().unwrap()),
            "statement": "probe counts are geometric: CV in [0.9, 1.1] and median/mean in [0.6, 0.8]",
            "cv": cv, "median_over_mean": mm,
            "verdict": if support { "SUPPORT" } else if falsified { "FALSIFIED" } else { "NOT_SUPPORTED" }}));
    }

    // L1: IC ns/probe on 4 KiB pages against the 8 GiB pointer chase on
    // 4 KiB pages, interleaved in the same run directories.
    let a_ns = median(&[run("a1").ns_per_probe, run("a2").ns_per_probe]);
    let a_chase = median(&[
        run("a1").chase["chase_ns_median"].as_f64().unwrap(),
        run("a2").chase["chase_ns_median"].as_f64().unwrap(),
    ]);
    let a_gather = median(&[
        run("a1").chase["gather_ns_median"].as_f64().unwrap(),
        run("a2").chase["gather_ns_median"].as_f64().unwrap(),
    ]);
    let l1 = a_ns / a_chase;
    out.push(json!({"id": "L1", "statement": "IC ns/probe (4 KiB pages) / dependent-chase ns/access over the table bytes",
        "ic_ns_per_probe": a_ns, "chase_ns": a_chase, "value": l1,
        "rule": "SUPPORT if in [1.0, 2.0]; FALSIFIED if < 0.5",
        "verdict": if (1.0..=2.0).contains(&l1) { "SUPPORT" } else if l1 < 0.5 { "FALSIFIED" } else { "NOT_SUPPORTED" }}));
    out.push(
        json!({"id": "L1g", "statement": "IC ns/probe / independent-gather ns/access (descriptive)",
        "gather_ns": a_gather, "value": a_ns / a_gather, "verdict": "DESCRIPTIVE"}),
    );

    // T1: 2 MiB pages, ABAB. r_k = B_k / A_k, drift d = |A2/A1 - 1|.
    let r1 = run("b1").ns_per_probe / run("a1").ns_per_probe;
    let r2 = run("b2").ns_per_probe / run("a2").ns_per_probe;
    let drift = (run("a2").ns_per_probe / run("a1").ns_per_probe - 1.0).abs();
    let table = |r: &IcRun| r.summary["root_table_allocated_bytes"].as_f64().unwrap();
    let huge_frac =
        |r: &IcRun| r.summary["anon_huge_pages_bytes"].as_f64().unwrap_or(0.0) / table(r);
    let applied = huge_frac(run("b1")) >= 0.9
        && huge_frac(run("b2")) >= 0.9
        && huge_frac(run("a1")) < 0.1
        && huge_frac(run("a2")) < 0.1;
    let chase_ratio = median(&[
        run("b1").chase["chase_ns_median"].as_f64().unwrap(),
        run("b2").chase["chase_ns_median"].as_f64().unwrap(),
    ]) / a_chase;
    let null_sees = chase_ratio <= 0.9;
    out.push(json!({"id": "T1.applied", "statement": "huge pages back >= 90% of the table in B runs and < 10% in A runs",
        "b_fraction": [huge_frac(run("b1")), huge_frac(run("b2"))],
        "a_fraction": [huge_frac(run("a1")), huge_frac(run("a2"))], "verdict": verdict(applied)}));
    out.push(json!({"id": "T1.null", "statement": "positive control: the 8 GiB chase on 2 MiB pages / on 4 KiB pages <= 0.9",
        "value": chase_ratio, "verdict": verdict(null_sees)}));
    let worst = r1.max(r2);
    let t1 = if !applied || !null_sees {
        "INSTRUMENT_CANNOT_SEE"
    } else if worst <= 0.67 {
        "SUPPORT"
    } else if (r1 - 1.0).abs() < 0.10 && (r2 - 1.0).abs() < 0.10 {
        "FALSIFIED"
    } else {
        "NOT_SUPPORTED"
    };
    out.push(json!({"id": "T1", "statement": "IC ns/probe on 2 MiB pages / on 4 KiB pages, max(r1, r2)",
        "r1": r1, "r2": r2, "drift": drift, "value": worst,
        "rule": "SUPPORT if <= 0.67; FALSIFIED if both |r - 1| < 0.10; gated on T1.applied and T1.null",
        "verdict": t1}));

    // R: the rho compute control, ns per walk step before and after.
    let rho: Vec<f64> = o
        .one("rho")
        .split(',')
        .map(|d| {
            let r = read_json(&Path::new(d).join("receipt.json"));
            r["online_phase_ns"]["walk_and_collision"].as_f64().unwrap()
                / r["outcome"]["walk_steps"].as_f64().unwrap()
        })
        .collect();
    let rho_spread =
        rho.iter().cloned().fold(f64::MIN, f64::max) / rho.iter().cloned().fold(f64::MAX, f64::min);
    out.push(json!({"id": "R", "statement": "tuned rho ns/step across the session (compute control, max/min)",
        "ns_per_step": rho, "value": rho_spread, "verdict": "DESCRIPTIVE"}));
    out.push(json!({"id": "D", "statement": "IC/rho online per unit on this host: ns/probe / ns/step, raw and over the chase latency",
        "ic_ns_per_probe_4k": a_ns, "rho_ns_per_step": median(&rho),
        "ratio_raw": a_ns / median(&rho), "ratio_over_chase": a_ns / median(&rho) / a_chase,
        "verdict": "DESCRIPTIVE"}));

    json!({"kind": "ic_probe_latency_claims", "schema": 1,
           "runs": runs.iter().map(|(k, r)| json!({"arm": k, "dir": r.dir, "ns_per_probe": r.ns_per_probe,
                "targets": r.targets.len()})).collect::<Vec<_>>(),
           "claims": out})
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geometric_sample_reads_geometric() {
        let mut r = SplitMix(7);
        let p: f64 = 1e-4;
        let xs: Vec<f64> = (0..4000)
            .map(|_| {
                let u = (r.next() >> 11) as f64 / (1u64 << 53) as f64;
                ((1.0 - u).ln() / (1.0 - p).ln()).ceil().max(1.0)
            })
            .collect();
        let s = probe_stats(&xs);
        assert!((s["cv"].as_f64().unwrap() - 1.0).abs() < 0.06);
        assert!((s["median_over_mean"].as_f64().unwrap() - 2f64.ln()).abs() < 0.04);
        assert!(s["ks_distance_vs_geometric"].as_f64().unwrap() < 0.03);
    }

    #[test]
    fn constant_sample_is_not_geometric() {
        let s = probe_stats(&[100.0; 64]);
        assert_eq!(s["cv"].as_f64().unwrap(), 0.0);
        assert!(s["ks_distance_vs_geometric"].as_f64().unwrap() > 0.3);
    }

    #[test]
    fn sizes_parse() {
        assert_eq!(parse_size("8Gi"), 8 << 30);
        assert_eq!(parse_size("16Mi"), 16 << 20);
        assert_eq!(parse_size("4096"), 4096);
    }
}
