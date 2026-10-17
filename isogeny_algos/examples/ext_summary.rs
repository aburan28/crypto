//! Summary of the external comparison (`external/run.sh`): reads every
//! `<dir>/round*/*.jsonl`, checks that all tools produced the same digests (products, squares,
//! x^p mod g, Phi_l) and the same group orders (SEA), and prints Markdown tables with the best
//! time per tool over all rounds and our time divided by each external tool's.
//!
//! `cargo run --release --example ext_summary -- results/p15-external > results/p15-external/SUMMARY.md`
//! (plus `--dat` to print `op p n ours flint pari` rows for the chart instead).
use isogeny_algos::json::Json;
use std::collections::{BTreeMap, BTreeSet};

fn tool_of(t: &str) -> &'static str {
    if t == "isogeny_algos" {
        "ours"
    } else if t.starts_with("flint") {
        "flint"
    } else if t.ends_with("-nosea") {
        "pari-nosea"
    } else {
        "pari"
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let dir = &args[0];
    let dat = args.iter().any(|a| a == "--dat");
    let mut recs: Vec<(String, Json)> = vec![];
    let mut rounds: Vec<_> = std::fs::read_dir(dir)
        .expect("dir")
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().starts_with("round"))
        .map(|e| e.path())
        .collect();
    rounds.sort();
    for r in &rounds {
        let mut files: Vec<_> = std::fs::read_dir(r)
            .unwrap()
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .collect();
        files.sort();
        for f in files {
            for line in std::fs::read_to_string(&f).unwrap().lines() {
                if let Ok(j) = Json::parse(line) {
                    recs.push((r.file_name().unwrap().to_string_lossy().into(), j));
                }
            }
        }
    }
    let s = |j: &Json, k: &str| j.get(k).and_then(|v| v.as_str()).unwrap_or("").to_string();
    let n = |j: &Json, k: &str| j.get(k).and_then(|v| v.as_f64()).unwrap_or(f64::NAN);
    let mut mismatches = 0usize;

    // polynomial workloads: (op, p, n) -> tool -> best_us; digests per key
    let mut best: BTreeMap<(String, String, u64), BTreeMap<&str, f64>> = BTreeMap::new();
    let mut digs: BTreeMap<(String, String, u64), BTreeSet<String>> = BTreeMap::new();
    let porder = |p: &str| match p {
        "p61" => 0,
        "p127" => 1,
        "p255" => 2,
        _ => 3,
    };
    for (_, j) in &recs {
        let op = s(j, "op");
        if !["mul", "sqr", "powx"].contains(&op.as_str()) {
            continue;
        }
        let key = (op, s(j, "p"), n(j, "n") as u64);
        let t = tool_of(&s(j, "tool"));
        let v = best
            .entry(key.clone())
            .or_default()
            .entry(t)
            .or_insert(f64::MAX);
        *v = v.min(n(j, "best_us"));
        digs.entry(key).or_default().insert(s(j, "digest"));
    }
    let mut keys: Vec<_> = best.keys().cloned().collect();
    keys.sort_by_key(|(op, p, nn)| (op.clone(), porder(p), *nn));
    if dat {
        for k in &keys {
            let b = &best[k];
            let g = |t: &str| b.get(t).copied().unwrap_or(f64::NAN);
            println!(
                "{} {} {} {:.2} {:.2} {:.2}",
                k.0,
                k.1,
                k.2,
                g("ours"),
                g("flint"),
                g("pari")
            );
        }
        return;
    }
    println!("# External comparison: isogeny_algos vs FLINT and PARI/GP\n");
    println!("Best time over all rounds per tool, in microseconds; ratio = ours / theirs (below 1: ours is faster). `same` = every tool and round printed the same digest.\n");
    for op in ["mul", "sqr", "powx"] {
        println!("## {op}\n");
        println!("| p | n | ours | FLINT | PARI | ours/FLINT | ours/PARI | outputs |");
        println!("|:--|--:|--:|--:|--:|--:|--:|:--|");
        for k in keys.iter().filter(|k| k.0 == op) {
            let b = &best[k];
            let g = |t: &str| b.get(t).copied().unwrap_or(f64::NAN);
            let same = digs[k].len() == 1;
            if !same {
                mismatches += 1;
            }
            println!(
                "| {} | {} | {:.1} | {:.1} | {:.1} | {:.2} | {:.2} | {} |",
                k.1,
                k.2,
                g("ours"),
                g("flint"),
                g("pari"),
                g("ours") / g("flint"),
                g("ours") / g("pari"),
                if same { "same" } else { "DIFFER" }
            );
        }
        println!();
    }

    // Phi_l mod p for all l <= max_ell
    println!("## Modular polynomials Phi_l mod p, all primes 3 <= l <= 89, p = 2^127 - 1\n");
    let mut phi_t: BTreeMap<&str, f64> = BTreeMap::new();
    let mut phi_d: BTreeMap<u64, BTreeSet<String>> = BTreeMap::new();
    for (_, j) in &recs {
        let op = s(j, "op");
        let t = tool_of(&s(j, "tool"));
        let label = match op.as_str() {
            "phi_all" => Some("ours: full Phi_l mod p (Hecke, compute_many)"),
            "phi_spec_all" => {
                Some("PARI: Phi_l(j, Y) mod p for one j (polmodular(l, 0, Mod(j, p)))")
            }
            "phi_full_all" => Some("PARI: full Phi_l over Z, reduced mod p (polmodular(l))"),
            _ => None,
        };
        if let Some(lb) = label {
            let v = phi_t.entry(lb).or_insert(f64::MAX);
            *v = v.min(n(j, "total_us"));
        }
        if op == "phi_digest" {
            let set = phi_d.entry(n(j, "ell") as u64).or_default();
            set.insert(s(j, "digest"));
            if t == "pari" {
                set.insert(s(j, "digest_full"));
            }
        }
    }
    println!("| computation | best total over rounds |");
    println!("|:--|--:|");
    for (lb, us) in &phi_t {
        println!("| {lb} | {:.0} ms |", us / 1e3);
    }
    let bad: Vec<_> = phi_d
        .iter()
        .filter(|(_, d)| d.len() != 1)
        .map(|(l, _)| *l)
        .collect();
    mismatches += bad.len();
    println!(
        "\nDigests Phi_l(987654321, 123456789) mod p: {} of {} levels agree across tools{}.\n",
        phi_d.len() - bad.len(),
        phi_d.len(),
        if bad.is_empty() {
            String::new()
        } else {
            format!(" (differ: {bad:?})")
        }
    );

    // SEA: per curve, median over rounds per tool; orders must agree
    for p in ["p61", "p127"] {
        println!("## SEA point counting, {p}\n");
        let mut per: BTreeMap<u64, BTreeMap<&str, Vec<f64>>> = BTreeMap::new();
        let mut ord: BTreeMap<u64, BTreeSet<String>> = BTreeMap::new();
        let mut pre = f64::MAX;
        for (_, j) in &recs {
            if s(j, "p") != p {
                continue;
            }
            match s(j, "op").as_str() {
                "sea" => {
                    let k = n(j, "k") as u64;
                    per.entry(k)
                        .or_default()
                        .entry(tool_of(&s(j, "tool")))
                        .or_default()
                        .push(n(j, "us"));
                    ord.entry(k).or_default().insert(s(j, "order"));
                }
                "sea_phi_precompute" => pre = pre.min(n(j, "total_us")),
                _ => {}
            }
        }
        let tools = ["ours", "pari", "pari-nosea"];
        let mut sums = [0.0f64; 3];
        let mut agree = 0;
        println!("| curve k | ours ms | PARI (seadata) ms | PARI (no seadata) ms | orders |");
        println!("|--:|--:|--:|--:|:--|");
        for (k, m) in &per {
            let mut row = vec![];
            for (i, t) in tools.iter().enumerate() {
                let v = m
                    .get(t)
                    .map(|v| v.iter().cloned().fold(f64::MAX, f64::min))
                    .unwrap_or(f64::NAN);
                sums[i] += v;
                row.push(v / 1e3);
            }
            let same = ord[k].len() == 1 && !ord[k].contains("none");
            agree += same as usize;
            println!(
                "| {k} | {:.1} | {:.1} | {:.1} | {} |",
                row[0],
                row[1],
                row[2],
                if same { "same" } else { "DIFFER" }
            );
        }
        mismatches += per.len() - agree;
        println!(
            "| **sum** | **{:.0}** | **{:.0}** | **{:.0}** | {agree}/{} agree |\n",
            sums[0] / 1e3,
            sums[1] / 1e3,
            sums[2] / 1e3,
            per.len()
        );
        println!(
            "Ours needs Phi_l for l up to {} computed once per field first: {:.0} ms (PARI with seadata reads precomputed polynomials from disk; without it, PARI computes what it needs per curve, inside the times above).\n",
            if p == "p61" { 61 } else { 89 },
            pre / 1e3
        );
    }
    println!("**Output mismatches: {mismatches}.**");
}
