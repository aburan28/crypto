//! Per-record comparison of two arms of `bench` records (`results/p14-bench/`): the minimum over
//! each arm's two repetitions, the geometric mean of base/new per group, and the largest changes
//! among records whose base time exceeds 5 ms.
//!
//! ```text
//! cargo run --release --example bench_compare -- results/p14-bench main round5-pre-bmss \
//!     sea phi bmss find csidh big kernel
//! ```
use isogeny_algos::json::Json;

fn load(dir: &str, g: &str, arm: &str, rep: u32) -> Vec<Json> {
    let path = format!("{dir}/{g}.{arm}.{rep}.jsonl");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("{path}: {e}"))
        .lines()
        .map(|l| Json::parse(l).expect("json line"))
        .collect()
}

/// The record's timing field: median_ns, ns or min_ns, else the first field ending in _ns.
fn tkey(r: &Json) -> Option<String> {
    for k in ["median_ns", "ns", "min_ns"] {
        if r.get(k).and_then(|v| v.as_f64()).is_some() {
            return Some(k.into());
        }
    }
    if let Json::Obj(kv) = r {
        for (k, v) in kv {
            if k.ends_with("_ns") && v.as_f64().is_some() {
                return Some(k.clone());
            }
        }
    }
    None
}

fn label(r: &Json) -> String {
    match r {
        Json::Obj(kv) => kv
            .iter()
            .filter(|(k, _)| k != "note" && k != "group")
            .filter_map(|(k, v)| v.as_str().map(|s| format!("{k}={s}")))
            .collect::<Vec<_>>()
            .join(" "),
        _ => String::new(),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (dir, base, new) = (&args[0], &args[1], &args[2]);
    for g in &args[3..] {
        let (b1, b2) = (load(dir, g, base, 1), load(dir, g, base, 2));
        let (n1, n2) = (load(dir, g, new, 1), load(dir, g, new, 2));
        assert!(
            b1.len() == b2.len() && b1.len() == n1.len() && b1.len() == n2.len(),
            "{g}: record counts differ"
        );
        let unver = n1
            .iter()
            .chain(&n2)
            .filter(|r| !matches!(r.get("verified"), Some(Json::Bool(true))))
            .count();
        // (base/new, base spread, base, new, label)
        let mut rows: Vec<(f64, f64, f64, f64, String)> = vec![];
        for i in 0..b1.len() {
            let Some(k) = tkey(&b1[i]) else { continue };
            let t = |r: &Json| r.get(&k).and_then(|v| v.as_f64()).unwrap_or(0.0);
            let (b, n) = (t(&b1[i]).min(t(&b2[i])), t(&n1[i]).min(t(&n2[i])));
            if b <= 0.0 || n <= 0.0 {
                continue;
            }
            let spread = t(&b1[i]).max(t(&b2[i])) / b;
            rows.push((b / n, spread, b, n, label(&b1[i])));
        }
        let gm = |v: &[&(f64, f64, f64, f64, String)]| {
            (v.iter().map(|r| r.0.ln()).sum::<f64>() / v.len() as f64).exp()
        };
        let all: Vec<_> = rows.iter().collect();
        let mut big: Vec<_> = rows.iter().filter(|r| r.2 > 5e6).collect();
        println!(
            "== {g}: {} records, new unverified={unver}, geomean {:.2}x; records >5ms: {}, geomean {:.2}x",
            rows.len(),
            gm(&all),
            big.len(),
            if big.is_empty() { f64::NAN } else { gm(&big) }
        );
        big.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        for r in big.iter().take(6) {
            println!(
                "   {:5.2}x (base spread {:.2})  {:9.2} -> {:9.2} ms  {}",
                r.0,
                r.1,
                r.2 / 1e6,
                r.3 / 1e6,
                r.4
            );
        }
        for r in big.iter().rev().take(3) {
            println!(
                "   low {:5.2}x (base spread {:.2})  {:9.2} -> {:9.2} ms  {}",
                r.0,
                r.1,
                r.2 / 1e6,
                r.3 / 1e6,
                r.4
            );
        }
    }
}
