//! Weak-isogeny-class coverage: GHS magic numbers over F_{2^N} and Smith's
//! genus-3 kernel over F_p, per curve and per isogeny class.
//!
//! ```text
//! cargo run --release --example weak_class_coverage -- \
//!     --out research/weak_isogeny_classes_20261013/data \
//!     --ghs 8,9,10,12,14,15,16,18,20,21,22,24 --walk-max 20 --smith 3,5,7,11
//! ```
//! Writes `ghs_N<N>.json`, `smith_chebotarev.json`, `smith_p<p>.json` and
//! `run_manifest.json` (arguments, timings, git revision).

use crypto_lib::cryptanalysis::{ghs_class_coverage, smith_kernel_coverage};
use std::{fs, path::PathBuf, time::Instant};

fn list(s: &str) -> Vec<u32> {
    s.split(',')
        .filter(|x| !x.is_empty())
        .map(|x| x.parse().expect("integer list"))
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .map(|i| args[i + 1].clone())
    };
    let out = PathBuf::from(get("--out").expect("--out DIR"));
    let ghs = list(&get("--ghs").unwrap_or_default());
    let walk_max: u32 = get("--walk-max").map(|s| s.parse().unwrap()).unwrap_or(0);
    let smith = list(&get("--smith").unwrap_or_default());
    fs::create_dir_all(&out).unwrap();
    let mut timings = Vec::new();
    for &n in &ghs {
        let t0 = Instant::now();
        let rep = ghs_class_coverage::field_report(n, n <= walk_max);
        let secs = t0.elapsed().as_secs_f64();
        fs::write(
            out.join(format!("ghs_N{n}.json")),
            serde_json::to_string_pretty(&rep).unwrap(),
        )
        .unwrap();
        eprintln!("ghs N={n}: {secs:.1}s (walk3={})", rep.walk3_built);
        timings.push(serde_json::json!({"run": format!("ghs_N{n}"), "seconds": secs}));
    }
    if !smith.is_empty() {
        let cheb = smith_kernel_coverage::chebotarev_report();
        fs::write(
            out.join("smith_chebotarev.json"),
            serde_json::to_string_pretty(&cheb).unwrap(),
        )
        .unwrap();
        for &p in &smith {
            let t0 = Instant::now();
            let rep = smith_kernel_coverage::prime_report(p, &cheb);
            let secs = t0.elapsed().as_secs_f64();
            fs::write(
                out.join(format!("smith_p{p}.json")),
                serde_json::to_string_pretty(&rep).unwrap(),
            )
            .unwrap();
            eprintln!("smith p={p}: {secs:.1}s");
            timings.push(serde_json::json!({"run": format!("smith_p{p}"), "seconds": secs}));
        }
    }
    let rev = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default();
    let manifest = serde_json::json!({
        "command": args, "git_head": rev, "timings": timings,
        "threads": std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1),
    });
    let name = if ghs.is_empty() {
        "run_manifest_smith.json"
    } else if smith.is_empty() {
        "run_manifest_ghs.json"
    } else {
        "run_manifest.json"
    };
    fs::write(
        out.join(name),
        serde_json::to_string_pretty(&manifest).unwrap(),
    )
    .unwrap();
}
