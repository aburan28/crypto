//! XL against MutantXL closure on committed refutation-degree draws.
//!
//! Every refutation degree the ladder studies record is XL's: the first
//! degree whose plain Macaulay row space contains `1`.  This replays the
//! committed unsatisfiable draws (same subspace basis, same `x(R)`, same
//! field) and scans each degree twice, plain and saturated under
//! MutantXL's fixed-degree step ([`closure_refutation_scan`]).
//!
//! ```sh
//! cargo run --release --example dreg_closure -- --max-cols 70000 \
//!     research/dreg_m4_chain_20261006/runs/cell-5-2-8.jsonl ...
//! ```
//!
//! One JSON line per unsatisfiable draw on stdout.  The scan runs from
//! degree 3 up to the draw's committed XL degree (its resolved degree, or
//! the degree its `≥D` bound was measured to), stops early once the closure
//! refutes, and builds no degree whose matrix could have more than
//! `--max-cols` columns (`C(N,≤D)`): the dense basis costs about the cube
//! of that.  Degrees not built are unmeasured, never refutations.

use crypto_lib::binary_ecc::F2mElement;
use crypto_lib::cryptanalysis::koblitz_bench::{closure_refutation_scan, ClosureDegree};
use crypto_lib::cryptanalysis::koblitz_groebner::{build_decomposition_system, FieldStructure};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use num_bigint::BigUint;
use serde_json::Value;

fn monomials_up_to(n_vars: u64, d: u32) -> f64 {
    let (mut total, mut binom) = (0f64, 1f64);
    for i in 0..=d as u64 {
        if i > n_vars {
            break;
        }
        total += binom;
        binom = binom * (n_vars - i) as f64 / (i + 1) as f64;
    }
    total
}

fn degree_json(d: &ClosureDegree) -> String {
    format!(
        r#"{{"degree":{},"xl_refuted":{},"closure_refuted":{},"saturated":{},"xl_rank":{},"closure_rank":{},"columns":{},"closure_products":{},"closure_word_ops":{},"secs":{:.3}}}"#,
        d.degree,
        d.xl_refuted,
        d.closure_refuted,
        d.saturated,
        d.xl_rank,
        d.closure_rank,
        d.columns,
        d.closure_products,
        d.closure_word_ops,
        d.secs
    )
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let flag = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let max_cols: f64 = flag("--max-cols").map_or(70_000.0, |v| v.parse().expect("--max-cols"));
    let skip = ["--max-cols"];
    let files: Vec<&String> = args
        .iter()
        .enumerate()
        .filter(|(i, a)| {
            !skip.contains(&a.as_str()) && !(*i > 0 && skip.contains(&args[i - 1].as_str()))
        })
        .map(|(_, a)| a)
        .collect();
    for file in files {
        let text = std::fs::read_to_string(file).unwrap_or_else(|e| panic!("{file}: {e}"));
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let v: Value = serde_json::from_str(line).expect("json");
            if v.get("control").is_some() || v["outcome"]["kind"] == "satisfiable" {
                continue;
            }
            let m = v.get("m").and_then(Value::as_u64).unwrap_or(3) as usize;
            let n = v["n"].as_u64().expect("n") as u32;
            let ell = v["ell"].as_u64().expect("ell") as usize;
            let irr = find_irreducible_sparse(n).expect("field");
            let st = FieldStructure::new(n, &irr);
            let b = F2mElement::one(n);
            let basis: Vec<F2mElement> = v["v_basis"]
                .as_array()
                .expect("v_basis")
                .iter()
                .map(|w| F2mElement::from_biguint(&BigUint::from(w.as_u64().unwrap()), n))
                .collect();
            let x_r = F2mElement::from_biguint(&BigUint::from(v["x_r"].as_u64().unwrap()), n);
            assert_eq!(basis.len(), ell);
            let sys = build_decomposition_system(&basis, &x_r, &b, m, &st).expect("system");
            assert_eq!(
                sys.n_vars as u64,
                v["n_vars"].as_u64().unwrap(),
                "{file}: rebuilt shape"
            );
            // The committed XL degree: resolved at `D`, or `≥D` (built to D−1).
            let committed = v["outcome"]["degree"].as_u64().expect("committed degree") as u32;
            let xl_known = v["outcome"]["kind"] == "resolved";
            let ceiling = if xl_known { committed } else { committed - 1 };
            let top = (3..=ceiling)
                .take_while(|&d| monomials_up_to(sys.n_vars as u64, d) <= max_cols)
                .last()
                .unwrap_or(2);
            let mut scan = Vec::new();
            for d in 3..=top {
                let one = closure_refutation_scan(&sys.equations, sys.n_vars, d, d);
                let Some(step) = one.into_iter().next() else {
                    break;
                };
                // XL identity: the plain basis refutes exactly at the
                // committed degree, never below it.
                assert_eq!(
                    step.xl_refuted,
                    xl_known && d >= committed,
                    "{file} draw {}: XL at degree {d} disagrees with the committed outcome",
                    v["draw"]
                );
                let done = step.closure_refuted;
                scan.push(step);
                if done {
                    break;
                }
            }
            let first =
                |f: fn(&ClosureDegree) -> bool| scan.iter().find(|d| f(d)).map(|d| d.degree);
            let show = |x: Option<u32>| x.map_or("null".to_string(), |d| d.to_string());
            let measured_to = scan.last().map(|d| d.degree).unwrap_or(2);
            println!(
                r#"{{"source":"{file}","cell":"{}","m":{m},"n":{n},"ell":{ell},"draw":{},"n_vars":{},"committed_outcome":{},"xl_degree":{},"closure_degree":{},"measured_to":{measured_to},"max_cols":{max_cols},"scan":[{}]}}"#,
                v["cell"].as_str().unwrap_or(""),
                v["draw"],
                sys.n_vars,
                v["outcome"],
                show(first(|d| d.xl_refuted)),
                show(first(|d| d.closure_refuted)),
                scan.iter().map(degree_json).collect::<Vec<_>>().join(",")
            );
            eprintln!(
                "m={m} ({n},{ell}) draw {}: committed {} | XL {} closure {} (measured to {measured_to})",
                v["draw"],
                v["outcome"],
                show(first(|d| d.xl_refuted)),
                show(first(|d| d.closure_refuted))
            );
        }
    }
    if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
        if let Some(l) = status.lines().find(|l| l.starts_with("VmHWM:")) {
            eprintln!("peak {l}");
        }
    }
}
