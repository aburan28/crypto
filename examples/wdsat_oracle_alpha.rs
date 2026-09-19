//! **The F₂ SAT oracle's exponent.** Protocol X7 in
//! `RESEARCH_ECC2K130_BEST.md`.
//!
//! For `m = 3` point decomposition over the factor base
//! `V = ⟨1, z, …, z^{l−1}⟩ ⊂ F_{2^n}`, the per-call cost of an oracle
//! is `Q = 2^{α l}` up to a constant, and the route-target inequality
//! says an oracle moves the index-calculus exponent only if `α < 2`;
//! rho parity at `n = 131` needs `α ≈ 0.38`. Pair enumeration sits at
//! `α = 2`. This runner measures `α` for WDSat, the published solver
//! for exactly this Weil-descended symmetrised `S₄` encoding, on a
//! ladder in `l`, against the in-repo pairs-and-solve oracle
//! (`semaev_decomp::decompose`, `O(2^{2l} l)` field operations) on the
//! same targets, same host, same run.
//!
//! Rejection instances (no decomposition) carry the fit: an attack
//! spends most calls rejecting, and an exhaustive answer has no
//! trajectory luck. Satisfiable instances are reported beside them.
//! Every WDSat answer is checked: a SAT witness must zero `f₃` and the
//! two oracles must agree on every target.
//!
//! ```bash
//! cargo run --release --example wdsat_oracle_alpha -- \
//!   --solver /tmp/wdsat-wide/wdsat_solver --max-l 9 --targets 6 \
//!   --out experiments/ecc2k130_wdsat_alpha_20260919
//! ```

use crypto_lib::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_semaev_s4::{wdsat_anf, weil_descend_s4};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{find_irreducible_sparse, points_with_x};
use crypto_lib::cryptanalysis::semaev_decomp::{decompose, eval_f3, Gf2};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::env;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

const TARGET_SEED: u64 = 0x5A7_0001;
/// `n ≈ 3l`, so a uniform target decomposes with probability near
/// `1/3!` on every rung and both classes occur.
const LADDER: &[(u32, u32)] = &[(15, 5), (19, 6), (21, 7), (24, 8), (27, 9), (30, 10), (33, 11)];
const M: u32 = 3;
/// Route-target scale: `l` at `n = 131`, `m = 3`, and the `α` parity needs.
const L_131: f64 = 45.0;
const ALPHA_PARITY: f64 = 0.38;

struct Args {
    solver: PathBuf,
    out: Option<PathBuf>,
    max_l: u32,
    targets: usize,
    timeout: Duration,
    jobs: usize,
}

fn parse_args() -> Args {
    let args: Vec<String> = env::args().collect();
    let value = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1).cloned())
    };
    Args {
        solver: PathBuf::from(value("--solver").expect("--solver PATH to wdsat_solver")),
        out: value("--out").map(PathBuf::from),
        max_l: value("--max-l").map_or(9, |v| v.parse().expect("--max-l")),
        targets: value("--targets").map_or(6, |v| v.parse().expect("--targets")),
        timeout: Duration::from_secs(value("--timeout").map_or(3600, |v| v.parse().expect("--timeout"))),
        jobs: value("--jobs").map_or(6, |v| v.parse().expect("--jobs")),
    }
}

fn sha256_file(path: &Path) -> String {
    let out = Command::new("sha256sum").arg(path).output().expect("sha256sum");
    String::from_utf8_lossy(&out.stdout)
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string()
}

/// A uniform nonzero abscissa that lifts to `K_0 : y² + xy = x³ + 1`.
fn draw_target(rng: &mut StdRng, n: u32, irr: &IrreduciblePoly) -> F2mElement {
    let curve = BinaryCurve {
        m: n,
        irreducible: irr.clone(),
        a: F2mElement::zero(n),
        b: F2mElement::one(n),
        generator: BinaryPoint::Infinity,
        order: BigUint::zero(),
        cofactor: BigUint::one(),
    };
    loop {
        let raw: u64 = rng.gen_range(1..(1u64 << n));
        let x = F2mElement::from_biguint(&BigUint::from(raw), n);
        if !points_with_x(&curve, &x).is_empty() {
            return x;
        }
    }
}

/// Nanoseconds per `F_{2^n}` multiplication on this host, a dependent
/// chain so nothing overlaps.
fn mul_ns(gf: &Gf2, n: u32) -> f64 {
    let reps = 1u64 << 22;
    let mut a = 0x9E37_79B9_7F4A_7C15u64 & ((1u64 << n) - 1) | 1;
    let b = 0x2545_F491_4F6C_DD1Du64 & ((1u64 << n) - 1) | 3;
    let t = Instant::now();
    for _ in 0..reps {
        a = gf.mul(a, b) | 1;
    }
    let ns = t.elapsed().as_nanos() as f64 / reps as f64;
    std::hint::black_box(a);
    ns
}

struct SolverRun {
    status: String,
    conflicts: Option<u64>,
    wall_s: f64,
    assignment: Option<String>,
    stdout_tail: String,
}

fn run_wdsat(solver: &Path, anf: &Path, n: u32, l: u32, timeout: Duration) -> SolverRun {
    let started = Instant::now();
    let mut child = Command::new(solver)
        .args(["-i", anf.to_str().unwrap(), "-n", &n.to_string(), "-l", &l.to_string(), "-m", &M.to_string(), "-b"])
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn wdsat");
    let mut stdout = child.stdout.take().unwrap();
    let reader = std::thread::spawn(move || {
        let mut s = String::new();
        stdout.read_to_string(&mut s).ok();
        s
    });
    let mut timed_out = false;
    loop {
        match child.try_wait().expect("wait") {
            Some(_) => break,
            None if started.elapsed() > timeout => {
                child.kill().ok();
                child.wait().ok();
                timed_out = true;
                break;
            }
            None => std::thread::sleep(Duration::from_millis(5)),
        }
    }
    let wall_s = started.elapsed().as_secs_f64();
    let out = reader.join().unwrap_or_default();
    let lines: Vec<&str> = out.lines().map(str::trim).filter(|s| !s.is_empty()).collect();
    if timed_out {
        return SolverRun { status: "TIMEOUT".into(), conflicts: None, wall_s, assignment: None, stdout_tail: out.chars().rev().take(200).collect::<String>().chars().rev().collect() };
    }
    let conflicts = lines.last().and_then(|s| s.parse::<u64>().ok());
    // WDSat prints `UNSAT` or the satisfying assignment as one line of
    // `9l − 3` bits, then the conflict count; there is no `SAT` line.
    let n_vars = (9 * l - 3) as usize;
    let assignment = lines
        .iter()
        .find(|s| s.len() == n_vars && s.bytes().all(|b| b == b'0' || b == b'1'))
        .map(|s| s.to_string());
    let status = if lines.iter().any(|s| *s == "UNSAT") {
        "UNSAT"
    } else if assignment.is_some() {
        "SAT"
    } else {
        "UNKNOWN"
    };
    SolverRun { status: status.into(), conflicts, wall_s, assignment, stdout_tail: lines.iter().rev().take(3).rev().map(|s| s.to_string()).collect::<Vec<_>>().join(" | ") }
}

fn bits_to_u64(bits: &[u8]) -> u64 {
    bits.iter().enumerate().fold(0u64, |acc, (j, &b)| acc | (u64::from(b == b'1') << j))
}

fn median(xs: &mut Vec<f64>) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let k = xs.len();
    Some(if k % 2 == 1 { xs[k / 2] } else { (xs[k / 2 - 1] + xs[k / 2]) / 2.0 })
}

fn slope(points: &[(f64, f64)]) -> Option<f64> {
    if points.len() < 4 {
        return None;
    }
    let k = points.len() as f64;
    let mx = points.iter().map(|p| p.0).sum::<f64>() / k;
    let my = points.iter().map(|p| p.1).sum::<f64>() / k;
    let sxx: f64 = points.iter().map(|p| (p.0 - mx).powi(2)).sum();
    let sxy: f64 = points.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    Some(sxy / sxx)
}

fn intercept(points: &[(f64, f64)], s: f64) -> f64 {
    let k = points.len() as f64;
    let mx = points.iter().map(|p| p.0).sum::<f64>() / k;
    let my = points.iter().map(|p| p.1).sum::<f64>() / k;
    my - s * mx
}

fn main() {
    let args = parse_args();
    let solver_sha = sha256_file(&args.solver);
    let hostname = fs::read_to_string("/etc/hostname").unwrap_or_default().trim().to_string();
    let code_hash = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string());
    let anf_dir = args.out.as_ref().map(|d| d.join("anf"));
    if let Some(d) = &anf_dir {
        fs::create_dir_all(d).expect("anf dir");
    }
    let tmp_dir = env::temp_dir().join(format!("wdsat-alpha-{}", std::process::id()));
    fs::create_dir_all(&tmp_dir).expect("tmp dir");

    println!(
        "wdsat oracle alpha: solver={} sha256={} max_l={} targets/rung={} timeout={}s jobs={} seed={TARGET_SEED:#x}",
        args.solver.display(),
        &solver_sha[..16],
        args.max_l,
        args.targets,
        args.timeout.as_secs(),
        args.jobs
    );
    println!(
        "{:>3} {:>3} {:>3} {:>7} {:>9} {:>11} {:>11} {:>11} {:>8}  {}",
        "n", "l", "t", "class", "confl", "confl/trip", "wdsat s", "pairs s", "ratio", "check"
    );

    let mut rungs = Vec::new();
    let mut rows = Vec::new();
    let mut all_ok = true;
    let mut rng = StdRng::seed_from_u64(TARGET_SEED);
    for &(n, l) in LADDER {
        if l > args.max_l {
            continue;
        }
        let irr = find_irreducible_sparse(n).expect("irreducible");
        let gf = Gf2::new(&irr);
        let ns_per_mul = mul_ns(&gf, n);
        let span = 1u64 << l;
        let triples = span * (span + 1) * (span + 2) / 6;
        let pairs = span * (span - 1) / 2;

        // Targets, instances and the pairs-and-solve arm, sequentially.
        struct Prep {
            t: usize,
            x_r: F2mElement,
            xr: u64,
            anf: PathBuf,
            pairs_s: f64,
            pairs_witness: Option<[u64; 3]>,
        }
        let mut preps = Vec::new();
        for t in 0..args.targets {
            let x_r = draw_target(&mut rng, n, &irr);
            let xr = gf.from_element(&x_r);
            let system = weil_descend_s4(n, l, &irr, &F2mElement::one(n), &x_r);
            let anf = wdsat_anf(&system);
            let name = format!("n{n}l{l}_t{t}.anf");
            let path = anf_dir.as_ref().unwrap_or(&tmp_dir).join(&name);
            fs::write(&path, &anf).expect("write anf");
            let started = Instant::now();
            let witness = decompose(xr, l, &gf);
            let pairs_s = started.elapsed().as_secs_f64();
            if let Some(w) = witness {
                assert_eq!(eval_f3(w[0], w[1], w[2], xr, &gf), 0, "pairs-and-solve witness must zero f3");
            }
            preps.push(Prep { t, x_r, xr, anf: path, pairs_s, pairs_witness: witness });
        }

        // WDSat arm, `jobs` processes at a time.
        let mut runs: Vec<Option<SolverRun>> = (0..preps.len()).map(|_| None).collect();
        for chunk in (0..preps.len()).collect::<Vec<_>>().chunks(args.jobs.max(1)) {
            let results: Vec<(usize, SolverRun)> = std::thread::scope(|s| {
                let handles: Vec<_> = chunk
                    .iter()
                    .map(|&i| {
                        let p = &preps[i];
                        let solver = args.solver.clone();
                        let timeout = args.timeout;
                        s.spawn(move || (i, run_wdsat(&solver, &p.anf, n, l, timeout)))
                    })
                    .collect();
                handles.into_iter().map(|h| h.join().unwrap()).collect()
            });
            for (i, r) in results {
                runs[i] = Some(r);
            }
        }

        let mut unsat_confl = Vec::new();
        let mut unsat_wall = Vec::new();
        let mut unsat_pairs = Vec::new();
        let mut sat_confl = Vec::new();
        let mut sat_wall = Vec::new();
        let mut sat_pairs = Vec::new();
        let mut rung_ok = true;
        let mut timeouts = 0usize;
        for (p, r) in preps.iter().zip(runs.into_iter()) {
            let r = r.unwrap();
            let expected = if p.pairs_witness.is_some() { "SAT" } else { "UNSAT" };
            let mut check = String::new();
            let mut ok = r.status == expected;
            if r.status == "TIMEOUT" {
                timeouts += 1;
                ok = false;
                check.push_str("timeout");
            }
            let mut wd_witness = None;
            if r.status == "SAT" {
                match &r.assignment {
                    Some(bits) if bits.len() >= 3 * l as usize => {
                        let b = bits.as_bytes();
                        let xs = [
                            bits_to_u64(&b[0..l as usize]),
                            bits_to_u64(&b[l as usize..2 * l as usize]),
                            bits_to_u64(&b[2 * l as usize..3 * l as usize]),
                        ];
                        let zero = eval_f3(xs[0], xs[1], xs[2], p.xr, &gf) == 0;
                        ok &= zero;
                        check.push_str(if zero { "witness zeroes f3" } else { "WITNESS FAILS f3" });
                        wd_witness = Some(xs);
                    }
                    _ => {
                        ok = false;
                        check.push_str("no assignment line");
                    }
                }
            } else if r.status == "UNSAT" {
                check.push_str(if ok { "agrees with exhaustive pairs" } else { "DISAGREES with pairs" });
            }
            if !ok {
                rung_ok = false;
                all_ok = false;
            }
            let ratio = r.wall_s / p.pairs_s.max(1e-9);
            println!(
                "{n:>3} {l:>3} {:>3} {:>7} {:>9} {:>11.3} {:>11.3} {:>11.5} {:>8.1}  {}",
                p.t,
                r.status,
                r.conflicts.map_or("-".to_string(), |c| c.to_string()),
                r.conflicts.map_or(f64::NAN, |c| c as f64 / triples as f64),
                r.wall_s,
                p.pairs_s,
                ratio,
                if check.is_empty() { r.stdout_tail.clone() } else { check.clone() }
            );
            if ok && r.status == "UNSAT" {
                unsat_confl.push(r.conflicts.unwrap() as f64);
                unsat_wall.push(r.wall_s);
                unsat_pairs.push(p.pairs_s);
            } else if ok && r.status == "SAT" {
                sat_confl.push(r.conflicts.unwrap() as f64);
                sat_wall.push(r.wall_s);
                sat_pairs.push(p.pairs_s);
            }
            rows.push(json!({
                "n": n, "l": l, "target_index": p.t,
                "target_x": p.x_r.to_biguint().to_string(),
                "anf": p.anf.file_name().unwrap().to_str(),
                "wdsat_status": r.status, "wdsat_conflicts": r.conflicts, "wdsat_wall_s": r.wall_s,
                "wdsat_witness": wd_witness,
                "pairs_status": expected, "pairs_wall_s": p.pairs_s, "pairs_witness": p.pairs_witness,
                "verified": ok,
            }));
        }
        let n_unsat = unsat_confl.len();
        let n_sat = sat_confl.len();
        let med_c = median(&mut unsat_confl.clone());
        let med_w = median(&mut unsat_wall.clone());
        let med_p = median(&mut unsat_pairs.clone());
        rungs.push(json!({
            "n": n, "l": l,
            "irreducible_low_terms": irr.low_terms,
            "factor_base_size": span, "sorted_triples": triples, "pairs": pairs,
            "ns_per_field_mul": ns_per_mul,
            "targets": preps.len(), "unsat": n_unsat, "sat": n_sat, "timeouts": timeouts,
            "complete": rung_ok && timeouts == 0,
            "unsat_median_conflicts": med_c,
            "unsat_median_conflicts_per_triple": med_c.map(|c| c / triples as f64),
            "unsat_median_wdsat_wall_s": med_w,
            "unsat_median_pairs_wall_s": med_p,
            "unsat_median_wdsat_mul_equivalents": med_w.map(|w| w * 1e9 / ns_per_mul),
            "unsat_median_pairs_mul_equivalents": med_p.map(|w| w * 1e9 / ns_per_mul),
            "sat_median_conflicts": median(&mut sat_confl.clone()),
            "sat_median_wdsat_wall_s": median(&mut sat_wall.clone()),
            "sat_median_pairs_wall_s": median(&mut sat_pairs.clone()),
        }));
    }

    // Fits on complete rungs with at least one rejection.
    let fit_rungs: Vec<&Value> = rungs
        .iter()
        .filter(|r| r["complete"] == true && r["unsat"].as_u64().unwrap_or(0) > 0)
        .collect();
    let pts = |key: &str| -> Vec<(f64, f64)> {
        fit_rungs
            .iter()
            .filter_map(|r| Some((r["l"].as_f64()?, r[key].as_f64()?.log2())))
            .collect()
    };
    let fit_of = |key: &str| -> Value {
        let p = pts(key);
        match slope(&p) {
            Some(s) => {
                let c = intercept(&p, s);
                json!({"alpha": s, "log2_intercept": c, "rungs": p.len(),
                       "log2_q_at_l45_extrapolated": s * L_131 + c})
            }
            None => json!({"alpha": null, "rungs": p.len()}),
        }
    };
    let fit = json!({
        "rule": "least squares of log2(median over rejection instances) on l over complete rungs; null below four rungs",
        "wdsat_conflicts": fit_of("unsat_median_conflicts"),
        "wdsat_wall": fit_of("unsat_median_wdsat_wall_s"),
        "wdsat_mul_equivalents": fit_of("unsat_median_wdsat_mul_equivalents"),
        "pairs_and_solve_wall": fit_of("unsat_median_pairs_wall_s"),
        "pairs_and_solve_mul_equivalents": fit_of("unsat_median_pairs_mul_equivalents"),
        "reference": {"enumeration_alpha": 2.0, "linear_oracle_alpha": 1.0, "parity_alpha_at_n131": ALPHA_PARITY, "l_at_n131": L_131},
    });
    let summary = json!({
        "protocol": "RESEARCH_ECC2K130_BEST.md X7",
        "class": "measurement",
        "unit": "per-call oracle cost: WDSat conflicts (hardware-free), wall seconds, and wall / measured ns-per-F_{2^n}-multiplication on this host",
        "solver": {"path": args.solver.to_str(), "sha256": solver_sha, "args": ["-n", "-l", "-m 3", "-b"]},
        "host": hostname, "code_hash": code_hash, "target_seed": TARGET_SEED,
        "targets_per_rung": args.targets, "timeout_s": args.timeout.as_secs(),
        "all_verified": all_ok,
        "fit": fit,
        "rungs": rungs,
        "rows": rows,
    });
    println!("{}", serde_json::to_string_pretty(&summary["fit"]).unwrap());
    println!("all_verified = {all_ok}");
    if let Some(dir) = &args.out {
        fs::create_dir_all(dir).expect("out dir");
        fs::write(dir.join("summary.json"), serde_json::to_vec_pretty(&summary).unwrap()).expect("summary");
    }
    fs::remove_dir_all(&tmp_dir).ok();
}
