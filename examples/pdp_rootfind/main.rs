//! `pdp_rootfind`: root-finding decomposition oracles for the PDP on both
//! tracks, measured against the frozen corpora and the committed F4 rows
//! (`research/pdp_root_finder_20261012/`).
//!
//! Binary (the PS-0 rungs, `research/pdp_oracle_speedup_20261011/baseline_v1/corpus/`):
//!
//! ```bash
//! cargo run --release --example pdp_rootfind -- bin-verify --l 8 --corpus corpus_l8.txt --arm rf-mitm
//! cargo run --release --example pdp_rootfind -- bin-count  --l 8 --corpus corpus_l8.txt --arm rf-mitm
//! cargo run --release --example pdp_rootfind -- bin-time   --l 11 --corpus corpus_l11.txt --arm rf-mitm --targets 16 --passes 3
//! ```
//!
//! Prime (Kummer cells rebuilt from committed PKM rows, or a fresh ladder):
//!
//! ```bash
//! cargo run --release --example pdp_rootfind -- prime-rows --rows runs/M3-t5.jsonl,runs/M3-t6.jsonl
//! cargo run --release --example pdp_rootfind -- prime-ladder --p 2013265921 --m 3 --t-min 4 --t-max 12 --arms pairs3,mitm3
//! ```
//!
//! A stage diagnostic (`AGENTS.md` §8): one decomposition per call, never an
//! end-to-end ECDLP.

// The PS harness's modules, shared by path: `pdp_speedup` uses every item,
// this driver only some.
#[allow(dead_code)]
#[path = "../pdp_speedup/arms.rs"]
mod arms;
#[allow(dead_code)]
#[path = "../pdp_speedup/common.rs"]
mod common;
#[allow(dead_code)]
#[path = "../pdp_speedup/corpus.rs"]
mod corpus;
#[allow(dead_code)]
#[path = "../pdp_speedup/group.rs"]
mod group;
#[allow(dead_code)]
#[path = "../pdp_speedup/kernel.rs"]
mod kernel;

mod analyze;
mod binrf;
mod fpset;
mod primerf;

use arms::Arm;
use common::*;
use corpus::*;
use crypto_lib::cryptanalysis::prime_fast::Pt;
use crypto_lib::cryptanalysis::semaev_decomp::decompose;
use group::{degenerate_target, ArtinSchreier};
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

fn append(out: &Option<String>, line: &str) {
    println!("{line}");
    if let Some(out) = out {
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(out)
            .unwrap();
        writeln!(f, "{line}").unwrap();
    }
}

// ── binary ──────────────────────────────────────────────────────────────

enum BinArm<'a> {
    Ps(Arm<'a>),
    Mitm(&'a Rung, binrf::Mitm<'a>),
    Pairs(&'a Rung, binrf::Pairs<'a>),
}

impl<'a> BinArm<'a> {
    fn build(name: &str, rung: &'a Rung, witnesses: bool) -> Self {
        match name {
            "rf-mitm" => {
                let solver = ArtinSchreier::new(&rung.gf);
                BinArm::Mitm(rung, binrf::Mitm::new(&rung.gf, &solver, rung.l, witnesses))
            }
            "rf-pairs" => {
                let solver = ArtinSchreier::new(&rung.gf);
                BinArm::Pairs(rung, binrf::Pairs::new(&rung.gf, &solver, rung.l))
            }
            other => BinArm::Ps(Arm::build(other, rung)),
        }
    }

    fn decide(&self, xr: u64, yr: u64) -> (bool, binrf::Ops) {
        match self {
            BinArm::Ps(a) => (a.decide(xr, yr), binrf::Ops::default()),
            BinArm::Mitm(r, m) => {
                if degenerate_target(&r.gf, r.l, xr) {
                    (decompose(xr, r.l, &r.gf).is_some(), binrf::Ops::default())
                } else {
                    let mut out = Vec::new();
                    m.run((xr, yr), false, &mut out)
                }
            }
            BinArm::Pairs(r, p) => {
                if degenerate_target(&r.gf, r.l, xr) {
                    (decompose(xr, r.l, &r.gf).is_some(), binrf::Ops::default())
                } else {
                    let mut out = Vec::new();
                    p.run((xr, yr), false, &mut out)
                }
            }
        }
    }

    fn all_roots(&self, xr: u64, yr: u64) -> Option<Vec<[u64; 3]>> {
        match self {
            BinArm::Ps(a) => a.all_roots(xr, yr),
            BinArm::Mitm(r, m) if !degenerate_target(&r.gf, r.l, xr) => {
                let mut out = Vec::new();
                m.run((xr, yr), true, &mut out);
                out.sort_unstable();
                out.dedup();
                Some(out)
            }
            BinArm::Pairs(r, p) if !degenerate_target(&r.gf, r.l, xr) => {
                let mut out = Vec::new();
                p.run((xr, yr), true, &mut out);
                out.sort_unstable();
                out.dedup();
                Some(out)
            }
            _ => None,
        }
    }

    fn table_bytes(&self) -> usize {
        match self {
            BinArm::Mitm(_, m) => m.table.bytes(),
            _ => 0,
        }
    }

    fn base_len(&self) -> usize {
        match self {
            BinArm::Mitm(_, m) => m.base.len(),
            BinArm::Pairs(_, p) => p.base.len(),
            _ => 0,
        }
    }
}

fn rung_of(args: &[String]) -> Rung {
    let l: u32 = arg(args, "--l").expect("--l");
    match arg::<u32>(args, "--n") {
        Some(n) => Rung::with(n, l),
        None => Rung::new(l),
    }
}

/// `bin-corpus --l L [--n N] --out F [--uniform U --planted P --threads T]`:
/// the PS-0 corpus generator on a rung with `n` chosen independently of `ℓ`.
fn cmd_bin_corpus(args: &[String]) {
    let rung = rung_of(args);
    let out = arg_s(args, "--out").expect("--out");
    let threads: usize = arg(args, "--threads").unwrap_or(3);
    let n_uniform: usize = arg(args, "--uniform").unwrap_or(N_UNIFORM);
    let n_planted: usize = arg(args, "--planted").unwrap_or(N_PLANTED);
    let t0 = Instant::now();
    let recs = build_records(&rung, threads, n_uniform, n_planted);
    let text = write_corpus(&rung, &recs);
    std::fs::write(&out, &text).unwrap();
    let found = recs.iter().filter(|r| !r.roots.is_empty()).count();
    println!(
        "{{\"cmd\":\"bin-corpus\",\"l\":{},\"n\":{},\"targets\":{},\"found\":{found},\"enumerate_secs\":{:.1},\"digest_verdicts\":\"{}\",\"digest_roots\":\"{}\"}}",
        rung.l, rung.n, recs.len(), t0.elapsed().as_secs_f64(), digest_verdicts(&recs), digest_roots(&recs)
    );
}

fn load(args: &[String]) -> (Rung, Loaded) {
    let l: u32 = arg(args, "--l").expect("--l");
    let path = arg_s(args, "--corpus").expect("--corpus");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let loaded = read_corpus(&text);
    let rung = rung_of(args);
    assert_eq!(loaded.l, l, "corpus is for a different rung");
    assert_eq!(loaded.modulus, rung.modulus(), "corpus modulus differs");
    (rung, loaded)
}

fn cmd_bin_verify(args: &[String]) {
    let (rung, loaded) = load(args);
    let arm_name = arg_s(args, "--arm").expect("--arm");
    let t0 = Instant::now();
    let arm = BinArm::build(&arm_name, &rung, true);
    let setup = t0.elapsed().as_secs_f64();
    let mut bad = 0usize;
    let mut replay = loaded.recs.clone();
    let mut roots_checked = 0usize;
    let t1 = Instant::now();
    for rec in replay.iter_mut() {
        let truth = !rec.roots.is_empty();
        let (got, _) = arm.decide(rec.t.xr, rec.t.yr);
        if got != truth {
            println!("FAIL target {} verdict {got} vs corpus {truth}", rec.t.idx);
            bad += 1;
        }
        if let Some(all) = arm.all_roots(rec.t.xr, rec.t.yr) {
            roots_checked += 1;
            if all != rec.roots {
                println!(
                    "FAIL target {} root set differs: arm {:?} corpus {:?}",
                    rec.t.idx, all, rec.roots
                );
                bad += 1;
            }
        }
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
    // Root-set digest from the arm's own root sets, where it lists them.
    let mut dr = String::from("-");
    if roots_checked == loaded.recs.len() {
        let mut own = loaded.recs.clone();
        for rec in own.iter_mut() {
            rec.roots = arm.all_roots(rec.t.xr, rec.t.yr).unwrap_or_default();
        }
        dr = digest_roots(&own);
        if dr != loaded.file_digest_roots {
            println!("FAIL root-set digest {dr} != corpus digest");
            bad += 1;
        }
    }
    println!(
        "{{\"cmd\":\"bin-verify\",\"arm\":\"{arm_name}\",\"l\":{},\"n\":{},\"targets\":{},\"root_sets_checked\":{roots_checked},\"setup_secs\":{setup:.3},\"secs\":{:.3},\"digest_verdicts\":\"{dv}\",\"digest_roots\":\"{dr}\",\"digest_verdicts_corpus\":\"{}\",\"digest_roots_corpus\":\"{}\",\"mismatches\":{bad},\"table_bytes\":{},\"base_points\":{}}}",
        rung.l,
        rung.n,
        loaded.recs.len(),
        t1.elapsed().as_secs_f64(),
        loaded.file_digest_verdicts,
        loaded.file_digest_roots,
        arm.table_bytes(),
        arm.base_len()
    );
    if bad != 0 {
        std::process::exit(1);
    }
}

fn cmd_bin_count(args: &[String]) {
    let (rung, loaded) = load(args);
    let arm_name = arg_s(args, "--arm").expect("--arm");
    let k: usize = arg(args, "--targets").unwrap_or(16);
    let arm = BinArm::build(&arm_name, &rung, false);
    let refuted: Vec<(u64, u64)> = loaded
        .recs
        .iter()
        .filter(|r| r.roots.is_empty())
        .map(|r| (r.t.xr, r.t.yr))
        .take(k)
        .collect();
    let mut tot = binrf::Ops::default();
    for &(xr, yr) in &refuted {
        let (found, ops) = arm.decide(xr, yr);
        assert!(!found);
        tot.add(ops);
    }
    let c = refuted.len() as f64;
    let build = match &arm {
        BinArm::Mitm(_, m) => m.table.build_ops,
        _ => binrf::Ops::default(),
    };
    println!(
        "{{\"cmd\":\"bin-count\",\"arm\":\"{arm_name}\",\"l\":{},\"n\":{},\"targets\":{},\"base_points\":{},\"per_target\":{{\"mul\":{:.1},\"sqr\":{:.1},\"inv\":{:.3},\"lookups\":{:.1},\"candidates\":{:.1}}},\"prepare\":{{\"mul\":{},\"sqr\":{},\"inv\":{},\"table_bytes\":{}}}}}",
        rung.l, rung.n, refuted.len(), arm.base_len(),
        tot.mul as f64 / c, tot.sqr as f64 / c, tot.inv as f64 / c, tot.lookups as f64 / c, tot.candidates as f64 / c,
        build.mul, build.sqr, build.inv, arm.table_bytes()
    );
}

fn cmd_bin_time(args: &[String]) {
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
    let ts = Instant::now();
    let arm = BinArm::build(&arm_name, &rung, false);
    let setup_secs = ts.elapsed().as_secs_f64();
    black_box(arm.decide(refuted[0].0, refuted[0].1));
    let t0 = Instant::now();
    let mut passes = 0usize;
    let mut tot = binrf::Ops::default();
    loop {
        for &(xr, yr) in &refuted {
            let (found, ops) = black_box(arm.decide(black_box(xr), black_box(yr)));
            assert!(!found);
            tot.add(ops);
        }
        passes += 1;
        if passes >= passes_min && t0.elapsed().as_secs_f64() >= min_secs {
            break;
        }
    }
    let secs = t0.elapsed().as_secs_f64();
    let calls = (refuted.len() * passes) as u64;
    let label = arg_s(args, "--label").unwrap_or_default();
    let line = format!(
        "{{\"cmd\":\"bin-time\",\"label\":\"{label}\",\"arm\":\"{arm_name}\",\"l\":{},\"n\":{},\"targets\":{},\"passes\":{passes},\"calls\":{calls},\"secs\":{secs:.6},\"setup_secs\":{setup_secs:.6},\"ns_per_target\":{:.1},\"mul_per_target\":{:.1},\"sqr_per_target\":{:.1},\"inv_per_target\":{:.3},\"lookups_per_target\":{:.1},\"candidates_per_target\":{:.1},\"table_bytes\":{},\"base_points\":{},\"kernel\":\"{}\"}}",
        rung.l, rung.n, refuted.len(),
        secs * 1e9 / calls as f64,
        tot.mul as f64 / calls as f64, tot.sqr as f64 / calls as f64, tot.inv as f64 / calls as f64,
        tot.lookups as f64 / calls as f64, tot.candidates as f64 / calls as f64,
        arm.table_bytes(), arm.base_len(), rung.gf.kernel_name()
    );
    append(&arg_s(args, "--out"), &line);
}

// ── prime ───────────────────────────────────────────────────────────────

struct PrimeArms<'a> {
    cell: &'a primerf::Cell,
    table: Option<primerf::PairSumTable>,
}

impl<'a> PrimeArms<'a> {
    fn decide(
        &self,
        arm: &str,
        r: Pt,
        all: bool,
        out: &mut Vec<primerf::Triple>,
    ) -> (bool, primerf::Ops) {
        match arm {
            "subtract2" => primerf::subtract2(self.cell, r, all, out),
            "pairs3" => primerf::pairs3(self.cell, r, all, out),
            "mitm3" => primerf::mitm3(self.cell, self.table.as_ref().expect("table"), r, all, out),
            "mitm4" => primerf::mitm4(self.cell, self.table.as_ref().expect("table"), r, all, out),
            other => panic!("unknown prime arm {other}"),
        }
    }
}

fn needs_table(arms: &[String]) -> bool {
    arms.iter().any(|a| a == "mitm3" || a == "mitm4")
}

fn median_f64(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn cmd_prime_rows(args: &[String]) {
    let paths: Vec<String> = arg_s(args, "--rows")
        .expect("--rows a.jsonl,b.jsonl")
        .split(',')
        .map(|s| s.to_string())
        .collect();
    let check_max_t: u32 = arg(args, "--check-max-t").unwrap_or(8);
    let brute_max_t: u32 = arg(args, "--brute-max-t").unwrap_or(9);
    let reps: usize = arg(args, "--reps").unwrap_or(5);
    let out = arg_s(args, "--out");
    for path in paths {
        let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
        for line in text.lines() {
            let row: serde_json::Value = match serde_json::from_str(line) {
                Ok(v) => v,
                Err(_) => continue,
            };
            let kind = row["kind"].as_str().unwrap_or("");
            // Only the tower systems themselves: the shape-matched nulls and
            // the generator-count ladders are different systems.
            if row["control"].as_str().unwrap_or("tower") != "tower"
                || row["g"].as_u64().unwrap_or(1) != 1
            {
                continue;
            }
            let (Some(p), Some(t), Some(m), Some(a), Some(b), Some(x_r)) = (
                row["p"].as_u64(),
                row["t"].as_u64(),
                row["m"].as_u64(),
                row["curve"]["a"].as_u64(),
                row["curve"]["b"].as_u64(),
                row["x_r"].as_u64(),
            ) else {
                continue;
            };
            let (t, m) = (t as u32, m as usize);
            let tower = &row["tower"];
            let cell = if let Some(v) = tower.get("v").and_then(|v| v.as_array()) {
                let vv: Vec<u64> = v.iter().map(|x| x.as_u64().unwrap()).collect();
                primerf::Cell::from_v(p, a, b, &vv)
            } else if kind == "kummer" {
                primerf::Cell::kummer(
                    p,
                    a,
                    b,
                    t,
                    tower["g"].as_u64().unwrap(),
                    tower["zeta"].as_u64().unwrap(),
                )
            } else if kind == "dickson" {
                // V = {w + 1/w : w = u·ζ^k}, as pkm_tower_check rebuilds it.
                let f = crypto_lib::cryptanalysis::prime_fast::Fp64::new(p);
                let u = f.to_mont(tower["u"].as_u64().unwrap());
                let z = f.to_mont(tower["zeta"].as_u64().unwrap());
                let mut w = u;
                let mut vv = Vec::with_capacity(1 << t);
                for _ in 0..(1u64 << t) {
                    vv.push(f.from_mont(f.add(w, f.inv(w))));
                    w = f.mul(w, z);
                }
                primerf::Cell::from_v(p, a, b, &vv)
            } else {
                continue;
            };
            let Some(r) = cell.lift(x_r) else {
                println!("{{\"row\":\"{path}\",\"skip\":\"x_r does not lift\"}}");
                continue;
            };
            let arms: Vec<String> = match m {
                2 => vec!["subtract2".into()],
                3 => vec!["pairs3".into(), "mitm3".into()],
                4 => vec!["mitm4".into()],
                _ => continue,
            };
            let t_table = Instant::now();
            let table = if needs_table(&arms) {
                Some(primerf::PairSumTable::build(&cell, true))
            } else {
                None
            };
            let table_secs = t_table.elapsed().as_secs_f64();
            let pa = PrimeArms { cell: &cell, table };
            let brute = if t <= brute_max_t || m == 2 {
                Some(primerf::brute_rational(&cell, r, m))
            } else {
                None
            };
            let alg = if m <= 3 && t <= check_max_t {
                primerf::algebraic_count(&cell, x_r, m)
            } else {
                None
            };
            let f4_inconsistent = row["inconsistent"].as_bool();
            // A verdict exists only when the F4 run finished: not timed out,
            // not stopped for size, no pairs left above its degree bound.
            let f4_finished = row["timed_out"].as_bool() == Some(false)
                && row["oversize"].as_bool() != Some(true)
                && row["pairs_above_bound"].as_u64().unwrap_or(0) == 0;
            let f4_muladds = row["muladds"].as_u64();
            let f4_ms = row["ms"].as_f64();
            let f4_d = row["degree_reached"].as_u64();
            let f4_width = row["max_cols"].as_u64();
            let mut arm_json = Vec::new();
            let mut all_found: Vec<Option<bool>> = Vec::new();
            for arm in &arms {
                let mut roots = Vec::new();
                let (found, ops) = pa.decide(arm, r, true, &mut roots);
                roots.sort_unstable();
                roots.dedup();
                all_found.push(Some(found));
                let brute_agree = brute.as_ref().map(|b| b == &roots);
                let mut times = Vec::new();
                for _ in 0..reps {
                    let mut o = Vec::new();
                    let t0 = Instant::now();
                    black_box(pa.decide(arm, black_box(r), false, &mut o));
                    times.push(t0.elapsed().as_secs_f64() * 1e9);
                }
                arm_json.push(format!(
                    "{{\"arm\":\"{arm}\",\"found\":{found},\"decompositions\":{},\"brute_agree\":{},\"mul\":{},\"sqr\":{},\"inv\":{},\"lookups\":{},\"candidates\":{},\"ns_median\":{:.0},\"ns_min\":{:.0}}}",
                    roots.len(),
                    brute_agree.map(|b| b.to_string()).unwrap_or("null".into()),
                    ops.mul, ops.sqr, ops.inv, ops.lookups, ops.candidates,
                    median_f64(&mut times.clone()),
                    times.iter().cloned().fold(f64::INFINITY, f64::min)
                ));
            }
            let table_bytes = pa.table.as_ref().map(|t| t.bytes()).unwrap_or(0);
            let table_ops = pa.table.as_ref().map(|t| t.build_ops).unwrap_or_default();
            let line = format!(
                "{{\"cmd\":\"prime-rows\",\"row\":\"{path}\",\"kind\":\"{kind}\",\"p\":{p},\"a\":{a},\"b\":{b},\"m\":{m},\"t\":{t},\"N\":{},\"x_r\":{x_r},\"v_size\":{},\"base_points\":{},\"f4\":{{\"inconsistent\":{},\"finished\":{f4_finished},\"muladds\":{},\"ms\":{},\"degree\":{},\"width\":{}}},\"algebraic_ordered_solutions\":{},\"brute_rational\":{},\"table\":{{\"secs\":{table_secs:.4},\"bytes\":{table_bytes},\"mul\":{},\"sqr\":{},\"inv\":{}}},\"arms\":[{}]}}",
                m as u32 * t,
                cell.v.len(),
                cell.base.len(),
                f4_inconsistent.map(|b| b.to_string()).unwrap_or("null".into()),
                f4_muladds.map(|v| v.to_string()).unwrap_or("null".into()),
                f4_ms.map(|v| format!("{v:.3}")).unwrap_or("null".into()),
                f4_d.map(|v| v.to_string()).unwrap_or("null".into()),
                f4_width.map(|v| v.to_string()).unwrap_or("null".into()),
                alg.map(|v| v.to_string()).unwrap_or("null".into()),
                brute.as_ref().map(|b| b.len().to_string()).unwrap_or("null".into()),
                table_ops.mul, table_ops.sqr, table_ops.inv,
                arm_json.join(",")
            );
            append(&out, &line);
        }
    }
}

fn mont_sqr_root_of_unity(
    f: &crypto_lib::cryptanalysis::prime_fast::Fp64,
    p: u64,
    t: u32,
    rng: &mut SplitMix64,
) -> u64 {
    assert!((p - 1).is_multiple_of(1u64 << t), "2^t must divide p − 1");
    loop {
        let h = rng.next() % (p - 2) + 2;
        let z = f.pow(f.to_mont(h), (p - 1) >> t);
        if t == 0 {
            return z;
        }
        let half = f.pow(z, 1u64 << (t - 1));
        if half != f.one {
            return z;
        }
    }
}

fn cmd_prime_ladder(args: &[String]) {
    let p: u64 = arg(args, "--p").unwrap_or(2013265921);
    let m: usize = arg(args, "--m").unwrap_or(3);
    let t_min: u32 = arg(args, "--t-min").unwrap_or(4);
    let t_max: u32 = arg(args, "--t-max").unwrap_or(10);
    let seed: u64 = arg(args, "--seed").unwrap_or(20261012);
    let n_uniform: usize = arg(args, "--uniform").unwrap_or(16);
    let n_planted: usize = arg(args, "--planted").unwrap_or(8);
    let arms: Vec<String> = arg_s(args, "--arms")
        .unwrap_or_else(|| "pairs3,mitm3".into())
        .split(',')
        .map(|s| s.to_string())
        .collect();
    let check_max_t: u32 = arg(args, "--check-max-t").unwrap_or(7);
    let alg_max_t: u32 = arg(args, "--alg-max-t").unwrap_or(7);
    let k_time: usize = arg(args, "--targets").unwrap_or(8);
    let passes_min: usize = arg(args, "--passes").unwrap_or(1);
    let min_secs: f64 = arg(args, "--min-secs").unwrap_or(0.0);
    let out = arg_s(args, "--out");
    let label = arg_s(args, "--label").unwrap_or_default();
    let f = crypto_lib::cryptanalysis::prime_fast::Fp64::new(p);
    for t in t_min..=t_max {
        let mut rng = SplitMix64(seed ^ ((t as u64) << 32) ^ (m as u64));
        // Curve y² = x³ + ax + b, nonsingular.
        let (a, b) = loop {
            let a = rng.next() % p;
            let b = rng.next() % p;
            let am = f.to_mont(a);
            let bm = f.to_mont(b);
            let disc = f.add(
                f.mul(f.to_mont(4), f.mul(f.sqr(am), am)),
                f.mul(f.to_mont(27), f.sqr(bm)),
            );
            if disc != 0 {
                break (a, b);
            }
        };
        let zeta = f.from_mont(mont_sqr_root_of_unity(&f, p, t, &mut rng));
        let g = rng.next() % (p - 1) + 1;
        let cell = primerf::Cell::kummer(p, a, b, t, g, zeta);
        // Targets: uniform abscissae that lift, then planted m-sums.
        let mut targets: Vec<(Pt, bool)> = Vec::new();
        while targets.len() < n_uniform {
            let x = rng.next() % p;
            if x == 0 {
                continue;
            }
            if let Some(pt) = cell.lift(x) {
                targets.push((pt, false));
            }
        }
        let mut planted = 0;
        while planted < n_planted {
            let mut acc: Option<Pt> = None;
            for _ in 0..m {
                let i = (rng.next() as usize) % cell.base.len();
                let pt = if rng.next() & 1 == 1 {
                    cell.curve.neg(cell.base[i])
                } else {
                    cell.base[i]
                };
                acc = cell.curve.add(acc, Some(pt));
            }
            if let Some(r) = acc {
                targets.push((r, true));
                planted += 1;
            }
        }
        let t_table = Instant::now();
        let table = if needs_table(&arms) {
            Some(primerf::PairSumTable::build(&cell, t <= check_max_t))
        } else {
            None
        };
        let table_secs = t_table.elapsed().as_secs_f64();
        let pa = PrimeArms { cell: &cell, table };
        // Exactness: every arm against the brute-force reference (t ≤ check_max_t),
        // the algebraic count (t ≤ alg_max_t, m ≤ 3), and against each other.
        let mut verdicts: Vec<Vec<bool>> = vec![Vec::new(); arms.len()];
        let mut mismatches = 0usize;
        let mut planted_found = 0usize;
        let mut uniform_found = 0usize;
        let mut alg_vs_rational = 0usize;
        for (ti, (r, is_planted)) in targets.iter().enumerate() {
            let brute = if t <= check_max_t {
                Some(primerf::brute_rational(&cell, *r, m))
            } else {
                None
            };
            let alg = if m <= 3 && t <= alg_max_t {
                primerf::algebraic_count(&cell, f.from_mont(r.x), m)
            } else {
                None
            };
            for (ai, arm) in arms.iter().enumerate() {
                let mut roots = Vec::new();
                let (found, _) = pa.decide(arm, *r, true, &mut roots);
                roots.sort_unstable();
                roots.dedup();
                verdicts[ai].push(found);
                if let Some(bf) = &brute {
                    // mitm witnesses are only stored when t ≤ check_max_t; compare sets then.
                    let comparable = arm != "mitm3" && arm != "mitm4"
                        || pa
                            .table
                            .as_ref()
                            .map(|t| t.witnesses.is_some())
                            .unwrap_or(false);
                    if comparable && bf != &roots {
                        mismatches += 1;
                        eprintln!(
                            "MISMATCH t={t} target={ti} arm={arm}: arm {roots:?} brute {bf:?}"
                        );
                    } else if (bf.is_empty()) == found {
                        mismatches += 1;
                        eprintln!(
                            "MISMATCH verdict t={t} target={ti} arm={arm}: found={found} brute={}",
                            bf.len()
                        );
                    }
                }
                if *is_planted && !found {
                    mismatches += 1;
                    eprintln!("MISMATCH planted not found t={t} target={ti} arm={arm}");
                }
                if ai == 0 {
                    if *is_planted && found {
                        planted_found += 1;
                    }
                    if !*is_planted && found {
                        uniform_found += 1;
                    }
                    if let (Some(alg), Some(bf)) = (alg, &brute) {
                        if (alg == 0) != bf.is_empty() {
                            alg_vs_rational += 1;
                        }
                    }
                }
            }
        }
        for ai in 1..arms.len() {
            if verdicts[ai] != verdicts[0] {
                mismatches += 1;
                eprintln!("MISMATCH between arms at t={t}");
            }
        }
        // Timing: the first k uniform targets in full-scan mode (every candidate
        // visited, as a refutation does; a hit does not stop the scan), so that
        // the cost is the refutation cost at every rung, including the rungs
        // where every uniform target decomposes.
        let refuted: Vec<Pt> = targets
            .iter()
            .filter(|(_, pl)| !pl)
            .map(|(r, _)| *r)
            .take(k_time)
            .collect();
        let timed_hits = targets
            .iter()
            .zip(&verdicts[0])
            .filter(|((_, pl), &v)| !pl && v)
            .take(k_time)
            .count();
        for arm in &arms {
            if refuted.is_empty() {
                break;
            }
            let mut o = Vec::new();
            black_box(pa.decide(arm, refuted[0], true, &mut o));
            let t0 = Instant::now();
            let mut passes = 0usize;
            let mut tot = primerf::Ops::default();
            loop {
                for r in &refuted {
                    o.clear();
                    let (_found, ops) = black_box(pa.decide(arm, black_box(*r), true, &mut o));
                    tot.add(ops);
                }
                passes += 1;
                if passes >= passes_min && t0.elapsed().as_secs_f64() >= min_secs {
                    break;
                }
            }
            let secs = t0.elapsed().as_secs_f64();
            let calls = (refuted.len() * passes) as f64;
            let table_bytes = pa.table.as_ref().map(|t| t.bytes()).unwrap_or(0);
            let table_ops = pa.table.as_ref().map(|t| t.build_ops).unwrap_or_default();
            let line = format!(
                "{{\"cmd\":\"prime-ladder\",\"label\":\"{label}\",\"arm\":\"{arm}\",\"p\":{p},\"m\":{m},\"t\":{t},\"N\":{},\"a\":{a},\"b\":{b},\"g\":{g},\"zeta\":{zeta},\"v_size\":{},\"base_points\":{},\"targets_uniform\":{n_uniform},\"targets_planted\":{n_planted},\"uniform_found\":{uniform_found},\"planted_found\":{planted_found},\"checked_brute\":{},\"checked_algebraic\":{},\"algebraic_vs_rational_disagreements\":{alg_vs_rational},\"mismatches\":{mismatches},\"timed_targets\":{},\"timed_mode\":\"full-scan\",\"timed_targets_with_hits\":{timed_hits},\"passes\":{passes},\"secs\":{secs:.6},\"ns_per_target\":{:.1},\"mul_per_target\":{:.1},\"sqr_per_target\":{:.1},\"inv_per_target\":{:.3},\"lookups_per_target\":{:.1},\"candidates_per_target\":{:.1},\"table\":{{\"secs\":{table_secs:.4},\"bytes\":{table_bytes},\"mul\":{},\"sqr\":{},\"inv\":{}}}}}",
                m as u32 * t,
                cell.v.len(),
                cell.base.len(),
                t <= check_max_t,
                m <= 3 && t <= alg_max_t,
                refuted.len(),
                secs * 1e9 / calls,
                tot.mul as f64 / calls, tot.sqr as f64 / calls, tot.inv as f64 / calls, tot.lookups as f64 / calls, tot.candidates as f64 / calls,
                table_ops.mul, table_ops.sqr, table_ops.inv
            );
            append(&out, &line);
        }
        if mismatches != 0 {
            eprintln!("t={t}: {mismatches} mismatches");
            std::process::exit(1);
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("bin-corpus") => cmd_bin_corpus(&args),
        Some("bin-verify") => cmd_bin_verify(&args),
        Some("bin-count") => cmd_bin_count(&args),
        Some("bin-time") => cmd_bin_time(&args),
        Some("prime-rows") => cmd_prime_rows(&args),
        Some("prime-ladder") => cmd_prime_ladder(&args),
        Some("analyze") => analyze::run(&args),
        _ => {
            eprintln!(
                "usage: pdp_rootfind bin-verify|bin-count|bin-time|prime-rows|prime-ladder ..."
            );
            std::process::exit(2);
        }
    }
}
