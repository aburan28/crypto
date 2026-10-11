//! lit2026 W5 probe for the inherited-F4 decomposition oracle.
//!
//! ```bash
//! f4_lane_probe profile   --a 1 --n 17 --dim 6 --targets 8 [--engines f4,f6p]
//! f4_lane_probe lanes     --a 1 --n 17 --dim 6 [--ks 2,4,8,16,32,64] [--reps 5]
//! f4_lane_probe lockstep  --a 1 --n 17 --dim 6 --targets 8 [--engines f4,f6p]
//! f4_lane_probe costmodel --a 1 --n 17 --dim 6 [--degrees 2,3,4]
//! ```
//!
//! - `profile`: where an arm's time goes — root builds against
//!   specialisations — with the stage unit, Macaulay rows generated and
//!   the stage units charged to root eliminations.
//! - `lanes`: the root phase alone, for `k` targets `[1]G … [k]G`: the `k`
//!   roots built one by one ([`ReducedBasis::from_system_with`], the
//!   production echelon root, or the RREF root under
//!   `KIC_F4_INHERIT_ROOT=rref`) against built together
//!   ([`ReducedBasis::from_systems_lanes`]), interleaved per repetition;
//!   each lane is checked to span its sequential root's row space (equal
//!   rows when the sequential root is RREF).
//! - `lockstep`: how much of the targets' splitting trees could run in
//!   lockstep — the nodes (assignments) that two or more targets reach —
//!   from the solver's node trace.  An upper bound for lane-batching below
//!   the root: lanes at one node must also share their pivot pattern.
//! - `costmodel`: root Macaulay matrix rows, columns and rank per degree,
//!   on the descended system and on a dense random system of the same
//!   shape, against the generic Boolean counts (Macaulay rows and columns,
//!   rank from the semi-regular Hilbert series).

use crypto_lib::cryptanalysis::inherited_f4::{
    f4_root_elimination_thread, f4_rows_generated_thread, macaulay_shape, ReducedBasis,
};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    f4_profile, f4_profile_reset, f4_word_ops_thread, inherited_basis_profile,
    inherited_basis_profile_reset, set_f4_inherited_basis_profile, set_node_trace, take_node_trace,
    FieldStructure, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, decomposition_solver_equations, groebner_decompose,
    groebner_decompose_f6_ic_pair_index, KoblitzCurve,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use num_bigint::BigUint;
use std::collections::HashMap;
use std::time::Instant;

fn flag(name: &str) -> Option<String> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1).cloned())
}

fn list(name: &str, default: &str) -> Vec<u32> {
    flag(name)
        .unwrap_or_else(|| default.to_string())
        .split(',')
        .map(|v| v.parse().expect("a number"))
        .collect()
}

/// User plus system CPU time of this process, in nanoseconds.
fn cpu_ns() -> u128 {
    // SAFETY: `getrusage` writes into the zeroed struct it is handed.
    let usage = unsafe {
        let mut usage: libc::rusage = std::mem::zeroed();
        libc::getrusage(libc::RUSAGE_SELF, &mut usage);
        usage
    };
    let tv = |t: libc::timeval| t.tv_sec as u128 * 1_000_000_000 + t.tv_usec as u128 * 1_000;
    tv(usage.ru_utime) + tv(usage.ru_stime)
}

fn binom(n: u64, k: u64) -> u128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut r: u128 = 1;
    for i in 0..k {
        r = r * (n - i) as u128 / (i + 1) as u128;
    }
    r
}

/// Coefficients `0..=top` of `(1 + t)^n / Π_j (1 + t^{d_j})`, truncated
/// before its first non-positive coefficient (the semi-regular Boolean
/// Hilbert series).
fn boolean_hilbert(n: u64, degrees: &[u32], top: usize) -> Vec<i128> {
    let mut h: Vec<i128> = (0..=top).map(|i| binom(n, i as u64) as i128).collect();
    for &d in degrees {
        // divide by (1 + t^d): h'_i = h_i − h'_{i−d}
        let d = d as usize;
        if d == 0 {
            continue;
        }
        for i in d..=top {
            h[i] -= h[i - d];
        }
    }
    if let Some(first) = h.iter().position(|&c| c <= 0) {
        for c in h.iter_mut().skip(first) {
            *c = 0;
        }
    }
    h
}

fn degree_of(p: &F2BoolPoly) -> u32 {
    p.terms
        .iter()
        .map(|t| t.mask.count_ones())
        .max()
        .unwrap_or(0)
}

fn occurring(system: &[F2BoolPoly]) -> u64 {
    system
        .iter()
        .flat_map(|p| p.terms.iter())
        .fold(0u64, |acc, t| acc | t.mask)
}

struct Cell {
    kc: KoblitzCurve,
    fb: crypto_lib::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    st: FieldStructure,
    m: usize,
    engine: SolverEngine,
}

impl Cell {
    fn target(&self, k: u64) -> crypto_lib::binary_ecc::BinaryPoint {
        self.kc.mul(self.kc.generator(), &BigUint::from(k))
    }
    fn system(&self, k: u64) -> (Vec<F2BoolPoly>, usize) {
        decomposition_solver_equations(
            &self.kc,
            &self.fb,
            &self.st,
            &self.target(k),
            self.m,
            self.engine,
        )
        .expect("a supported target")
    }
}

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_default();
    let a: u8 = flag("--a").and_then(|v| v.parse().ok()).unwrap_or(0);
    let n: u32 = flag("--n").and_then(|v| v.parse().ok()).unwrap_or(17);
    let dim: u32 = flag("--dim").and_then(|v| v.parse().ok()).unwrap_or(6);
    let m: usize = flag("--m").and_then(|v| v.parse().ok()).unwrap_or(3);
    let targets: u64 = flag("--targets").and_then(|v| v.parse().ok()).unwrap_or(8);
    let budget: usize = flag("--budget")
        .and_then(|v| v.parse().ok())
        .unwrap_or(20_000);
    let reps: usize = flag("--reps").and_then(|v| v.parse().ok()).unwrap_or(5);
    let engines = flag("--engines").unwrap_or_else(|| "f4,f6p".to_string());
    let kc = KoblitzCurve::new(a, n).expect("curve has a usable subgroup");
    let fb = build_standard_subspace_factor_base(&kc, dim).expect("factor base");
    let index_of = fb.index_map();
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let engine = SolverEngine::InheritedF4 { max_degree: 3 };
    let cell = Cell {
        kc,
        fb,
        st,
        m,
        engine,
    };
    let decide = |arm: &str, k: u64| {
        let target = cell.target(k);
        match arm {
            "f4" => groebner_decompose(
                &cell.kc, &cell.fb, &index_of, &cell.st, &target, m, engine, budget,
            ),
            "f6p" => groebner_decompose_f6_ic_pair_index(
                &cell.kc, &cell.fb, &index_of, &cell.st, &target, m, engine, budget,
            ),
            other => panic!("unknown engine {other}; use f4 or f6p"),
        }
    };
    let cell_id = format!("K{a}/2^{n} dim {dim}");
    match mode.as_str() {
        "profile" => {
            set_f4_inherited_basis_profile(true);
            for arm in engines.split(',') {
                inherited_basis_profile_reset();
                f4_profile_reset();
                let x0 = f4_word_ops_thread();
                let r0 = f4_rows_generated_thread();
                let e0 = f4_root_elimination_thread();
                let c0 = cpu_ns();
                let t0 = Instant::now();
                for k in 1..=targets {
                    let _ = decide(arm, k);
                }
                let ms = t0.elapsed().as_secs_f64() * 1e3;
                let cpu = (cpu_ns() - c0) as f64 / 1e6;
                let p = inherited_basis_profile();
                let f = f4_profile();
                let r = f4_rows_generated_thread();
                println!(
                    "{{\"cell\":\"{cell_id}\",\"arm\":\"{arm}\",\"targets\":{targets},\"ms\":{ms:.3},\"cpu_ms\":{cpu:.3},\"word_xors\":{},\"root_elim_xors\":{},\"rows_generated\":{},\"row_words_generated\":{},\"root_calls\":{},\"root_ms\":{:.3},\"specialise_calls\":{},\"specialise_ms\":{:.3},\"basis_specialise_ms\":{:.3},\"build_ms\":{:.3},\"reduce_ms\":{:.3},\"calls\":{},\"rank_sum\":{},\"cols_sum\":{}}}",
                    f4_word_ops_thread() - x0,
                    f4_root_elimination_thread() - e0,
                    r.0 - r0.0,
                    r.1 - r0.1,
                    p.root_calls,
                    p.root_ns as f64 / 1e6,
                    p.specialise_calls,
                    p.specialise_total_ns as f64 / 1e6,
                    p.basis_specialise_ns as f64 / 1e6,
                    f.build_ns as f64 / 1e6,
                    f.reduce_ns as f64 / 1e6,
                    f.calls,
                    f.rows,
                    f.cols,
                );
            }
        }
        "lanes" => {
            let ks = list("--ks", "2,4,8,16,32,64");
            let kmax = *ks.iter().max().unwrap() as u64;
            let systems: Vec<Vec<F2BoolPoly>> = (1..=kmax).map(|k| cell.system(k).0).collect();
            let n_vars = cell.system(1).1;
            let degree = 3;
            for &k in &ks {
                let k = k as usize;
                let slices: Vec<&[F2BoolPoly]> =
                    systems[..k].iter().map(|s| s.as_slice()).collect();
                let mut seq_ms = Vec::new();
                let mut lane_ms = Vec::new();
                let mut seq_cpu = Vec::new();
                let mut lane_cpu = Vec::new();
                let mut seq_ops = 0u64;
                let mut lane_stats = None;
                for _ in 0..reps {
                    // sequential
                    let e0 = f4_root_elimination_thread();
                    let c0 = cpu_ns();
                    let t0 = Instant::now();
                    let seq: Vec<ReducedBasis> = slices
                        .iter()
                        .map(|s| {
                            ReducedBasis::from_system_with(s, n_vars, degree, true)
                                .expect("fits")
                                .0
                        })
                        .collect();
                    seq_ms.push(t0.elapsed().as_secs_f64() * 1e3);
                    seq_cpu.push((cpu_ns() - c0) as f64 / 1e6);
                    seq_ops = f4_root_elimination_thread() - e0;
                    // lanes
                    let c0 = cpu_ns();
                    let t0 = Instant::now();
                    let (lanes, stats) =
                        ReducedBasis::from_systems_lanes(&slices, n_vars, degree, true);
                    lane_ms.push(t0.elapsed().as_secs_f64() * 1e3);
                    lane_cpu.push((cpu_ns() - c0) as f64 / 1e6);
                    // identity: each lane spans its sequential root's row space
                    for (i, (s, l)) in seq.iter().zip(&lanes).enumerate() {
                        let l = &l.as_ref().expect("lane fits").0;
                        assert!(
                            l.same_row_space(s),
                            "{cell_id}: lane {i} of {k} differs from its sequential root"
                        );
                    }
                    lane_stats = Some(stats);
                }
                let s = lane_stats.unwrap();
                let med = |v: &mut Vec<f64>| {
                    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    (v[v.len() / 2], v[0], v[v.len() - 1])
                };
                let (sm, smin, smax) = med(&mut seq_ms);
                let (lm, lmin, lmax) = med(&mut lane_ms);
                let (sc, _, _) = med(&mut seq_cpu);
                let (lc, _, _) = med(&mut lane_cpu);
                println!(
                    "{{\"cell\":\"{cell_id}\",\"k\":{k},\"reps\":{reps},\"seq_ms\":[{sm:.3},{smin:.3},{smax:.3}],\"lane_ms\":[{lm:.3},{lmin:.3},{lmax:.3}],\"seq_cpu_ms\":{sc:.3},\"lane_cpu_ms\":{lc:.3},\"seq_root_xors\":{seq_ops},\"lane_word_equivalent\":{},\"lane_word_ops\":{},\"lane_bits\":{},\"rows\":{},\"columns\":{},\"pivot_columns\":{},\"pivot_passes\":{},\"shared_columns\":{}}}",
                    s.word_equivalent,
                    s.lane_word_ops,
                    s.lane_bits,
                    s.rows,
                    s.columns,
                    s.pivot_columns,
                    s.pivot_passes,
                    s.shared_columns,
                );
            }
        }
        "lockstep" => {
            for arm in engines.split(',') {
                // node key → (targets reaching it, work)
                let mut nodes: HashMap<(u64, u64), (u32, u64)> = HashMap::new();
                let mut visits = 0u64;
                let mut work_total = 0u64;
                for k in 1..=targets {
                    set_node_trace(true);
                    let _ = decide(arm, k);
                    let end = f4_word_ops_thread();
                    let events = take_node_trace();
                    set_node_trace(false);
                    let mut seen: HashMap<(u64, u64), u64> = HashMap::new();
                    for (i, e) in events.iter().enumerate() {
                        let next = events.get(i + 1).map_or(end, |n| n.word_ops);
                        *seen.entry((e.assigned, e.values)).or_default() += next - e.word_ops;
                        visits += 1;
                    }
                    for (key, work) in seen {
                        let entry = nodes.entry(key).or_default();
                        entry.0 += 1;
                        entry.1 += work;
                        work_total += work;
                    }
                }
                let distinct = nodes.len() as u64;
                let shared_nodes = nodes.values().filter(|v| v.0 >= 2).count() as u64;
                let shared_work: u64 = nodes.values().filter(|v| v.0 >= 2).map(|v| v.1).sum();
                let lane_visits: u64 = nodes.values().map(|v| v.0 as u64).sum();
                // Work a perfect lockstep would leave: one lane-batched pass
                // per distinct node, costing what its most expensive lane cost
                // is not known; bound it by work / lanes at the node.
                let mut hist = [0u64; 9];
                for v in nodes.values() {
                    hist[(v.0 as usize).min(8)] += 1;
                }
                println!(
                    "{{\"cell\":\"{cell_id}\",\"arm\":\"{arm}\",\"targets\":{targets},\"node_visits\":{visits},\"distinct_nodes\":{distinct},\"distinct_node_target_pairs\":{lane_visits},\"mean_lanes_per_node\":{:.4},\"nodes_shared_by_2plus\":{shared_nodes},\"work_total\":{work_total},\"work_at_shared_nodes\":{shared_work},\"work_share_shared\":{:.4},\"lanes_histogram\":{:?}}}",
                    lane_visits as f64 / distinct.max(1) as f64,
                    shared_work as f64 / work_total.max(1) as f64,
                    hist,
                );
            }
        }
        "costmodel" => {
            let degrees = list("--degrees", "2,3,4");
            let (system, n_vars) = cell.system(1);
            let system: Vec<F2BoolPoly> = system.into_iter().filter(|p| !p.is_zero()).collect();
            let occ = occurring(&system);
            let n_occ = occ.count_ones() as u64;
            let gen_degrees: Vec<u32> = system.iter().map(degree_of).collect();
            let mut by_degree: HashMap<u32, usize> = HashMap::new();
            for &d in &gen_degrees {
                *by_degree.entry(d).or_default() += 1;
            }
            let mean_support: f64 = system
                .iter()
                .map(|p| occurring(std::slice::from_ref(p)).count_ones() as f64)
                .sum::<f64>()
                / system.len() as f64;
            let terms: usize = system.iter().map(|p| p.terms.len()).sum();
            // A dense random system of the same shape: every monomial of
            // degree ≤ d_j over the occurring variables with probability 1/2.
            let vars: Vec<u32> = (0..64).filter(|v| occ & (1u64 << v) != 0).collect();
            let mut seed = 0x9E37_79B9_7F4A_7C15u64 ^ u64::from(n);
            let mut next = || {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                seed
            };
            let dense: Vec<F2BoolPoly> = gen_degrees
                .iter()
                .map(|&d| {
                    let mut monos = Vec::new();
                    let mut stack: Vec<(u64, usize, u32)> = vec![(0, 0, 0)];
                    while let Some((mask, from, deg)) = stack.pop() {
                        if next() & 1 == 1 || (mask == 0) {
                            monos.push(F2BoolMono::from_mask(mask));
                        }
                        if deg < d {
                            for (j, &v) in vars.iter().enumerate().skip(from) {
                                stack.push((mask | (1u64 << v), j + 1, deg + 1));
                            }
                        }
                    }
                    F2BoolPoly::from_monos(monos, n_vars)
                })
                .collect();
            let top = *degrees.iter().max().unwrap() as usize;
            let hilbert = boolean_hilbert(n_occ, &gen_degrees, top);
            for &d in &degrees {
                let cols_generic: u128 = (0..=d as u64).map(|i| binom(n_occ, i)).sum();
                let rows_generic: u128 = gen_degrees
                    .iter()
                    .filter(|&&g| g <= d)
                    .map(|&g| (0..=(d - g) as u64).map(|i| binom(n_occ, i)).sum::<u128>())
                    .sum();
                let staircase: i128 = hilbert[..=d as usize].iter().sum();
                let rank_generic = (cols_generic as i128 - staircase).min(rows_generic as i128);
                let shape =
                    |s: &[F2BoolPoly], local: bool| match macaulay_shape(s, n_vars, d, local) {
                        Some((r, c, k)) => format!("[{r},{c},{k}]"),
                        None => "null".to_string(),
                    };
                println!(
                    "{{\"cell\":\"{cell_id}\",\"degree\":{d},\"n_vars\":{n_vars},\"n_occurring\":{n_occ},\"generators\":{},\"generator_degrees\":{:?},\"mean_support\":{mean_support:.2},\"terms\":{terms},\"generic\":[{rows_generic},{cols_generic},{rank_generic}],\"descended_support_local\":{},\"descended_occurring\":{},\"dense_support_local\":{},\"dense_occurring\":{}}}",
                    system.len(),
                    {
                        let mut v: Vec<(u32, usize)> = by_degree.iter().map(|(a, b)| (*a, *b)).collect();
                        v.sort_unstable();
                        v
                    },
                    shape(&system, true),
                    shape(&system, false),
                    shape(&dense, true),
                    shape(&dense, false),
                );
            }
        }
        other => panic!("unknown mode {other:?}; use profile, lanes, lockstep or costmodel"),
    }
}
