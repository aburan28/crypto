//! **Volcano-level panel at arity m on a walked Koblitz isogeny class.**
//!
//! A variant of `koblitz_presentation_sweep` that (a) takes the summand count
//! `m` as an argument so the Gröbner step is non-trivial (m = 3 builds the
//! Weil-descended `S_4` system), and (b) measures an explicit list of class
//! members, so a caller can stratify by the per-prime volcano position
//! `(v_ℓ(f_E))_ℓ` read off the ℓ-labelled edges of
//! `experiments/koblitz_isogeny_class_walk.json`.  Everything else (full-group
//! probes, fixed probe budget per cell, BSGS-checked transported logarithm,
//! append-per-cell checkpointing) is the sweep's.
//!
//! Run: `cargo run --release --example koblitz_level_panel <n> <a2> <l> <m> <probes> <seeds,..|mono> <a6,a6,..|all>`
//! `seeds` = `mono` uses the monomial subspace `V = ⟨1, z, …, z^{l−1}⟩`
//! (the only one affordable at m = 3); otherwise each seed draws a random V.
//! Cells are appended to `$KOBLITZ_PRES_CHECKPOINT/levelpanel_<n>_<a2>_l<l>_m<m>.jsonl`
//! (default `experiments/`) and skipped on rerun.  The per-ℓ level vector is
//! joined on `a6` by the analysis script; this binary records `a6` and the
//! breadth-first depth only.

use std::collections::{BTreeMap, HashSet, VecDeque};
use std::io::Write;

use crypto_lib::cryptanalysis::koblitz_isogeny_cost::*;
use rayon::prelude::*;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.len() < 7 {
        eprintln!("usage: koblitz_level_panel <n> <a2> <l> <m> <probes> <seeds,..|mono> <a6,..|all>");
        std::process::exit(2);
    }
    let n: u32 = a[0].parse().unwrap();
    let a2: u8 = a[1].parse().unwrap();
    let l: u32 = a[2].parse().unwrap();
    let m: usize = a[3].parse().unwrap();
    let probes: usize = a[4].parse().unwrap();
    let seeds: Vec<Option<u64>> = if a[5] == "mono" {
        vec![None]
    } else {
        a[5].split(',').map(|s| Some(s.parse().unwrap())).collect()
    };

    let text = std::fs::read_to_string("experiments/koblitz_isogeny_class_walk.json").unwrap();
    let v: serde_json::Value = serde_json::from_str(&text).unwrap();
    let case = v["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["n"].as_u64() == Some(n as u64) && c["a2"].as_u64() == Some(a2 as u64))
        .expect("class walked");
    let mut adj: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
    for e in case["edge_list"].as_array().unwrap() {
        let (f, t) = (e["from"].as_u64().unwrap(), e["to"].as_u64().unwrap());
        adj.entry(f).or_default().push(t);
    }
    let mut depth: BTreeMap<u64, u32> = BTreeMap::new();
    let mut q = VecDeque::from([(1u64, 0u32)]);
    while let Some((x, d)) = q.pop_front() {
        if depth.contains_key(&x) {
            continue;
        }
        depth.insert(x, d);
        for &y in adj.get(&x).into_iter().flatten() {
            q.push_back((y, d + 1));
        }
    }
    let members: Vec<u64> = if a[6] == "all" {
        depth.keys().copied().collect()
    } else {
        a[6].split(',').map(|s| s.parse().unwrap()).collect()
    };
    for &x in &members {
        assert!(depth.contains_key(&x), "a6 = {x} is not a walked member of the class");
    }

    let dir = std::env::var("KOBLITZ_PRES_CHECKPOINT").unwrap_or_else(|_| "experiments".into());
    let path = std::path::PathBuf::from(dir).join(format!("levelpanel_{n}_{a2}_l{l}_m{m}.jsonl"));
    let done: HashSet<(u64, i64)> = std::fs::read_to_string(&path)
        .unwrap_or_default()
        .lines()
        .filter_map(|line| serde_json::from_str::<serde_json::Value>(line).ok())
        .map(|r| (r["a6"].as_u64().unwrap(), r["basis_seed"].as_i64().unwrap_or(-1)))
        .collect();
    let mut cells: Vec<(u64, Option<u64>)> = members
        .iter()
        .flat_map(|&x| seeds.iter().map(move |&s| (x, s)))
        .filter(|&(x, s)| !done.contains(&(x, s.map(|v| v as i64).unwrap_or(-1))))
        .collect();
    {
        use rand::seq::SliceRandom;
        use rand::SeedableRng;
        cells.shuffle(&mut rand::rngs::StdRng::seed_from_u64(DEFAULT_SEED ^ n as u64 ^ (m as u64) << 8));
    }
    eprintln!("n={n} a2={a2} l={l} m={m}: {} members x {} V, {} cells left", members.len(), seeds.len(), cells.len());

    let ffd_targets: u32 = std::env::var("KOBLITZ_PRES_FFD_TARGETS").ok().and_then(|v| v.parse().ok()).unwrap_or(0);
    let threads: usize = std::env::var("KOBLITZ_PANEL_THREADS").ok().and_then(|v| v.parse().ok()).unwrap_or(3);
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();
    let irr = field_for(n).unwrap();
    let order = koblitz_family_order(n, a2) as u64;
    let file = std::sync::Mutex::new(std::fs::OpenOptions::new().create(true).append(true).open(&path).unwrap());
    let finished = std::sync::atomic::AtomicUsize::new(0);
    let total = cells.len();
    let t0 = std::time::Instant::now();
    cells.par_iter().with_max_len(1).for_each(|&(a6, seed)| {
        let opts = IcCostOptions {
            l,
            m,
            yield_probes: probes,
            max_trials: 0,
            ffd_targets,
            ffd_d_max: 6,
            known_order: Some(order),
            full_group_probes: true,
            basis_seed: seed,
            geometric_basis: false,
            // m = 3 calls run hundreds of reductions; give each decomposition
            // room so a cell is a measurement, not a budget artefact.
            node_budget: 2_000_000,
            ..Default::default()
        };
        let t1 = std::time::Instant::now();
        let line = match measure_member(n, &irr, a2, a6, &opts) {
            Ok(r) => serde_json::json!({
                "n": n, "a2": a2, "l": l, "m": m, "a6": a6, "basis_seed": seed, "depth": depth[&a6],
                "probes": r.trials, "relations": r.relations, "factor_base_points": r.factor_base_points,
                "unknowns": r.unknowns, "groebner_ns": r.groebner_ns, "groebner_calls": r.groebner_calls,
                "reductions": r.reductions, "infeasible_branches": r.infeasible_branches,
                "refutations": r.refutations, "unknown_probes": r.unknown_probes,
                "first_fall_hist": r.first_fall_hist, "d_star_hist": r.d_star_hist,
                "inconsistent": r.inconsistent_relations, "independent_relations": r.independent_relations,
                "bsgs_ok": r.bsgs_log == Some(r.planted), "r": r.r, "cofactor": r.cofactor,
                "cell_wall_s": t1.elapsed().as_secs_f64(),
            }),
            Err(e) => serde_json::json!({"n": n, "a2": a2, "l": l, "m": m, "a6": a6, "basis_seed": seed,
                                         "depth": depth[&a6], "skipped": format!("{e:?}")}),
        };
        writeln!(file.lock().unwrap(), "{line}").unwrap();
        let k = finished.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
        eprintln!("   {k}/{total} cells, {:.0}s", t0.elapsed().as_secs_f64());
    });
}
