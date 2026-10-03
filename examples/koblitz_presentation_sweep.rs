//! **Curve effect vs presentation effect in index-calculus cost.**
//!
//! Every (class member, factor-base subspace V) cell of a walked isogeny
//! class is measured with full-group probes at a fixed probe budget, so a
//! random-effects analysis can split cost variation into a curve part and a
//! presentation part.  Design and pre-registered criteria:
//! `research/notes/koblitz-isogeny/presentation-vs-curve-design-20261003.md`.
//!
//! Run: `cargo run --release --example koblitz_presentation_sweep <n> <a2> <l> <probes> <sample|all> <seeds,..>`
//! Members come from `experiments/koblitz_isogeny_class_walk.json`.
//! Cells are appended to `$KOBLITZ_PRES_CHECKPOINT/pres_<n>_<a2>_l<l>.jsonl`
//! (default `experiments/`) and skipped on rerun.

use std::collections::{BTreeMap, HashSet, VecDeque};
use std::io::Write;

use crypto_lib::cryptanalysis::koblitz_isogeny_cost::*;
use rayon::prelude::*;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let n: u32 = a[0].parse().unwrap();
    let a2: u8 = a[1].parse().unwrap();
    let l: u32 = a[2].parse().unwrap();
    let probes: usize = a[3].parse().unwrap();
    let sample: Option<usize> = a[4].parse().ok();
    let seeds: Vec<u64> = a[5].split(',').map(|s| s.parse().unwrap()).collect();

    // members and their depth below K_a, from the walk's edges
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
    let mix = |x: u64| {
        let mut h = x ^ DEFAULT_SEED;
        h ^= h >> 33;
        h = h.wrapping_mul(0xFF51_AFD7_ED55_8CCD);
        h ^ (h >> 33)
    };
    let mut members: Vec<u64> = depth.keys().copied().filter(|&x| x != 1).collect();
    members.sort_by_key(|&x| mix(x));
    let mut members: Vec<u64> = std::iter::once(1).chain(members).collect();
    if let Some(k) = sample {
        members.truncate(k);
    }

    let dir = std::env::var("KOBLITZ_PRES_CHECKPOINT").unwrap_or_else(|_| "experiments".into());
    let path = std::path::PathBuf::from(dir).join(format!("pres_{n}_{a2}_l{l}.jsonl"));
    let done: HashSet<(u64, u64)> = std::fs::read_to_string(&path)
        .unwrap_or_default()
        .lines()
        .filter_map(|line| serde_json::from_str::<serde_json::Value>(line).ok())
        .map(|r| (r["a6"].as_u64().unwrap(), r["basis_seed"].as_u64().unwrap()))
        .collect();
    let cells: Vec<(u64, u64)> = members
        .iter()
        .flat_map(|&m| seeds.iter().map(move |&s| (m, s)))
        .filter(|c| !done.contains(c))
        .collect();
    eprintln!(
        "n={n} a2={a2} l={l}: {} members x {} V, {} cells left",
        members.len(),
        seeds.len(),
        cells.len()
    );

    let irr = field_for(n).unwrap();
    let order = koblitz_family_order(n, a2) as u64;
    let file = std::sync::Mutex::new(
        std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .unwrap(),
    );
    let finished = std::sync::atomic::AtomicUsize::new(0);
    let total = cells.len();
    let t0 = std::time::Instant::now();
    cells.par_iter().for_each(|&(a6, seed)| {
        let opts = IcCostOptions {
            l,
            m: 2,
            yield_probes: probes,
            max_trials: 0,
            ffd_targets: 4,
            ffd_d_max: 6,
            known_order: Some(order),
            full_group_probes: true,
            basis_seed: Some(seed),
            ..Default::default()
        };
        let line = match measure_member(n, &irr, a2, a6, &opts) {
            Ok(r) => serde_json::json!({
                "n": n, "a2": a2, "l": l, "a6": a6, "basis_seed": seed, "depth": depth[&a6],
                "probes": r.trials, "relations": r.relations, "factor_base_points": r.factor_base_points,
                "unknowns": r.unknowns, "groebner_ns": r.groebner_ns, "groebner_calls": r.groebner_calls,
                "reductions": r.reductions, "first_fall_hist": r.first_fall_hist,
                "d_star_hist": r.d_star_hist, "inconsistent": r.inconsistent_relations,
                "bsgs_ok": r.bsgs_log == Some(r.planted), "r": r.r, "cofactor": r.cofactor,
            }),
            Err(e) => serde_json::json!({"n": n, "a2": a2, "l": l, "a6": a6, "basis_seed": seed,
                                         "depth": depth[&a6], "skipped": format!("{e:?}")}),
        };
        writeln!(file.lock().unwrap(), "{line}").unwrap();
        let k = finished.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
        if k % 32 == 0 || k == total {
            eprintln!("   {k}/{total} cells, {:.0}s", t0.elapsed().as_secs_f64());
        }
    });
}
