//! Full pricing: every native counter of a solve charged in the unit.
//!
//! The ordinary methods price native work with the ratios pinned in
//! `docs/ic/calibration.json`.  That table has no field multiplication,
//! squaring or inversion, and some curves (`K_0 / GF(2^19)`) have no row
//! at all, so a Frobenius-orbit factor base's construction and the
//! canonicalisations of the Frobenius walks are counted under
//! `*_uncharged` names and the total is a lower bound.
//!
//! The `*_full` methods run the same solver unchanged and then price
//! every counter from one table, `docs/ecbench/full_prices.json`, which
//! `ecbench price-units` measures on one host for every unit and curve
//! at once.  A counter the rules below do not know stays unpriced, so the
//! record still says it is a lower bound rather than silently dropping
//! work.  The pinned table and the methods that read it are untouched,
//! so every committed session replays as before.

use std::collections::BTreeMap;
use std::time::Instant;

use rand::{rngs::StdRng, Rng, SeedableRng};
use serde_json::{json, Value};

use crate::cryptanalysis::ecbench::canonical::sha256_hex;
use crate::cryptanalysis::ecbench::methods::{PhaseRecord, SolveReport};
use crate::cryptanalysis::ic_boundary::{
    calibrate_binary_instance, fast_map, BinaryInstance, FastMap,
};

/// The measured table, compiled in so that a record's price is a
/// function of its counts and this file alone.
pub const FULL_PRICES: &str = include_str!("../../../docs/ecbench/full_prices.json");

/// The method a `*_full` id runs before repricing.
pub fn base_of(id: &str) -> Option<&'static str> {
    match id {
        "ic.pipeline_full" => Some("ic.pipeline"),
        "rho.signed_frobenius_full" => Some("rho.signed_frobenius"),
        "rho.signed_frobenius_strong_full" => Some("rho.signed_frobenius_strong"),
        _ => None,
    }
}

/// The units of the table, as ratios to one affine addition.
pub const UNITS: &[&str] = &[
    "double",
    "as_solve",
    "lookup",
    "row_op",
    "frobenius",
    "canon",
    "field_mul",
    "field_sqr",
    "field_inv",
    "word_xor",
];

/// How one counter is charged: work already inside the phase's
/// additions and doublings (or a count of outcomes, not of work), or a
/// sum of units.
enum Rule {
    InGroupOps,
    Units(&'static [(&'static str, f64)]),
}

/// Every counter the IC pipeline (compact-orbit-scan and
/// koblitz-standard-subspace bases; mitm-frobenius-counted and
/// pdp3-koblitz rootfind oracles; walk, stride and stride-renew
/// targets; incremental Gauss) and the two Frobenius rho walks emit.
fn rule(counter: &str) -> Option<Rule> {
    use Rule::*;
    Some(match counter {
        // Counts of outcomes or of work already charged as additions.
        "cofactor_projected_points"
        | "multiplier_points"
        | "pair_table_representatives"
        | "pair_table_kind"
        | "pair_table_points"
        | "pair_table_build_adds"
        | "derived_counts"
        | "abscissae_scanned"
        | "direct_relations_skipped"
        | "frobfold_mismatches"
        | "lift_failures"
        | "relations"
        | "repeated_targets_skipped"
        | "s4_hits"
        | "trials"
        | "unliftable_systems"
        | "walk_jumps"
        | "walk_pool_offsets"
        | "walk_pool_ops"
        | "walk_pool_restarts"
        | "walk_restarts"
        | "walk_steps"
        | "stride_renewals"
        | "stride_restarts"
        | "stride_steps"
        | "columns"
        | "pinned_by_repeated_row"
        | "rank"
        | "repeated_column_rows"
        | "rows"
        | "single_column_rows"
        | "two_column_rows"
        // rho.signed_frobenius
        | "cycle_escape_doublings"
        | "cycles_length_2"
        | "cycles_length_4"
        | "cycles_other_length"
        | "distinguished_point_bits"
        | "jump_coefficient_bits"
        | "jumps"
        | "look_ahead_additions"
        | "setup_additions"
        | "setup_doublings"
        | "start_additions"
        | "steps"
        | "useless_collisions"
        | "verification_operations"
        | "walk_operations"
        | "walks"
        | "walks_abandoned_in_a_cycle"
        | "walks_capped"
        // rho.signed_frobenius_strong (scalar multiplications are in its
        // phase charge at 1.5·log2 r additions each)
        | "capped_walks"
        | "failed_collisions"
        | "fruitless_walks"
        | "scalar_multiplications"
        | "table_entries"
        | "wasted_merges" => InGroupOps,
        "as_solves" | "legacy_as_solves_uncharged" => Units(&[("as_solve", 1.0)]),
        // Hash probes and inserts; the scan's per-abscissa loop and the
        // per-lift bookkeeping are charged one probe each.
        "lookups"
        | "target_guard_probes"
        | "hash_probes_uncharged"
        | "factor_point_inserts_uncharged"
        | "pair_table_entries"
        | "abscissae_scanned_uncharged"
        | "point_lifts_uncharged"
        | "table_inserts_uncharged"
        | "table_queries_uncharged"
        | "distinguished_points" => Units(&[("lookup", 1.0)]),
        "row_ops" => Units(&[("row_op", 1.0)]),
        "frobenius_maps" => Units(&[("frobenius", 1.0)]),
        "canonicalisations" | "canonicalisations_uncharged" => Units(&[("canon", 1.0)]),
        "inversions" => Units(&[("field_inv", 1.0)]),
        "lift_field_squares_uncharged" | "min_key_field_squares_uncharged" => {
            Units(&[("field_sqr", 1.0)])
        }
        "lift_field_multiplies_uncharged" => Units(&[("field_mul", 1.0)]),
        // y² + xy = x³ + ax² + b: two squarings, two multiplications.
        "curve_equation_checks_uncharged" => Units(&[("field_mul", 2.0), ("field_sqr", 2.0)]),
        // −(x, y) = (x, x + y): one word XOR.
        "point_negations_uncharged" | "partition_hashes_uncharged" => {
            Units(&[("word_xor", 1.0)])
        }
        // Writing the base: one word store per eight bytes.
        "factor_base_vec_bytes_floor_uncharged" => Units(&[("word_xor", 0.125)]),
        _ => return None,
    })
}

/// The table row for one curve, as unit → ratio to an addition.
pub fn row(slug: &str) -> Result<BTreeMap<String, f64>, String> {
    let doc: Value = serde_json::from_str(FULL_PRICES).map_err(|e| e.to_string())?;
    let r = doc
        .get("instances")
        .and_then(|i| i.get(slug))
        .and_then(Value::as_object)
        .ok_or_else(|| {
            format!("docs/ecbench/full_prices.json has no row for {slug}; measure it with `ecbench price-units`")
        })?;
    let mut out = BTreeMap::new();
    for u in UNITS {
        let v = r
            .get(*u)
            .and_then(Value::as_f64)
            .ok_or_else(|| format!("full_prices.json row {slug} lacks `{u}`"))?;
        out.insert((*u).to_string(), v);
    }
    Ok(out)
}

/// Price one counter; `None` when no rule covers it.
fn price(counter: &str, count: u64, units: &BTreeMap<String, f64>) -> Option<f64> {
    match rule(counter)? {
        Rule::InGroupOps => Some(0.0),
        Rule::Units(us) => Some(
            us.iter()
                .map(|(u, k)| k * units[*u] * count as f64)
                .sum::<f64>(),
        ),
    }
}

/// Reprice a report of the base method in place.  `ic` reports are
/// rebuilt per phase (additions and doublings at one each, the
/// convention every ecbench row uses, plus every native counter); rho
/// reports keep their group-operation phases and gain a `natives` phase.
pub fn reprice(rep: &mut SolveReport, family: &str, slug: &str) -> Result<(), String> {
    let units = row(slug)?;
    let mut by_counter: BTreeMap<String, f64> = BTreeMap::new();
    let mut unknown: Vec<String> = Vec::new();
    let mut charge = |name: &str, count: u64, by: &mut BTreeMap<String, f64>| -> f64 {
        if count == 0 {
            return 0.0;
        }
        match price(name, count, &units) {
            Some(g) => {
                *by.entry(name.to_string()).or_default() += g;
                g
            }
            None => {
                let label = if name.ends_with("_uncharged") {
                    name.to_string()
                } else {
                    format!("{name}_uncharged")
                };
                if !unknown.contains(&label) {
                    unknown.push(label);
                }
                0.0
            }
        }
    };
    // Unpriced items that are not native counters (solver wall time,
    // counts lost at a cap) stay what they were.
    let mut kept: Vec<String> = rep
        .unpriced
        .iter()
        .filter(|u| u.starts_with("solver_") || u.starts_with("counts_lost"))
        .cloned()
        .collect();
    if family == "ic" {
        for p in rep.phases.iter_mut() {
            let mut g = (p.adds + p.doubles) as f64;
            for (name, &count) in &p.native {
                g += charge(name, count, &mut by_counter);
            }
            p.gae = g;
        }
        // The pipeline's shape counters (rows, rank, trials) are not work.
        rep.total_gae = rep.phases.iter().map(|p| p.gae).sum();
    } else {
        let mut natives = PhaseRecord {
            name: "natives".into(),
            ..Default::default()
        };
        let mut g = 0.0;
        for (name, &count) in &rep.counters {
            let before = by_counter.get(name).copied().unwrap_or(0.0);
            g += charge(name, count, &mut by_counter);
            if by_counter.get(name).copied().unwrap_or(0.0) > before {
                natives.native.insert(name.clone(), count);
            }
        }
        natives.gae = g;
        rep.total_gae += g;
        rep.phases.push(natives);
    }
    kept.extend(unknown);
    kept.sort();
    kept.dedup();
    rep.unpriced = kept;
    let table_sha = sha256_hex(FULL_PRICES.as_bytes());
    if let Value::Object(d) = &mut rep.detail {
        d.insert(
            "full_pricing".into(),
            json!({
                "table": "docs/ecbench/full_prices.json",
                "table_sha256": table_sha,
                "units": units,
                "gae_by_counter": by_counter,
            }),
        );
    } else {
        let old = std::mem::take(&mut rep.detail);
        rep.detail = json!({"base": old, "full_pricing": {
            "table": "docs/ecbench/full_prices.json",
            "table_sha256": table_sha,
            "units": units,
            "gae_by_counter": by_counter,
        }});
    }
    Ok(())
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// Nanoseconds per field multiplication, squaring and inversion, and per
/// probe of the pipeline's own hash map type ([`fast_map`]) holding `2^20`
/// entries, larger than any table of the panel and past the L2 cache, so
/// the probe is priced at or above what a run's table costs (half the
/// probes hit).
fn field_and_lookup_ns(inst: &BinaryInstance, seed: u64) -> BTreeMap<&'static str, f64> {
    let gf = &inst.gf;
    let mut rng = StdRng::seed_from_u64(seed);
    let mask = (1u64 << inst.n) - 1;
    let xs: Vec<u64> = (0..1024)
        .map(|_| (rng.gen::<u64>() & mask).max(1))
        .collect();
    let samples = 400_000u64;
    let mut out = BTreeMap::new();
    let start = Instant::now();
    let mut acc = xs[0];
    for i in 0..samples {
        acc = gf.mul(acc, xs[(i % 1024) as usize]) ^ 1;
    }
    out.insert(
        "field_mul",
        start.elapsed().as_nanos() as f64 / samples as f64,
    );
    std::hint::black_box(acc);
    let start = Instant::now();
    let mut acc = xs[1];
    for i in 0..samples {
        acc = gf.sqr(acc) ^ xs[(i % 1024) as usize];
    }
    out.insert(
        "field_sqr",
        start.elapsed().as_nanos() as f64 / samples as f64,
    );
    std::hint::black_box(acc);
    let inv_samples = samples / 4;
    let start = Instant::now();
    let mut acc = 0u64;
    for i in 0..inv_samples {
        acc ^= gf.inv(xs[(i % 1024) as usize] ^ (acc & 1)).max(1);
    }
    out.insert(
        "field_inv",
        start.elapsed().as_nanos() as f64 / inv_samples as f64,
    );
    std::hint::black_box(acc);
    let mut map: FastMap<u32> = fast_map(1 << 20);
    for i in 0..(1u32 << 20) {
        map.insert(rng.gen::<u64>(), i);
    }
    let keys: Vec<u64> = map.keys().copied().step_by(97).take(1 << 13).collect();
    let probes: Vec<u64> = (0..(1usize << 14))
        .map(|i| {
            if i % 2 == 0 {
                keys[(i / 2) % keys.len()]
            } else {
                rng.gen::<u64>()
            }
        })
        .collect();
    let start = Instant::now();
    let mut hits = 0u64;
    for i in 0..samples {
        if map.contains_key(&probes[(i as usize) % probes.len()]) {
            hits += 1;
        }
    }
    out.insert("lookup", start.elapsed().as_nanos() as f64 / samples as f64);
    std::hint::black_box(hits);
    out
}

/// Measure every unit for one curve `reps` times on this host and return
/// the median ratio to an addition per unit, with the per-rep ratios.
pub fn measure(inst: &BinaryInstance, reps: usize) -> Result<Value, String> {
    let mut ratios: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
    let mut add_ns = Vec::new();
    for k in 0..reps.max(1) {
        let c = calibrate_binary_instance(inst);
        let add = c.ns_per_add;
        if add <= 0.0 {
            return Err("addition measured at zero nanoseconds".into());
        }
        add_ns.push(add);
        let canon = c
            .ns_per_canon
            .ok_or("no Frobenius canonicaliser on this curve (not Koblitz?)")?;
        let mut push = |u: &'static str, ns: f64| ratios.entry(u).or_default().push(ns / add);
        push("double", c.ns_per_double);
        push("as_solve", c.ns_per_as_solve.unwrap_or(f64::NAN));
        push("row_op", c.ns_per_row_op);
        push("frobenius", c.ns_per_frobenius.unwrap_or(f64::NAN));
        push("canon", canon);
        push("word_xor", c.ns_per_word_xor.unwrap_or(f64::NAN));
        for (u, ns) in field_and_lookup_ns(inst, 0xF0 + k as u64) {
            push(u, ns);
        }
    }
    let mut row = serde_json::Map::new();
    let mut per_rep = serde_json::Map::new();
    for (u, v) in &ratios {
        if v.iter().any(|x| !x.is_finite()) {
            return Err(format!("unit {u} did not measure"));
        }
        row.insert((*u).into(), json!((median(v.clone()) * 1e6).round() / 1e6));
        per_rep.insert((*u).into(), json!(v));
    }
    Ok(json!({
        "row": row,
        "reps": reps,
        "ns_per_add_median": median(add_ns.clone()),
        "ns_per_add": add_ns,
        "ratios_per_rep": per_rep,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_full_id_has_a_base_and_the_table_parses() {
        for id in [
            "ic.pipeline_full",
            "rho.signed_frobenius_full",
            "rho.signed_frobenius_strong_full",
        ] {
            assert!(base_of(id).is_some());
        }
        let doc: Value = serde_json::from_str(FULL_PRICES).unwrap();
        for (slug, _) in doc["instances"].as_object().unwrap() {
            let r = row(slug).unwrap();
            assert_eq!(r.len(), UNITS.len());
            assert!(
                r.values().all(|v| *v > 0.0),
                "{slug}: a unit priced at zero"
            );
        }
    }

    #[test]
    fn an_unknown_counter_stays_unpriced() {
        let doc: Value = serde_json::from_str(FULL_PRICES).unwrap();
        let slug = doc["instances"]
            .as_object()
            .unwrap()
            .keys()
            .next()
            .unwrap()
            .clone();
        let mut rep = SolveReport {
            recovered: None,
            exhausted: false,
            phases: vec![PhaseRecord {
                name: "relations".into(),
                adds: 10,
                doubles: 2,
                native: [
                    ("lookups".to_string(), 4u64),
                    ("mystery_work".to_string(), 3),
                ]
                .into_iter()
                .collect(),
                ..Default::default()
            }],
            total_gae: 0.0,
            automorphisms_used: 2,
            counters: BTreeMap::new(),
            unpriced: vec!["lookups_uncharged".into()],
            deterministic: true,
            nondeterminism: vec![],
            solve_wall_ns: 0,
            factor_base: None,
            detail: json!({}),
            online: None,
            online_error: None,
            solver: None,
            field_ops: None,
        };
        reprice(&mut rep, "ic", &slug).unwrap();
        assert_eq!(rep.unpriced, vec!["mystery_work_uncharged".to_string()]);
        let lookup = row(&slug).unwrap()["lookup"];
        assert!((rep.total_gae - (12.0 + 4.0 * lookup)).abs() < 1e-9);
    }
}
