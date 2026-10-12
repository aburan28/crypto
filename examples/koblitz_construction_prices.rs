//! Ledger §21: what the Koblitz workflow's constructions cost, one by one.
//!
//! Ledger §20 put the workflow's constructions — everything it builds
//! around selection, the table, collection, the linear algebra and the
//! descent — at 7–70% of `S`, on clocks that lump several constructors
//! together.  This prices each constructor alone, on one thread, in the
//! thread's unit (one batched affine addition, `add_many` over 1,024
//! points) measured in the same process:
//!
//! - the curve (`KoblitzCurve::new`, the group order factorised);
//! - the base's projected signed-orbit map, as
//!   `projected_signed_orbit_count` builds it;
//! - the point index (`FrobeniusFactorBase::index_map`);
//! - the decomposition check (`m_can_decompose`);
//! - the field structure (`FieldStructure::new`);
//! - each consumer the workflow builds: `ColumnCoverage::new`,
//!   `RelationCollector::with_pair_table`, `FactorBaseLogSolver::new`.
//!
//! Each is timed `rounds` times and reported as the median in units, and
//! per base point.  The first call of a cached construction and a later
//! one are timed apart, so a cache's effect is visible in one run.
//!
//!     cargo run --release --example koblitz_construction_prices -- <params.json> [rounds]
use crypto_lib::cryptanalysis::koblitz_fast::{BatchScratch, FastCurve, FastPoint};
use crypto_lib::cryptanalysis::koblitz_groebner::FieldStructure;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_subgroup_orbit_factor_base_with_cost, projected_signed_orbit_count, ColumnCoverage,
    DecompositionStrategy, FactorBaseLogSolver, KoblitzCurve, KoblitzIcOptions, PairSumTable,
    RelationCollector,
};
use serde_json::{json, Value};
use std::time::Instant;

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

fn main() {
    let path = std::env::args().nth(1).expect("a workflow parameter file");
    let rounds: usize = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);
    let p: Value =
        serde_json::from_str(&std::fs::read_to_string(&path).expect("params")).expect("json");
    let n = p["curve"]["degree"].as_u64().expect("degree") as u32;
    let a = p["curve"]["curve_a"].as_u64().unwrap_or(0) as u8;
    let spec = &p["factor_base"]["spec"];
    let points = spec["points"].as_u64().expect("points") as usize;
    let seed = spec["seed"].as_u64().expect("seed");
    let m = p["summands"].as_u64().unwrap_or(3) as usize;

    let kc = KoblitzCurve::new(a, n).expect("a Koblitz curve");
    let fc = FastCurve::new(&kc.curve).expect("single-word arithmetic");
    let (fb, _) = build_subgroup_orbit_factor_base_with_cost(&kc, seed, points).expect("a base");
    let pair = PairSumTable::build_folded_within(&kc, &fb, PairSumTable::DEFAULT_BYTE_BUDGET)
        .expect("a folded table");
    let opts = KoblitzIcOptions {
        m,
        strategy: DecompositionStrategy::PairTable,
        collapse_negation: true,
        collapse_projected_orbits: true,
        allow_direct_relation: false,
        ..KoblitzIcOptions::default()
    };

    // The unit.
    let g = fc.lift(kc.generator());
    let batch: Vec<FastPoint> = (0..1024u64)
        .map(|i| fc.mul(g, &num_bigint::BigUint::from(1 + i * 7919)))
        .collect();
    let mut scratch = BatchScratch::default();
    let mut out = Vec::with_capacity(1024);
    let mut unit = || {
        let t = Instant::now();
        for _ in 0..1000 {
            out.clear();
            fc.add_many(g, &batch, &mut out, &mut scratch);
        }
        t.elapsed().as_nanos() as f64 / 1_024_000.0
    };
    unit();

    // Run with RAYON_NUM_THREADS=1: the constructors' parallel loops are
    // priced as one thread runs them.
    let mut report = serde_json::Map::new();
    let mut time = |name: &str, f: &mut dyn FnMut()| {
        let mut units = Vec::with_capacity(rounds);
        let mut first = None;
        for _ in 0..rounds {
            let before = unit();
            let t = Instant::now();
            f();
            let ns = t.elapsed().as_nanos() as f64;
            let after = unit();
            let u = ns / ((before + after) / 2.0);
            first.get_or_insert(u);
            units.push(u);
        }
        let med = median(&mut units);
        report.insert(
            name.to_string(),
            json!({"median_units": med, "first_units": first, "per_point": med / fb.points.len() as f64}),
        );
    };
    time("curve", &mut || {
        std::hint::black_box(KoblitzCurve::new(a, n));
    });
    time("orbit_map", &mut || {
        std::hint::black_box(projected_signed_orbit_count(&kc, &fb));
    });
    time("index_map", &mut || {
        std::hint::black_box(fb.index_map());
    });
    time("m_can_decompose", &mut || {
        std::hint::black_box(fb.m_can_decompose(&kc, m));
    });
    time("field_structure", &mut || {
        std::hint::black_box(FieldStructure::new(kc.n, &kc.curve.irreducible));
    });
    time("column_coverage", &mut || {
        std::hint::black_box(ColumnCoverage::new(&kc, &fb));
    });
    time("relation_collector", &mut || {
        std::hint::black_box(
            RelationCollector::with_pair_table(&kc, &fb, &opts, Some(&pair)).is_some(),
        );
    });
    time("log_solver", &mut || {
        std::hint::black_box(FactorBaseLogSolver::new(&kc, &fb, &opts).is_some());
    });
    let out = json!({
        "params": path, "n": n, "a": a, "points": fb.points.len(), "signed_orbits": fb.signed_orbits.len(),
        "cofactor_bits": kc.cofactor.bits(), "rounds": rounds,
        "constructors": Value::Object(report),
    });
    println!("{}", serde_json::to_string_pretty(&out).expect("json"));
}
