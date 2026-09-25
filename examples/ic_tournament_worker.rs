//! Complete, public-synthetic IC/rho jobs for the candidate tournament.
//! Uses the production research library; JSON stdin contains no known scalar.
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_factor_base_search::FactorBaseSpec;
use crypto_lib::cryptanalysis::koblitz_groebner::SolverEngine;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    koblitz_signed_frobenius_rho_with_preparation, point_key, points_with_x, DecompositionStrategy,
    FactorBaseLogSolver, IndividualLogSolver, KoblitzCurve, KoblitzIcOptions,
    KoblitzSignedRhoOptions, LinearAlgebra, PairSumTable, RelationCollector, RelationWorkUnit,
};
use crypto_lib::cryptanalysis::koblitz_sparse_la::SparseSolveOptions;
use crypto_lib::cryptanalysis::semaev_sat::XorEncoding;
use num_bigint::BigUint;
use serde::Deserialize;
use serde_json::{json, Value};
use std::io::Read;
use std::time::Instant;

/// The amd64 Valgrind client-request ABI, as specified in valgrind.h.
/// Outside Valgrind the rotation preamble and self-exchange have no effect.
/// DUMP_STATS_AT writes the interval then resets counters; collection stays on.
#[inline(never)]
fn dump(label: &'static [u8]) {
    assert_eq!(label.last(), Some(&0));
    #[cfg(target_arch = "x86_64")]
    unsafe {
        let args = [0x4354_0003usize, label.as_ptr() as usize, 0, 0, 0, 0];
        std::arch::asm!(
            "rol rdi, 3", "rol rdi, 13", "rol rdi, 61", "rol rdi, 51",
            "xchg rbx, rbx",
            in("rax") args.as_ptr(), inout("rdx") 0usize => _, out("rdi") _,
            options(nostack)
        );
    }
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Job {
    mode: String,
    degree: u32,
    curve_a: u8,
    #[serde(default)]
    target_seeds: Vec<u64>,
    #[serde(default)]
    public_targets: Option<Vec<[String; 2]>>,
    algorithm_seed: u64,
    factor_base: FactorBaseSpec,
    config: Config,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, deny_unknown_fields)]
struct Config {
    solver: String,
    linear_algebra: String,
    batch_trials: u64,
    max_trials: u64,
    summands: usize,
    collection_window: Option<usize>,
    sparse: SparseSolveOptions,
    factor_base_orbits: Option<usize>,
    factor_base_cube_root: bool,
    factor_base: Option<FactorBaseSpec>,
    groebner_degree: u32,
    node_budget: usize,
    conflict_budget: u64,
    rho_parallel_walks: usize,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            solver: "pair_table".into(),
            linear_algebra: "sparse".into(),
            batch_trials: 64,
            max_trials: 4096,
            summands: 3,
            collection_window: None,
            sparse: SparseSolveOptions::default(),
            factor_base_orbits: None,
            factor_base_cube_root: false,
            factor_base: None,
            groebner_degree: 3,
            node_budget: 4096,
            conflict_budget: 100_000,
            rho_parallel_walks: 32,
        }
    }
}

fn point(p: &BinaryPoint) -> Value {
    match p {
        BinaryPoint::Infinity => Value::Null,
        BinaryPoint::Affine { x, y } => {
            json!([x.to_biguint().to_string(), y.to_biguint().to_string()])
        }
    }
}

// None means no valid-range scalar was available to replay; this is not a
// successful check and must not claim that scalar replay was performed.
fn scalar_replay(c: &KoblitzCurve, q: &BinaryPoint, log: Option<&BigUint>) -> Option<bool> {
    log.filter(|d| *d < &c.subgroup_order)
        .map(|d| c.mul(c.generator(), d) == *q)
}

// Strict public coordinates: do not silently reduce or repair supplied inputs.
fn public_point(c: &KoblitzCurve, encoded: &[String; 2]) -> Result<BinaryPoint, String> {
    let coordinate = |text: &str| -> Result<F2mElement, String> {
        if text.is_empty()
            || !text.bytes().all(|b| b.is_ascii_digit())
            || (text.len() > 1 && text.starts_with('0'))
        {
            return Err("noncanonical public coordinate".into());
        }
        let value = text
            .parse::<BigUint>()
            .map_err(|_| "invalid public coordinate")?;
        if value.bits() > u64::from(c.n) {
            return Err("public coordinate outside field".into());
        }
        Ok(F2mElement::from_biguint(&value, c.n))
    };
    let q = BinaryPoint::Affine {
        x: coordinate(&encoded[0])?,
        y: coordinate(&encoded[1])?,
    };
    if !c.curve.is_on_curve(&q) || c.mul(&q, &c.subgroup_order) != BinaryPoint::Infinity {
        return Err("public point is outside the declared subgroup".into());
    }
    Ok(q)
}

// Same public hash-to-curve domain as ic workflow. No target scalar is created.
fn target(c: &KoblitzCurve, seed: u64) -> Result<BinaryPoint, String> {
    for counter in 0u64..1_000_000 {
        let mut h = blake3::Hasher::new();
        h.update(b"ic-workflow-public-target-v1\0");
        h.update(&c.n.to_le_bytes());
        h.update(&[c.a]);
        h.update(&c.k.to_le_bytes());
        h.update(&c.b_index.to_le_bytes());
        h.update(&seed.to_le_bytes());
        h.update(&counter.to_le_bytes());
        let digest = h.finalize();
        let x =
            u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap()) & ((1u64 << c.n) - 1);
        let x = F2mElement::from_biguint(&BigUint::from(x), c.n);
        let mut lifts = points_with_x(&c.curve, &x);
        lifts.sort_by_key(point_key);
        if lifts.is_empty() {
            continue;
        }
        let q = c.mul(
            &lifts[usize::from(digest.as_bytes()[8] & 1) % lifts.len()],
            &c.cofactor,
        );
        if q != BinaryPoint::Infinity {
            return Ok(q);
        }
    }
    Err("hash-to-curve exhausted".into())
}

fn run(job: &Job) -> Result<Value, String> {
    if !(5..=31).contains(&job.degree) || job.degree.is_multiple_of(2) || job.curve_a > 1 {
        return Err("worker accepts bounded odd-degree Koblitz fixtures (5..31)".into());
    }
    if job.mode != "fixture" && job.public_targets.as_ref().map(Vec::len) != Some(1) {
        return Err(
            "measured jobs require exactly one supplied public point; use fixture mode first"
                .into(),
        );
    }
    let cfg = &job.config;
    if !(1..=256).contains(&cfg.rho_parallel_walks)
        || !(2..=4).contains(&cfg.groebner_degree)
        || cfg.node_budget == 0
        || cfg.conflict_budget == 0
    {
        return Err("invalid solver or rho limits".into());
    }
    if cfg.batch_trials == 0
        || cfg.batch_trials > 4096
        || cfg.max_trials > 65536
        || cfg.max_trials < cfg.batch_trials
        || !(2..=4).contains(&cfg.summands)
    {
        return Err("invalid collection limits".into());
    }
    let start = Instant::now();
    let c = KoblitzCurve::new(job.curve_a, job.degree).ok_or("curve has no usable subgroup")?;
    let targets = match &job.public_targets {
        Some(encoded) => {
            if encoded.is_empty()
                || encoded.len() > 100
                || (!job.target_seeds.is_empty() && job.target_seeds.len() != encoded.len())
            {
                return Err("invalid public-point/seed counts".into());
            }
            encoded
                .iter()
                .map(|p| public_point(&c, p))
                .collect::<Result<Vec<_>, _>>()?
        }
        None => {
            if job.mode != "fixture" || job.target_seeds.is_empty() || job.target_seeds.len() > 100
            {
                return Err("fixture mode requires 1..100 target seeds".into());
            }
            job.target_seeds
                .iter()
                .map(|&s| target(&c, s))
                .collect::<Result<Vec<_>, _>>()?
        }
    };
    let target_seeds: Vec<Option<u64>> = if job.target_seeds.is_empty() {
        vec![None; targets.len()]
    } else {
        job.target_seeds.iter().copied().map(Some).collect()
    };
    let metadata = json!({
        "degree": c.n, "curve_a": c.a,
        "irreducible": {"degree": c.curve.irreducible.degree, "low_terms": c.curve.irreducible.low_terms},
        "subgroup_order": c.subgroup_order.to_string(), "group_order": c.group_order.to_string(),
        "cofactor": c.cofactor.to_string(), "lambda": c.lambda.to_string(),
        "generator": point(c.generator()), "targets": targets.iter().map(point).collect::<Vec<_>>(),
        "target_seeds": target_seeds, "target_scalar_constructed": false
    });
    dump(b"curve_and_targets\0");
    if job.mode == "fixture" {
        return Ok(json!({"schema_version":1,"status":"fixture","fixture":metadata}));
    }
    if job.mode == "rho" {
        let q = &targets[0];
        let options = KoblitzSignedRhoOptions {
            seed: job.algorithm_seed,
            max_iterations_per_restart: cfg.max_trials,
            parallel_walks: cfg.rho_parallel_walks,
            ..Default::default()
        };
        let mut online_start = None;
        let mut field_kernel = None;
        let answer = koblitz_signed_frobenius_rho_with_preparation(
            &c,
            q,
            &options,
            &mut |fc| {
                field_kernel = Some(fc.field.kernel_name());
                online_start = Some(Instant::now());
            },
            &mut |_| {},
        )
        .ok_or("rho has no prepared packed backend")?;
        dump(b"rho_solve\0");
        let replay = scalar_replay(&c, q, answer.recovered_log.as_ref());
        let verified = replay == Some(true);
        let online_ns = online_start
            .expect("prepared rho boundary")
            .elapsed()
            .as_nanos();
        dump(b"final_verification\0");
        let solutions = vec![
            json!({"index":0,"recovered":answer.recovered_log.as_ref().map(ToString::to_string),
            "verified":answer.verified,"iterations":answer.iterations,"restarts":answer.restarts_attempted,
            "walk_group_additions":answer.charges.walk_group_additions,"effective_walks":answer.parallel_walks}),
        ];
        return Ok(
            json!({"schema_version":1,"mode":"rho","status":if verified {"complete"} else {"incomplete"},
            "fixture":metadata,"solutions":solutions,"automorphism_order":2*c.n,
            "elapsed_seconds":start.elapsed().as_secs_f64(),"online_timing_schema":1,
            "online_wall_ns":online_ns,"target_input":"supplied_public_point",
            "reusable_setup_excluded":true,"scalar_replay_included":replay.is_some(),"field_kernel":field_kernel}),
        );
    }
    if job.mode != "ic" {
        return Err("unknown mode".into());
    }
    let strategy = match cfg.solver.as_str() {
        "pair_table" => DecompositionStrategy::PairTable,
        "enumerate" => DecompositionStrategy::Enumerate,
        "f4" | "f5" | "inherited_f4" => DecompositionStrategy::Groebner,
        "sat_xor" | "sat_cnf" => DecompositionStrategy::Sat,
        _ => return Err("unsupported decomposition backend".into()),
    };
    let la = match cfg.linear_algebra.as_str() {
        "dense" => LinearAlgebra::Dense,
        "sparse" => LinearAlgebra::Sparse(cfg.sparse),
        _ => return Err("unknown linear algebra mode".into()),
    };
    let mut opts = KoblitzIcOptions {
        m: cfg.summands,
        seed: job.algorithm_seed,
        max_trials: cfg.max_trials as usize,
        strategy,
        linear_algebra: la,
        node_budget: cfg.node_budget,
        collection_window: cfg.collection_window,
        allow_direct_relation: false,
        ..KoblitzIcOptions::default()
    };
    opts.engine = match cfg.solver.as_str() {
        "f5" => SolverEngine::MatrixF5 {
            max_degree: cfg.groebner_degree,
        },
        "inherited_f4" => SolverEngine::InheritedF4 {
            max_degree: cfg.groebner_degree,
        },
        _ => SolverEngine::MatrixF4 {
            max_degree: cfg.groebner_degree,
        },
    };
    opts.sat_options.conflict_budget = cfg.conflict_budget;
    opts.sat_options.encoding = if cfg.solver == "sat_cnf" {
        XorEncoding::Cnf
    } else {
        XorEncoding::Native
    };
    if cfg
        .factor_base_orbits
        .is_some_and(|n| !(1..=8).contains(&n))
        || (cfg.factor_base_orbits.is_some() && cfg.factor_base_cube_root)
        || (cfg.factor_base.is_some()
            && (cfg.factor_base_orbits.is_some() || cfg.factor_base_cube_root))
    {
        return Err("invalid factor-base policy".into());
    }
    let effective_base = if cfg.factor_base_orbits.is_some() || cfg.factor_base_cube_root {
        let FactorBaseSpec::SubgroupOrbits { seed, .. } = &job.factor_base else {
            return Err("orbit policy requires a subgroup-orbit base recipe".into());
        };
        let points = if let Some(orbits) = cfg.factor_base_orbits {
            2 * c.n as usize * orbits
        } else {
            let r = c.subgroup_order.to_u64_digits()[0];
            let mut b = 1u64;
            while b * b * b < r / 2 {
                b += 1;
            }
            (b as usize).max(2 * c.n as usize)
        };
        FactorBaseSpec::SubgroupOrbits {
            seed: *seed,
            points,
        }
    } else {
        cfg.factor_base
            .clone()
            .unwrap_or_else(|| job.factor_base.clone())
    };
    let fb = effective_base.materialize(&c)?;
    let pair = if strategy == DecompositionStrategy::PairTable {
        Some(PairSumTable::build(&c, &fb).ok_or("pair table unavailable")?)
    } else {
        None
    };
    let collector = RelationCollector::with_pair_table(&c, &fb, &opts, pair.as_ref())
        .ok_or("unsupported decomposition")?;
    let mut system = FactorBaseLogSolver::new(&c, &fb, &opts).ok_or("no projected columns")?;
    let expected_columns = system.columns();
    dump(b"factor_base_and_tables\0");
    let mut relations = Vec::new();
    let mut trials = 0;
    let mut collection_reports = Vec::new();
    let mut outcome = None;
    while trials < cfg.max_trials && outcome.is_none() {
        let count = cfg.batch_trials.min(cfg.max_trials - trials);
        let (rows, collection_report) = collector.collect_observed(RelationWorkUnit {
            seed: job.algorithm_seed,
            start: trials,
            count,
        });
        trials += count;
        dump(b"collection_and_decomposition\0");
        system.push(&rows);
        collection_reports.push(collection_report);
        outcome = system.try_solve();
        relations.extend(rows);
        dump(b"verify_filter_and_linear_algebra\0");
    }
    let base = fb.points.iter().map(point).collect::<Vec<_>>();
    let mut report = system.report();
    report.trials = trials as usize;
    let Some((table, mut solved)) = outcome else {
        return Ok(
            json!({"schema_version":1,"mode":"ic","status":"incomplete","fixture":metadata,
            "trials":trials,"columns":expected_columns,"relations":relations,
            "factor_base":base,"rejected_relations":report.rejected_relations,
            "query_schema_version":1,"collection_reports":collection_reports,
            "log_table_report":report,"solve_attempts":report.solve_attempts,
            "effective_factor_base":effective_base,"summands":cfg.summands,
            "online_timing_schema":1,"online_wall_ns":null,"target_input":"supplied_public_point",
            "reusable_setup_excluded":true,"scalar_replay_included":false}),
        );
    };
    solved.trials = trials as usize;
    if solved.rejected_relations != 0 || !table.verify(&c) {
        return Ok(
            json!({"schema_version":1,"mode":"ic","status":"invalid_certificate",
            "fixture":metadata,"factor_base":base,"relations":relations,"trials":trials,
            "query_schema_version":1,"collection_reports":collection_reports,
            "log_table_report":solved,"reason":"invalid relation or factor-base log certificate"}),
        );
    }
    let columns = table
        .columns
        .iter()
        .map(|(p, l)| json!({"point":point(p),"log":l.to_string()}))
        .collect::<Vec<_>>();
    dump(b"log_certification\0");
    let solver = IndividualLogSolver::new(&c, &fb, &table, &opts, pair.as_ref())
        .ok_or("descent setup failed")?;
    let online_start = Instant::now();
    let q = &targets[0];
    let answer = solver.solve_observed(q);
    dump(b"individual_log\0");
    let replay = scalar_replay(&c, q, answer.log.as_ref());
    let verified = replay == Some(true);
    let online_ns = online_start.elapsed().as_nanos();
    dump(b"final_verification\0");
    let solutions = vec![
        json!({"index":0,"recovered":answer.log.as_ref().map(ToString::to_string),
        "trials":answer.trials,"relation":answer.relation,"attempts":answer.attempts}),
    ];
    Ok(
        json!({"schema_version":1,"mode":"ic","status":if verified{"complete"}else{"incomplete"},
        "fixture":metadata,"solutions":solutions,"factor_base":base,"relations":relations,
        "column_logs":columns,"columns":expected_columns,"trials":trials,"solve_attempts":solved.solve_attempts,
        "query_schema_version":1,"collection_reports":collection_reports,"log_table_report":solved,
        "online_timing_schema":1,"online_wall_ns":online_ns,"target_input":"supplied_public_point",
        "reusable_setup_excluded":true,"scalar_replay_included":replay.is_some(),
        "accepted_relations":solved.relations,"duplicate_relations":solved.duplicate_relations,
        "rejected_relations":solved.rejected_relations,"summands":cfg.summands,
        "effective_factor_base":effective_base,
        "sparse_report":solved.sparse_report,"elapsed_seconds":start.elapsed().as_secs_f64()}),
    )
}

fn main() {
    let mut input = String::new();
    std::io::stdin()
        .take(1_048_577)
        .read_to_string(&mut input)
        .expect("stdin");
    let result = if input.len() > 1_048_576 {
        Err("job too large".into())
    } else {
        serde_json::from_str::<Job>(&input)
            .map_err(|e| e.to_string())
            .and_then(|job| {
                dump(b"startup_and_input\0");
                run(&job)
            })
    };
    let report = result
        .unwrap_or_else(|reason| json!({"schema_version":1,"status":"error","reason":reason}));
    let success = matches!(report["status"].as_str(), Some("complete" | "fixture"));
    println!("{}", report);
    if !success {
        std::process::exit(2);
    }
}

#[cfg(test)]
mod public_input_tests {
    use super::*;

    fn job(mode: &str) -> Job {
        serde_json::from_value(json!({"mode":mode,"degree":9,"curve_a":0,
            "target_seeds":[2026092555],"algorithm_seed":2026092555,
            "factor_base":{"kind":"factor","index":0},
            "config":{"solver":"pair_table","linear_algebra":"dense","summands":2,
                "batch_trials":8,"max_trials":256}}))
        .unwrap()
    }

    fn supplied(mode: &str) -> Job {
        let fixture = run(&job("fixture")).unwrap();
        let mut measured = job(mode);
        measured.public_targets =
            Some(serde_json::from_value(fixture["fixture"]["targets"].clone()).unwrap());
        measured
    }

    #[test]
    fn public_input_rejects_seed_only_batch_and_noncanonical_points() {
        for mode in ["ic", "rho"] {
            assert!(run(&job(mode))
                .unwrap_err()
                .contains("exactly one supplied"));
            let original = supplied(mode);
            let mut batch = original.clone();
            let point = batch.public_targets.as_ref().unwrap()[0].clone();
            batch.public_targets.as_mut().unwrap().push(point);
            assert!(run(&batch).unwrap_err().contains("exactly one supplied"));
            for bad in ["", "00", "01", "+1", "-1", " 1", "512"] {
                let mut changed = original.clone();
                changed.public_targets.as_mut().unwrap()[0][0] = bad.into();
                assert!(run(&changed).is_err(), "accepted {bad:?}");
            }
            let mut changed = original.clone();
            changed.public_targets = Some(vec![["0".into(), "0".into()]]); // off curve
            assert!(run(&changed)
                .unwrap_err()
                .contains("outside the declared subgroup"));
            changed.public_targets = Some(vec![["0".into(), "1".into()]]); // order-two point
            assert!(run(&changed)
                .unwrap_err()
                .contains("outside the declared subgroup"));
        }
    }

    #[test]
    fn public_input_seeds_are_provenance_and_online_interval_includes_replay() {
        for mode in ["ic", "rho"] {
            let original = supplied(mode);
            let first = run(&original).unwrap();
            assert_eq!(first["status"], "complete");
            let online = first["online_wall_ns"].as_u64().unwrap();
            assert!(online > 0);
            assert!(online as f64 <= first["elapsed_seconds"].as_f64().unwrap() * 1e9);
            assert_eq!(first["target_input"], "supplied_public_point");
            assert_eq!(first["reusable_setup_excluded"], true);
            assert_eq!(first["scalar_replay_included"], true);
            if mode == "rho" {
                assert!(matches!(
                    first["field_kernel"].as_str(),
                    Some("portable" | "pclmulqdq" | "pmull")
                ));
            }
            for seeds in [Vec::new(), vec![1]] {
                let mut changed = original.clone();
                changed.target_seeds = seeds;
                let after = run(&changed).unwrap();
                assert_eq!(after["fixture"]["targets"], first["fixture"]["targets"]);
                assert_eq!(after["solutions"], first["solutions"]);
            }
        }
    }

    #[test]
    fn public_input_absent_or_out_of_range_scalar_does_not_claim_replay() {
        let c = KoblitzCurve::new(0, 9).unwrap();
        let q = c.generator();
        assert_eq!(scalar_replay(&c, q, None), None);
        assert_eq!(scalar_replay(&c, q, Some(&c.subgroup_order)), None);
        assert_eq!(scalar_replay(&c, q, Some(&BigUint::from(1u32))), Some(true));
        assert_eq!(
            scalar_replay(&c, q, Some(&BigUint::from(2u32))),
            Some(false)
        );
    }

    #[test]
    fn public_input_failed_preparation_has_no_online_interval() {
        let mut job = supplied("ic");
        job.config.max_trials = 1;
        job.config.batch_trials = 1;
        let report = run(&job).unwrap();
        assert_eq!(report["status"], "incomplete");
        assert_eq!(report["trials"], 1);
        assert!(report["online_wall_ns"].is_null());
        assert_eq!(report["scalar_replay_included"], false);
        assert_eq!(
            report["collection_reports"][0]["attempts"]
                .as_array()
                .unwrap()
                .len(),
            1
        );
    }
}
