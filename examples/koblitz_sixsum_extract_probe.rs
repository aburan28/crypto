//! Higher-arity relation extraction over the Frobenius pair table: the
//! 4-sum baseline, kernel A (iterate one table), kernel B (S4 per probe,
//! row order), kernel B2 (S4 per probe, diagonal stride) and the 5-sum
//! control, measured per cell `(n, K, kernel)` in `S₃`/`S₄` solves,
//! first-summand probes, targets and wall per verified relation, against
//! the models `r/(n²K²)` and `45·r/(n⁴K⁴)` and the structural
//! dense/sparse expectations.  Prints one JSON document.
//!
//! ```text
//! cargo run --release --example koblitz_sixsum_extract_probe -- \
//!     --n=23,31,41 --a=0 --columns=5,10,20 --kernels=4sum,a,b,b2,c \
//!     --relations=32 --seed=20261012 --max-targets=4096 \
//!     --max-s3=400000000 --max-s4=50000000 --cell-wall-s=300
//! ```
//!
//! Public synthetic fixtures only (`[a]G` targets on the registered Koblitz
//! rungs); every relation counted has been re-added as a group equation.
//! No key-recovery, sub-rho or deployed-curve claim is made or implied.

use std::collections::BTreeMap;
use std::time::{Duration, Instant};

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use crypto_lib::cryptanalysis::ic_boundary::{
    koblitz_instance, BinaryGroup, BinaryInstance, CountedGroup, FactorBase, GroupOps,
};
use crypto_lib::cryptanalysis::ic_framework::plugins::CompactOrbitScanBase;
use crypto_lib::cryptanalysis::ic_framework::stages::{FactorBaseBuilder, InstanceCtx, Params};
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::cryptanalysis::koblitz_sixsum_extract::{
    Budget, BuildStats, Kernel, ProbeCounters, ScanEnd, SixSumTable,
};

struct Args {
    ns: Vec<u32>,
    a: u8,
    columns: Vec<usize>,
    kernels: Vec<Kernel>,
    relations: u64,
    seed: u64,
    max_targets: u64,
    max_s3: u64,
    max_s4: u64,
    cell_wall: Duration,
    early_abort: bool,
    raw_x_cap: u64,
    keep_relations: usize,
}

fn parse_list<T: std::str::FromStr>(raw: &str, what: &str) -> Vec<T> {
    raw.split(',')
        .filter(|s| !s.trim().is_empty())
        .map(|s| {
            s.trim()
                .parse::<T>()
                .unwrap_or_else(|_| panic!("bad {what} entry `{s}`"))
        })
        .collect()
}

fn parse_args() -> Args {
    let mut args = Args {
        ns: vec![23, 31, 41],
        a: 0,
        columns: vec![5, 10, 20],
        kernels: vec![
            Kernel::FourSum,
            Kernel::IterateTable,
            Kernel::S4PerProbe,
            Kernel::S4Strided,
            Kernel::FiveSumControl,
        ],
        relations: 32,
        seed: 20_261_012,
        max_targets: 4096,
        max_s3: 400_000_000,
        max_s4: 50_000_000,
        cell_wall: Duration::from_secs(300),
        early_abort: true,
        raw_x_cap: 100_000_000,
        keep_relations: 4,
    };
    for arg in std::env::args().skip(1) {
        let (key, value) = arg
            .trim_start_matches("--")
            .split_once('=')
            .unwrap_or((arg.trim_start_matches("--"), "1"));
        match key {
            "n" => args.ns = parse_list(value, "n"),
            "a" => args.a = value.parse().expect("a is 0 or 1"),
            "columns" | "k" => args.columns = parse_list(value, "columns"),
            "kernels" | "kernel" => {
                args.kernels = value
                    .split(',')
                    .filter(|s| !s.trim().is_empty())
                    .map(|s| Kernel::parse(s).unwrap_or_else(|e| panic!("{e}")))
                    .collect()
            }
            "relations" => args.relations = value.parse().expect("relations"),
            "seed" => args.seed = value.parse().expect("seed"),
            "max-targets" => args.max_targets = value.parse().expect("max-targets"),
            "max-s3" => args.max_s3 = value.parse().expect("max-s3"),
            "max-s4" => args.max_s4 = value.parse().expect("max-s4"),
            "cell-wall-s" => {
                args.cell_wall = Duration::from_secs_f64(value.parse().expect("cell-wall-s"))
            }
            "early-abort" => args.early_abort = matches!(value, "1" | "true" | "yes"),
            "raw-x-cap" => args.raw_x_cap = value.parse().expect("raw-x-cap"),
            "keep-relations" => args.keep_relations = value.parse().expect("keep-relations"),
            "help" | "h" => {
                eprintln!(
                    "--n=23,31,41 --a=0 --columns=5,10,20 --kernels=4sum,a,b,b2,c --relations=32 --seed=S --max-targets=N --max-s3=N --max-s4=N --cell-wall-s=S --early-abort=1 --raw-x-cap=N --keep-relations=N"
                );
                std::process::exit(0);
            }
            other => panic!("unknown argument `{other}`"),
        }
    }
    args
}

#[derive(Serialize, Clone, Copy)]
struct Models {
    /// `r / (n² K²)`: the feasibility note's 4-sum fit (C = 1).
    r_over_n2k2: f64,
    /// `r / (4 n² K²)`: the structural per-`S₃`-solve estimate with
    /// `K²n` table states — two roots per solve, each landing in the
    /// `K²n²` covered abscissae out of `r/2`.
    r_over_4n2k2: f64,
    /// `45 r / (n⁴ K⁴)`: the note's §4.7 6-sum model.
    model_45r_over_n4k4: f64,
    /// `r / (4 · abscissa_states)`: the same structural estimate from
    /// the table actually built (`abscissa_states = states × n`).
    structural_s3_solves_per_relation: f64,
    /// `r / (8 · abscissa_states)`: kernel B's quartic has four roots.
    structural_s4_solves_per_relation: f64,
}

#[derive(Serialize)]
struct SupportHistogram(BTreeMap<usize, u64>);

#[derive(Serialize)]
struct Cell {
    n: u32,
    a: u8,
    columns: usize,
    kernel: &'static str,
    arity: usize,
    r: u64,
    cofactor: u64,
    base_points: usize,
    table_states: usize,
    abscissa_states: usize,
    early_abort: bool,
    counters: ProbeCounters,
    relations: u64,
    /// Relations whose support after cancellation equals the arity.
    full_support_relations: u64,
    support_histogram: SupportHistogram,
    s3_solves_per_relation: Option<f64>,
    s4_solves_per_relation: Option<f64>,
    /// `S₃` solves plus `S₄` solves: one count of root solves.
    solves_per_relation: Option<f64>,
    targets_per_relation: Option<f64>,
    wall_s_per_relation: Option<f64>,
    point_adds_per_relation: Option<f64>,
    lookups_per_relation: Option<f64>,
    /// The note's unit: first-summand probes per relation (`S₃` solves
    /// for the 4-sum, outer states `s₁` for kernel A, base points for
    /// kernel C, `S₄` solves for kernel B).
    outer_probes_per_relation: Option<f64>,
    /// Structural expectation of `outer_probes_per_relation` in the
    /// **dense regime** (several relations per target, so the first
    /// hit lands inside one target) from the table actually built
    /// (`A` = abscissa states): `r/(4A)` for the 4-sum, `r/(4A²)` for
    /// kernels A and C (one inner scan per outer probe, hit probability
    /// `4A²/r`), `r/(8A)` for kernel B.
    structural_outer_probes_per_relation: f64,
    /// Expected relations of this arity per target: `B^m/(m!·r)` for
    /// `B` base points and arity `m` (unordered `m`-sums).
    expected_relations_per_target: f64,
    /// Structural expectation of targets per relation,
    /// `1/(1 − exp(−expected_relations_per_target))`.
    structural_targets_per_relation: f64,
    /// Solves in one *hit-or-miss scan unit*: the unit whose hits are
    /// clustered together.  `A` (`S₃`) for the 4-sum (one target
    /// scan), `A` for kernels A and C (one inner 4-sum scan of `R − s₁`
    /// or `R − P_e`), `A` (`S₄`) for kernel B (one row, `s₁` fixed,
    /// which is the 4-sum scan of `R ± s₁`), `A²/2` for kernel B2 (the
    /// whole pair square, de-clustered by the diagonal walk).
    scan_unit_solves: f64,
    /// Expected relations in one scan unit: `λ₄ = B⁴/(24r)` for the
    /// 4-sum and kernels A, C; `2λ₄` for kernel B's row; `B⁶/(720r)`
    /// for kernel B2.
    scan_unit_lambda: f64,
    /// Dense-regime structural solves per relation: `r/(4A)` (`S₃`) or
    /// `r/(8A)` (`S₄`): the first hit of a scan unit with several.
    structural_dense_solves_per_relation: f64,
    /// Sparse-regime structural solves per relation (fewer than one
    /// relation per scan unit, so most units run to completion):
    /// `scan_unit_solves / scan_unit_lambda`.
    structural_sparse_solves_per_relation: f64,
    /// Hits of one relation inside its scan unit: 6 for every kernel
    /// whose unit is a 4-sum scan (three pairings, two states each),
    /// 45 for kernel B2 (the unordered pair choices of a 6-sum).
    scan_unit_multiplicity: f64,
    /// The interpolating structural expectation,
    /// `scan_unit_solves × (e^{−λ}/(1 − e^{−λ}) + 1/(μλ + 1))`: empty
    /// units paid in full plus the first of `μλ` uniformly placed hits.
    /// Tends to `sparse` for `λ ≪ 1` and to `dense` for `λ ≫ 1`.
    structural_solves_per_relation: f64,
    models: Models,
    /// Measured solves per relation over `r/(n²K²)`.
    fit_c_solves: Option<f64>,
    /// Measured targets per relation over `45r/(n⁴K⁴)`.
    targets_over_model_45: Option<f64>,
    /// Measured outer probes per relation over `45r/(n⁴K⁴)` (the
    /// note's own unit).
    outer_probes_over_model_45: Option<f64>,
    scan_ends: BTreeMap<String, u64>,
    stop_reason: String,
    sample_relations: Vec<crypto_lib::cryptanalysis::koblitz_sixsum_extract::Relation>,
}

#[derive(Serialize)]
struct RungRecord {
    n: u32,
    a: u8,
    name: String,
    r: u64,
    cofactor: u64,
    group_order: u64,
    bases: Vec<BaseRecord>,
}

#[derive(Serialize)]
struct BaseRecord {
    columns: usize,
    points: usize,
    table: BuildStats,
    base_wall_s: f64,
}

#[derive(Serialize)]
struct Report {
    tool: &'static str,
    version: u32,
    date_utc: String,
    args: BTreeMap<String, String>,
    note: &'static str,
    rungs: Vec<RungRecord>,
    cells: Vec<Cell>,
    total_wall_s: f64,
}

fn build_base(inst: &BinaryInstance, columns: usize, raw_x_cap: u64) -> Result<FactorBase<FastPoint>, String> {
    let group = BinaryGroup(&inst.fast);
    let ctx = InstanceCtx {
        group: &group,
        generator: inst.generator,
        target: group.identity(),
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: inst.name.clone(),
        field_degree: Some(inst.n),
    };
    let mut params = Params::default();
    params.set("columns", columns.to_string());
    params.set("raw_x_cap", raw_x_cap.to_string());
    let mut ops = GroupOps::default();
    CompactOrbitScanBase { instance: inst }.build(&ctx, &params, &mut ops)
}

fn ratio(num: f64, den: u64) -> Option<f64> {
    (den > 0).then(|| num / den as f64)
}

fn main() {
    let args = parse_args();
    let started = Instant::now();
    let mut cells = Vec::new();
    let mut rungs = Vec::new();
    for &n in &args.ns {
        let Some(inst) = koblitz_instance(args.a, n) else {
            eprintln!("# no registered K_{} rung at n = {n}; skipped", args.a);
            continue;
        };
        eprintln!(
            "# {} r = {} ({:.1} bits), cofactor {}",
            inst.name,
            inst.r,
            (inst.r as f64).log2(),
            inst.cofactor
        );
        let mut rung = RungRecord {
            n,
            a: args.a,
            name: inst.name.clone(),
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            bases: Vec::new(),
        };
        let group = BinaryGroup(&inst.fast);
        for &columns in &args.columns {
            let base_start = Instant::now();
            let fb = match build_base(&inst, columns, args.raw_x_cap) {
                Ok(fb) => fb,
                Err(e) => {
                    eprintln!("# n = {n}, K = {columns}: base failed: {e}");
                    continue;
                }
            };
            let base_wall_s = base_start.elapsed().as_secs_f64();
            let table = match SixSumTable::build(&inst, &fb) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("# n = {n}, K = {columns}: table failed: {e}");
                    continue;
                }
            };
            eprintln!(
                "# n = {n}, K = {columns}: {} points, {} table states ({} abscissa states), build {:.2} s",
                fb.points.len(),
                table.len(),
                table.abscissa_states(),
                table.build.wall_ns as f64 * 1e-9
            );
            rung.bases.push(BaseRecord {
                columns,
                points: fb.points.len(),
                table: table.build.clone(),
                base_wall_s,
            });
            let nf = n as f64;
            let kf = columns as f64;
            let rf = inst.r as f64;
            let models = Models {
                r_over_n2k2: rf / (nf * nf * kf * kf),
                r_over_4n2k2: rf / (4.0 * nf * nf * kf * kf),
                model_45r_over_n4k4: 45.0 * rf / (nf.powi(4) * kf.powi(4)),
                structural_s3_solves_per_relation: rf / (4.0 * table.abscissa_states() as f64),
                structural_s4_solves_per_relation: rf / (8.0 * table.abscissa_states() as f64),
            };
            for &kernel in &args.kernels {
                let cell_start = Instant::now();
                let deadline = cell_start + args.cell_wall;
                let mut ctr = ProbeCounters::default();
                let budget = Budget {
                    early_abort: args.early_abort,
                    max_s3_solves: args.max_s3,
                    max_s4_solves: args.max_s4,
                    deadline: Some(deadline),
                };
                let mut rng = StdRng::seed_from_u64(
                    args.seed ^ ((n as u64) << 40) ^ ((columns as u64) << 20) ^ kernel.arity() as u64,
                );
                let mut relations = 0u64;
                let mut full_support = 0u64;
                let mut support_hist: BTreeMap<usize, u64> = BTreeMap::new();
                let mut scan_ends: BTreeMap<String, u64> = BTreeMap::new();
                let mut samples = Vec::new();
                let mut target_ops = GroupOps::default();
                let stop_reason;
                let mut t = 0u64;
                loop {
                    if relations >= args.relations {
                        stop_reason = "relations".to_string();
                        break;
                    }
                    if t >= args.max_targets {
                        stop_reason = "max_targets".to_string();
                        break;
                    }
                    if ctr.s3_solves >= args.max_s3 {
                        stop_reason = "max_s3".to_string();
                        break;
                    }
                    if ctr.s4_solves >= args.max_s4 {
                        stop_reason = "max_s4".to_string();
                        break;
                    }
                    if Instant::now() >= deadline {
                        stop_reason = "cell_wall".to_string();
                        break;
                    }
                    let scalar = rng.gen_range(1..inst.r);
                    let target = group.mul(&mut target_ops, inst.generator, scalar);
                    let start = rng.gen::<u64>();
                    let (rels, end) =
                        table.extract(&inst, &fb, target, kernel, start, &budget, &mut ctr);
                    *scan_ends.entry(format!("{end:?}")).or_insert(0) += 1;
                    for rel in rels {
                        relations += 1;
                        if rel.support == rel.arity {
                            full_support += 1;
                        }
                        *support_hist.entry(rel.support).or_insert(0) += 1;
                        if samples.len() < args.keep_relations {
                            samples.push(rel);
                        }
                    }
                    t += 1;
                    if end == ScanEnd::Budget && (ctr.s3_solves >= args.max_s3 || ctr.s4_solves >= args.max_s4) {
                        stop_reason = if ctr.s3_solves >= args.max_s3 { "max_s3" } else { "max_s4" }.to_string();
                        break;
                    }
                }
                let wall = cell_start.elapsed().as_secs_f64();
                let solves = (ctr.s3_solves + ctr.s4_solves) as f64;
                let targets = ctr.targets;
                let abs_states = table.abscissa_states() as f64;
                let base_f = fb.points.len() as f64;
                let m = kernel.arity() as i32;
                let m_fact: f64 = (1..=m).map(|v| v as f64).product();
                let rel_per_target = base_f.powi(m) / (m_fact * rf);
                let lambda4 = base_f.powi(4) / (24.0 * rf);
                // (outer probes/relation in the dense regime, solves of
                // one hit-or-miss scan unit, expected relations found by
                // one scan unit, dense solves/relation).
                let (structural_outer, scan_solves, unit_lambda, dense) = match kernel {
                    // One target scan; it has a 4-sum with rate λ₄.
                    Kernel::FourSum => (rf / (4.0 * abs_states), abs_states, lambda4, rf / (4.0 * abs_states)),
                    // One inner 4-sum scan of R − s₁ per outer probe.
                    Kernel::IterateTable => (
                        rf / (4.0 * abs_states * abs_states),
                        abs_states,
                        lambda4,
                        rf / (4.0 * abs_states),
                    ),
                    // One inner 4-sum scan of R − P_e per outer probe.
                    Kernel::FiveSumControl => (
                        rf / (4.0 * abs_states * abs_states),
                        abs_states,
                        lambda4,
                        rf / (4.0 * abs_states),
                    ),
                    // One row (s₁ fixed) tests R − s₁ and R + s₁ at once.
                    Kernel::S4PerProbe => (
                        rf / (8.0 * abs_states),
                        abs_states,
                        2.0 * lambda4,
                        rf / (8.0 * abs_states),
                    ),
                    // The whole pair square, de-clustered: the unit is
                    // the target and its 6-sum count.
                    Kernel::S4Strided => (
                        rf / (8.0 * abs_states),
                        abs_states * abs_states / 2.0,
                        rel_per_target,
                        rf / (8.0 * abs_states),
                    ),
                };
                let structural_targets = 1.0 / (1.0 - (-rel_per_target).exp());
                let sparse = scan_solves / unit_lambda;
                let mu = if kernel == Kernel::S4Strided { 45.0 } else { 6.0 };
                let interp = scan_solves
                    * ((-unit_lambda).exp() / (1.0 - (-unit_lambda).exp())
                        + 1.0 / (mu * unit_lambda + 1.0));
                let cell = Cell {
                    n,
                    a: args.a,
                    columns,
                    kernel: kernel.name(),
                    arity: kernel.arity(),
                    r: inst.r,
                    cofactor: inst.cofactor,
                    base_points: fb.points.len(),
                    table_states: table.len(),
                    abscissa_states: table.abscissa_states(),
                    early_abort: args.early_abort,
                    counters: ctr,
                    relations,
                    full_support_relations: full_support,
                    support_histogram: SupportHistogram(support_hist),
                    s3_solves_per_relation: ratio(ctr.s3_solves as f64, relations),
                    s4_solves_per_relation: ratio(ctr.s4_solves as f64, relations),
                    solves_per_relation: ratio(solves, relations),
                    targets_per_relation: ratio(targets as f64, relations),
                    wall_s_per_relation: ratio(wall, relations),
                    point_adds_per_relation: ratio(ctr.point_adds as f64, relations),
                    lookups_per_relation: ratio(ctr.lookups as f64, relations),
                    outer_probes_per_relation: ratio(ctr.outer_probes as f64, relations),
                    structural_outer_probes_per_relation: structural_outer,
                    expected_relations_per_target: rel_per_target,
                    structural_targets_per_relation: structural_targets,
                    scan_unit_solves: scan_solves,
                    scan_unit_lambda: unit_lambda,
                    structural_dense_solves_per_relation: dense,
                    structural_sparse_solves_per_relation: sparse,
                    scan_unit_multiplicity: mu,
                    structural_solves_per_relation: interp,
                    fit_c_solves: ratio(solves, relations).map(|v| v / models.r_over_n2k2),
                    targets_over_model_45: ratio(targets as f64, relations)
                        .map(|v| v / models.model_45r_over_n4k4),
                    outer_probes_over_model_45: ratio(ctr.outer_probes as f64, relations)
                        .map(|v| v / models.model_45r_over_n4k4),
                    models,
                    scan_ends,
                    stop_reason,
                    sample_relations: samples,
                };
                eprintln!(
                    "# n={n} K={columns} {:<22} rel={:>3} targets={:>5} S3={:>11} S4={:>10} S3/rel={:>12} S4/rel={:>10} outer/rel={:>10} tgt/rel={:>8} wall/rel={:>9} ({:.1}s, {})",
                    kernel.name(),
                    relations,
                    targets,
                    ctr.s3_solves,
                    ctr.s4_solves,
                    cell.s3_solves_per_relation.map(|v| format!("{v:.3e}")).unwrap_or("-".into()),
                    cell.s4_solves_per_relation.map(|v| format!("{v:.3e}")).unwrap_or("-".into()),
                    cell.outer_probes_per_relation.map(|v| format!("{v:.3e}")).unwrap_or("-".into()),
                    cell.targets_per_relation.map(|v| format!("{v:.3}")).unwrap_or("-".into()),
                    cell.wall_s_per_relation.map(|v| format!("{v:.3e}")).unwrap_or("-".into()),
                    wall,
                    cell.stop_reason
                );
                cells.push(cell);
            }
        }
        rungs.push(rung);
    }
    let mut argmap = BTreeMap::new();
    argmap.insert("n".into(), format!("{:?}", args.ns));
    argmap.insert("a".into(), args.a.to_string());
    argmap.insert("columns".into(), format!("{:?}", args.columns));
    argmap.insert(
        "kernels".into(),
        args.kernels.iter().map(|k| k.name()).collect::<Vec<_>>().join(","),
    );
    argmap.insert("relations".into(), args.relations.to_string());
    argmap.insert("seed".into(), args.seed.to_string());
    argmap.insert("max_targets".into(), args.max_targets.to_string());
    argmap.insert("max_s3".into(), args.max_s3.to_string());
    argmap.insert("max_s4".into(), args.max_s4.to_string());
    argmap.insert("cell_wall_s".into(), args.cell_wall.as_secs_f64().to_string());
    argmap.insert("early_abort".into(), args.early_abort.to_string());
    argmap.insert("raw_x_cap".into(), args.raw_x_cap.to_string());
    let report = Report {
        tool: "koblitz_sixsum_extract_probe",
        version: 3,
        date_utc: chrono_free_utc_now(),
        args: argmap,
        note: "Public synthetic [a]G targets on registered Koblitz rungs; every relation re-added as a group equation (fast and general arithmetic) before it counts. Units: S3/S4 root solves, targets, wall. No key-recovery, sub-rho or deployed-curve claim.",
        rungs,
        cells,
        total_wall_s: started.elapsed().as_secs_f64(),
    };
    println!("{}", serde_json::to_string_pretty(&report).expect("json"));
}

/// ISO-8601 UTC timestamp without a chrono dependency.
fn chrono_free_utc_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Civil-from-days (Howard Hinnant's algorithm).
    let days = (secs / 86_400) as i64;
    let rem = secs % 86_400;
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!(
        "{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z",
        rem / 3600,
        (rem % 3600) / 60,
        rem % 60
    )
}
