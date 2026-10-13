//! Lead A2 — cyclotomic-sparse Frobenius-stable factor bases against the
//! Hamming-weight ONB control: census, natural 2-sum / 3-sum yields, and
//! the solving degree / F4 wall / WDSat wall of the Weil-descended
//! systems with the two membership encodings.
//!
//! ```bash
//! CARGO_TARGET_DIR=target cargo run --release --example a2_cyclotomic_base -- \
//!     --out research/ic_techniques_20261012/a2_cyclotomic_and_delta \
//!     --targets 16 --yield-targets 64 --cell-timeout-secs 150
//! ```
//!
//! `WDSAT_BINARY=/path/to/wdsat_solver` runs the WDSat leg on an existing
//! binary; `WDSAT_SRC=/path/to/WDSat/src` compiles one sized for these
//! systems into `/tmp` first (WDSat's capacities are compile-time
//! constants).  Without either the WDSat leg is reported as unavailable.
//!
//! The solver cells run in child processes (`--cell JSON`) so a cell that
//! exceeds its wall budget is killed and reported as such instead of
//! stalling the sweep; per-target records already printed are kept.
//!
//! Public synthetic fixtures only: Koblitz `K_a` at `n ∈ {11, 23}`,
//! targets `[k]G` from a seeded stream, no logarithm computed.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::cyclotomic_factor_base::{
    build_cyclotomic_factor_base, build_matched_onb_weight_factor_base, build_membership_system,
    index_map, measure_system, seeded_targets, BaseCensus, BaseFamily, CensusedBase,
    CyclotomicField, MembershipEncoding, MembershipSystem, SystemMeasurement,
};
use crypto_lib::cryptanalysis::koblitz_factor_base_search::{
    census_from_witnesses, TargetSet, Witness, YieldCensus,
};
use crypto_lib::cryptanalysis::koblitz_groebner::{FieldStructure, SolverEngine};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    cofactor_project_factor_base, KoblitzCurve, PairSumTable,
};
use crypto_lib::cryptanalysis::wdsat_oracle::{
    format_anf, parse_wdsat_model, run_wdsat, suggest_capacity, validate_anf, AnfRow,
    WdsatSolveOptions,
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

// ── Arguments ──────────────────────────────────────────────────────

struct Args {
    out: PathBuf,
    targets: usize,
    yield_targets: usize,
    seed: u64,
    cell_timeout: Duration,
    wdsat_timeout: Duration,
    node_budget: usize,
    wdsat_binary: Option<PathBuf>,
    wdsat_src: Option<PathBuf>,
    cell: Option<String>,
    skip_n23_solver: bool,
    /// Macaulay solving-degree probe: `system degree + extra` on systems
    /// of at most 24 unknowns, the system degree alone above that.
    macaulay_extra: u32,
    engine: String,
    /// Witnesses kept per target in the yield census (coverage stays exact).
    witness_cap: usize,
    /// Which solver cells to run: `all`, `n11`, `n23`, or `none`.
    cells: String,
    /// Output file suffix (`a2_results_<label>.json`).
    label: Option<String>,
}

fn parse_args() -> Args {
    let mut args = Args {
        out: PathBuf::from("research/ic_techniques_20261012/a2_cyclotomic_and_delta"),
        targets: 16,
        yield_targets: 64,
        seed: 20261012,
        cell_timeout: Duration::from_secs(150),
        wdsat_timeout: Duration::from_secs(20),
        node_budget: 1024,
        wdsat_binary: std::env::var_os("WDSAT_BINARY").map(PathBuf::from),
        wdsat_src: std::env::var_os("WDSAT_SRC").map(PathBuf::from),
        cell: None,
        skip_n23_solver: false,
        macaulay_extra: 1,
        engine: "inherited".into(),
        witness_cap: 256,
        cells: "all".into(),
        label: None,
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut value = || it.next().expect("flag value");
        match flag.as_str() {
            "--out" => args.out = PathBuf::from(value()),
            "--targets" => args.targets = value().parse().expect("targets"),
            "--yield-targets" => args.yield_targets = value().parse().expect("yield targets"),
            "--seed" => args.seed = value().parse().expect("seed"),
            "--cell-timeout-secs" => {
                args.cell_timeout = Duration::from_secs(value().parse().expect("secs"))
            }
            "--wdsat-timeout-secs" => {
                args.wdsat_timeout = Duration::from_secs(value().parse().expect("secs"))
            }
            "--node-budget" => args.node_budget = value().parse().expect("budget"),
            "--wdsat" => args.wdsat_binary = Some(PathBuf::from(value())),
            "--wdsat-src" => args.wdsat_src = Some(PathBuf::from(value())),
            "--cell" => args.cell = Some(value()),
            "--skip-n23-solver" => args.skip_n23_solver = true,
            "--macaulay-extra" => args.macaulay_extra = value().parse().expect("extra degrees"),
            "--engine" => args.engine = value(),
            "--witness-cap" => args.witness_cap = value().parse().expect("cap"),
            "--cells" => args.cells = value(),
            "--label" => args.label = Some(value()),
            other => panic!("unknown flag {other}"),
        }
    }
    args
}

// ── Cells ──────────────────────────────────────────────────────────

/// One solver cell: a curve, a family, a summand count and a target list.
#[derive(Clone, Debug, Serialize, Deserialize)]
struct CellSpec {
    a: u8,
    n: u32,
    /// Cyclotomic weight; the ONB arm is matched to this weight's columns.
    w: usize,
    m: usize,
    /// `"cyclotomic"` or `"onb_matched"`.
    arm: String,
    targets: usize,
    seed: u64,
    node_budget: usize,
    wdsat_binary: Option<PathBuf>,
    wdsat_timeout_secs: u64,
    macaulay_extra: u32,
    engine: String,
}

fn engine_of(name: &str) -> SolverEngine {
    match name {
        "matrix" | "matrix_f4" => SolverEngine::MatrixF4 { max_degree: 3 },
        "inherited" | "inherited_f4" => SolverEngine::InheritedF4 { max_degree: 3 },
        "f5" | "matrix_f5" => SolverEngine::MatrixF5 { max_degree: 3 },
        other => panic!("unknown engine {other}"),
    }
}

impl CellSpec {
    fn label(&self) -> String {
        format!("K{}/n{}/w{}/m{}/{}", self.a, self.n, self.w, self.m, self.arm)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct WdsatRecord {
    wall_ms: f64,
    verdict: String,
    /// The returned model satisfies every ANF row (checked here, not
    /// trusted from the solver).
    #[serde(default)]
    model_valid: bool,
    lifted: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct TargetLine {
    measurement: SystemMeasurement,
    wdsat: Option<WdsatRecord>,
}

fn curve(a: u8, n: u32) -> KoblitzCurve {
    KoblitzCurve::new(a, n).unwrap_or_else(|| panic!("K_{a} / F_2^{n} has no usable subgroup"))
}

/// The base of a cell, plus the encoding its system uses.
fn cell_base(kc: &KoblitzCurve, cf: &CyclotomicField, spec: &CellSpec) -> (CensusedBase, MembershipEncoding) {
    let cyc = build_cyclotomic_factor_base(kc, cf, spec.w).expect("cyclotomic base");
    match spec.arm.as_str() {
        "cyclotomic" => (cyc, MembershipEncoding::CyclotomicUnity { w: spec.w }),
        "onb_matched" => {
            let columns = cyc.base.signed_orbits.len().max(1);
            let onb = build_matched_onb_weight_factor_base(kc, cf, columns).expect("matched ONB base");
            let BaseFamily::OnbWeight { w } = onb.census.family else {
                unreachable!()
            };
            (onb, MembershipEncoding::OnbWeight { w })
        }
        other => panic!("unknown arm {other}"),
    }
}

/// Highest Macaulay degree the solving-degree probe builds: the system
/// degree plus `extra` on small systems (≤ 24 unknowns), the system degree
/// alone above that, where the sparse profile already self-limits.
fn d_max_for(n_vars: usize, degree: u32, extra: u32) -> u32 {
    if n_vars <= 24 {
        degree + extra
    } else {
        degree
    }
}

fn system_for(
    cf: &CyclotomicField,
    st: &FieldStructure,
    kc: &KoblitzCurve,
    target: &BinaryPoint,
    encoding: MembershipEncoding,
    m: usize,
) -> Option<MembershipSystem> {
    let BinaryPoint::Affine { x, .. } = target else {
        return None;
    };
    build_membership_system(cf, st, &kc.curve.b, x, encoding, m)
}

fn wdsat_verdict(stdout: &str, n_vars: usize) -> (String, Option<Vec<bool>>) {
    for line in stdout.lines() {
        let value = line.trim();
        if matches!(value, "UNSAT" | "UNSAT on XORGAUSS init") {
            return ("unsat".into(), None);
        }
        if !value.is_empty() {
            if let Some(model) = parse_wdsat_model(value, n_vars) {
                return ("model".into(), Some(model));
            }
        }
    }
    ("unknown".into(), None)
}

/// Child mode: run one cell, one JSON line per target on stdout.
fn run_cell_child(spec: &CellSpec) {
    let kc = curve(spec.a, spec.n);
    let cf = CyclotomicField::for_curve(&kc).expect("cyclotomic field");
    let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
    let (built, encoding) = cell_base(&kc, &cf, spec);
    let index_of = index_map(&built.base);
    let targets = seeded_targets(&kc, spec.targets, spec.seed);
    let stdout = std::io::stdout();
    println!(
        "{}",
        json!({"census": built.census, "encoding": encoding}).to_string()
    );
    for (i, (_, target)) in targets.iter().enumerate() {
        let Some(sys) = system_for(&cf, &st, &kc, target, encoding, spec.m) else {
            println!(
                "{}",
                json!({"target_index": i, "unsupported": "system exceeds the 64-variable engine"})
            );
            continue;
        };
        let d_max = d_max_for(sys.n_vars, sys.degree, spec.macaulay_extra);
        let measurement = measure_system(
            &kc,
            &built.base,
            &index_of,
            &sys,
            target,
            i,
            d_max,
            engine_of(&spec.engine),
            spec.node_budget,
        );
        let wdsat = spec.wdsat_binary.as_ref().map(|binary| {
            let rows: Vec<AnfRow> = sys.equations.iter().map(AnfRow::from_poly).collect();
            let anf = format_anf(sys.n_vars as u32, &rows);
            let branch = sys.branch_vars();
            let options = WdsatSolveOptions {
                binary: binary.clone(),
                work_dir: Some(std::env::temp_dir().join("a2-wdsat-runs")),
                timeout: Duration::from_secs(spec.wdsat_timeout_secs),
                keep_anf: None,
            };
            let started = Instant::now();
            let debug_dir = std::env::var_os("A2_WDSAT_DEBUG").map(|_| {
                let dir = std::env::temp_dir().join("a2-wdsat-debug");
                let _ = std::fs::create_dir_all(&dir);
                dir
            });
            let stem = format!("{}-{i}", spec.label().replace('/', "_"));
            if let Some(dir) = &debug_dir {
                let _ = std::fs::write(dir.join(format!("{stem}.anf")), &anf);
            }
            match run_wdsat(&anf, sys.n_vars as u32, &branch, &options) {
                Ok((out, err, elapsed)) => {
                    if let Some(dir) = &debug_dir {
                        let _ = std::fs::write(dir.join(format!("{stem}.out")), &out);
                        let _ = std::fs::write(dir.join(format!("{stem}.err")), &err);
                    }
                    if err.contains("not enough")
                        || out.contains("FATAL")
                        || err.contains("FATAL")
                        || err.contains("Assertion")
                        || err.contains("assert")
                    {
                        return WdsatRecord {
                            wall_ms: elapsed.as_secs_f64() * 1e3,
                            verdict: "capacity".into(),
                            model_valid: false,
                            lifted: false,
                        };
                    }
                    let (verdict, model) = wdsat_verdict(&out, sys.n_vars);
                    let model_valid = model.as_ref().is_some_and(|model| validate_anf(&anf, model));
                    let lifted = model_valid
                        && model.as_ref().is_some_and(|model| {
                            let root = model
                                .iter()
                                .enumerate()
                                .fold(0u64, |acc, (i, &b)| acc | ((b as u64) << i));
                            let xs: Vec<_> =
                                (0..sys.m).map(|i| sys.summand_x(root, i, kc.n)).collect();
                            lift_ok(&kc, &built, &index_of, &xs, target)
                        });
                    WdsatRecord {
                        wall_ms: elapsed.as_secs_f64() * 1e3,
                        verdict,
                        model_valid,
                        lifted,
                    }
                }
                Err(e) => WdsatRecord {
                    wall_ms: started.elapsed().as_secs_f64() * 1e3,
                    verdict: if e.contains("timed out") {
                        "timeout".into()
                    } else {
                        format!("error: {e}")
                    },
                    model_valid: false,
                    lifted: false,
                },
            }
        });
        let line = TargetLine { measurement, wdsat };
        let mut lock = stdout.lock();
        writeln!(lock, "{}", serde_json::to_string(&line).unwrap()).unwrap();
        lock.flush().unwrap();
    }
}

fn lift_ok(
    kc: &KoblitzCurve,
    built: &CensusedBase,
    index_of: &HashMap<(num_bigint::BigUint, num_bigint::BigUint), usize>,
    xs: &[crypto_lib::binary_ecc::F2mElement],
    target: &BinaryPoint,
) -> bool {
    // Same walk as `lift_candidate`, over the public point set.
    fn walk(
        kc: &KoblitzCurve,
        built: &CensusedBase,
        index_of: &HashMap<(num_bigint::BigUint, num_bigint::BigUint), usize>,
        xs: &[crypto_lib::binary_ecc::F2mElement],
        depth: usize,
        acc: &BinaryPoint,
        target: &BinaryPoint,
    ) -> bool {
        if depth == xs.len() {
            return acc == target;
        }
        for p in crypto_lib::cryptanalysis::koblitz_index_calculus::points_with_x(
            &kc.curve, &xs[depth],
        ) {
            let key = crypto_lib::cryptanalysis::koblitz_index_calculus::point_key(&p);
            let Some(&idx) = index_of.get(&key) else {
                continue;
            };
            let next = kc.add(acc, &built.base.points[idx]);
            if walk(kc, built, index_of, xs, depth + 1, &next, target) {
                return true;
            }
        }
        false
    }
    walk(kc, built, index_of, xs, 0, &BinaryPoint::Infinity, target)
}

/// Parent side: spawn the cell as a child with a wall budget.
#[derive(Clone, Debug, Serialize)]
struct CellResult {
    spec: CellSpec,
    census: Option<Value>,
    encoding: Option<Value>,
    lines: Vec<TargetLine>,
    unsupported: usize,
    timed_out: bool,
    wall_ms: f64,
    summary: Value,
}

fn summarise(lines: &[TargetLine]) -> Value {
    if lines.is_empty() {
        return json!({"targets": 0});
    }
    let n = lines.len() as f64;
    let mean = |f: fn(&TargetLine) -> f64| lines.iter().map(f).sum::<f64>() / n;
    let mut f4: Vec<f64> = lines.iter().map(|l| l.measurement.f4_ms).collect();
    f4.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let degrees: Vec<Option<u32>> = lines
        .iter()
        .map(|l| l.measurement.macaulay_solving_degree)
        .collect();
    let resolved: Vec<u32> = degrees.iter().flatten().copied().collect();
    let wd: Vec<&WdsatRecord> = lines.iter().filter_map(|l| l.wdsat.as_ref()).collect();
    json!({
        "targets": lines.len(),
        "n_vars": lines[0].measurement.n_vars,
        "equations": lines[0].measurement.equations,
        "system_degree": lines[0].measurement.system_degree,
        "macaulay_solving_degree_resolved": resolved.len(),
        "macaulay_solving_degree_min": resolved.iter().min(),
        "macaulay_solving_degree_max": resolved.iter().max(),
        "macaulay_degree_reached_max": lines.iter().map(|l| l.measurement.macaulay_degree_reached).max(),
        "macaulay_mean_ms": mean(|l| l.measurement.macaulay_ms),
        "f4_mean_ms": mean(|l| l.measurement.f4_ms),
        "f4_median_ms": f4[f4.len() / 2],
        "f4_mean_word_xors": mean(|l| l.measurement.f4_word_xors as f64),
        "f4_max_degree_built": lines.iter().map(|l| l.measurement.f4_max_degree_built).max(),
        "f4_mean_splits": mean(|l| l.measurement.f4_splits as f64),
        "f4_exhausted": lines.iter().filter(|l| l.measurement.f4_exhausted).count(),
        "boolean_roots_total": lines.iter().map(|l| l.measurement.boolean_roots).sum::<usize>(),
        "lifted_decompositions": lines.iter().filter(|l| l.measurement.lifted_decomposition).count(),
        "wdsat_runs": wd.len(),
        "wdsat_mean_ms": if wd.is_empty() { Value::Null } else { json!(wd.iter().map(|r| r.wall_ms).sum::<f64>() / wd.len() as f64) },
        "wdsat_verdicts": wd.iter().fold(HashMap::<String, usize>::new(), |mut acc, r| { *acc.entry(r.verdict.clone()).or_default() += 1; acc }),
        "wdsat_models_valid": wd.iter().filter(|r| r.model_valid).count(),
        "wdsat_lifted": wd.iter().filter(|r| r.lifted).count(),
    })
}

fn run_cell_parent(spec: &CellSpec, timeout: Duration) -> CellResult {
    let started = Instant::now();
    let exe = std::env::current_exe().expect("current exe");
    let mut child = Command::new(exe)
        .arg("--cell")
        .arg(serde_json::to_string(spec).unwrap())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .expect("spawn cell child");
    let stdout = child.stdout.take().expect("child stdout");
    let (tx, rx) = std::sync::mpsc::channel::<String>();
    std::thread::spawn(move || {
        for line in BufReader::new(stdout).lines().map_while(Result::ok) {
            if tx.send(line).is_err() {
                break;
            }
        }
    });
    let deadline = started + timeout;
    let mut census = None;
    let mut encoding = None;
    let mut lines = Vec::new();
    let mut unsupported = 0;
    let mut timed_out = false;
    loop {
        match rx.recv_timeout(Duration::from_millis(200)) {
            Ok(line) => {
                if let Ok(v) = serde_json::from_str::<Value>(&line) {
                    if v.get("census").is_some() {
                        census = v.get("census").cloned();
                        encoding = v.get("encoding").cloned();
                    } else if v.get("unsupported").is_some() {
                        unsupported += 1;
                    } else if let Ok(t) = serde_json::from_value::<TargetLine>(v) {
                        eprintln!(
                            "    [{}] target {} F4 {:.0} ms deg {:?}/{} roots {} lifted {} wdsat {:?}",
                            spec.label(),
                            t.measurement.target_index,
                            t.measurement.f4_ms,
                            t.measurement.macaulay_solving_degree,
                            t.measurement.macaulay_degree_reached,
                            t.measurement.boolean_roots,
                            t.measurement.lifted_decomposition,
                            t.wdsat.as_ref().map(|w| (&w.verdict, w.wall_ms as u64))
                        );
                        lines.push(t);
                    }
                }
            }
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {}
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => break,
        }
        if let Ok(Some(_)) = child.try_wait() {
            // Drain whatever is left.
            while let Ok(line) = rx.try_recv() {
                if let Ok(t) = serde_json::from_str::<TargetLine>(&line) {
                    lines.push(t);
                }
            }
            break;
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            let _ = child.wait();
            timed_out = true;
            break;
        }
    }
    let _ = child.wait();
    let summary = summarise(&lines);
    CellResult {
        spec: spec.clone(),
        census,
        encoding,
        lines,
        unsupported,
        timed_out,
        wall_ms: started.elapsed().as_secs_f64() * 1e3,
        summary,
    }
}

// ── WDSat build ────────────────────────────────────────────────────

/// What one system needs of WDSat's compile-time capacities.
#[derive(Clone, Copy, Debug, Default, Serialize)]
struct WdsatNeed {
    n_vars: u32,
    /// Distinct monomials of degree ≥ 2 (each becomes a linearised id).
    nonlinear_monomials: u32,
    max_degree: u32,
    equations: u32,
    /// Terms in the longest equation.
    max_terms: u32,
    /// Terms over all equations (the parse buffer).
    total_terms: u32,
}

impl WdsatNeed {
    fn of(sys: &MembershipSystem) -> Self {
        let mut distinct = std::collections::HashSet::new();
        let mut max_degree = 1u32;
        let mut max_terms = 0u32;
        let mut total_terms = 0u32;
        for eq in &sys.equations {
            let mut terms = 0u32;
            for mono in &eq.terms {
                terms += 1;
                let d = mono.mask.count_ones();
                if d >= 2 {
                    distinct.insert(mono.mask);
                }
                max_degree = max_degree.max(d.max(1));
            }
            max_terms = max_terms.max(terms);
            total_terms += terms;
        }
        Self {
            n_vars: sys.n_vars as u32,
            nonlinear_monomials: distinct.len() as u32,
            max_degree,
            equations: sys.equations.len() as u32,
            max_terms,
            total_terms,
        }
    }

    fn max(&self, other: &Self) -> Self {
        Self {
            n_vars: self.n_vars.max(other.n_vars),
            nonlinear_monomials: self.nonlinear_monomials.max(other.nonlinear_monomials),
            max_degree: self.max_degree.max(other.max_degree),
            equations: self.equations.max(other.equations),
            max_terms: self.max_terms.max(other.max_terms),
            total_terms: self.total_terms.max(other.total_terms),
        }
    }

    // The `config.h` constants, with a quarter of slack on the counts
    // that vary with the target's constants.
    fn max_anf_id(&self) -> u64 {
        (self.n_vars as u64 + 1).max(8)
    }
    fn max_degree_plus_one(&self) -> u64 {
        (self.max_degree as u64 + 1).max(3)
    }
    fn max_id(&self) -> u64 {
        (self.n_vars as u64 + self.nonlinear_monomials as u64 * 5 / 4 + 64).max(64)
    }
    fn max_eq(&self) -> u64 {
        (self.nonlinear_monomials as u64 * 5 / 4 * (self.max_degree as u64 + 1) + 64).max(64)
    }
    fn max_eq_size(&self) -> u64 {
        (self.max_degree as u64 + 2).max(4)
    }
    fn max_xeq(&self) -> u64 {
        (self.equations as u64 + 64).max(64)
    }
    fn max_xeq_size(&self) -> u64 {
        (self.max_terms as u64 * 5 / 4 + 64).max(64)
    }
    fn max_buffer(&self) -> u64 {
        (self.total_terms as u64 * (self.max_degree as u64 + 2) * 5 / 4 + 5000).max(5000)
    }

    /// Bytes of the static arrays WDSat declares with these constants
    /// (`xorgauss.c`, `xorset.c`, `dimacs.c`, `cnf.c`), the quantity the
    /// arm64 linker and the host's memory see.
    fn static_bytes(&self) -> u64 {
        let id_size = self.max_id() + 1;
        let sz_gauss = id_size / 64 + 1;
        let anf = self.max_anf_id();
        let eqv_history = anf * id_size * sz_gauss * 8;
        let eqv = id_size * sz_gauss * 8 * 2;
        let monomials_to_column = anf * id_size * (self.max_degree_plus_one() - 1) * 8;
        let xorset_history = 2 * id_size * self.max_xeq_size() * 8;
        let xor_equations = 2 * self.max_xeq() * self.max_xeq_size() * 8;
        let boolean_equations = self.max_eq() * self.max_eq_size() * 8;
        let buffers = 2 * self.max_buffer() * 8;
        let small = anf * id_size * 8 * 3;
        eqv_history
            + eqv
            + monomials_to_column
            + xorset_history
            + xor_equations
            + boolean_equations
            + buffers
            + small
    }

    fn config_h(&self) -> String {
        format!(
            "#define __XG_ENHANCED__\n#ifdef __XG_ENHANCED__\n#define __MAX_ANF_ID__ {}\n#define __MAX_DEGREE__ {}\n#endif\n#define __MAX_ID__ {}\n#define __MAX_BUFFER_SIZE__ {}\n#define __MAX_EQ__ {}\n#define __MAX_EQ_SIZE__ {}\n#define __MAX_XEQ__ {}\n#define __MAX_XEQ_SIZE__ {}\n",
            self.max_anf_id(),
            self.max_degree_plus_one(),
            self.max_id(),
            self.max_buffer(),
            self.max_eq(),
            self.max_eq_size(),
            self.max_xeq(),
            self.max_xeq_size()
        )
    }
}

/// Compile a WDSat sized for the largest system of the cells, if a source
/// tree was given.
fn prepare_wdsat(args: &Args, cells: &[CellSpec]) -> (Option<PathBuf>, Value) {
    if let Some(binary) = &args.wdsat_binary {
        return (
            Some(binary.clone()),
            json!({"status": "provided", "binary": binary}),
        );
    }
    let Some(src) = &args.wdsat_src else {
        return (
            None,
            json!({"status": "unavailable", "reason": "neither WDSAT_BINARY nor WDSAT_SRC set and no wdsat_solver on PATH"}),
        );
    };
    if !src.join("wdsat.c").is_file() {
        return (
            None,
            json!({"status": "unavailable", "reason": format!("{} has no wdsat.c", src.display())}),
        );
    }
    // Size the capacities from the first two targets of every cell.  The
    // library's `write_config_for_anf` sizes `__MAX_ID__` and
    // `__MAX_XEQ_SIZE__` from the *total* XOR-atom count, which on these
    // dense Weil-descent systems makes WDSat's `xorset_history` and
    // `xorgauss_equivalency_history` arrays (quadratic in `__MAX_ID__`)
    // run to gigabytes and the arm64 link fail (ADRP out of range).  The
    // capacities here come from the monomial structure itself: distinct
    // monomials, the longest equation, the equation count.
    let mut sized = Vec::new();
    let mut per_cell: Vec<(String, WdsatNeed)> = Vec::new();
    for spec in cells {
        let kc = curve(spec.a, spec.n);
        let cf = CyclotomicField::for_curve(&kc).expect("field");
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let (_, encoding) = cell_base(&kc, &cf, spec);
        let targets = seeded_targets(&kc, 2, spec.seed);
        let mut need = WdsatNeed::default();
        let mut atoms = 0u32;
        for (_, target) in &targets {
            let Some(sys) = system_for(&cf, &st, &kc, target, encoding, spec.m) else {
                continue;
            };
            need = need.max(&WdsatNeed::of(&sys));
            let rows: Vec<AnfRow> = sys.equations.iter().map(AnfRow::from_poly).collect();
            atoms = atoms.max(suggest_capacity(&format_anf(sys.n_vars as u32, &rows)).xor_atom_count);
        }
        if need.n_vars == 0 {
            continue;
        }
        sized.push(json!({"cell": spec.label(), "xor_atoms": atoms, "need": need, "static_bytes_alone": need.static_bytes()}));
        per_cell.push((spec.label(), need));
    }
    if per_cell.is_empty() {
        return (None, json!({"status": "unavailable", "reason": "no cell produced a system"}));
    }
    // Union the cells that fit under the static-memory cap, smallest first,
    // so one oversized cell costs only itself.
    const STATIC_CAP: u64 = 2_500_000_000;
    per_cell.sort_by_key(|(_, need)| need.static_bytes());
    let mut union = WdsatNeed::default();
    let mut covered = Vec::new();
    let mut over_capacity = Vec::new();
    for (label, need) in &per_cell {
        let candidate = union.max(need);
        if candidate.static_bytes() <= STATIC_CAP {
            union = candidate;
            covered.push(label.clone());
        } else {
            over_capacity.push(label.clone());
        }
    }
    let build_dir = std::env::temp_dir().join(format!("a2-wdsat-build-{}", std::process::id()));
    let src_copy = build_dir.join("src");
    std::fs::create_dir_all(&src_copy).expect("build dir");
    for entry in std::fs::read_dir(src).expect("read WDSat src") {
        let entry = entry.expect("entry");
        let name = entry.file_name();
        let name_str = name.to_string_lossy();
        if name_str.ends_with(".c") || name_str.ends_with(".h") {
            std::fs::copy(entry.path(), src_copy.join(&name)).expect("copy source");
        }
    }
    std::fs::write(src_copy.join("config.h"), union.config_h()).expect("config.h");
    let cap = &union;
    let sources: Vec<PathBuf> = std::fs::read_dir(&src_copy)
        .expect("read copy")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|e| e == "c"))
        .collect();
    let binary = build_dir.join("wdsat_solver");
    let output = Command::new("cc")
        .arg("-O3")
        .arg("-w")
        .args(&sources)
        .arg("-o")
        .arg(&binary)
        .arg("-lm")
        .output();
    match output {
        Ok(out) if out.status.success() && binary.is_file() => (
            Some(binary.clone()),
            json!({"status": "built", "binary": binary, "source": src, "config_h": cap.config_h(), "static_bytes": cap.static_bytes(), "covered_cells": covered, "over_capacity_cells": over_capacity, "cells": sized}),
        ),
        Ok(out) => (
            None,
            json!({"status": "build_failed", "stderr": String::from_utf8_lossy(&out.stderr).chars().take(2000).collect::<String>(), "cells": sized}),
        ),
        Err(e) => (None, json!({"status": "build_failed", "error": e.to_string()})),
    }
}

// ── Census and yields (parent) ─────────────────────────────────────

#[derive(Clone, Debug, Serialize)]
struct YieldRow {
    family: BaseFamily,
    /// `false`: the raw rational points above the set; `true`: their
    /// non-identity cofactor images `[h]P` (the subgroup-usable base the
    /// relation collector would actually sum over).
    projected: bool,
    columns: usize,
    points: usize,
    m: usize,
    census: YieldCensus,
    /// `coverage / columns`: the falsifier's yield per column.
    coverage_per_column: f64,
    table_tier: String,
    table_ms: f64,
}

/// [`witness_list`] with at most `cap` witnesses kept per target.
/// Coverage (does the target decompose at all) is exact either way; the
/// witness counts are censored at `cap`, which on the oversaturated
/// `w = 3` bases (≈ 10⁵ 3-sum witnesses per target) is what keeps the
/// census inside the run budget.
fn capped_witness_list(
    kc: &KoblitzCurve,
    fb: &crypto_lib::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    table: &PairSumTable,
    targets: &TargetSet,
    m: usize,
    cap: usize,
) -> Vec<Witness> {
    let mut out = Vec::new();
    for (t, target) in targets.points.iter().enumerate() {
        let mut kept = 0usize;
        table.witnesses(kc, fb, target, m, &mut |idxs| {
            let mut orbits: Vec<u32> = idxs
                .iter()
                .map(|&i| fb.signed_orbit_of[i].0 as u32)
                .collect();
            orbits.sort_unstable();
            orbits.dedup();
            out.push(Witness {
                target: t as u32,
                orbits,
            });
            kept += 1;
            kept < cap
        });
    }
    out
}

fn yields_of(
    kc: &KoblitzCurve,
    family: BaseFamily,
    base: &crypto_lib::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    projected: bool,
    targets: &TargetSet,
    witness_cap: usize,
) -> Vec<YieldRow> {
    let mut rows = Vec::new();
    if base.points.len() < 4 {
        return rows;
    }
    let started = Instant::now();
    let Some(table) = PairSumTable::build(kc, base) else {
        return rows;
    };
    let table_ms = started.elapsed().as_secs_f64() * 1e3;
    for m in [2usize, 3] {
        let started = Instant::now();
        let witnesses = capped_witness_list(kc, base, &table, targets, m, witness_cap);
        let census = census_from_witnesses(
            &witnesses,
            targets,
            m,
            base.signed_orbits.len(),
            base.signed_orbits.len(),
            0,
            started.elapsed().as_secs_f64() * 1e3,
        );
        let columns = base.signed_orbits.len().max(1);
        rows.push(YieldRow {
            family,
            projected,
            columns,
            points: base.points.len(),
            m,
            coverage_per_column: census.coverage / columns as f64,
            census,
            table_tier: table.tier().to_string(),
            table_ms,
        });
    }
    rows
}

/// Yields of a base on its raw points and on its cofactor projection.
fn yields(
    kc: &KoblitzCurve,
    built: &CensusedBase,
    targets: &TargetSet,
    witness_cap: usize,
) -> Vec<YieldRow> {
    let mut rows = yields_of(kc, built.census.family, &built.base, false, targets, witness_cap);
    if let Ok(projected) = cofactor_project_factor_base(kc, &built.base) {
        rows.extend(yields_of(kc, built.census.family, &projected, true, targets, witness_cap));
    }
    rows
}

fn main() {
    let args = parse_args();
    if let Some(cell) = &args.cell {
        let spec: CellSpec = serde_json::from_str(cell).expect("cell spec");
        run_cell_child(&spec);
        return;
    }
    let sweep_started = Instant::now();
    std::fs::create_dir_all(&args.out).expect("output dir");

    // 1. Census on every rung with the structure, both models.
    eprintln!("== census ==");
    let mut census_rows: Vec<BaseCensus> = Vec::new();
    let mut matched_pairs: Vec<Value> = Vec::new();
    let rungs: [(u8, u32, &[usize]); 3] = [(1, 11, &[1, 2, 3]), (0, 23, &[1, 2, 3]), (1, 23, &[1, 2, 3])];
    let mut bases: HashMap<(u8, u32, usize), (CensusedBase, Option<CensusedBase>)> = HashMap::new();
    for (a, n, ws) in rungs {
        let kc = curve(a, n);
        let cf = CyclotomicField::for_curve(&kc).expect("cyclotomic field");
        eprintln!(
            "K_{a}/F_2^{n}: p = {}, γ = {:#x}, #E = {}, r = {}, h = {}",
            cf.p,
            cf.gamma.raw_bits()[0],
            kc.group_order,
            kc.subgroup_order,
            kc.cofactor
        );
        for &w in ws {
            let cyc = match build_cyclotomic_factor_base(&kc, &cf, w) {
                Ok(b) => b,
                Err(e) => {
                    eprintln!("  w = {w}: {e}");
                    continue;
                }
            };
            eprintln!(
                "  cyclotomic w={w}: |C_w|={} (C(p,w)={}) stable={} points={} in_subgroup={} orbits={} signed={} projected_cols={} [{:.0} ms]",
                cyc.census.abscissae, cyc.census.predicted_abscissae, cyc.census.frobenius_stable,
                cyc.census.points, cyc.census.points_in_subgroup, cyc.census.orbits,
                cyc.census.signed_orbits, cyc.census.projected_columns, cyc.census.build_ms
            );
            census_rows.push(cyc.census.clone());
            let columns = cyc.base.signed_orbits.len();
            let onb = if columns > 0 {
                match build_matched_onb_weight_factor_base(&kc, &cf, columns) {
                    Ok(onb) => {
                        eprintln!(
                            "  matched ONB {:?}: abscissae={} points={} in_subgroup={} signed={} projected_cols={} retained={:?}",
                            onb.census.family, onb.census.abscissae, onb.census.points,
                            onb.census.points_in_subgroup, onb.census.signed_orbits,
                            onb.census.projected_columns, onb.census.retained_orbits
                        );
                        census_rows.push(onb.census.clone());
                        matched_pairs.push(json!({"a": a, "n": n, "w": w, "cyclotomic_columns": columns, "onb": onb.census}));
                        Some(onb)
                    }
                    Err(e) => {
                        eprintln!("  matched ONB: {e}");
                        None
                    }
                }
            } else {
                None
            };
            bases.insert((a, n, w), (cyc, onb));
        }
    }

    // 2. Natural 2-sum / 3-sum yields at n = 23 (and the n = 11 toy).
    eprintln!("== yields ==");
    let mut yield_rows: Vec<Value> = Vec::new();
    let yield_rungs: &[(u8, u32)] = if args.yield_targets == 0 {
        eprintln!("  (skipped: --yield-targets 0)");
        &[]
    } else {
        &[(1u8, 11u32), (0, 23), (1, 23)]
    };
    for &(a, n) in yield_rungs {
        let kc = curve(a, n);
        let targets = TargetSet::new(&kc, args.yield_targets, 0, args.seed);
        for w in [2usize, 3] {
            let Some((cyc, onb)) = bases.get(&(a, n, w)) else {
                continue;
            };
            for built in std::iter::once(cyc).chain(onb.iter()) {
                for row in yields(&kc, built, &targets, args.witness_cap) {
                    eprintln!(
                        "  K_{a}/n{n} w={w} {:?}{} m={}: coverage {:.4} ({}/{}), mean witnesses {:.3}, per column {:.3e}, columns {}, points {} [{} table]",
                        row.family, if row.projected { " [projected]" } else { "" }, row.m, row.census.coverage, row.census.covered, row.census.targets,
                        row.census.mean_witnesses, row.coverage_per_column, row.columns, row.points, row.table_tier
                    );
                    let mut v = serde_json::to_value(&row).unwrap();
                    v["a"] = json!(a);
                    v["n"] = json!(n);
                    v["cyclotomic_w"] = json!(w);
                    // Drop the per-orbit witness vector from the report.
                    if let Some(c) = v.get_mut("census") {
                        c.as_object_mut().unwrap().remove("orbit_witnesses");
                    }
                    yield_rows.push(v);
                }
            }
        }
    }
    drop(bases);

    // 3. Solver cells.
    let mut cells: Vec<CellSpec> = Vec::new();
    let mk = |a: u8, n: u32, w: usize, m: usize, arm: &str| CellSpec {
        a,
        n,
        w,
        m,
        arm: arm.into(),
        targets: args.targets,
        seed: args.seed,
        node_budget: args.node_budget,
        wdsat_binary: None,
        wdsat_timeout_secs: args.wdsat_timeout.as_secs(),
        macaulay_extra: args.macaulay_extra,
        engine: args.engine.clone(),
    };
    let want_n11 = matches!(args.cells.as_str(), "all" | "n11");
    let want_n23 = matches!(args.cells.as_str(), "all" | "n23") && !args.skip_n23_solver;
    if want_n11 {
        for arm in ["cyclotomic", "onb_matched"] {
            cells.push(mk(1, 11, 1, 2, arm));
        }
        for arm in ["cyclotomic", "onb_matched"] {
            cells.push(mk(1, 11, 2, 2, arm));
        }
        for arm in ["cyclotomic", "onb_matched"] {
            cells.push(mk(1, 11, 1, 3, arm));
        }
    }
    if want_n23 {
        for arm in ["cyclotomic", "onb_matched"] {
            cells.push(mk(0, 23, 1, 2, arm));
        }
        for arm in ["cyclotomic", "onb_matched"] {
            cells.push(mk(1, 23, 1, 2, arm));
        }
        // w = 2 at n = 23 needs 92 unknowns: reported as unsupported by the child.
        cells.push(mk(0, 23, 2, 2, "cyclotomic"));
    }
    let (wdsat_binary, wdsat_status) = if cells.is_empty() {
        (None, json!({"status": "unavailable", "reason": "no solver cells requested"}))
    } else {
        prepare_wdsat(&args, &cells)
    };
    eprintln!("== wdsat: {} ==", wdsat_status["status"]);
    for cell in &mut cells {
        cell.wdsat_binary = wdsat_binary.clone();
    }
    eprintln!("== solver cells ({} × {} targets, {} s each) ==", cells.len(), args.targets, args.cell_timeout.as_secs());
    let mut results = Vec::new();
    for spec in &cells {
        eprintln!("-- {}", spec.label());
        let result = run_cell_parent(spec, args.cell_timeout);
        eprintln!(
            "   {} targets, timed_out={}, unsupported={}, {}",
            result.lines.len(),
            result.timed_out,
            result.unsupported,
            result.summary
        );
        results.push(result);
    }

    let report = json!({
        "lead": "A2 cyclotomic-sparse Frobenius-stable factor bases",
        "generated_utc": chrono_like_now(),
        "args": {
            "targets": args.targets, "yield_targets": args.yield_targets, "seed": args.seed,
            "cell_timeout_secs": args.cell_timeout.as_secs(), "wdsat_timeout_secs": args.wdsat_timeout.as_secs(),
            "node_budget": args.node_budget, "solver_engine": args.engine,
            "solver_engine_note": "Boolean F4 ladder from the system degree (max_degree 3 is the floor; a degree-d system builds degree d), splitting on the engine's rule when the algebra stalls",
            "macaulay_extra": args.macaulay_extra,
            "witness_cap_per_target": args.witness_cap,
        },
        "host": {"os": std::env::consts::OS, "arch": std::env::consts::ARCH},
        "wdsat": wdsat_status,
        "census": census_rows,
        "matched_pairs": matched_pairs,
        "yields": yield_rows,
        "solver_cells": results,
        "wall_ms": sweep_started.elapsed().as_secs_f64() * 1e3,
    });
    let suffix = args.label.as_ref().map(|l| format!("_{l}")).unwrap_or_default();
    let path = args.out.join(format!("a2_results{suffix}.json"));
    std::fs::write(&path, serde_json::to_string_pretty(&report).unwrap()).expect("write results");
    eprintln!("wrote {} after {:.1} s", path.display(), sweep_started.elapsed().as_secs_f64());
    write_summary_table(&args.out.join(format!("a2_summary{suffix}.md")), &report);
}

fn chrono_like_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    format!("unix:{secs}")
}

fn write_summary_table(path: &Path, report: &Value) {
    let mut s = String::new();
    s.push_str("# A2 summary (generated)\n\n## Census\n\n| a | n | family | predicted | abscissae | stable | with points | points | in subgroup | orbits | signed | projected cols | retained |\n|--:|--:|:--|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|\n");
    for c in report["census"].as_array().unwrap() {
        s.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            c["a"], c["n"], c["family"]["family"].as_str().unwrap_or("?").to_string() + &format!("(w={})", c["family"]["w"]),
            c["predicted_abscissae"], c["abscissae"], c["frobenius_stable"], c["abscissae_with_points"], c["points"],
            c["points_in_subgroup"], c["orbits"], c["signed_orbits"], c["projected_columns"], c["retained_orbits"]
        ));
    }
    s.push_str("\n## Yields\n\n| a | n | cyc w | family | points | m | columns | points | coverage | covered/targets | mean witnesses | coverage per column |\n|--:|--:|--:|:--|:--|--:|--:|--:|--:|:--|--:|--:|\n");
    for y in report["yields"].as_array().unwrap() {
        s.push_str(&format!(
            "| {} | {} | {} | {}(w={}) | {} | {} | {} | {} | {:.4} | {}/{} | {:.3} | {:.3e} |\n",
            y["a"], y["n"], y["cyclotomic_w"], y["family"]["family"].as_str().unwrap_or("?"), y["family"]["w"],
            if y["projected"].as_bool().unwrap_or(false) { "projected" } else { "raw" }, y["m"], y["columns"], y["points"],
            y["census"]["coverage"].as_f64().unwrap_or(0.0), y["census"]["covered"], y["census"]["targets"],
            y["census"]["mean_witnesses"].as_f64().unwrap_or(0.0), y["coverage_per_column"].as_f64().unwrap_or(0.0)
        ));
    }
    s.push_str("\n## Solver cells\n\n| cell | encoding | vars | eqs | sys deg | Macaulay solving deg (resolved/targets) | F4 mean ms | F4 median ms | F4 mean word XORs | F4 max deg | mean splits | lifted | WDSat runs | WDSat mean ms | WDSat verdicts | timed out |\n|:--|:--|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|\n");
    for c in report["solver_cells"].as_array().unwrap() {
        let sm = &c["summary"];
        let enc = c["encoding"].as_object().map(|e| format!("{}(w={})", e.get("encoding").and_then(|v| v.as_str()).unwrap_or("?"), e.get("w").cloned().unwrap_or(Value::Null))).unwrap_or_else(|| "?".into());
        s.push_str(&format!(
            "| {} | {} | {} | {} | {} | {}..{} ({}/{}) reached {} | {:.1} | {:.1} | {:.3e} | {} | {:.2} | {} | {} | {} | {} | {} |\n",
            c["spec"]["a"].as_u64().map(|a| format!("K{a}/n{}/w{}/m{}/{}", c["spec"]["n"], c["spec"]["w"], c["spec"]["m"], c["spec"]["arm"].as_str().unwrap_or("?"))).unwrap_or_default(),
            enc, sm["n_vars"], sm["equations"], sm["system_degree"],
            sm["macaulay_solving_degree_min"], sm["macaulay_solving_degree_max"], sm["macaulay_solving_degree_resolved"], sm["targets"], sm["macaulay_degree_reached_max"],
            sm["f4_mean_ms"].as_f64().unwrap_or(0.0), sm["f4_median_ms"].as_f64().unwrap_or(0.0), sm["f4_mean_word_xors"].as_f64().unwrap_or(0.0),
            sm["f4_max_degree_built"], sm["f4_mean_splits"].as_f64().unwrap_or(0.0), sm["lifted_decompositions"],
            sm["wdsat_runs"], sm["wdsat_mean_ms"], sm["wdsat_verdicts"], c["timed_out"]
        ));
    }
    std::fs::write(path, s).expect("write summary");
}
