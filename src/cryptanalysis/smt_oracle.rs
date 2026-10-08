//! SMT-LIB 2 point-decomposition oracle: cvc5, Z3, or any SMT-LIB solver
//! on the Weil-descended Semaev systems the native CDCL+XOR path
//! ([`crate::cryptanalysis::koblitz_index_calculus::sat_decompose`]) and
//! the WDSat oracle ([`crate::cryptanalysis::wdsat_oracle`]) already solve.
//!
//! The instance is the same `F_2` polynomial system: `m` summand
//! abscissae confined to an `ℓ`-dimensional factor-base subspace, the
//! summation polynomial descended to `n` Boolean equations, plus the
//! group trace row and, optionally, the degree-`D` Macaulay rows the
//! native path adds.  What changes is only who searches it, so a run of
//! this oracle against the native solver is a measurement of the solver,
//! not of the problem.
//!
//! Two encodings are emitted, because SMT solvers route them to different
//! engines:
//!
//! - [`SmtEncoding::Bool`]: `QF_UF` with Boolean constants, `xor` and
//!   `and`.  cvc5 and Z3 hand this to their propositional core; neither
//!   has native parity reasoning, so a dense `xor` row is Tseitin-split
//!   much as [`crate::cryptanalysis::semaev_sat::XorEncoding::Cnf`] does.
//! - [`SmtEncoding::BitVec1`]: `QF_BV` over `(_ BitVec 1)` with `bvxor`
//!   and `bvand`.  Both solvers bit-blast this through their bit-vector
//!   preprocessors (cvc5 to CaDiCaL by default, Z3 to its own SAT core),
//!   which is where any word-level simplification they own would act.
//!
//! Models are enumerated with blocking assertions, as the native path
//! does, so a final `unsat` is a proven refutation and every accepted
//! model is checked against the original equations and lifted through
//! the group before it counts.  Solver-reported conflict and decision
//! counts are harvested when the solver prints them; they are the
//! solver's own counters, not this repository's counted unit, and are
//! recorded as a practicality note beside wall time.
//!
//! Scope: a stage diagnostic of the decomposition oracle.  Nothing here
//! prices relation collection, linear algebra or the whole method, so no
//! `S`, speed or `vs_rho` figure follows from it.

use crate::binary_ecc::{BinaryPoint, F2mElement};
use crate::cryptanalysis::koblitz_groebner::{matrix_f4_f2, DecompositionSystem, FieldStructure};
use crate::cryptanalysis::koblitz_index_calculus::{
    absolute_trace_bit, lift_candidate, FrobeniusFactorBase, KoblitzCurve,
    SatDecompositionStats,
};
use crate::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use crate::cryptanalysis::wdsat_oracle::AnfRow;
use num_bigint::BigUint;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

/// How the Boolean system is written in SMT-LIB 2.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmtEncoding {
    /// `QF_UF`: Boolean constants, `xor` rows, `and` monomials.
    Bool,
    /// `QF_BV`: one-bit bit-vectors, `bvxor` rows, `bvand` monomials.
    BitVec1,
}

impl SmtEncoding {
    /// The `set-logic` name.
    pub fn logic(self) -> &'static str {
        match self {
            Self::Bool => "QF_UF",
            Self::BitVec1 => "QF_BV",
        }
    }

    /// Parse the names used on command lines and in reports.
    pub fn parse(name: &str) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "bool" | "qf_uf" => Some(Self::Bool),
            "bv" | "bv1" | "bitvec1" | "qf_bv" => Some(Self::BitVec1),
            _ => None,
        }
    }

    /// The name reports use.
    pub fn name(self) -> &'static str {
        match self {
            Self::Bool => "bool",
            Self::BitVec1 => "bv1",
        }
    }
}

/// Variable name of the zero-based solver variable `i`.
fn var_name(i: u32) -> String {
    format!("x{i}")
}

/// One monomial as an SMT term.
fn monomial_term(vars: &[u32], encoding: SmtEncoding) -> String {
    match vars {
        [] => match encoding {
            SmtEncoding::Bool => "true".into(),
            SmtEncoding::BitVec1 => "#b1".into(),
        },
        [v] => var_name(*v),
        _ => {
            let op = match encoding {
                SmtEncoding::Bool => "and",
                SmtEncoding::BitVec1 => "bvand",
            };
            let mut s = format!("({op}");
            for v in vars {
                s.push(' ');
                s.push_str(&var_name(*v));
            }
            s.push(')');
            s
        }
    }
}

/// One ANF row `Σ monomials + constant = 0` as an SMT assertion, or
/// `None` when the row is the trivial `0 = 0`.
fn row_assertion(row: &AnfRow, encoding: SmtEncoding) -> Option<String> {
    let terms: Vec<String> = row
        .monomials
        .iter()
        .map(|m| monomial_term(m, encoding))
        .collect();
    // Σ monomials = constant (over F_2).
    let (xor_op, one, zero) = match encoding {
        SmtEncoding::Bool => ("xor", "true", "false"),
        SmtEncoding::BitVec1 => ("bvxor", "#b1", "#b0"),
    };
    let rhs = if row.constant { one } else { zero };
    let lhs = match terms.len() {
        0 => return row.constant.then(|| "(assert false)".to_string()),
        1 => terms[0].clone(),
        _ => {
            let mut s = format!("({xor_op}");
            for t in &terms {
                s.push(' ');
                s.push_str(t);
            }
            s.push(')');
            s
        }
    };
    Some(match encoding {
        SmtEncoding::Bool if row.constant => format!("(assert {lhs})"),
        SmtEncoding::Bool => format!("(assert (not {lhs}))"),
        SmtEncoding::BitVec1 => format!("(assert (= {lhs} {rhs}))"),
    })
}

/// A blocking assertion excluding one full assignment.
fn blocking_assertion(model: &[bool], encoding: SmtEncoding) -> String {
    let mut s = String::from("(assert (not (and");
    for (i, bit) in model.iter().enumerate() {
        let v = var_name(i as u32);
        match (encoding, *bit) {
            (SmtEncoding::Bool, true) => write!(s, " {v}").unwrap(),
            (SmtEncoding::Bool, false) => write!(s, " (not {v})").unwrap(),
            (SmtEncoding::BitVec1, true) => write!(s, " (= {v} #b1)").unwrap(),
            (SmtEncoding::BitVec1, false) => write!(s, " (= {v} #b0)").unwrap(),
        }
    }
    s.push_str(")))");
    s
}

/// Write the system as a complete SMT-LIB 2 script: declarations, one
/// assertion per row, any blocking assertions, `check-sat`, and a
/// `get-value` over every variable so the model is read back in order.
pub fn format_smtlib(
    n_vars: u32,
    rows: &[AnfRow],
    blocked: &[Vec<bool>],
    encoding: SmtEncoding,
) -> String {
    let mut s = String::new();
    s.push_str("(set-option :produce-models true)\n");
    writeln!(s, "(set-logic {})", encoding.logic()).unwrap();
    let sort = match encoding {
        SmtEncoding::Bool => "Bool",
        SmtEncoding::BitVec1 => "(_ BitVec 1)",
    };
    for i in 0..n_vars {
        writeln!(s, "(declare-const {} {sort})", var_name(i)).unwrap();
    }
    for row in rows {
        if let Some(a) = row_assertion(row, encoding) {
            s.push_str(&a);
            s.push('\n');
        }
    }
    for model in blocked {
        s.push_str(&blocking_assertion(model, encoding));
        s.push('\n');
    }
    s.push_str("(check-sat)\n(get-value (");
    for i in 0..n_vars {
        if i > 0 {
            s.push(' ');
        }
        s.push_str(&var_name(i));
    }
    s.push_str("))\n(exit)\n");
    s
}

/// Which solver is being driven; fixes the command line and where the
/// solver prints its counters.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SmtSolverKind {
    Cvc5,
    Z3,
    /// Any SMT-LIB 2 solver given its full argument list.
    Generic,
}

impl SmtSolverKind {
    /// Parse the names used on command lines and in reports.
    pub fn parse(name: &str) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "cvc5" => Some(Self::Cvc5),
            "z3" => Some(Self::Z3),
            "generic" | "smtlib" => Some(Self::Generic),
            _ => None,
        }
    }

    /// The name reports use.
    pub fn name(self) -> &'static str {
        match self {
            Self::Cvc5 => "cvc5",
            Self::Z3 => "z3",
            Self::Generic => "generic",
        }
    }

    /// Default arguments placed before the script path.  Statistics are
    /// requested so conflict and decision counts can be harvested.
    pub fn default_args(self) -> Vec<String> {
        match self {
            Self::Cvc5 => ["--lang", "smt2", "--stats-internal"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            Self::Z3 => ["-smt2", "-st"].iter().map(|s| s.to_string()).collect(),
            Self::Generic => Vec::new(),
        }
    }
}

/// How to invoke an external SMT solver.
#[derive(Clone, Debug)]
pub struct SmtSolveOptions {
    /// Which solver, for defaults and counter parsing.
    pub kind: SmtSolverKind,
    /// Path to the solver executable.
    pub binary: PathBuf,
    /// Arguments placed before the script path; `None` uses the kind's defaults.
    pub args: Option<Vec<String>>,
    /// Working directory for the script and the child process.
    pub work_dir: Option<PathBuf>,
    /// Soft wall-clock budget per solver call.
    pub timeout: Duration,
    /// Keep the generated script under this path instead of a temporary file
    /// (each enumeration round appends `.k` for its blocking depth).
    pub keep_script: Option<PathBuf>,
    /// Boolean or one-bit bit-vector encoding.
    pub encoding: SmtEncoding,
    /// Add the Kosters–Yeo trace row, as the native path does by default.
    pub trace_constraint: bool,
    /// Add the degree-`D` Macaulay rows the native path adds for `Some(D)`.
    pub macaulay_degree: Option<u32>,
    /// Models to examine before giving up, as `sat_decompose`'s `max_models`.
    pub max_models: usize,
}

impl SmtSolveOptions {
    /// Defaults matching `sat_decompose(.., max_models = 64, macaulay_degree = Some(2))`.
    pub fn new(kind: SmtSolverKind, binary: impl Into<PathBuf>) -> Self {
        Self {
            kind,
            binary: binary.into(),
            args: None,
            work_dir: None,
            timeout: Duration::from_secs(600),
            keep_script: None,
            encoding: SmtEncoding::Bool,
            trace_constraint: true,
            macaulay_degree: Some(2),
            max_models: 64,
        }
    }
}

/// What one solver invocation returned.
#[derive(Clone, Debug, Default)]
pub struct SmtRun {
    pub stdout: String,
    pub stderr: String,
    pub elapsed: Duration,
    /// Solver-reported conflicts, when printed.
    pub conflicts: Option<u64>,
    /// Solver-reported decisions, when printed.
    pub decisions: Option<u64>,
}

/// The verdict on the first `sat`/`unsat`/`unknown` line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SmtVerdict {
    Sat(Vec<bool>),
    Unsat,
    Unknown,
}

/// Read the `get-value` answer: every `(xK v)` pair with `v` in
/// `true/false/#b1/#b0`.  Returns `None` unless all `n_vars` are present.
pub fn parse_smt_model(stdout: &str, n_vars: usize) -> Option<Vec<bool>> {
    let flat = stdout.replace(['(', ')'], " ");
    let tokens: Vec<&str> = flat.split_whitespace().collect();
    let mut model = vec![None; n_vars];
    let mut i = 0;
    while i + 1 < tokens.len() {
        if let Some(idx) = tokens[i].strip_prefix('x').and_then(|t| t.parse::<usize>().ok()) {
            let value = match tokens[i + 1] {
                "true" | "#b1" => Some(true),
                "false" | "#b0" => Some(false),
                _ => None,
            };
            if let Some(v) = value {
                if idx < n_vars {
                    model[idx] = Some(v);
                }
                i += 2;
                continue;
            }
        }
        i += 1;
    }
    model.into_iter().collect()
}

/// Classify solver output.
pub fn parse_smt_stdout(stdout: &str, n_vars: usize) -> SmtVerdict {
    for line in stdout.lines() {
        match line.trim() {
            "sat" => {
                return match parse_smt_model(stdout, n_vars) {
                    Some(model) => SmtVerdict::Sat(model),
                    None => SmtVerdict::Unknown,
                }
            }
            "unsat" => return SmtVerdict::Unsat,
            "unknown" | "timeout" => return SmtVerdict::Unknown,
            _ => {}
        }
    }
    SmtVerdict::Unknown
}

/// Pull `conflicts`/`decisions` counters out of whatever the solver printed.
///
/// Z3 `-st` prints `:conflicts N` and `:decisions N`, omitting counters
/// that are zero.  cvc5 1.4 drives CaDiCaL and exposes no conflict count;
/// `--stats-internal` prints `resource::steps::resource = { ...,
/// DecisionStep: N, ... }`, which is taken as its decision count.  Any
/// `*conflicts = N` or `*decisions = N` line from another solver or
/// version is accepted too.  The last matching number wins, which is the
/// final cumulative count in every known format.
pub fn harvest_counters(text: &str) -> (Option<u64>, Option<u64>) {
    let mut conflicts = None;
    let mut decisions = None;
    let flat = text.replace(['(', ')', ',', '=', '{', '}'], " ");
    let tokens: Vec<&str> = flat.split_whitespace().collect();
    for w in tokens.windows(2) {
        let key = w[0].trim_matches(':').to_ascii_lowercase();
        let Ok(value) = w[1].parse::<u64>() else {
            continue;
        };
        if key.ends_with("conflicts") {
            conflicts = Some(value);
        } else if key.ends_with("decisions") || key == "decisionstep" {
            decisions = Some(value);
        }
    }
    (conflicts, decisions)
}

/// Run the solver once on a script.
pub fn run_smt(script: &str, round: usize, options: &SmtSolveOptions) -> Result<SmtRun, String> {
    let path = match &options.keep_script {
        Some(path) => {
            let path = if round == 0 {
                path.clone()
            } else {
                path.with_extension(format!("{round}.smt2"))
            };
            if let Some(parent) = path.parent() {
                fs::create_dir_all(parent).map_err(|e| format!("create script parent: {e}"))?;
            }
            path
        }
        None => {
            let dir = options.work_dir.clone().unwrap_or_else(std::env::temp_dir);
            fs::create_dir_all(&dir).map_err(|e| format!("create work dir: {e}"))?;
            dir.join(format!("smt-pdp-{}-{round}.smt2", std::process::id()))
        }
    };
    fs::write(&path, script).map_err(|e| format!("write script {}: {e}", path.display()))?;

    let args = options
        .args
        .clone()
        .unwrap_or_else(|| options.kind.default_args());
    let started = Instant::now();
    let mut command = Command::new(&options.binary);
    command
        .args(&args)
        .arg(&path)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    if let Some(dir) = &options.work_dir {
        command.current_dir(dir);
    }
    let mut child = command
        .spawn()
        .map_err(|e| format!("spawn {} {}: {e}", options.kind.name(), options.binary.display()))?;
    let deadline = started + options.timeout;
    let timed_out = loop {
        match child.try_wait() {
            Ok(Some(_)) => break false,
            Ok(None) if Instant::now() >= deadline => {
                let _ = child.kill();
                break true;
            }
            Ok(None) => std::thread::sleep(Duration::from_millis(5)),
            Err(e) => return Err(format!("wait {}: {e}", options.kind.name())),
        }
    };
    let output = child
        .wait_with_output()
        .map_err(|e| format!("collect {} output: {e}", options.kind.name()))?;
    let elapsed = started.elapsed();
    if options.keep_script.is_none() {
        let _ = fs::remove_file(&path);
    }
    if timed_out {
        return Err(format!(
            "{} timed out after {} ms",
            options.kind.name(),
            options.timeout.as_millis()
        ));
    }
    let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    // Solvers exit non-zero on `unknown` or on a parse error; the verdict
    // line decides, so only a missing verdict is an invocation failure.
    if !output.status.success()
        && !stdout
            .lines()
            .any(|l| matches!(l.trim(), "sat" | "unsat" | "unknown"))
    {
        return Err(format!(
            "{} exited with status {}: {}",
            options.kind.name(),
            output.status,
            stderr.lines().take(3).collect::<Vec<_>>().join(" | ")
        ));
    }
    let (mut conflicts, mut decisions) = harvest_counters(&format!("{stdout}\n{stderr}"));
    // Z3 prints its statistics block only with `-st` and leaves out every
    // zero counter, so a block without `:conflicts` means none occurred.
    if options.kind == SmtSolverKind::Z3 && stdout.contains("(:") {
        conflicts.get_or_insert(0);
        decisions.get_or_insert(0);
    }
    // cvc5 1.4 exposes no SAT conflict count (its `*Conflict` keys are
    // theory-inference tallies), and on the bit-vector path the search
    // runs inside CaDiCaL, outside the `DecisionStep` resource counter.
    if options.kind == SmtSolverKind::Cvc5 {
        conflicts = None;
        if options.encoding == SmtEncoding::BitVec1 {
            decisions = None;
        }
    }
    Ok(SmtRun {
        stdout,
        stderr,
        elapsed,
        conflicts,
        decisions,
    })
}

/// Per-call accounting beyond [`SatDecompositionStats`].
#[derive(Clone, Debug, Default)]
pub struct SmtDecompositionReport {
    /// Solver calls (one per enumeration round).
    pub calls: usize,
    /// Wall time inside the solver, summed over rounds.
    pub solver_wall: Duration,
    /// Solver-reported conflicts summed over rounds, when every round printed one.
    pub conflicts: Option<u64>,
    /// Solver-reported decisions summed over rounds, when every round printed one.
    pub decisions: Option<u64>,
    /// Boolean equations handed to the solver (system, trace row, Macaulay rows).
    pub rows: usize,
    /// Solver variables.
    pub n_vars: usize,
    /// First invocation error, if the attempt ended in one.
    pub error: Option<String>,
}

/// Build the Boolean rows exactly as the native path does: the descended
/// system, the trace row when requested, then any Macaulay rows.
fn build_rows(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    sys: &DecompositionSystem,
    x_r: &F2mElement,
    m: usize,
    options: &SmtSolveOptions,
    stats: &mut SatDecompositionStats,
) -> Vec<F2BoolPoly> {
    let mut equations = sys.equations.clone();
    if options.trace_constraint {
        // Kosters–Yeo: Σ Tr(x_i) = Tr(x_R) + (m+1)·Tr(a) for rational points.
        let rhs = absolute_trace_bit(x_r, kc.n, &kc.curve.irreducible)
            ^ (m.is_multiple_of(2) && absolute_trace_bit(&kc.curve.a, kc.n, &kc.curve.irreducible));
        let mut terms = Vec::new();
        for (j, basis_element) in fb.subspace_basis.iter().enumerate() {
            if absolute_trace_bit(basis_element, kc.n, &kc.curve.irreducible) {
                for i in 0..m {
                    terms.push(F2BoolMono::var((i * fb.subspace_basis.len() + j) as u32));
                }
            }
        }
        if rhs {
            terms.push(F2BoolMono::one());
        }
        equations.push(F2BoolPoly::from_monos(terms, sys.n_vars));
    }
    if let Some(d) = options.macaulay_degree {
        if let Some(rows) = matrix_f4_f2(&sys.equations, sys.n_vars, d) {
            stats.implied_rows = rows.len();
            equations.extend(rows);
        }
    }
    equations
}

fn model_to_u64(model: &[bool]) -> Option<u64> {
    if model.len() > 64 {
        return None;
    }
    Some(
        model
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .fold(0u64, |acc, (i, _)| acc | (1u64 << i)),
    )
}

/// Solve one Semaev decomposition instance with an external SMT solver,
/// returning the same shape as `sat_decompose` plus the solver report.
///
/// Semantics follow the native path: models are enumerated with blocking
/// assertions up to `max_models`; a model is accepted only if it
/// satisfies the original equations and lifts to factor-base points that
/// sum to the target; a final `unsat` is a refutation; a timeout, an
/// invocation failure, `unknown`, or a spurious model ends the attempt
/// as exhausted.
pub fn smt_decompose_detailed(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    options: &SmtSolveOptions,
) -> (Option<Vec<usize>>, SatDecompositionStats, SmtDecompositionReport) {
    let mut stats = SatDecompositionStats::default();
    let mut report = SmtDecompositionReport::default();
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => {
            stats.exhausted = true;
            stats.unsupported = true;
            return (None, stats, report);
        }
    };
    let sys = match crate::cryptanalysis::polynomial_reuse::build_decomposition_system_reusing(
        &fb.subspace_basis,
        &x_r,
        &kc.curve.b,
        m,
        st,
    ) {
        Some(sys) => sys,
        None => {
            stats.exhausted = true;
            stats.unsupported = true;
            return (None, stats, report);
        }
    };
    let equations = build_rows(kc, fb, &sys, &x_r, m, options, &mut stats);
    let rows: Vec<AnfRow> = equations.iter().map(AnfRow::from_poly).collect();
    report.rows = rows.len();
    report.n_vars = sys.n_vars;

    let mut blocked: Vec<Vec<bool>> = Vec::new();
    let mut conflicts = Some(0u64);
    let mut decisions = Some(0u64);
    loop {
        if blocked.len() >= options.max_models {
            stats.exhausted = true;
            break;
        }
        let script = format_smtlib(sys.n_vars as u32, &rows, &blocked, options.encoding);
        stats.solver_calls += 1;
        report.calls += 1;
        let run = match run_smt(&script, blocked.len(), options) {
            Ok(run) => run,
            Err(e) => {
                report.error = Some(e);
                stats.exhausted = true;
                break;
            }
        };
        report.solver_wall += run.elapsed;
        conflicts = match (conflicts, run.conflicts) {
            (Some(a), Some(b)) => Some(a + b),
            _ => None,
        };
        decisions = match (decisions, run.decisions) {
            (Some(a), Some(b)) => Some(a + b),
            _ => None,
        };
        match parse_smt_stdout(&run.stdout, sys.n_vars) {
            SmtVerdict::Unsat => {
                stats.refuted = true;
                break;
            }
            SmtVerdict::Unknown => {
                stats.exhausted = true;
                break;
            }
            SmtVerdict::Sat(model) => {
                stats.models += 1;
                let Some(root) = model_to_u64(&model) else {
                    stats.exhausted = true;
                    break;
                };
                if !equations.iter().all(|e| e.eval(root) == 0) {
                    // Fail closed: a model that violates the rows handed to
                    // the solver is an encoding or parsing bug, not a near miss.
                    stats.spurious += 1;
                    stats.exhausted = true;
                    break;
                }
                let xs: Vec<F2mElement> = (0..m)
                    .map(|i| sys.summand_x(&fb.subspace_basis, root, i, kc.n))
                    .collect();
                if let Some(idxs) = lift_candidate(kc, fb, index_of, &xs, target) {
                    report.conflicts = conflicts;
                    report.decisions = decisions;
                    stats.conflicts = conflicts.unwrap_or(0);
                    return (Some(idxs), stats, report);
                }
                blocked.push(model);
            }
        }
    }
    report.conflicts = conflicts;
    report.decisions = decisions;
    stats.conflicts = conflicts.unwrap_or(0);
    (None, stats, report)
}

/// [`smt_decompose_detailed`] with the oracle-API return shape of
/// `sat_decompose` and `wdsat_decompose`.
pub fn smt_decompose(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
    options: &SmtSolveOptions,
) -> (Option<Vec<usize>>, SatDecompositionStats) {
    let (out, stats, _) = smt_decompose_detailed(kc, fb, index_of, st, target, m, options);
    (out, stats)
}

/// Locate a solver binary from `SMT_<KIND>_BIN` (e.g. `SMT_CVC5_BIN`), then `PATH`.
pub fn find_solver(kind: SmtSolverKind) -> Option<PathBuf> {
    let env_key = format!("SMT_{}_BIN", kind.name().to_ascii_uppercase());
    if let Ok(path) = std::env::var(&env_key) {
        let path = PathBuf::from(path);
        if path.is_file() {
            return Some(path);
        }
    }
    let exe = kind.name();
    std::env::var_os("PATH").and_then(|paths| {
        std::env::split_paths(&paths)
            .map(|dir| dir.join(exe))
            .find(|p| p.is_file())
    })
}

/// Evaluate ANF rows on a model; true when every row vanishes.
pub fn rows_satisfied(rows: &[AnfRow], model: &[bool]) -> bool {
    rows.iter().all(|row| {
        let mut acc = row.constant;
        for mono in &row.monomials {
            if mono.iter().all(|&v| model.get(v as usize).copied().unwrap_or(false)) {
                acc = !acc;
            }
        }
        !acc
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(monos: &[&[u32]], constant: bool) -> AnfRow {
        AnfRow {
            monomials: monos.iter().map(|m| m.to_vec()).collect(),
            constant,
        }
    }

    /// Every variable assignment of a small system, for brute-force checks.
    fn brute_force(rows: &[AnfRow], n_vars: usize) -> Vec<Vec<bool>> {
        (0u64..(1 << n_vars))
            .map(|code| (0..n_vars).map(|i| (code >> i) & 1 == 1).collect::<Vec<bool>>())
            .filter(|m| rows_satisfied(rows, m))
            .collect()
    }

    #[test]
    fn bool_script_has_one_assertion_per_nontrivial_row() {
        let rows = vec![
            row(&[&[0], &[1, 2]], true),
            row(&[], false), // 0 = 0, dropped
            row(&[&[2]], false),
        ];
        let script = format_smtlib(3, &rows, &[], SmtEncoding::Bool);
        assert!(script.contains("(set-logic QF_UF)"));
        assert_eq!(script.matches("(declare-const").count(), 3);
        assert_eq!(script.matches("(assert").count(), 2);
        assert!(script.contains("(assert (xor x0 (and x1 x2)))"));
        assert!(script.contains("(assert (not x2))"));
        assert!(script.contains("(get-value (x0 x1 x2))"));
    }

    #[test]
    fn bv_script_uses_one_bit_vectors() {
        let rows = vec![row(&[&[0], &[1, 2]], true), row(&[&[2]], false)];
        let script = format_smtlib(3, &rows, &[], SmtEncoding::BitVec1);
        assert!(script.contains("(set-logic QF_BV)"));
        assert!(script.contains("(declare-const x0 (_ BitVec 1))"));
        assert!(script.contains("(assert (= (bvxor x0 (bvand x1 x2)) #b1))"));
        assert!(script.contains("(assert (= x2 #b0))"));
    }

    #[test]
    fn inconsistent_constant_row_asserts_false() {
        let rows = vec![row(&[], true)];
        let script = format_smtlib(1, &rows, &[], SmtEncoding::Bool);
        assert!(script.contains("(assert false)"));
    }

    #[test]
    fn blocking_assertion_excludes_exactly_that_model() {
        let model = vec![true, false, true];
        let b = blocking_assertion(&model, SmtEncoding::Bool);
        assert_eq!(b, "(assert (not (and x0 (not x1) x2)))");
        let b = blocking_assertion(&model, SmtEncoding::BitVec1);
        assert_eq!(b, "(assert (not (and (= x0 #b1) (= x1 #b0) (= x2 #b1))))");
    }

    #[test]
    fn model_parser_reads_both_value_syntaxes_in_any_order() {
        let out = "sat\n((x2 true)\n (x0 false) (x1 true))\n";
        assert_eq!(parse_smt_model(out, 3), Some(vec![false, true, true]));
        let out = "sat\n((x0 #b1) (x1 #b0))\n";
        assert_eq!(parse_smt_model(out, 2), Some(vec![true, false]));
        assert_eq!(parse_smt_model("sat\n((x0 true))\n", 2), None);
    }

    #[test]
    fn verdict_follows_first_status_line() {
        assert_eq!(parse_smt_stdout("unsat\n", 2), SmtVerdict::Unsat);
        assert_eq!(parse_smt_stdout("unknown\n", 2), SmtVerdict::Unknown);
        assert_eq!(parse_smt_stdout("", 2), SmtVerdict::Unknown);
        assert_eq!(
            parse_smt_stdout("sat\n((x0 true) (x1 false))\n", 2),
            SmtVerdict::Sat(vec![true, false])
        );
        // `sat` without a readable model is not a model.
        assert_eq!(parse_smt_stdout("sat\n", 2), SmtVerdict::Unknown);
    }

    #[test]
    fn counters_are_harvested_from_z3_and_cvc5_formats() {
        let z3 = "sat\n((x0 true))\n(:conflicts 12\n :decisions 34\n :max-memory 5.1)\n";
        assert_eq!(harvest_counters(z3), (Some(12), Some(34)));
        let cvc5 = "sat\nresource::steps::resource = { BvSatStep: 1, CnfStep: 5, DecisionStep: 9, RewriteStep: 124 }\n";
        assert_eq!(harvest_counters(cvc5), (None, Some(9)));
        let other = "sat\nSatSolver::conflicts = 7\nSatSolver::decisions = 9\n";
        assert_eq!(harvest_counters(other), (Some(7), Some(9)));
        assert_eq!(harvest_counters("sat\n"), (None, None));
    }

    #[test]
    fn rows_satisfied_matches_anf_semantics() {
        // x0 + x1·x2 + 1 = 0  and  x2 = 0  ⇒  x0 = 1, x2 = 0, x1 free.
        let rows = vec![row(&[&[0], &[1, 2]], true), row(&[&[2]], false)];
        let sols = brute_force(&rows, 3);
        assert_eq!(sols.len(), 2);
        assert!(sols.iter().all(|m| m[0] && !m[2]));
    }

    /// With a real solver on `PATH` or `SMT_<KIND>_BIN`, the script round-trips:
    /// the solver's model satisfies the rows, blocking enumerates every
    /// solution, and the final answer is `unsat`.
    fn solver_round_trip(kind: SmtSolverKind, encoding: SmtEncoding) {
        let Some(binary) = find_solver(kind) else {
            eprintln!("{} not available; skipping", kind.name());
            return;
        };
        let rows = vec![
            row(&[&[0], &[1, 2]], true),
            row(&[&[2], &[3]], false),
            row(&[&[1, 3], &[0, 3]], false),
        ];
        let expected = brute_force(&rows, 4);
        assert!(!expected.is_empty());
        let mut options = SmtSolveOptions::new(kind, binary);
        options.encoding = encoding;
        options.timeout = Duration::from_secs(60);
        let mut blocked = Vec::new();
        loop {
            let script = format_smtlib(4, &rows, &blocked, encoding);
            let run = run_smt(&script, blocked.len(), &options).expect("solver runs");
            match parse_smt_stdout(&run.stdout, 4) {
                SmtVerdict::Sat(model) => {
                    assert!(rows_satisfied(&rows, &model), "solver model violates rows");
                    assert!(!blocked.contains(&model), "blocked model returned again");
                    blocked.push(model);
                    assert!(blocked.len() <= expected.len(), "more models than solutions");
                }
                SmtVerdict::Unsat => break,
                SmtVerdict::Unknown => panic!("solver answered unknown: {}", run.stderr),
            }
        }
        blocked.sort();
        let mut expected = expected;
        expected.sort();
        assert_eq!(blocked, expected, "enumeration must equal brute force");
    }

    #[test]
    fn cvc5_round_trip_bool() {
        solver_round_trip(SmtSolverKind::Cvc5, SmtEncoding::Bool);
    }

    #[test]
    fn cvc5_round_trip_bv1() {
        solver_round_trip(SmtSolverKind::Cvc5, SmtEncoding::BitVec1);
    }

    #[test]
    fn z3_round_trip_bool() {
        solver_round_trip(SmtSolverKind::Z3, SmtEncoding::Bool);
    }

    #[test]
    fn z3_round_trip_bv1() {
        solver_round_trip(SmtSolverKind::Z3, SmtEncoding::BitVec1);
    }
}
