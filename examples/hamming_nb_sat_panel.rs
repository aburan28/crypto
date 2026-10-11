//! Hamming-weight factor base in a type-II optimal normal basis, decided by
//! off-the-shelf SAT solvers: the native replacement for the Python
//! `ecc2k130/codegen` route and its measured `n = 9, 11, 23` rungs.
//!
//! Factor base `F_w = {P : wt_ONB(x(P)) ≤ w}` on `K_0 : y² + xy = x³ + 1`
//! over `F_{2^n}` in the type-II ONB (Frobenius-stable at every `n`, no
//! storage, membership by popcount).  The decomposition question for
//! `m = 2` is `S₃(x₁, x₂, x_R) = 0`; for `m = 3` the chain
//! `S₃(x₁, x₂, u) = 0`, `S₃(u, x₃, x_R) = 0` with `u ∈ F_{2^n}` free.
//! Squaring is a coordinate permutation and the ONB product has at most
//! `2n − 1` terms per coordinate, so each `S₃` link is `≤ 4n²` AND gates and
//! `O(n)` parity rows; the weight bound is a Sinz sequential counter and the
//! summand symmetry is broken lexicographically.
//!
//! Every model is lifted through the group: a model whose abscissae do not
//! all lift, or whose signed sums miss `±R`, is blocked and the solver is
//! re-run, so a reported relation is a verified point relation and a final
//! UNSAT is a refutation of the algebraic system (not merely of the point
//! question).  An exact pair-table oracle decides every target first and
//! supplies the truth and the enumeration floor.
//!
//! ```text
//! cargo run --release --example hamming_nb_sat_panel -- \
//!     --n 11 --w 3 --m 3 --targets 6 --planted 2 --seed 20261011 \
//!     --solvers cms,cadical --timeout-s 300 --out results/n11-w3-m3.jsonl
//! ```
//!
//! Stage diagnostic: no `S`, no end-to-end cost, no `vs_rho` figure.

use crypto_lib::cryptanalysis::koblitz_rotated_chain::{OnbCurve, OnbPoint, Type2Onb};
use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::io::Write as _;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

// ---------------------------------------------------------------- CNF builder

type Lit = i64; // DIMACS literal, never 0

#[derive(Default)]
struct Cnf {
    vars: i64,
    clauses: Vec<Vec<Lit>>,
    /// Parity rows: XOR of literals == rhs.
    xors: Vec<(Vec<Lit>, bool)>,
    ands: usize,
}

impl Cnf {
    fn var(&mut self) -> Lit {
        self.vars += 1;
        self.vars
    }
    fn clause(&mut self, c: Vec<Lit>) {
        self.clauses.push(c);
    }
    fn and(&mut self, a: Lit, b: Lit) -> Lit {
        let g = self.var();
        self.clause(vec![-g, a]);
        self.clause(vec![-g, b]);
        self.clause(vec![g, -a, -b]);
        self.ands += 1;
        g
    }
    fn or(&mut self, a: Lit, b: Lit) -> Lit {
        -self.and(-a, -b)
    }
    fn xor2(&mut self, a: Lit, b: Lit) -> Lit {
        let g = self.var();
        self.xors.push((vec![g, a, b], false));
        g
    }
    /// Fresh variable `t` with `t = XOR(lits) ⊕ constant`.
    fn define_parity(&mut self, lits: &[Lit], constant: bool) -> Lit {
        let t = self.var();
        let mut row = lits.to_vec();
        row.push(t);
        self.xors.push((row, constant));
        t
    }
    /// Sinz sequential counter: at most `w` of `xs` true.
    fn at_most(&mut self, xs: &[Lit], w: usize) {
        let n = xs.len();
        if w >= n {
            return;
        }
        if w == 0 {
            for &x in xs {
                self.clause(vec![-x]);
            }
            return;
        }
        let s: Vec<Vec<Lit>> = (0..n - 1)
            .map(|_| (0..w).map(|_| self.var()).collect())
            .collect();
        self.clause(vec![-xs[0], s[0][0]]);
        for j in 1..w {
            self.clause(vec![-s[0][j]]);
        }
        for i in 1..n - 1 {
            self.clause(vec![-xs[i], s[i][0]]);
            self.clause(vec![-s[i - 1][0], s[i][0]]);
            for j in 1..w {
                self.clause(vec![-xs[i], -s[i - 1][j - 1], s[i][j]]);
                self.clause(vec![-s[i - 1][j], s[i][j]]);
            }
            self.clause(vec![-xs[i], -s[i - 1][w - 1]]);
        }
        self.clause(vec![-xs[n - 1], -s[n - 2][w - 1]]);
    }
    /// `a ≤ b` as unsigned integers with bit `n−1` most significant.
    fn lex_le(&mut self, a: &[Lit], b: &[Lit]) {
        let mut gt: Option<Lit> = None;
        let mut eq: Option<Lit> = None;
        for i in (0..a.len()).rev() {
            let a_gt_b = self.and(a[i], -b[i]);
            let term = match eq {
                Some(e) => self.and(e, a_gt_b),
                None => a_gt_b,
            };
            gt = Some(match gt {
                Some(g) => self.or(g, term),
                None => term,
            });
            let same = -self.xor2(a[i], b[i]);
            eq = Some(match eq {
                Some(e) => self.and(e, same),
                None => same,
            });
        }
        if let Some(g) = gt {
            self.clause(vec![-g]);
        }
    }
    fn write_dimacs(&self, path: &PathBuf, native_xor: bool) -> (usize, usize) {
        let mut out = String::new();
        let mut clauses = self.clauses.clone();
        let mut extra_vars = 0i64;
        let mut xor_lines = Vec::new();
        if native_xor {
            for (row, rhs) in &self.xors {
                // CryptoMiniSat: `x` line asserts XOR(lits) = true; negate one literal for rhs=false.
                let mut r = row.clone();
                if !*rhs {
                    r[0] = -r[0];
                }
                xor_lines.push(r);
            }
        } else {
            // Tseitin-split each parity row into 3-ary XORs.
            for (row, rhs) in &self.xors {
                let mut acc: Vec<Lit> = row.clone();
                while acc.len() > 3 {
                    let a = acc.pop().unwrap();
                    let b = acc.pop().unwrap();
                    extra_vars += 1;
                    let t = self.vars + extra_vars;
                    // t = a xor b
                    clauses.push(vec![-t, a, b]);
                    clauses.push(vec![-t, -a, -b]);
                    clauses.push(vec![t, -a, b]);
                    clauses.push(vec![t, a, -b]);
                    acc.push(t);
                }
                // XOR(acc) = rhs, |acc| in 1..=3: forbid every assignment of the wrong parity.
                let k = acc.len();
                for mask in 0..(1u32 << k) {
                    let parity = mask.count_ones() % 2 == 1;
                    if parity != *rhs {
                        let c: Vec<Lit> = (0..k)
                            .map(|i| if mask >> i & 1 == 1 { -acc[i] } else { acc[i] })
                            .collect();
                        clauses.push(c);
                    }
                }
            }
        }
        let nvars = self.vars + extra_vars;
        let _ = writeln!(out, "p cnf {} {}", nvars, clauses.len() + xor_lines.len());
        for c in &clauses {
            for l in c {
                let _ = write!(out, "{l} ");
            }
            out.push_str("0\n");
        }
        for r in &xor_lines {
            out.push('x');
            for l in r {
                let _ = write!(out, " {l}");
            }
            out.push_str(" 0\n");
        }
        fs::write(path, out).expect("write dimacs");
        (nvars as usize, clauses.len() + xor_lines.len())
    }
}

// -------------------------------------------------- symbolic ONB arithmetic

/// A Boolean affine form: XOR of literals (positive variables) plus a constant.
#[derive(Clone, Default)]
struct Lin {
    vars: Vec<Lit>,
    c: bool,
}

impl Lin {
    fn constant(c: bool) -> Self {
        Lin { vars: vec![], c }
    }
    fn var(v: Lit) -> Self {
        Lin {
            vars: vec![v],
            c: false,
        }
    }
    fn add(&self, o: &Lin) -> Lin {
        let mut vars = self.vars.clone();
        for &v in &o.vars {
            if let Some(pos) = vars.iter().position(|&x| x == v) {
                vars.swap_remove(pos);
            } else {
                vars.push(v);
            }
        }
        Lin {
            vars,
            c: self.c ^ o.c,
        }
    }
}

struct SymField {
    f: Type2Onb,
    n: usize,
    /// `perm[i]` = coordinate index that coordinate `i` moves to under Frobenius.
    perm: Vec<usize>,
}

impl SymField {
    fn new(f: Type2Onb) -> Self {
        let n = f.degree();
        let perm = (0..n)
            .map(|i| {
                let img = f.frobenius(1u128 << i, 1);
                img.trailing_zeros() as usize
            })
            .collect();
        SymField { f, n, perm }
    }
    fn fold(&self, index: usize) -> usize {
        let ring = 2 * self.n + 1;
        let index = index % ring;
        index.min(ring - index)
    }
    fn constant(&self, v: u128) -> Vec<Lin> {
        (0..self.n)
            .map(|i| Lin::constant(v >> i & 1 == 1))
            .collect()
    }
    fn inputs(&self, cnf: &mut Cnf) -> Vec<Lin> {
        (0..self.n).map(|_| Lin::var(cnf.var())).collect()
    }
    fn add(&self, a: &[Lin], b: &[Lin]) -> Vec<Lin> {
        a.iter().zip(b).map(|(x, y)| x.add(y)).collect()
    }
    fn square(&self, a: &[Lin]) -> Vec<Lin> {
        let mut out = vec![Lin::default(); self.n];
        for (i, l) in a.iter().enumerate() {
            out[self.perm[i]] = l.clone();
        }
        out
    }
    /// Materialise every coordinate as a single literal (or a constant).
    fn materialise(&self, cnf: &mut Cnf, a: &[Lin]) -> Vec<Option<Lit>> {
        a.iter()
            .map(|l| {
                if l.vars.is_empty() {
                    None
                } else if l.vars.len() == 1 && !l.c {
                    Some(l.vars[0])
                } else {
                    Some(cnf.define_parity(&l.vars, l.c))
                }
            })
            .collect()
    }
    /// ONB product.  Coordinates are 1-indexed in the fold rule, as in `Type2Onb::mul`.
    fn mul(&self, cnf: &mut Cnf, a: &[Lin], b: &[Lin]) -> Vec<Lin> {
        let ma = self.materialise(cnf, a);
        let mb = self.materialise(cnf, b);
        let mut out = vec![Lin::default(); self.n];
        for i0 in 0..self.n {
            let a_bit = (ma[i0], a[i0].c && a[i0].vars.is_empty());
            for j0 in 0..self.n {
                let b_bit = (mb[j0], b[j0].c && b[j0].vars.is_empty());
                // term = a_i * b_j as a Lin
                let term: Lin = match (a_bit, b_bit) {
                    ((None, false), _) | (_, (None, false)) => continue, // zero
                    ((None, true), (None, true)) => Lin::constant(true),
                    ((None, true), (Some(v), _)) | ((Some(v), _), (None, true)) => Lin::var(v),
                    ((Some(u), _), (Some(v), _)) => Lin::var(cnf.and(u, v)),
                };
                let i = i0 + 1;
                let j = j0 + 1;
                let s = self.fold(i + j);
                if s != 0 {
                    out[s - 1] = out[s - 1].add(&term);
                }
                let d = self.fold(i.abs_diff(j));
                if d != 0 {
                    out[d - 1] = out[d - 1].add(&term);
                }
            }
        }
        out
    }
    /// `S₃(a, b, c) = (ab + ac + bc)² + abc + 1` for `K_0` (`b_curve = 1`).
    fn s3(&self, cnf: &mut Cnf, a: &[Lin], b: &[Lin], c: &[Lin]) -> Vec<Lin> {
        let ab = self.mul(cnf, a, b);
        let ac = self.mul(cnf, a, c);
        let bc = self.mul(cnf, b, c);
        let sum = self.add(&self.add(&ab, &ac), &bc);
        let sq = self.square(&sum);
        let abc = self.mul(cnf, &ab, c);
        let one = self.constant(self.f.one());
        self.add(&self.add(&sq, &abc), &one)
    }
    fn assert_zero(&self, cnf: &mut Cnf, a: &[Lin]) {
        for l in a {
            if l.vars.is_empty() {
                assert!(!l.c, "constant 1 = 0: inconsistent constraint");
                continue;
            }
            cnf.xors.push((l.vars.clone(), l.c));
        }
    }
}

// ------------------------------------------------------------- the instance

struct Instance {
    cnf: Cnf,
    /// Input literals of each summand abscissa, coordinate order.
    summands: Vec<Vec<Lit>>,
}

fn build_instance(sf: &SymField, m: usize, w: usize, x_r: u128, symmetry: bool) -> Instance {
    let mut cnf = Cnf::default();
    let xs: Vec<Vec<Lin>> = (0..m).map(|_| sf.inputs(&mut cnf)).collect();
    let target = sf.constant(x_r);
    match m {
        2 => {
            let s = sf.s3(&mut cnf, &xs[0], &xs[1], &target);
            sf.assert_zero(&mut cnf, &s);
        }
        3 => {
            let u = sf.inputs(&mut cnf);
            let s1 = sf.s3(&mut cnf, &xs[0], &xs[1], &u);
            sf.assert_zero(&mut cnf, &s1);
            let s2 = sf.s3(&mut cnf, &u, &xs[2], &target);
            sf.assert_zero(&mut cnf, &s2);
        }
        _ => panic!("m must be 2 or 3"),
    }
    let summands: Vec<Vec<Lit>> = xs
        .iter()
        .map(|x| x.iter().map(|l| l.vars[0]).collect())
        .collect();
    for s in &summands {
        cnf.at_most(s, w);
    }
    if symmetry {
        for i in 1..m {
            cnf.lex_le(&summands[i - 1], &summands[i]);
        }
    }
    Instance { cnf, summands }
}

// ----------------------------------------------------------- reference oracle

struct Reference {
    base: Vec<OnbPoint>,
    /// x(P_i + P_j) for i < j  ->  number of unordered pairs with that abscissa.
    pair_x: HashMap<u128, u32>,
}

fn weight_le(x: u128, w: usize) -> bool {
    (x.count_ones() as usize) <= w
}

fn enumerate_weight(n: usize, w: usize) -> Vec<u128> {
    let mut out = vec![0u128];
    let mut cur: Vec<u128> = vec![0];
    for _ in 0..w {
        let mut next = Vec::new();
        for &v in &cur {
            let top = if v == 0 {
                0
            } else {
                128 - v.leading_zeros() as usize
            };
            for b in top..n {
                next.push(v | (1u128 << b));
            }
        }
        out.extend(next.iter().copied());
        cur = next;
    }
    out
}

impl Reference {
    fn new(curve: OnbCurve, n: usize, w: usize) -> Self {
        let mut base = Vec::new();
        for x in enumerate_weight(n, w) {
            for p in curve.points_with_x(x) {
                base.push(p);
            }
        }
        let mut pair_x = HashMap::new();
        for i in 0..base.len() {
            for j in i + 1..base.len() {
                let s = curve.add(base[i], base[j]);
                if !s.infinity {
                    *pair_x.entry(s.x).or_insert(0) += 1;
                }
            }
        }
        Reference { base, pair_x }
    }
    /// Number of decompositions `R = Σ P_i` with `P_i ∈ F` (points, unordered, distinct).
    fn count(&self, curve: OnbCurve, m: usize, r: OnbPoint, w: usize) -> u64 {
        match m {
            2 => {
                let mut c = 0u64;
                for &p in &self.base {
                    let q = curve.add(r, curve.neg(p));
                    if !q.infinity && weight_le(q.x, w) && q != p {
                        c += 1;
                    }
                }
                c / 2
            }
            3 => {
                let mut c = 0u64;
                for &p in &self.base {
                    let q = curve.add(r, curve.neg(p));
                    if !q.infinity {
                        c += *self.pair_x.get(&q.x).unwrap_or(&0) as u64;
                    }
                }
                // Each unordered triple is counted once per choice of the singled-out
                // summand (3 times) and the pair table counts both Q and −Q (the pair
                // sum must equal Q exactly; −Q pairs are a different triple summing to a
                // different point unless …). Keep the raw count: it is an upper bound used
                // only for the decomposable / not-decomposable verdict.
                c
            }
            _ => unreachable!(),
        }
    }
    /// Exact verdict: does any multiset of `m` factor-base points sum to `r`?
    fn decomposable(&self, curve: OnbCurve, m: usize, r: OnbPoint, w: usize) -> bool {
        match m {
            2 => self.base.iter().any(|&p| {
                let q = curve.add(r, curve.neg(p));
                !q.infinity && weight_le(q.x, w)
            }),
            3 => self.base.iter().any(|&p| {
                let q = curve.add(r, curve.neg(p));
                !q.infinity && {
                    // need P_a + P_b = q exactly
                    self.base.iter().any(|&a| {
                        let b = curve.add(q, curve.neg(a));
                        !b.infinity && weight_le(b.x, w)
                    })
                }
            }),
            _ => unreachable!(),
        }
    }
}

// ------------------------------------------------------------------ solvers

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Solver {
    Cms,
    Cadical,
}

impl Solver {
    fn name(self) -> &'static str {
        match self {
            Solver::Cms => "cryptominisat5",
            Solver::Cadical => "cadical",
        }
    }
    fn native_xor(self) -> bool {
        matches!(self, Solver::Cms)
    }
}

struct SolveOutcome {
    status: &'static str, // sat | unsat | timeout | error
    model: Vec<bool>,
    conflicts: Option<u64>,
    wall: Duration,
}

fn run_solver(solver: Solver, path: &PathBuf, timeout: Duration, nvars: usize) -> SolveOutcome {
    let start = Instant::now();
    let mut cmd = Command::new("timeout");
    cmd.arg(format!("{}", timeout.as_secs() + 1));
    match solver {
        Solver::Cms => {
            cmd.arg("cryptominisat5").arg("--verb").arg("1").arg(path);
        }
        Solver::Cadical => {
            cmd.arg("cadical")
                .arg("-t")
                .arg(format!("{}", timeout.as_secs()))
                .arg(path);
        }
    }
    let output = cmd
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("spawn solver");
    let wall = start.elapsed();
    let text = String::from_utf8_lossy(&output.stdout).to_string();
    let mut status = "error";
    let mut model = vec![false; nvars + 1];
    let mut conflicts = None;
    for line in text.lines() {
        if line.starts_with("s SATISFIABLE") {
            status = "sat";
        } else if line.starts_with("s UNSATISFIABLE") {
            status = "unsat";
        } else if line.starts_with("s UNKNOWN") || line.starts_with("c UNKNOWN") {
            status = "timeout";
        } else if let Some(rest) = line.strip_prefix("v ") {
            for tok in rest.split_whitespace() {
                let v: i64 = tok.parse().unwrap_or(0);
                if v > 0 && (v as usize) < model.len() {
                    model[v as usize] = true;
                }
            }
        } else if line.starts_with("c conflicts") || line.starts_with("c   conflicts") {
            // cadical: "c conflicts:   1234   56.7  per second"; cms: "c conflicts  : 1234 ..."
            let digits: String = line
                .split(':')
                .nth(1)
                .unwrap_or("")
                .trim()
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect();
            if let Ok(v) = digits.parse::<u64>() {
                conflicts = Some(v);
            }
        }
    }
    if status == "error" && wall >= timeout {
        status = "timeout";
    }
    SolveOutcome {
        status,
        model,
        conflicts,
        wall,
    }
}

// --------------------------------------------------------------------- main

struct Args {
    n: usize,
    w: usize,
    m: usize,
    targets: usize,
    planted: usize,
    seed: u64,
    solvers: Vec<Solver>,
    timeout: Duration,
    max_models: usize,
    symmetry: bool,
    out: Option<String>,
    work: PathBuf,
}

fn parse_args() -> Args {
    let mut a = Args {
        n: 11,
        w: 3,
        m: 3,
        targets: 6,
        planted: 2,
        seed: 20261011,
        solvers: vec![Solver::Cms, Solver::Cadical],
        timeout: Duration::from_secs(300),
        max_models: 64,
        symmetry: true,
        out: None,
        work: std::env::temp_dir().join("hamming_nb_sat_panel"),
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let v = argv.get(i + 1).cloned().unwrap_or_default();
        match argv[i].as_str() {
            "--n" => a.n = v.parse().unwrap(),
            "--w" => a.w = v.parse().unwrap(),
            "--m" => a.m = v.parse().unwrap(),
            "--targets" => a.targets = v.parse().unwrap(),
            "--planted" => a.planted = v.parse().unwrap(),
            "--seed" => a.seed = v.parse().unwrap(),
            "--timeout-s" => a.timeout = Duration::from_secs(v.parse().unwrap()),
            "--max-models" => a.max_models = v.parse().unwrap(),
            "--no-symmetry" => {
                a.symmetry = false;
                i += 1;
                continue;
            }
            "--solvers" => {
                a.solvers = v
                    .split(',')
                    .map(|s| match s {
                        "cms" => Solver::Cms,
                        "cadical" => Solver::Cadical,
                        _ => panic!("unknown solver {s}"),
                    })
                    .collect()
            }
            "--out" => a.out = Some(v),
            "--work" => a.work = PathBuf::from(v),
            other => panic!("unknown argument {other}"),
        }
        i += 2;
    }
    a
}

/// SplitMix64, so the panel needs no crate beyond the library.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
}

fn random_point(curve: OnbCurve, n: usize, rng: &mut Rng) -> OnbPoint {
    loop {
        let x = (rng.next() as u128) & ((1u128 << n) - 1);
        let pts = curve.points_with_x(x);
        if !pts.is_empty() {
            return pts[(rng.next() & 1) as usize % pts.len()];
        }
    }
}

fn verify_model(
    curve: OnbCurve,
    m: usize,
    xs: &[u128],
    r: OnbPoint,
    w: usize,
) -> Option<Vec<OnbPoint>> {
    let mut lifts: Vec<Vec<OnbPoint>> = Vec::new();
    for &x in xs {
        if !weight_le(x, w) {
            return None;
        }
        let p = curve.points_with_x(x);
        if p.is_empty() {
            return None;
        }
        lifts.push(p);
    }
    // Try every choice of sign (point) per summand.
    let choices = 1usize << m;
    for c in 0..choices {
        let mut sum = OnbPoint::INFINITY;
        let mut pts = Vec::new();
        for i in 0..m {
            let p = lifts[i][(c >> i & 1) % lifts[i].len()];
            pts.push(p);
            sum = curve.add(sum, p);
        }
        if sum == r {
            return Some(pts);
        }
    }
    None
}

fn main() {
    let args = parse_args();
    let field = Type2Onb::new(args.n).expect("type-II ONB");
    let curve = OnbCurve::new(field);
    let sf = SymField::new(field);
    fs::create_dir_all(&args.work).unwrap();
    let mut out: Box<dyn std::io::Write> = match &args.out {
        Some(p) => Box::new(
            fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(p)
                .unwrap(),
        ),
        None => Box::new(std::io::stdout()),
    };
    let t0 = Instant::now();
    let reference = Reference::new(curve, args.n, args.w);
    let ref_build = t0.elapsed();
    let base_size = reference.base.len();
    let _ = writeln!(
        out,
        "{{\"kind\":\"header\",\"n\":{},\"w\":{},\"m\":{},\"base_points\":{},\"pair_table_entries\":{},\"ref_build_ms\":{},\"targets\":{},\"planted\":{},\"seed\":{},\"timeout_s\":{},\"max_models\":{},\"symmetry\":{},\"solvers\":{:?}}}",
        args.n, args.w, args.m, base_size, reference.pair_x.len(), ref_build.as_millis(), args.targets, args.planted, args.seed, args.timeout.as_secs(), args.max_models, args.symmetry,
        args.solvers.iter().map(|s| s.name()).collect::<Vec<_>>()
    );
    let mut rng = Rng(args.seed);
    let mut targets: Vec<(&'static str, OnbPoint)> = Vec::new();
    for _ in 0..args.targets {
        targets.push(("uniform", random_point(curve, args.n, &mut rng)));
    }
    for _ in 0..args.planted {
        let mut sum = OnbPoint::INFINITY;
        for _ in 0..args.m {
            let p = reference.base[(rng.next() as usize) % base_size];
            sum = curve.add(sum, p);
        }
        if sum.infinity {
            continue;
        }
        targets.push(("planted", sum));
    }
    for (ti, (stratum, r)) in targets.iter().enumerate() {
        let t = Instant::now();
        let truth = reference.decomposable(curve, args.m, *r, args.w);
        let count = reference.count(curve, args.m, *r, args.w);
        let ref_wall = t.elapsed();
        let inst = build_instance(&sf, args.m, args.w, r.x, args.symmetry);
        let floor = match args.m {
            2 => base_size as u64,
            _ => (base_size as u64) * (base_size as u64 - 1) / 2,
        };
        for &solver in &args.solvers {
            let path = args.work.join(format!(
                "n{}-w{}-m{}-t{}-{}.cnf",
                args.n,
                args.w,
                args.m,
                ti,
                solver.name()
            ));
            let mut cnf_clauses = inst.cnf.clauses.clone();
            let mut cnf = Cnf {
                vars: inst.cnf.vars,
                clauses: std::mem::take(&mut cnf_clauses),
                xors: inst.cnf.xors.clone(),
                ands: inst.cnf.ands,
            };
            let mut total_wall = Duration::ZERO;
            let mut total_conflicts = 0u64;
            let mut conflicts_known = true;
            let mut blocked = 0usize;
            let mut relation: Option<Vec<u128>> = None;
            let mut calls = 0usize;
            let (nvars, nclauses) = cnf.write_dimacs(&path, solver.native_xor());
            let verdict = loop {
                if total_wall >= args.timeout || blocked >= args.max_models {
                    break if blocked >= args.max_models {
                        "model_budget"
                    } else {
                        "timeout"
                    };
                }
                let remaining = args.timeout - total_wall;
                let o = run_solver(solver, &path, remaining, nvars);
                calls += 1;
                total_wall += o.wall;
                match o.conflicts {
                    Some(c) => total_conflicts += c,
                    None => conflicts_known = false,
                }
                match o.status {
                    "unsat" => break "unsat",
                    "sat" => {
                        let xs: Vec<u128> = inst
                            .summands
                            .iter()
                            .map(|lits| {
                                let mut v = 0u128;
                                for (i, &l) in lits.iter().enumerate() {
                                    if o.model[l as usize] {
                                        v |= 1u128 << i;
                                    }
                                }
                                v
                            })
                            .collect();
                        if verify_model(curve, args.m, &xs, *r, args.w).is_some() {
                            relation = Some(xs);
                            break "found";
                        }
                        // Block this abscissa assignment and continue.
                        blocked += 1;
                        let mut clause = Vec::new();
                        for lits in &inst.summands {
                            for &l in lits {
                                clause.push(if o.model[l as usize] { -l } else { l });
                            }
                        }
                        cnf.clauses.push(clause);
                        cnf.write_dimacs(&path, solver.native_xor());
                    }
                    "timeout" => break "timeout",
                    _ => break "error",
                }
            };
            let agrees = match verdict {
                "found" => truth,
                "unsat" => !truth,
                _ => true,
            };
            let _ = writeln!(
                out,
                "{{\"kind\":\"cell\",\"n\":{},\"w\":{},\"m\":{},\"target\":{},\"stratum\":\"{}\",\"x_r\":\"0x{:x}\",\"truth_decomposable\":{},\"ref_count\":{},\"ref_wall_us\":{},\"enumeration_floor\":{},\"solver\":\"{}\",\"vars\":{},\"clauses\":{},\"and_gates\":{},\"parity_rows\":{},\"verdict\":\"{}\",\"agrees\":{},\"calls\":{},\"blocked_models\":{},\"conflicts\":{},\"wall_ms\":{},\"relation\":{}}}",
                args.n, args.w, args.m, ti, stratum, r.x, truth, count, ref_wall.as_micros(), floor, solver.name(), nvars, nclauses, inst.cnf.ands, inst.cnf.xors.len(), verdict, agrees, calls, blocked,
                if conflicts_known { total_conflicts.to_string() } else { "null".to_string() },
                total_wall.as_millis(),
                relation.map(|xs| format!("[{}]", xs.iter().map(|x| format!("\"0x{x:x}\"")).collect::<Vec<_>>().join(","))).unwrap_or_else(|| "null".to_string())
            );
            if !agrees {
                eprintln!("DISAGREEMENT at target {ti} solver {}", solver.name());
            }
        }
    }
}
