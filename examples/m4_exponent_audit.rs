//! The **`m = 4` exponent audit** ladder of
//! `research/ic_m4_exponent_audit_20260928/PREREGISTRATION.md`
//! (survey: `research/ic_candidate_tournament_20260915/campaign_20260916/DECOMPOSITION-SURVEY.md` §5).
//!
//! One process measures one cell `(K_a, n)` on one arm and prints one JSON
//! line per target on stdout (and to `--out FILE`, which is never
//! overwritten).  The cell's factor base is `F = {P : x(P) ∈ V}` with `V` a
//! uniformly random `ℓ`-dimensional `F_2`-subspace, `ℓ = round(n/4)`, drawn
//! exactly as `examples/dreg_ladder.rs` draws it (`random_subspace_basis`,
//! reproduced below verbatim so this file also builds against the frozen
//! source commit, which predates that library function).  Targets are
//! uniform (natural, unplanted) points `[k]G` of the prime-order subgroup.
//!
//! Arms (`--arm`):
//!
//! - `semaev`    — the chained `m = 4` Semaev system through
//!   `koblitz_index_calculus::groebner_decompose`, the oracle
//!   `examples/groebner_stage_bench.rs` measures, with its
//!   20,000-node budget.  Metric: 64-bit word XORs in the
//!   Macaulay eliminations for that target (deterministic).
//!   Each target's ground truth (does an `m = 4` decomposition
//!   over `F` exist?) comes from exhaustive enumeration,
//!   uncharged, and is compared with the oracle's verdict.
//! - `enumerate` — the enumeration null: non-decreasing index tuples,
//!   depth first, stopping at the first decomposition, exactly
//!   the library's `decompose`.  Metric: point additions.
//! - `null`      — the random-system null: per target, a random Boolean
//!   system (`koblitz_bench::random_control_system`) with the
//!   target's own Semaev system's unknown count, equation
//!   count, degree and mean terms per equation, solved by the
//!   same engine, options and node budget with every root
//!   rejected (a full-tree search, the refutation analogue).
//!   Metric: word XORs.
//! - `degree`    — the secondary check: exact Boolean solution count of the
//!   target's chained system, then, on systems with none, the
//!   refutation degree (`solving_degree`, natural layout) up to
//!   `--d-max`.
//!
//! ```sh
//! KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete \
//!   cargo run --release --example m4_exponent_audit -- \
//!     --arm semaev --a 1 --n 13 --targets 16 --seed 20260928 --out cell.jsonl
//! ```
//!
//! Stage diagnostic only (AGENTS.md §8): no relation is collected, no
//! logarithm is computed, and nothing here is an end-to-end ECDLP cost.

use crypto_lib::binary_ecc::curve::point_neg;
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_bench::random_control_system;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, f4_profile, f4_profile_reset, solve_boolean_system_filtered,
    solving_degree, split_rule_default, system_degree, FieldStructure, SolveOptions, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_explicit_frobenius_orbit_factor_base, groebner_decompose, point_key, points_with_x,
    span_f2, FactorBaseDomain, FrobeniusFactorBase, KoblitzCurve,
};
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

const M: usize = 4;

/// Verbatim copy of `koblitz_bench::random_subspace_basis` (HEAD), so the
/// subspace draw is the one `dreg_ladder` makes on either source tree.
fn random_subspace_basis(n: u32, ell: usize, rng: &mut StdRng) -> Vec<F2mElement> {
    assert!(
        ell >= 1 && (ell as u32) < n && n < 64,
        "need 1 ≤ ℓ < n < 64"
    );
    let mask = (1u64 << n) - 1;
    let mut echelon: Vec<u64> = Vec::new();
    let mut basis = Vec::with_capacity(ell);
    while basis.len() < ell {
        let v = rng.gen::<u64>() & mask;
        let r = echelon.iter().fold(v, |r, &e| r.min(r ^ e));
        if r == 0 {
            continue;
        }
        echelon.push(r);
        echelon.sort_unstable_by(|a, b| b.cmp(a));
        basis.push(F2mElement::from_biguint(&BigUint::from(v), n));
    }
    basis
}

fn word(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

/// `F = {P : x(P) ∈ span(basis)}` as a [`FrobeniusFactorBase`] whose orbit
/// fields carry negation only (the subspace is not Frobenius-stable).  Built
/// by editing the public fields of a one-orbit base before any consumer has
/// used it, so no library change is needed on either source tree.  The
/// oracle reads only `subspace_basis`, `points` and the caller's index map.
fn subspace_factor_base(kc: &KoblitzCurve, basis: &[F2mElement]) -> Option<FrobeniusFactorBase> {
    let n = kc.n;
    let mut fb = build_explicit_frobenius_orbit_factor_base(kc, &[F2mElement::zero(n)])?;
    let subspace = span_f2(basis, n);
    if subspace.len() != 1usize << basis.len() {
        return None;
    }
    let mut points = Vec::new();
    for x in &subspace {
        points.extend(points_with_x(&kc.curve, x));
    }
    if points.is_empty() {
        return None;
    }
    let index: HashMap<(BigUint, BigUint), usize> = points
        .iter()
        .enumerate()
        .map(|(i, p)| (point_key(p), i))
        .collect();
    let mut signed_orbits: Vec<Vec<usize>> = Vec::new();
    let mut signed_orbit_of = vec![(usize::MAX, 0u32, false); points.len()];
    for i in 0..points.len() {
        if signed_orbit_of[i].0 != usize::MAX {
            continue;
        }
        let o = signed_orbits.len();
        let j = *index.get(&point_key(&point_neg(&points[i])))?;
        signed_orbit_of[i] = (o, 0, false);
        let mut members = vec![i];
        if j != i {
            signed_orbit_of[j] = (o, 0, true);
            members.push(j);
        }
        signed_orbits.push(members);
    }
    fb.domain = FactorBaseDomain::LinearSubspace;
    fb.ell = basis.len() as u32;
    fb.f_j = 0;
    fb.linearised_exponents = Vec::new();
    fb.subspace = subspace;
    fb.subspace_basis = basis.to_vec();
    fb.orbits = (0..points.len()).map(|i| vec![i]).collect();
    fb.orbit_of = (0..points.len()).map(|i| (i, 0)).collect();
    fb.points = points;
    fb.signed_orbits = signed_orbits;
    fb.signed_orbit_of = signed_orbit_of;
    Some(fb)
}

/// Whether `m` summands whose cofactor classes `[r]P` range over `F`'s can
/// sum to `O` — the admissibility test, computed directly.
fn admissible(kc: &KoblitzCurve, fb: &FrobeniusFactorBase, m: usize) -> bool {
    let mut classes: Vec<BinaryPoint> = Vec::new();
    for p in &fb.points {
        let c = kc.mul(p, &kc.subgroup_order);
        if !classes.contains(&c) {
            classes.push(c);
        }
    }
    let mut layer: Vec<BinaryPoint> = vec![BinaryPoint::Infinity];
    for _ in 0..m {
        let mut next: Vec<BinaryPoint> = Vec::new();
        for a in &layer {
            for c in &classes {
                let s = kc.add(a, c);
                if !next.contains(&s) {
                    next.push(s);
                }
            }
        }
        layer = next;
    }
    layer.contains(&BinaryPoint::Infinity)
}

/// The library's exhaustive `decompose`, with every point addition counted.
struct Enumerator<'a> {
    kc: &'a KoblitzCurve,
    neg: Vec<BinaryPoint>,
    index: &'a HashMap<(BigUint, BigUint), usize>,
    adds: u64,
}

impl Enumerator<'_> {
    fn decompose(&mut self, target: &BinaryPoint, m: usize, start: usize) -> Option<Vec<usize>> {
        if m == 0 {
            return (*target == BinaryPoint::Infinity).then(Vec::new);
        }
        if m == 1 {
            let idx = *self.index.get(&point_key(target))?;
            return (idx >= start).then(|| vec![idx]);
        }
        for i in start..self.neg.len() {
            let rest = self.kc.add(target, &self.neg[i]);
            self.adds += 1;
            if let Some(mut tail) = self.decompose(&rest, m - 1, i) {
                let mut out = vec![i];
                out.append(&mut tail);
                return Some(out);
            }
        }
        None
    }
}

/// Exact number of Boolean roots of the chained `m = 4` system:
/// `(x₁,x₂,x₃,x₄) ∈ V⁴`, `e₁, e₂ ∈ F_{2^n}` with `S₃(x₁,x₂,e₁) = 0`,
/// `S₃(e₁,x₃,e₂) = 0`, `S₃(e₂,x₄,x_R) = 0`, where
/// `S₃(a,c,d) = (a+c)²d² + acd + (ac)² + b` (the chain
/// `build_decomposition_system` descends; roots ↔ assignments).
/// Cost `O(2^{2ℓ+n})` field multiplications.
fn chained_s3_m4_solution_count(kc: &KoblitzCurve, basis: &[F2mElement], x_r: u64) -> u64 {
    let n = kc.n;
    let gf = Gf2::new(&kc.curve.irreducible);
    let bb = word(&kc.curve.b);
    let words: Vec<u64> = basis.iter().map(word).collect();
    let span: Vec<u64> = (0..1u64 << words.len())
        .map(|c| {
            (0..words.len())
                .filter(|t| (c >> t) & 1 == 1)
                .fold(0, |acc, t| acc ^ words[t])
        })
        .collect();
    let size = 1usize << n;
    let squares: Vec<u64> = (0..size as u64).map(|u| gf.sqr(u)).collect();
    let s3 = |a: u64, c: u64, d: u64| {
        let ac = gf.mul(a, c);
        gf.mul(gf.sqr(a ^ c), gf.sqr(d)) ^ gf.mul(ac, d) ^ gf.sqr(ac) ^ bb
    };
    // For fixed (a, c), all d with S₃(a, c, d) = 0, by scanning d.
    let roots_in_third = |a: u64, c: u64, out: &mut Vec<u64>| {
        out.clear();
        let (lead, mid) = (gf.sqr(a ^ c), gf.mul(a, c));
        let tail = gf.sqr(mid) ^ bb;
        for (d, &dd) in squares.iter().enumerate() {
            let d = d as u64;
            if gf.mul(lead, dd) ^ gf.mul(mid, d) ^ tail == 0 {
                out.push(d);
            }
        }
    };
    // P1[e1] = #{(x1, x2) ∈ V² : S₃(x1, x2, e1) = 0}.
    let mut p1 = vec![0u64; size];
    let mut buf = Vec::new();
    for &x1 in &span {
        for &x2 in &span {
            roots_in_third(x1, x2, &mut buf);
            for &e1 in &buf {
                p1[e1 as usize] += 1;
            }
        }
    }
    // C3[e2] = #{x4 ∈ V : S₃(e2, x4, x_R) = 0}; S₃ is symmetric, so the
    // e2 with C3 > 0 are the roots in the third slot of S₃(x4, x_R, ·).
    let mut c3: HashMap<u64, u64> = HashMap::new();
    for &x4 in &span {
        roots_in_third(x4, x_r, &mut buf);
        for &e2 in &buf {
            *c3.entry(e2).or_insert(0) += 1;
        }
    }
    let mut total = 0u64;
    let mut e2s: Vec<(u64, u64)> = c3.into_iter().collect();
    e2s.sort_unstable();
    for (e2, c) in e2s {
        for &x3 in &span {
            // e1 with S₃(e1, x3, e2) = 0 = roots in the third slot of S₃(x3, e2, ·).
            roots_in_third(x3, e2, &mut buf);
            for &e1 in &buf {
                debug_assert_eq!(s3(e1, x3, e2), 0);
                total += p1[e1 as usize] * c;
            }
        }
    }
    total
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let flag = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let num = |name: &str, default: u64| {
        flag(name)
            .map(|v| {
                v.parse::<u64>()
                    .unwrap_or_else(|_| panic!("{name} takes an integer"))
            })
            .unwrap_or(default)
    };
    let arm = flag("--arm").expect("--arm semaev|enumerate|null|degree");
    assert!(
        matches!(arm.as_str(), "semaev" | "enumerate" | "null" | "degree"),
        "unknown arm {arm}"
    );
    let a = num("--a", 0) as u8;
    let n = num("--n", 9) as u32;
    let ell = flag("--ell")
        .map(|v| v.parse::<usize>().expect("--ell"))
        .unwrap_or(((n + 2) / 4) as usize); // round(n/4); n odd, so no ties
    let targets = num("--targets", 16) as u32;
    let first = num("--first", 0) as u32;
    let seed = num("--seed", 20_260_928);
    let node_budget = num("--node-budget", 20_000) as usize;
    let d_max = num("--d-max", 5) as u32;
    // Degree arm: measure the refutation degree on at most this many
    // Boolean-unsatisfiable targets of the cell (the first ones in order).
    let max_unsat = num("--max-unsat", 4) as usize;
    let mut unsat_measured = 0usize;
    let label = flag("--label").unwrap_or_else(|| "unlabelled".into());
    let mut out_file = flag("--out").map(|p| {
        assert!(
            !std::path::Path::new(&p).exists(),
            "{p} exists; never overwritten"
        );
        if let Some(dir) = std::path::Path::new(&p).parent() {
            std::fs::create_dir_all(dir).expect("create output directory");
        }
        std::fs::File::create(&p).expect("create --out file")
    });

    let kc = KoblitzCurve::new(a, n).unwrap_or_else(|| panic!("no curve K_{a}/2^{n}"));
    let r = kc.subgroup_order.clone();
    let r64: u64 = r.iter_u64_digits().next().unwrap_or(0);
    assert!(r.bits() <= 63 && r64 > 2, "subgroup order out of range");
    let st = FieldStructure::new(n, &kc.curve.irreducible);
    let engine = SolverEngine::default();
    let policy = serde_json::json!({
        "KIC_CHAIN_ORDER": std::env::var("KIC_CHAIN_ORDER").ok(),
        "KIC_LINEAR_ELIM": std::env::var("KIC_LINEAR_ELIM").ok(),
        "KIC_F4_DROP": std::env::var("KIC_F4_DROP").ok(),
        "KIC_F4_MULTIPLIERS": std::env::var("KIC_F4_MULTIPLIERS").ok(),
        "F4_F2_MAX_ROWS": std::env::var("F4_F2_MAX_ROWS").ok(),
        "F4_F2_MAX_COLS": std::env::var("F4_F2_MAX_COLS").ok(),
    });

    // The cell: dreg_ladder's per-cell seed with the curve folded in; the
    // subspace is the first admissible draw from it.
    let cell_seed = seed ^ (u64::from(a) << 48) ^ (u64::from(n) << 40) ^ ((ell as u64) << 32);
    let mut rng = StdRng::seed_from_u64(cell_seed);
    let mut draws = 0u32;
    let (basis, fb) = loop {
        draws += 1;
        assert!(draws <= 64, "no admissible subspace in 64 draws");
        let basis = random_subspace_basis(n, ell, &mut rng);
        if let Some(fb) = subspace_factor_base(&kc, &basis) {
            if admissible(&kc, &fb, M) {
                break (basis, fb);
            }
        }
    };
    let index_of = fb.index_map();
    let v_basis: Vec<u64> = basis.iter().map(word).collect();
    let neg: Vec<BinaryPoint> = fb.points.iter().map(point_neg).collect();
    let cell = format!("K{a}n{n}l{ell}");
    eprintln!(
        "{cell} arm={arm} label={label}: |F| = {} points, V = {v_basis:?}, draws = {draws}, engine = {:?}",
        fb.points.len(),
        engine.effective()
    );

    let g = kc.generator().clone();
    let mut emit = |line: String| {
        println!("{line}");
        if let Some(f) = out_file.as_mut() {
            writeln!(f, "{line}").expect("write --out");
            f.flush().ok();
        }
        std::io::stdout().flush().ok();
    };

    for t in first..first + targets {
        // Uniform natural target in the prime-order subgroup, its own seed.
        let mut trng = StdRng::seed_from_u64(
            cell_seed
                .wrapping_mul(0x9E37_79B9_7F4A_7C15)
                .wrapping_add(0x7A26_0000 + u64::from(t)),
        );
        let k = 1 + trng.gen::<u64>() % (r64 - 1);
        let target = kc.mul(&g, &BigUint::from(k));
        let x_r = match &target {
            BinaryPoint::Affine { x, .. } => x.clone(),
            BinaryPoint::Infinity => unreachable!("k in [1, r-1]"),
        };
        let sys = build_decomposition_system(&fb.subspace_basis, &x_r, &kc.curve.b, M, &st)
            .expect("system fits");
        let deg = system_degree(&sys.equations);
        let terms =
            sys.equations.iter().map(|e| e.terms.len()).sum::<usize>() / sys.equations.len().max(1);
        let head = format!(
            r#""label":"{label}","arm":"{arm}","cell":"{cell}","a":{a},"n":{n},"ell":{ell},"m":{M},"seed":{seed},"cell_seed":{cell_seed},"subspace_draws":{draws},"v_basis":{v_basis:?},"fb_points":{},"target":{t},"k":{k},"x_r":{},"n_vars":{},"n_eqs":{},"degree":{deg},"terms_per_eq":{terms},"node_budget":{node_budget},"policy":{policy}"#,
            fb.points.len(),
            word(&x_r),
            sys.n_vars,
            sys.equations.len(),
        );
        let started = Instant::now();
        let body = match arm.as_str() {
            "semaev" => {
                // Ground truth first, uncharged.
                let truth = Enumerator {
                    kc: &kc,
                    neg: neg.clone(),
                    index: &index_of,
                    adds: 0,
                }
                .decompose(&target, M, 0)
                .is_some();
                let started = Instant::now();
                f4_profile_reset();
                let (found, stats) =
                    groebner_decompose(&kc, &fb, &index_of, &st, &target, M, engine, node_budget);
                let p = f4_profile();
                let secs = started.elapsed().as_secs_f64();
                let verified = found.as_ref().map(|idxs| {
                    idxs.len() == M
                        && idxs
                            .iter()
                            .fold(BinaryPoint::Infinity, |acc, &i| kc.add(&acc, &fb.points[i]))
                            == target
                });
                let verdict = match (&found, stats.exhausted) {
                    (Some(_), _) => "satisfiable",
                    (None, false) => "refuted",
                    (None, true) => "censored",
                };
                let agree = match verdict {
                    "satisfiable" => Some(truth),
                    "refuted" => Some(!truth),
                    _ => None,
                };
                format!(
                    r#""verdict":"{verdict}","truth_decomposable":{truth},"agree":{},"verified":{},"word_ops":{},"f4_calls":{},"f4_oversize":{},"matrix_rows":{},"matrix_cols":{},"reductions":{},"infeasible_branches":{},"splits":{},"oversize":{},"eliminated":{},"max_degree_built":{},"exhausted":{},"secs":{secs:.4}"#,
                    agree.map_or("null".into(), |b| b.to_string()),
                    verified.map_or("null".into(), |b| b.to_string()),
                    p.word_ops,
                    p.calls,
                    p.oversize,
                    p.rows,
                    p.cols,
                    stats.reductions,
                    stats.infeasible_branches,
                    stats.splits,
                    stats.oversize,
                    stats.eliminated,
                    stats.max_degree_built,
                    stats.exhausted,
                )
            }
            "enumerate" => {
                let mut en = Enumerator {
                    kc: &kc,
                    neg: neg.clone(),
                    index: &index_of,
                    adds: 0,
                };
                let found = en.decompose(&target, M, 0);
                let secs = started.elapsed().as_secs_f64();
                format!(
                    r#""verdict":"{}","group_adds":{},"secs":{secs:.4}"#,
                    if found.is_some() {
                        "satisfiable"
                    } else {
                        "refuted"
                    },
                    en.adds
                )
            }
            "null" => {
                let null_seed = cell_seed ^ 0x0C01_7201_u64.wrapping_mul(u64::from(t) + 1);
                let polys =
                    random_control_system(sys.n_vars, sys.equations.len(), deg, terms, null_seed);
                let opts = SolveOptions {
                    engine,
                    max_solutions: usize::MAX,
                    node_budget,
                    split_rule: split_rule_default(),
                };
                let mut roots = 0u64;
                f4_profile_reset();
                let (_, stats) = solve_boolean_system_filtered(&polys, sys.n_vars, &opts, |_| {
                    roots += 1;
                    false
                });
                let p = f4_profile();
                let secs = started.elapsed().as_secs_f64();
                format!(
                    r#""verdict":"{}","null_seed":{null_seed},"roots_seen":{roots},"word_ops":{},"f4_calls":{},"f4_oversize":{},"reductions":{},"infeasible_branches":{},"splits":{},"oversize":{},"exhausted":{},"secs":{secs:.4}"#,
                    if stats.exhausted {
                        "censored"
                    } else {
                        "complete"
                    },
                    p.word_ops,
                    p.calls,
                    p.oversize,
                    stats.reductions,
                    stats.infeasible_branches,
                    stats.splits,
                    stats.oversize,
                    stats.exhausted,
                )
            }
            "degree" => {
                let solutions = chained_s3_m4_solution_count(&kc, &fb.subspace_basis, word(&x_r));
                // Instrument check: the engine's own full-tree root count on
                // the same system must equal the exact count when complete.
                let opts = SolveOptions {
                    engine,
                    max_solutions: usize::MAX,
                    node_budget,
                    split_rule: split_rule_default(),
                };
                let mut solver_roots = 0u64;
                let (_, cstats) =
                    solve_boolean_system_filtered(&sys.equations, sys.n_vars, &opts, |_| {
                        solver_roots += 1;
                        false
                    });
                let count_check = if cstats.exhausted {
                    "null".to_string()
                } else {
                    (solver_roots == solutions).to_string()
                };
                let outcome = if solutions > 0 {
                    r#"{"kind":"satisfiable"}"#.to_string()
                } else if unsat_measured >= max_unsat {
                    r#"{"kind":"not_measured"}"#.to_string()
                } else {
                    unsat_measured += 1;
                    let (d, profs) = solving_degree(&sys.equations, sys.n_vars, d_max);
                    let built = profs.last().map(|p| p.degree);
                    match d {
                        Some(degree) => format!(
                            r#"{{"kind":"resolved","degree":{degree},"refuted":{}}}"#,
                            profs.last().map(|p| p.refuted).unwrap_or(false)
                        ),
                        None if built == Some(d_max) => {
                            format!(r#"{{"kind":"at_least","degree":{}}}"#, d_max + 1)
                        }
                        None => format!(
                            r#"{{"kind":"caps_hit","built":{}}}"#,
                            built.map_or("null".into(), |b| b.to_string())
                        ),
                    }
                };
                format!(
                    r#""boolean_solutions":{solutions},"solver_roots":{solver_roots},"solver_exhausted":{},"count_check":{count_check},"d_max":{d_max},"outcome":{outcome},"secs":{:.4}"#,
                    cstats.exhausted,
                    started.elapsed().as_secs_f64()
                )
            }
            _ => unreachable!(),
        };
        emit(format!("{{{head},{body}}}"));
    }
}
