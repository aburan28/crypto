//! Complete bounded generic Boolean solves with matched search semantics.
#[allow(dead_code)]
mod algebra {
    include!("kernel.rs");
}
use std::hint::black_box;
use std::time::Instant;

type System = Vec<Vec<u64>>;
include!("quadratic.rs");
include!("packed.rs");
include!("basis.rs");
include!("affine.rs");
include!("compact.rs");
include!("linearization.rs");
include!("syndrome.rs");
include!("leaf.rs");
#[derive(Clone, Debug, PartialEq, Eq)]
enum Outcome {
    Sat(u64),
    Unsat,
    Unknown(&'static str),
}
#[derive(Default, Clone, Debug, PartialEq, Eq)]
struct Logical {
    nodes: u64,
    kernel_calls: u64,
    decisions: u64,
    forced: u64,
    affine_eliminated: u64,
    derived_rows: u64,
    enumeration_points: u64,
    enumeration_batches: u64,
    enumeration_leaves: u64,
    specialized_terms: u64,
    source_rows: u64,
    source_columns: u64,
    max_depth: usize,
}
#[derive(Default)]
struct Profile {
    kernel_ns: u128,
    substitution_ns: u128,
    enumeration_ns: u128,
    kernel_calls_by_active: Vec<u64>,
    flat_calls: u64,
}
struct Solver {
    n: u8,
    limit: u64,
    session: Option<algebra::Session>,
    logical: Logical,
    profile: Profile,
    trace: u64,
    small_degree_two: bool,
    merge_specialization: bool,
}
impl Solver {
    fn event(&mut self, value: u64) {
        self.trace = (self.trace ^ value).wrapping_mul(0x100000001b3);
    }
    fn specialize(&mut self, system: System, mask: u64, values: u64) -> System {
        let zero = mask & !values;
        let mut out = Vec::new();
        for mut poly in system {
            self.logical.specialized_terms += poly.len() as u64;
            let terms = if self.merge_specialization {
                poly.retain(|m| m & zero == 0);
                let ones = mask & values;
                if ones.count_ones() == 1 {
                    specialize_one_in_place(&mut poly, ones);
                } else if ones != 0 {
                    for m in &mut poly {
                        *m &= !ones;
                    }
                    algebra::canonical(&mut poly);
                }
                poly
            } else {
                let mut terms: Vec<_> = poly
                    .into_iter()
                    .filter(|m| m & zero == 0)
                    .map(|m| m & !mask)
                    .collect();
                algebra::canonical(&mut terms);
                terms
            };
            if !terms.is_empty() && !out.contains(&terms) {
                out.push(terms);
            }
        }
        out
    }
    fn assign_from_affine(&mut self, rows: &[u64]) -> Result<(u64, u64), ()> {
        let constant = 1u64 << self.n;
        let mut mask = 0;
        let mut values = 0;
        for &row in rows {
            if row == constant {
                return Err(());
            }
            let variables = row & (constant - 1);
            if variables.count_ones() == 1 {
                let value = if row & constant != 0 { variables } else { 0 };
                if mask & variables != 0 && (values & variables) != value {
                    return Err(());
                }
                mask |= variables;
                values |= value;
            }
        }
        Ok((mask, values))
    }
    fn visit(
        &mut self,
        mut system: System,
        mut known: u64,
        mut values: u64,
        depth: usize,
    ) -> Outcome {
        if self.logical.nodes == self.limit {
            self.event(0xcaf);
            return Outcome::Unknown("NODE_CAP");
        }
        self.logical.nodes += 1;
        self.logical.max_depth = self.logical.max_depth.max(depth);
        self.event(0x100);
        self.event(known);
        self.event(values);
        self.event(system.len() as u64);
        for poly in &system {
            self.event(poly.len() as u64);
            for &m in poly {
                self.event(m);
            }
        }
        loop {
            if system.is_empty() {
                self.event(0x501);
                self.event(values);
                return Outcome::Sat(values);
            }
            if system.iter().any(|p| p.len() == 1 && p[0] == 0) {
                self.event(0x500);
                return Outcome::Unsat;
            }
            let affine: Vec<_> = system
                .iter()
                .filter(|p| p.iter().all(|m| m.count_ones() <= 1))
                .map(|p| {
                    p.iter()
                        .fold(0, |v, &m| v ^ if m == 0 { 1u64 << self.n } else { m })
                })
                .collect();
            let all_affine = affine.len() == system.len();
            let reduced = algebra::affine_rref(self.n, affine.into_iter());
            let (mask, ones) = match self.assign_from_affine(&reduced) {
                Ok(x) => x,
                Err(()) => {
                    self.event(0x502);
                    return Outcome::Unsat;
                }
            };
            if all_affine {
                // Free variables are zero. In RREF every pivot depends only on
                // free variables and the constant, so read its constant value.
                let constant = 1u64 << self.n;
                for row in reduced {
                    let variable = row & (constant - 1);
                    if variable != 0 {
                        let pivot = variable & variable.wrapping_neg();
                        if row & constant != 0 {
                            values |= pivot;
                        }
                    }
                }
                self.event(0x503);
                self.event(values);
                return Outcome::Sat(values);
            }
            if mask != 0 {
                self.logical.forced += mask.count_ones() as u64;
                self.event(0x200);
                self.event(mask);
                self.event(ones);
                known |= mask;
                values |= ones;
                system = self.specialize(system, mask, ones);
                continue;
            }
            if self.session.is_some() {
                let active = system.iter().flatten().fold(0, |m, &term| m | term);
                let width = active.count_ones() as usize;
                if self.small_degree_two && width > 10 {
                    break;
                }
                self.logical.kernel_calls += 1;
                self.profile.kernel_calls_by_active[width] += 1;
                let tick = Instant::now();
                let (original, reply) = if self.small_degree_two {
                    self.session.as_mut().unwrap().run_degree(self.n, 2, system)
                } else {
                    self.session.as_mut().unwrap().run(self.n, system)
                };
                self.profile.kernel_ns += tick.elapsed().as_nanos();
                system = original;
                let reply = match reply {
                    Ok(r) => r,
                    Err(reason) => {
                        self.event(0x600);
                        return Outcome::Unknown(reason);
                    }
                };
                self.profile.flat_calls += u64::from(reply.used_flat);
                self.logical.source_rows += reply.source_rows as u64;
                self.logical.source_columns += reply.source_columns as u64;
                self.event(0x300);
                self.event(reply.tail.len() as u64);
                for &row in &reply.tail {
                    self.event(row);
                }
                let (mask, ones) = match self.assign_from_affine(&reply.tail) {
                    Ok(x) => x,
                    Err(()) => {
                        self.event(0x504);
                        return Outcome::Unsat;
                    }
                };
                if mask != 0 {
                    self.logical.forced += mask.count_ones() as u64;
                    self.event(0x201);
                    self.event(mask);
                    self.event(ones);
                    known |= mask;
                    values |= ones;
                    system = self.specialize(system, mask, ones);
                    continue;
                }
            }
            break;
        }
        // Most frequent occurring variable; ties go to the lowest index.
        let mut counts = vec![0usize; self.n as usize];
        for poly in &system {
            for &term in poly {
                let mut m = term;
                while m != 0 {
                    counts[m.trailing_zeros() as usize] += 1;
                    m &= m - 1;
                }
            }
        }
        let variable = (0..self.n as usize)
            .filter(|&i| counts[i] > 0)
            .max_by_key(|&i| (counts[i], std::cmp::Reverse(i)))
            .unwrap();
        let bit = 1u64 << variable;
        self.logical.decisions += 1;
        self.event(0x400);
        self.event(variable as u64);
        let left = self.specialize(system.clone(), bit, 0);
        match self.visit(left, known | bit, values, depth + 1) {
            Outcome::Unsat => {
                let right = self.specialize(system, bit, bit);
                self.visit(right, known | bit, values | bit, depth + 1)
            }
            result => result,
        }
    }
}

fn specialize_one_in_place(poly: &mut Vec<u64>, variable: u64) {
    let mut removed = Vec::new();
    let mut keep = 0;
    for i in 0..poly.len() {
        let m = poly[i];
        if m & variable != 0 {
            removed.push(m & !variable);
        } else {
            poly[keep] = m;
            keep += 1;
        }
    }
    poly.truncate(keep);
    if removed.is_empty() {
        return;
    }
    // Both sublists retain monomial order: removing one common factor preserves
    // relative total degrees and the symmetric differences defining DegRevLex.
    let total = keep + removed.len();
    poly.resize(total, 0);
    let (mut i, mut j, mut write) = (keep, removed.len(), total);
    while i > 0 && j > 0 {
        match algebra::mono_order(poly[i - 1], removed[j - 1]) {
            std::cmp::Ordering::Greater => {
                write -= 1;
                i -= 1;
                poly[write] = poly[i];
            }
            std::cmp::Ordering::Less => {
                write -= 1;
                j -= 1;
                poly[write] = removed[j];
            }
            std::cmp::Ordering::Equal => {
                i -= 1;
                j -= 1;
            }
        }
    }
    while i > 0 {
        write -= 1;
        i -= 1;
        poly[write] = poly[i];
    }
    while j > 0 {
        write -= 1;
        j -= 1;
        poly[write] = removed[j];
    }
    poly.copy_within(write..total, 0);
    poly.truncate(total - write);
}
struct Solved {
    outcome: Outcome,
    logical: Logical,
    profile: Profile,
    trace: u64,
}
fn solve(system: &System, n: u8, arm: &str, limit: u64) -> Solved {
    if arm.starts_with("packed_gray") {
        return solve_leaf(
            system,
            n,
            limit,
            if arm.contains("12") { 12 } else { 8 },
            arm.ends_with("simd"),
        );
    }
    if arm == "gray_scalar" || arm == "gray_simd" {
        return solve_gray(system, n, limit, arm == "gray_simd");
    }
    if [
        "affine_basis_list",
        "affine_basis_fast",
        "affine_sl_basis_list",
        "affine_sl_basis_fast",
    ]
    .contains(&arm)
    {
        return solve_affine_basis_only(
            system,
            n,
            limit,
            arm.contains("_sl_"),
            arm.ends_with("_fast"),
        );
    }
    if arm == "affine_fast" {
        return solve_affine_fast(system, n, limit);
    }
    if arm == "affine_sl_fast" {
        return solve_affine_linearized_fast(system, n, limit);
    }
    if arm == "affine_sl_list" || arm == "affine_sl_compact" {
        return solve_affine_linearized(system, n, limit, arm == "affine_sl_compact");
    }
    if arm == "affine_compact" {
        return solve_affine_compact(system, n, limit);
    }
    if arm == "affine_list" || arm == "affine_wide" {
        return solve_affine(system, n, limit, arm == "affine_wide");
    }
    if arm == "tail_list" || arm == "tail_wide" {
        return solve_tail_basis(system, n, limit, arm == "tail_wide");
    }
    if arm == "basis_list" || arm == "basis_wide" {
        return solve_basis(system, n, limit, arm == "basis_wide");
    }
    if arm == "packed_state" {
        return solve_packed(system, n, limit);
    }
    if arm == "quadratic_state" {
        return solve_quadratic(system, n, limit);
    }
    let session = match arm {
        "search" | "merge_search" => None,
        "flat" => Some(algebra::Session::new(algebra::Method::Flat)),
        "bucket" => Some(algebra::Session::new(algebra::Method::Bucket)),
        "hybrid" => Some(algebra::Session::new(algebra::Method::Hybrid)),
        "small_flat" => Some(algebra::Session::new(algebra::Method::Flat)),
        "word_tail" => Some(algebra::Session::new(algebra::Method::Word2)),
        _ => panic!("arm"),
    };
    let mut solver = Solver {
        n,
        limit,
        session,
        logical: Logical::default(),
        profile: Profile {
            kernel_calls_by_active: vec![0; n as usize + 1],
            ..Profile::default()
        },
        trace: 0xcbf29ce484222325,
        small_degree_two: matches!(arm, "small_flat" | "word_tail"),
        merge_specialization: arm == "merge_search",
    };
    let outcome = solver.visit(system.clone(), 0, 0, 0);
    Solved {
        outcome,
        logical: solver.logical,
        profile: solver.profile,
        trace: solver.trace,
    }
}
fn satisfies(system: &System, assignment: u64) -> bool {
    system
        .iter()
        .all(|p| p.iter().fold(false, |v, &m| v ^ (assignment & m == m)) == false)
}
#[cfg(test)]
fn brute(system: &System, n: u8) -> Outcome {
    match (0..(1u64 << n)).find(|&a| satisfies(system, a)) {
        Some(a) => Outcome::Sat(a),
        None => Outcome::Unsat,
    }
}
fn next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}
fn fixture(n: u8, seed: u64, family: &str) -> (System, Option<u64>) {
    assert!(["planted", "cross_planted", "unplanted"].contains(&family));
    let mut state = seed;
    let witness = next(&mut state) & ((1u64 << n) - 1);
    let mut system = Vec::new();
    for _ in 0..n + 2 {
        let mut terms = std::collections::BTreeSet::new();
        while terms.len() < 12 {
            let a = next(&mut state) % u64::from(n);
            let b = next(&mut state) % u64::from(n);
            if a != b {
                terms.insert((1u64 << a) | (1u64 << b));
            }
        }
        for _ in 0..4 {
            terms.insert(1u64 << (next(&mut state) % u64::from(n)));
        }
        let mut poly: Vec<_> = terms.into_iter().collect();
        let constant = if family == "unplanted" {
            next(&mut state) & 1 != 0
        } else {
            poly.iter().fold(false, |v, &m| v ^ (witness & m == m))
        };
        if constant {
            poly.push(0);
        }
        algebra::canonical(&mut poly);
        system.push(poly);
    }
    if family == "cross_planted" {
        let bit = 1u64 << (seed % u64::from(n));
        let mut second = system[0].clone();
        second.push(bit);
        if witness & bit != 0 {
            second.push(0);
        }
        algebra::canonical(&mut second);
        system[1] = second;
    }
    let planted = if family == "unplanted" {
        None
    } else {
        Some(witness)
    };
    if let Some(w) = planted {
        assert!(satisfies(&system, w));
    }
    (system, planted)
}
fn status(outcome: &Outcome) -> &'static str {
    match outcome {
        Outcome::Sat(_) => "SAT",
        Outcome::Unsat => "UNSAT",
        Outcome::Unknown(_) => "UNKNOWN",
    }
}
const ALL_ARMS: [&str; 30] = [
    "search",
    "flat",
    "bucket",
    "hybrid",
    "small_flat",
    "word_tail",
    "merge_search",
    "quadratic_state",
    "packed_state",
    "basis_list",
    "basis_wide",
    "tail_list",
    "tail_wide",
    "affine_list",
    "affine_wide",
    "affine_compact",
    "affine_sl_list",
    "affine_sl_compact",
    "affine_fast",
    "affine_sl_fast",
    "affine_basis_list",
    "affine_basis_fast",
    "affine_sl_basis_list",
    "affine_sl_basis_fast",
    "gray_scalar",
    "gray_simd",
    "packed_gray8_scalar",
    "packed_gray8_simd",
    "packed_gray12_scalar",
    "packed_gray12_simd",
];
fn benchmark_arms(mode: Option<&str>) -> &'static [&'static str] {
    let count = match mode {
        None => 4,
        Some("with-word") => 6,
        Some("with-merge") => 7,
        Some("with-quadratic") => 8,
        Some("with-packed") => 9,
        Some("with-basis") => 11,
        Some("with-tail") => 13,
        Some("with-affine") => 15,
        Some("with-compact") => 16,
        Some("with-products") => 18,
        Some("with-fast") => 20,
        Some("with-basis-only") => 24,
        Some("with-gray") => 26,
        Some("with-leaf") => 30,
        _ => panic!("unknown benchmark mode"),
    };
    &ALL_ARMS[..count]
}
fn policy_group(arm: &str) -> usize {
    match arm {
        "flat" | "bucket" | "hybrid" => 0,
        "small_flat" | "word_tail" => 1,
        "search" | "merge_search" | "quadratic_state" | "packed_state" => 2,
        "basis_list" | "basis_wide" => 3,
        "tail_list" | "tail_wide" => 4,
        "affine_list" | "affine_wide" | "affine_compact" | "affine_fast" => 5,
        "affine_sl_list" | "affine_sl_compact" | "affine_sl_fast" => 6,
        "affine_basis_list" | "affine_basis_fast" => 7,
        "affine_sl_basis_list" | "affine_sl_basis_fast" => 8,
        "gray_scalar" | "gray_simd" => 9,
        "packed_gray8_scalar" | "packed_gray8_simd" => 10,
        "packed_gray12_scalar" | "packed_gray12_simd" => 11,
        _ => panic!("unknown inference policy"),
    }
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert!(
        [6, 7].contains(&args.len()),
        "worker N SEED FAMILY REPETITIONS NODE_CAP [with-word|with-merge|with-quadratic|with-packed|with-basis|with-tail|with-affine|with-compact|with-products]"
    );
    let mode = args.get(6).map(String::as_str);
    let n = args[1].parse::<u8>().unwrap();
    let seed = args[2].parse::<u64>().unwrap();
    let family = &args[3];
    let reps = args[4].parse::<usize>().unwrap();
    let limit = args[5].parse::<u64>().unwrap();
    assert!([12, 16, 20, 24].contains(&n) && (1..=30).contains(&reps) && limit <= 200_000);
    let (system, planted) = fixture(n, seed, family);
    let reference = solve(&system, n, "search", limit);
    let planted_json = planted.map_or("null".into(), |x| x.to_string());
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"family\":\"{family}\",\"polys\":{:?},\"planted_witness\":{planted_json},\"search_reference\":\"{}\"}}",system,status(&reference.outcome));
    let names = benchmark_arms(mode);
    let with_merge = names.contains(&"merge_search");
    if names.contains(&"packed_state") {
        assert!(system
            .iter()
            .all(|p| p.iter().filter(|m| m.count_ones() == 2).count() + (n as usize) < 64));
    }
    let mut expected_kernel: [Option<(Outcome, Logical, u64)>; 12] = [
        None, None, None, None, None, None, None, None, None, None, None, None,
    ];
    for rep in 0..reps {
        for order in 0..names.len() {
            let arm = names[(rep + order) % names.len()];
            let start = Instant::now();
            let got = solve(black_box(&system), n, arm, limit);
            let solve_ns = start.elapsed().as_nanos();
            let tick = Instant::now();
            match &got.outcome {
                Outcome::Sat(model) => {
                    assert!(satisfies(&system, *model));
                    assert!(*model < (1u64 << n));
                }
                Outcome::Unsat => {
                    assert!(planted.is_none());
                    if !matches!(reference.outcome, Outcome::Unknown(_)) {
                        assert_eq!(reference.outcome, Outcome::Unsat);
                    }
                }
                Outcome::Unknown(_) => (),
            }
            if !matches!(got.outcome, Outcome::Unknown(_))
                && !matches!(reference.outcome, Outcome::Unknown(_))
            {
                assert_eq!(status(&got.outcome), status(&reference.outcome));
            }
            if arm != "search" || with_merge {
                let signature = (got.outcome.clone(), got.logical.clone(), got.trace);
                let group = policy_group(arm);
                if let Some(expected) = &expected_kernel[group] {
                    assert_eq!(&signature, expected);
                } else {
                    expected_kernel[group] = Some(signature);
                }
            }
            let verified = matches!(got.outcome, Outcome::Sat(_))
                || (got.outcome == Outcome::Unsat && reference.outcome == Outcome::Unsat);
            let validation_ns = tick.elapsed().as_nanos();
            let total_ns = start.elapsed().as_nanos();
            let model = match got.outcome {
                Outcome::Sat(m) => m.to_string(),
                _ => "null".into(),
            };
            let reason = match got.outcome {
                Outcome::Unknown(r) => format!("\"{r}\""),
                _ => "null".into(),
            };
            let outcome = status(&got.outcome);
            let l = got.logical;
            let profile = got.profile;
            let substitution_ns = if arm.starts_with("affine_") {
                profile.substitution_ns.to_string()
            } else {
                "null".into()
            };
            let enumeration_ns = if arm.starts_with("gray_") || arm.starts_with("packed_gray") {
                profile.enumeration_ns.to_string()
            } else {
                "null".into()
            };
            let trace = got.trace;
            println!("{{\"type\":\"sample\",\"rep\":{rep},\"order\":{order},\"variant\":\"{arm}\",\"outcome\":\"{outcome}\",\"model\":{model},\"reason\":{reason},\"verified\":{verified},\"solve_ns\":{solve_ns},\"validation_ns\":{validation_ns},\"total_ns\":{total_ns},\"substitution_ns\":{substitution_ns},\"affine_eliminated\":{},\"derived_rows\":{},\"enumeration_points\":{},\"enumeration_batches\":{},\"enumeration_leaves\":{},\"enumeration_ns\":{enumeration_ns},\"kernel_ns\":{},\"nodes\":{},\"kernel_calls\":{},\"decisions\":{},\"forced\":{},\"specialized_terms\":{},\"source_rows\":{},\"source_columns\":{},\"max_depth\":{},\"trace\":{trace},\"flat_calls\":{},\"kernel_calls_by_active\":{:?}}}",l.affine_eliminated,l.derived_rows,l.enumeration_points,l.enumeration_batches,l.enumeration_leaves,profile.kernel_ns,l.nodes,l.kernel_calls,l.decisions,l.forced,l.specialized_terms,l.source_rows,l.source_columns,l.max_depth,profile.flat_calls,profile.kernel_calls_by_active);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exhaustive_small_systems_agree_with_boolean_enumeration() {
        for mask in 0..256u16 {
            let mut poly: Vec<_> = (0u64..8).filter(|m| mask & (1 << m) != 0).collect();
            algebra::canonical(&mut poly);
            for extra in [vec![], vec![0], vec![1], vec![1, 0], vec![2, 4]] {
                let mut second = extra;
                algebra::canonical(&mut second);
                let system = vec![poly.clone(), second];
                let expected = brute(&system, 3);
                for arm in [
                    "search",
                    "flat",
                    "bucket",
                    "hybrid",
                    "small_flat",
                    "word_tail",
                    "merge_search",
                    "quadratic_state",
                    "packed_state",
                    "basis_list",
                    "basis_wide",
                    "tail_list",
                    "tail_wide",
                    "affine_list",
                    "affine_wide",
                    "affine_compact",
                    "affine_sl_list",
                    "affine_sl_compact",
                ] {
                    let got = solve(&system, 3, arm, 1000);
                    assert_eq!(status(&got.outcome), status(&expected));
                    if let Outcome::Sat(model) = got.outcome {
                        assert!(satisfies(&system, model));
                    }
                }
            }
        }
    }
    #[test]
    fn cap_is_unknown_not_unsat() {
        let (system, _) = fixture(12, 17, "planted");
        for arm in [
            "search",
            "flat",
            "bucket",
            "hybrid",
            "small_flat",
            "word_tail",
            "merge_search",
            "quadratic_state",
            "packed_state",
            "basis_list",
            "basis_wide",
            "tail_list",
            "tail_wide",
            "affine_list",
            "affine_wide",
            "affine_compact",
            "affine_sl_list",
            "affine_sl_compact",
        ] {
            assert_eq!(
                solve(&system, 12, arm, 0).outcome,
                Outcome::Unknown("NODE_CAP")
            );
        }
    }
    #[test]
    fn kernel_variants_have_identical_traces_and_verified_models() {
        for family in ["planted", "cross_planted", "unplanted"] {
            let (system, _) = fixture(12, 31337, family);
            let a = solve(&system, 12, "flat", 20000);
            for arm in ["bucket", "hybrid"] {
                let b = solve(&system, 12, arm, 20000);
                assert_eq!(a.outcome, b.outcome);
                assert_eq!(a.logical, b.logical);
                assert_eq!(a.trace, b.trace);
            }
            assert!(!matches!(a.outcome, Outcome::Unknown(_)));
            assert_eq!(status(&a.outcome), status(&brute(&system, 12)));
        }
    }
    #[test]
    fn specialization_and_affine_completion_preserve_solutions() {
        let mut p = vec![0, 1, 2, 3];
        algebra::canonical(&mut p);
        let system = vec![p];
        let mut solver = Solver {
            n: 2,
            limit: 10,
            session: None,
            logical: Logical::default(),
            profile: Profile::default(),
            trace: 0,
            small_degree_two: false,
            merge_specialization: false,
        };
        for bit in [1u64, 2] {
            for value in [0, bit] {
                let reduced = solver.specialize(system.clone(), bit, value);
                for a in 0..4 {
                    if a & bit == value {
                        assert_eq!(satisfies(&system, a), satisfies(&reduced, a & !bit));
                    }
                }
            }
        }
        let mut a = vec![0, 1, 2];
        algebra::canonical(&mut a);
        let system = vec![a, vec![2, 4]];
        for arm in [
            "search",
            "flat",
            "bucket",
            "hybrid",
            "small_flat",
            "word_tail",
            "merge_search",
            "quadratic_state",
            "packed_state",
            "basis_list",
            "basis_wide",
            "tail_list",
            "tail_wide",
            "affine_list",
            "affine_wide",
            "affine_compact",
            "affine_sl_list",
            "affine_sl_compact",
        ] {
            let answer = solve(&system, 3, arm, 100);
            if let Outcome::Sat(model) = answer.outcome {
                assert!(satisfies(&system, model));
            } else {
                panic!("linear system should be satisfiable");
            }
        }
    }
    #[test]
    fn word_degree_two_matches_flat_on_embedded_active_variables() {
        for n in [3, 12, 20, 24] {
            let active: Vec<_> = (0..n.min(8)).map(|i| i * n / n.min(8)).collect();
            let mut state = 31337;
            for _ in 0..32 {
                let mut system = Vec::new();
                for p in 0..12 {
                    let mut terms = Vec::new();
                    for _ in 0..12 {
                        let a = active[next(&mut state) as usize % active.len()];
                        let b = active[next(&mut state) as usize % active.len()];
                        terms.push(if p % 3 == 0 {
                            1u64 << a
                        } else {
                            (1u64 << a) | (1u64 << b)
                        });
                    }
                    if next(&mut state) & 1 != 0 {
                        terms.push(0);
                    }
                    algebra::canonical(&mut terms);
                    system.push(terms);
                }
                let (_, a) =
                    algebra::Session::new(algebra::Method::Flat).run_degree(n, 2, system.clone());
                let (_, b) = algebra::Session::new(algebra::Method::Word2).run_degree(n, 2, system);
                let (a, b) = (a.unwrap(), b.unwrap());
                assert_eq!(a.tail, b.tail);
                assert_eq!(a.source_rows, b.source_rows);
                assert_eq!(a.source_columns, b.source_columns);
            }
        }
        let (_, answer) =
            algebra::Session::new(algebra::Method::Word2).run_degree(12, 2, vec![vec![0]]);
        assert_eq!(answer.unwrap().tail, vec![1 << 12]);
    }
    #[test]
    fn selective_policies_match_each_other_on_complete_search() {
        for family in ["planted", "cross_planted", "unplanted"] {
            let (system, _) = fixture(16, 31337, family);
            let a = solve(&system, 16, "small_flat", 200000);
            let b = solve(&system, 16, "word_tail", 200000);
            assert_eq!(a.outcome, b.outcome);
            assert_eq!(a.logical, b.logical);
            assert_eq!(a.trace, b.trace);
            assert!(a
                .profile
                .kernel_calls_by_active
                .iter()
                .enumerate()
                .all(|(n, &count)| n <= 10 || count == 0));
            assert_eq!(status(&a.outcome), status(&brute(&system, 16)));
        }
    }
    #[test]
    fn ordered_specialization_matches_sorting_for_every_small_assignment() {
        let make = |fast| Solver {
            n: 3,
            limit: 1000,
            session: None,
            logical: Logical::default(),
            profile: Profile::default(),
            trace: 0,
            small_degree_two: false,
            merge_specialization: fast,
        };
        for selected in 0..256u16 {
            let mut poly: Vec<_> = (0u64..8).filter(|m| selected & (1 << m) != 0).collect();
            algebra::canonical(&mut poly);
            let system = vec![poly.clone(), poly, vec![1, 0]];
            for mask in 0u64..8 {
                for values in 0u64..8 {
                    if values & !mask != 0 {
                        continue;
                    }
                    let mut slow = make(false);
                    let mut fast = make(true);
                    let a = slow.specialize(system.clone(), mask, values);
                    let b = fast.specialize(system.clone(), mask, values);
                    assert_eq!(a, b);
                    assert_eq!(
                        slow.logical.specialized_terms,
                        fast.logical.specialized_terms
                    );
                    for point in 0u64..8 {
                        if point & mask == values {
                            assert_eq!(satisfies(&system, point), satisfies(&b, point & !mask));
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn merge_search_keeps_the_complete_search_trace() {
        for family in ["planted", "cross_planted", "unplanted"] {
            let (system, _) = fixture(16, 31337, family);
            let a = solve(&system, 16, "search", 200000);
            let b = solve(&system, 16, "merge_search", 200000);
            assert_eq!(a.outcome, b.outcome);
            assert_eq!(a.logical, b.logical);
            assert_eq!(a.trace, b.trace);
        }
    }
    #[test]
    fn quadratic_state_matches_every_small_partial_assignment() {
        let terms: Vec<_> = (0u64..8).filter(|m| m.count_ones() <= 2).collect();
        for selected in 0..(1usize << terms.len()) {
            let mut poly: Vec<_> = terms
                .iter()
                .enumerate()
                .filter(|(i, _)| selected & (1 << i) != 0)
                .map(|(_, &m)| m)
                .collect();
            algebra::canonical(&mut poly);
            let system = vec![poly.clone(), poly, vec![], vec![1, 0]];
            let (model, state) = QuadraticModel::compile(&system, 3).unwrap();
            assert_eq!(model.materialize(&state), system);
            for mask in 0u64..8 {
                for values in 0u64..8 {
                    if values & !mask != 0 {
                        continue;
                    }
                    let mut reference = Solver {
                        n: 3,
                        limit: 1000,
                        session: None,
                        logical: Logical::default(),
                        profile: Profile::default(),
                        trace: 0,
                        small_degree_two: false,
                        merge_specialization: true,
                    };
                    let expected = reference.specialize(system.clone(), mask, values);
                    let mut logical = Logical::default();
                    let got = model.specialize(state.clone(), mask, values, &mut logical);
                    assert_eq!(model.materialize(&got), expected);
                    assert_eq!(
                        logical.specialized_terms,
                        reference.logical.specialized_terms
                    );
                    for next_mask in 0u64..8 {
                        let next_values = next_mask & !values;
                        let expected_next =
                            reference.specialize(expected.clone(), next_mask, next_values);
                        let actual_next = model.specialize(
                            got.clone(),
                            next_mask,
                            next_values,
                            &mut Logical::default(),
                        );
                        assert_eq!(model.materialize(&actual_next), expected_next);
                    }
                }
            }
        }
    }
    #[test]
    fn quadratic_trace_matches_search_after_deduplication_and_backtracking() {
        for n in [12, 16] {
            for family in ["planted", "cross_planted", "unplanted"] {
                let (system, _) = fixture(n, 31337, family);
                let a = solve(&system, n, "merge_search", 200000);
                let b = solve(&system, n, "quadratic_state", 200000);
                assert_eq!(a.outcome, b.outcome);
                assert_eq!(a.logical, b.logical);
                assert_eq!(a.trace, b.trace);
            }
        }
        let system = vec![vec![7, 1, 0]];
        assert!(QuadraticModel::compile(&system, 3).is_none());
        let a = solve(&system, 3, "merge_search", 1000);
        let b = solve(&system, 3, "quadratic_state", 1000);
        assert_eq!(a.outcome, b.outcome);
        assert_eq!(a.trace, b.trace);
    }
}

#[cfg(test)]
mod packed_tests {
    use super::*;
    #[test]
    fn stack_affine_matches_every_small_row_space() {
        for subset in 0u64..(1u64 << 16) {
            let input: Vec<_> = (0..16).filter(|r| subset & (1u64 << r) != 0).collect();
            let expected = algebra::affine_rref(3, input.iter().copied());
            for reversed in [false, true] {
                let mut got = StackAffine::new();
                let mut ordered = input.clone();
                if reversed {
                    ordered.reverse();
                }
                for row in ordered {
                    got.insert(row);
                }
                got.reduce();
                let actual: Vec<_> = got.rows.iter().copied().filter(|&r| r != 0).collect();
                assert_eq!(actual, expected);
            }
        }
        let mut got = StackAffine::new();
        let rows = [(1u64 << 36) | 1, 1, (1u64 << 35) | 1, 0];
        for &r in &rows {
            got.insert(r);
        }
        got.reduce();
        assert_eq!(
            got.rows
                .iter()
                .copied()
                .filter(|&r| r != 0)
                .collect::<Vec<_>>(),
            algebra::affine_rref(36, rows.into_iter())
        );
    }
    #[test]
    fn local_coordinate_rows_preserve_all_small_specializations() {
        let terms: Vec<_> = (0u64..8).filter(|m| m.count_ones() <= 2).collect();
        for chosen in 0usize..128 {
            let mut p: Vec<_> = terms
                .iter()
                .enumerate()
                .filter(|(i, _)| chosen & (1 << i) != 0)
                .map(|(_, &m)| m)
                .collect();
            algebra::canonical(&mut p);
            let system = vec![p.clone(), vec![3, 6, 1, 0], vec![], p, vec![1, 0]];
            let (reference, original) = QuadraticModel::compile(&system, 3).unwrap();
            let (packed, initial) = PackedModel::compile(&system, 3).unwrap();
            assert_eq!(packed.materialize(&initial), system);
            for mask in 0..8 {
                for values in (0..8).filter(|v| v & !mask == 0) {
                    let mut a = Logical::default();
                    let mut b = Logical::default();
                    let expected = reference.specialize(original.clone(), mask, values, &mut a);
                    let actual = packed.specialize(initial.clone(), mask, values, &mut b);
                    assert_eq!(
                        reference.materialize(&expected),
                        packed.materialize(&actual)
                    );
                    assert_eq!(a, b);
                    for next_mask in 0..8 {
                        for next_values in (0..8).filter(|v| v & !next_mask == 0) {
                            let mut a = Logical::default();
                            let mut b = Logical::default();
                            let left = reference.specialize(
                                expected.clone(),
                                next_mask,
                                next_values,
                                &mut a,
                            );
                            let right =
                                packed.specialize(actual.clone(), next_mask, next_values, &mut b);
                            assert_eq!(reference.materialize(&left), packed.materialize(&right));
                            assert_eq!(a, b);
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn packed_solver_preserves_trace_and_capacity_fallback() {
        for n in [12, 16] {
            for seed in [17, 937, 20261013, 1703939] {
                for family in ["planted", "cross_planted", "unplanted"] {
                    let (system, _) = fixture(n, seed, family);
                    let a = solve_quadratic(&system, n, 200000);
                    let b = solve_packed(&system, n, 200000);
                    assert_eq!(a.outcome, b.outcome);
                    assert_eq!(a.logical, b.logical);
                    assert_eq!(a.trace, b.trace);
                }
            }
        }
        let mut p: Vec<u64> = (0..256).filter(|m: &u64| m.count_ones() == 2).collect();
        algebra::canonical(&mut p);
        let system = vec![p];
        assert!(PackedModel::compile(&system, 36).is_none());
        let a = solve_quadratic(&system, 36, 1000);
        let b = solve_packed(&system, 36, 1000);
        assert_eq!(a.outcome, b.outcome);
        assert_eq!(a.logical, b.logical);
        assert_eq!(a.trace, b.trace);
    }
}

#[cfg(test)]
mod basis_tests {
    use super::*;
    include!("basis_tests.rs");
}

#[cfg(test)]
mod affine_tests {
    use super::*;
    include!("affine_tests.rs");
}

#[cfg(test)]
mod linearization_tests {
    use super::*;
    include!("linearization_tests.rs");
}

#[cfg(test)]
mod syndrome_tests {
    use super::*;
    include!("syndrome_tests.rs");
}
