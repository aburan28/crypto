//! Complete bounded generic Boolean solves with matched search semantics.
#[allow(dead_code)]
mod algebra {
    include!("kernel.rs");
}
use std::hint::black_box;
use std::time::Instant;

type System = Vec<Vec<u64>>;
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
    specialized_terms: u64,
    source_rows: u64,
    source_columns: u64,
    max_depth: usize,
}
#[derive(Default)]
struct Profile {
    kernel_ns: u128,
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
}
impl Solver {
    fn event(&mut self, value: u64) {
        self.trace = (self.trace ^ value).wrapping_mul(0x100000001b3);
    }
    fn specialize(&mut self, system: System, mask: u64, values: u64) -> System {
        let zero = mask & !values;
        let mut out = Vec::new();
        for poly in system {
            self.logical.specialized_terms += poly.len() as u64;
            let mut terms: Vec<_> = poly
                .into_iter()
                .filter(|m| m & zero == 0)
                .map(|m| m & !mask)
                .collect();
            algebra::canonical(&mut terms);
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
                self.logical.kernel_calls += 1;
                self.profile.kernel_calls_by_active[width] += 1;
                let tick = Instant::now();
                let (original, reply) = self.session.as_mut().unwrap().run(self.n, system);
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
struct Solved {
    outcome: Outcome,
    logical: Logical,
    profile: Profile,
    trace: u64,
}
fn solve(system: &System, n: u8, arm: &str, limit: u64) -> Solved {
    let session = match arm {
        "search" => None,
        "flat" => Some(algebra::Session::new(algebra::Method::Flat)),
        "bucket" => Some(algebra::Session::new(algebra::Method::Bucket)),
        "hybrid" => Some(algebra::Session::new(algebra::Method::Hybrid)),
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
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 6, "worker N SEED FAMILY REPETITIONS NODE_CAP");
    let n = args[1].parse::<u8>().unwrap();
    let seed = args[2].parse::<u64>().unwrap();
    let family = &args[3];
    let reps = args[4].parse::<usize>().unwrap();
    let limit = args[5].parse::<u64>().unwrap();
    assert!([12, 16, 20, 24].contains(&n) && (1..=8).contains(&reps) && limit <= 200_000);
    let (system, planted) = fixture(n, seed, family);
    let reference = solve(&system, n, "search", limit);
    let planted_json = planted.map_or("null".into(), |x| x.to_string());
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"family\":\"{family}\",\"polys\":{:?},\"planted_witness\":{planted_json},\"search_reference\":\"{}\"}}",system,status(&reference.outcome));
    let names = ["search", "flat", "bucket", "hybrid"];
    let mut expected_kernel: Option<(Outcome, Logical, u64)> = None;
    for rep in 0..reps {
        for order in 0..4 {
            let arm = names[(rep + order) % 4];
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
            if arm != "search" {
                let signature = (got.outcome.clone(), got.logical.clone(), got.trace);
                if let Some(expected) = &expected_kernel {
                    assert_eq!(&signature, expected);
                } else {
                    expected_kernel = Some(signature);
                }
            }
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
            let trace = got.trace;
            println!("{{\"type\":\"sample\",\"rep\":{rep},\"order\":{order},\"variant\":\"{arm}\",\"outcome\":\"{outcome}\",\"model\":{model},\"reason\":{reason},\"solve_ns\":{solve_ns},\"validation_ns\":{validation_ns},\"total_ns\":{total_ns},\"kernel_ns\":{},\"nodes\":{},\"kernel_calls\":{},\"decisions\":{},\"forced\":{},\"specialized_terms\":{},\"source_rows\":{},\"source_columns\":{},\"max_depth\":{},\"trace\":{trace},\"flat_calls\":{},\"kernel_calls_by_active\":{:?}}}",profile.kernel_ns,l.nodes,l.kernel_calls,l.decisions,l.forced,l.specialized_terms,l.source_rows,l.source_columns,l.max_depth,profile.flat_calls,profile.kernel_calls_by_active);
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
                for arm in ["search", "flat", "bucket", "hybrid"] {
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
        for arm in ["search", "flat", "bucket", "hybrid"] {
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
        for arm in ["search", "flat", "bucket", "hybrid"] {
            let answer = solve(&system, 3, arm, 100);
            if let Outcome::Sat(model) = answer.outcome {
                assert!(satisfies(&system, model));
            } else {
                panic!("linear system should be satisfiable");
            }
        }
    }
}
