//! Versioned width-generic IC: unknown orbit logs, exact witnesses and rank.
//! This is a bounded enumerated pair-table method, not an algebraic shortcut.
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use serde_json::json;

use super::canonical::{sha256_hex, short_id};
use super::claw::{ClawGroup, NarrowClaw, WideClaw};
use super::methods::{FactorBaseFacts, PhaseRecord, ResolvedMethod, SolveReport};
use super::workload::{CurveFacts, WideInstance};
use crate::cryptanalysis::ic_boundary::{BinaryInstance, GroupOps};
use crate::cryptanalysis::ic_measurement::{self as measurement, Phase};
use crate::cryptanalysis::koblitz_fast::FastPoint;
use crate::cryptanalysis::koblitz_strong_rho::{RawPointG, RhoScalar};

type Key = (bool, u128, u128);
type Counts = BTreeMap<String, u64>;

fn count(c: &mut Counts, name: &str, by: u64) {
    *c.entry(name.into()).or_default() += by;
}

#[derive(Clone, Copy)]
struct Shape {
    columns: usize,
    raw_x_cap: u64,
    max_trials: u64,
    max_pair_visits: u64,
    folded: bool,
}

impl Shape {
    fn parse(m: &ResolvedMethod) -> Result<Self, String> {
        let number = |k: &str| -> Result<u64, String> {
            m.params
                .get(k)
                .ok_or_else(|| format!("missing {k}"))?
                .parse()
                .map_err(|_| format!("invalid integer {k}"))
        };
        let s = Self {
            columns: usize::try_from(number("columns")?).map_err(|e| e.to_string())?,
            raw_x_cap: number("raw_x_cap")?,
            max_trials: number("max_trials")?,
            max_pair_visits: number("max_pair_visits")?,
            folded: match m.params.get("table_mode").map(String::as_str) {
                Some("orbit") => true,
                Some("direct") => false,
                _ => return Err("table_mode must be orbit or direct".into()),
            },
        };
        if !(1..=64).contains(&s.columns)
            || s.raw_x_cap == 0
            || s.raw_x_cap > 1_000_000
            || s.max_trials > 1_000_000
            || s.max_pair_visits == 0
            || s.max_pair_visits > 2_000_000
        {
            return Err("orbit_rank_v1 requires columns1..64, raw_x_cap1..1000000, max_trials0..1000000 and max_pair_visits1..2000000".into());
        }
        Ok(s)
    }
}

#[derive(Clone, Copy, Debug)]
struct Term<S> {
    column: usize,
    coefficient: S,
}

#[derive(Clone, Copy)]
struct PointTerm<P, S> {
    point: P,
    term: Term<S>,
}

#[derive(Clone, Copy)]
struct Pair<P, S> {
    point: P,
    terms: [Term<S>; 2],
    canonical_shift: u32,
}

/// Exact modular RREF. Scalars and field elements have independent widths.
struct Rank<S: RhoScalar> {
    modulus: S,
    rows: BTreeMap<usize, Vec<S>>,
    columns: usize,
}

impl<S: RhoScalar> Rank<S> {
    fn new(modulus: S, columns: usize) -> Self {
        Self {
            modulus,
            rows: BTreeMap::new(),
            columns,
        }
    }

    fn subtract(modulus: S, a: &mut [S], b: &[S], f: S, c: &mut Counts) {
        for (x, &y) in a.iter_mut().zip(b) {
            count(c, "matrix_entries_scanned_uncharged", 1);
            if y != S::zero() {
                *x = S::sub_mod(*x, S::mul_mod(f, y, modulus), modulus);
                count(c, "matrix_multiply_subtracts_uncharged", 1);
            }
        }
    }

    fn insert(&mut self, mut row: Vec<S>, c: &mut Counts) -> Result<(), String> {
        for (&pivot, existing) in &self.rows {
            let f = row[pivot];
            if f != S::zero() {
                Self::subtract(self.modulus, &mut row, existing, f, c);
            }
        }
        let Some(pivot) = (0..self.columns).find(|&i| row[i] != S::zero()) else {
            if row[self.columns] != S::zero() {
                return Err("inconsistent independently point-checked relation matrix".into());
            }
            count(c, "dependent_rows", 1);
            return Ok(());
        };
        count(c, "scalar_inversions_uncharged", 1);
        let inv = S::inverse_mod(row[pivot], self.modulus).ok_or("noninvertible pivot")?;
        for x in &mut row {
            *x = S::mul_mod(*x, inv, self.modulus);
            count(c, "matrix_normalizations_uncharged", 1);
        }
        for existing in self.rows.values_mut() {
            let f = existing[pivot];
            if f != S::zero() {
                Self::subtract(self.modulus, existing, &row, f, c);
            }
        }
        self.rows.insert(pivot, row);
        count(c, "independent_rows", 1);
        Ok(())
    }

    fn pinned(&self, column: usize) -> Option<S> {
        let row = self.rows.get(&column)?;
        (0..self.columns)
            .all(|i| i == column || row[i] == S::zero())
            .then_some(row[self.columns])
    }
}

fn phase(name: &str, ops: GroupOps) -> PhaseRecord {
    PhaseRecord {
        name: name.into(),
        adds: ops.adds,
        doubles: ops.doubles,
        scalar_mults: ops.scalar_mults,
        gae: ops.gae(),
        ..Default::default()
    }
}

fn narrow_key(p: FastPoint) -> Key {
    (p.infinity, u128::from(p.x), u128::from(p.y))
}

fn wide_key(p: RawPointG<u128>) -> Key {
    match p {
        RawPointG::Infinity => (true, 0, 0),
        RawPointG::Affine { x, y } => (false, x, y),
    }
}

pub fn solve_narrow(
    m: &ResolvedMethod,
    inst: &BinaryInstance,
    curve: &CurveFacts,
    target: FastPoint,
    seed: u64,
) -> Result<SolveReport, String> {
    let g = NarrowClaw::new(inst).ok_or("orbit_rank_v1 needs a prime-degree Koblitz curve")?;
    run(
        &g,
        inst.cofactor,
        |x| inst.points_with_x(x),
        narrow_key,
        m,
        curve,
        target,
        seed,
    )
}

pub fn solve_wide(
    m: &ResolvedMethod,
    inst: &WideInstance,
    curve: &CurveFacts,
    target: RawPointG<u128>,
    seed: u64,
) -> Result<SolveReport, String> {
    let g = WideClaw::new(inst);
    run(
        &g,
        inst.kc.cofactor,
        |x| {
            inst.kc
                .points_with_x(u128::from(x))
                .into_iter()
                .map(|(x, y)| RawPointG::Affine { x, y })
                .collect()
        },
        wide_key,
        m,
        curve,
        target,
        seed,
    )
}

/// The same collector at either width. No planted scalar or factor-base logs.
#[allow(clippy::too_many_arguments)]
fn run<G: ClawGroup>(
    g: &G,
    cofactor: G::S,
    lift: impl Fn(u64) -> Vec<G::P>,
    key: impl Fn(G::P) -> Key,
    method: &ResolvedMethod,
    curve: &CurveFacts,
    target: G::P,
    seed: u64,
) -> Result<SolveReport, String> {
    let shape = Shape::parse(method)?;
    let started = Instant::now();
    let (r, n) = (g.r(), g.n());
    let mut c = Counts::new();
    let mut base_ops = GroupOps::default();
    let mut table_ops = GroupOps::default();
    let mut query_ops = GroupOps::default();
    let mut pdp_ops = GroupOps::default();
    let mut check_ops = GroupOps::default();
    let mut columns = Vec::new();
    let mut classes = BTreeSet::new();
    measurement::mark(Phase::FactorBase);
    for x in 1..=shape.raw_x_cap {
        count(&mut c, "raw_x_lifts_uncharged", 1);
        let Some(p) = lift(x).into_iter().next() else {
            continue;
        };
        let p = g.counted_mul(&mut base_ops, p, cofactor);
        if g.is_identity(&p) {
            continue;
        }
        count(&mut c, "canonicalisations_uncharged", 1);
        let class = g.class(&p).0;
        if classes.insert(class) {
            columns.push(p);
            if columns.len() == shape.columns {
                break;
            }
        }
    }
    if columns.len() != shape.columns {
        return Err(format!(
            "raw-x cap found {} of {} requested orbit columns",
            columns.len(),
            shape.columns
        ));
    }
    let mut powers = Vec::with_capacity(n as usize);
    let mut power = G::S::one();
    for _ in 0..n {
        powers.push(power);
        power = G::S::mul_mod(power, g.lambda(), r);
        count(&mut c, "scalar_products_uncharged", 1);
    }
    if power != G::S::one() {
        return Err("Frobenius eigenvalue fails lambda^n=1".into());
    }
    let mut points = Vec::new();
    let mut point_keys = BTreeSet::new();
    for (column, &p) in columns.iter().enumerate() {
        for (t, &coefficient) in powers.iter().enumerate() {
            let q = g.frob(t as u32, p);
            count(&mut c, "frobenius_maps_uncharged", 1);
            for (point, coefficient) in [
                (q, coefficient),
                (g.neg(q), G::S::sub_mod(G::S::zero(), coefficient, r)),
            ] {
                if point_keys.insert(key(point)) {
                    points.push(PointTerm {
                        point,
                        term: Term {
                            column,
                            coefficient,
                        },
                    });
                }
            }
        }
    }
    c.insert("factor_base_columns".into(), columns.len() as u64);
    c.insert("factor_base_signed_points".into(), points.len() as u64);
    // The base identity excludes the oracle choice and target-dependent seed.
    let inventory: Vec<_> = point_keys
        .iter()
        .map(|&(inf, x, y)| json!([inf, format!("0x{x:x}"), format!("0x{y:x}")]))
        .collect();
    let points_sha256 = sha256_hex(super::canonical::canonical(&json!(inventory))?.as_bytes());
    let params = BTreeMap::from([
        ("columns".into(), shape.columns.to_string()),
        ("raw_x_cap".into(), shape.raw_x_cap.to_string()),
        (
            "key_encoding".into(),
            "canonical-json-sorted-full-point-tuples-hex-v1".into(),
        ),
    ]);
    let identity = json!({"schema":"ecbench.factor_base/v1", "curve":curve.slug,
        "family":"projected-raw-x-orbits-v1", "params":params,
        "columns":columns.len(), "signed_points":points.len(), "points_sha256":points_sha256});
    let (fb_id, fb_sha256) = short_id("FB1", &identity)?;
    let base = FactorBaseFacts { fb_id, fb_sha256, family:"projected-raw-x-orbits-v1".into(),
        params, description:"first raw x lifts, cofactor projected, distinct signed Frobenius orbits; logarithms unknown".into(),
        signed_points:points.len() as u64, abscissae:point_keys.iter().map(|k|k.1).collect::<BTreeSet<_>>().len() as u64,
        columns:columns.len() as u64, dimension:None, points_sha256 };

    measurement::mark(Phase::Precompute);
    let mut quotient = BTreeMap::<u128, Pair<G::P, G::S>>::new();
    let mut direct = BTreeMap::<Key, Pair<G::P, G::S>>::new();
    let mut visits = 0;
    let mut table_complete = true;
    if shape.folded {
        'pairs: for (i, &p) in columns.iter().enumerate() {
            for (j, &q) in columns.iter().enumerate() {
                for (t, &lambda) in powers.iter().enumerate() {
                    let f = g.frob(t as u32, q);
                    count(&mut c, "frobenius_maps_uncharged", 1);
                    for negative in [false, true] {
                        if visits == shape.max_pair_visits {
                            table_complete = false;
                            break 'pairs;
                        }
                        visits += 1;
                        let sum = g.add(&mut table_ops, p, if negative { g.neg(f) } else { f });
                        count(&mut c, "canonicalisations_uncharged", 1);
                        let (class, canonical_shift) = g.class(&sum);
                        let coefficient = if negative {
                            G::S::sub_mod(G::S::zero(), lambda, r)
                        } else {
                            lambda
                        };
                        quotient.entry(class).or_insert(Pair {
                            point: sum,
                            terms: [
                                Term {
                                    column: i,
                                    coefficient: G::S::one(),
                                },
                                Term {
                                    column: j,
                                    coefficient,
                                },
                            ],
                            canonical_shift,
                        });
                        count(&mut c, "table_inserts_uncharged", 1);
                    }
                }
            }
        }
    } else {
        'pairs: for (i, p) in points.iter().enumerate() {
            for q in &points[i..] {
                if visits == shape.max_pair_visits {
                    table_complete = false;
                    break 'pairs;
                }
                visits += 1;
                let sum = g.add(&mut table_ops, p.point, q.point);
                direct.entry(key(sum)).or_insert(Pair {
                    point: sum,
                    terms: [p.term, q.term],
                    canonical_shift: 0,
                });
                count(&mut c, "table_inserts_uncharged", 1);
            }
        }
    }
    c.insert("pair_visits".into(), visits);
    c.insert(
        "table_entries".into(),
        (quotient.len() + direct.len()) as u64,
    );
    c.insert("table_complete".into(), u64::from(table_complete));
    let mut rank = Rank::new(r, columns.len() + 1);
    let mut recovered = None;
    let mut witnesses = Vec::new();
    if table_complete {
        measurement::begin_online(Phase::TargetQuery);
        'queries: for trial in 0..shape.max_trials {
            measurement::mark(Phase::TargetQuery);
            let a = g.seed_scalar(seed, 2 * trial);
            let b = g.seed_scalar(seed, 2 * trial + 1);
            let ag = g.counted_mul(&mut query_ops, g.generator(), a);
            let bq = g.counted_mul(&mut query_ops, target, b);
            let query = g.add(&mut query_ops, ag, bq);
            count(&mut c, "trials", 1);
            measurement::mark(Phase::TargetPdp);
            for third in &points {
                let residual = g.add(&mut pdp_ops, query, g.neg(third.point));
                count(&mut c, "table_probes_uncharged", 1);
                let mut terms = vec![third.term];
                if !g.is_identity(&residual) {
                    let pair_terms = if shape.folded {
                        count(&mut c, "canonicalisations_uncharged", 1);
                        let (class, shift) = g.class(&residual);
                        let Some(pair) = quotient.get(&class) else {
                            continue;
                        };
                        let rotation = (pair.canonical_shift + n - shift) % n;
                        let p = g.frob(rotation, pair.point);
                        count(&mut c, "frobenius_maps_uncharged", 1);
                        let factor = if p == residual {
                            powers[rotation as usize]
                        } else if g.neg(p) == residual {
                            G::S::sub_mod(G::S::zero(), powers[rotation as usize], r)
                        } else {
                            return Err("quotient witness failed full-point lift".into());
                        };
                        pair.terms.map(|term| Term {
                            column: term.column,
                            coefficient: G::S::mul_mod(factor, term.coefficient, r),
                        })
                    } else {
                        let Some(pair) = direct.get(&key(residual)) else {
                            continue;
                        };
                        pair.terms
                    };
                    terms.extend(pair_terms);
                }
                measurement::mark(Phase::TargetRelationCheck);
                let mut sum = g.identity();
                let mut row = vec![G::S::zero(); columns.len() + 2];
                for term in &terms {
                    let p = g.counted_mul(&mut check_ops, columns[term.column], term.coefficient);
                    sum = g.add(&mut check_ops, sum, p);
                    row[term.column] = G::S::add_mod(row[term.column], term.coefficient, r);
                    count(&mut c, "scalar_products_uncharged", 1);
                }
                if sum != query {
                    return Err("relation failed independent full-point sum".into());
                }
                count(&mut c, "point_checked_relations", 1);
                row[columns.len()] = G::S::sub_mod(G::S::zero(), b, r);
                row[columns.len() + 1] = a;
                if witnesses.len() < 16 {
                    witnesses.push(json!({"trial":trial,"a":a.to_string(),"b":b.to_string(),
                        "query":key(query), "terms":terms.iter().map(|t|
                        json!({"column":t.column,"coefficient":t.coefficient.to_string()})).collect::<Vec<_>>() }));
                }
                measurement::mark(Phase::TargetDescent);
                rank.insert(row, &mut c)?;
                if let Some(d) = rank.pinned(columns.len()) {
                    recovered = Some(d.to_u128());
                    break 'queries;
                }
                break;
            }
        }
        measurement::end_online();
    }
    c.insert("rank".into(), rank.rows.len() as u64);
    let phases = vec![
        phase("factor_base", base_ops),
        phase("oracle_setup", table_ops),
        phase("target_query", query_ops),
        phase("target_PDP", pdp_ops),
        phase("target_relation_check", check_ops),
        phase("linear_algebra", GroupOps::default()),
    ];
    let total_gae = phases.iter().map(|p| p.gae).sum();
    let unpriced = c
        .iter()
        .filter(|(k, v)| k.ends_with("_uncharged") && **v > 0)
        .map(|(k, _)| k.clone())
        .collect();
    Ok(SolveReport {
        recovered,
        exhausted: recovered.is_none(),
        phases,
        total_gae,
        automorphisms_used: if shape.folded { 2 * n } else { 1 },
        counters: c,
        unpriced,
        deterministic: true,
        solve_wall_ns: started.elapsed().as_nanos() as u64,
        factor_base: Some(base),
        detail: json!({"schema":"ecbench.orbit-rank/v1", "table_mode":if shape.folded {"orbit"} else {"direct"},
            "base_columns":columns.iter().map(|&p|key(p)).collect::<Vec<_>>(),"relation_witnesses":witnesses,
            "known_base_logs":false,"table_complete":table_complete,
            "stop":if !table_complete {"pair_visit_cap"} else if recovered.is_some() {"target_pinned"} else {"query_cap"},
            "cost_scope":"counted group operations; native work unpriced; no full-cost or wall speedup"}),
        ..Default::default()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn method(mode: &str, columns: u64, trials: u64) -> ResolvedMethod {
        super::super::methods::resolve(&super::super::methods::MethodSpec {
            id: "ic.orbit_rank_v1".into(),
            params: BTreeMap::from([
                ("columns".into(), columns.to_string()),
                ("raw_x_cap".into(), "4096".into()),
                ("max_trials".into(), trials.to_string()),
                ("table_mode".into(), mode.into()),
                ("max_pair_visits".into(), "1000000".into()),
            ]),
        })
        .unwrap()
    }

    #[test]
    fn exact_and_quotient_collectors_recover_and_replay_same_points() {
        use super::super::workload::{CurveSpec, Instance};
        use crate::cryptanalysis::ic_boundary::{BinaryGroup, CountedGroup};
        let spec = CurveSpec::Koblitz { a: 0, n: 13 };
        let inst = spec.build().unwrap();
        let facts = inst.facts(&spec);
        let Instance::Binary(ref i) = inst else {
            unreachable!()
        };
        let g = BinaryGroup(&i.fast);
        for scalar in [1, 123, 321] {
            let target = g.mul(&mut GroupOps::default(), i.generator, scalar);
            let direct = solve_narrow(&method("direct", 2, 512), i, &facts, target, 7351).unwrap();
            let orbit = solve_narrow(&method("orbit", 2, 512), i, &facts, target, 7351).unwrap();
            assert_eq!(direct.recovered, Some(u128::from(scalar)));
            assert_eq!(orbit.recovered, direct.recovered);
            assert_eq!(direct.factor_base, orbit.factor_base);
            let replay = solve_narrow(&method("orbit", 2, 512), i, &facts, target, 7351).unwrap();
            assert_eq!(orbit.counters, replay.counters);
            assert_eq!(orbit.phases, replay.phases);
            assert_eq!(orbit.detail, replay.detail);
        }
    }

    #[test]
    fn same_collector_recovers_with_wide_field_and_scalar_types() {
        use super::super::workload::{CurveSpec, Instance};
        use crate::cryptanalysis::{koblitz_strong_rho::WideStrongRho, wide_gf2m::WideGf2};
        let spec = CurveSpec::Koblitz { a: 0, n: 13 };
        let inst = spec.build().unwrap();
        let facts = inst.facts(&spec);
        let Instance::Binary(ref i) = inst else {
            unreachable!()
        };
        let narrow = NarrowClaw::new(i).unwrap();
        let rho = WideStrongRho::from_parts(
            WideGf2::new(&i.irreducible),
            u128::from(i.a),
            u128::from(i.r),
            u128::from(narrow.lambda()),
            RawPointG::Affine {
                x: u128::from(i.generator.x),
                y: u128::from(i.generator.y),
            },
        );
        let wide = WideClaw::from_parts(&rho, u128::from(i.r), u128::from(narrow.lambda()));
        let target = wide.counted_mul(&mut GroupOps::default(), wide.generator(), 123);
        let rep = run(
            &wide,
            u128::from(i.cofactor),
            |x| {
                i.points_with_x(x)
                    .into_iter()
                    .map(|p| RawPointG::Affine {
                        x: u128::from(p.x),
                        y: u128::from(p.y),
                    })
                    .collect()
            },
            wide_key,
            &method("orbit", 2, 512),
            &facts,
            target,
            7351,
        )
        .unwrap();
        assert_eq!(rep.recovered, Some(123));
    }

    #[test]
    fn deterministic_pair_cap_preserves_incomplete_setup_work() {
        use super::super::workload::{CurveSpec, Instance};
        let spec = CurveSpec::Koblitz { a: 0, n: 13 };
        let inst = spec.build().unwrap();
        let facts = inst.facts(&spec);
        let Instance::Binary(ref i) = inst else {
            unreachable!()
        };
        let mut m = method("orbit", 2, 512);
        m.params.insert("max_pair_visits".into(), "3".into());
        let rep = solve_narrow(&m, i, &facts, i.generator, 7351).unwrap();
        assert!(rep.exhausted);
        assert_eq!(rep.counters["pair_visits"], 3);
        assert_eq!(rep.detail["stop"], "pair_visit_cap");
        assert!(rep.total_gae > 0.);
    }

    #[test]
    fn exact_n83_witness_lifts_with_coefficients_past_u64() {
        use super::super::workload::{CurveSpec, Instance};
        let spec = CurveSpec::Koblitz { a: 0, n: 83 };
        let inst = spec.build().unwrap();
        let facts = inst.facts(&spec);
        let Instance::Wide(ref i) = inst else {
            unreachable!()
        };
        let g = WideClaw::new(i);
        let seed = 7351;
        let (a, b) = (g.seed_scalar(seed, 0), g.seed_scalar(seed, 1));
        let desired = <u128 as RhoScalar>::add_mod(i.kc.lambda, 2, i.kc.r);
        let scalar = <u128 as RhoScalar>::mul_mod(
            <u128 as RhoScalar>::sub_mod(desired, a, i.kc.r),
            <u128 as RhoScalar>::inverse_mod(b, i.kc.r).unwrap(),
            i.kc.r,
        );
        let target = g.counted_mul(&mut GroupOps::default(), g.generator(), scalar);
        let rep = solve_wide(&method("orbit", 2, 1), i, &facts, target, seed).unwrap();
        assert_eq!(rep.counters["point_checked_relations"], 1);
        assert_eq!(rep.counters["rank"], 1);
        assert!(rep.detail["relation_witnesses"][0]["terms"]
            .as_array()
            .unwrap()
            .iter()
            .any(
                |t| t["coefficient"].as_str().unwrap().parse::<u128>().unwrap()
                    > u128::from(u64::MAX)
            ));
        assert!(rep.exhausted); // A planted relation is not a full target recovery.
    }

    #[test]
    fn wide_rank_keeps_modulus_and_coefficients_above_one_word() {
        let p = 2_417_851_639_230_796_216_685_689u128;
        let mut rank = Rank::new(p, 2);
        let (x, y) = ((1u128 << 75) + 39, (1u128 << 73) + 7);
        let rhs = <u128 as RhoScalar>::add_mod(x, y, p);
        let mut counts = Counts::new();
        rank.insert(vec![1, 1, rhs], &mut counts).unwrap();
        assert_eq!(rank.pinned(1), None);
        rank.insert(vec![1, 0, x], &mut counts).unwrap();
        assert_eq!(rank.pinned(0), Some(x));
        assert_eq!(rank.pinned(1), Some(y));
        assert!(rank.insert(vec![0, 0, 1], &mut counts).is_err());
    }
}
