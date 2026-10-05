//! The breadth-first walk over the `ℓ`-isogeny graphs of a prime-field
//! curve's `F_p`-isogeny class, and the records it writes.

use std::collections::HashMap;
use std::time::Instant;

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use rayon::prelude::*;

use super::curve::{self, Model, OrderAudit};
use super::field::{is_probable_prime, Fe, Field};
use super::kernel::{self, EdgeError, Isogeny};
use super::modpoly::{self, ModPoly};
use super::poly::{self, Poly};
use super::record::{self, V};
use super::traits;
use crate::cryptanalysis::curve_id;
use crate::ecc::curve::CurveParams;

/// Trial-division bound for the discriminant and twist factorisations.
pub const TRIAL_BOUND: u64 = 1 << 20;
/// Embedding degrees searched.
pub const EMBEDDING_SEARCH: u64 = 1000;

/// The walk's root: a registered curve with its generator.
#[derive(Clone, Debug)]
pub struct StartCurve {
    /// Standard name (`P-256`), or a free label for a custom curve.
    pub name: String,
    /// EC1 alias tag: the standard name lower-cased, else `fp`.
    pub tag: String,
    pub p: BigUint,
    pub a: BigUint,
    pub b: BigUint,
    /// `#E(F_p)`.
    pub order: BigUint,
    pub subgroup_order: BigUint,
    pub cofactor: BigUint,
    pub gx: BigUint,
    pub gy: BigUint,
}

impl StartCurve {
    pub fn from_params(c: &CurveParams, standard: bool) -> Self {
        let tag = if standard {
            c.name
                .to_ascii_lowercase()
                .replace(|ch: char| !ch.is_ascii_alphanumeric(), "")
        } else {
            "fp".into()
        };
        let h = BigUint::from(c.h);
        StartCurve {
            name: c.name.to_string(),
            tag,
            p: c.p.clone(),
            a: c.a.clone(),
            b: c.b.clone(),
            order: &c.n * &h,
            subgroup_order: c.n.clone(),
            cofactor: h,
            gx: c.gx.clone(),
            gy: c.gy.clone(),
        }
    }

    pub fn p256() -> Self {
        Self::from_params(&CurveParams::p256(), true)
    }

    pub fn p192() -> Self {
        Self::from_params(&CurveParams::p192(), true)
    }

    pub fn p224() -> Self {
        Self::from_params(&CurveParams::p224(), true)
    }
}

#[derive(Clone, Debug)]
pub struct WalkConfig {
    /// Odd primes `ℓ` whose graphs are walked.
    pub primes: Vec<u64>,
    /// Stop expanding once this many curves are known.
    pub max_curves: usize,
    /// Do not expand curves at this depth or deeper.
    pub max_depth: usize,
    /// Root-finding seed (results do not depend on it).
    pub seed: u64,
    /// Points per order audit.
    pub audit_points: usize,
    /// Primes requested but not walked, with the reason (an Atkin prime
    /// has no rational `ℓ`-isogeny).
    pub skipped: Vec<(u64, String)>,
}

impl Default for WalkConfig {
    fn default() -> Self {
        WalkConfig {
            primes: vec![3, 5, 7, 11, 13],
            max_curves: 64,
            max_depth: usize::MAX,
            seed: 1,
            audit_points: 1,
            skipped: Vec::new(),
        }
    }
}

/// How `ℓ` splits in `Z[π]`, from `(D_π / ℓ)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EllKind {
    Elkies,
    Atkin,
    Divides,
}

impl EllKind {
    fn name(self) -> &'static str {
        match self {
            EllKind::Elkies => "elkies",
            EllKind::Atkin => "atkin",
            EllKind::Divides => "divides_discriminant",
        }
    }
}

/// Invariants of the whole `F_p`-isogeny class: they hold for every curve
/// the walk can reach.
#[derive(Clone, Debug)]
pub struct ClassInfo {
    pub trace: BigInt,
    /// `D_π = t² − 4p`.
    pub disc: BigInt,
    pub disc_factors: Vec<(u64, u32)>,
    pub disc_cofactor: BigUint,
    pub disc_cofactor_status: &'static str,
    /// `D_π` squarefree and `≡ 1 (mod 4)` given the cofactor's primality:
    /// then `Z[π]` is maximal and `End(E) = Z[π]` for every curve.
    pub disc_fundamental_if_cofactor_prime: bool,
    pub order_prime: bool,
    pub twist_order: BigUint,
    pub twist_factors: Vec<(u64, u32)>,
    pub twist_cofactor: BigUint,
    pub twist_cofactor_status: &'static str,
    pub embedding_degree: Option<u64>,
    pub ells: Vec<(u64, EllKind, u32)>,
}

fn trial_factor(n: &BigUint, bound: u64) -> (Vec<(u64, u32)>, BigUint) {
    let mut n = n.clone();
    let mut out = Vec::new();
    let mut sieve = vec![true; bound as usize + 1];
    for q in 2..=bound {
        if !sieve[q as usize] {
            continue;
        }
        let mut m = q * q;
        while m <= bound {
            sieve[m as usize] = false;
            m += q;
        }
        let mut e = 0;
        while (&n % q).is_zero() {
            n /= q;
            e += 1;
        }
        if e > 0 {
            out.push((q, e));
        }
    }
    (out, n)
}

fn cofactor_status(c: &BigUint) -> &'static str {
    if c.is_one() {
        "one"
    } else if is_probable_prime(c) {
        "probable_prime"
    } else if c.sqrt().pow(2) == *c {
        "perfect_square"
    } else {
        "composite_unfactored"
    }
}

impl ClassInfo {
    pub fn compute(start: &StartCurve, primes: &[u64]) -> Self {
        let p = BigInt::from(start.p.clone());
        let trace = &p + 1 - BigInt::from(start.order.clone());
        let disc: BigInt = &trace * &trace - &p * 4u32;
        let (disc_factors, disc_cofactor) = trial_factor(disc.magnitude(), TRIAL_BOUND);
        let disc_cofactor_status = cofactor_status(&disc_cofactor);
        let disc_mod4 = disc.mod_floor(&BigInt::from(4));
        let disc_fundamental_if_cofactor_prime = disc_factors.iter().all(|(_, e)| *e == 1)
            && matches!(disc_cofactor_status, "one" | "probable_prime")
            && disc_mod4 == BigInt::one();
        let twist_order = (&start.p << 1) + 2u8 - &start.order;
        let (twist_factors, twist_cofactor) = trial_factor(&twist_order, TRIAL_BOUND);
        let twist_cofactor_status = cofactor_status(&twist_cofactor);
        let r = &start.subgroup_order;
        let embedding_degree = if r.is_zero() {
            None
        } else {
            let pr = &start.p % r;
            let mut acc = pr.clone();
            (1..=EMBEDDING_SEARCH).find(|_| {
                let hit = acc.is_one();
                acc = (&acc * &pr) % r;
                hit
            })
        };
        let ells = primes
            .iter()
            .map(|&l| {
                let v = disc_factors
                    .iter()
                    .find(|(q, _)| *q == l)
                    .map_or(0, |(_, e)| *e);
                let dm = disc.mod_floor(&BigInt::from(l)).to_u64().unwrap();
                let kind = if dm == 0 {
                    EllKind::Divides
                } else if BigUint::from(dm)
                    .modpow(&BigUint::from((l - 1) / 2), &BigUint::from(l))
                    .is_one()
                {
                    EllKind::Elkies
                } else {
                    EllKind::Atkin
                };
                (l, kind, v)
            })
            .collect();
        ClassInfo {
            trace,
            disc,
            disc_factors,
            disc_cofactor,
            disc_cofactor_status,
            disc_fundamental_if_cofactor_prime,
            order_prime: is_probable_prime(&start.order),
            twist_order,
            twist_factors,
            twist_cofactor,
            twist_cofactor_status,
            embedding_degree,
            ells,
        }
    }

    /// `v_ℓ(D_π)` for a walked `ℓ`.
    pub fn v_ell(&self, ell: u64) -> u32 {
        self.ells.iter().find(|e| e.0 == ell).map_or(0, |e| e.2)
    }

    /// The `ℓ`-volcano's depth `v_ℓ(f_π) = ⌊v_ℓ(D_π)/2⌋` (odd `ℓ`).
    pub fn depth(&self, ell: u64) -> u32 {
        self.v_ell(ell) / 2
    }

    /// Every curve sits at level 0 of a depth-0 `ℓ`-volcano when
    /// `v_ℓ(D_π) ≤ 1` (odd `ℓ`): `ℓ` does not divide `[O_K : Z[π]]`.
    pub fn level_zero_proved(&self, ell: u64) -> bool {
        ell % 2 == 1 && self.v_ell(ell) <= 1
    }

    /// Rational `ℓ`-isogenies a curve has when the volcano has depth 0.
    pub fn expected_neighbours(&self, ell: u64) -> Option<usize> {
        if !self.level_zero_proved(ell) {
            return None;
        }
        self.ells.iter().find(|e| e.0 == ell).map(|e| match e.1 {
            EllKind::Elkies => 2,
            EllKind::Atkin => 0,
            EllKind::Divides => 1,
        })
    }

    pub fn record(&self) -> V {
        let factors = |fs: &[(u64, u32)]| {
            V::Seq(
                fs.iter()
                    .map(|(q, e)| V::Seq(vec![V::int(q), V::int(e)]))
                    .collect(),
            )
        };
        V::map(vec![
            ("trace", V::bigi(&self.trace)),
            ("frobenius_discriminant", V::bigi(&self.disc)),
            (
                "frobenius_discriminant_trial_factors",
                factors(&self.disc_factors),
            ),
            ("frobenius_discriminant_trial_bound", V::int(TRIAL_BOUND)),
            (
                "frobenius_discriminant_cofactor",
                V::big(&self.disc_cofactor),
            ),
            (
                "frobenius_discriminant_cofactor_status",
                V::s(self.disc_cofactor_status),
            ),
            (
                "frobenius_discriminant_fundamental_if_cofactor_prime",
                V::Bool(self.disc_fundamental_if_cofactor_prime),
            ),
            ("order_probable_prime", V::Bool(self.order_prime)),
            ("twist_order", V::big(&self.twist_order)),
            ("twist_order_trial_factors", factors(&self.twist_factors)),
            ("twist_order_cofactor", V::big(&self.twist_cofactor)),
            (
                "twist_order_cofactor_status",
                V::s(self.twist_cofactor_status),
            ),
            (
                "embedding_degree",
                V::opt(self.embedding_degree.map(V::int)),
            ),
            ("embedding_degree_searched_to", V::int(EMBEDDING_SEARCH)),
            (
                "ell",
                V::Seq(
                    self.ells
                        .iter()
                        .map(|(l, k, v)| {
                            V::map(vec![
                                ("ell", V::int(l)),
                                ("splitting", V::s(k.name())),
                                ("v_ell_discriminant", V::int(v)),
                                ("volcano_depth", V::int(self.depth(*l))),
                                ("level_zero_proved", V::Bool(self.level_zero_proved(*l))),
                                (
                                    "expected_rational_neighbours",
                                    V::opt(self.expected_neighbours(*l).map(V::int)),
                                ),
                            ])
                        })
                        .collect(),
                ),
            ),
        ])
    }
}

/// A curve the walk found.
#[derive(Clone, Debug)]
pub struct Node {
    pub j: Fe,
    pub model: Model,
    pub depth: usize,
    /// The BFS-tree edge into this node.
    pub parent_edge: Option<usize>,
    pub expanded: bool,
    pub generator: (Fe, Fe),
    pub audit: Option<OrderAudit>,
    /// `(ℓ, rational roots of Φ_ℓ(j, Y))`, once expanded.
    pub roots: Vec<(u64, usize)>,
    /// Proved `(ℓ, level, proof)`: `v_ℓ` of the conductor of `End(E)`.
    pub levels: Vec<(u64, u32, String)>,
}

/// A verified edge `source → target` of degree `ℓ`.
#[derive(Clone, Debug)]
pub struct Edge {
    pub ell: u64,
    /// `h` when every curve is proved to sit at level 0, else `x`.
    pub dir: char,
    pub source: usize,
    pub target: usize,
    pub kernel: Poly,
    pub velu: Model,
    /// `u²` of the isomorphism from Vélu's codomain onto the target model.
    pub u2: Fe,
}

/// A root of `Φ_ℓ(j, Y)` that did not yield a certified edge.
#[derive(Clone, Debug)]
pub struct Failure {
    pub source: usize,
    pub ell: u64,
    pub j2: Fe,
    pub error: EdgeError,
}

#[derive(Clone, Debug, Default)]
pub struct Stats {
    pub phi_ms: Vec<(u64, u128, usize)>,
    pub walk_ms: u128,
    pub audit_ms: u128,
    pub edges_beyond_cap: usize,
    /// `(ℓ, roots of Φ_ℓ(j, Y)) → count of expanded curves`.
    pub root_histogram: HashMap<(u64, usize), usize>,
}

pub struct Walk {
    pub field: Field,
    pub start: StartCurve,
    pub config: WalkConfig,
    pub class: ClassInfo,
    pub phis: Vec<ModPoly>,
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
    pub failures: Vec<Failure>,
    pub stats: Stats,
    /// [`traits::class_audits`] once run.
    pub class_audits: Option<V>,
}

type Expansion = Vec<(u64, usize, Vec<(Fe, Result<Isogeny, EdgeError>)>)>;

fn expand(f: &Field, phis: &[ModPoly], m: &Model, seed: u64) -> Expansion {
    let j = m.j(f).expect("nonsingular");
    phis.iter()
        .map(|phi| {
            let roots = poly::roots(f, &phi.at_x(f, &j), seed ^ phi.ell);
            let n = roots.len();
            let out = roots
                .into_iter()
                .map(|j2| (j2, kernel::elkies_isogeny(f, m, phi, &j2)))
                .collect();
            (phi.ell, n, out)
        })
        .collect()
}

impl Walk {
    pub fn new(start: StartCurve, config: WalkConfig) -> Result<Self, String> {
        let field = Field::new(&start.p).ok_or("p is not an odd prime below 2^256")?;
        for &l in &config.primes {
            if l < 3 || l % 2 == 0 || !is_probable_prime(&BigUint::from(l)) {
                return Err(format!(
                    "ℓ = {l}: the walk takes odd primes (ℓ = 2 needs rational 2-torsion; see the class record)"
                ));
            }
        }
        let root = Model {
            a: field.from_big(&start.a),
            b: field.from_big(&start.b),
        };
        let (gx, gy) = (field.from_big(&start.gx), field.from_big(&start.gy));
        if !root.on_curve(&field, &gx, &gy) {
            return Err("the generator is not on the curve".into());
        }
        if curve::scalar_mul(&field, &root, &gx, &gy, &start.subgroup_order).is_some() {
            return Err("[r]G is not the identity".into());
        }
        let j = root.j(&field).ok_or("singular curve")?;
        let class = ClassInfo::compute(&start, &config.primes);
        let mut stats = Stats::default();
        let mut phis = Vec::new();
        let built: Vec<Result<(ModPoly, u128), String>> = config
            .primes
            .par_iter()
            .map(|&l| {
                let t = Instant::now();
                modpoly::modular_polynomial(&field, l).map(|phi| (phi, t.elapsed().as_millis()))
            })
            .collect();
        for b in built {
            let (phi, ms) = b?;
            stats.phi_ms.push((phi.ell, ms, phi.checks_passed));
            phis.push(phi);
        }
        Ok(Walk {
            field,
            start,
            config,
            class,
            phis,
            nodes: vec![Node {
                j,
                model: root,
                depth: 0,
                parent_edge: None,
                expanded: false,
                generator: (gx, gy),
                audit: None,
                roots: Vec::new(),
                levels: Vec::new(),
            }],
            edges: Vec::new(),
            failures: Vec::new(),
            stats,
            class_audits: None,
        })
    }

    pub fn run(&mut self) {
        let t = Instant::now();
        let f = &self.field;
        let mut index: HashMap<Fe, usize> = HashMap::new();
        index.insert(self.nodes[0].j, 0);
        let mut next = 0usize;
        let chunk = rayon::current_num_threads().max(1) * 4;
        while next < self.nodes.len() && self.nodes.len() < self.config.max_curves {
            let depth = self.nodes[next].depth;
            if depth >= self.config.max_depth {
                break;
            }
            let end = (next + chunk).min(self.nodes.len()).min(
                // Stay within one depth so the order is breadth-first.
                (next..self.nodes.len())
                    .find(|&i| self.nodes[i].depth != depth)
                    .unwrap_or(self.nodes.len()),
            );
            let models: Vec<Model> = self.nodes[next..end].iter().map(|n| n.model).collect();
            let phis = &self.phis;
            let seed = self.config.seed;
            let results: Vec<Expansion> = models
                .par_iter()
                .map(|m| expand(f, phis, m, seed))
                .collect();
            for (offset, exp) in results.into_iter().enumerate() {
                let src = next + offset;
                if self.nodes.len() >= self.config.max_curves {
                    break;
                }
                self.nodes[src].expanded = true;
                for (ell, nroots, found) in exp {
                    *self.stats.root_histogram.entry((ell, nroots)).or_default() += 1;
                    self.nodes[src].roots.push((ell, nroots));
                    for (j2, res) in found {
                        let iso = match res {
                            Ok(iso) => iso,
                            Err(error) => {
                                self.failures.push(Failure {
                                    source: src,
                                    ell,
                                    j2,
                                    error,
                                });
                                continue;
                            }
                        };
                        let dst = match index.get(&j2) {
                            Some(&d) => d,
                            None => {
                                if self.nodes.len() >= self.config.max_curves {
                                    self.stats.edges_beyond_cap += 1;
                                    continue;
                                }
                                let Some((model, _)) = curve::canonical_model(f, &iso.codomain)
                                else {
                                    self.failures.push(Failure {
                                        source: src,
                                        ell,
                                        j2,
                                        error: EdgeError::SpecialJ,
                                    });
                                    continue;
                                };
                                let d = self.nodes.len();
                                let generator =
                                    curve::deterministic_generator(f, &model, &self.start.cofactor);
                                self.nodes.push(Node {
                                    j: j2,
                                    model,
                                    depth: self.nodes[src].depth + 1,
                                    parent_edge: Some(self.edges.len()),
                                    expanded: false,
                                    generator,
                                    audit: None,
                                    roots: Vec::new(),
                                    levels: Vec::new(),
                                });
                                index.insert(j2, d);
                                d
                            }
                        };
                        match kernel::link_to_target(f, &iso.codomain, &self.nodes[dst].model) {
                            Ok(u2) => self.edges.push(Edge {
                                ell,
                                dir: 'x',
                                source: src,
                                target: dst,
                                kernel: iso.kernel,
                                velu: iso.codomain,
                                u2,
                            }),
                            Err(error) => self.failures.push(Failure {
                                source: src,
                                ell,
                                j2,
                                error,
                            }),
                        }
                    }
                }
            }
            next = end;
        }
        self.resolve_levels();
        self.stats.walk_ms = t.elapsed().as_millis();
        let t = Instant::now();
        let f = &self.field;
        let n = self.start.order.clone();
        let prime = self.class.order_prime;
        let pts = self.config.audit_points;
        let audits: Vec<OrderAudit> = self
            .nodes
            .par_iter()
            .map(|node| curve::audit_order(f, &node.model, &n, prime, pts, 0))
            .collect();
        for (node, a) in self.nodes.iter_mut().zip(audits) {
            node.audit = Some(a);
        }
        self.stats.audit_ms = t.elapsed().as_millis();
    }
}

impl Walk {
    /// Prove each curve's level in each walked `ℓ`-volcano where the
    /// evidence allows, and orient the edges by the levels.
    ///
    /// The `ℓ`-volcano's depth is `v_ℓ(f_π) = ⌊v_ℓ(t² − 4p)/2⌋` for odd `ℓ`
    /// (`v_ℓ(D_K) ≤ 1`).  Depth 0: every curve is at level 0.  Depth 1: a
    /// curve with `ℓ + 1` rational `ℓ`-isogenies is on the surface, one with
    /// a single rational `ℓ`-isogeny is on the floor, and the one neighbour
    /// of a floor curve is on the surface.  Deeper volcanoes are left
    /// unresolved.  An edge is `d`, `u` or `h` when both levels are proved
    /// (level change `+1`, `−1`, `0`), else `x`.
    fn resolve_levels(&mut self) {
        let primes = self.config.primes.clone();
        for &l in &primes {
            let v = self.class.v_ell(l);
            let depth = self.class.depth(l);
            let mut lv: Vec<Option<(u32, String)>> = vec![None; self.nodes.len()];
            if depth == 0 {
                for x in lv.iter_mut() {
                    *x = Some((
                        0,
                        format!("v_{l}(t^2-4p) = {v} <= 1: the {l}-volcano has depth 0"),
                    ));
                }
            } else if depth == 1 {
                for (i, n) in self.nodes.iter().enumerate() {
                    if let Some(&(_, r)) = n.roots.iter().find(|x| x.0 == l) {
                        if r as u64 == l + 1 {
                            lv[i] = Some((0, format!("v_{l}(t^2-4p) = {v}: depth 1; {r} rational {l}-isogenies, so surface")));
                        } else if r == 1 {
                            lv[i] = Some((1, format!("v_{l}(t^2-4p) = {v}: depth 1; one rational {l}-isogeny, so floor")));
                        }
                    }
                }
                for e in self.edges.iter().filter(|e| e.ell == l) {
                    for (a, b) in [(e.source, e.target), (e.target, e.source)] {
                        if lv[a].as_ref().map(|x| x.0) == Some(1) && lv[b].is_none() {
                            lv[b] = Some((
                                0,
                                format!("depth 1: {l}-isogenous to a floor curve, so surface"),
                            ));
                        }
                    }
                }
            }
            for (i, x) in lv.iter().enumerate() {
                if let Some((level, proof)) = x {
                    self.nodes[i].levels.push((l, *level, proof.clone()));
                }
            }
            for e in self.edges.iter_mut().filter(|e| e.ell == l) {
                e.dir = match (lv[e.source].as_ref(), lv[e.target].as_ref()) {
                    (Some(a), Some(b)) => match b.0 as i64 - a.0 as i64 {
                        1 => 'd',
                        -1 => 'u',
                        0 => 'h',
                        _ => 'x',
                    },
                    _ => 'x',
                };
            }
        }
    }

    /// Run the class audits ([`traits::class_audits`]) on the root, and
    /// re-run them on the first `sample` walked curves to check that their
    /// verdicts are those of the root.
    pub fn run_class_audits(&mut self, sample: usize) {
        let f = &self.field;
        let picked: Vec<_> = self
            .nodes
            .iter()
            .skip(1)
            .take(sample)
            .map(|n| {
                (
                    f.to_big(&n.model.a),
                    f.to_big(&n.model.b),
                    f.to_big(&n.generator.0),
                    f.to_big(&n.generator.1),
                )
            })
            .collect();
        self.class_audits = Some(traits::class_audits(&self.start, &picked));
    }

    /// `<EC1>V<ℓ>L<level>…` over the proved levels, `ℓ` ascending.
    pub fn position_alias(&self, i: usize, ec1: &str) -> String {
        let mut lv = self.nodes[i].levels.clone();
        lv.sort();
        let segs: String = lv.iter().map(|(l, v, _)| format!("V{l}L{v}")).collect();
        format!("{ec1}{segs}")
    }
}

/// Identities and display names of one node, computed once for output.
pub struct NodeIds {
    pub slug: String,
    pub icv1: String,
    pub ec1: String,
    pub uid: String,
    pub field: V,
    pub curve: V,
    pub is_root: bool,
}

impl Walk {
    fn fe(&self, x: &Fe) -> BigUint {
        self.field.to_big(x)
    }

    pub fn node_ids(&self, i: usize) -> NodeIds {
        let node = &self.nodes[i];
        let (a, b) = (self.fe(&node.model.a), self.fe(&node.model.b));
        let id = curve_id::prime(&self.start.p, &a, &b, &self.start.order)
            .expect("nonsingular curve with order in the Hasse interval");
        let (gx, gy) = (self.fe(&node.generator.0), self.fe(&node.generator.1));
        let (field, curve) = record::ec1_records(
            &self.start.p,
            &a,
            &b,
            &self.start.subgroup_order,
            &self.start.cofactor,
            (&gx, &gy),
        );
        let tag = if i == 0 {
            self.start.tag.as_str()
        } else {
            "fp"
        };
        let (ec1, uid) = record::ec1_identity(&self.start.p, &field, &curve, tag);
        NodeIds {
            slug: id.slug,
            icv1: id.icv1,
            ec1,
            uid,
            field,
            curve,
            is_root: i == 0,
        }
    }

    fn poly_digits(&self, h: &Poly) -> V {
        V::Seq(h.iter().map(|c| V::s(self.fe(c).to_string())).collect())
    }

    /// The canonical record an edge's id and the routes through it hash.
    pub fn edge_record(&self, e: &Edge, ids: &[NodeIds]) -> V {
        V::map(vec![
            ("degree", V::int(e.ell)),
            ("direction", V::s(e.dir.to_string())),
            ("separable", V::Bool(true)),
            ("source_curve_uid", V::s(ids[e.source].uid.clone())),
            ("target_curve_uid", V::s(ids[e.target].uid.clone())),
            ("kernel_polynomial", self.poly_digits(&e.kernel)),
            (
                "velu_codomain_a4_a6",
                V::Seq(vec![
                    V::s(self.fe(&e.velu.a).to_string()),
                    V::s(self.fe(&e.velu.b).to_string()),
                ]),
            ),
            ("isomorphism_u2", V::s(self.fe(&e.u2).to_string())),
        ])
    }

    pub fn edge_id(&self, e: &Edge, ids: &[NodeIds]) -> String {
        let sha = record::sha256_hex(&self.edge_record(e, ids).canonical_json());
        format!("e{}{}1_{}", e.ell, e.dir, &sha[..12])
    }

    /// `IW1…` for an ordered list of edges (`VOLCANO_NAMING.md`).
    pub fn route_id(&self, path: &[usize], ids: &[NodeIds]) -> String {
        let edges: Vec<&Edge> = path.iter().map(|&i| &self.edges[i]).collect();
        let mut segs = String::new();
        let mut k = 0;
        while k < edges.len() {
            let mut run = 1;
            while k + run < edges.len()
                && edges[k + run].ell == edges[k].ell
                && edges[k + run].dir == edges[k].dir
            {
                run += 1;
            }
            segs.push_str(&format!("E{}{}{}", edges[k].ell, edges[k].dir, run));
            k += run;
        }
        let rec = V::map(vec![
            ("source_curve_uid", V::s(ids[edges[0].source].uid.clone())),
            (
                "target_curve_uid",
                V::s(ids[edges[edges.len() - 1].target].uid.clone()),
            ),
            (
                "intermediate_curve_uids",
                V::Seq(
                    edges[..edges.len() - 1]
                        .iter()
                        .map(|e| V::s(ids[e.target].uid.clone()))
                        .collect(),
                ),
            ),
            (
                "edges",
                V::Seq(edges.iter().map(|e| self.edge_record(e, ids)).collect()),
            ),
        ]);
        let sha = record::sha256_hex(&rec.canonical_json());
        format!("IW1{segs}h{}", &sha[..12])
    }

    /// BFS-tree edges from the root to node `i`, in walk order.
    pub fn root_path(&self, i: usize) -> Vec<usize> {
        let mut path = Vec::new();
        let mut v = i;
        while let Some(e) = self.nodes[v].parent_edge {
            path.push(e);
            v = self.edges[e].source;
        }
        path.reverse();
        path
    }

    fn equation(&self, m: &Model) -> String {
        format!("y^2=x^3+{}*x+{}", self.fe(&m.a), self.fe(&m.b))
    }

    pub fn curve_ref(&self, i: usize) -> String {
        format!("icwalk/{}/n{:06}", self.start.tag, i)
    }

    /// The curves in the `docs/curves/ic/curves.yaml` format.
    pub fn curves_yaml(&self, ids: &[NodeIds], routes: &RouteIndex) -> String {
        let detectors = traits::default_detectors();
        let detectors = &detectors;
        let f = &self.field;
        let class = &self.class;
        let mut curves = Vec::new();
        for (i, node) in self.nodes.iter().enumerate() {
            let id = &ids[i];
            let (a, b) = (self.fe(&node.model.a), self.fe(&node.model.b));
            let order_status = match node.audit {
                Some(OrderAudit::ProvedPrime) => "proved",
                Some(OrderAudit::Consistent) => "unproved_in_this_registry",
                _ => "unknown",
            };
            let end_disc = if class.disc_fundamental_if_cofactor_prime {
                V::map(vec![
                    ("value", V::bigi(&class.disc)),
                    ("status", V::s("unproved_in_this_registry")),
                    (
                        "proof_ref",
                        V::s("walk.json#class: D_pi is squarefree and 1 mod 4 if its trial-division cofactor is prime; the cofactor passes Miller-Rabin only"),
                    ),
                ])
            } else {
                V::map(vec![("value", V::Null), ("status", V::s("unknown"))])
            };
            let levels: Vec<(String, V)> = node
                .levels
                .iter()
                .map(|(l, level, proof)| {
                    (
                        l.to_string(),
                        V::map(vec![
                            ("level", V::int(level)),
                            ("status", V::s("proved")),
                            ("proof_ref", V::s(proof.clone())),
                        ]),
                    )
                })
                .collect();
            let degree_counts: Vec<(String, V)> = self
                .config
                .primes
                .iter()
                .map(|&l| {
                    let c = routes.out_by_ell.get(&(i, l)).copied().unwrap_or(0);
                    (l.to_string(), V::int(c))
                })
                .collect();
            let icv1_status = if id.is_root {
                "registered_representation"
            } else {
                "unresolved"
            };
            let label = if id.is_root {
                format!("{} (walk root)", self.start.name)
            } else {
                format!("{}-isogenous curve, walk node {i}", self.start.name)
            };
            let mut traits: Vec<(String, V)> = vec![
                (
                    "ordinary",
                    V::map(vec![
                        ("value", V::Bool(true)),
                        ("status", V::s("derived_from_model")),
                    ]),
                ),
                (
                    "j_invariant",
                    V::map(vec![
                        ("value", V::big(&self.fe(&node.j))),
                        ("status", V::s("derived_from_model")),
                    ]),
                ),
                ("endomorphism_discriminant", end_disc),
                (
                    "volcano_component",
                    V::map(vec![("value", V::Null), ("status", V::s("unknown"))]),
                ),
                (
                    "volcano_total_depth",
                    V::map(vec![("value", V::Null), ("status", V::s("unmeasured"))]),
                ),
                (
                    "group_order",
                    V::map(vec![
                        ("value", V::big(&self.start.order)),
                        ("status", V::s(order_status)),
                    ]),
                ),
                (
                    "trace",
                    V::map(vec![
                        ("value", V::bigi(&class.trace)),
                        ("status", V::s("derived_from_model")),
                    ]),
                ),
            ]
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect();
            let ctx = traits::CurveCtx {
                field: f,
                model: &node.model,
                generator: node.generator,
                start: &self.start,
            };
            for d in detectors {
                traits.push((
                    d.name().to_string(),
                    V::map(vec![
                        ("value", d.detect(&ctx)),
                        ("status", V::s(d.status())),
                    ]),
                ));
            }
            curves.push((
                id.slug.clone(),
                V::map(vec![
                    ("label", V::s(label)),
                    (
                        "curve_tag",
                        V::s(if id.is_root {
                            self.start.tag.clone()
                        } else {
                            "fp".into()
                        }),
                    ),
                    ("curve_ref", V::s(self.curve_ref(i))),
                    ("curve_id", V::s(id.ec1.clone())),
                    ("curve_uid", V::s(id.uid.clone())),
                    ("position_alias", V::s(self.position_alias(i, &id.ec1))),
                    (
                        "icv1_identity",
                        V::map(vec![
                            ("slug", V::s(id.slug.clone())),
                            ("full", V::s(id.icv1.clone())),
                            ("status", V::s(icv1_status)),
                        ]),
                    ),
                    (
                        "icv1_registration",
                        V::s(if id.is_root {
                            "registered in docs/curves/registry.json"
                        } else {
                            "computed by curve_id::prime; not registered"
                        }),
                    ),
                    ("trait_status", V::Map(traits)),
                    ("factor_base_refs", V::Seq(vec![])),
                    ("factor_base_link_status", V::s("not_reconciled")),
                    ("field", id.field.clone()),
                    ("curve", id.curve.clone()),
                    ("equation", V::s(self.equation(&node.model))),
                    ("a", V::big(&a)),
                    ("b", V::big(&b)),
                    (
                        "endomorphism",
                        V::map(vec![
                            ("endomorphism_order_conductor", V::Null),
                            ("frobenius_order_conductor", V::Null),
                            ("frobenius_discriminant", V::bigi(&class.disc)),
                            (
                                "volcano_levels",
                                if levels.is_empty() {
                                    V::Null
                                } else {
                                    V::Map(levels)
                                },
                            ),
                            ("status", V::s("partial")),
                        ]),
                    ),
                    (
                        "isogeny_routes",
                        V::map(vec![
                            (
                                "incoming",
                                V::Seq(
                                    routes
                                        .incoming
                                        .get(&i)
                                        .cloned()
                                        .unwrap_or_default()
                                        .into_iter()
                                        .map(V::s)
                                        .collect(),
                                ),
                            ),
                            (
                                "outgoing",
                                V::Seq(
                                    routes
                                        .outgoing
                                        .get(&i)
                                        .cloned()
                                        .unwrap_or_default()
                                        .into_iter()
                                        .map(V::s)
                                        .collect(),
                                ),
                            ),
                        ]),
                    ),
                    ("root_route", V::opt(routes.root.get(&i).cloned().map(V::s))),
                    (
                        "walk",
                        V::map(vec![
                            ("node", V::int(i)),
                            ("depth", V::int(node.depth)),
                            ("expanded", V::Bool(node.expanded)),
                            ("rational_neighbours_by_ell", V::Map(degree_counts)),
                        ]),
                    ),
                    ("evidence_ref", V::s("isogeny_routes.json")),
                ]),
            ));
        }
        let doc = V::map(vec![
            ("schema_version", V::int(1)),
            (
                "identity_rule",
                V::s("sha256_sorted_key_compact_utf8_json_of_field_and_curve"),
            ),
            ("curves", V::Map(curves)),
        ]);
        format!(
            "# Isogeny walk from {} written by `isogeny_walk` (src/bin/isogeny_walk.rs).\n\
             # Format: docs/curves/ic/curves.yaml and its schema.  `field` and `curve` are\n\
             # exactly the EC1 preimage, so curve_uid is recomputable from this file.\n\
             # Models: {}\n# Generators: {}\n\
             # ICV1 slugs of walked curves are computed, not registered (AGENTS.md §11:\n\
             # register a slug before citing it in prose).\n{}",
            self.start.name,
            curve::CANONICAL_RULE,
            curve::GENERATOR_RULE,
            doc.yaml()
        )
    }

    /// `isogeny_routes.json` in cryptanalysis's catalog format.
    pub fn routes_json(&self, ids: &[NodeIds], routes: &RouteIndex) -> V {
        let nodes = self
            .nodes
            .iter()
            .enumerate()
            .map(|(i, n)| {
                let levels: Vec<(String, V)> = n
                    .levels
                    .iter()
                    .map(|(l, level, _)| (l.to_string(), V::int(level)))
                    .collect();
                V::map(vec![
                    ("ref", V::s(self.curve_ref(i))),
                    ("icv1_slug", V::s(ids[i].slug.clone())),
                    ("curve_id", V::s(ids[i].ec1.clone())),
                    ("curve_uid", V::s(ids[i].uid.clone())),
                    ("field_characteristic", V::s(self.start.p.to_string())),
                    ("equation", V::s(self.equation(&n.model))),
                    (
                        "coefficients_a1_a2_a3_a4_a6",
                        V::Seq(vec![
                            V::int(0),
                            V::int(0),
                            V::int(0),
                            V::s(self.fe(&n.model.a).to_string()),
                            V::s(self.fe(&n.model.b).to_string()),
                        ]),
                    ),
                    ("j_invariant", V::s(self.fe(&n.j).to_string())),
                    (
                        "generator_G",
                        V::Seq(vec![
                            V::s(self.fe(&n.generator.0).to_string()),
                            V::s(self.fe(&n.generator.1).to_string()),
                        ]),
                    ),
                    ("group_order", V::s(self.start.order.to_string())),
                    (
                        "subgroup_order",
                        V::s(self.start.subgroup_order.to_string()),
                    ),
                    ("cofactor", V::s(self.start.cofactor.to_string())),
                    ("depth", V::int(n.depth)),
                    ("expanded", V::Bool(n.expanded)),
                    ("proved_volcano_levels", V::Map(levels)),
                    ("position_alias", V::s(self.position_alias(i, &ids[i].ec1))),
                ])
            })
            .collect();
        let edges = self
            .edges
            .iter()
            .map(|e| {
                let mut rec = match self.edge_record(e, ids) {
                    V::Map(m) => m,
                    _ => unreachable!(),
                };
                rec.insert(0, ("id".into(), V::s(self.edge_id(e, ids))));
                rec.insert(1, ("status".into(), V::s("verified")));
                rec.insert(2, ("source_curve_ref".into(), V::s(self.curve_ref(e.source))));
                rec.insert(3, ("target_curve_ref".into(), V::s(self.curve_ref(e.target))));
                rec.push((
                    "kernel_certificate".into(),
                    V::s("kernel::verify_kernel: squarefree, divides psi_ell, closed under [g], Velu codomain F_p-isomorphic to the target"),
                ));
                V::Map(rec)
            })
            .collect();
        V::map(vec![
            ("schema_version", V::int(1)),
            ("curve_nodes", V::Seq(nodes)),
            ("edges", V::Seq(edges)),
            ("routes", V::Seq(routes.records.clone())),
        ])
    }

    /// Run summary: configuration, class invariants, counts and timings.
    pub fn summary(&self, ids: &[NodeIds]) -> V {
        let mut hist: Vec<_> = self.stats.root_histogram.iter().collect();
        hist.sort();
        let failures: Vec<V> = self
            .failures
            .iter()
            .map(|fl| {
                V::map(vec![
                    ("source", V::s(self.curve_ref(fl.source))),
                    ("ell", V::int(fl.ell)),
                    ("j2", V::s(self.fe(&fl.j2).to_string())),
                    ("error", V::s(format!("{:?}", fl.error))),
                ])
            })
            .collect();
        let qr: Vec<u32> = self
            .nodes
            .iter()
            .map(|n| curve::qr_prefix_64(&self.field, &n.model))
            .collect();
        let a3 = self
            .nodes
            .iter()
            .filter(|n| !curve::a_minus_3_models(&self.field, &n.model).is_empty())
            .count();
        let audits = |want: OrderAudit| self.nodes.iter().filter(|n| n.audit == Some(want)).count();
        V::map(vec![
            ("schema", V::s("isogeny-walk-summary/v1")),
            ("tool", V::s("src/bin/isogeny_walk.rs")),
            ("root", V::map(vec![
                ("name", V::s(self.start.name.clone())),
                ("icv1_slug", V::s(ids[0].slug.clone())),
                ("curve_id", V::s(ids[0].ec1.clone())),
                ("curve_uid", V::s(ids[0].uid.clone())),
            ])),
            ("config", V::map(vec![
                ("primes", V::Seq(self.config.primes.iter().map(V::int).collect())),
                ("max_curves", V::int(self.config.max_curves)),
                ("max_depth", if self.config.max_depth == usize::MAX { V::Null } else { V::int(self.config.max_depth) }),
                ("seed", V::int(self.config.seed)),
                ("audit_points", V::int(self.config.audit_points)),
                ("skipped_primes", V::Seq(self.config.skipped.iter().map(|(l, why)| V::map(vec![
                    ("ell", V::int(l)),
                    ("reason", V::s(why.clone())),
                ])).collect())),
                ("threads", V::int(rayon::current_num_threads())),
            ])),
            ("rules", V::map(vec![
                ("canonical_model", V::s(curve::CANONICAL_RULE)),
                ("generator", V::s(curve::GENERATOR_RULE)),
                ("edge_direction", V::s("d/u/h when both endpoint levels are proved (depth 0: all level 0; depth 1: ell+1 rational ell-isogenies = surface, one = floor); x otherwise")),
            ])),
            ("class", self.class.record()),
            ("modular_polynomials", V::Seq(self.stats.phi_ms.iter().map(|(l, ms, checks)| V::map(vec![
                ("ell", V::int(l)),
                ("build_ms", V::int(ms)),
                ("vanishing_checks_passed", V::int(checks)),
            ])).collect())),
            ("counts", V::map(vec![
                ("curves", V::int(self.nodes.len())),
                ("curves_expanded", V::int(self.nodes.iter().filter(|n| n.expanded).count())),
                ("max_depth_reached", V::int(self.nodes.iter().map(|n| n.depth).max().unwrap_or(0))),
                ("edges_verified", V::int(self.edges.len())),
                ("edges_beyond_cap", V::int(self.stats.edges_beyond_cap)),
                ("failures", V::int(self.failures.len())),
                ("order_proved_prime", V::int(audits(OrderAudit::ProvedPrime))),
                ("order_consistent_unproved", V::int(audits(OrderAudit::Consistent))),
                ("order_refuted", V::int(audits(OrderAudit::Refuted))),
                ("a_minus_3_models", V::int(a3)),
                ("qr_prefix_64_min", V::opt(qr.iter().min().map(V::int))),
                ("qr_prefix_64_max", V::opt(qr.iter().max().map(V::int))),
            ])),
            ("qr_prefix_64_histogram", {
                let mut h = std::collections::BTreeMap::new();
                for q in &qr {
                    *h.entry(*q).or_insert(0usize) += 1;
                }
                V::Seq(h.into_iter().map(|(q, c)| V::Seq(vec![V::int(q), V::int(c)])).collect())
            }),
            ("proved_levels", V::Seq(self.config.primes.iter().map(|&l| {
                let mut by = std::collections::BTreeMap::new();
                let mut unresolved = 0usize;
                for n in &self.nodes {
                    match n.levels.iter().find(|x| x.0 == l) {
                        Some((_, level, _)) => *by.entry(*level).or_insert(0usize) += 1,
                        None => unresolved += 1,
                    }
                }
                V::map(vec![
                    ("ell", V::int(l)),
                    ("volcano_depth", V::int(self.class.depth(l))),
                    ("curves_by_level", V::Seq(by.into_iter().map(|(lv, c)| V::Seq(vec![V::int(lv), V::int(c)])).collect())),
                    ("curves_unresolved", V::int(unresolved)),
                ])
            }).collect())),
            ("roots_per_expanded_curve", V::Seq(hist.iter().map(|((l, r), c)| V::map(vec![
                ("ell", V::int(l)),
                ("roots", V::int(r)),
                ("curves", V::int(c)),
            ])).collect())),
            ("failures", V::Seq(failures)),
            ("class_audits", self.class_audits.clone().unwrap_or(V::Null)),
            ("timing_practicality_note", V::map(vec![
                ("walk_ms", V::int(self.stats.walk_ms)),
                ("audit_ms", V::int(self.stats.audit_ms)),
            ])),
        ])
    }
}

/// Route ids per node, built once from the edges.
pub struct RouteIndex {
    pub incoming: HashMap<usize, Vec<String>>,
    pub outgoing: HashMap<usize, Vec<String>>,
    pub out_by_ell: HashMap<(usize, u64), usize>,
    pub root: HashMap<usize, String>,
    pub records: Vec<V>,
}

impl RouteIndex {
    pub fn build(w: &Walk, ids: &[NodeIds]) -> Self {
        let mut idx = RouteIndex {
            incoming: HashMap::new(),
            outgoing: HashMap::new(),
            out_by_ell: HashMap::new(),
            root: HashMap::new(),
            records: Vec::new(),
        };
        let mut seen = std::collections::HashSet::new();
        let mut push = |idx: &mut RouteIndex, path: &[usize], rid: &str| {
            if seen.insert(rid.to_string()) {
                let first = &w.edges[path[0]];
                let last = &w.edges[path[path.len() - 1]];
                idx.records.push(V::map(vec![
                    ("id", V::s(rid)),
                    ("status", V::s("verified")),
                    ("source_curve_ref", V::s(w.curve_ref(first.source))),
                    ("target_curve_ref", V::s(w.curve_ref(last.target))),
                    (
                        "edge_ids",
                        V::Seq(
                            path.iter()
                                .map(|&e| V::s(w.edge_id(&w.edges[e], ids)))
                                .collect(),
                        ),
                    ),
                    ("search_prime", V::Null),
                ]));
            }
        };
        for (k, e) in w.edges.iter().enumerate() {
            let rid = w.route_id(&[k], ids);
            idx.outgoing.entry(e.source).or_default().push(rid.clone());
            idx.incoming.entry(e.target).or_default().push(rid.clone());
            *idx.out_by_ell.entry((e.source, e.ell)).or_default() += 1;
            push(&mut idx, &[k], &rid);
        }
        for i in 1..w.nodes.len() {
            let path = w.root_path(i);
            if path.is_empty() {
                continue;
            }
            let rid = w.route_id(&path, ids);
            idx.root.insert(i, rid.clone());
            push(&mut idx, &path, &rid);
        }
        idx
    }
}

/// Parse a decimal string field.
fn dec(v: &serde_json::Value, what: &str) -> Result<BigUint, String> {
    v.as_str()
        .and_then(|s| BigUint::parse_bytes(s.as_bytes(), 10))
        .ok_or_else(|| format!("{what}: expected a decimal string"))
}

/// Re-verify a walk's `isogeny_routes.json` from its records alone.
///
/// - Every node: its model, ICV1 slug, EC1 identity, generator and an order
///   audit; depth-0 level claims against `t² − 4p`.
/// - Every edge: its kernel certificate ([`kernel::verify_kernel`]), Vélu
///   codomain, isomorphism to the target, its direction against the
///   endpoints' recorded levels, and its id against the hash of its record.
/// - Every route: its `IW1` id against the hash of its ordered edges.
///
/// Shares no code with the construction beyond the field, the polynomial
/// arithmetic, the identities and the certificate; it never touches `Φ_ℓ`.
/// Depth-1 level claims rest on root counts the walk observed and are
/// checked here only for consistency with the edge directions.
pub fn verify_routes(
    routes: &serde_json::Value,
    start: &StartCurve,
    audit_points: usize,
) -> Result<(usize, usize), String> {
    let f = Field::new(&start.p).ok_or("bad field")?;
    let nodes = routes["curve_nodes"].as_array().ok_or("curve_nodes")?;
    let n = &start.order;
    let n_prime = is_probable_prime(n);
    let primes: Vec<u64> = nodes
        .iter()
        .flat_map(|nd| {
            nd["proved_volcano_levels"]
                .as_object()
                .map(|m| {
                    m.keys()
                        .filter_map(|k| k.parse().ok())
                        .collect::<Vec<u64>>()
                })
                .unwrap_or_default()
        })
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    let class = ClassInfo::compute(start, &primes);
    type Checked = (String, Model, String, HashMap<u64, u64>);
    let checked: Vec<Result<Checked, String>> = nodes
        .par_iter()
        .enumerate()
        .map(|(i, node)| {
            let r = node["ref"].as_str().ok_or("ref")?.to_string();
            let co = node["coefficients_a1_a2_a3_a4_a6"]
                .as_array()
                .ok_or("coefficients")?;
            let (a, b) = (dec(&co[3], "a4")?, dec(&co[4], "a6")?);
            let m = Model {
                a: f.from_big(&a),
                b: f.from_big(&b),
            };
            if i == 0 && (a != &start.a % &start.p || b != &start.b % &start.p) {
                return Err("the root is not the registered model".into());
            }
            let id = curve_id::prime(&start.p, &a, &b, n).ok_or("identity")?;
            if node["icv1_slug"].as_str() != Some(id.slug.as_str()) {
                return Err(format!("{r}: slug differs from the recomputed {}", id.slug));
            }
            let g = node["generator_G"].as_array().ok_or("generator")?;
            let (gx, gy) = (dec(&g[0], "gx")?, dec(&g[1], "gy")?);
            let (field, curve) = record::ec1_records(
                &start.p,
                &a,
                &b,
                &start.subgroup_order,
                &start.cofactor,
                (&gx, &gy),
            );
            let tag = if i == 0 { start.tag.as_str() } else { "fp" };
            let (ec1, uid) = record::ec1_identity(&start.p, &field, &curve, tag);
            if node["curve_uid"].as_str() != Some(uid.as_str())
                || node["curve_id"].as_str() != Some(ec1.as_str())
            {
                return Err(format!(
                    "{r}: EC1 identity differs from the recomputed {ec1}"
                ));
            }
            let (fx, fy) = (f.from_big(&gx), f.from_big(&gy));
            if !m.on_curve(&f, &fx, &fy)
                || curve::scalar_mul(&f, &m, &fx, &fy, &start.subgroup_order).is_some()
            {
                return Err(format!("{r}: generator is not a point of order dividing r"));
            }
            if curve::audit_order(&f, &m, n, n_prime, audit_points, 7) == OrderAudit::Refuted {
                return Err(format!("{r}: a point has [n]P != O"));
            }
            let mut levels = HashMap::new();
            if let Some(map) = node["proved_volcano_levels"].as_object() {
                for (k, v) in map {
                    let l: u64 = k.parse().map_err(|_| format!("{r}: level key {k}"))?;
                    let lv = v.as_u64().ok_or(format!("{r}: level of {k}"))?;
                    if class.depth(l) == 0 && lv != 0 {
                        return Err(format!("{r}: level {lv} claimed in a depth-0 {l}-volcano"));
                    }
                    if lv > u64::from(class.depth(l)) {
                        return Err(format!("{r}: level {lv} below the {l}-volcano's floor"));
                    }
                    levels.insert(l, lv);
                }
            }
            Ok((r, m, uid, levels))
        })
        .collect();
    let mut models: HashMap<String, (Model, String, HashMap<u64, u64>)> = HashMap::new();
    for c in checked {
        let (r, m, uid, levels) = c?;
        models.insert(r, (m, uid, levels));
    }
    let edges = routes["edges"].as_array().ok_or("edges")?;
    let edge_records: Vec<Result<(String, u64, char, V), String>> = edges
        .par_iter()
        .map(|e| {
            let id = e["id"].as_str().unwrap_or("?").to_string();
            let ell = e["degree"].as_u64().ok_or("degree")?;
            let dir = e["direction"]
                .as_str()
                .and_then(|d| d.chars().next())
                .ok_or("direction")?;
            let (src, su, sl) = models
                .get(e["source_curve_ref"].as_str().unwrap_or(""))
                .ok_or("source")?;
            let (dst, du, dl) = models
                .get(e["target_curve_ref"].as_str().unwrap_or(""))
                .ok_or("target")?;
            if e["source_curve_uid"].as_str() != Some(su.as_str())
                || e["target_curve_uid"].as_str() != Some(du.as_str())
            {
                return Err(format!("{id}: endpoint uid mismatch"));
            }
            let want_dir = match (sl.get(&ell), dl.get(&ell)) {
                (Some(a), Some(b)) => match *b as i64 - *a as i64 {
                    1 => 'd',
                    -1 => 'u',
                    0 => 'h',
                    _ => return Err(format!("{id}: endpoint levels differ by more than one")),
                },
                _ => 'x',
            };
            if dir != want_dir {
                return Err(format!(
                    "{id}: direction {dir} but the endpoint levels give {want_dir}"
                ));
            }
            let coeffs = e["kernel_polynomial"].as_array().ok_or("kernel")?;
            let h: Vec<Fe> = coeffs
                .iter()
                .map(|c| dec(c, "kernel coefficient").map(|v| f.from_big(&v)))
                .collect::<Result<_, _>>()?;
            let velu =
                kernel::verify_kernel(&f, src, ell, &h).map_err(|err| format!("{id}: {err:?}"))?;
            let u2 =
                kernel::link_to_target(&f, &velu, dst).map_err(|err| format!("{id}: {err:?}"))?;
            let rec = e["velu_codomain_a4_a6"].as_array().ok_or("velu")?;
            if f.from_big(&dec(&rec[0], "velu a")?) != velu.a
                || f.from_big(&dec(&rec[1], "velu b")?) != velu.b
                || f.from_big(&dec(&e["isomorphism_u2"], "u2")?) != u2
            {
                return Err(format!("{id}: recorded codomain or isomorphism differs"));
            }
            let fe = |x: &Fe| V::s(f.to_big(x).to_string());
            let record = V::map(vec![
                ("degree", V::int(ell)),
                ("direction", V::s(dir.to_string())),
                ("separable", V::Bool(true)),
                ("source_curve_uid", V::s(su.clone())),
                ("target_curve_uid", V::s(du.clone())),
                ("kernel_polynomial", V::Seq(h.iter().map(fe).collect())),
                (
                    "velu_codomain_a4_a6",
                    V::Seq(vec![fe(&velu.a), fe(&velu.b)]),
                ),
                ("isomorphism_u2", fe(&u2)),
            ]);
            let sha = record::sha256_hex(&record.canonical_json());
            if id != format!("e{ell}{dir}1_{}", &sha[..12]) {
                return Err(format!("{id}: id is not the hash of the edge record"));
            }
            Ok((id, ell, dir, record))
        })
        .collect();
    let mut by_id: HashMap<String, (u64, char, V)> = HashMap::new();
    for r in edge_records {
        let (id, ell, dir, rec) = r?;
        by_id.insert(id, (ell, dir, rec));
    }
    for route in routes["routes"].as_array().ok_or("routes")? {
        let rid = route["id"].as_str().unwrap_or("?");
        let path: Vec<&(u64, char, V)> = route["edge_ids"]
            .as_array()
            .ok_or("edge_ids")?
            .iter()
            .map(|e| {
                by_id
                    .get(e.as_str().unwrap_or(""))
                    .ok_or(format!("{rid}: unknown edge"))
            })
            .collect::<Result<_, _>>()?;
        if path.is_empty() {
            return Err(format!("{rid}: empty route"));
        }
        let uid = |rec: &V, key: &str| -> V {
            match rec {
                V::Map(m) => m
                    .iter()
                    .find(|(k, _)| k == key)
                    .map(|(_, v)| v.clone())
                    .unwrap_or(V::Null),
                _ => V::Null,
            }
        };
        for w in path.windows(2) {
            if uid(&w[0].2, "target_curve_uid") != uid(&w[1].2, "source_curve_uid") {
                return Err(format!("{rid}: edges do not chain"));
            }
        }
        let mut segs = String::new();
        let mut k = 0;
        while k < path.len() {
            let mut run = 1;
            while k + run < path.len()
                && path[k + run].0 == path[k].0
                && path[k + run].1 == path[k].1
            {
                run += 1;
            }
            segs.push_str(&format!("E{}{}{}", path[k].0, path[k].1, run));
            k += run;
        }
        let rec = V::map(vec![
            ("source_curve_uid", uid(&path[0].2, "source_curve_uid")),
            (
                "target_curve_uid",
                uid(&path[path.len() - 1].2, "target_curve_uid"),
            ),
            (
                "intermediate_curve_uids",
                V::Seq(
                    path[..path.len() - 1]
                        .iter()
                        .map(|e| uid(&e.2, "target_curve_uid"))
                        .collect(),
                ),
            ),
            ("edges", V::Seq(path.iter().map(|e| e.2.clone()).collect())),
        ]);
        let sha = record::sha256_hex(&rec.canonical_json());
        if rid != format!("IW1{segs}h{}", &sha[..12]) {
            return Err(format!("{rid}: id is not the hash of the route record"));
        }
    }
    Ok((nodes.len(), edges.len()))
}
