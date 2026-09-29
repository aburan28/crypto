//! # Endomorphism-invariant factor bases for index calculus.
//!
//! The binary Koblitz pipeline folds its factor base by the Frobenius:
//! a base closed under `π` needs one unknown per signed orbit, because
//! `log π(P) = λ_π · log P`.  The same identity holds for **any**
//! endomorphism `φ` of the curve — `φ` acts on the prime-order subgroup
//! as multiplication by an eigenvalue `λ_φ` — so a base closed under `φ`
//! folds the same way.  This module is that fold, written once over the
//! framework's [`CountedGroup`] so that every family below runs through
//! the same code and reports in the same unit:
//!
//! | type | example | `deg φ` | `ord_r(λ_φ)` | invariant base exists? |
//! |:--|:--|--:|:--|:--|
//! | A. automorphism | `j = 0`: `(x, y) ↦ (ζx, y)`; `j = 1728`: `(x, y) ↦ (−x, iy)`; negation | 1 | 3, 4, 2 | yes: any abscissa set closed under `x ↦ ζx` (`x ↦ −x`) |
//! | B. Frobenius-type | Koblitz `τ`; GLS `ψ = τ_u ∘ π_p ∘ τ_u⁻¹` on the twist over `F_{p²}`; `π_p` on a subfield curve over `F_{p³}` (`subfield_fp3`) | `p` | `n`; 4; 3 | yes: a `π`-stable subspace; the line `x ∈ u·s·F_p` (`gls_fp2`); the eigenline `x ∈ s^k·F_p` |
//! | A × B. composite | `ζ` and `ψ` on a `j = 0` twist over `F_{p²}`; `ι` and `ψ` on a `j = 1728` twist; `ζ` and `π` on a `j = 0` subfield curve | — | 12; 4; 3 | yes, and the fold is the **order of the subgroup of `(Z/rZ)^*` the eigenvalues generate**: 12 on the `j = 0` twist, but 4 on the `j = 1728` twist (`ι = ±ψ` on `⟨G⟩`) and 3 on the `j = 0` subfield curve (`ζ = π^{±1}` on `⟨G⟩`) |
//! | C. small-degree CM | `D = −7, −8`: degree 2 by Vélu; `1 + i` on `j = 1728`; `D = −11` and `√−3` on `j = 0`: degree 3 by Vélu | 2, 3 | large | **no** inside a prime-order subgroup (below) |
//!
//! ## The one fact the fold rests on, and its converse
//!
//! Let `F` be a factor base, `⟨G⟩` the subgroup of prime order `r`, `h`
//! the cofactor, and `φ` an endomorphism with `φ(P) = [λ]P` on `⟨G⟩`.
//! Every relation is written over `[h]P` (the relation loop multiplies
//! by `h`), and `[h]φ(P) = φ([h]P) = [λ][h]P`, so
//!
//! ```text
//!     log [h]φ(P) = λ · log [h]P          for every P ∈ F.
//! ```
//!
//! If `φ(F) ⊆ F` the unknowns of an orbit `{P, φP, φ²P, …}` are one
//! unknown times `1, λ, λ², …`, and the base has one column per orbit.
//! The orbit is finite exactly when `λ` has finite order modulo `r` on
//! the points of `F` that matter, i.e. when **`ord_r(λ)` is small**: an
//! automorphism of order `w` has `λ^w = 1`; the Frobenius of `E(F_{q^n})`
//! has `λ^n = 1`.  A degree-2 endomorphism such as the `D = −7` map has
//! `λ² − λ + 2 ≡ 0 (mod r)`, and `λ` then has order dividing `r − 1`
//! with no reason to be small; the orbit of any `[h]P ≠ O` under it has
//! `ord_r(λ)` points, so a `φ`-stable base inside `⟨G⟩` is either
//! nearly all of `⟨G⟩` or empty.  That is why type C gives no invariant
//! base, and the module *measures* it rather than asserting it:
//! [`eigenvalue_order`] reports `ord_r(λ)` and [`endomorphism_overlap`]
//! counts how many base points a type-C map keeps in the base.
//!
//! ## What is generic and what is not
//!
//! [`Endomorphism`] is the plug point: a name, a degree, the eigenvalue,
//! and the map on points.  [`fold_by_endomorphisms`] takes any list of
//! them, closes a seed set under the maps (or refuses a set that is not
//! closed), walks the orbits, checks that every two paths to one point
//! agree on the eigenvalue product, and returns a framework
//! [`FactorBase`] with one column per orbit.  It is the same function
//! for the prime-field automorphisms here, for the GLS map in
//! [`crate::cryptanalysis::gls_fp2`], and for whatever a later family
//! adds; the control (the same point set folded by negation alone) is
//! the same function with a shorter generator list.
//!
//! [`verify_endomorphism`] is what every constructor's output goes
//! through before it is used: the map sends points to points, is
//! additive, and acts as `[λ]` on `⟨G⟩`.  A constructor that gets a
//! formula wrong is caught here, not in a wrong logarithm three stages
//! later.
//!
//! ## What this module claims
//!
//! The fold is exact: a run over a folded base recovers the planted
//! logarithm or reports that it did not, and its column count is the
//! orbit count.  The instance generators produce curves whose group
//! order is certified on random points against the CM candidate orders
//! (Cornacchia), never assumed.
//!
//! ## What it does not claim
//!
//! Nothing about speed.  A fold by an automorphism group of order `w`
//! divides the columns by `w/2` against the negation fold, and the
//! relation count moves with its floor — the `j = 0` quotient of
//! `RESEARCH_GLV_INDEX_CALCULUS.md` measured exactly that (`3.0×`,
//! class *engineering*).  Rho takes the same `√w`.  The evaluation plan
//! (`RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`) states the boundaries.
//!
//! Companion to `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`.

use std::collections::VecDeque;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

use crate::cryptanalysis::ic_boundary::{
    CountedGroup, FactorBase, GroupOps, PrimeCurve, PrimeInstance, PrimePoint,
};
use crate::cryptanalysis::residual_walk::{inv_mod, is_prime_u64, pow_mod, sqrt_mod};

// ── Arithmetic modulo a word-sized modulus ─────────────────────────

/// `a·b mod m` without overflow, `m < 2^63`.
pub fn mulm(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}
pub fn addm(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + b as u128) % m as u128) as u64
}
pub fn subm(a: u64, b: u64, m: u64) -> u64 {
    addm(a % m, m - b % m, m)
}
fn negm(a: u64, m: u64) -> u64 {
    subm(0, a, m)
}

/// The multiplicative order of `lambda` modulo the prime `r`, by
/// factoring `r − 1` with trial division.  `r < 2^62`.
pub fn eigenvalue_order(lambda: u64, r: u64) -> u64 {
    assert!(r > 1 && !lambda.is_multiple_of(r), "order of a non-unit");
    let lambda = lambda % r;
    let mut order = r - 1;
    let mut n = r - 1;
    let mut f = 2u64;
    let mut primes = Vec::new();
    while f * f <= n {
        if n.is_multiple_of(f) {
            primes.push(f);
            while n.is_multiple_of(f) {
                n /= f;
            }
        }
        f += if f == 2 { 1 } else { 2 };
    }
    if n > 1 {
        primes.push(n);
    }
    for q in primes {
        while order.is_multiple_of(q) && pow_mod(lambda, order / q, r) == 1 {
            order /= q;
        }
    }
    order
}

/// Trial-division factorisation, `n < 2^62`, small enough for the
/// instance generators (which need the largest prime factor of a group
/// order of at most about `2^40`).
pub fn factor_u64(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut f = 2u64;
    while f * f <= n {
        if n.is_multiple_of(f) {
            let mut e = 0;
            while n.is_multiple_of(f) {
                n /= f;
                e += 1;
            }
            out.push((f, e));
        }
        f += if f == 2 { 1 } else { 2 };
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

// ── The plug point ─────────────────────────────────────────────────

/// An endomorphism of the curve as the fold sees it: the map on points
/// and the scalar it acts as on the prime-order subgroup.
///
/// `eigenvalue` must satisfy `apply(G) = [eigenvalue]G`; the fold also
/// uses it on points outside `⟨G⟩`, which is sound because relations
/// are written over `[h]P` and `φ` commutes with `[h]`.
/// [`verify_endomorphism`] checks the contract on random points.
pub trait Endomorphism<G: CountedGroup> {
    fn name(&self) -> String;
    /// The degree of the map: 1 for an automorphism, `p` for a
    /// Frobenius-type map, 2 or 3 for the small CM endomorphisms.
    fn degree(&self) -> u64;
    /// `λ` modulo `r` with `φ(P) = [λ]P` for `P ∈ ⟨G⟩`.
    fn eigenvalue(&self) -> u64;
    fn apply(&self, g: &G, p: G::Elt) -> G::Elt;
}

/// `P ↦ −P`, eigenvalue `r − 1`.  Every base is folded by it; it is a
/// generator like any other so that the control base (negation alone)
/// and the folded base (negation plus the family's maps) are the same
/// call with different generator lists.
#[derive(Clone, Copy, Debug)]
pub struct Negation {
    pub r: u64,
}

impl<G: CountedGroup> Endomorphism<G> for Negation {
    fn name(&self) -> String {
        "negation".into()
    }
    fn degree(&self) -> u64 {
        1
    }
    fn eigenvalue(&self) -> u64 {
        self.r - 1
    }
    fn apply(&self, g: &G, p: G::Elt) -> G::Elt {
        g.neg(p)
    }
}

/// What [`verify_endomorphism`] checked and found.
#[derive(Clone, Debug, Serialize)]
pub struct EndomorphismCheck {
    pub name: String,
    pub degree: u64,
    pub eigenvalue: u64,
    /// `ord_r(λ)`: the orbit length the fold would see on `⟨G⟩`.
    pub eigenvalue_order: u64,
    pub samples: u64,
    /// `φ([k]G) = [λ][k]G` on every sample.
    pub eigenvalue_holds: bool,
    /// `φ(P + Q) = φ(P) + φ(Q)` on every sample pair.
    pub additive: bool,
}

/// Check that `endo` is what it says: acts as `[λ]` on `samples` random
/// multiples of the generator and is additive on random pairs.  Any
/// failure is an error naming the check; a constructor's output goes
/// through this before it is used.
pub fn verify_endomorphism<G: CountedGroup>(
    g: &G,
    generator: G::Elt,
    r: u64,
    endo: &dyn Endomorphism<G>,
    samples: u64,
    seed: u64,
) -> Result<EndomorphismCheck, String> {
    let mut rng = StdRng::seed_from_u64(seed ^ 0x454E_444F);
    let mut ops = GroupOps::default();
    let lambda = endo.eigenvalue();
    if lambda >= r {
        return Err(format!(
            "{}: eigenvalue {lambda} is not reduced modulo r = {r}",
            endo.name()
        ));
    }
    let mut eigenvalue_holds = true;
    let mut additive = true;
    for _ in 0..samples {
        let k = rng.gen_range(1..r);
        let l = rng.gen_range(1..r);
        let p = g.mul(&mut ops, generator, k);
        let q = g.mul(&mut ops, generator, l);
        let phi_p = endo.apply(g, p);
        let expect = g.mul(&mut ops, p, lambda);
        if phi_p != expect {
            eigenvalue_holds = false;
        }
        let lhs = endo.apply(g, g.add(&mut ops, p, q));
        let rhs = g.add(&mut ops, phi_p, endo.apply(g, q));
        if lhs != rhs {
            additive = false;
        }
    }
    let check = EndomorphismCheck {
        name: endo.name(),
        degree: endo.degree(),
        eigenvalue: lambda,
        eigenvalue_order: eigenvalue_order(lambda, r),
        samples,
        eigenvalue_holds,
        additive,
    };
    if !eigenvalue_holds {
        return Err(format!(
            "{}: φ(P) ≠ [{lambda}]P on a random multiple of G",
            endo.name()
        ));
    }
    if !additive {
        return Err(format!("{}: φ is not additive", endo.name()));
    }
    Ok(check)
}

// ── The fold ───────────────────────────────────────────────────────

/// Whether a seed set that is not closed under the maps is completed
/// or refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Closure {
    /// Add every image until the set is closed.  This is how an
    /// invariant base is *built* from a seed (the smallest abscissae,
    /// say): the closure of `B` seeds under `⟨ζ⟩` is `3B` abscissae.
    Close,
    /// Refuse a set with an image outside it.  This is how a base that
    /// is invariant *by construction* (a `π`-stable subspace, the GLS
    /// line) is checked rather than trusted.
    Strict,
}

/// What the fold did, beside the base it built.
#[derive(Clone, Debug, Serialize)]
pub struct FoldReport {
    pub generators: Vec<String>,
    pub seed_points: usize,
    /// Points added by [`Closure::Close`]; zero under `Strict`.
    pub points_added: usize,
    pub points: usize,
    pub orbits: usize,
    pub smallest_orbit: usize,
    pub largest_orbit: usize,
    /// Points per orbit: the fold, measured.
    pub points_per_orbit: f64,
    /// Seed points dropped because `[h]P = O`: they carry no logarithm
    /// and would make the eigenvalue consistency check meaningless.
    pub torsion_dropped: usize,
    pub endomorphism_maps: u64,
}

/// Close `seed` under `gens` (or check that it is closed), walk the
/// orbits, and return a base with one column per orbit and coefficient
/// `Π λ` along the path from the orbit's representative.
///
/// Every point reached twice is checked to carry the same coefficient
/// by both paths; a disagreement means the generators do not act as
/// their eigenvalues say on this point — which is what happens when a
/// generator of infinite order on `⟨G⟩` is handed in, and also the
/// error a wrong eigenvalue produces — and is reported, not folded.  An
/// orbit that grows past `max_orbit` is the other symptom of an
/// infinite-order generator and is refused the same way.
///
/// `h` is the cofactor: seed points with `[h]P = O` are dropped (they
/// have no logarithm to fold).  The base is also required to be closed
/// under negation, which every generator list that includes
/// [`Negation`] guarantees.
#[allow(clippy::too_many_arguments)]
pub fn fold_by_endomorphisms<G: CountedGroup>(
    g: &G,
    r: u64,
    h: u64,
    seed: Vec<G::Elt>,
    gens: &[&dyn Endomorphism<G>],
    closure: Closure,
    max_orbit: usize,
    key_of: impl Fn(&G::Elt) -> u64,
    abscissa_of: impl Fn(&G::Elt) -> u64,
    description: String,
) -> Result<(FactorBase<G::Elt>, FoldReport), String> {
    use std::collections::HashMap;
    let mut ops = GroupOps::default();
    let mut points: Vec<G::Elt> = Vec::with_capacity(seed.len());
    let mut index: HashMap<u64, usize> = HashMap::with_capacity(seed.len() * 2);
    let mut torsion_dropped = 0usize;
    let seed_points = seed.len();
    for p in seed {
        if g.is_identity(&p) {
            continue;
        }
        if h > 1 && g.is_identity(&g.mul(&mut ops, p, h)) {
            torsion_dropped += 1;
            continue;
        }
        if let std::collections::hash_map::Entry::Vacant(e) = index.entry(key_of(&p)) {
            e.insert(points.len());
            points.push(p);
        }
    }
    let mut col_of: Vec<Option<usize>> = vec![None; points.len()];
    let mut coef_of: Vec<u64> = vec![0; points.len()];
    let mut orbits = 0usize;
    let mut smallest = usize::MAX;
    let mut largest = 0usize;
    let mut maps = 0u64;
    let mut points_added = 0usize;
    let mut i = 0usize;
    while i < points.len() {
        if col_of[i].is_some() {
            i += 1;
            continue;
        }
        let orbit = orbits;
        orbits += 1;
        col_of[i] = Some(orbit);
        coef_of[i] = 1;
        let mut size = 0usize;
        let mut queue = VecDeque::from([i]);
        while let Some(j) = queue.pop_front() {
            size += 1;
            if size > max_orbit {
                return Err(format!(
                    "orbit {orbit} exceeds {max_orbit} points: a generator has no finite \
                     order on this base (eigenvalue orders {:?}); it cannot fold",
                    gens.iter()
                        .map(|e| eigenvalue_order(e.eigenvalue(), r))
                        .collect::<Vec<_>>()
                ));
            }
            let pt = points[j];
            let c = coef_of[j];
            for e in gens {
                let img = e.apply(g, pt);
                maps += 1;
                if g.is_identity(&img) {
                    return Err(format!(
                        "{} maps a base point to the identity; a base point in its kernel \
                         cannot be folded",
                        e.name()
                    ));
                }
                let coef = mulm(c, e.eigenvalue(), r);
                let k = key_of(&img);
                match index.get(&k).copied() {
                    Some(t) => match col_of[t] {
                        Some(o) => {
                            if o != orbit || coef_of[t] != coef {
                                return Err(format!(
                                    "{} reaches point {t} with coefficient {coef} but it already \
                                     carries {} in orbit {o} (this orbit is {orbit}): the \
                                     generators do not act as their eigenvalues on this base",
                                    e.name(),
                                    coef_of[t]
                                ));
                            }
                        }
                        None => {
                            col_of[t] = Some(orbit);
                            coef_of[t] = coef;
                            queue.push_back(t);
                        }
                    },
                    None => match closure {
                        Closure::Strict => {
                            return Err(format!(
                                "{} maps a base point outside the base; the base is not \
                                 invariant",
                                e.name()
                            ))
                        }
                        Closure::Close => {
                            let t = points.len();
                            points.push(img);
                            index.insert(k, t);
                            col_of.push(Some(orbit));
                            coef_of.push(coef);
                            points_added += 1;
                            queue.push_back(t);
                        }
                    },
                }
            }
        }
        smallest = smallest.min(size);
        largest = largest.max(size);
        i += 1;
    }
    let col_of: Vec<usize> = col_of.into_iter().map(|c| c.expect("assigned")).collect();
    let mut fb = FactorBase::from_column_map(
        description,
        points,
        col_of,
        coef_of,
        orbits,
        &key_of,
        |p| key_of(&g.neg(*p)),
        &abscissa_of,
    )?;
    fb.cost.count("endomorphism_maps", maps);
    fb.cost.count("torsion_dropped", torsion_dropped as u64);
    fb.cost.group_ops.merge(ops);
    let report = FoldReport {
        generators: gens.iter().map(|e| e.name()).collect(),
        seed_points,
        points_added,
        points: fb.points.len(),
        orbits,
        smallest_orbit: if orbits == 0 { 0 } else { smallest },
        largest_orbit: largest,
        points_per_orbit: fb.points.len() as f64 / orbits.max(1) as f64,
        torsion_dropped,
        endomorphism_maps: maps,
    };
    Ok((fb, report))
}

/// How many points of a base an endomorphism keeps inside it, and the
/// orbit length it would need: the measurement behind "type C gives no
/// invariant base".  For an automorphism of the base every image is in
/// it; for a degree-2 CM map the count is the chance overlap of two
/// sets of size `|F|` in a group of size `#E`.
#[derive(Clone, Debug, Serialize)]
pub struct OverlapReport {
    pub name: String,
    pub degree: u64,
    pub eigenvalue: u64,
    pub eigenvalue_order: u64,
    pub base_points: usize,
    /// Base points `P` with `φ(P)` also in the base.
    pub images_in_base: usize,
    /// `images_in_base / base_points`.
    pub fraction: f64,
    /// `|F| / #E`: the fraction a random map would keep.
    pub chance_fraction: f64,
}

pub fn endomorphism_overlap<G: CountedGroup>(
    g: &G,
    fb: &FactorBase<G::Elt>,
    endo: &dyn Endomorphism<G>,
    r: u64,
    group_order: u64,
) -> OverlapReport {
    let mut inside = 0usize;
    for p in &fb.points {
        let img = endo.apply(g, *p);
        if !g.is_identity(&img) && fb.index_of_key(g.key(&img)).is_some() {
            inside += 1;
        }
    }
    OverlapReport {
        name: endo.name(),
        degree: endo.degree(),
        eigenvalue: endo.eigenvalue(),
        eigenvalue_order: eigenvalue_order(endo.eigenvalue(), r),
        base_points: fb.points.len(),
        images_in_base: inside,
        fraction: inside as f64 / fb.points.len().max(1) as f64,
        chance_fraction: fb.points.len() as f64 / group_order as f64,
    }
}

// ── Prime-field automorphisms (type A) ─────────────────────────────

/// `(x, y) ↦ (cx·x, cy·y)`: the shape of every automorphism of a short
/// Weierstrass curve over `F_p`.  `j = 0` has `cx = ζ` (order 3),
/// `j = 1728` has `(cx, cy) = (−1, i)` (order 4).
#[derive(Clone, Debug, Serialize)]
pub struct DiagonalAutomorphism {
    pub label: String,
    pub cx: u64,
    pub cy: u64,
    pub eigenvalue: u64,
    /// The order of the map on the curve.
    pub order: u32,
}

impl Endomorphism<PrimeCurve> for DiagonalAutomorphism {
    fn name(&self) -> String {
        self.label.clone()
    }
    fn degree(&self) -> u64 {
        1
    }
    fn eigenvalue(&self) -> u64 {
        self.eigenvalue
    }
    fn apply(&self, g: &PrimeCurve, p: PrimePoint) -> PrimePoint {
        if p.infinity {
            return p;
        }
        PrimePoint::affine(mulm(self.cx, p.x, g.p), mulm(self.cy, p.y, g.p))
    }
}

/// An element of exact order `n` in `F_p^*`, when `n | p − 1`.
fn root_of_unity(n: u64, p: u64) -> Option<u64> {
    if !(p - 1).is_multiple_of(n) {
        return None;
    }
    let cofactor = (p - 1) / n;
    // Check every proper divisor of n: exact order.
    let divisors: Vec<u64> = factor_u64(n).into_iter().map(|(q, _)| n / q).collect();
    (2..p)
        .map(|b| pow_mod(b, cofactor, p))
        .find(|&z| z != 1 && divisors.iter().all(|&d| pow_mod(z, d, p) != 1))
}

/// The eigenvalue on `⟨G⟩` of a diagonal map, among the roots of its
/// characteristic polynomial: `ψ(G) = [λ]G` decides between them.
fn eigenvalue_among(
    inst: &PrimeInstance,
    candidates: &[u64],
    image_of_g: PrimePoint,
) -> Option<u64> {
    let mut ops = GroupOps::default();
    let g = inst.generator_point();
    candidates
        .iter()
        .copied()
        .find(|&l| inst.curve.mul(&mut ops, g, l) == image_of_g)
}

/// `ψ(x, y) = (ζx, y)` on `y² = x³ + b`, `p ≡ 1 (mod 3)`: the order-3
/// automorphism, with `λ² + λ + 1 ≡ 0 (mod r)` decided on `G`.
pub fn j0_automorphism(inst: &PrimeInstance) -> Result<DiagonalAutomorphism, String> {
    let (p, r) = (inst.curve.p, inst.r);
    if inst.curve.a != 0 {
        return Err("not a j = 0 curve (a ≠ 0)".into());
    }
    let zeta = root_of_unity(3, p).ok_or("p ≢ 1 (mod 3): the curve is supersingular here")?;
    let s = sqrt_mod(r - 3, r).ok_or("−3 is not a square mod r; λ is undefined")?;
    let half = inv_mod(2, r);
    let cands = [
        mulm(subm(s, 1, r), half, r),
        mulm(subm(r - s, 1, r), half, r),
    ];
    let g = inst.generator_point();
    let psi_g = PrimePoint::affine(mulm(zeta, g.x, p), g.y);
    let eigenvalue = eigenvalue_among(inst, &cands, psi_g)
        .ok_or("neither root of λ² + λ + 1 is the eigenvalue of ψ on G")?;
    Ok(DiagonalAutomorphism {
        label: "zeta3".into(),
        cx: zeta,
        cy: 1,
        eigenvalue,
        order: 3,
    })
}

/// `ι(x, y) = (−x, iy)` on `y² = x³ + ax`, `p ≡ 1 (mod 4)`: the order-4
/// automorphism (`ι² = −1`), with `λ² ≡ −1 (mod r)` decided on `G`.
pub fn j1728_automorphism(inst: &PrimeInstance) -> Result<DiagonalAutomorphism, String> {
    let (p, r) = (inst.curve.p, inst.r);
    if inst.curve.b != 0 {
        return Err("not a j = 1728 curve (b ≠ 0)".into());
    }
    let i = root_of_unity(4, p).ok_or("p ≢ 1 (mod 4): the curve is supersingular here")?;
    let s = sqrt_mod(r - 1, r).ok_or("−1 is not a square mod r; λ is undefined")?;
    let g = inst.generator_point();
    let iota_g = PrimePoint::affine(negm(g.x, p), mulm(i, g.y, p));
    let eigenvalue = eigenvalue_among(inst, &[s, r - s], iota_g)
        .ok_or("neither square root of −1 is the eigenvalue of ι on G")?;
    Ok(DiagonalAutomorphism {
        label: "iota4".into(),
        cx: p - 1,
        cy: i,
        eigenvalue,
        order: 4,
    })
}

/// Which automorphism group a prime-field base is folded by.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum AutomorphismGroup {
    /// Whatever the curve has: `μ₆` on `j = 0`, `μ₄` on `j = 1728`,
    /// `{±1}` otherwise.
    Auto,
    /// `{±1}` only: the control on any curve.
    Negation,
    J0,
    J1728,
}

impl AutomorphismGroup {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "auto" => Ok(Self::Auto),
            "negation" => Ok(Self::Negation),
            "j0" => Ok(Self::J0),
            "j1728" => Ok(Self::J1728),
            other => Err(format!(
                "unknown automorphism group `{other}`; try auto, negation, j0 or j1728"
            )),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Negation => "negation",
            Self::J0 => "j0",
            Self::J1728 => "j1728",
        }
    }
}

/// The generators of the chosen automorphism group on this instance:
/// negation always, plus the curve's extra automorphism when asked for
/// and present.  `Auto` on a curve with neither is negation alone, and
/// the returned list says so.
pub fn automorphism_generators(
    inst: &PrimeInstance,
    group: AutomorphismGroup,
) -> Result<Vec<Box<dyn Endomorphism<PrimeCurve>>>, String> {
    let mut gens: Vec<Box<dyn Endomorphism<PrimeCurve>>> = vec![Box::new(Negation { r: inst.r })];
    let p = inst.curve.p;
    let want_j0 = match group {
        AutomorphismGroup::J0 => true,
        AutomorphismGroup::Auto => inst.curve.a == 0 && p % 3 == 1,
        _ => false,
    };
    let want_j1728 = match group {
        AutomorphismGroup::J1728 => true,
        AutomorphismGroup::Auto => inst.curve.b == 0 && p % 4 == 1,
        _ => false,
    };
    if want_j0 {
        gens.push(Box::new(j0_automorphism(inst)?));
    }
    if want_j1728 {
        gens.push(Box::new(j1728_automorphism(inst)?));
    }
    Ok(gens)
}

/// The `size` smallest abscissae carrying a point, both signs: the seed
/// every prime-field fold starts from (the same seed as
/// `prime-abscissa`, so the control and the folded base share it).
pub fn smallest_abscissa_points(inst: &PrimeInstance, size: usize) -> Vec<PrimePoint> {
    let curve = &inst.curve;
    let mut out = Vec::with_capacity(2 * size);
    let mut abscissae = 0usize;
    let mut x = 1u64;
    while abscissae < size && x < curve.p {
        let rhs = curve.rhs(x);
        if curve.legendre(rhs) >= 0 {
            if let Some(y) = curve.sqrt(rhs) {
                let p = PrimePoint::affine(x, y);
                out.push(p);
                let q = curve.neg(p);
                if q != p {
                    out.push(q);
                }
                abscissae += 1;
            }
        }
        x += 1;
    }
    out
}

/// A prime-field base: the closure of the `size` smallest abscissae
/// under `group`, folded by `group` (or, with `fold = false`, by
/// negation alone: the control with the same points).
pub fn glv_orbit_base(
    inst: &PrimeInstance,
    size: usize,
    group: AutomorphismGroup,
    fold: bool,
) -> Result<(FactorBase<PrimePoint>, FoldReport), String> {
    let start = std::time::Instant::now();
    let gens = automorphism_generators(inst, group)?;
    let refs: Vec<&dyn Endomorphism<PrimeCurve>> = gens.iter().map(|b| b.as_ref()).collect();
    let seed = smallest_abscissa_points(inst, size);
    let curve = &inst.curve;
    // Close under the whole group first, so that the control folds the
    // same point set.
    let (closed, _) = fold_by_endomorphisms(
        curve,
        inst.r,
        inst.cofactor,
        seed,
        &refs,
        Closure::Close,
        64,
        |p| curve.key(p),
        |p| p.x,
        String::new(),
    )?;
    let fold_refs: Vec<&dyn Endomorphism<PrimeCurve>> =
        if fold { refs.clone() } else { vec![refs[0]] };
    let names: Vec<String> = fold_refs.iter().map(|e| e.name()).collect();
    let description = format!(
        "closure of the {size} smallest abscissae under {}, folded by {}",
        refs.iter()
            .map(|e| e.name())
            .collect::<Vec<_>>()
            .join(" and "),
        names.join(" and ")
    );
    let (mut fb, report) = fold_by_endomorphisms(
        curve,
        inst.r,
        inst.cofactor,
        closed.points,
        &fold_refs,
        Closure::Strict,
        64,
        |p| curve.key(p),
        |p| p.x,
        description,
    )?;
    fb.cost.wall_ns = start.elapsed().as_nanos() as u64;
    fb.cost.count(
        "abscissae_scanned",
        fb.abscissa_list.last().copied().unwrap_or(0),
    );
    fb.cost.count("sqrt_solves", size as u64);
    Ok((fb, report))
}

// ── CM instance generation ─────────────────────────────────────────

/// The curve families the generator builds, by CM discriminant.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum CmFamily {
    /// `y² = x³ + b`, `D = −3`, `Aut = μ₆`.
    J0,
    /// `y² = x³ + ax`, `D = −4`, `Aut = μ₄`.
    J1728,
    /// `j = −3375`, `D = −7`: a degree-2 endomorphism, `φ² − φ + 2 = 0`.
    D7,
    /// `j = 8000`, `D = −8`: a degree-2 endomorphism, `φ² + 2 = 0`.
    D8,
    /// `j = −32768`, `D = −11`: a degree-3 endomorphism, `φ² − φ + 3 = 0`.
    D11,
    /// A random `(a, b)`: negation only (the control family).
    Generic,
}

impl CmFamily {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "j0" => Ok(Self::J0),
            "j1728" => Ok(Self::J1728),
            "d7" => Ok(Self::D7),
            "d8" => Ok(Self::D8),
            "d11" => Ok(Self::D11),
            "generic" => Ok(Self::Generic),
            other => Err(format!(
                "unknown family `{other}`; try j0, j1728, d7, d8, d11 or generic"
            )),
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::J0 => "j0",
            Self::J1728 => "j1728",
            Self::D7 => "d7",
            Self::D8 => "d8",
            Self::D11 => "d11",
            Self::Generic => "generic",
        }
    }
}

fn random_prime(bits: u32, rng: &mut StdRng) -> u64 {
    loop {
        let c = rng.gen_range((1u64 << (bits - 1))..(1u64 << bits)) | 1;
        if is_prime_u64(c) {
            return c;
        }
    }
}

fn isqrt(n: u64) -> u64 {
    let mut s = (n as f64).sqrt() as u64;
    while s * s > n {
        s -= 1;
    }
    while (s + 1) * (s + 1) <= n {
        s += 1;
    }
    s
}

/// `(t, v)` with `4p = t² + |D| v²`, `t, v > 0`, if any: the trace of
/// the CM curves of discriminant `D` up to sign and unit.  A scan over
/// `v`, fine for `p < 2^40`.
pub fn cornacchia_4p(p: u64, abs_d: u64) -> Option<(u64, u64)> {
    let four_p = 4 * p;
    let mut v = 1u64;
    while abs_d * v * v < four_p {
        let rest = four_p - abs_d * v * v;
        let t = isqrt(rest);
        if t * t == rest && t > 0 {
            return Some((t, v));
        }
        v += 1;
    }
    None
}

/// The candidate group orders `p + 1 − t'` of the curves with the
/// family's `j`-invariant over `F_p`, from Cornacchia's `(t, v)`.
pub fn cm_candidate_orders(family: CmFamily, p: u64) -> Option<Vec<u64>> {
    let n = |t: i128| -> u64 { (p as i128 + 1 - t) as u64 };
    match family {
        CmFamily::J0 => {
            // 4p = t² + 3v²; the six traces are ±t, ±(t + 3v)/2, ±(t − 3v)/2.
            let (t, v) = cornacchia_4p(p, 3)?;
            let (t, v) = (t as i128, v as i128);
            let mut out = Vec::new();
            for tr in [t, (t + 3 * v) / 2, (t - 3 * v) / 2] {
                if (t + 3 * v) % 2 != 0 {
                    continue;
                }
                out.push(n(tr));
                out.push(n(-tr));
            }
            Some(out)
        }
        CmFamily::J1728 => {
            // p = u² + v²; the four traces are ±2u, ±2v.
            let (t, v) = cornacchia_4p(p, 4)?;
            // 4p = t² + 4v² ⇒ t even, p = (t/2)² + v².
            let (u, v) = ((t / 2) as i128, v as i128);
            Some(vec![n(2 * u), n(-2 * u), n(2 * v), n(-2 * v)])
        }
        CmFamily::D7 => {
            let (t, _) = cornacchia_4p(p, 7)?;
            Some(vec![n(t as i128), n(-(t as i128))])
        }
        CmFamily::D8 => {
            let (t, _) = cornacchia_4p(p, 8)?;
            Some(vec![n(t as i128), n(-(t as i128))])
        }
        CmFamily::D11 => {
            let (t, _) = cornacchia_4p(p, 11)?;
            Some(vec![n(t as i128), n(-(t as i128))])
        }
        CmFamily::Generic => None,
    }
}

/// The `j`-invariant's Weierstrass model `y² = x³ + 3kx + 2k`,
/// `k = j / (1728 − j)`, for `j ∉ {0, 1728}`.
fn curve_with_j(j: i64, p: u64) -> Option<(u64, u64)> {
    let jm = j.rem_euclid(p as i64) as u64;
    let denom = subm(1728 % p, jm, p);
    if denom == 0 || jm == 0 {
        return None;
    }
    let k = mulm(jm, inv_mod(denom, p), p);
    Some((mulm(3, k, p), mulm(2, k, p)))
}

/// Which of the candidate orders annihilates `points`; `None` if none
/// or more than one does.
fn certified_order(curve: &PrimeCurve, candidates: &[u64], points: &[PrimePoint]) -> Option<u64> {
    let mut ops = GroupOps::default();
    let mut found = None;
    for &n in candidates {
        if points.iter().all(|&pt| curve.mul(&mut ops, pt, n).infinity) {
            if found.is_some() && found != Some(n) {
                return None;
            }
            found = Some(n);
        }
    }
    found
}

fn random_point(curve: &PrimeCurve, rng: &mut StdRng) -> PrimePoint {
    loop {
        let x = rng.gen_range(1..curve.p);
        let rhs = curve.rhs(x);
        if curve.legendre(rhs) >= 0 {
            if let Some(y) = curve.sqrt(rhs) {
                if y != 0 {
                    return PrimePoint::affine(x, y);
                }
            }
        }
    }
}

/// **A curve of the family at about `bits` bits, with its group order
/// certified.**  The prime is random (constrained so the family's
/// endomorphism is rational), the curve is a random twist of the
/// family's model, the order is the unique CM candidate that kills
/// three random points, and `r` is its largest prime factor with the
/// cofactor at most `max_cofactor`.  Deterministic in `seed`.
///
/// `J0` is searched to prime order (`max_cofactor = 1` is
/// satisfiable); `J1728`, `D7` and `D8` always have a rational 2-torsion
/// point (the kernel of their degree-2 map), so their cofactor is at
/// least 2 and `max_cofactor` must allow it.  `Generic` counts points
/// in `O(p)` and is meant for `bits ≤ 26`.
pub fn generate_cm_instance(
    family: CmFamily,
    bits: u32,
    seed: u64,
    max_cofactor: u64,
) -> Result<PrimeInstance, String> {
    if !(10..=40).contains(&bits) {
        return Err(format!("bits = {bits} outside 10..=40"));
    }
    let mut rng = StdRng::seed_from_u64(seed ^ 0x434D_5F49_4E53_5400 ^ (bits as u64));
    for _attempt in 0..100_000u32 {
        let p = random_prime(bits, &mut rng);
        let admissible = match family {
            CmFamily::J0 => p % 3 == 1,
            CmFamily::J1728 => p % 4 == 1,
            CmFamily::D7 => pow_mod(p - 7, (p - 1) / 2, p) == 1,
            CmFamily::D8 => pow_mod(p - 8, (p - 1) / 2, p) == 1,
            CmFamily::D11 => pow_mod(p - 11, (p - 1) / 2, p) == 1,
            CmFamily::Generic => bits <= 26,
        };
        if !admissible {
            if family == CmFamily::Generic {
                return Err("generic instances count points in O(p): bits ≤ 26".into());
            }
            continue;
        }
        let (a, b) = match family {
            CmFamily::J0 => (0, rng.gen_range(1..p)),
            CmFamily::J1728 => (rng.gen_range(1..p), 0),
            CmFamily::D7 | CmFamily::D8 | CmFamily::D11 => {
                let j = match family {
                    CmFamily::D7 => -3375,
                    CmFamily::D8 => 8000,
                    _ => -32768,
                };
                let Some((a0, b0)) = curve_with_j(j, p) else {
                    continue;
                };
                let d = rng.gen_range(1..p);
                (
                    mulm(a0, mulm(d, d, p), p),
                    mulm(b0, mulm(mulm(d, d, p), d, p), p),
                )
            }
            CmFamily::Generic => (rng.gen_range(1..p), rng.gen_range(1..p)),
        };
        let curve = PrimeCurve { p, a, b };
        let disc = addm(
            mulm(4, mulm(mulm(a, a, p), a, p), p),
            mulm(27, mulm(b, b, p), p),
            p,
        );
        if disc == 0 {
            continue;
        }
        let order = match family {
            CmFamily::Generic => curve.point_count(),
            _ => {
                let Some(cands) = cm_candidate_orders(family, p) else {
                    continue;
                };
                let pts: Vec<PrimePoint> = (0..3).map(|_| random_point(&curve, &mut rng)).collect();
                let Some(n) = certified_order(&curve, &cands, &pts) else {
                    continue;
                };
                n
            }
        };
        let factors = factor_u64(order);
        let Some(&(r, _)) = factors.last() else {
            continue;
        };
        let h = order / r;
        if h > max_cofactor || r < 64 {
            continue;
        }
        // A generator of the prime-order subgroup.
        let mut ops = GroupOps::default();
        let g = loop {
            let pt = random_point(&curve, &mut rng);
            let g = curve.mul(&mut ops, pt, h);
            if !g.infinity {
                break g;
            }
        };
        if !curve.mul(&mut ops, g, r).infinity {
            return Err(format!(
                "certified order {order} is wrong on {}: [r]G ≠ O",
                curve.p
            ));
        }
        return Ok(PrimeInstance {
            name: format!("{}-{bits}bit-p{p}", family.name()),
            curve,
            group_order: order,
            r,
            cofactor: h,
            generator: (g.x, g.y),
        });
    }
    Err(format!(
        "no {} instance found at {bits} bits with cofactor ≤ {max_cofactor}",
        family.name()
    ))
}

// ── Polynomials over F_p, enough to find roots of a cubic ──────────

/// Coefficients low to high, no trailing zeros.
fn poly_trim(mut v: Vec<u64>) -> Vec<u64> {
    while v.last() == Some(&0) {
        v.pop();
    }
    v
}

fn poly_mul(a: &[u64], b: &[u64], p: u64, muls: &mut u64) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            out[i + j] = addm(out[i + j], mulm(x, y, p), p);
            *muls += 1;
        }
    }
    poly_trim(out)
}

fn poly_rem(a: &[u64], m: &[u64], p: u64, muls: &mut u64) -> Vec<u64> {
    let mut a = a.to_vec();
    let dm = m.len() - 1;
    let inv_lead = inv_mod(m[dm], p);
    while a.len() > dm {
        let da = a.len() - 1;
        let c = mulm(a[da], inv_lead, p);
        *muls += 1;
        if c != 0 {
            for k in 0..=dm {
                a[da - dm + k] = subm(a[da - dm + k], mulm(c, m[k], p), p);
            }
            *muls += dm as u64 + 1;
        }
        a.pop();
        a = poly_trim(a);
    }
    a
}

fn poly_gcd(a: &[u64], b: &[u64], p: u64, muls: &mut u64) -> Vec<u64> {
    let (mut a, mut b) = (poly_trim(a.to_vec()), poly_trim(b.to_vec()));
    while !b.is_empty() {
        let r = poly_rem(&a, &b, p, muls);
        a = b;
        b = r;
    }
    if let Some(&lead) = a.last() {
        let inv = inv_mod(lead, p);
        for c in &mut a {
            *c = mulm(*c, inv, p);
        }
        *muls += a.len() as u64;
    }
    a
}

/// `base^e mod m`.
fn poly_powmod(base: &[u64], mut e: u64, m: &[u64], p: u64, muls: &mut u64) -> Vec<u64> {
    let mut acc = vec![1u64];
    let mut b = poly_rem(base, m, p, muls);
    while e > 0 {
        if e & 1 == 1 {
            acc = poly_rem(&poly_mul(&acc, &b, p, muls), m, p, muls);
        }
        b = poly_rem(&poly_mul(&b, &b, p, muls), m, p, muls);
        e >>= 1;
    }
    acc
}

/// The roots in `F_p` of a polynomial (coefficients low to high), by
/// `gcd(x^p − x, f)` and Cantor–Zassenhaus splitting.  Deterministic
/// in `seed`.
pub fn poly_roots_fp(f: &[u64], p: u64, seed: u64) -> Vec<u64> {
    let mut muls = 0u64;
    poly_roots_fp_counted(f, p, seed, &mut muls)
}

/// [`poly_roots_fp`] with its `F_p` multiplications added to `muls`.
pub fn poly_roots_fp_counted(f: &[u64], p: u64, seed: u64, muls: &mut u64) -> Vec<u64> {
    let f = poly_trim(f.to_vec());
    if f.len() <= 1 {
        return Vec::new();
    }
    if f.len() == 2 {
        return vec![negm(mulm(f[0], inv_mod(f[1], p), p), p)];
    }
    // x^p − x mod f
    let xp = poly_powmod(&[0, 1], p, &f, p, muls);
    let mut xp_minus_x = xp;
    if xp_minus_x.len() < 2 {
        xp_minus_x.resize(2, 0);
    }
    xp_minus_x[1] = subm(xp_minus_x[1], 1, p);
    let g = poly_gcd(&f, &poly_trim(xp_minus_x), p, muls);
    let mut rng = StdRng::seed_from_u64(seed ^ 0x524F_4F54);
    let mut out = Vec::new();
    let mut stack = vec![g];
    while let Some(h) = stack.pop() {
        match h.len() {
            0 | 1 => {}
            2 => out.push(negm(mulm(h[0], inv_mod(h[1], p), p), p)),
            _ => {
                // Split h with gcd((x + δ)^((p−1)/2) − 1, h).
                let delta = rng.gen_range(0..p);
                let mut s = poly_powmod(&[delta, 1], (p - 1) / 2, &h, p, muls);
                if s.is_empty() {
                    s.push(0);
                }
                s[0] = subm(s[0], 1, p);
                let d = poly_gcd(&h, &poly_trim(s), p, muls);
                if d.len() <= 1 || d.len() == h.len() {
                    stack.push(h);
                    continue;
                }
                let mut q = h.clone();
                // h / d by long division (exact).
                let mut quotient = vec![0u64; h.len() - d.len() + 1];
                let inv_lead = inv_mod(*d.last().unwrap(), p);
                while q.len() >= d.len() && !q.is_empty() {
                    let dq = q.len() - 1;
                    let c = mulm(q[dq], inv_lead, p);
                    quotient[dq - (d.len() - 1)] = c;
                    for k in 0..d.len() {
                        q[dq - (d.len() - 1) + k] =
                            subm(q[dq - (d.len() - 1) + k], mulm(c, d[k], p), p);
                    }
                    *muls += d.len() as u64 + 1;
                    q.pop();
                    q = poly_trim(q);
                }
                stack.push(d);
                stack.push(poly_trim(quotient));
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

// ── Degree-2 CM endomorphisms by Vélu (type C) ─────────────────────

/// A degree-2 endomorphism: the 2-isogeny with kernel `(x0, 0)` by
/// Vélu's formulas, composed with the isomorphism `(X, Y) ↦ (u²X, u³Y)`
/// from its codomain back to the curve.  Exists over `F_p` when the
/// codomain is `F_p`-isomorphic to the curve — on a CM curve whose
/// order has an element of norm 2 (`D = −7, −8`, and `1 + i` on
/// `D = −4`).
#[derive(Clone, Debug, Serialize)]
pub struct VeluDegree2 {
    pub x0: u64,
    /// Vélu's `t = 3x0² + a`.
    pub t: u64,
    pub u2: u64,
    pub u3: u64,
    pub eigenvalue: u64,
    /// The trace of the endomorphism: `φ² − trace·φ + 2 = 0`.
    pub trace: i64,
}

impl Endomorphism<PrimeCurve> for VeluDegree2 {
    fn name(&self) -> String {
        format!("velu2[x0={}, trace={}]", self.x0, self.trace)
    }
    fn degree(&self) -> u64 {
        2
    }
    fn eigenvalue(&self) -> u64 {
        self.eigenvalue
    }
    fn apply(&self, g: &PrimeCurve, pt: PrimePoint) -> PrimePoint {
        if pt.infinity || pt.x == self.x0 {
            return PrimePoint::INFINITY;
        }
        let p = g.p;
        let d = subm(pt.x, self.x0, p);
        let inv = inv_mod(d, p);
        let inv2 = mulm(inv, inv, p);
        let xx = addm(pt.x, mulm(self.t, inv, p), p);
        let yy = mulm(pt.y, subm(1, mulm(self.t, inv2, p), p), p);
        PrimePoint::affine(mulm(self.u2, xx, p), mulm(self.u3, yy, p))
    }
}

/// The rational isomorphisms `(X, Y) ↦ (u²X, u³Y)` from
/// `y² = x³ + a₂x + b₂` to `y² = x³ + ax + b`, as `(u², u³)` pairs:
/// none when the curves are not `F_p`-isomorphic (a twist, or a
/// different `j`), otherwise both signs of `u`.
pub fn isomorphism_scalings(p: u64, a: u64, b: u64, a2: u64, b2: u64) -> Vec<(u64, u64)> {
    let u2_candidates: Vec<u64> = if a != 0 && b != 0 {
        if a2 == 0 || b2 == 0 {
            return Vec::new();
        }
        vec![mulm(mulm(b, a2, p), inv_mod(mulm(a, b2, p), p), p)]
    } else if b == 0 {
        if a2 == 0 || b2 != 0 {
            return Vec::new();
        }
        // u⁴ = a / a₂
        poly_roots_fp(&[negm(mulm(a, inv_mod(a2, p), p), p), 0, 1], p, 3)
    } else {
        if b2 == 0 || a2 != 0 {
            return Vec::new();
        }
        // u⁶ = b / b₂
        poly_roots_fp(&[negm(mulm(b, inv_mod(b2, p), p), p), 0, 0, 1], p, 3)
    };
    let mut out = Vec::new();
    for u2 in u2_candidates {
        if mulm(mulm(u2, u2, p), a2, p) != a || mulm(mulm(mulm(u2, u2, p), u2, p), b2, p) != b {
            continue;
        }
        let Some(u) = sqrt_mod(u2, p) else {
            continue; // the isomorphism is over F_{p²}: a twist, not rational
        };
        for uu in [u, p - u] {
            out.push((u2, mulm(u2, uu, p)));
        }
    }
    out
}

/// The eigenvalue of a degree-`d` endomorphism on `⟨G⟩` from its image
/// of `G`, among the roots of `λ² − tλ + d` for `t² < 4d`: `(trace, λ)`.
fn eigenvalue_of_degree(
    inst: &PrimeInstance,
    image_of_g: PrimePoint,
    d: u64,
) -> Option<(i64, u64)> {
    let r = inst.r;
    let g = inst.generator_point();
    let mut ops = GroupOps::default();
    let half = inv_mod(2, r);
    let tmax = ((4 * d) as f64).sqrt() as i64;
    for t in -tmax..=tmax {
        if (t * t) as u64 >= 4 * d {
            continue;
        }
        let tm = t.rem_euclid(r as i64) as u64;
        let disc = subm(mulm(tm, tm, r), (4 * d) % r, r);
        let Some(sq) = sqrt_mod(disc, r) else {
            continue;
        };
        for l in [
            mulm(addm(tm, sq, r), half, r),
            mulm(subm(tm, sq, r), half, r),
        ] {
            if inst.curve.mul(&mut ops, g, l) == image_of_g {
                return Some((t, l));
            }
        }
    }
    None
}

/// Every rational degree-2 endomorphism of the instance's curve found
/// through its rational 2-torsion points, each verified on `G` with an
/// eigenvalue among the roots of `λ² − tλ + 2` for `|t| ≤ 2`.  Empty
/// on a curve with no rational 2-torsion, or whose 2-isogenous curves
/// are not isomorphic to it over `F_p` (every generic curve).
pub fn velu_degree2_endomorphisms(inst: &PrimeInstance) -> Vec<VeluDegree2> {
    let curve = &inst.curve;
    let (p, a, b) = (curve.p, curve.a, curve.b);
    let mut out = Vec::new();
    let roots = poly_roots_fp(&[b, a, 0, 1], p, 7);
    let g = inst.generator_point();
    for x0 in roots {
        let t = addm(mulm(3, mulm(x0, x0, p), p), a, p);
        let w = mulm(x0, t, p);
        let a2 = subm(a, mulm(5, t, p), p);
        let b2 = subm(b, mulm(7, w, p), p);
        for (u2, u3) in isomorphism_scalings(p, a, b, a2, b2) {
            let mut endo = VeluDegree2 {
                x0,
                t,
                u2,
                u3,
                eigenvalue: 0,
                trace: 0,
            };
            let img = endo.apply(curve, g);
            if !curve.is_on_curve(img) {
                continue;
            }
            if let Some((tr, l)) = eigenvalue_of_degree(inst, img, 2) {
                endo.eigenvalue = l;
                endo.trace = tr;
                out.push(endo);
            }
        }
    }
    out
}

// ── Degree-3 CM endomorphisms by Vélu (type C) ─────────────────────

/// A degree-3 endomorphism: the 3-isogeny with kernel `{O, ±Q}`,
/// `Q = (x0, y0)` with `x0 ∈ F_p` a root of the 3-division polynomial
/// (`y0` may lie in `F_{p²}`; only `y0²` enters the formulas), by
/// Vélu's formulas for an odd-order kernel, composed with the
/// isomorphism back to the curve.  Exists over `F_p` when the CM order
/// has an element of norm 3: `D = −11` (`φ² − φ + 3 = 0`) and `√−3` on
/// `j = 0` (`φ² + 3 = 0`, kernel `x0 = 0`).
#[derive(Clone, Debug, Serialize)]
pub struct VeluDegree3 {
    pub x0: u64,
    /// `y0² = x0³ + a x0 + b`.
    pub y0_sq: u64,
    /// Vélu's `g^x = 3x0² + a`, `t = 2g^x`, `u = 4y0²`.
    pub gx: u64,
    pub t: u64,
    pub u: u64,
    pub u2: u64,
    pub u3: u64,
    pub eigenvalue: u64,
    /// `φ² − trace·φ + 3 = 0`.
    pub trace: i64,
}

impl Endomorphism<PrimeCurve> for VeluDegree3 {
    fn name(&self) -> String {
        format!("velu3[x0={}, trace={}]", self.x0, self.trace)
    }
    fn degree(&self) -> u64 {
        3
    }
    fn eigenvalue(&self) -> u64 {
        self.eigenvalue
    }
    fn apply(&self, g: &PrimeCurve, pt: PrimePoint) -> PrimePoint {
        if pt.infinity || pt.x == self.x0 {
            return PrimePoint::INFINITY;
        }
        let p = g.p;
        let inv = inv_mod(subm(pt.x, self.x0, p), p);
        let inv2 = mulm(inv, inv, p);
        let inv3 = mulm(inv2, inv, p);
        let xx = addm(
            addm(pt.x, mulm(self.t, inv, p), p),
            mulm(self.u, inv2, p),
            p,
        );
        let factor = subm(
            subm(1, mulm(2, mulm(self.u, inv3, p), p), p),
            mulm(2, mulm(self.gx, inv2, p), p),
            p,
        );
        let yy = mulm(pt.y, factor, p);
        PrimePoint::affine(mulm(self.u2, xx, p), mulm(self.u3, yy, p))
    }
}

/// Every rational degree-3 endomorphism of the instance's curve found
/// through the `F_p`-rational roots of its 3-division polynomial
/// `3x⁴ + 6ax² + 12bx − a²`, each verified on `G`.
pub fn velu_degree3_endomorphisms(inst: &PrimeInstance) -> Vec<VeluDegree3> {
    let curve = &inst.curve;
    let (p, a, b) = (curve.p, curve.a, curve.b);
    let mut out = Vec::new();
    let psi3 = [negm(mulm(a, a, p), p), mulm(12, b, p), mulm(6, a, p), 0, 3];
    let g = inst.generator_point();
    for x0 in poly_roots_fp(&psi3, p, 11) {
        let y0_sq = curve.rhs(x0);
        let gx = addm(mulm(3, mulm(x0, x0, p), p), a, p);
        let t = mulm(2, gx, p);
        let u = mulm(4, y0_sq, p);
        let w = addm(u, mulm(x0, t, p), p);
        let a2 = subm(a, mulm(5, t, p), p);
        let b2 = subm(b, mulm(7, w, p), p);
        for (u2, u3) in isomorphism_scalings(p, a, b, a2, b2) {
            let mut endo = VeluDegree3 {
                x0,
                y0_sq,
                gx,
                t,
                u,
                u2,
                u3,
                eigenvalue: 0,
                trace: 0,
            };
            let img = endo.apply(curve, g);
            if !curve.is_on_curve(img) {
                continue;
            }
            if let Some((tr, l)) = eigenvalue_of_degree(inst, img, 3) {
                endo.eigenvalue = l;
                endo.trace = tr;
                out.push(endo);
            }
        }
    }
    out
}

// ── Rho folded by the same group (E6) ──────────────────────────────

/// Classes of a Pollard walk under the finite group generated by a
/// list of endomorphisms: the representative of a point is the
/// least-keyed point of its orbit, with the scalar `μ` (a product of
/// eigenvalues) such that `rep = [μ]P`.  This is the matched rho
/// reference for a base folded by the same generators — the `√A`
/// that AGENTS.md §1 says a generic algorithm already takes.
pub struct EndomorphismClasses<'a, G: CountedGroup> {
    pub gens: Vec<&'a dyn Endomorphism<G>>,
    pub r: u64,
    /// The orbit length of the generator: the `A` of the floor.
    pub group_order: u32,
    max_orbit: usize,
}

impl<'a, G: CountedGroup> EndomorphismClasses<'a, G> {
    /// Build the classes and measure the group's order on `generator`.
    pub fn new(
        g: &G,
        generator: G::Elt,
        r: u64,
        gens: Vec<&'a dyn Endomorphism<G>>,
    ) -> Result<Self, String> {
        let mut classes = Self {
            gens,
            r,
            group_order: 1,
            max_orbit: 256,
        };
        let orbit = classes.orbit(g, generator)?;
        classes.group_order = orbit.len() as u32;
        Ok(classes)
    }

    /// Every `(Q, μ)` with `Q = [μ]P` in the orbit of `P`.
    fn orbit(&self, g: &G, p: G::Elt) -> Result<Vec<(G::Elt, u64)>, String> {
        let mut out: Vec<(G::Elt, u64)> = vec![(p, 1)];
        let mut i = 0usize;
        while i < out.len() {
            let (pt, c) = out[i];
            for e in &self.gens {
                let img = e.apply(g, pt);
                let coef = mulm(c, e.eigenvalue(), self.r);
                if let Some(&(_, c2)) = out.iter().find(|(q, _)| *q == img) {
                    if c2 != coef {
                        return Err(format!(
                            "{} reaches a point with two eigenvalue products",
                            e.name()
                        ));
                    }
                } else {
                    out.push((img, coef));
                    if out.len() > self.max_orbit {
                        return Err(
                            "the orbit does not close: a generator has infinite order".into()
                        );
                    }
                }
            }
            i += 1;
        }
        Ok(out)
    }
}

impl<G: CountedGroup> crate::cryptanalysis::ic_boundary::RhoClasses<G>
    for EndomorphismClasses<'_, G>
{
    fn automorphisms(&self) -> u32 {
        self.group_order
    }
    fn canon(&self, g: &G, p: G::Elt) -> (G::Elt, u64) {
        let orbit = match self.orbit(g, p) {
            Ok(o) => o,
            Err(_) => return (p, 1),
        };
        orbit
            .into_iter()
            .min_by_key(|(q, _)| g.key(q))
            .expect("the orbit holds p")
    }
    fn method(&self) -> &'static str {
        "endomorphism-folded r-adding walk: least-keyed orbit representative, look-ahead, short-cycle escape by doubling, distinguished points, stride starts"
    }
}

/// A counted rho whose classes are the orbits of `gens`: the matched
/// reference for a base folded by the same generators.
pub fn rho_reference_folded<G: CountedGroup>(
    g: &G,
    generator: G::Elt,
    target: G::Elt,
    r: u64,
    seed: u64,
    max_steps: u64,
    gens: &[&dyn Endomorphism<G>],
) -> Result<crate::cryptanalysis::ic_boundary::RhoResult, String> {
    let classes = EndomorphismClasses::new(g, generator, r, gens.to_vec())?;
    Ok(crate::cryptanalysis::ic_boundary::rho_walk_with(
        g,
        &classes,
        generator,
        target,
        r,
        seed,
        max_steps,
        crate::cryptanalysis::ic_boundary::RhoWalk::negation(),
    ))
}

// ── Tests ──────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::{
        collect_and_solve, rho_reference, RestartPool, TargetSource,
    };
    use crate::cryptanalysis::ic_framework::plugins::SubtractOracle;
    use crate::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx, Params};

    fn ctx_for(inst: &PrimeInstance, planted: u64) -> InstanceCtx<'_, PrimeCurve> {
        let mut ops = GroupOps::default();
        let g = inst.generator_point();
        InstanceCtx {
            group: &inst.curve,
            generator: g,
            target: inst.curve.mul(&mut ops, g, planted),
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        }
    }

    /// Run the shared relation loop with the `subtract` oracle over a
    /// base and return the recovered logarithm.
    fn solve_over(inst: &PrimeInstance, fb: &FactorBase<PrimePoint>, planted: u64) -> Option<u64> {
        let ctx = ctx_for(inst, planted);
        let mut oracle = SubtractOracle;
        oracle
            .prepare(&ctx, fb, &Params::default(), &mut GroupOps::default())
            .unwrap();
        let out = collect_and_solve(
            &inst.curve,
            ctx.generator,
            ctx.target,
            inst.r,
            inst.cofactor,
            fb,
            11,
            4_000_000,
            TargetSource::Walk,
            RestartPool::Lazy,
            |ops, ctr, pt| oracle.decompose(&ctx, fb, ops, ctr, pt),
        );
        assert!(out.verified || out.recovered.is_none());
        out.recovered
    }

    #[test]
    fn eigenvalue_order_is_the_multiplicative_order() {
        // 2 has order 10 mod 11; 3 has order 5; 10 has order 2.
        assert_eq!(eigenvalue_order(2, 11), 10);
        assert_eq!(eigenvalue_order(3, 11), 5);
        assert_eq!(eigenvalue_order(10, 11), 2);
        assert_eq!(eigenvalue_order(1, 11), 1);
    }

    #[test]
    fn cornacchia_finds_the_cm_representation() {
        // 4·31 = 124 = 11² + 3·1² (D = −3); 4·13 = 52 = 6² + 4·2² (D = −4).
        assert_eq!(cornacchia_4p(31, 3), Some((11, 1)));
        let (t, v) = cornacchia_4p(13, 4).unwrap();
        assert_eq!(t * t + 4 * v * v, 52);
        // p ≡ 2 (mod 3) has no D = −3 representation.
        assert_eq!(cornacchia_4p(11, 3), None);
    }

    #[test]
    fn poly_roots_finds_every_root_of_a_split_cubic() {
        // (x − 2)(x − 5)(x − 9) over F_101.
        let p = 101;
        let mut muls = 0u64;
        let f = poly_mul(
            &poly_mul(&[p - 2, 1], &[p - 5, 1], p, &mut muls),
            &[p - 9, 1],
            p,
            &mut muls,
        );
        assert_eq!(poly_roots_fp(&f, p, 1), vec![2, 5, 9]);
        // x² + 1 over F_7 (7 ≡ 3 mod 4) has no root.
        assert!(poly_roots_fp(&[1, 0, 1], 7, 1).is_empty());
        // x² + 1 over F_13 has roots 5 and 8.
        assert_eq!(poly_roots_fp(&[1, 0, 1], 13, 1), vec![5, 8]);
    }

    #[test]
    fn j0_instance_has_a_verified_order_three_automorphism() {
        let inst = generate_cm_instance(CmFamily::J0, 16, 1, 1).unwrap();
        assert_eq!(inst.cofactor, 1);
        assert_eq!(inst.curve.a, 0);
        assert_eq!(inst.curve.p % 3, 1);
        let psi = j0_automorphism(&inst).unwrap();
        let check =
            verify_endomorphism(&inst.curve, inst.generator_point(), inst.r, &psi, 50, 1).unwrap();
        assert_eq!(check.eigenvalue_order, 3);
        assert_eq!(
            addm(
                addm(
                    mulm(psi.eigenvalue, psi.eigenvalue, inst.r),
                    psi.eigenvalue,
                    inst.r
                ),
                1,
                inst.r
            ),
            0,
            "λ² + λ + 1 ≡ 0"
        );
    }

    #[test]
    fn j1728_instance_has_a_verified_order_four_automorphism() {
        let inst = generate_cm_instance(CmFamily::J1728, 16, 2, 4).unwrap();
        assert!(
            inst.cofactor == 2 || inst.cofactor == 4,
            "(0, 0) is rational 2-torsion"
        );
        assert_eq!(inst.curve.b, 0);
        assert_eq!(inst.curve.p % 4, 1);
        let iota = j1728_automorphism(&inst).unwrap();
        let check =
            verify_endomorphism(&inst.curve, inst.generator_point(), inst.r, &iota, 50, 2).unwrap();
        assert_eq!(check.eigenvalue_order, 4);
        assert_eq!(
            mulm(iota.eigenvalue, iota.eigenvalue, inst.r),
            inst.r - 1,
            "λ² ≡ −1"
        );
    }

    #[test]
    fn a_wrong_eigenvalue_is_caught_by_verification() {
        let inst = generate_cm_instance(CmFamily::J0, 14, 3, 1).unwrap();
        let mut psi = j0_automorphism(&inst).unwrap();
        psi.eigenvalue = addm(psi.eigenvalue, 1, inst.r);
        let err = verify_endomorphism(&inst.curve, inst.generator_point(), inst.r, &psi, 5, 3)
            .unwrap_err();
        assert!(err.contains("≠"), "{err}");
    }

    /// **The headline invariant of the fold**: `6` signed points per
    /// column on `j = 0`, `4` on `j = 1728`, `2` on a generic curve,
    /// and the control on the same points has `2` everywhere.
    #[test]
    fn the_fold_divides_the_columns_by_the_automorphism_group() {
        for (family, seed, per_column, h) in [
            (CmFamily::J0, 5, 6.0, 1),
            (CmFamily::J1728, 6, 4.0, 4),
            (CmFamily::Generic, 7, 2.0, 1),
        ] {
            let inst = generate_cm_instance(family, 14, seed, h).unwrap();
            let (folded, rep) = glv_orbit_base(&inst, 12, AutomorphismGroup::Auto, true).unwrap();
            let (control, crep) =
                glv_orbit_base(&inst, 12, AutomorphismGroup::Auto, false).unwrap();
            assert_eq!(
                folded.points.len(),
                control.points.len(),
                "{family:?}: same point set"
            );
            assert!(
                (rep.points_per_orbit - per_column).abs() < 1e-9,
                "{family:?}: {} points per column, expected {per_column}",
                rep.points_per_orbit
            );
            assert!((crep.points_per_orbit - 2.0).abs() < 1e-9);
            assert_eq!(folded.columns * (per_column as usize) / 2, control.columns);
            // Every coefficient is a unit, and every column has a
            // representative with coefficient 1.
            let mut has_rep = vec![false; folded.columns];
            for (i, &c) in folded.coef_of.iter().enumerate() {
                assert!(c != 0 && c < inst.r);
                if c == 1 {
                    has_rep[folded.col_of[i]] = true;
                }
            }
            assert!(
                has_rep.iter().all(|&b| b),
                "{family:?}: a column without representative"
            );
        }
    }

    /// A base that is not invariant is refused under `Strict`, and the
    /// error names the map.
    #[test]
    fn strict_closure_refuses_a_non_invariant_base() {
        let inst = generate_cm_instance(CmFamily::J0, 14, 8, 1).unwrap();
        let psi = j0_automorphism(&inst).unwrap();
        let neg = Negation { r: inst.r };
        let seed = smallest_abscissa_points(&inst, 10);
        let gens: Vec<&dyn Endomorphism<PrimeCurve>> = vec![&neg, &psi];
        let curve = &inst.curve;
        let err = fold_by_endomorphisms(
            curve,
            inst.r,
            1,
            seed,
            &gens,
            Closure::Strict,
            64,
            |p| curve.key(p),
            |p| p.x,
            String::new(),
        )
        .err()
        .expect("a non-invariant base is refused");
        assert!(err.contains("not invariant"), "{err}");
    }

    /// **The fold's coefficients are right**: over a folded base the
    /// shared relation loop recovers the planted logarithm on every
    /// family, and the control over the same points recovers the same
    /// one.  A wrong `λ^k` would make the matrix solve a different
    /// problem and the verification would say so.
    #[test]
    fn a_folded_base_recovers_the_planted_logarithm() {
        for (family, seed, h) in [
            (CmFamily::J0, 21, 1),
            (CmFamily::J1728, 22, 4),
            (CmFamily::Generic, 23, 1),
        ] {
            let inst = generate_cm_instance(family, 15, seed, h).unwrap();
            let planted = 12_345 % inst.r;
            let (folded, _) = glv_orbit_base(&inst, 24, AutomorphismGroup::Auto, true).unwrap();
            let (control, _) = glv_orbit_base(&inst, 24, AutomorphismGroup::Auto, false).unwrap();
            let a = solve_over(&inst, &folded, planted);
            let b = solve_over(&inst, &control, planted);
            assert_eq!(a, Some(planted), "{family:?}: folded base");
            assert_eq!(b, Some(planted), "{family:?}: control base");
        }
    }

    /// The fold works on a cofactor curve too: points of `E` outside
    /// `⟨G⟩` fold by the same eigenvalue because relations are written
    /// over `[h]P`.
    #[test]
    fn the_fold_is_sound_on_a_cofactor_curve() {
        let inst = generate_cm_instance(CmFamily::D7, 15, 31, 8).unwrap();
        assert!(inst.cofactor >= 2, "D = −7 has rational 2-torsion");
        let planted = 777 % inst.r;
        let (fb, rep) = glv_orbit_base(&inst, 24, AutomorphismGroup::Auto, true).unwrap();
        assert_eq!(rep.generators, vec!["negation".to_string()]);
        assert_eq!(solve_over(&inst, &fb, planted), Some(planted));
    }

    /// **Type C, measured.**  The `D = −7` and `D = −8` curves carry a
    /// rational degree-2 endomorphism: Vélu finds it, verification
    /// accepts it, its eigenvalue satisfies `λ² − tλ + 2 = 0`, its
    /// order modulo `r` is large, and it keeps only a chance fraction
    /// of a base inside the base — so the fold refuses it.
    #[test]
    fn degree_two_cm_endomorphisms_exist_verify_and_do_not_fold() {
        for (family, seed, traces) in [(CmFamily::D7, 41, vec![1, -1]), (CmFamily::D8, 42, vec![0])]
        {
            let inst = generate_cm_instance(family, 16, seed, 8).unwrap();
            let endos = velu_degree2_endomorphisms(&inst);
            assert!(
                !endos.is_empty(),
                "{family:?}: no degree-2 endomorphism found"
            );
            let phi = &endos[0];
            assert!(
                traces.contains(&phi.trace),
                "{family:?}: trace {}",
                phi.trace
            );
            let check =
                verify_endomorphism(&inst.curve, inst.generator_point(), inst.r, phi, 30, seed)
                    .unwrap();
            let (l, r) = (phi.eigenvalue, inst.r);
            let tm = phi.trace.rem_euclid(r as i64) as u64;
            assert_eq!(addm(subm(mulm(l, l, r), mulm(tm, l, r), r), 2, r), 0);
            assert!(
                check.eigenvalue_order > 64,
                "{family:?}: ord_r(λ) = {} is small",
                check.eigenvalue_order
            );
            let (fb, _) = glv_orbit_base(&inst, 32, AutomorphismGroup::Negation, true).unwrap();
            let overlap = endomorphism_overlap(&inst.curve, &fb, phi, inst.r, inst.group_order);
            assert!(
                overlap.fraction < 0.25,
                "{family:?}: {} of the base maps into the base",
                overlap.fraction
            );
            let neg = Negation { r: inst.r };
            let gens: Vec<&dyn Endomorphism<PrimeCurve>> = vec![&neg, phi];
            let curve = &inst.curve;
            let err = fold_by_endomorphisms(
                curve,
                inst.r,
                inst.cofactor,
                fb.points.clone(),
                &gens,
                Closure::Close,
                256,
                |p| curve.key(p),
                |p| p.x,
                String::new(),
            )
            .err()
            .expect("a degree-2 map cannot fold");
            assert!(
                err.contains("finite order")
                    || err.contains("kernel")
                    || err.contains("eigenvalues"),
                "{err}"
            );
        }
    }

    /// `1 + i` on `j = 1728` is a degree-2 endomorphism of a curve that
    /// also has the order-4 automorphism: type C living beside type A.
    #[test]
    fn one_plus_i_is_a_degree_two_endomorphism_of_a_j1728_curve() {
        // Need rational 2-torsion: y² = x³ + ax always has (0, 0).
        let inst = generate_cm_instance(CmFamily::J1728, 15, 51, 4).unwrap();
        let endos = velu_degree2_endomorphisms(&inst);
        assert!(!endos.is_empty(), "no degree-2 map through (0, 0)");
        let phi = &endos[0];
        assert_eq!(phi.trace.abs(), 2, "1 ± i has trace ±2");
        verify_endomorphism(&inst.curve, inst.generator_point(), inst.r, phi, 20, 51).unwrap();
    }

    /// **E4b.**  `D = −11` carries a rational degree-3 endomorphism with
    /// `φ² ∓ φ + 3 = 0`, and every `j = 0` curve carries `√−3` (kernel
    /// `x0 = 0`, `φ² + 3 = 0`); Vélu finds both, verification accepts
    /// them, and their eigenvalue orders are large.
    #[test]
    fn degree_three_cm_endomorphisms_exist_and_verify() {
        let inst = generate_cm_instance(CmFamily::D11, 16, 71, 16).unwrap();
        let endos = velu_degree3_endomorphisms(&inst);
        assert!(!endos.is_empty(), "D = −11: no degree-3 endomorphism found");
        let phi = &endos[0];
        assert_eq!(phi.trace.abs(), 1, "trace of (1 ± √−11)/2");
        let check =
            verify_endomorphism(&inst.curve, inst.generator_point(), inst.r, phi, 30, 71).unwrap();
        let (l, r) = (phi.eigenvalue, inst.r);
        let tm = phi.trace.rem_euclid(r as i64) as u64;
        assert_eq!(addm(subm(mulm(l, l, r), mulm(tm, l, r), r), 3, r), 0);
        assert!(
            check.eigenvalue_order > 64,
            "ord = {}",
            check.eigenvalue_order
        );

        let inst = generate_cm_instance(CmFamily::J0, 16, 72, 4).unwrap();
        let endos = velu_degree3_endomorphisms(&inst);
        assert!(!endos.is_empty(), "j = 0: √−3 not found");
        let phi = endos.iter().find(|e| e.x0 == 0).expect("kernel at x = 0");
        assert_eq!(phi.trace, 0, "φ² = −3");
        verify_endomorphism(&inst.curve, inst.generator_point(), inst.r, phi, 30, 72).unwrap();
        assert_eq!(
            mulm(phi.eigenvalue, phi.eigenvalue, inst.r),
            inst.r - 3,
            "λ² ≡ −3"
        );
        assert!(eigenvalue_order(phi.eigenvalue, inst.r) > 64);
    }

    /// **E6.**  The rho folded by `⟨−1, ζ⟩` walks on classes of six
    /// points, verifies its answer, and takes fewer steps than the
    /// negation walk on the same instance and targets.
    #[test]
    fn the_folded_rho_walks_on_the_group_s_classes() {
        use crate::cryptanalysis::ic_boundary::rho_reference_negation;
        let inst = generate_cm_instance(CmFamily::J0, 18, 81, 1).unwrap();
        let gens = automorphism_generators(&inst, AutomorphismGroup::Auto).unwrap();
        let refs: Vec<&dyn Endomorphism<PrimeCurve>> = gens.iter().map(|b| b.as_ref()).collect();
        let g = inst.generator_point();
        let mut ops = GroupOps::default();
        let (mut folded_steps, mut neg_steps) = (0u64, 0u64);
        for k in 0..6u64 {
            let planted = 1 + (k * 7919) % (inst.r - 1);
            let q = inst.curve.mul(&mut ops, g, planted);
            let f =
                rho_reference_folded(&inst.curve, g, q, inst.r, 100 + k, 1 << 30, &refs).unwrap();
            assert_eq!(f.automorphisms, 6);
            assert!(f.verified, "the folded walk must verify its answer");
            assert_eq!(f.recovered, Some(planted));
            let n = rho_reference_negation(&inst.curve, g, q, inst.r, 100 + k, 1 << 30);
            assert!(n.verified);
            folded_steps += f.steps;
            neg_steps += n.steps;
        }
        assert!(
            folded_steps < neg_steps,
            "folded {folded_steps} steps against negation {neg_steps}"
        );
    }

    /// The counted rho reference runs on a generated instance, so the
    /// bench can put every family's row against its own rho.
    #[test]
    fn rho_reference_runs_on_a_generated_instance() {
        let inst = generate_cm_instance(CmFamily::J0, 14, 61, 1).unwrap();
        let mut ops = GroupOps::default();
        let g = inst.generator_point();
        let q = inst.curve.mul(&mut ops, g, 99);
        let res = rho_reference(&inst.curve, g, q, inst.r, 5, 1 << 22);
        assert!(res.s > 0.0);
    }
}
