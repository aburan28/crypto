//! # A pair table over orbit representatives (E12).
//!
//! The meet-in-the-middle oracle of the framework stores every pair sum
//! `P_i + P_j` of a base, keyed by the sum's abscissa so that the
//! negation is folded out (`PairTable::build_negation_folded`): about
//! `|F|²/4` additions and entries.  When the base is closed under a
//! finite group `H` of endomorphisms — the fold of
//! `research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`
//! — the group acts on pairs diagonally, `h(P + Q) = h(P) + h(Q)`, so one
//! entry per `H`-orbit of pairs is enough: the pairs `(rep_c, P_j)` with
//! `rep_c` the first point of column `c` and `col(P_j) ≥ c` meet every
//! orbit, which is `|F|²/(2w)` additions for `|H| = w`, a factor `w/2`
//! below the negation table.  This is the odd-characteristic,
//! any-group form of `FrobeniusPairTable` (Koblitz).
//!
//! The table is keyed by an **orbit invariant** of the sum: a number
//! equal on two points exactly when they lie in one `H`-orbit
//! ([`OrbitKey`]).  For `⟨−1, π⟩` on a subfield curve over `F_{p³}` it
//! is the minimal polynomial of `x` over `F_p` ([`SubfieldOrbitKey`],
//! about a dozen `F_p` multiplications); for `⟨−1, ζ⟩` on a `j = 0`
//! curve over `F_p` it is `x³`, for `⟨−1, ι⟩` on `j = 1728` it is `x²`
//! ([`PrimePowerKey`]).  A probe pays the invariant once; on a hit it
//! adds the stored pair, finds the group element taking that sum to the
//! target by walking the sum's orbit, reads the two image points off
//! the base's column map, and checks with one more addition that they
//! sum to the target.  Hits are rare, so the probe's cost is the
//! invariant: the table is `w/2` smaller and `w/2` cheaper to build, and
//! each probe dearer by the invariant.
//!
//! [`OrbitMitmOracle`] is the oracle over it, `m = 2` (one probe) or
//! `m = 3` (one probe per `R − P_i`), the same loop as `decompose_mitm`.

use std::time::Instant;

use crate::cryptanalysis::ext_curve::{ExtCurve, Fp3};
use crate::cryptanalysis::glv_invariant_base::{mulm, EndomorphismClasses};
use crate::cryptanalysis::ic_boundary::{
    fast_map, CountedGroup, FactorBase, FastMap, GroupOps, OracleCounters, PrimeCurve,
};
use crate::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx, Params};

/// An invariant of the orbits of a finite group acting on a curve.
pub trait OrbitKey<G: CountedGroup> {
    fn name(&self) -> &'static str;
    /// Equal on two points exactly when one is the image of the other
    /// under the group; never 0 for an affine point.
    fn key(&self, g: &G, p: &G::Elt) -> u64;
    /// `F_p` multiplications one [`OrbitKey::key`] costs that the group's
    /// own arithmetic does not count (the prime-field curve counts none;
    /// `ExtCurve` counts its field's multiplications itself).
    fn uncounted_muls(&self) -> u64;
    /// The same for one application of a group generator on a hit.
    fn uncounted_muls_per_map(&self) -> u64;
}

/// `⟨−1, π⟩` on a subfield curve `E / F_p` over `F_{p³}`: the orbit of
/// `P` has abscissae `{x, x^p, x^{p²}}`, the roots of the minimal
/// polynomial `T³ − σ₁T² + σ₂T − σ₃` of `x`, so `(σ₁, σ₂, σ₃)` is the
/// invariant.  For `x = a + b s + c s²` with `s³ = ν`: `σ₁ = Tr x = 3a`,
/// `σ₂ = (σ₁² − Tr x²)/2 = 3(a² − ν b c)`, `σ₃ = N x = a³ + ν b³ + ν² c³
/// − 3ν a b c` — twelve `F_p` multiplications, uncounted by the group
/// and charged by the oracle.
pub struct SubfieldOrbitKey;

impl OrbitKey<ExtCurve<Fp3>> for SubfieldOrbitKey {
    fn name(&self) -> &'static str {
        "minimal polynomial of x over F_p"
    }
    fn key(&self, g: &ExtCurve<Fp3>, pt: &<ExtCurve<Fp3> as CountedGroup>::Elt) -> u64 {
        if pt.infinity {
            return 0;
        }
        let p = g.f.p;
        let nu = g.f.nu;
        let [a, b, c] = pt.x;
        let a2 = mulm(a, a, p);
        let bc = mulm(b, c, p);
        let nubc = mulm(nu, bc, p);
        let s1 = (3 * a) % p;
        let s2 = (3 * ((a2 + p - nubc) % p)) % p;
        let a3 = mulm(a2, a, p);
        let b3 = mulm(mulm(b, b, p), b, p);
        let c3 = mulm(mulm(c, c, p), c, p);
        let nu_b3 = mulm(nu, b3, p);
        let nu2_c3 = mulm(mulm(nu, nu, p), c3, p);
        let abc3nu = (3 * mulm(a, nubc, p)) % p;
        let s3 = ((a3 + nu_b3) % p + (nu2_c3 + p - abc3nu) % p) % p;
        1 + (s1 * p + s2) * p + s3
    }
    fn uncounted_muls(&self) -> u64 {
        12
    }
    fn uncounted_muls_per_map(&self) -> u64 {
        0
    }
}

/// `⟨−1, ζ⟩` on `j = 0` (`x ↦ βx`, `β³ = 1`): key `x³`; `⟨−1, ι⟩` on
/// `j = 1728` (`x ↦ −x`): key `x²`; negation alone: key `x`.  The orbit's
/// abscissae are exactly the `e`-th roots of `x^e`.
pub struct PrimePowerKey {
    pub exponent: u32,
}

impl OrbitKey<PrimeCurve> for PrimePowerKey {
    fn name(&self) -> &'static str {
        match self.exponent {
            1 => "x",
            2 => "x²",
            3 => "x³",
            _ => "x^e",
        }
    }
    fn key(&self, g: &PrimeCurve, pt: &<PrimeCurve as CountedGroup>::Elt) -> u64 {
        if pt.infinity {
            return 0;
        }
        let mut v = pt.x % g.p;
        for _ in 1..self.exponent {
            v = mulm(v, pt.x, g.p);
        }
        v + 1
    }
    fn uncounted_muls(&self) -> u64 {
        u64::from(self.exponent.saturating_sub(1))
    }
    fn uncounted_muls_per_map(&self) -> u64 {
        1
    }
}

/// The table: one entry per `H`-orbit of pairs, keyed by the orbit
/// invariant of the sum.
pub struct OrbitPairTable {
    map: FastMap<(u32, u32)>,
    /// The base indices of each column, for reading `h(P_i)` off the
    /// column map on a hit.
    members: Vec<Vec<usize>>,
    pub entries: u64,
    pub representatives: u64,
    pub build_ops: GroupOps,
    pub build_keys: u64,
    pub build_wall_ns: u64,
}

impl OrbitPairTable {
    /// Build over `fb`, whose columns must be the orbits of the group
    /// `key` is an invariant of.
    pub fn build<G: CountedGroup, K: OrbitKey<G>>(g: &G, fb: &FactorBase<G::Elt>, key: &K) -> Self {
        let start = Instant::now();
        let n = fb.points.len();
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); fb.columns];
        for i in 0..n {
            members[fb.col_of[i]].push(i);
        }
        let mut ops = GroupOps::default();
        let mut keys = 0u64;
        let mut map: FastMap<(u32, u32)> = fast_map(n * n / (2 * fb.columns.max(1)).max(1) + 1);
        let mut representatives = 0u64;
        for (c, col) in members.iter().enumerate() {
            let Some(&rep) = col.first() else {
                continue;
            };
            representatives += 1;
            for j in 0..n {
                if fb.col_of[j] < c {
                    continue;
                }
                let s = g.add(&mut ops, fb.points[rep], fb.points[j]);
                if g.is_identity(&s) {
                    continue;
                }
                keys += 1;
                map.entry(key.key(g, &s)).or_insert((rep as u32, j as u32));
            }
        }
        Self {
            entries: map.len() as u64,
            representatives,
            map,
            members,
            build_ops: ops,
            build_keys: keys,
            build_wall_ns: start.elapsed().as_nanos() as u64,
        }
    }

    /// The base point `h(P_i)` for the group element `h` whose eigenvalue
    /// is `lambda`: the member of `P_i`'s column carrying `λ·coef(P_i)`.
    fn image<E>(&self, fb: &FactorBase<E>, i: usize, lambda: u64, r: u64) -> Option<usize> {
        let want = mulm(lambda, fb.coef_of[i], r);
        self.members[fb.col_of[i]]
            .iter()
            .copied()
            .find(|&j| fb.coef_of[j] == want)
    }

    /// Indices `(a, b)` with `P_a + P_b = target`, if some group image of
    /// a stored pair sums to it.
    #[allow(clippy::too_many_arguments)]
    pub fn probe<G: CountedGroup, K: OrbitKey<G>>(
        &self,
        g: &G,
        fb: &FactorBase<G::Elt>,
        key: &K,
        classes: &EndomorphismClasses<'_, G>,
        ops: &mut GroupOps,
        ctr: &mut OracleCounters,
        stats: &mut OrbitProbeStats,
        target: G::Elt,
    ) -> Option<(usize, usize)> {
        ctr.canonicalisations += 1;
        stats.keys += 1;
        let k = key.key(g, &target);
        ctr.lookups += 1;
        let &(i, j) = self.map.get(&k)?;
        let (i, j) = (i as usize, j as usize);
        stats.hits += 1;
        let s = g.add(ops, fb.points[i], fb.points[j]);
        // The group element taking the stored sum to the target.
        let orbit = classes.orbit(g, s).ok()?;
        stats.maps += (orbit.len() * classes.gens.len()) as u64;
        let Some(&(_, lambda)) = orbit.iter().find(|(q, _)| *q == target) else {
            stats.key_collisions += 1;
            return None;
        };
        let (a, b) = (
            self.image(fb, i, lambda, classes.r)?,
            self.image(fb, j, lambda, classes.r)?,
        );
        if g.add(ops, fb.points[a], fb.points[b]) != target {
            stats.image_mismatches += 1;
            return None;
        }
        Some((a, b))
    }
}

/// What the probes did beyond the framework's counters.
#[derive(Clone, Copy, Debug, Default, serde::Serialize)]
pub struct OrbitProbeStats {
    pub keys: u64,
    pub hits: u64,
    /// Generator applications walking a stored sum's orbit on a hit.
    pub maps: u64,
    /// Hits whose key matched but whose orbit did not hold the target;
    /// zero for an exact invariant.
    pub key_collisions: u64,
    /// Hits whose read-off images did not sum to the target; a wrong
    /// column map would show here, so it must stay zero.
    pub image_mismatches: u64,
}

/// Meet in the middle over an [`OrbitPairTable`].
pub struct OrbitMitmOracle<'a, G: CountedGroup, K: OrbitKey<G>> {
    m: u32,
    key: K,
    classes: EndomorphismClasses<'a, G>,
    table: Option<OrbitPairTable>,
    pub stats: OrbitProbeStats,
}

impl<'a, G: CountedGroup, K: OrbitKey<G>> OrbitMitmOracle<'a, G, K> {
    pub fn new(m: u32, key: K, classes: EndomorphismClasses<'a, G>) -> Self {
        Self {
            m,
            key,
            classes,
            table: None,
            stats: OrbitProbeStats::default(),
        }
    }

    pub fn table(&self) -> Option<&OrbitPairTable> {
        self.table.as_ref()
    }

    /// `F_p` multiplications the group did not count: the probes' and the
    /// build's invariants, and the generator maps on hits.
    pub fn uncounted_muls(&self) -> u64 {
        let build = self.table.as_ref().map(|t| t.build_keys).unwrap_or(0);
        (build + self.stats.keys) * self.key.uncounted_muls()
            + self.stats.maps * self.key.uncounted_muls_per_map()
    }

    pub fn build_uncounted_muls(&self) -> u64 {
        self.table.as_ref().map(|t| t.build_keys).unwrap_or(0) * self.key.uncounted_muls()
    }
}

impl<G: CountedGroup, K: OrbitKey<G>> DecompositionOracle<G> for OrbitMitmOracle<'_, G, K> {
    fn name(&self) -> &str {
        "orbit-mitm"
    }

    fn summands(&self) -> u32 {
        self.m
    }

    fn describe(&self, _params: &Params) -> String {
        format!(
            "a table of pair sums over the base's orbit representatives (one entry per orbit of pairs under a group of order {}), keyed by {}, probed once per target (m = 2) or once per R − P (m = 3)",
            self.classes.group_order,
            self.key.name()
        )
    }

    fn prepare(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        _params: &Params,
        ops: &mut GroupOps,
    ) -> Result<(), String> {
        let table = OrbitPairTable::build(ctx.group, fb, &self.key);
        ops.merge(table.build_ops);
        self.table = Some(table);
        Ok(())
    }

    fn decompose(
        &mut self,
        ctx: &InstanceCtx<G>,
        fb: &FactorBase<G::Elt>,
        ops: &mut GroupOps,
        counters: &mut OracleCounters,
        point: G::Elt,
    ) -> Option<Vec<usize>> {
        let table = self.table.as_ref()?;
        let g = ctx.group;
        match self.m {
            2 => table
                .probe(
                    g,
                    fb,
                    &self.key,
                    &self.classes,
                    ops,
                    counters,
                    &mut self.stats,
                    point,
                )
                .map(|(a, b)| vec![a, b]),
            3 => {
                for i in 0..fb.points.len() {
                    let s = g.add(ops, point, g.neg(fb.points[i]));
                    if g.is_identity(&s) {
                        continue;
                    }
                    if let Some((a, b)) = table.probe(
                        g,
                        fb,
                        &self.key,
                        &self.classes,
                        ops,
                        counters,
                        &mut self.stats,
                        s,
                    ) {
                        return Some(vec![i, a, b]);
                    }
                }
                None
            }
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::glv_invariant_base::{
        automorphism_generators, generate_cm_instance, glv_orbit_base, AutomorphismGroup, CmFamily,
        Endomorphism, Negation,
    };
    use crate::cryptanalysis::ic_framework::plugins::MitmOracle;
    use crate::cryptanalysis::subfield_fp3::{
        generate_subfield_instance, subfield_line_base, SubfieldFold,
    };

    /// The orbit table is `w/2` smaller than the negation table, and the
    /// two oracles decompose the same targets, every returned triple
    /// summing to its target.
    #[test]
    fn the_orbit_table_agrees_with_the_negation_table_on_a_subfield_line() {
        let inst = generate_subfield_instance(7, 2, false, 8).unwrap();
        let (fb, _) = subfield_line_base(&inst, SubfieldFold::Frobenius).unwrap();
        let mut ops = GroupOps::default();
        let target = inst.curve.mul(&mut ops, inst.generator, 3);
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: inst.generator,
            target,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: Some(3),
        };
        let neg = Negation { r: inst.r };
        let gens: Vec<&dyn Endomorphism<_>> = vec![&neg, &inst.frobenius];
        let classes = EndomorphismClasses::new(&inst.curve, inst.generator, inst.r, gens).unwrap();
        assert_eq!(classes.group_order, 6);
        let mut orbit = OrbitMitmOracle::new(3, SubfieldOrbitKey, classes);
        orbit
            .prepare(&ctx, &fb, &Params::default(), &mut ops)
            .unwrap();
        let mut plain = MitmOracle::new(3);
        let mut params = Params::default();
        params.set("negation_folded", "1");
        let mut plain_ops = GroupOps::default();
        plain.prepare(&ctx, &fb, &params, &mut plain_ops).unwrap();
        let t = orbit.table().unwrap();
        assert!(
            (t.build_ops.adds as f64) * 2.4 < plain_ops.adds as f64,
            "orbit build {} against negation build {}",
            t.build_ops.adds,
            plain_ops.adds
        );
        let (mut ca, mut cb) = (OracleCounters::default(), OracleCounters::default());
        let mut hits = 0;
        for k in 2..150u64 {
            let pt = inst.curve.mul(&mut ops, inst.generator, k);
            let a = orbit.decompose(&ctx, &fb, &mut ops, &mut ca, pt);
            let b = plain.decompose(&ctx, &fb, &mut ops, &mut cb, pt);
            assert_eq!(a.is_some(), b.is_some(), "k = {k}");
            if let Some(idx) = a {
                let sum = idx.iter().fold(inst.curve.identity(), |acc, &i| {
                    inst.curve.add(&mut ops, acc, fb.points[i])
                });
                assert_eq!(sum, pt);
                hits += 1;
            }
        }
        assert!(hits > 0);
        assert_eq!(orbit.stats.key_collisions, 0);
        assert_eq!(orbit.stats.image_mismatches, 0);
    }

    /// The closed-form key is the minimal polynomial of `x`: equal on the
    /// six points of an orbit, and it agrees with the field's own trace
    /// and norm.
    #[test]
    fn the_subfield_key_is_the_minimal_polynomial_and_constant_on_orbits() {
        use crate::cryptanalysis::ext_curve::{random_point, ExtField};
        use rand::SeedableRng;
        let inst = generate_subfield_instance(8, 1, false, 8).unwrap();
        let f = &inst.curve.f;
        let p = f.p;
        let mut rng = rand::rngs::StdRng::seed_from_u64(4);
        for _ in 0..50 {
            let pt = random_point(&inst.curve, &mut rng);
            let k = SubfieldOrbitKey.key(&inst.curve, &pt);
            let tr = (3 * pt.x[0]) % p;
            let n = f.norm(pt.x);
            assert_eq!((k - 1) / (p * p), tr);
            assert_eq!((k - 1) % p, n);
            let fp = crate::cryptanalysis::ext_curve::ExtPoint::affine(f.frob(pt.x), f.frob(pt.y));
            assert_eq!(SubfieldOrbitKey.key(&inst.curve, &fp), k);
            assert_eq!(SubfieldOrbitKey.key(&inst.curve, &inst.curve.neg(pt)), k);
        }
    }

    /// The same on `F_p` with `m = 2`: `j = 0` keyed by `x³`, `j = 1728`
    /// by `x²`.
    #[test]
    fn the_orbit_table_agrees_with_the_negation_table_on_prime_field_automorphism_curves() {
        for (family, exponent, w) in [(CmFamily::J0, 3u32, 6u32), (CmFamily::J1728, 2, 4)] {
            let inst = generate_cm_instance(family, 18, 3, 8).unwrap();
            let (fb, _) = glv_orbit_base(&inst, 40, AutomorphismGroup::Auto, true).unwrap();
            let gens = automorphism_generators(&inst, AutomorphismGroup::Auto).unwrap();
            let refs: Vec<&dyn Endomorphism<_>> = gens.iter().map(|b| b.as_ref()).collect();
            let g = inst.generator_point();
            let classes = EndomorphismClasses::new(&inst.curve, g, inst.r, refs).unwrap();
            assert_eq!(classes.group_order, w);
            let mut ops = GroupOps::default();
            let target = inst.curve.mul(&mut ops, g, 5);
            let ctx = InstanceCtx {
                group: &inst.curve,
                generator: g,
                target,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: None,
            };
            let mut orbit = OrbitMitmOracle::new(2, PrimePowerKey { exponent }, classes);
            orbit
                .prepare(&ctx, &fb, &Params::default(), &mut ops)
                .unwrap();
            let mut plain = MitmOracle::new(2);
            let mut params = Params::default();
            params.set("negation_folded", "1");
            plain.prepare(&ctx, &fb, &params, &mut ops).unwrap();
            let (mut ca, mut cb) = (OracleCounters::default(), OracleCounters::default());
            let mut hits = 0;
            for k in 2..2000u64 {
                let pt = inst.curve.mul(&mut ops, g, k);
                let a = orbit.decompose(&ctx, &fb, &mut ops, &mut ca, pt);
                let b = plain.decompose(&ctx, &fb, &mut ops, &mut cb, pt);
                assert_eq!(a.is_some(), b.is_some(), "{family:?} k = {k}");
                if let Some(idx) = a {
                    assert_eq!(
                        inst.curve
                            .add(&mut ops, fb.points[idx[0]], fb.points[idx[1]]),
                        pt
                    );
                    hits += 1;
                }
            }
            assert!(hits > 0, "{family:?}");
            assert_eq!(orbit.stats.key_collisions, 0);
            assert_eq!(orbit.stats.image_mismatches, 0);
        }
    }
}
