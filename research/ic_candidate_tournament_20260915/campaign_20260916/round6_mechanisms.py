"""Exact source changes for round 0006, each applied in isolation to the
round-0005 winner (`combined_descent`) and tested for equivalence before any
measurement. Every transform asserts its anchor text occurs exactly once so a
drifted source fails loudly instead of silently building the wrong candidate."""
import re

IC = 'src/cryptanalysis/koblitz_index_calculus.rs'
GF = 'src/cryptanalysis/semaev_decomp.rs'


def _replace_once(text, needle, replacement):
    assert text.count(needle) == 1, f'anchor must occur exactly once: {needle[:60]!r}'
    return text.replace(needle, replacement, 1)


# ── Field inversion (semaev_decomp.rs) ───────────────────────────────

INV_DOC = '''    /// `a^{-1}` by Fermat: `a^(2^n − 2)`.  Zero maps to zero.
    ///
    /// Still `n` squarings and `n` multiplications — which is why the
    /// code above it goes to some length to need only one per *row* of
    /// the pair loop rather than one per polynomial division.
'''
INV_HEAD = INV_DOC + '''    pub fn inv(&self, a: u64) -> u64 {
        if a == 0 {
            return 0;
        }
        let mut result = 1u64;
'''

INV_TEST = '''
#[cfg(test)]
mod tournament_inverse_equivalence {
    use super::*;
    use crate::cryptanalysis::koblitz_index_calculus::KoblitzCurve;

    #[test]
    fn inverse_matches_fermat_and_multiplies_to_one() {
        let mut checked = 0u64;
        for (n, a) in [(5, 0), (7, 1), (13, 0), (17, 1), (19, 0), (19, 1), (23, 0), (31, 1)] {
            let Some(kc) = KoblitzCurve::new(a, n) else { continue };
            let f = Gf2::new(&kc.curve.irreducible);
            let mask = (1u64 << n) - 1;
            let count = if n <= 13 { 1u64 << n } else { 1u64 << 16 };
            for i in 0..count {
                let x = if n <= 13 { i } else { (i.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ (i >> 7)) & mask };
                let inv = f.inv(x);
                assert_eq!(inv, f.inv_fermat(x), "degree {n}, a={a}, x={x}");
                if x == 0 {
                    assert_eq!(inv, 0);
                } else {
                    assert_eq!(f.mul(x, inv), 1, "degree {n}, a={a}, x={x}");
                }
                checked += 1;
            }
        }
        assert!(checked > 100_000);
    }
}
'''


def _inversion(text, body):
    """Install `body` as `inv` and keep the original as `inv_fermat` for the test."""
    text = _replace_once(text, INV_HEAD, body + INV_DOC + '''    /// Retained as the reference for the equivalence test of the inversion above.
    pub fn inv_fermat(&self, a: u64) -> u64 {
        if a == 0 {
            return 0;
        }
        let mut result = 1u64;
''')
    return text + INV_TEST


def it_inv(text):
    return _inversion(text, '''    /// `a^{-1}` as the same power `a^(2^n − 2)` Fermat uses, evaluated with
    /// the Itoh–Tsujii addition chain: `β_k = a^(2^k − 1)`,
    /// `β_{2k} = β_k^(2^k)·β_k`, `β_{k+1} = β_k²·a` along the binary expansion
    /// of `n − 1`, then one squaring. The same `n − 1` squarings, but
    /// `⌊log₂(n−1)⌋ + popcount(n−1) − 1` multiplications instead of `n − 2`,
    /// and no data-dependent branch.
    pub fn inv(&self, a: u64) -> u64 {
        if a == 0 {
            return 0;
        }
        let e = self.n - 1;
        if e == 0 {
            return a;
        }
        let top = 31 - e.leading_zeros();
        let mut beta = a;
        let mut k = 1u32;
        for bit in (0..top).rev() {
            beta = self.mul(self.sqr_k(beta, k), beta);
            k *= 2;
            if (e >> bit) & 1 == 1 {
                beta = self.mul(self.sqr(beta), a);
                k += 1;
            }
        }
        debug_assert_eq!(k, e);
        self.sqr(beta)
    }

''')


def euclid_inv(text):
    return _inversion(text, '''    /// `a^{-1}` by the binary extended Euclidean algorithm over `F_2[z]`
    /// (Hankerson–Menezes–Vanstone, Algorithm 2.48). The invariants
    /// `g1·a ≡ u` and `g2·a ≡ v (mod irr)` hold until `u = 1`; `irr`
    /// carries its leading `z^n` bit, so `v` starts at the modulus itself.
    pub fn inv(&self, a: u64) -> u64 {
        if a == 0 {
            return 0;
        }
        let (mut u, mut v) = (a, self.irr);
        let (mut g1, mut g2) = (1u64, 0u64);
        while u != 1 {
            let mut j = v.leading_zeros() as i32 - u.leading_zeros() as i32;
            if j < 0 {
                std::mem::swap(&mut u, &mut v);
                std::mem::swap(&mut g1, &mut g2);
                j = -j;
            }
            u ^= v << j;
            g1 ^= g2 << j;
        }
        g1
    }

''')


# ── One single-word curve per build / per verified batch ─────────────

def fast_curve_once(text):
    text = _replace_once(text, '''fn factor_base_points_with_x(curve: &BinaryCurve, x: &F2mElement) -> Vec<BinaryPoint> {
    if curve.m % 2 == 0 {
        return points_with_x(curve, x);
    }
    let Some(fc) = FastCurve::new(curve) else {
        return points_with_x(curve, x);
    };
    let f = &fc.field;
''', '''/// The single-word curve the factor-base lift uses, when it applies.
fn factor_base_fast_curve(curve: &BinaryCurve) -> Option<FastCurve> {
    if curve.m % 2 == 0 {
        return None;
    }
    FastCurve::new(curve)
}

fn factor_base_points_with_x(curve: &BinaryCurve, x: &F2mElement) -> Vec<BinaryPoint> {
    match factor_base_fast_curve(curve) {
        Some(fc) => factor_base_points_with_fast(&fc, curve, x),
        None => points_with_x(curve, x),
    }
}

/// [`factor_base_points_with_x`] with the single-word curve built once by
/// the caller instead of once per abscissa.
fn factor_base_points_with_fast(fc: &FastCurve, curve: &BinaryCurve, x: &F2mElement) -> Vec<BinaryPoint> {
    let f = &fc.field;
''')
    text = _replace_once(text, '''    let mut points: Vec<BinaryPoint> = Vec::new();
    for x in &subspace {
        for p in factor_base_points_with_x(&kc.curve, x) {
            points.push(p);
        }
    }
''', '''    let mut points: Vec<BinaryPoint> = Vec::new();
    let lift_curve = factor_base_fast_curve(&kc.curve);
    for x in &subspace {
        let lifts = match &lift_curve {
            Some(fc) => factor_base_points_with_fast(fc, &kc.curve, x),
            None => points_with_x(&kc.curve, x),
        };
        for p in lifts {
            points.push(p);
        }
    }
''')
    text = _replace_once(text, '''            let Some(point) = factor_base_points_with_x(&kc.curve, &x).into_iter().next() else {
                continue;
            };
''', '''            let lifts = if kc.curve.m % 2 == 0 {
                points_with_x(&kc.curve, &x)
            } else {
                factor_base_points_with_fast(&curve, &kc.curve, &x)
            };
            let Some(point) = lifts.into_iter().next() else {
                continue;
            };
''')
    text = _replace_once(text, '''pub fn verify_collected_relation(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    m: usize,
    rel: &CollectedRelation,
) -> bool {
    let r_u64 = kc
''', '''pub fn verify_collected_relation(
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    m: usize,
    rel: &CollectedRelation,
) -> bool {
    verify_collected_relation_with(FastCurve::new(&kc.curve).as_ref(), kc, fb, m, rel)
}

/// [`verify_collected_relation`] with the single-word curve built once by
/// the caller for a whole batch instead of once per relation.
pub fn verify_collected_relation_with(
    fast: Option<&FastCurve>,
    kc: &KoblitzCurve,
    fb: &FrobeniusFactorBase,
    m: usize,
    rel: &CollectedRelation,
) -> bool {
    let r_u64 = kc
''')
    text = _replace_once(text, '''    if let Some(fc) = FastCurve::new(&kc.curve) {
        let sum = rel.points.iter().fold(FastPoint::INFINITY,
''', '''    if let Some(fc) = fast {
        let sum = rel.points.iter().fold(FastPoint::INFINITY,
''')
    text = _replace_once(text, '''        let verdicts: Vec<bool> = relations
            .par_iter()
            .map(|rel| verify_collected_relation(self.kc, self.fb, self.opts.m, rel))
            .collect();
''', '''        let fast = FastCurve::new(&self.kc.curve);
        let verdicts: Vec<bool> = relations
            .par_iter()
            .map(|rel| verify_collected_relation_with(fast.as_ref(), self.kc, self.fb, self.opts.m, rel))
            .collect();
''')
    return text + '''
#[cfg(test)]
mod tournament_fast_curve_once_equivalence {
    use super::*;

    #[test]
    fn shared_curve_lift_and_verification_match_per_call_construction() {
        let mut checked = 0usize;
        let mut relations = 0usize;
        for (n, a) in [(5, 0), (7, 1), (13, 0), (17, 1), (19, 0), (19, 1), (23, 0), (31, 1)] {
            let Some(kc) = KoblitzCurve::new(a, n) else { continue };
            let fc = FastCurve::new(&kc.curve).unwrap();
            let mask = (1u64 << n) - 1;
            for i in 0..2048u64 {
                let raw = i.wrapping_mul(0x9e37_79b9_7f4a_7c15) & mask;
                let x = F2mElement::from_biguint(&BigUint::from(raw), n);
                assert_eq!(factor_base_points_with_fast(&fc, &kc.curve, &x), points_with_x(&kc.curve, &x));
                checked += 1;
            }
            // The smallest cells cannot reach 6n subgroup points; the lift check above still ran.
            let Ok(fb) = build_subgroup_orbit_factor_base(&kc, 43, 6 * n as usize) else { continue };
            let r = kc.subgroup_order.to_u64_digits()[0];
            let g = fc.lift(kc.generator());
            for i in 0..64usize {
                let points = vec![i % fb.points.len(), (7 * i + 1) % fb.points.len(), (13 * i + 2) % fb.points.len()];
                let sum = points.iter().fold(FastPoint::INFINITY, |s, &j| fc.add(s, fc.lift(&fb.points[j])));
                // A true relation, when the sum is a generator multiple, and a forged coefficient.
                let mut candidates = vec![CollectedRelation { trial: i as u64, a: 1 + (i as u64 % (r - 1)), points: points.clone() }];
                if !sum.infinity {
                    let mut acc = g;
                    for coefficient in 1..2048u64 {
                        if acc == sum {
                            candidates.push(CollectedRelation { trial: i as u64, a: coefficient, points: points.clone() });
                            break;
                        }
                        acc = fc.add(acc, g);
                    }
                }
                for rel in &candidates {
                    assert_eq!(verify_collected_relation_with(Some(&fc), &kc, &fb, 3, rel),
                               verify_collected_relation(&kc, &fb, 3, rel));
                    relations += 1;
                }
            }
        }
        assert!(checked > 10_000 && relations > 200, "{checked} lifts, {relations} relations");
    }
}
'''


# ── Lazy symbolic field structure ────────────────────────────────────

FIELD_ACCESSOR = '''    /// The symbolic field structure, built on first use: only the general
    /// (non-pair-table) decomposition path needs it.
    fn field_structure(&self) -> &FieldStructure {
        self.field
            .get_or_init(|| FieldStructure::new(self.kc.n, &self.kc.curve.irreducible))
    }

'''


def lazy_field(text):
    needle = '\n    field: FieldStructure,\n'
    assert text.count(needle) == 2, text.count(needle)
    text = text.replace(needle, '\n    field: std::sync::OnceLock<FieldStructure>,\n')
    needle = '            field: FieldStructure::new(kc.n, &kc.curve.irreducible),\n'
    assert text.count(needle) == 2, text.count(needle)
    text = text.replace(needle, '            field: std::sync::OnceLock::new(),\n')
    uses = re.findall(r'&self\.field,', text)
    assert len(uses) == 3, uses
    text = text.replace('&self.field,', 'self.field_structure(),')
    for anchor in ("impl<'a> RelationCollector<'a> {\n", "impl<'a> IndividualLogSolver<'a> {\n"):
        text = _replace_once(text, anchor, anchor + FIELD_ACCESSOR)
    return text


# ── Frobenius eigenvalue powers computed once per projected map ──────

def lambda_table(text):
    text = _replace_once(text, '''    orbit_of: Vec<Option<(usize, u32, bool)>>,
    representatives: Vec<BinaryPoint>,
}
''', '''    orbit_of: Vec<Option<(usize, u32, bool)>>,
    representatives: Vec<BinaryPoint>,
    /// `λ^k mod r` for `k < n`: the coefficient of a summand `π^k(rep)`,
    /// computed once here rather than by a modular exponentiation per summand.
    lambda_powers: Vec<BigUint>,
}

/// `λ^0, …, λ^(n−1) mod r` by repeated multiplication; equal to
/// `λ.modpow(k, r)` for every `k < n`.
fn lambda_powers(kc: &KoblitzCurve) -> Vec<BigUint> {
    let r = &kc.subgroup_order;
    let mut powers = Vec::with_capacity(kc.n as usize);
    let mut current = BigUint::one() % r;
    for _ in 0..kc.n {
        powers.push(current.clone());
        current = (&current * &kc.lambda) % r;
    }
    powers
}
''')
    text = _replace_once(text, '''    ProjectedSignedOrbitMap {
        orbit_of,
        representatives,
    }
''', '''    ProjectedSignedOrbitMap {
        orbit_of,
        representatives,
        lambda_powers: lambda_powers(kc),
    }
''')
    text = _replace_once(text, '''    Some(ProjectedSignedOrbitMap {
        orbit_of,
        representatives: representatives.into_iter().map(|p| fc.lower(p)).collect(),
    })
''', '''    Some(ProjectedSignedOrbitMap {
        orbit_of,
        representatives: representatives.into_iter().map(|p| fc.lower(p)).collect(),
        lambda_powers: lambda_powers(kc),
    })
''')
    text = _replace_once(text, '''        let mut coeff = kc.lambda.modpow(&BigUint::from(k), r);
''', '''        let mut coeff = match projected_orbits {
            Some(projected) if (k as usize) < projected.lambda_powers.len() => {
                projected.lambda_powers[k as usize].clone()
            }
            _ => kc.lambda.modpow(&BigUint::from(k), r),
        };
''')
    return text + '''
#[cfg(test)]
mod tournament_lambda_table_equivalence {
    use super::*;

    #[test]
    fn lambda_powers_match_modular_exponentiation() {
        let mut checked = 0usize;
        for (n, a) in [(5, 0), (7, 1), (13, 0), (17, 1), (19, 0), (19, 1), (23, 0), (31, 1)] {
            let Some(kc) = KoblitzCurve::new(a, n) else { continue };
            let powers = lambda_powers(&kc);
            assert_eq!(powers.len(), n as usize);
            for (k, power) in powers.iter().enumerate() {
                assert_eq!(*power, kc.lambda.modpow(&BigUint::from(k), &kc.subgroup_order), "degree {n}, k={k}");
                checked += 1;
            }
        }
        assert!(checked > 100);
    }
}
'''


# ── Orbit tables in single-word arithmetic ───────────────────────────

ORBIT_BLOCK_START = '    // Index points for the orbit walk and for relation lookups.\n'
ORBIT_BLOCK_END = '''    if signed_orbits.iter().map(Vec::len).sum::<usize>() != points.len() {
        return None;
    }
'''


def fast_orbits(text):
    start = text.index(ORBIT_BLOCK_START)
    assert text.count(ORBIT_BLOCK_START) == 1
    end = text.index(ORBIT_BLOCK_END, start) + len(ORBIT_BLOCK_END)
    original_block = text[start:end]
    general = original_block.replace('points.len()', 'points.len()')
    dispatch = '''    let (orbits, orbit_of, signed_orbits, signed_orbit_of) = match FastCurve::new(&kc.curve) {
        Some(fc) => orbit_structure_fast(kc, &fc, &points)?,
        None => orbit_structure_general(kc, &points)?,
    };
'''
    text = text[:start] + dispatch + text[end:]
    helpers = '''
/// The Frobenius orbit table, its inverse map, the signed-Frobenius orbit
/// table and its inverse map of a point list.
type OrbitStructure = (Vec<Vec<usize>>, Vec<(usize, u32)>, Vec<Vec<usize>>, Vec<(usize, u32, bool)>);

/// The orbit tables in the general arithmetic: the original construction.
fn orbit_structure_general(kc: &KoblitzCurve, points: &[BinaryPoint]) -> Option<OrbitStructure> {
''' + general + '''    Some((orbits, orbit_of, signed_orbits, signed_orbit_of))
}

/// [`orbit_structure_general`] in single-word arithmetic: the same walks in
/// the same order, keyed by [`FastPoint::pack`], which is injective on
/// affine points and maps `O` to the key no affine point uses.
fn orbit_structure_fast(kc: &KoblitzCurve, fc: &FastCurve, points: &[BinaryPoint]) -> Option<OrbitStructure> {
    let lifted: Vec<FastPoint> = points.iter().map(|p| fc.lift(p)).collect();
    let mut index_of: HashMap<u64, usize> = HashMap::with_capacity(points.len());
    for (i, p) in lifted.iter().enumerate() {
        index_of.insert(p.pack(), i);
    }

    let mut orbit_of: Vec<(usize, u32)> = vec![(usize::MAX, 0); points.len()];
    let mut orbits: Vec<Vec<usize>> = Vec::new();
    for start in 0..points.len() {
        if orbit_of[start].0 != usize::MAX {
            continue;
        }
        let o = orbits.len();
        let mut cycle = Vec::new();
        let mut cur = lifted[start];
        let mut k = 0u32;
        loop {
            let idx = *index_of.get(&cur.pack())?;
            if orbit_of[idx].0 != usize::MAX {
                break;
            }
            orbit_of[idx] = (o, k);
            cycle.push(idx);
            cur = fc.frobenius_k(cur, kc.k);
            k += 1;
        }
        orbits.push(cycle);
    }

    let mut signed_orbit_of = vec![(usize::MAX, 0, false); points.len()];
    let mut signed_orbits = Vec::new();
    for start in 0..points.len() {
        if signed_orbit_of[start].0 != usize::MAX {
            continue;
        }
        let signed_orbit = signed_orbits.len();
        let mut members = Vec::new();
        let mut current = lifted[start];
        for k in 0..kc.n {
            for (negated, point) in [(false, current), (true, fc.neg(current))] {
                let index = *index_of.get(&point.pack())?;
                if signed_orbit_of[index].0 == usize::MAX {
                    signed_orbit_of[index] = (signed_orbit, k, negated);
                    members.push(index);
                } else if signed_orbit_of[index].0 != signed_orbit {
                    return None;
                }
            }
            current = fc.frobenius_k(current, kc.k);
        }
        if current != lifted[start] {
            return None;
        }
        signed_orbits.push(members);
    }
    if signed_orbits.iter().map(Vec::len).sum::<usize>() != points.len() {
        return None;
    }
    Some((orbits, orbit_of, signed_orbits, signed_orbit_of))
}

#[cfg(test)]
mod tournament_fast_orbits_equivalence {
    use super::*;

    #[test]
    fn single_word_orbit_tables_match_general_tables() {
        let mut checked = 0usize;
        for (n, a) in [(5, 0), (7, 1), (13, 0), (17, 1), (19, 0), (19, 1), (23, 0), (31, 1)] {
            let Some(kc) = KoblitzCurve::new(a, n) else { continue };
            let fc = FastCurve::new(&kc.curve).unwrap();
            for (seed, size) in [(43u64, 6 * n as usize), (7, 2 * n as usize), (1, 10 * n as usize)] {
                let Ok(fb) = build_subgroup_orbit_factor_base(&kc, seed, size) else { continue };
                let general = orbit_structure_general(&kc, &fb.points).unwrap();
                let fast = orbit_structure_fast(&kc, &fc, &fb.points).unwrap();
                assert_eq!(general, fast, "degree {n}, a={a}, seed {seed}");
                assert_eq!(fb.orbits, general.0);
                assert_eq!(fb.orbit_of, general.1);
                assert_eq!(fb.signed_orbits, general.2);
                assert_eq!(fb.signed_orbit_of, general.3);
                // A list that is not Frobenius-closed must fail identically.
                let truncated = &fb.points[..fb.points.len() - 1];
                assert_eq!(orbit_structure_general(&kc, truncated).is_none(), orbit_structure_fast(&kc, &fc, truncated).is_none());
                checked += fb.points.len();
            }
        }
        assert!(checked > 1000);
    }
}
'''
    # Place the helpers directly before `point_key`, which follows the function.
    return _replace_once(text, '''/// Hashable identity of a point — the key of
/// [`FrobeniusFactorBase::index_map`].
pub fn point_key(''', helpers + '''
/// Hashable identity of a point — the key of
/// [`FrobeniusFactorBase::index_map`].
pub fn point_key(''')


# ── Descent certificate: the index-calculus admission requirement ─────
#
# Applied to the baseline and therefore to every candidate. Each recovered
# logarithm reports the one relation `[a]G + [b]Q = Σ P_i` it was derived from,
# so the independent checker can confirm that every target's scalar is an
# index-calculus consequence of the factor-base logs, not a generic collision.

WORKER = 'examples/ic_tournament_worker.rs'


def descent_certificate(text):
    text = _replace_once(text, """pub struct IndividualLogReport {
    /// `[a]G + [b]Q` probes drawn before one descended.
    pub trials: usize,
    /// The recovered logarithm, if the descent succeeded and verified.
    pub log: Option<BigUint>,
}
""", """pub struct IndividualLogReport {
    /// `[a]G + [b]Q` probes drawn before one descended.
    pub trials: usize,
    /// The recovered logarithm, if the descent succeeded and verified.
    pub log: Option<BigUint>,
    /// The relation the logarithm was derived from: the certificate that it
    /// is an index-calculus consequence of the factor-base logs.
    pub relation: Option<DescentRelation>,
}

/// The single relation `[a]G + [b]Q = Σ_i P_i` over factor-base points
/// `P_i` (by index) that an individual-logarithm descent recovered its
/// answer from. An empty summand list is the degenerate `[a]G + [b]Q = O`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DescentRelation {
    pub a: u64,
    pub b: u64,
    pub points: Vec<usize>,
}
""")
    text = _replace_once(text, """                if self.kc.mul(self.kc.generator(), &d) == *q { return Some(Some(d)); }
            }
        } else if let Some(idxs) = pair.decompose_fast(initial, m) {
            if let Some(d) = self.logarithm_from(q, &idxs, a0, b) { return Some(Some(d)); }
        }
""", """                if self.kc.mul(self.kc.generator(), &d) == *q {
                    report.relation = Some(DescentRelation { a: a0, b, points: Vec::new() });
                    return Some(Some(d));
                }
            }
        } else if let Some(idxs) = pair.decompose_fast(initial, m) {
            if let Some(d) = self.logarithm_from(q, &idxs, a0, b) {
                report.relation = Some(DescentRelation { a: a0, b, points: idxs });
                return Some(Some(d));
            }
        }
""")
    text = _replace_once(text, """                        if self.kc.mul(self.kc.generator(), &d) == *q {
                            return Some(Some(d));
                        }
""", """                        if self.kc.mul(self.kc.generator(), &d) == *q {
                            report.relation = Some(DescentRelation { a: *a, b: *b, points: Vec::new() });
                            return Some(Some(d));
                        }
""")
    text = _replace_once(text, """                if let Some(d) = self.logarithm_from(q, &idxs, *a, *b) {
                    return Some(Some(d));
                }
""", """                if let Some(d) = self.logarithm_from(q, &idxs, *a, *b) {
                    report.relation = Some(DescentRelation { a: *a, b: *b, points: idxs });
                    return Some(Some(d));
                }
""")
    text = _replace_once(text, """                    if let Some(d) = solve_for_d(&BigUint::from(a), &BigUint::from(b), r) {
                        if self.verify_log(q, &d) {
                            report.log = Some(d.clone());
                            return Some((d, report));
                        }
                    }
                    continue;
""", """                    if let Some(d) = solve_for_d(&BigUint::from(a), &BigUint::from(b), r) {
                        if self.verify_log(q, &d) {
                            report.log = Some(d.clone());
                            report.relation = Some(DescentRelation { a, b, points: Vec::new() });
                            return Some((d, report));
                        }
                    }
                    continue;
""")
    text = _replace_once(text, """            let a = BigUint::from(a);
            let b = BigUint::from(b);
            let relation = relation_from_decomposition_with_mode(
                kc,
                self.fb,
                &idxs,
""", """            let (probe_a, probe_b) = (a, b);
            let a = BigUint::from(a);
            let b = BigUint::from(b);
            let relation = relation_from_decomposition_with_mode(
                kc,
                self.fb,
                &idxs,
""")
    text = _replace_once(text, """            let d = (numerator * hb_inv) % r;
            if self.verify_log(q, &d) {
                report.log = Some(d.clone());
                return Some((d, report));
            }
""", """            let d = (numerator * hb_inv) % r;
            if self.verify_log(q, &d) {
                report.log = Some(d.clone());
                report.relation = Some(DescentRelation { a: probe_a, b: probe_b, points: idxs.clone() });
                return Some((d, report));
            }
""")
    return text


def worker_certificate(text):
    return _replace_once(text, """        solutions.push(json!({"index":i,"recovered":answer.as_ref().map(|(d,_)|d.to_string()),
            "trials":answer.as_ref().map(|(_,r)|r.trials)}));
""", """        solutions.push(json!({"index":i,"recovered":answer.as_ref().map(|(d,_)|d.to_string()),
            "trials":answer.as_ref().map(|(_,r)|r.trials),
            "relation":answer.as_ref().and_then(|(_,r)|r.relation.as_ref())
                .map(|rel|json!({"a":rel.a,"b":rel.b,"points":rel.points}))}));
""")


BASELINE = {IC: [descent_certificate], WORKER: [worker_certificate]}

CANDIDATES = {
    'it_inv': {GF: [it_inv]},
    'euclid_inv': {GF: [euclid_inv]},
    'fast_curve_once': {IC: [fast_curve_once]},
    'lazy_field': {IC: [lazy_field]},
    'lambda_table': {IC: [lambda_table]},
    'fast_orbits': {IC: [fast_orbits]},
    'combined': {GF: [it_inv], IC: [fast_curve_once, lazy_field, lambda_table, fast_orbits]},
}

HYPOTHESES = {
    'it_inv': 'Field inversion by the Itoh–Tsujii chain: the same power as Fermat with 6 instead of 21 multiplications at degree 23; inversion is 41% of descent and 18% of setup instructions.',
    'euclid_inv': 'Field inversion by binary extended Euclid over F_2[z]: fewer instructions than Fermat, at the price of data-dependent branches.',
    'fast_curve_once': 'Build the single-word curve (and its 256-entry reduction tables) once per factor-base build and once per verified relation batch instead of once per abscissa and per relation.',
    'lazy_field': 'Do not build the symbolic FieldStructure the pair-table path never reads; it is constructed twice per job at 1.6M instructions each.',
    'lambda_table': 'Compute the n Frobenius eigenvalue powers once per projected orbit map instead of one BigUint modular exponentiation per relation summand.',
    'fast_orbits': 'Walk the Frobenius and signed-Frobenius orbits of the factor base in single-word arithmetic with packed keys, replacing big-integer squarings and hash keys.',
    'combined': 'Measure it_inv, fast_curve_once, lazy_field, lambda_table and fast_orbits together, every check retained.',
}
