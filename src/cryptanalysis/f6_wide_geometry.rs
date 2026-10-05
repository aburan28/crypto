//! Exact four-summand geometric closure for two-word binary curves.
//!
//! This is the wide-field pair closure component of F6-IC. It is a
//! meet-in-the-middle method, not a new Gröbner basis algorithm. It keeps
//! field coordinates in two words and amortizes inversions across a fixed
//! point's additions. The index is target-independent; `solve4` does only
//! target-dependent work. Every returned witness is checked in the curve
//! group. A miss is exact for the supplied point list.

use std::collections::HashMap;

use crate::binary_ecc::curve::{point_add, point_neg};
use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
enum PointKey {
    Infinity,
    Affine(u128, u128),
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
enum PointXKey {
    Infinity,
    Affine(u128),
}

fn packed(element: &F2mElement) -> Option<u128> {
    let bits = element.raw_bits();
    if bits.len() > 2 {
        return None;
    }
    Some(u128::from(bits[0]) | (u128::from(*bits.get(1).unwrap_or(&0)) << 64))
}

fn key(point: &BinaryPoint) -> Option<PointKey> {
    match point {
        BinaryPoint::Infinity => Some(PointKey::Infinity),
        BinaryPoint::Affine { x, y } => Some(PointKey::Affine(packed(x)?, packed(y)?)),
    }
}

fn x_key(point: &BinaryPoint) -> Option<PointXKey> {
    match point {
        BinaryPoint::Infinity => Some(PointXKey::Infinity),
        BinaryPoint::Affine { x, .. } => Some(PointXKey::Affine(packed(x)?)),
    }
}

/// Add `fixed` to a batch with one field inversion for all ordinary
/// distinct-x additions. Exceptional cases use the reference group law.
pub fn batch_add_fixed(
    curve: &BinaryCurve,
    fixed: &BinaryPoint,
    others: &[BinaryPoint],
) -> Vec<BinaryPoint> {
    let BinaryPoint::Affine { x: x1, y: y1 } = fixed else {
        return others.to_vec();
    };
    let irr = &curve.irreducible;
    let mut result = vec![BinaryPoint::Infinity; others.len()];
    let mut positions = Vec::new();
    let mut denominators = Vec::new();
    for (i, other) in others.iter().enumerate() {
        match other {
            BinaryPoint::Affine { x: x2, .. } if x1 != x2 => {
                positions.push(i);
                denominators.push(x1.add(x2));
            }
            _ => result[i] = point_add(curve, fixed, other),
        }
    }
    if denominators.is_empty() {
        return result;
    }
    let mut prefix = Vec::with_capacity(denominators.len() + 1);
    prefix.push(F2mElement::one(curve.m));
    for denominator in &denominators {
        prefix.push(prefix.last().unwrap().mul(denominator, irr));
    }
    let mut inverse = prefix
        .last()
        .unwrap()
        .flt_inverse(irr)
        .expect("nonzero product");
    for j in (0..positions.len()).rev() {
        let i = positions[j];
        let denominator_inverse = inverse.mul(&prefix[j], irr);
        inverse = inverse.mul(&denominators[j], irr);
        let BinaryPoint::Affine { x: x2, y: y2 } = &others[i] else {
            unreachable!();
        };
        let lambda = y1.add(y2).mul(&denominator_inverse, irr);
        let x3 = lambda
            .square(irr)
            .add(&lambda)
            .add(x1)
            .add(x2)
            .add(&curve.a);
        let y3 = lambda.mul(&x1.add(&x3), irr).add(&x3).add(y1);
        result[i] = BinaryPoint::Affine { x: x3, y: y3 };
    }
    result
}

/// Compute `fixed + other` and `fixed - other` using the same batch of
/// denominator inverses. On a binary curve the two operands have the same
/// x-coordinate, so their ordinary-addition denominators coincide.
pub fn batch_add_fixed_both_signs(
    curve: &BinaryCurve,
    fixed: &BinaryPoint,
    others: &[BinaryPoint],
) -> (Vec<BinaryPoint>, Vec<BinaryPoint>) {
    let BinaryPoint::Affine { x: x1, y: y1 } = fixed else {
        return (others.to_vec(), others.iter().map(point_neg).collect());
    };
    let irr = &curve.irreducible;
    let mut plus = vec![BinaryPoint::Infinity; others.len()];
    let mut minus = vec![BinaryPoint::Infinity; others.len()];
    let mut positions = Vec::new();
    let mut denominators = Vec::new();
    for (i, other) in others.iter().enumerate() {
        match other {
            BinaryPoint::Affine { x: x2, .. } if x1 != x2 => {
                positions.push(i);
                denominators.push(x1.add(x2));
            }
            _ => {
                plus[i] = point_add(curve, fixed, other);
                minus[i] = point_add(curve, fixed, &point_neg(other));
            }
        }
    }
    if denominators.is_empty() {
        return (plus, minus);
    }
    let mut prefix = Vec::with_capacity(denominators.len() + 1);
    prefix.push(F2mElement::one(curve.m));
    for denominator in &denominators {
        prefix.push(prefix.last().unwrap().mul(denominator, irr));
    }
    let mut inverse = prefix
        .last()
        .unwrap()
        .flt_inverse(irr)
        .expect("nonzero product");
    for j in (0..positions.len()).rev() {
        let i = positions[j];
        let denominator_inverse = inverse.mul(&prefix[j], irr);
        inverse = inverse.mul(&denominators[j], irr);
        let BinaryPoint::Affine { x: x2, y: y2 } = &others[i] else {
            unreachable!();
        };
        let lambda_plus = y1.add(y2).mul(&denominator_inverse, irr);
        let lambda_minus = lambda_plus.add(&x2.mul(&denominator_inverse, irr));
        let finish = |lambda: &F2mElement| {
            let x3 = lambda.square(irr).add(lambda).add(x1).add(x2).add(&curve.a);
            let y3 = lambda.mul(&x1.add(&x3), irr).add(&x3).add(y1);
            BinaryPoint::Affine { x: x3, y: y3 }
        };
        plus[i] = finish(&lambda_plus);
        minus[i] = finish(&lambda_minus);
    }
    (plus, minus)
}

/// Query only the x keys of both sums. The y-coordinate is needed solely
/// when a key is present in the signed pair index, so the caller can defer
/// that calculation to the rare candidate hit.
fn batch_x_keys_fixed_both_signs(
    curve: &BinaryCurve,
    fixed: &BinaryPoint,
    others: &[BinaryPoint],
) -> Option<Vec<(PointXKey, PointXKey)>> {
    let BinaryPoint::Affine { x: x1, y: y1 } = fixed else {
        return others
            .iter()
            .map(|point| {
                let x = x_key(point)?;
                Some((x, x))
            })
            .collect();
    };
    let irr = &curve.irreducible;
    let mut keys = vec![(PointXKey::Infinity, PointXKey::Infinity); others.len()];
    let mut positions = Vec::new();
    let mut denominators = Vec::new();
    for (i, other) in others.iter().enumerate() {
        match other {
            BinaryPoint::Affine { x: x2, .. } if x1 != x2 => {
                positions.push(i);
                denominators.push(x1.add(x2));
            }
            _ => {
                let plus = point_add(curve, fixed, other);
                let minus = point_add(curve, fixed, &point_neg(other));
                keys[i] = (x_key(&plus)?, x_key(&minus)?);
            }
        }
    }
    if denominators.is_empty() {
        return Some(keys);
    }
    let mut prefix = Vec::with_capacity(denominators.len() + 1);
    prefix.push(F2mElement::one(curve.m));
    for denominator in &denominators {
        prefix.push(prefix.last().unwrap().mul(denominator, irr));
    }
    let mut inverse = prefix.last()?.flt_inverse(irr)?;
    for j in (0..positions.len()).rev() {
        let i = positions[j];
        let denominator_inverse = inverse.mul(&prefix[j], irr);
        inverse = inverse.mul(&denominators[j], irr);
        let BinaryPoint::Affine { x: x2, y: y2 } = &others[i] else {
            unreachable!();
        };
        let lambda = y1.add(y2).mul(&denominator_inverse, irr);
        let delta = x2.mul(&denominator_inverse, irr);
        let x_plus = lambda
            .square(irr)
            .add(&lambda)
            .add(x1)
            .add(x2)
            .add(&curve.a);
        // lambda_minus = lambda_plus + delta in characteristic two.
        let x_minus = x_plus.add(&delta.square(irr)).add(&delta);
        keys[i] = (
            PointXKey::Affine(packed(&x_plus)?),
            PointXKey::Affine(packed(&x_minus)?),
        );
    }
    Some(keys)
}

/// The memory cap counts unordered pairs, including repeats. This cap is
/// independent of the target, and a failed construction allocates no table.
pub struct F6WidePairIndex {
    curve: BinaryCurve,
    points: Vec<BinaryPoint>,
    sums: Vec<(BinaryPoint, (usize, usize))>,
    // A group sum determines the residual exactly. One representative pair
    // per sum suffices because repeated factor-base points are allowed.
    lookup: HashMap<PointKey, (usize, usize)>,
}

impl F6WidePairIndex {
    pub fn new(curve: &BinaryCurve, points: &[BinaryPoint], max_pairs: usize) -> Option<Self> {
        if curve.m > 128 || points.is_empty() || points.iter().any(|p| !curve.is_on_curve(p)) {
            return None;
        }
        let pair_count = points.len().checked_mul(points.len().checked_add(1)?)? / 2;
        if pair_count > max_pairs {
            return None;
        }
        let mut sums = Vec::with_capacity(pair_count);
        let mut lookup: HashMap<PointKey, (usize, usize)> = HashMap::with_capacity(pair_count);
        for i in 0..points.len() {
            let row = batch_add_fixed(curve, &points[i], &points[i..]);
            for (offset, sum) in row.into_iter().enumerate() {
                let pair = (i, i + offset);
                lookup.entry(key(&sum)?).or_insert(pair);
                sums.push((sum, pair));
            }
        }
        Some(Self {
            curve: curve.clone(),
            points: points.to_vec(),
            sums,
            lookup,
        })
    }

    pub fn pair_count(&self) -> usize {
        self.sums.len()
    }

    /// Exact four-summand search. The first match is returned in stable
    /// unordered-pair order; repetitions of a factor-base point are allowed.
    pub fn solve4(&self, target: &BinaryPoint) -> Option<[usize; 4]> {
        self.solve4_impl(target, true)
    }

    /// Reference path for paired stage diagnostics. It uses the same pair
    /// index and lookup order, with one affine inversion per residual.
    pub fn solve4_scalar_reference(&self, target: &BinaryPoint) -> Option<[usize; 4]> {
        self.solve4_impl(target, false)
    }

    fn solve4_impl(&self, target: &BinaryPoint, batched: bool) -> Option<[usize; 4]> {
        if !self.curve.is_on_curve(target) {
            return None;
        }
        const BATCH: usize = 256;
        for chunk in self.sums.chunks(BATCH) {
            let negated: Vec<_> = chunk.iter().map(|(sum, _)| point_neg(sum)).collect();
            let residuals = if batched {
                batch_add_fixed(&self.curve, target, &negated)
            } else {
                negated
                    .iter()
                    .map(|point| point_add(&self.curve, target, point))
                    .collect()
            };
            for (residual, (_, (i, j))) in residuals.iter().zip(chunk) {
                if let Some(&(k, l)) = self.lookup.get(&key(residual)?) {
                    let indices = [*i, *j, k, l];
                    let sum = indices.iter().fold(BinaryPoint::Infinity, |acc, &index| {
                        point_add(&self.curve, &acc, &self.points[index])
                    });
                    if sum == *target {
                        return Some(indices);
                    }
                }
            }
        }
        None
    }
}

struct SignedPairSum {
    point: BinaryPoint,
    pair: (usize, usize),
    neg_pair: (usize, usize),
}

/// Exact sign-quotient pair index for a distinct, negation-closed base.
/// Each stored pair sum also represents its negative via the negated source
/// points. A query checks both signs and verifies every returned witness.
pub struct F6SignedPairIndex {
    curve: BinaryCurve,
    points: Vec<BinaryPoint>,
    sums: Vec<SignedPairSum>,
    lookup: HashMap<PointXKey, usize>,
    x_filter: Vec<u64>,
    pair_count: usize,
}

impl F6SignedPairIndex {
    const X_FILTER_BITS: usize = 25;
    const X_FILTER_MASK: usize = (1 << Self::X_FILTER_BITS) - 1;

    pub fn new(curve: &BinaryCurve, points: &[BinaryPoint], max_pairs: usize) -> Option<Self> {
        if curve.m > 128 || points.is_empty() || points.iter().any(|p| !curve.is_on_curve(p)) {
            return None;
        }
        let pair_count = points.len().checked_mul(points.len().checked_add(1)?)? / 2;
        if pair_count > max_pairs {
            return None;
        }
        let point_index: HashMap<_, _> = points
            .iter()
            .enumerate()
            .map(|(i, point)| Some((key(point)?, i)))
            .collect::<Option<_>>()?;
        if point_index.len() != points.len() {
            return None;
        }
        let negatives: Vec<usize> = points
            .iter()
            .map(|point| point_index.get(&key(&point_neg(point))?).copied())
            .collect::<Option<_>>()?;
        let capacity = pair_count.div_ceil(2);
        let mut sums = Vec::with_capacity(capacity);
        let mut lookup = HashMap::with_capacity(capacity);
        let mut x_filter = vec![0u64; 1 << (Self::X_FILTER_BITS - 6)];
        for i in 0..points.len() {
            let row = batch_add_fixed(curve, &points[i], &points[i..]);
            for (offset, sum) in row.into_iter().enumerate() {
                let x = x_key(&sum)?;
                if let std::collections::hash_map::Entry::Vacant(entry) = lookup.entry(x) {
                    if let PointXKey::Affine(value) = x {
                        let bit = (value as usize) & Self::X_FILTER_MASK;
                        x_filter[bit >> 6] |= 1u64 << (bit & 63);
                    }
                    let j = i + offset;
                    entry.insert(sums.len());
                    sums.push(SignedPairSum {
                        point: sum,
                        pair: (i, j),
                        neg_pair: (negatives[i], negatives[j]),
                    });
                }
            }
        }
        Some(Self {
            curve: curve.clone(),
            points: points.to_vec(),
            sums,
            lookup,
            x_filter,
            pair_count,
        })
    }

    pub fn pair_count(&self) -> usize {
        self.pair_count
    }

    pub fn signed_sum_count(&self) -> usize {
        self.sums.len()
    }

    fn may_contain_x(&self, x: PointXKey) -> bool {
        match x {
            PointXKey::Infinity => true,
            PointXKey::Affine(value) => {
                let bit = (value as usize) & Self::X_FILTER_MASK;
                self.x_filter[bit >> 6] & (1u64 << (bit & 63)) != 0
            }
        }
    }

    fn lookup_pair(&self, residual: &BinaryPoint) -> Option<(usize, usize)> {
        let index = *self.lookup.get(&x_key(residual)?)?;
        self.pair_at(index, residual)
    }

    fn pair_at(&self, index: usize, residual: &BinaryPoint) -> Option<(usize, usize)> {
        let entry = &self.sums[index];
        if *residual == entry.point {
            Some(entry.pair)
        } else if *residual == point_neg(&entry.point) {
            Some(entry.neg_pair)
        } else {
            None
        }
    }

    fn verify(
        &self,
        target: &BinaryPoint,
        first: (usize, usize),
        second: (usize, usize),
    ) -> Option<[usize; 4]> {
        let indices = [first.0, first.1, second.0, second.1];
        let sum = indices.iter().fold(BinaryPoint::Infinity, |acc, &index| {
            point_add(&self.curve, &acc, &self.points[index])
        });
        (sum == *target).then_some(indices)
    }

    pub fn solve4(&self, target: &BinaryPoint) -> Option<[usize; 4]> {
        if !self.curve.is_on_curve(target) {
            return None;
        }
        const BATCH: usize = 256;
        for chunk in self.sums.chunks(BATCH) {
            let points: Vec<_> = chunk.iter().map(|entry| entry.point.clone()).collect();
            let (plus, minus) = batch_add_fixed_both_signs(&self.curve, target, &points);
            for ((entry, plus_residual), minus_residual) in
                chunk.iter().zip(plus.iter()).zip(minus.iter())
            {
                if let Some(pair) = self.lookup_pair(minus_residual) {
                    if let Some(witness) = self.verify(target, entry.pair, pair) {
                        return Some(witness);
                    }
                }
                if let Some(pair) = self.lookup_pair(plus_residual) {
                    if let Some(witness) = self.verify(target, entry.neg_pair, pair) {
                        return Some(witness);
                    }
                }
            }
        }
        None
    }

    /// Exact four-summand search with y-coordinate work deferred until an
    /// x-key match. The search order and verification are the same as
    /// `solve4`; only the bulk residual representation changes.
    pub fn solve4_xonly(&self, target: &BinaryPoint) -> Option<[usize; 4]> {
        if !self.curve.is_on_curve(target) {
            return None;
        }
        const BATCH: usize = 256;
        for chunk in self.sums.chunks(BATCH) {
            let points: Vec<_> = chunk.iter().map(|entry| entry.point.clone()).collect();
            let keys = batch_x_keys_fixed_both_signs(&self.curve, target, &points)?;
            for (entry, (plus_key, minus_key)) in chunk.iter().zip(keys) {
                let minus_hit = self
                    .may_contain_x(minus_key)
                    .then(|| self.lookup.get(&minus_key).copied())
                    .flatten();
                if let Some(index) = minus_hit {
                    let residual = point_add(&self.curve, target, &point_neg(&entry.point));
                    if let Some(pair) = self.pair_at(index, &residual) {
                        if let Some(witness) = self.verify(target, entry.pair, pair) {
                            return Some(witness);
                        }
                    }
                }
                let plus_hit = self
                    .may_contain_x(plus_key)
                    .then(|| self.lookup.get(&plus_key).copied())
                    .flatten();
                if let Some(index) = plus_hit {
                    let residual = point_add(&self.curve, target, &entry.point);
                    if let Some(pair) = self.pair_at(index, &residual) {
                        if let Some(witness) = self.verify(target, entry.neg_pair, pair) {
                            return Some(witness);
                        }
                    }
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::curve::scalar_mul;
    use crate::binary_ecc::IrreduciblePoly;
    use crate::cryptanalysis::koblitz_index_calculus::{
        build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
    };
    use num_bigint::BigUint;

    fn n83_curve() -> BinaryCurve {
        let m = 83;
        let parse = |s: &str| BigUint::parse_bytes(s.as_bytes(), 16).unwrap();
        BinaryCurve {
            m,
            irreducible: IrreduciblePoly {
                degree: m,
                low_terms: vec![0, 1, 2, 45],
            },
            a: F2mElement::one(m),
            b: F2mElement::one(m),
            generator: BinaryPoint::Affine {
                x: F2mElement::from_biguint(&parse("68a212cfe19a809fe0598"), m),
                y: F2mElement::from_biguint(&parse("244a245ea0b17d8cc8297"), m),
            },
            order: BigUint::from(8_569_786_107_849_059u64),
            cofactor: BigUint::from(1_128_547_018u64),
        }
    }

    #[test]
    fn n83_batch_add_matches_reference_with_exceptions() {
        let curve = n83_curve();
        assert!(curve.is_on_curve(&curve.generator));
        let points: Vec<_> = (1u32..=12)
            .map(|i| scalar_mul(&curve, &curve.generator, &BigUint::from(i)))
            .collect();
        let fixed = &points[3];
        let mut inputs = points.clone();
        inputs.push(point_neg(fixed));
        inputs.push(BinaryPoint::Infinity);
        assert_eq!(
            batch_add_fixed(&curve, fixed, &inputs),
            inputs
                .iter()
                .map(|p| point_add(&curve, fixed, p))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn n83_signed_batch_add_matches_reference_with_exceptions() {
        let curve = n83_curve();
        let points: Vec<_> = (1u32..=12)
            .map(|i| scalar_mul(&curve, &curve.generator, &BigUint::from(i)))
            .collect();
        let fixed = &points[3];
        let mut inputs = points.clone();
        inputs.push(point_neg(fixed));
        inputs.push(BinaryPoint::Infinity);
        let (plus, minus) = batch_add_fixed_both_signs(&curve, fixed, &inputs);
        let keys = batch_x_keys_fixed_both_signs(&curve, fixed, &inputs).unwrap();
        for ((plus_key, minus_key), (plus_point, minus_point)) in
            keys.iter().zip(plus.iter().zip(minus.iter()))
        {
            assert_eq!(*plus_key, x_key(plus_point).unwrap());
            assert_eq!(*minus_key, x_key(minus_point).unwrap());
        }
        assert_eq!(
            plus,
            inputs
                .iter()
                .map(|p| point_add(&curve, fixed, p))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            minus,
            inputs
                .iter()
                .map(|p| point_add(&curve, fixed, &point_neg(p)))
                .collect::<Vec<_>>()
        );
        let (plus_inf, minus_inf) =
            batch_add_fixed_both_signs(&curve, &BinaryPoint::Infinity, &inputs);
        let keys_inf =
            batch_x_keys_fixed_both_signs(&curve, &BinaryPoint::Infinity, &inputs).unwrap();
        for ((plus_key, minus_key), (plus_point, minus_point)) in
            keys_inf.iter().zip(plus_inf.iter().zip(minus_inf.iter()))
        {
            assert_eq!(*plus_key, x_key(plus_point).unwrap());
            assert_eq!(*minus_key, x_key(minus_point).unwrap());
        }
        assert_eq!(plus_inf, inputs);
        assert_eq!(minus_inf, inputs.iter().map(point_neg).collect::<Vec<_>>());
    }

    #[test]
    fn n83_signed_pair_index_matches_exhaustive_small_base() {
        let curve = n83_curve();
        let positive: Vec<_> = (1u32..=4)
            .map(|i| scalar_mul(&curve, &curve.generator, &BigUint::from(i)))
            .collect();
        assert!(F6SignedPairIndex::new(&curve, &positive, 10).is_none());
        let mut points = positive.clone();
        points.extend(positive.iter().map(point_neg));
        let index = F6SignedPairIndex::new(&curve, &points, 36).unwrap();
        assert_eq!(index.pair_count(), 36);
        assert!(index.signed_sum_count() <= 36);
        assert!(index.lookup.keys().all(|&x| index.may_contain_x(x)));
        let mut four_sums = std::collections::HashSet::new();
        for a in &points {
            for b in &points {
                let ab = point_add(&curve, a, b);
                for c in &points {
                    let abc = point_add(&curve, &ab, c);
                    for d in &points {
                        four_sums.insert(key(&point_add(&curve, &abc, d)).unwrap());
                    }
                }
            }
        }
        for scalar in 0u32..=40 {
            let target = scalar_mul(&curve, &curve.generator, &BigUint::from(scalar));
            let expected = four_sums.contains(&key(&target).unwrap());
            let witness = index.solve4(&target);
            assert_eq!(witness.is_some(), expected, "target scalar {scalar}");
            let xonly_witness = index.solve4_xonly(&target);
            assert_eq!(
                xonly_witness.is_some(),
                expected,
                "x-only target scalar {scalar}"
            );
            if let Some(indices) = xonly_witness {
                let replay = indices.iter().fold(BinaryPoint::Infinity, |acc, &i| {
                    point_add(&curve, &acc, &points[i])
                });
                assert_eq!(replay, target);
            }
            if let Some(indices) = witness {
                let replay = indices.iter().fold(BinaryPoint::Infinity, |acc, &i| {
                    point_add(&curve, &acc, &points[i])
                });
                assert_eq!(replay, target);
            }
        }
    }

    #[test]
    fn n83_pair_closure_matches_exhaustive_group_search() {
        let curve = n83_curve();
        let points: Vec<_> = (1u32..=8)
            .map(|i| scalar_mul(&curve, &curve.generator, &BigUint::from(i)))
            .collect();
        let index = F6WidePairIndex::new(&curve, &points, 36).unwrap();
        assert_eq!(index.pair_count(), 36);
        assert!(F6WidePairIndex::new(&curve, &points, 35).is_none());
        let target = scalar_mul(&curve, &curve.generator, &BigUint::from(23u32));
        let witness = index.solve4(&target).unwrap();
        let sum = witness.iter().fold(BinaryPoint::Infinity, |acc, &i| {
            point_add(&curve, &acc, &points[i])
        });
        assert_eq!(sum, target);
        let missing = scalar_mul(&curve, &curve.generator, &BigUint::from(33u32));
        assert!(index.solve4(&missing).is_none());
    }

    #[test]
    fn n83_registered_curve_and_subgroup_usable_base_close_exactly() {
        let kc = KoblitzCurve::known_n83_k1().unwrap();
        assert_eq!(kc.label(), "icv1-f2m83-t6151469093347-cdcc5432");
        let parent = build_standard_subspace_factor_base(&kc, 5).unwrap();
        let base = cofactor_project_factor_base(&kc, &parent).unwrap();
        assert_eq!(parent.points.len(), 37);
        assert_eq!(base.points.len(), 36);
        assert_eq!(base.signed_orbits.len(), 18);
        for point in &base.points {
            assert_eq!(kc.mul(point, &kc.subgroup_order), BinaryPoint::Infinity);
        }
        let index = F6WidePairIndex::new(&kc.curve, &base.points, 1024).unwrap();
        let planted = base.points[..4]
            .iter()
            .fold(BinaryPoint::Infinity, |sum, point| kc.add(&sum, point));
        let witness = index.solve4(&planted).unwrap();
        let replay = witness.iter().fold(BinaryPoint::Infinity, |sum, &i| {
            kc.add(&sum, &base.points[i])
        });
        assert_eq!(replay, planted);
    }

    #[test]
    fn n83_k0_confidence_gate_generator_and_frobenius_match() {
        let kc = KoblitzCurve::known_n83_k0().unwrap();
        assert_eq!(kc.label(), "icv1-f2m83-tm6151469093347-debefd74");
        assert_eq!(kc.cofactor, BigUint::from(4u32));
        assert_eq!(
            kc.mul(kc.generator(), &kc.subgroup_order),
            BinaryPoint::Infinity
        );
        assert_eq!(
            kc.frobenius(kc.generator()),
            kc.mul(kc.generator(), &kc.lambda)
        );
    }
}
