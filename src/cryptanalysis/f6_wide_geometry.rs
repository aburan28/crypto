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

fn stored_point(point: &BinaryPoint) -> Option<Option<(u128, u128)>> {
    match point {
        BinaryPoint::Infinity => Some(None),
        BinaryPoint::Affine { x, y } => Some(Some((packed(x)?, packed(y)?))),
    }
}

fn restore_point(point: Option<(u128, u128)>, m: u32) -> BinaryPoint {
    let Some((x, y)) = point else {
        return BinaryPoint::Infinity;
    };
    BinaryPoint::Affine {
        x: F2mElement::from_words(&[x as u64, (x >> 64) as u64], m),
        y: F2mElement::from_words(&[y as u64, (y >> 64) as u64], m),
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

#[cfg(target_arch = "aarch64")]
mod packed83 {
    use super::*;
    use num_bigint::BigUint;

    pub(super) const MASK: u128 = (1u128 << 83) - 1;

    #[inline(always)]
    fn fold(value: u128) -> u128 {
        let high = value >> 83;
        (value & MASK) ^ high ^ (high << 1) ^ (high << 2) ^ (high << 45)
    }

    #[inline(always)]
    fn reduce(low: u128, high: u128) -> u128 {
        let upper = (low >> 83) ^ (high << 45);
        let first = (low & MASK) ^ upper ^ (upper << 1) ^ (upper << 2) ^ (upper << 45);
        let result = fold(fold(first));
        debug_assert_eq!(result & !MASK, 0);
        result
    }

    #[target_feature(enable = "aes")]
    unsafe fn clmul(a: u64, b: u64) -> u128 {
        std::arch::aarch64::vmull_p64(a, b)
    }

    /// Carry-less multiply modulo the pinned degree-83 polynomial.
    ///
    /// # Safety
    /// The caller must have checked ARM64 AES/PMULL support.
    #[target_feature(enable = "aes")]
    pub(super) unsafe fn mul(a: u128, b: u128) -> u128 {
        debug_assert_eq!((a | b) & !MASK, 0);
        let a0 = a as u64;
        let a1 = (a >> 64) as u64;
        let b0 = b as u64;
        let b1 = (b >> 64) as u64;
        let low_product = unsafe { clmul(a0, b0) };
        let cross = unsafe { clmul(a0, b1) ^ clmul(a1, b0) };
        let high_product = unsafe { clmul(a1, b1) };
        reduce(low_product ^ (cross << 64), (cross >> 64) ^ high_product)
    }

    /// Field squaring uses only two carry-less limb products.
    ///
    /// # Safety
    /// The caller must have checked ARM64 AES/PMULL support.
    #[target_feature(enable = "aes")]
    pub(super) unsafe fn square(a: u128) -> u128 {
        debug_assert_eq!(a & !MASK, 0);
        reduce(unsafe { clmul(a as u64, a as u64) }, unsafe {
            clmul((a >> 64) as u64, (a >> 64) as u64)
        })
    }

    /// # Safety
    /// The caller must have checked ARM64 AES/PMULL support and the pinned
    /// degree-83 irreducible polynomial.
    #[target_feature(enable = "aes")]
    pub(super) unsafe fn batch_x_keys(
        curve: &BinaryCurve,
        fixed: &BinaryPoint,
        others: &[BinaryPoint],
    ) -> Option<Vec<(PointXKey, PointXKey)>> {
        let BinaryPoint::Affine { x: x1, y: y1 } = fixed else {
            return batch_x_keys_fixed_both_signs(curve, fixed, others);
        };
        let x1 = packed(x1)?;
        let y1 = packed(y1)?;
        let a = packed(&curve.a)?;
        let mut keys = vec![(PointXKey::Infinity, PointXKey::Infinity); others.len()];
        let mut positions = Vec::new();
        let mut denominators = Vec::new();
        for (i, other) in others.iter().enumerate() {
            match other {
                BinaryPoint::Affine { x: x2, .. } if x1 != packed(x2)? => {
                    positions.push(i);
                    denominators.push(x1 ^ packed(x2)?);
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
        prefix.push(1u128);
        for &denominator in &denominators {
            prefix.push(unsafe { mul(*prefix.last()?, denominator) });
        }
        let product = F2mElement::from_biguint(&BigUint::from(*prefix.last()?), 83);
        let mut inverse = packed(&product.flt_inverse(&curve.irreducible)?)?;
        for j in (0..positions.len()).rev() {
            let i = positions[j];
            let denominator_inverse = unsafe { mul(inverse, prefix[j]) };
            inverse = unsafe { mul(inverse, denominators[j]) };
            let BinaryPoint::Affine { x: x2, y: y2 } = &others[i] else {
                unreachable!();
            };
            let x2 = packed(x2)?;
            let lambda = unsafe { mul(y1 ^ packed(y2)?, denominator_inverse) };
            let delta = unsafe { mul(x2, denominator_inverse) };
            let x_plus = unsafe { square(lambda) } ^ lambda ^ x1 ^ x2 ^ a;
            let x_minus = x_plus ^ unsafe { square(delta) } ^ delta;
            keys[i] = (PointXKey::Affine(x_plus), PointXKey::Affine(x_minus));
        }
        Some(keys)
    }

    /// The signed index stores affine coordinates as raw field words. This
    /// query path avoids rebuilding heap-backed points for every batch.
    /// # Safety
    /// The caller must check the pinned polynomial and ARM64 AES/PMULL.
    #[target_feature(enable = "aes")]
    pub(super) unsafe fn batch_x_keys_stored(
        curve: &BinaryCurve,
        fixed: &BinaryPoint,
        others: &[SignedPairSum],
    ) -> Option<Vec<(PointXKey, PointXKey)>> {
        let BinaryPoint::Affine { x: x1, y: y1 } = fixed else {
            return others
                .iter()
                .map(|entry| {
                    let x = entry
                        .point
                        .map_or(PointXKey::Infinity, |(x, _)| PointXKey::Affine(x));
                    Some((x, x))
                })
                .collect();
        };
        let x1 = packed(x1)?;
        let y1 = packed(y1)?;
        let a = packed(&curve.a)?;
        let mut keys = vec![(PointXKey::Infinity, PointXKey::Infinity); others.len()];
        let mut positions = Vec::new();
        let mut denominators = Vec::new();
        for (i, entry) in others.iter().enumerate() {
            match entry.point {
                Some((x2, _)) if x1 != x2 => {
                    positions.push(i);
                    denominators.push(x1 ^ x2);
                }
                _ => {
                    let other = restore_point(entry.point, 83);
                    let plus = point_add(curve, fixed, &other);
                    let minus = point_add(curve, fixed, &point_neg(&other));
                    keys[i] = (x_key(&plus)?, x_key(&minus)?);
                }
            }
        }
        if denominators.is_empty() {
            return Some(keys);
        }
        let mut prefix = Vec::with_capacity(denominators.len() + 1);
        prefix.push(1u128);
        for &denominator in &denominators {
            prefix.push(unsafe { mul(*prefix.last()?, denominator) });
        }
        let product = *prefix.last()?;
        let product = F2mElement::from_words(&[product as u64, (product >> 64) as u64], 83);
        let mut inverse = packed(&product.flt_inverse(&curve.irreducible)?)?;
        for j in (0..positions.len()).rev() {
            let i = positions[j];
            let denominator_inverse = unsafe { mul(inverse, prefix[j]) };
            inverse = unsafe { mul(inverse, denominators[j]) };
            let (x2, y2) = others[i].point.expect("ordinary addition");
            let lambda = unsafe { mul(y1 ^ y2, denominator_inverse) };
            let delta = unsafe { mul(x2, denominator_inverse) };
            let x_plus = unsafe { square(lambda) } ^ lambda ^ x1 ^ x2 ^ a;
            let x_minus = x_plus ^ unsafe { square(delta) } ^ delta;
            keys[i] = (PointXKey::Affine(x_plus), PointXKey::Affine(x_minus));
        }
        Some(keys)
    }

    /// # Safety
    /// The caller must have checked ARM64 AES/PMULL support and the pinned
    /// degree-83 irreducible polynomial.
    #[target_feature(enable = "aes")]
    pub(super) unsafe fn batch_add_fixed(
        curve: &BinaryCurve,
        fixed: &BinaryPoint,
        others: &[BinaryPoint],
    ) -> Option<Vec<BinaryPoint>> {
        let BinaryPoint::Affine { x: x1, y: y1 } = fixed else {
            return Some(others.to_vec());
        };
        let x1 = packed(x1)?;
        let y1 = packed(y1)?;
        let a = packed(&curve.a)?;
        let mut result = vec![BinaryPoint::Infinity; others.len()];
        let mut positions = Vec::new();
        let mut denominators = Vec::new();
        for (i, other) in others.iter().enumerate() {
            match other {
                BinaryPoint::Affine { x: x2, .. } if x1 != packed(x2)? => {
                    positions.push(i);
                    denominators.push(x1 ^ packed(x2)?);
                }
                _ => result[i] = point_add(curve, fixed, other),
            }
        }
        if denominators.is_empty() {
            return Some(result);
        }
        let mut prefix = Vec::with_capacity(denominators.len() + 1);
        prefix.push(1u128);
        for &denominator in &denominators {
            prefix.push(unsafe { mul(*prefix.last()?, denominator) });
        }
        let product = *prefix.last()?;
        let product = F2mElement::from_words(&[product as u64, (product >> 64) as u64], 83);
        let mut inverse = packed(&product.flt_inverse(&curve.irreducible)?)?;
        for j in (0..positions.len()).rev() {
            let i = positions[j];
            let denominator_inverse = unsafe { mul(inverse, prefix[j]) };
            inverse = unsafe { mul(inverse, denominators[j]) };
            let BinaryPoint::Affine { x: x2, y: y2 } = &others[i] else {
                unreachable!();
            };
            let x2 = packed(x2)?;
            let lambda = unsafe { mul(y1 ^ packed(y2)?, denominator_inverse) };
            let x3 = unsafe { square(lambda) } ^ lambda ^ x1 ^ x2 ^ a;
            let y3 = unsafe { mul(lambda, x1 ^ x3) } ^ x3 ^ y1;
            result[i] = BinaryPoint::Affine {
                x: F2mElement::from_words(&[x3 as u64, (x3 >> 64) as u64], 83),
                y: F2mElement::from_words(&[y3 as u64, (y3 >> 64) as u64], 83),
            };
        }
        Some(result)
    }
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
    // No heap allocations per sum. `None` is the identity and is distinct
    // from every affine point even when the curve admits (0, 0).
    point: Option<(u128, u128)>,
    pair: (u32, u32),
    neg_pair: (u32, u32),
}

/// Exact sign-quotient pair index for a distinct, negation-closed base.
/// Each stored pair sum also represents its negative via the negated source
/// points. A query checks both signs and verifies every returned witness.
pub struct F6SignedPairIndex {
    curve: BinaryCurve,
    points: Vec<BinaryPoint>,
    sums: Vec<SignedPairSum>,
    lookup: HashMap<u128, usize>,
    infinity_index: Option<usize>,
    pair_count: usize,
}

impl F6SignedPairIndex {
    pub fn new(curve: &BinaryCurve, points: &[BinaryPoint], max_pairs: usize) -> Option<Self> {
        if curve.m > 128 || points.is_empty() || points.iter().any(|p| !curve.is_on_curve(p)) {
            return None;
        }
        u32::try_from(points.len()).ok()?;
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
        let negatives: Vec<u32> = points
            .iter()
            .map(|point| u32::try_from(*point_index.get(&key(&point_neg(point))?)?).ok())
            .collect::<Option<_>>()?;
        let capacity = pair_count.div_ceil(2);
        let mut sums = Vec::with_capacity(capacity);
        let mut lookup = HashMap::with_capacity(capacity);
        let mut infinity_index = None;
        #[cfg(target_arch = "aarch64")]
        let use_packed = curve.m == 83
            && curve.irreducible.degree == 83
            && curve.irreducible.low_terms == [0, 1, 2, 45]
            && std::arch::is_aarch64_feature_detected!("aes");
        for i in 0..points.len() {
            #[cfg(target_arch = "aarch64")]
            let row = if use_packed {
                // SAFETY: the exact polynomial and ARM64 AES/PMULL feature
                // were checked once above; all inputs were curve-checked.
                unsafe { packed83::batch_add_fixed(curve, &points[i], &points[i..]) }?
            } else {
                batch_add_fixed(curve, &points[i], &points[i..])
            };
            #[cfg(not(target_arch = "aarch64"))]
            let row = batch_add_fixed(curve, &points[i], &points[i..]);
            for (offset, sum) in row.into_iter().enumerate() {
                let first_of_x_orbit = match x_key(&sum)? {
                    PointXKey::Infinity => {
                        if infinity_index.is_none() {
                            infinity_index = Some(sums.len());
                            true
                        } else {
                            false
                        }
                    }
                    PointXKey::Affine(x) => {
                        if let std::collections::hash_map::Entry::Vacant(entry) = lookup.entry(x) {
                            entry.insert(sums.len());
                            true
                        } else {
                            false
                        }
                    }
                };
                if first_of_x_orbit {
                    let j = i + offset;
                    sums.push(SignedPairSum {
                        point: stored_point(&sum)?,
                        pair: (i as u32, j as u32),
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
            infinity_index,
            pair_count,
        })
    }

    pub fn pair_count(&self) -> usize {
        self.pair_count
    }

    pub fn signed_sum_count(&self) -> usize {
        self.sums.len()
    }

    fn lookup_index(&self, x: PointXKey) -> Option<usize> {
        match x {
            PointXKey::Infinity => self.infinity_index,
            PointXKey::Affine(x) => self.lookup.get(&x).copied(),
        }
    }

    fn lookup_pair(&self, residual: &BinaryPoint) -> Option<(usize, usize)> {
        let index = self.lookup_index(x_key(residual)?)?;
        self.pair_at(index, residual)
    }

    fn pair_at(&self, index: usize, residual: &BinaryPoint) -> Option<(usize, usize)> {
        let entry = &self.sums[index];
        let residual = stored_point(residual)?;
        if residual == entry.point {
            Some((entry.pair.0 as usize, entry.pair.1 as usize))
        } else if residual == entry.point.map(|(x, y)| (x, x ^ y)) {
            Some((entry.neg_pair.0 as usize, entry.neg_pair.1 as usize))
        } else {
            None
        }
    }

    fn verify(
        &self,
        target: &BinaryPoint,
        first: (u32, u32),
        second: (usize, usize),
    ) -> Option<[usize; 4]> {
        let indices = [first.0 as usize, first.1 as usize, second.0, second.1];
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
            let points: Vec<_> = chunk
                .iter()
                .map(|entry| restore_point(entry.point, self.curve.m))
                .collect();
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
            let points: Vec<_> = chunk
                .iter()
                .map(|entry| restore_point(entry.point, self.curve.m))
                .collect();
            let keys = batch_x_keys_fixed_both_signs(&self.curve, target, &points)?;
            for (entry, (plus_key, minus_key)) in chunk.iter().zip(keys) {
                if let Some(index) = self.lookup_index(minus_key) {
                    let point = restore_point(entry.point, self.curve.m);
                    let residual = point_add(&self.curve, target, &point_neg(&point));
                    if let Some(pair) = self.pair_at(index, &residual) {
                        if let Some(witness) = self.verify(target, entry.pair, pair) {
                            return Some(witness);
                        }
                    }
                }
                if let Some(index) = self.lookup_index(plus_key) {
                    let point = restore_point(entry.point, self.curve.m);
                    let residual = point_add(&self.curve, target, &point);
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

    /// ARM64 PMULL residual arithmetic for the pinned degree-83 field.
    /// All other curves and CPUs use the exact x-only reference path.
    pub fn solve4_xonly_pmull83(&self, target: &BinaryPoint) -> Option<[usize; 4]> {
        #[cfg(not(target_arch = "aarch64"))]
        {
            self.solve4_xonly(target)
        }
        #[cfg(target_arch = "aarch64")]
        {
            if self.curve.m != 83
                || self.curve.irreducible.degree != 83
                || self.curve.irreducible.low_terms != [0, 1, 2, 45]
                || !std::arch::is_aarch64_feature_detected!("aes")
            {
                return self.solve4_xonly(target);
            }
            if !self.curve.is_on_curve(target) {
                return None;
            }
            const BATCH: usize = 256;
            for chunk in self.sums.chunks(BATCH) {
                // SAFETY: the field polynomial and ARM64 AES/PMULL feature
                // were checked above. The result is verified in the group.
                let keys = unsafe { packed83::batch_x_keys_stored(&self.curve, target, chunk) }?;
                for (entry, (plus_key, minus_key)) in chunk.iter().zip(keys) {
                    if let Some(index) = self.lookup_index(minus_key) {
                        let point = restore_point(entry.point, self.curve.m);
                        let residual = point_add(&self.curve, target, &point_neg(&point));
                        if let Some(pair) = self.pair_at(index, &residual) {
                            if let Some(witness) = self.verify(target, entry.pair, pair) {
                                return Some(witness);
                            }
                        }
                    }
                    if let Some(index) = self.lookup_index(plus_key) {
                        let point = restore_point(entry.point, self.curve.m);
                        let residual = point_add(&self.curve, target, &point);
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

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn n83_pmull_field_matches_reference() {
        if !std::arch::is_aarch64_feature_detected!("aes") {
            return;
        }
        let curve = n83_curve();
        let check = |a: u128, b: u128| {
            let lhs = F2mElement::from_biguint(&BigUint::from(a), 83);
            let rhs = F2mElement::from_biguint(&BigUint::from(b), 83);
            // SAFETY: ARM64 AES/PMULL support was checked above.
            assert_eq!(
                unsafe { packed83::mul(a, b) },
                packed(&lhs.mul(&rhs, &curve.irreducible)).unwrap()
            );
            assert_eq!(
                unsafe { packed83::square(a) },
                packed(&lhs.square(&curve.irreducible)).unwrap()
            );
        };
        let edges = [
            0,
            1,
            1 << 45,
            1 << 64,
            1 << 82,
            (1 << 64) - 1,
            packed83::MASK,
        ];
        for a in edges {
            for b in edges {
                check(a, b);
            }
        }
        let mut state = 0x93cd_b69f_47e2_180d_a54b_8c31_f739_62e1u128;
        for _ in 0..10_000 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let a = state & packed83::MASK;
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            check(a, state & packed83::MASK);
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn n83_pmull_batch_keys_match_reference_with_exceptions() {
        if !std::arch::is_aarch64_feature_detected!("aes") {
            return;
        }
        let curve = n83_curve();
        let points: Vec<_> = (1u32..=12)
            .map(|i| scalar_mul(&curve, &curve.generator, &BigUint::from(i)))
            .collect();
        let fixed = &points[3];
        let mut inputs = points.clone();
        inputs.push(point_neg(fixed));
        inputs.push(BinaryPoint::Infinity);
        inputs.push(fixed.clone());
        let expected = batch_x_keys_fixed_both_signs(&curve, fixed, &inputs).unwrap();
        // SAFETY: ARM64 AES/PMULL support and the pinned polynomial were checked.
        let observed = unsafe { packed83::batch_x_keys(&curve, fixed, &inputs) }.unwrap();
        assert_eq!(observed, expected);
        let stored: Vec<_> = inputs
            .iter()
            .map(|point| SignedPairSum {
                point: stored_point(point).unwrap(),
                pair: (0, 0),
                neg_pair: (0, 0),
            })
            .collect();
        let observed_stored =
            unsafe { packed83::batch_x_keys_stored(&curve, fixed, &stored) }.unwrap();
        assert_eq!(observed_stored, expected);
        let expected_inf =
            batch_x_keys_fixed_both_signs(&curve, &BinaryPoint::Infinity, &inputs).unwrap();
        let observed_inf =
            unsafe { packed83::batch_x_keys(&curve, &BinaryPoint::Infinity, &inputs) }.unwrap();
        assert_eq!(observed_inf, expected_inf);
        let observed_stored_inf =
            unsafe { packed83::batch_x_keys_stored(&curve, &BinaryPoint::Infinity, &stored) }
                .unwrap();
        assert_eq!(observed_stored_inf, expected_inf);
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
        #[cfg(target_arch = "aarch64")]
        if std::arch::is_aarch64_feature_detected!("aes") {
            // SAFETY: the pinned polynomial and CPU feature were checked.
            let packed = unsafe { packed83::batch_add_fixed(&curve, fixed, &inputs) }.unwrap();
            assert_eq!(packed, batch_add_fixed(&curve, fixed, &inputs));
            let packed_inf =
                unsafe { packed83::batch_add_fixed(&curve, &BinaryPoint::Infinity, &inputs) }
                    .unwrap();
            assert_eq!(packed_inf, inputs);
        }
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn n83_pmull_pair_build_matches_k0_reference() {
        if !std::arch::is_aarch64_feature_detected!("aes") {
            return;
        }
        let kc = KoblitzCurve::known_n83_k0().unwrap();
        let parent = build_standard_subspace_factor_base(&kc, 8).unwrap();
        let base = cofactor_project_factor_base(&kc, &parent).unwrap();
        assert_eq!(base.points.len(), 258);
        for i in [0usize, 1, 2, 7, 64, 127, 257] {
            let reference = batch_add_fixed(&kc.curve, &base.points[i], &base.points[i..]);
            // SAFETY: the pinned polynomial and CPU feature were checked.
            let candidate =
                unsafe { packed83::batch_add_fixed(&kc.curve, &base.points[i], &base.points[i..]) }
                    .unwrap();
            assert_eq!(candidate, reference, "base row {i}");
        }
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
            let pmull_witness = index.solve4_xonly_pmull83(&target);
            assert_eq!(
                pmull_witness.is_some(),
                expected,
                "PMULL target scalar {scalar}"
            );
            if let Some(indices) = pmull_witness {
                let replay = indices.iter().fold(BinaryPoint::Infinity, |acc, &i| {
                    point_add(&curve, &acc, &points[i])
                });
                assert_eq!(replay, target);
            }
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
