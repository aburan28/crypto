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

/// The memory cap counts unordered pairs, including repeats. This cap is
/// independent of the target, and a failed construction allocates no table.
pub struct F6WidePairIndex {
    curve: BinaryCurve,
    points: Vec<BinaryPoint>,
    sums: Vec<(BinaryPoint, (usize, usize))>,
    lookup: HashMap<PointKey, Vec<(usize, usize)>>,
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
        let mut lookup: HashMap<PointKey, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..points.len() {
            let row = batch_add_fixed(curve, &points[i], &points[i..]);
            for (offset, sum) in row.into_iter().enumerate() {
                let pair = (i, i + offset);
                lookup.entry(key(&sum)?).or_default().push(pair);
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
                if let Some(pairs) = self.lookup.get(&key(residual)?) {
                    for &(k, l) in pairs {
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
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::curve::scalar_mul;
    use crate::binary_ecc::IrreduciblePoly;
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
}
