//! Complete fixed-base signed point-sum lookup for small binary curves.
//!
//! This is an exact group-law oracle. It enumerates all signed multisets of
//! at most three factors, then meets two such halves for at most six factors
//! (or a two-factor left half for the five-factor control). It makes no
//! assumption about the distribution of sums or the availability of an
//! endomorphism on the target curve.

use std::collections::HashMap;

use super::koblitz_fast::{BatchScratch, FastCurve, FastPoint};

const QUERY_BATCH: usize = 4096;

#[derive(Clone, Copy, Debug)]
struct HalfSum {
    indices: [u8; 3],
    len: u8,
}

/// Counts logical point-addition requests, including batch additions whose
/// operands are degenerate. These are *not* field-operation equivalents.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize)]
pub struct MitmCounts {
    pub raw_half_entries: u64,
    pub distinct_half_sums: u64,
    pub build_adds: u64,
    pub query_adds: u64,
    pub query_batches: u64,
    pub lookups: u64,
    pub witness_adds: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize)]
pub struct SignedWitness {
    /// Net coefficients in the order of the supplied base points.
    pub coefficients: Vec<i8>,
    /// Indices into `[+B0,-B0,+B1,-B1,...]`, before cancellation.
    pub signed_indices: Vec<u8>,
    pub left_len: u8,
    pub right_len: u8,
}

impl SignedWitness {
    pub fn l1_norm(&self) -> usize {
        self.coefficients
            .iter()
            .map(|&coefficient| coefficient.unsigned_abs() as usize)
            .sum()
    }
}

/// A complete signed `3+3` lookup over one fixed, ordered factor base.
pub struct NativeSignedMitm<'a> {
    curve: &'a FastCurve,
    signed: Vec<FastPoint>,
    halves: Vec<HalfSum>,
    negative_halves: Vec<FastPoint>,
    up_to_two: usize,
    by_point: HashMap<FastPoint, usize>,
    build_counts: MitmCounts,
}

impl<'a> NativeSignedMitm<'a> {
    pub fn build(curve: &'a FastCurve, base: &[FastPoint]) -> Result<Self, String> {
        if base.is_empty() || base.len() > 127 {
            return Err("base must contain 1..=127 sign classes".into());
        }
        for (i, &point) in base.iter().enumerate() {
            if point.infinity || !curve.is_on_curve(point) {
                return Err(format!("base point {i} is infinity or off curve"));
            }
            if base[..i].iter().any(|other| other.x == point.x) {
                return Err(format!("base point {i} duplicates a sign class"));
            }
        }
        let mut signed = Vec::with_capacity(2 * base.len());
        for &point in base {
            signed.push(point);
            signed.push(curve.neg(point));
        }
        let mut table = Self {
            curve,
            signed,
            halves: Vec::new(),
            negative_halves: Vec::new(),
            up_to_two: 0,
            by_point: HashMap::new(),
            build_counts: MitmCounts::default(),
        };
        table.insert(FastPoint::INFINITY, [0; 3], 0);
        for i in 0..table.signed.len() {
            table.insert(table.signed[i], [i as u8, 0, 0], 1);
        }
        let mut pairs = Vec::new();
        for i in 0..table.signed.len() {
            for j in i..table.signed.len() {
                let point = curve.add(table.signed[i], table.signed[j]);
                table.build_counts.build_adds += 1;
                table.insert(point, [i as u8, j as u8, 0], 2);
                pairs.push((i, j, point));
            }
        }
        table.up_to_two = table.halves.len();
        for (i, j, pair) in pairs {
            for k in j..table.signed.len() {
                let point = curve.add(pair, table.signed[k]);
                table.build_counts.build_adds += 1;
                table.insert(point, [i as u8, j as u8, k as u8], 3);
            }
        }
        table.build_counts.distinct_half_sums = table.halves.len() as u64;
        Ok(table)
    }

    fn insert(&mut self, point: FastPoint, indices: [u8; 3], len: u8) {
        self.build_counts.raw_half_entries += 1;
        if self.by_point.contains_key(&point) {
            return;
        }
        self.by_point.insert(point, self.halves.len());
        self.halves.push(HalfSum { indices, len });
        self.negative_halves.push(self.curve.neg(point));
    }

    pub fn build_counts(&self) -> MitmCounts {
        self.build_counts
    }

    pub fn base_size(&self) -> usize {
        self.signed.len() / 2
    }

    /// `left_max=2` checks all sums of at most five signed factors;
    /// `left_max=3` checks all sums of at most six. A `None` result has
    /// examined the complete corresponding left-hand set.
    pub fn query(
        &self,
        target: FastPoint,
        left_max: u8,
    ) -> Result<(Option<SignedWitness>, MitmCounts), String> {
        if !matches!(left_max, 2 | 3) || !self.curve.is_on_curve(target) {
            return Err("query needs a curve point and left_max=2 or 3".into());
        }
        let limit = if left_max == 2 {
            self.up_to_two
        } else {
            self.halves.len()
        };
        let mut counts = MitmCounts::default();
        let mut out = Vec::with_capacity(QUERY_BATCH);
        let mut scratch = BatchScratch::default();
        for start in (0..limit).step_by(QUERY_BATCH) {
            let end = (start + QUERY_BATCH).min(limit);
            out.clear();
            self.curve.add_many(
                target,
                &self.negative_halves[start..end],
                &mut out,
                &mut scratch,
            );
            if out.len() != end - start {
                return Err("batched addition returned the wrong length".into());
            }
            counts.query_adds += out.len() as u64;
            counts.query_batches += 1;
            for (offset, &remainder) in out.iter().enumerate() {
                counts.lookups += 1;
                let Some(&right) = self.by_point.get(&remainder) else {
                    continue;
                };
                let left = self.halves[start + offset];
                let right = self.halves[right];
                let mut signed_indices = Vec::with_capacity((left.len + right.len) as usize);
                signed_indices.extend_from_slice(&left.indices[..left.len as usize]);
                signed_indices.extend_from_slice(&right.indices[..right.len as usize]);
                let mut sum = FastPoint::INFINITY;
                let mut coefficients = vec![0i8; self.base_size()];
                for &index in &signed_indices {
                    let index = index as usize;
                    sum = self.curve.add(sum, self.signed[index]);
                    counts.witness_adds += 1;
                    coefficients[index / 2] += if index.is_multiple_of(2) { 1 } else { -1 };
                }
                if sum != target {
                    return Err("point-sum witness failed full group-law replay".into());
                }
                let witness = SignedWitness {
                    coefficients,
                    signed_indices,
                    left_len: left.len,
                    right_len: right.len,
                };
                if witness.l1_norm() > (left_max + 3) as usize {
                    return Err("witness exceeds claimed signed arity".into());
                }
                return Ok((Some(witness), counts));
            }
        }
        Ok((None, counts))
    }

    /// Independently reconstruct a reported coefficient vector from the
    /// supplied base. This does not use the half-sum table.
    pub fn verify_coefficients(&self, target: FastPoint, coefficients: &[i8]) -> bool {
        if coefficients.len() != self.base_size() {
            return false;
        }
        let mut sum = FastPoint::INFINITY;
        for (i, &coefficient) in coefficients.iter().enumerate() {
            let point = if coefficient < 0 {
                self.signed[2 * i + 1]
            } else {
                self.signed[2 * i]
            };
            for _ in 0..coefficient.unsigned_abs() {
                sum = self.curve.add(sum, point);
            }
        }
        sum == target
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;
    use crate::binary_ecc::IrreduciblePoly;
    use crate::cryptanalysis::binary_velu::Curve;

    #[test]
    fn small_field_exhaustive_membership_agrees_with_direct_sums() {
        let curve = Curve::new(8, &IrreduciblePoly::deg_8(), 0, 1).expect("small curve");
        let mut base = Vec::new();
        for x in 1..256 {
            if let Some(point) = curve.points_with_x(x).into_iter().next() {
                base.push(point);
            }
            if base.len() == 3 {
                break;
            }
        }
        assert_eq!(base.len(), 3);
        let oracle = NativeSignedMitm::build(&curve.fast, &base).expect("table");
        assert_eq!(oracle.build_counts().raw_half_entries, 84);
        let signed: Vec<_> = base
            .iter()
            .flat_map(|&point| [point, curve.fast.neg(point)])
            .collect();
        let mut reachable = vec![HashSet::from([FastPoint::INFINITY])];
        for _ in 0..6 {
            let next = reachable.last().unwrap();
            let mut extended = next.clone();
            for &sum in next {
                for &point in &signed {
                    extended.insert(curve.fast.add(sum, point));
                }
            }
            reachable.push(extended);
        }
        for x in 0..256 {
            for target in curve.points_with_x(x) {
                for (arity, expected) in [(2, &reachable[5]), (3, &reachable[6])] {
                    let (witness, _) = oracle.query(target, arity).expect("query");
                    assert_eq!(witness.is_some(), expected.contains(&target));
                    if let Some(witness) = witness {
                        assert!(oracle.verify_coefficients(target, &witness.coefficients));
                        assert!(witness.l1_norm() <= (arity + 3) as usize);
                    }
                }
            }
        }
    }
}
