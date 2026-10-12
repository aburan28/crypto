//! Deterministic polynomial-basis isomorphisms between two representations
//! of the same binary field. The n37 cold-control and archived degree-73
//! descent use distinct irreducible polynomials; their coordinate words
//! cannot be compared until this map is applied to both coordinates.

use crate::binary_ecc::{F2mElement, F2mPoly, IrreduciblePoly};
use crate::cryptanalysis::binary_isogeny::find_roots_in_f2m;
use crate::cryptanalysis::koblitz_fast::{FastCurve, FastPoint};
use crate::cryptanalysis::semaev_decomp::Gf2;

/// An isomorphism from a source polynomial basis to a target polynomial
/// basis, with a checked inverse. The image of the source generator fixes
/// one of the `n` Frobenius-conjugate isomorphisms and must be pinned when
/// experiments on separate hosts need identical coordinate words.
#[derive(Clone, Debug)]
pub struct BinaryFieldBasis {
    degree: u32,
    mask: u64,
    generator_image: u64,
    forward_columns: Vec<u64>,
    inverse_columns: Vec<u64>,
}

impl BinaryFieldBasis {
    /// Find every target-field root of the source modulus, choose the least
    /// coordinate word, and verify that its powers form a full-rank basis.
    /// This is target-independent setup; charge it to setup in a comparison.
    pub fn discover(source: &IrreduciblePoly, target: &IrreduciblePoly) -> Result<Self, String> {
        let n = checked_degree(source, target)?;
        if source.low_terms == target.low_terms {
            // In GF(2), x reduces to the constant term of the modulus.
            let image = if n == 1 {
                u64::from(source.low_terms.contains(&0))
            } else {
                2
            };
            return Self::from_generator_image(source, target, image);
        }
        let mut coeffs = vec![F2mElement::zero(n); n as usize + 1];
        coeffs[n as usize] = F2mElement::one(n);
        for &bit in &source.low_terms {
            if bit >= n {
                return Err("source modulus has an out-of-range term".into());
            }
            coeffs[bit as usize] = F2mElement::one(n);
        }
        let polynomial = F2mPoly::from_coeffs(coeffs, n);
        let mut roots: Vec<u64> = find_roots_in_f2m(&polynomial, n, target)
            .iter()
            .map(|root| {
                root.to_biguint()
                    .to_u64_digits()
                    .first()
                    .copied()
                    .unwrap_or(0)
            })
            .collect();
        roots.sort_unstable();
        roots.dedup();
        if roots.len() != n as usize {
            return Err(format!(
                "source modulus has {} target-field roots; expected {n}",
                roots.len()
            ));
        }
        Self::from_generator_image(source, target, roots[0])
    }

    /// Construct from a frozen image of the source polynomial generator.
    /// Rejects a non-root, a dependent power basis, or inconsistent field
    /// dimensions instead of silently reinterpreting coordinates.
    pub fn from_generator_image(
        source: &IrreduciblePoly,
        target: &IrreduciblePoly,
        generator_image: u64,
    ) -> Result<Self, String> {
        let n = checked_degree(source, target)?;
        let mask = (1u64 << n) - 1;
        if generator_image & !mask != 0 {
            return Err("generator image is outside target field".into());
        }
        let field = Gf2::new(target);
        let mut powers = Vec::with_capacity(n as usize + 1);
        let mut power = 1u64;
        for _ in 0..=n {
            powers.push(power);
            power = field.mul(power, generator_image);
        }
        let mut residual = powers[n as usize];
        for &bit in &source.low_terms {
            if bit >= n {
                return Err("source modulus has an out-of-range term".into());
            }
            residual ^= powers[bit as usize];
        }
        if residual != 0 {
            return Err("generator image is not a root of source modulus".into());
        }
        let forward_columns = powers[..n as usize].to_vec();
        let inverse_rows =
            invert_f2(&forward_columns, n).ok_or("generator image has a dependent power basis")?;
        let inverse_columns = (0..n)
            .map(|bit| {
                (0..n).fold(0u64, |word, row| {
                    word | (((inverse_rows[row as usize] >> bit) & 1) << row)
                })
            })
            .collect();
        Ok(Self {
            degree: n,
            mask,
            generator_image,
            forward_columns,
            inverse_columns,
        })
    }

    pub fn degree(&self) -> u32 {
        self.degree
    }

    pub fn generator_image(&self) -> u64 {
        self.generator_image
    }

    /// Source polynomial-basis word to target polynomial-basis word.
    pub fn map_element(&self, value: u64) -> Option<u64> {
        if value & !self.mask != 0 {
            return None;
        }
        Some(apply_columns(value, &self.forward_columns))
    }

    /// Target polynomial-basis word to source polynomial-basis word.
    pub fn inverse_element(&self, value: u64) -> Option<u64> {
        if value & !self.mask != 0 {
            return None;
        }
        Some(apply_columns(value, &self.inverse_columns))
    }

    /// Apply the same field map to both coordinates; infinity stays infinity.
    pub fn map_point(&self, point: FastPoint) -> Option<FastPoint> {
        if point.infinity {
            return Some(FastPoint::INFINITY);
        }
        Some(FastPoint::affine(
            self.map_element(point.x)?,
            self.map_element(point.y)?,
        ))
    }

    pub fn inverse_point(&self, point: FastPoint) -> Option<FastPoint> {
        if point.infinity {
            return Some(FastPoint::INFINITY);
        }
        Some(FastPoint::affine(
            self.inverse_element(point.x)?,
            self.inverse_element(point.y)?,
        ))
    }
}

fn checked_degree(source: &IrreduciblePoly, target: &IrreduciblePoly) -> Result<u32, String> {
    let n = source.degree;
    if n == 0 || n > FastCurve::MAX_DEGREE || target.degree != n {
        return Err("field degrees differ or exceed FastPoint capacity".into());
    }
    for modulus in [source, target] {
        let mut seen = 0u64;
        for &bit in &modulus.low_terms {
            if bit >= n || (seen >> bit) & 1 == 1 {
                return Err("field modulus has invalid or duplicate low terms".into());
            }
            seen |= 1u64 << bit;
        }
    }
    Ok(n)
}

fn apply_columns(mut value: u64, columns: &[u64]) -> u64 {
    let mut result = 0;
    while value != 0 {
        let bit = value.trailing_zeros() as usize;
        result ^= columns[bit];
        value &= value - 1;
    }
    result
}

/// Return the rows of the inverse of a binary matrix stored as columns.
/// Keep this local: frozen benchmark workflows overlay their older
/// `koblitz_fast` implementation, where its equivalent helper is private.
fn invert_f2(columns: &[u64], n: u32) -> Option<Vec<u64>> {
    let mut matrix: Vec<u64> = (0..n)
        .map(|row| {
            (0..n).fold(0u64, |bits, col| {
                bits | (((columns[col as usize] >> row) & 1) << col)
            })
        })
        .collect();
    let mut inverse: Vec<u64> = (0..n).map(|row| 1u64 << row).collect();
    for col in 0..n as usize {
        let pivot = (col..n as usize).find(|&row| (matrix[row] >> col) & 1 == 1)?;
        matrix.swap(col, pivot);
        inverse.swap(col, pivot);
        for row in 0..n as usize {
            if row != col && (matrix[row] >> col) & 1 == 1 {
                matrix[row] ^= matrix[col];
                inverse[row] ^= inverse[col];
            }
        }
    }
    Some(inverse)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::binary_velu::Curve;

    #[test]
    fn both_directions_preserve_small_field_arithmetic() {
        let source = IrreduciblePoly::deg_8();
        let target = IrreduciblePoly {
            degree: 8,
            low_terms: vec![0, 2, 3, 4],
        };
        let bridge = BinaryFieldBasis::discover(&source, &target).expect("isomorphism");
        let a = Gf2::new(&source);
        let b = Gf2::new(&target);
        assert_eq!(bridge.map_element(0), Some(0));
        assert_eq!(bridge.map_element(1), Some(1));
        assert_eq!(bridge.map_element(2), Some(bridge.generator_image()));
        assert_eq!(bridge.map_element(256), None);
        assert!(BinaryFieldBasis::from_generator_image(&source, &target, 0).is_err());
        for x in 0..256u64 {
            assert_eq!(
                bridge.inverse_element(bridge.map_element(x).unwrap()),
                Some(x)
            );
            assert_eq!(
                bridge.map_element(a.sqr(x)),
                Some(b.sqr(bridge.map_element(x).unwrap()))
            );
        }
        for x in (0..256u64).step_by(7) {
            for y in (0..256u64).step_by(11) {
                assert_eq!(
                    bridge.map_element(a.mul(x, y)),
                    Some(b.mul(
                        bridge.map_element(x).unwrap(),
                        bridge.map_element(y).unwrap()
                    ))
                );
            }
        }
        let identity = BinaryFieldBasis::discover(&source, &source).expect("identity");
        for x in 0..256u64 {
            assert_eq!(identity.map_element(x), Some(x));
        }
        let one_bit = IrreduciblePoly {
            degree: 1,
            low_terms: vec![0],
        };
        let one_bit_identity = BinaryFieldBasis::discover(&one_bit, &one_bit).expect("GF(2)");
        assert_eq!(one_bit_identity.map_element(0), Some(0));
        assert_eq!(one_bit_identity.map_element(1), Some(1));
        let one_bit_x_modulus = IrreduciblePoly {
            degree: 1,
            low_terms: vec![],
        };
        let one_bit_x_identity =
            BinaryFieldBasis::discover(&one_bit_x_modulus, &one_bit_x_modulus).expect("GF(2)");
        assert_eq!(one_bit_x_identity.generator_image(), 0);
        assert_eq!(one_bit_x_identity.map_element(1), Some(1));
    }

    #[test]
    fn frozen_n37_public_targets_cross_without_logs() {
        let frozen: serde_json::Value = serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/research/notes/ecc2k130/disjoint_cold_v2_20261001/FROZEN.json"
        )))
        .expect("frozen campaign");
        let spec = &frozen["specs"]["n37_L1024"];
        assert_eq!(spec["K"], 42);
        assert_eq!(spec["L"], 1024);
        let low_terms: Vec<u32> = serde_json::from_value(spec["field_modulus_low_terms"].clone())
            .expect("control modulus");
        assert_eq!(low_terms, vec![0, 1, 4, 6]);
        let source_modulus = IrreduciblePoly {
            degree: 37,
            low_terms,
        };
        let archive_modulus = IrreduciblePoly {
            degree: 37,
            low_terms: vec![0, 1, 2, 3, 4, 5],
        };
        let bridge = BinaryFieldBasis::discover(&source_modulus, &archive_modulus)
            .expect("frozen-to-archive basis bridge");
        assert_eq!(bridge.generator_image(), 10_156_182_909);
        let replay = BinaryFieldBasis::from_generator_image(
            &source_modulus,
            &archive_modulus,
            10_156_182_909,
        )
        .expect("frozen root reconstructs the same bridge");
        let source_field = Gf2::new(&source_modulus);
        let archive_field = Gf2::new(&archive_modulus);
        let mask = (1u64 << 37) - 1;
        let mut state = 0x37_2026_1002u64;
        for _ in 0..256 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let x = state & mask;
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let y = state & mask;
            let mapped_x = bridge.map_element(x).unwrap();
            let mapped_y = bridge.map_element(y).unwrap();
            assert_eq!(bridge.inverse_element(mapped_x), Some(x));
            assert_eq!(
                bridge.map_element(source_field.sqr(x)),
                Some(archive_field.sqr(mapped_x))
            );
            assert_eq!(
                bridge.map_element(source_field.mul(x, y)),
                Some(archive_field.mul(mapped_x, mapped_y))
            );
        }
        let r = spec["subgroup_order"].as_u64().expect("subgroup order");
        assert_eq!(r, 230_603_167);
        let order = 137_439_487_532;
        let source = Curve::with_order(37, &source_modulus, 0, 1, order).expect("source model");
        let archive = Curve::with_order(37, &archive_modulus, 0, 1, order).expect("archive model");
        let generator_words: [u64; 2] =
            serde_json::from_value(spec["generator"].clone()).expect("frozen generator");
        let generator = FastPoint::affine(generator_words[0], generator_words[1]);
        assert!(source.fast.is_on_curve(generator));
        assert!(source.fast.mul_u64(generator, r).infinity);
        let mapped_generator = bridge.map_point(generator).expect("mapped generator");
        assert!(archive.fast.is_on_curve(mapped_generator));
        assert!(archive.fast.mul_u64(mapped_generator, r).infinity);
        assert_eq!(bridge.inverse_point(mapped_generator), Some(generator));

        let points = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b00.points.jsonl"
        ));
        let mut previous = None;
        let mut count = 0;
        for line in points.lines() {
            let words: [u64; 2] = serde_json::from_str(line).expect("frozen public point");
            let q = FastPoint::affine(words[0], words[1]);
            assert!(source.fast.is_on_curve(q));
            let image = bridge.map_point(q).expect("mapped public point");
            assert_eq!(replay.map_point(q), Some(image));
            assert!(archive.fast.is_on_curve(image));
            assert_eq!(bridge.inverse_point(image), Some(q));
            assert!(archive.fast.mul_u64(image, r).infinity);
            if count < 16 {
                assert_eq!(
                    bridge.map_point(source.fast.mul_u64(q, 7)),
                    Some(archive.fast.mul_u64(image, 7))
                );
                if let Some(p) = previous {
                    assert_eq!(
                        bridge.map_point(source.fast.add(p, q)),
                        Some(archive.fast.add(bridge.map_point(p).unwrap(), image))
                    );
                }
            }
            previous = Some(q);
            count += 1;
        }
        assert_eq!(count, 1024);
    }
}
