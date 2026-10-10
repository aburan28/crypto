//! Full-width, native-XOR SAT circuits for chained binary Semaev S3.
//!
//! A chain of `m-1` S3 equations is a necessary x-coordinate condition for
//! an `m`-point decomposition. Intermediate sums can be the identity, which
//! has no x-coordinate, so each possible intermediate-identity pattern is a
//! separate finite SAT instance. The caller must check every model against
//! the original equations and lift its summands in the curve group. A capped
//! pattern is inconclusive, not a refutation.

use num_bigint::BigUint;

use crate::binary_ecc::{F2mElement, IrreduciblePoly};
use crate::cryptanalysis::binary_semaev::binary_semaev_s3;
use crate::cryptanalysis::sat::{Lit, Solver};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChainBuildError {
    Unsupported,
    VariableCap { required: u64, maximum: u32 },
    DomainClauseCap { required: usize, maximum: usize },
}

/// Worst-case SAT variable count for the field-product circuit, before a
/// finite-coordinate domain. The bound is exact for the all-affine pattern.
pub fn chain_s3_max_variables(n: u32, m: usize) -> Option<u64> {
    if !(2..=6).contains(&m) || !(2..=127).contains(&n) {
        return None;
    }
    let n = u64::from(n);
    let m = m as u64;
    // Summands, m-2 internal outputs and target; one sum block and three
    // reduced products per S3. Squaring is a linear map in characteristic 2.
    Some(n * (2 * m - 1) + 4 * n * (m - 1) + 3 * n * n * (m - 1))
}

/// Exact number of missing-prefix clauses for `m` identical finite x-domains.
/// This is computed without constructing a solver or enumerating the domain's
/// complement, and is a memory preflight rather than a runtime estimate.
pub fn finite_domain_clause_count(codes: &[BigUint], n: u32, m: usize) -> Option<usize> {
    if n == 0 || m == 0 || codes.iter().any(|code| code.bits() > u64::from(n)) {
        return None;
    }
    fn count(codes: &[BigUint], bit: usize) -> Option<usize> {
        if codes.is_empty() {
            return Some(1);
        }
        if bit == 0 {
            return Some(0);
        }
        let position = bit - 1;
        let split = codes.partition_point(|code| !code.bit(position as u64));
        count(&codes[..split], position)?.checked_add(count(&codes[split..], position)?)
    }
    let mut sorted = codes.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    count(&sorted, n as usize)?.checked_mul(m)
}

fn reduced_monomial(degree: usize, irr: &IrreduciblePoly) -> Vec<usize> {
    let n = irr.degree as usize;
    let mut bits = vec![false; 2 * n - 1];
    bits[degree] = true;
    for high in (n..=degree).rev() {
        if bits[high] {
            bits[high] = false;
            for &low in &irr.low_terms {
                bits[high - n + low as usize] ^= true;
            }
        }
    }
    (0..n).filter(|&bit| bits[bit]).collect()
}

struct Circuit {
    n: usize,
    next: u32,
    reductions: Vec<Vec<usize>>,
    square_inputs: Vec<Vec<usize>>,
    ands: Vec<(u32, u32, u32)>,
    xors: Vec<(Vec<u32>, bool)>,
    units: Vec<(u32, bool)>,
    impossible: bool,
}

impl Circuit {
    fn new(irr: &IrreduciblePoly) -> Self {
        let n = irr.degree as usize;
        let reductions: Vec<_> = (0..2 * n - 1)
            .map(|degree| reduced_monomial(degree, irr))
            .collect();
        let mut square_inputs = vec![Vec::new(); n];
        for bit in 0..n {
            for &output in &reductions[2 * bit] {
                square_inputs[output].push(bit);
            }
        }
        Self {
            n,
            next: 1,
            reductions,
            square_inputs,
            ands: Vec::new(),
            xors: Vec::new(),
            units: Vec::new(),
            impossible: false,
        }
    }

    fn wire(&mut self) -> u32 {
        let wire = self.next;
        self.next += 1;
        wire
    }

    fn block(&mut self) -> Vec<u32> {
        (0..self.n).map(|_| self.wire()).collect()
    }

    fn equal(&mut self, left: &[u32], right: &[u32]) {
        for (&a, &b) in left.iter().zip(right) {
            self.xors.push((vec![a, b], false));
        }
    }

    fn product(&mut self, left: &[u32], right: &[u32]) -> Vec<u32> {
        let mut rows = vec![Vec::new(); self.n];
        for (i, &a) in left.iter().enumerate() {
            for (j, &b) in right.iter().enumerate() {
                let gate = self.wire();
                self.ands.push((gate, a, b));
                for &output in &self.reductions[i + j] {
                    rows[output].push(gate);
                }
            }
        }
        let output = self.block();
        for (bit, mut row) in rows.into_iter().enumerate() {
            row.push(output[bit]);
            self.xors.push((row, false));
        }
        output
    }

    fn s3(&mut self, x: &[u32], y: &[u32], z: &[u32], b: &F2mElement) {
        let sum = self.block();
        for bit in 0..self.n {
            self.xors.push((vec![sum[bit], x[bit], y[bit]], false));
        }
        let xy = self.product(x, y);
        let sum_z = self.product(&sum, z);
        let xyz = self.product(&xy, z);
        let b_bits = b.to_biguint();
        for output in 0..self.n {
            let mut row = vec![xyz[output]];
            for &input in &self.square_inputs[output] {
                row.extend([xy[input], sum_z[input]]);
            }
            self.xors.push((row, b_bits.bit(output as u64)));
        }
    }

    fn finish(self) -> Solver {
        let mut solver = Solver::new(self.next - 1);
        for (gate, left, right) in self.ands {
            let (gate, left, right) = (gate as Lit, left as Lit, right as Lit);
            solver.add_clause(vec![-gate, left]);
            solver.add_clause(vec![-gate, right]);
            solver.add_clause(vec![gate, -left, -right]);
        }
        for (row, rhs) in self.xors {
            solver.add_xor(&row, rhs);
        }
        for (wire, value) in self.units {
            solver.add_clause(vec![if value { wire as Lit } else { -(wire as Lit) }]);
        }
        if self.impossible {
            solver.add_clause(vec![]);
        }
        solver
    }
}

/// An x-coordinate circuit for one intermediate-identity pattern. Summand
/// blocks occupy SAT variables `1..=m*n`, followed by internal and target
/// blocks. The product gates and all S3 rows use native XOR constraints.
pub struct ChainedS3Encoding {
    pub solver: Solver,
    pub n_and_gates: usize,
    pub n_s3_nodes: usize,
    n: u32,
    m: usize,
    identity_mask: u32,
    blocks: Vec<Vec<u32>>,
    target_x: Option<F2mElement>,
    b: F2mElement,
    irr: IrreduciblePoly,
    domain_installed: bool,
}

impl ChainedS3Encoding {
    pub fn build(
        m: usize,
        irr: &IrreduciblePoly,
        b: &F2mElement,
        target_x: Option<&F2mElement>,
        identity_mask: u32,
        max_variables: u32,
    ) -> Result<Self, ChainBuildError> {
        let n = irr.degree;
        let Some(required) = chain_s3_max_variables(n, m) else {
            return Err(ChainBuildError::Unsupported);
        };
        if max_variables == 0 || required > u64::from(max_variables) {
            return Err(ChainBuildError::VariableCap {
                required,
                maximum: max_variables,
            });
        }
        if b.m_value() != n
            || target_x.is_some_and(|x| x.m_value() != n)
            || identity_mask >= (1u32 << (m - 2))
            || !irr.low_terms.contains(&0)
            || irr.low_terms.iter().any(|&term| term >= n)
            || {
                let mut terms = irr.low_terms.clone();
                terms.sort_unstable();
                terms.dedup();
                terms.len() != irr.low_terms.len()
            }
        {
            return Err(ChainBuildError::Unsupported);
        }
        let mut circuit = Circuit::new(irr);
        let blocks: Vec<Vec<u32>> = (0..2 * m - 1).map(|_| circuit.block()).collect();
        let target = &blocks[2 * m - 2];
        let target_bits = target_x.map(F2mElement::to_biguint);
        for (bit, &wire) in target.iter().enumerate() {
            circuit.units.push((
                wire,
                target_bits.as_ref().is_some_and(|x| x.bit(bit as u64)),
            ));
        }

        let mut n_s3_nodes = 0;
        let mut previous = &blocks[0];
        for step in 1..m {
            let output = if step == m - 1 {
                target
            } else {
                &blocks[m + step - 1]
            };
            let previous_infinity = step > 1 && (identity_mask & (1 << (step - 2))) != 0;
            let output_infinity = if step == m - 1 {
                target_x.is_none()
            } else {
                (identity_mask & (1 << (step - 1))) != 0
            };
            if step < m - 1 && output_infinity {
                for &wire in output {
                    circuit.units.push((wire, false));
                }
            }
            match (previous_infinity, output_infinity) {
                (true, true) => circuit.impossible = true, // O + affine != O
                (true, false) => circuit.equal(&blocks[step], output),
                (false, true) => circuit.equal(previous, &blocks[step]),
                (false, false) => {
                    circuit.s3(previous, &blocks[step], output, b);
                    n_s3_nodes += 1;
                }
            }
            previous = output;
        }
        let n_and_gates = circuit.ands.len();
        let mut solver = circuit.finish();
        solver.set_branch_priority(&(1..=(m as u32 * n)).collect::<Vec<_>>());
        Ok(Self {
            solver,
            n_and_gates,
            n_s3_nodes,
            n,
            m,
            identity_mask,
            blocks,
            target_x: target_x.cloned(),
            b: b.clone(),
            irr: irr.clone(),
            domain_installed: false,
        })
    }

    /// Add the exact finite set of rational factor-base x-coordinates to
    /// every summand. No constraint is placed on affine intermediate sums.
    pub fn constrain_summands(
        &mut self,
        codes: &[BigUint],
        max_domain_clauses: usize,
    ) -> Result<usize, ChainBuildError> {
        if self.domain_installed {
            return Err(ChainBuildError::Unsupported);
        }
        let required = finite_domain_clause_count(codes, self.n, self.m)
            .ok_or(ChainBuildError::Unsupported)?;
        if required > max_domain_clauses {
            return Err(ChainBuildError::DomainClauseCap {
                required,
                maximum: max_domain_clauses,
            });
        }
        let mut codes = codes.to_vec();
        codes.sort_unstable();
        codes.dedup();
        for block in &self.blocks[..self.m] {
            add_coordinate_domain(&mut self.solver, block, &codes, self.n as usize);
        }
        self.domain_installed = true;
        Ok(required)
    }

    fn decode_block(&self, block: &[u32], model: &[bool]) -> F2mElement {
        let set: Vec<u32> = block
            .iter()
            .enumerate()
            .filter_map(|(bit, &wire)| model[wire as usize - 1].then_some(bit as u32))
            .collect();
        F2mElement::from_bit_positions(&set, self.n)
    }

    pub fn decode_summands(&self) -> Vec<F2mElement> {
        let model = self.solver.model();
        self.blocks[..self.m]
            .iter()
            .map(|block| self.decode_block(block, &model))
            .collect()
    }

    /// Re-evaluate every active S3 over the original field and each
    /// intermediate-identity equality. Call this after a SAT result.
    pub fn verify_model(&self) -> bool {
        let model = self.solver.model();
        let values: Vec<_> = self
            .blocks
            .iter()
            .map(|block| self.decode_block(block, &model))
            .collect();
        if self.target_x.as_ref() != Some(&values[2 * self.m - 2])
            && !(self.target_x.is_none() && values[2 * self.m - 2].is_zero())
        {
            return false;
        }
        let mut previous = &values[0];
        for step in 1..self.m {
            let output = if step == self.m - 1 {
                &values[2 * self.m - 2]
            } else {
                &values[self.m + step - 1]
            };
            let previous_infinity = step > 1 && (self.identity_mask & (1 << (step - 2))) != 0;
            let output_infinity = if step == self.m - 1 {
                self.target_x.is_none()
            } else {
                (self.identity_mask & (1 << (step - 1))) != 0
            };
            let valid = match (previous_infinity, output_infinity) {
                (true, true) => false,
                (true, false) => output == &values[step],
                (false, true) => previous == &values[step],
                (false, false) => {
                    binary_semaev_s3(previous, &values[step], output, &self.b, &self.irr).is_zero()
                }
            };
            if !valid || (step < self.m - 1 && output_infinity && !output.is_zero()) {
                return false;
            }
            previous = output;
        }
        true
    }
}

fn add_coordinate_domain(solver: &mut Solver, block: &[u32], codes: &[BigUint], bit: usize) {
    fn visit(
        solver: &mut Solver,
        block: &[u32],
        codes: &[BigUint],
        bit: usize,
        prefix: &mut Vec<Lit>,
    ) {
        if codes.is_empty() {
            solver.add_clause(prefix.clone());
            return;
        }
        if bit == 0 {
            return;
        }
        let position = bit - 1;
        let split = codes.partition_point(|code| !code.bit(position as u64));
        let wire = block[position] as Lit;
        prefix.push(wire);
        visit(solver, block, &codes[..split], position, prefix);
        *prefix.last_mut().expect("pushed") = -wire;
        visit(solver, block, &codes[split..], position, prefix);
        prefix.pop();
    }
    visit(solver, block, codes, bit, &mut Vec::new());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::binary_ecc::BinaryPoint;
    use crate::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
    use crate::cryptanalysis::sat::SolveResult;

    #[test]
    fn product_gates_match_field_multiplication_and_reject_a_wrong_output() {
        for n in [7u32, 9] {
            let irr = match n {
                7 => IrreduciblePoly {
                    degree: 7,
                    low_terms: vec![0, 1],
                },
                _ => IrreduciblePoly {
                    degree: 9,
                    low_terms: vec![0, 4],
                },
            };
            for (a, b) in [(3u32, 7u32), (65, 31), (127, 255)] {
                let left = F2mElement::from_biguint(&BigUint::from(a), n);
                let right = F2mElement::from_biguint(&BigUint::from(b), n);
                let expected = left.mul(&right, &irr).to_biguint();
                let mut circuit = Circuit::new(&irr);
                let x = circuit.block();
                let y = circuit.block();
                let z = circuit.product(&x, &y);
                for (bit, &wire) in x.iter().enumerate() {
                    circuit
                        .units
                        .push((wire, left.to_biguint().bit(bit as u64)));
                }
                for (bit, &wire) in y.iter().enumerate() {
                    circuit
                        .units
                        .push((wire, right.to_biguint().bit(bit as u64)));
                }
                let mut solver = circuit.finish();
                assert_eq!(solver.solve(), SolveResult::Sat);
                let model = solver.model();
                let actual =
                    BigUint::from(z.iter().enumerate().fold(0u32, |bits, (bit, &wire)| {
                        bits | (u32::from(model[wire as usize - 1]) << bit)
                    }));
                assert_eq!(actual, expected);
                solver.reset_search();
                solver.add_clause(vec![if expected.bit(0) {
                    -(z[0] as Lit)
                } else {
                    z[0] as Lit
                }]);
                assert_eq!(solver.solve(), SolveResult::Unsat);
            }
        }
    }

    #[test]
    fn product_gates_preserve_bits_above_one_word_in_the_pinned_field() {
        let n = 83u32;
        let irr = IrreduciblePoly {
            degree: n,
            low_terms: vec![0, 1, 2, 45],
        };
        let left_bits =
            (BigUint::from(1u8) << 82) | (BigUint::from(1u8) << 64) | BigUint::from(9u8);
        let right_bits =
            (BigUint::from(1u8) << 81) | (BigUint::from(1u8) << 65) | BigUint::from(3u8);
        let left = F2mElement::from_biguint(&left_bits, n);
        let right = F2mElement::from_biguint(&right_bits, n);
        let expected = left.mul(&right, &irr).to_biguint();
        let mut circuit = Circuit::new(&irr);
        let x = circuit.block();
        let y = circuit.block();
        let z = circuit.product(&x, &y);
        for (bit, &wire) in x.iter().enumerate() {
            circuit.units.push((wire, left_bits.bit(bit as u64)));
        }
        for (bit, &wire) in y.iter().enumerate() {
            circuit.units.push((wire, right_bits.bit(bit as u64)));
        }
        let mut solver = circuit.finish();
        assert_eq!(solver.solve(), SolveResult::Sat);
        let model = solver.model();
        for (bit, &wire) in z.iter().enumerate() {
            assert_eq!(model[wire as usize - 1], expected.bit(bit as u64));
        }
    }

    #[test]
    fn s3_circuit_agrees_with_original_field_polynomial_on_roots_and_nonroots() {
        let n = 7u32;
        let irr = IrreduciblePoly {
            degree: n,
            low_terms: vec![0, 1],
        };
        let b = F2mElement::one(n);
        let mut cases = Vec::new();
        let (mut roots, mut nonroots) = (0, 0);
        for x in 0..8u32 {
            for y in 0..8u32 {
                for z in 0..128u32 {
                    let values =
                        [x, y, z].map(|word| F2mElement::from_biguint(&BigUint::from(word), n));
                    let is_root =
                        binary_semaev_s3(&values[0], &values[1], &values[2], &b, &irr).is_zero();
                    if (is_root && roots < 4) || (!is_root && nonroots < 4) {
                        cases.push((values, is_root));
                        if is_root {
                            roots += 1;
                        } else {
                            nonroots += 1;
                        }
                    }
                    if roots == 4 && nonroots == 4 {
                        break;
                    }
                }
                if roots == 4 && nonroots == 4 {
                    break;
                }
            }
            if roots == 4 && nonroots == 4 {
                break;
            }
        }
        assert_eq!((roots, nonroots), (4, 4));
        for (values, is_root) in cases {
            let mut circuit = Circuit::new(&irr);
            let blocks = [circuit.block(), circuit.block(), circuit.block()];
            circuit.s3(&blocks[0], &blocks[1], &blocks[2], &b);
            for (block, value) in blocks.iter().zip(values.iter()) {
                let bits = value.to_biguint();
                for (bit, &wire) in block.iter().enumerate() {
                    circuit.units.push((wire, bits.bit(bit as u64)));
                }
            }
            let mut solver = circuit.finish();
            assert_eq!(solver.solve() == SolveResult::Sat, is_root);
        }
    }

    #[test]
    fn source_bound_covers_degree_83_without_constructing_a_model() {
        assert_eq!(chain_s3_max_variables(83, 5), Some(84_743));
        assert_eq!(chain_s3_max_variables(83, 6), Some(105_908));
        assert_eq!(chain_s3_max_variables(83, 7), None);
    }

    #[test]
    fn intermediate_and_final_identity_patterns_preserve_real_group_sums() {
        let kc = KoblitzCurve::new(0, 7).unwrap();
        let point = kc.generator();
        let BinaryPoint::Affine { x: point_x, .. } = point else {
            panic!("generator is affine");
        };
        let three = kc.mul(point, &BigUint::from(3u32));
        let BinaryPoint::Affine { x: three_x, .. } = three else {
            panic!("three times generator is affine");
        };
        // P + (-P) + P + P + P = 3P. The first internal sum is O.
        let mut affine = ChainedS3Encoding::build(
            5,
            &kc.curve.irreducible,
            &kc.curve.b,
            Some(&three_x),
            0b001,
            2_000,
        )
        .unwrap();
        affine
            .constrain_summands(&[point_x.to_biguint()], 1_000)
            .unwrap();
        assert_eq!(affine.solver.solve(), SolveResult::Sat);
        assert!(affine.verify_model());
        assert!(affine.decode_summands().iter().all(|x| x == point_x));

        // Three P,-P pairs sum to O. Internal sums after 2 and 4 points
        // are O, while the intermediate sums after 3 and 5 are affine.
        let mut identity =
            ChainedS3Encoding::build(6, &kc.curve.irreducible, &kc.curve.b, None, 0b0101, 2_000)
                .unwrap();
        identity
            .constrain_summands(&[point_x.to_biguint()], 1_000)
            .unwrap();
        assert_eq!(identity.solver.solve(), SolveResult::Sat);
        assert!(identity.verify_model());
        assert!(identity.decode_summands().iter().all(|x| x == point_x));
    }

    #[test]
    fn finite_domain_clause_preflight_is_exact_and_rejects_before_installation() {
        let codes = [BigUint::from(0u32), BigUint::from(1u32)];
        assert_eq!(finite_domain_clause_count(&codes, 3, 5), Some(10));
        let kc = KoblitzCurve::new(0, 7).unwrap();
        let mut enc =
            ChainedS3Encoding::build(5, &kc.curve.irreducible, &kc.curve.b, None, 0, 2_000)
                .unwrap();
        let before = enc.solver.n_clauses();
        let required = finite_domain_clause_count(&codes, 7, 5).unwrap();
        assert_eq!(
            enc.constrain_summands(&codes, required - 1),
            Err(ChainBuildError::DomainClauseCap {
                required,
                maximum: required - 1,
            })
        );
        assert_eq!(enc.solver.n_clauses(), before);
        assert_eq!(enc.constrain_summands(&codes, required), Ok(required));
        assert_eq!(enc.solver.n_clauses(), before + required);
    }
}
