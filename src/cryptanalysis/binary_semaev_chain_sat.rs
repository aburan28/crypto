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

/// Exact additions for a signed-Frobenius x-domain represented by distinct
/// orbit representatives. Each summand chooses one representative and one
/// phase in `0..n`; sign does not alter x. The root point and every phase
/// transform are tied to the original summand bits, so this is an exact
/// finite-domain constraint rather than a post-solve filter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrobeniusOrbitDomainSize {
    pub variables: u64,
    pub clauses: usize,
    pub xor_rows: usize,
}

pub fn frobenius_orbit_domain_size(
    n: u32,
    m: usize,
    representatives: usize,
) -> Option<FrobeniusOrbitDomainSize> {
    if !(2..=127).contains(&n) || m == 0 || representatives < 2 {
        return None;
    }
    let bits = (u32::BITS - (n - 1).leading_zeros()) as usize;
    let n = n as usize;
    let vars_per = representatives
        .checked_mul(2)?
        .checked_sub(1)?
        .checked_add(n)?
        .checked_add(bits)?
        .checked_add(2usize.checked_mul(n)?.checked_mul(bits)?)?;
    let clauses_per = representatives
        .checked_mul(n)?
        .checked_add(3usize.checked_mul(representatives)?.checked_sub(3)?)?
        .checked_add((1usize << bits).checked_sub(n)?)?
        .checked_add(4usize.checked_mul(n)?.checked_mul(bits)?)?;
    Some(FrobeniusOrbitDomainSize {
        variables: u64::try_from(vars_per.checked_mul(m)?).ok()?,
        clauses: clauses_per.checked_mul(m)?,
        xor_rows: bits.checked_add(1)?.checked_mul(n)?.checked_mul(m)?,
    })
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
    first_pair_exchangeable: bool,
    first_pair_ordered: bool,
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
            first_pair_exchangeable: false,
            first_pair_ordered: false,
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
        self.first_pair_exchangeable = true;
        Ok(required)
    }

    /// Install the exact union of all binary-Frobenius phases of the given
    /// representative x-coordinates. One-hot representative selection and
    /// a bounded binary phase are linked to each summand by native XOR rows.
    /// All caps are checked before the solver is mutated.
    pub fn constrain_frobenius_orbit_summands(
        &mut self,
        representatives: &[BigUint],
        max_variables: u32,
        max_domain_clauses: usize,
    ) -> Result<FrobeniusOrbitDomainSize, ChainBuildError> {
        if self.domain_installed || representatives.iter().any(|x| x.bits() > u64::from(self.n)) {
            return Err(ChainBuildError::Unsupported);
        }
        let mut sorted = representatives.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() != representatives.len() {
            return Err(ChainBuildError::Unsupported);
        }
        let required = frobenius_orbit_domain_size(self.n, self.m, sorted.len())
            .ok_or(ChainBuildError::Unsupported)?;
        let total = u64::from(self.solver.n_vars())
            .checked_add(required.variables)
            .ok_or(ChainBuildError::Unsupported)?;
        if total > u64::from(max_variables) || total > i32::MAX as u64 {
            return Err(ChainBuildError::VariableCap {
                required: total,
                maximum: max_variables,
            });
        }
        if required.clauses > max_domain_clauses {
            return Err(ChainBuildError::DomainClauseCap {
                required: required.clauses,
                maximum: max_domain_clauses,
            });
        }
        for block in &self.blocks[..self.m] {
            add_frobenius_orbit_domain(&mut self.solver, block, &sorted, &self.irr);
        }
        self.domain_installed = true;
        self.first_pair_exchangeable = true;
        Ok(required)
    }

    /// Give each summand its own exact signed-Frobenius x-domain. This is the
    /// circuit prerequisite for small-base summands plus one or two residual
    /// domains; a caller must still lift and verify any partial point equation.
    /// Caps and every domain are checked before mutating the solver.
    pub fn constrain_frobenius_orbit_domains(
        &mut self,
        representatives_by_summand: &[Vec<BigUint>],
        max_variables: u32,
        max_domain_clauses: usize,
    ) -> Result<FrobeniusOrbitDomainSize, ChainBuildError> {
        if self.domain_installed || representatives_by_summand.len() != self.m {
            return Err(ChainBuildError::Unsupported);
        }
        let mut sorted_domains = Vec::with_capacity(self.m);
        let mut required = FrobeniusOrbitDomainSize {
            variables: 0,
            clauses: 0,
            xor_rows: 0,
        };
        for representatives in representatives_by_summand {
            if representatives.iter().any(|x| x.bits() > u64::from(self.n)) {
                return Err(ChainBuildError::Unsupported);
            }
            let mut sorted = representatives.clone();
            sorted.sort_unstable();
            sorted.dedup();
            if sorted.len() != representatives.len() {
                return Err(ChainBuildError::Unsupported);
            }
            let size = frobenius_orbit_domain_size(self.n, 1, sorted.len())
                .ok_or(ChainBuildError::Unsupported)?;
            required.variables = required
                .variables
                .checked_add(size.variables)
                .ok_or(ChainBuildError::Unsupported)?;
            required.clauses = required
                .clauses
                .checked_add(size.clauses)
                .ok_or(ChainBuildError::Unsupported)?;
            required.xor_rows = required
                .xor_rows
                .checked_add(size.xor_rows)
                .ok_or(ChainBuildError::Unsupported)?;
            sorted_domains.push(sorted);
        }
        let total = u64::from(self.solver.n_vars())
            .checked_add(required.variables)
            .ok_or(ChainBuildError::Unsupported)?;
        if total > u64::from(max_variables) || total > i32::MAX as u64 {
            return Err(ChainBuildError::VariableCap {
                required: total,
                maximum: max_variables,
            });
        }
        if required.clauses > max_domain_clauses {
            return Err(ChainBuildError::DomainClauseCap {
                required: required.clauses,
                maximum: max_domain_clauses,
            });
        }
        self.first_pair_exchangeable = sorted_domains[0] == sorted_domains[1];
        for (block, representatives) in self.blocks[..self.m].iter().zip(&sorted_domains) {
            add_frobenius_orbit_domain(&mut self.solver, block, representatives, &self.irr);
        }
        self.domain_installed = true;
        Ok(required)
    }

    /// Keep one representative of the first summand-pair swap. The first
    /// S3 node is symmetric in these two x-coordinates, and every later
    /// node depends only on its output. This is optional so ordered and
    /// unrestricted searches remain separate sweep arms.
    pub fn order_first_pair(&mut self, max_variables: u32) -> Result<(), ChainBuildError> {
        if self.first_pair_ordered || !self.first_pair_exchangeable {
            return Err(ChainBuildError::Unsupported);
        }
        let extra = self.n - 1;
        let required = u64::from(self.solver.n_vars()) + u64::from(extra);
        if required > u64::from(max_variables) || required > i32::MAX as u64 {
            return Err(ChainBuildError::VariableCap {
                required,
                maximum: max_variables,
            });
        }
        add_lex_leq(&mut self.solver, &self.blocks[0], &self.blocks[1]);
        self.first_pair_ordered = true;
        Ok(())
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

fn add_frobenius_orbit_domain(
    solver: &mut Solver,
    block: &[u32],
    representatives: &[BigUint],
    irr: &IrreduciblePoly,
) {
    fn take(next: &mut u32, count: usize) -> Vec<u32> {
        let start = *next;
        *next += count as u32;
        (start..*next).collect()
    }

    let n = block.len();
    let r = representatives.len();
    let bits = (u32::BITS - ((n as u32) - 1).leading_zeros()) as usize;
    let size = frobenius_orbit_domain_size(n as u32, 1, r).expect("preflighted domain");
    let allocated = solver.add_vars(size.variables as u32);
    let mut next = *allocated.start();
    let selectors = take(&mut next, r);
    let prefix = take(&mut next, r - 1);
    let root = take(&mut next, n);
    let phase = take(&mut next, bits);

    // Sequential exactly-one selection. No quadratic pairwise clauses.
    solver.add_clause(selectors.iter().map(|&v| v as Lit).collect());
    solver.add_clause(vec![-(selectors[0] as Lit), prefix[0] as Lit]);
    for i in 1..r - 1 {
        solver.add_clause(vec![-(selectors[i] as Lit), prefix[i] as Lit]);
        solver.add_clause(vec![-(prefix[i - 1] as Lit), prefix[i] as Lit]);
        solver.add_clause(vec![-(selectors[i] as Lit), -(prefix[i - 1] as Lit)]);
    }
    solver.add_clause(vec![-(selectors[r - 1] as Lit), -(prefix[r - 2] as Lit)]);
    for (&selector, rep) in selectors.iter().zip(representatives) {
        for (bit, &wire) in root.iter().enumerate() {
            solver.add_clause(vec![
                -(selector as Lit),
                if rep.bit(bit as u64) {
                    wire as Lit
                } else {
                    -(wire as Lit)
                },
            ]);
        }
    }

    // Bits encode an exponent in 0..n-1; every other 7-bit value at N83
    // is explicitly forbidden. Compose powers 1,2,4,... of Frobenius.
    for value in n..(1usize << bits) {
        solver.add_clause(
            phase
                .iter()
                .enumerate()
                .map(|(bit, &wire)| {
                    if (value >> bit) & 1 == 1 {
                        -(wire as Lit)
                    } else {
                        wire as Lit
                    }
                })
                .collect(),
        );
    }
    let mut current = root;
    for (power, &choose) in phase.iter().enumerate() {
        let shifted = take(&mut next, n);
        let output = take(&mut next, n);
        let mut rows = vec![Vec::new(); n];
        for (input, &wire) in current.iter().enumerate() {
            let basis = F2mElement::from_bit_positions(&[input as u32], n as u32);
            let image = basis.square_k_times(1u32 << power, irr).to_biguint();
            for (bit, row) in rows.iter_mut().enumerate() {
                if image.bit(bit as u64) {
                    row.push(wire);
                }
            }
        }
        for bit in 0..n {
            rows[bit].push(shifted[bit]);
            solver.add_xor(&rows[bit], false);
            let (input, transformed, result) =
                (current[bit] as Lit, shifted[bit] as Lit, output[bit] as Lit);
            let choose = choose as Lit;
            solver.add_clause(vec![choose, -input, result]);
            solver.add_clause(vec![choose, input, -result]);
            solver.add_clause(vec![-choose, -transformed, result]);
            solver.add_clause(vec![-choose, transformed, -result]);
        }
        current = output;
    }
    for (&output, &summand) in current.iter().zip(block) {
        solver.add_xor(&[output, summand], false);
    }
    debug_assert_eq!(next - 1, *allocated.end());
}

/// Enforce unsigned `a <= b`, with the highest field bit most significant.
/// Equality forces each next prefix wire true; no other assignment is lost.
fn add_lex_leq(solver: &mut Solver, a: &[u32], b: &[u32]) {
    debug_assert_eq!(a.len(), b.len());
    let top = a.len() - 1;
    let aux: Vec<Lit> = solver
        .add_vars(top as u32)
        .map(|wire| wire as Lit)
        .collect();
    solver.add_clause(vec![-(a[top] as Lit), b[top] as Lit]);
    for bit in (0..top).rev() {
        let above = if bit + 1 < top {
            Some(aux[bit + 1])
        } else {
            None
        };
        for (left, right) in [
            (a[bit + 1] as Lit, b[bit + 1] as Lit),
            (-(a[bit + 1] as Lit), -(b[bit + 1] as Lit)),
        ] {
            let mut clause = vec![left, right, aux[bit]];
            if let Some(prefix) = above {
                clause.push(-prefix);
            }
            solver.add_clause(clause);
        }
        solver.add_clause(vec![-aux[bit], -(a[bit] as Lit), b[bit] as Lit]);
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
    fn first_pair_lex_constraint_keeps_exactly_one_swap_orientation() {
        for width in 1..=4usize {
            let mut solver = Solver::new((2 * width) as u32);
            let a: Vec<u32> = (1..=width as u32).collect();
            let b: Vec<u32> = ((width + 1) as u32..=(2 * width) as u32).collect();
            add_lex_leq(&mut solver, &a, &b);
            let mut seen = std::collections::HashSet::new();
            while solver.solve() == SolveResult::Sat {
                let model = solver.model();
                let code = |first: usize| {
                    (0..width).fold(0u32, |word, bit| {
                        word | (u32::from(model[first + bit]) << bit)
                    })
                };
                let pair = (code(0), code(width));
                assert!(pair.0 <= pair.1);
                assert!(seen.insert(pair));
                solver.reset_search();
                solver.add_clause(
                    (0..2 * width)
                        .map(|bit| {
                            if model[bit] {
                                -((bit + 1) as Lit)
                            } else {
                                (bit + 1) as Lit
                            }
                        })
                        .collect(),
                );
            }
            let choices = 1usize << width;
            assert_eq!(seen.len(), choices * (choices + 1) / 2);
        }
    }

    #[test]
    fn first_pair_order_cap_is_checked_before_mutation() {
        let irr = IrreduciblePoly {
            degree: 7,
            low_terms: vec![0, 1],
        };
        let mut encoding =
            ChainedS3Encoding::build(5, &irr, &F2mElement::one(7), None, 0, 2_000).unwrap();
        encoding
            .constrain_summands(&[BigUint::from(1u8), BigUint::from(2u8)], 1_000)
            .unwrap();
        let variables = encoding.solver.n_vars();
        let clauses = encoding.solver.n_clauses();
        assert!(matches!(
            encoding.order_first_pair(variables + 5),
            Err(ChainBuildError::VariableCap { required, .. }) if required == u64::from(variables + 6)
        ));
        assert_eq!(encoding.solver.n_vars(), variables);
        assert_eq!(encoding.solver.n_clauses(), clauses);
        encoding.order_first_pair(variables + 6).unwrap();
        assert_eq!(encoding.solver.n_vars(), variables + 6);
        // Root assignments from the finite domain may remove clauses.
        assert!((clauses + 1..=clauses + 19).contains(&encoding.solver.n_clauses()));
        assert!(matches!(
            encoding.order_first_pair(variables + 12),
            Err(ChainBuildError::Unsupported)
        ));
    }

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

    #[test]
    fn frobenius_orbit_domain_accepts_exactly_the_small_field_closure() {
        let irr = IrreduciblePoly {
            degree: 3,
            low_terms: vec![0, 1],
        };
        let representatives = [BigUint::from(1u8), BigUint::from(2u8)];
        let mut closure = std::collections::BTreeSet::new();
        for representative in &representatives {
            let mut value = F2mElement::from_biguint(representative, 3);
            for _ in 0..3 {
                closure.insert(value.to_biguint());
                value = value.square(&irr);
            }
        }
        assert_eq!(closure.len(), 4);
        for x in 0..8u8 {
            let mut solver = Solver::new(3);
            let before = (solver.n_vars(), solver.n_clauses(), solver.n_xors());
            add_frobenius_orbit_domain(&mut solver, &[1, 2, 3], &representatives, &irr);
            let size = frobenius_orbit_domain_size(3, 1, 2).unwrap();
            assert_eq!(solver.n_vars() as u64 - before.0 as u64, size.variables);
            assert_eq!(solver.n_clauses() - before.1, size.clauses);
            assert_eq!(solver.n_xors() - before.2, size.xor_rows);
            for bit in 0..3 {
                solver.add_clause(vec![if x & (1 << bit) != 0 {
                    (bit + 1) as Lit
                } else {
                    -((bit + 1) as Lit)
                }]);
            }
            assert_eq!(
                solver.solve() == SolveResult::Sat,
                closure.contains(&BigUint::from(x)),
                "x={x}"
            );
        }
    }

    #[test]
    fn frobenius_orbit_preflight_caps_before_mutation_and_counts_n83() {
        assert_eq!(
            frobenius_orbit_domain_size(83, 5, 2_001),
            Some(FrobeniusOrbitDomainSize {
                variables: 26_265,
                clauses: 872_260,
                xor_rows: 3_320,
            })
        );
        let kc = KoblitzCurve::new(0, 7).unwrap();
        let mut enc =
            ChainedS3Encoding::build(5, &kc.curve.irreducible, &kc.curve.b, None, 0, 2_000)
                .unwrap();
        let reps = [BigUint::from(1u8), BigUint::from(2u8)];
        let size = frobenius_orbit_domain_size(7, 5, reps.len()).unwrap();
        let before = (
            enc.solver.n_vars(),
            enc.solver.n_clauses(),
            enc.solver.n_xors(),
        );
        assert_eq!(
            enc.constrain_frobenius_orbit_summands(&reps, before.0, size.clauses),
            Err(ChainBuildError::VariableCap {
                required: before.0 as u64 + size.variables,
                maximum: before.0,
            })
        );
        assert_eq!(
            enc.constrain_frobenius_orbit_summands(
                &reps,
                before.0 + size.variables as u32,
                size.clauses - 1,
            ),
            Err(ChainBuildError::DomainClauseCap {
                required: size.clauses,
                maximum: size.clauses - 1,
            })
        );
        assert_eq!(
            (
                enc.solver.n_vars(),
                enc.solver.n_clauses(),
                enc.solver.n_xors(),
            ),
            before
        );
        assert_eq!(
            enc.constrain_frobenius_orbit_summands(
                &reps,
                before.0 + size.variables as u32,
                size.clauses,
            ),
            Ok(size)
        );
        assert_eq!(enc.solver.n_vars() as u64 - before.0 as u64, size.variables);
        assert_eq!(enc.solver.n_clauses() - before.1, size.clauses);
        assert_eq!(enc.solver.n_xors() - before.2, size.xor_rows);
    }

    #[test]
    fn heterogeneous_orbit_domains_preflight_and_guard_pair_symmetry() {
        let kc = KoblitzCurve::new(0, 7).unwrap();
        let BinaryPoint::Affine { x, .. } = kc.generator() else {
            panic!("generator is affine");
        };
        let point_x = x.to_biguint();
        let alternatives: Vec<_> = (0u8..128)
            .map(BigUint::from)
            .filter(|candidate| candidate != &point_x)
            .take(3)
            .collect();
        let domains = vec![
            vec![point_x.clone(), alternatives[0].clone()],
            vec![
                point_x.clone(),
                alternatives[1].clone(),
                alternatives[2].clone(),
            ],
        ];
        let small = frobenius_orbit_domain_size(7, 1, 2).unwrap();
        let large = frobenius_orbit_domain_size(7, 1, 3).unwrap();
        let size = FrobeniusOrbitDomainSize {
            variables: small.variables + large.variables,
            clauses: small.clauses + large.clauses,
            xor_rows: small.xor_rows + large.xor_rows,
        };
        let mut enc =
            ChainedS3Encoding::build(2, &kc.curve.irreducible, &kc.curve.b, None, 0, 2_000)
                .unwrap();
        let before = (
            enc.solver.n_vars(),
            enc.solver.n_clauses(),
            enc.solver.n_xors(),
        );
        assert_eq!(
            enc.constrain_frobenius_orbit_domains(&domains, before.0, size.clauses),
            Err(ChainBuildError::VariableCap {
                required: u64::from(before.0) + size.variables,
                maximum: before.0,
            })
        );
        assert_eq!(
            enc.constrain_frobenius_orbit_domains(
                &domains,
                before.0 + size.variables as u32,
                size.clauses - 1,
            ),
            Err(ChainBuildError::DomainClauseCap {
                required: size.clauses,
                maximum: size.clauses - 1,
            })
        );
        assert_eq!(
            (
                enc.solver.n_vars(),
                enc.solver.n_clauses(),
                enc.solver.n_xors()
            ),
            before
        );
        assert_eq!(
            enc.constrain_frobenius_orbit_domains(
                &domains,
                before.0 + size.variables as u32,
                size.clauses,
            ),
            Ok(size)
        );
        assert_eq!(
            enc.solver.n_vars() as u64 - u64::from(before.0),
            size.variables
        );
        assert_eq!(enc.solver.n_clauses() - before.1, size.clauses);
        assert_eq!(enc.solver.n_xors() - before.2, size.xor_rows);
        assert_eq!(
            enc.order_first_pair(2_000),
            Err(ChainBuildError::Unsupported)
        );
        assert_eq!(enc.solver.solve(), SolveResult::Sat);
        assert!(enc.verify_model());
        let found = enc.decode_summands();
        assert_eq!(found[0], found[1]);

        let mut equal =
            ChainedS3Encoding::build(2, &kc.curve.irreducible, &kc.curve.b, None, 0, 2_000)
                .unwrap();
        let equal_domains = vec![
            domains[0].clone(),
            domains[0].iter().rev().cloned().collect(),
        ];
        let equal_size = frobenius_orbit_domain_size(7, 2, 2).unwrap();
        equal
            .constrain_frobenius_orbit_domains(&equal_domains, 2_000, equal_size.clauses)
            .unwrap();
        equal.order_first_pair(2_000).unwrap();
        assert_eq!(equal.solver.solve(), SolveResult::Sat);
        assert!(equal.verify_model());
    }
}
