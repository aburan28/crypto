//! Factored n83 S3 circuit and independently checked SAT model.
//! Protocol: research/f6_n83_factored_s3_sat_20261005/PROTOCOL.md.
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::time::Instant;

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, KoblitzCurve,
};
use crypto_lib::cryptanalysis::wide_groebner::{TwoWordFieldStructure, WideFieldTable};
use crypto_lib::cryptanalysis::wide_sixsum::{Mono512, System512};
use num_bigint::BigUint;
use serde_json::json;

const N: usize = 83;
const INPUTS: usize = 339;
const PLANTED: [usize; 5] = [0, 2, 4, 6, 8];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Bit {
    Const(bool),
    Lit(i32),
}

impl Bit {
    fn not(self) -> Self {
        match self {
            Self::Const(value) => Self::Const(!value),
            Self::Lit(lit) => Self::Lit(-lit),
        }
    }
}

#[derive(Clone, Copy)]
enum Kind {
    Xor,
    And,
}

#[derive(Clone, Copy)]
struct Gate {
    kind: Kind,
    a: i32,
    b: i32,
    out: i32,
}

struct Circuit {
    next_var: i32,
    gates: Vec<Gate>,
    units: Vec<Bit>,
    outputs: Vec<Bit>,
    xors: usize,
    ands: usize,
}

impl Circuit {
    fn new() -> Self {
        Self {
            next_var: INPUTS as i32 + 1,
            gates: Vec::new(),
            units: Vec::new(),
            outputs: Vec::new(),
            xors: 0,
            ands: 0,
        }
    }

    fn input(index: usize) -> Bit {
        assert!(index < INPUTS);
        Bit::Lit(index as i32 + 1)
    }

    fn gate(&mut self, kind: Kind, a: i32, b: i32) -> Bit {
        let out = self.next_var;
        self.next_var = self.next_var.checked_add(1).expect("SAT variable limit");
        self.gates.push(Gate { kind, a, b, out });
        match kind {
            Kind::Xor => self.xors += 1,
            Kind::And => self.ands += 1,
        }
        Bit::Lit(out)
    }

    fn xor(&mut self, a: Bit, b: Bit) -> Bit {
        match (a, b) {
            (Bit::Const(false), other) | (other, Bit::Const(false)) => other,
            (Bit::Const(true), other) | (other, Bit::Const(true)) => other.not(),
            (Bit::Lit(x), Bit::Lit(y)) if x == y => Bit::Const(false),
            (Bit::Lit(x), Bit::Lit(y)) if x == -y => Bit::Const(true),
            (Bit::Lit(x), Bit::Lit(y)) => self.gate(Kind::Xor, x, y),
        }
    }

    fn and(&mut self, a: Bit, b: Bit) -> Bit {
        match (a, b) {
            (Bit::Const(false), _) | (_, Bit::Const(false)) => Bit::Const(false),
            (Bit::Const(true), other) | (other, Bit::Const(true)) => other,
            (Bit::Lit(x), Bit::Lit(y)) if x == y => Bit::Lit(x),
            (Bit::Lit(x), Bit::Lit(y)) if x == -y => Bit::Const(false),
            (Bit::Lit(x), Bit::Lit(y)) => self.gate(Kind::And, x, y),
        }
    }

    fn xor_many(&mut self, mut bits: Vec<Bit>) -> Bit {
        if bits.is_empty() {
            return Bit::Const(false);
        }
        while bits.len() > 1 {
            let mut next = Vec::with_capacity(bits.len().div_ceil(2));
            for pair in bits.chunks(2) {
                next.push(if pair.len() == 2 {
                    self.xor(pair[0], pair[1])
                } else {
                    pair[0]
                });
            }
            bits = next;
        }
        bits[0]
    }

    fn field_add(&mut self, a: &[Bit], b: &[Bit]) -> Vec<Bit> {
        (0..N).map(|i| self.xor(a[i], b[i])).collect()
    }

    fn field_square(&mut self, a: &[Bit], table: &impl WideFieldTable) -> Vec<Bit> {
        let mut terms = vec![Vec::new(); N];
        for (i, &bit) in a.iter().enumerate() {
            if bit == Bit::Const(false) {
                continue;
            }
            let mut mask = table.square_bits(i);
            while mask != 0 {
                let output = mask.trailing_zeros() as usize;
                mask &= mask - 1;
                terms[output].push(bit);
            }
        }
        terms.into_iter().map(|bits| self.xor_many(bits)).collect()
    }

    fn field_mul(&mut self, a: &[Bit], b: &[Bit], table: &impl WideFieldTable) -> Vec<Bit> {
        let mut terms = vec![Vec::new(); N];
        for (i, &left) in a.iter().enumerate() {
            if left == Bit::Const(false) {
                continue;
            }
            for (j, &right) in b.iter().enumerate() {
                if right == Bit::Const(false) {
                    continue;
                }
                let product = self.and(left, right);
                if product == Bit::Const(false) {
                    continue;
                }
                let mut mask = table.product_bits(i, j);
                while mask != 0 {
                    let output = mask.trailing_zeros() as usize;
                    mask &= mask - 1;
                    terms[output].push(product);
                }
            }
        }
        terms.into_iter().map(|bits| self.xor_many(bits)).collect()
    }

    fn s3(&mut self, x: &[Bit], y: &[Bit], z: &[Bit], b: &[Bit], table: &impl WideFieldTable) {
        let sum = self.field_add(x, y);
        let sum_sq = self.field_square(&sum, table);
        let z_sq = self.field_square(z, table);
        let first = self.field_mul(&sum_sq, &z_sq, table);
        let xy = self.field_mul(x, y, table);
        let second = self.field_mul(&xy, z, table);
        let third = self.field_square(&xy, table);
        for i in 0..N {
            let a = self.xor(first[i], second[i]);
            let a = self.xor(a, third[i]);
            let output = self.xor(a, b[i]);
            self.outputs.push(output);
        }
    }

    fn value(bit: Bit, values: &[bool]) -> bool {
        match bit {
            Bit::Const(value) => value,
            Bit::Lit(lit) if lit > 0 => values[lit as usize],
            Bit::Lit(lit) => !values[(-lit) as usize],
        }
    }

    fn evaluate(&self, inputs: &Mono512) -> Vec<bool> {
        let mut values = vec![false; self.next_var as usize];
        for i in 0..INPUTS {
            values[i + 1] = inputs.0[i / 64] >> (i % 64) & 1 != 0;
        }
        for gate in &self.gates {
            let a = Self::value(Bit::Lit(gate.a), &values);
            let b = Self::value(Bit::Lit(gate.b), &values);
            values[gate.out as usize] = match gate.kind {
                Kind::Xor => a ^ b,
                Kind::And => a & b,
            };
        }
        self.outputs
            .iter()
            .map(|&bit| Self::value(bit, &values))
            .collect()
    }

    fn write_cnf(&self, path: &str) -> std::io::Result<(usize, usize, u64)> {
        let unit = self
            .units
            .iter()
            .filter(|&&constraint| constraint != Bit::Const(true))
            .count();
        let clauses = self.xors * 4 + self.ands * 3 + unit;
        let unit_literals = self
            .units
            .iter()
            .filter(|&&constraint| matches!(constraint, Bit::Lit(_)))
            .count();
        let literals = (self.xors * 12 + self.ands * 7 + unit_literals) as u64;
        let mut out = BufWriter::new(File::create(path)?);
        writeln!(out, "p cnf {} {clauses}", self.next_var - 1)?;
        for gate in &self.gates {
            let (a, b, o) = (gate.a, gate.b, gate.out);
            match gate.kind {
                Kind::Xor => {
                    writeln!(out, "{} {} {} 0", -a, -b, -o)?;
                    writeln!(out, "{a} {b} {} 0", -o)?;
                    writeln!(out, "{a} {} {o} 0", -b)?;
                    writeln!(out, "{} {b} {o} 0", -a)?;
                }
                Kind::And => {
                    writeln!(out, "{} {} {o} 0", -a, -b)?;
                    writeln!(out, "{a} {} 0", -o)?;
                    writeln!(out, "{b} {} 0", -o)?;
                }
            }
        }
        for &bit in &self.units {
            match bit {
                Bit::Const(true) => {}
                Bit::Const(false) => writeln!(out, "0")?,
                Bit::Lit(lit) => writeln!(out, "{lit} 0")?,
            }
        }
        out.flush()?;
        Ok((clauses, unit, literals))
    }

    fn write_xcnf(&self, path: &str) -> std::io::Result<(usize, usize, u64)> {
        let unit = self
            .units
            .iter()
            .filter(|&&constraint| constraint != Bit::Const(true))
            .count();
        let unit_literals = self
            .units
            .iter()
            .filter(|&&constraint| matches!(constraint, Bit::Lit(_)))
            .count();
        let constraints = self.xors + self.ands * 3 + unit;
        let literals = (self.xors * 3 + self.ands * 7 + unit_literals) as u64;
        let mut out = BufWriter::new(File::create(path)?);
        writeln!(out, "p cnf {} {constraints}", self.next_var - 1)?;
        for gate in &self.gates {
            let (a, b, o) = (gate.a, gate.b, gate.out);
            match gate.kind {
                Kind::Xor => {
                    // Extended DIMACS uses odd parity for all-positive
                    // literals. Negating the first literal makes it even.
                    let parity = (a < 0) ^ (b < 0);
                    let first = if parity { a.abs() } else { -a.abs() };
                    writeln!(out, "x {first} {} {o} 0", b.abs())?;
                }
                Kind::And => {
                    writeln!(out, "{} {} {o} 0", -a, -b)?;
                    writeln!(out, "{a} {} 0", -o)?;
                    writeln!(out, "{b} {} 0", -o)?;
                }
            }
        }
        for &bit in &self.units {
            match bit {
                Bit::Const(true) => {}
                Bit::Const(false) => writeln!(out, "0")?,
                Bit::Lit(lit) => writeln!(out, "{lit} 0")?,
            }
        }
        out.flush()?;
        Ok((constraints, unit, literals))
    }
}

fn bits_of(element: &F2mElement) -> Vec<Bit> {
    let words = element.raw_bits();
    (0..N)
        .map(|i| Bit::Const(words[i / 64] >> (i % 64) & 1 != 0))
        .collect()
}

fn x_of(point: &BinaryPoint) -> &F2mElement {
    match point {
        BinaryPoint::Affine { x, .. } => x,
        BinaryPoint::Infinity => panic!("affine point required"),
    }
}

fn input_assignment(points: &[BinaryPoint], indices: &[usize; 5], kc: &KoblitzCurve) -> Mono512 {
    let mut assignment = Mono512::default();
    let mut prefix = BinaryPoint::Infinity;
    for (summand, &index) in indices.iter().enumerate() {
        let point = &points[index];
        let words = x_of(point).raw_bits();
        assert_eq!(words[0] >> 18, 0);
        assert!(words.get(1).is_none_or(|&word| word == 0));
        for bit in 0..18 {
            if words[0] >> bit & 1 != 0 {
                let variable = summand * 18 + bit;
                assignment.0[variable / 64] |= 1u64 << (variable % 64);
            }
        }
        prefix = kc.add(&prefix, point);
        if (1..=3).contains(&summand) {
            let words = x_of(&prefix).raw_bits();
            let start = 90 + (summand - 1) * N;
            for bit in 0..N {
                if words[bit / 64] >> (bit % 64) & 1 != 0 {
                    let variable = start + bit;
                    assignment.0[variable / 64] |= 1u64 << (variable % 64);
                }
            }
        }
    }
    assignment
}

fn point_json(point: &BinaryPoint) -> serde_json::Value {
    match point {
        BinaryPoint::Infinity => json!({"infinity":true}),
        BinaryPoint::Affine { x, y } => json!({
            "x":x.to_biguint().to_str_radix(16), "y":y.to_biguint().to_str_radix(16)
        }),
    }
}

fn target_for(kc: &KoblitzCurve, points: &[BinaryPoint], mode: &str, offset: usize) -> BinaryPoint {
    if mode == "planted" {
        return PLANTED
            .iter()
            .fold(BinaryPoint::Infinity, |sum, &i| kc.add(&sum, &points[i]));
    }
    assert_eq!(mode, "ordinary");
    assert!(offset < 4);
    let public = BinaryPoint::Affine {
        x: F2mElement::from_hex("355fb5df7a905f16921eb", 83),
        y: F2mElement::from_hex("5900a390f42d290f1bbe", 83),
    };
    assert!(kc.curve.is_on_curve(&public));
    assert_eq!(kc.mul(&public, &kc.subgroup_order), BinaryPoint::Infinity);
    let four = BigUint::from(4u32);
    let inv_four = (&kc.subgroup_order * BigUint::from(3u32) + BigUint::from(1u32)) / &four;
    let preimage = kc.mul(&public, &inv_four);
    assert_eq!(kc.mul(&preimage, &four), public);
    let one = F2mElement::one(83);
    let zero = F2mElement::zero(83);
    let torsion = [
        BinaryPoint::Infinity,
        BinaryPoint::Affine {
            x: zero.clone(),
            y: one.clone(),
        },
        BinaryPoint::Affine {
            x: one.clone(),
            y: zero,
        },
        BinaryPoint::Affine {
            x: one.clone(),
            y: one,
        },
    ];
    for point in &torsion {
        assert!(kc.curve.is_on_curve(point));
        assert_eq!(kc.mul(point, &four), BinaryPoint::Infinity);
    }
    kc.add(&preimage, &torsion[offset])
}

fn construct_circuit(
    basis: &[F2mElement],
    target: &F2mElement,
    b: &F2mElement,
    table: &impl WideFieldTable,
) -> Circuit {
    let mut circuit = Circuit::new();
    let mut source = Vec::with_capacity(5);
    for summand in 0..5 {
        let mut field = Vec::with_capacity(N);
        for output in 0..N {
            let terms: Vec<_> = basis
                .iter()
                .enumerate()
                .filter(|(_, element)| {
                    let words = element.raw_bits();
                    words[output / 64] >> (output % 64) & 1 != 0
                })
                .map(|(bit, _)| Circuit::input(summand * 18 + bit))
                .collect();
            field.push(circuit.xor_many(terms));
        }
        source.push(field);
    }
    let intermediate: Vec<Vec<Bit>> = (0..3)
        .map(|index| {
            (0..N)
                .map(|bit| Circuit::input(90 + index * N + bit))
                .collect()
        })
        .collect();
    let target = bits_of(target);
    let b = bits_of(b);
    circuit.s3(&source[0], &source[1], &intermediate[0], &b, table);
    circuit.s3(&intermediate[0], &source[2], &intermediate[1], &b, table);
    circuit.s3(&intermediate[1], &source[3], &intermediate[2], &b, table);
    circuit.s3(&intermediate[2], &source[4], &target, &b, table);
    assert_eq!(circuit.outputs.len(), 332);
    circuit
        .units
        .extend(circuit.outputs.iter().copied().map(Bit::not));
    circuit
}

fn compare_expanded(circuit: &Circuit, system: &System512, planted: Mono512) {
    let mut cases = vec![planted];
    let mut state = 0xf6_83_5a_7c_1020_2605u64;
    for _ in 0..8 {
        let mut assignment = Mono512::default();
        for i in 0..INPUTS {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            if state & 1 != 0 {
                assignment.0[i / 64] |= 1u64 << (i % 64);
            }
        }
        cases.push(assignment);
    }
    for (case, assignment) in cases.iter().enumerate() {
        let values = circuit.evaluate(assignment);
        for (i, polynomial) in system.equations.iter().enumerate() {
            assert_eq!(
                values[i],
                polynomial.eval(assignment),
                "case {case}, S3 bit {i}"
            );
        }
    }
}

fn parse_model(path: &str) -> Option<Vec<bool>> {
    let contents = fs::read_to_string(path).expect("solver output");
    if !contents.lines().any(|line| line.trim() == "s SATISFIABLE") {
        return None;
    }
    let mut values = vec![None; INPUTS + 1];
    for line in contents.lines() {
        if !line.starts_with("v ") {
            continue;
        }
        for token in line[2..].split_whitespace() {
            let lit: i32 = token.parse().expect("Kissat model literal");
            if lit == 0 {
                continue;
            }
            let variable = lit.unsigned_abs() as usize;
            if variable <= INPUTS {
                values[variable] = Some(lit > 0);
            }
        }
    }
    Some(
        values
            .into_iter()
            .map(|value| value.unwrap_or(false))
            .collect(),
    )
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        args.len() == 5 || args.len() == 6,
        "usage: probe emit|emit_xor|verify planted|ordinary OFFSET PATH [all|source|none]"
    );
    let action = args[1].as_str();
    let mode = args[2].as_str();
    let offset: usize = args[3].parse().expect("integer torsion offset");
    let path = &args[4];
    assert!(action == "emit" || action == "emit_xor" || action == "verify");
    assert!(args.len() == 5 || action == "emit_xor");
    assert!(mode == "planted" || mode == "ordinary");
    assert!(offset < 4);
    let start = Instant::now();
    let kc = KoblitzCurve::known_n83_k0().expect("pinned K0 curve");
    assert_eq!(kc.label(), "icv1-f2m83-tm6151469093347-debefd74");
    let base = build_standard_subspace_factor_base(&kc, 18).expect("dimension-18 base");
    assert_eq!(base.points.len(), 261_447);
    assert_eq!(base.subspace_basis.len(), 18);
    let table = TwoWordFieldStructure::new(83, &kc.curve.irreducible).expect("field table");
    let target = target_for(&kc, &base.points, mode, offset);
    let setup_ns = start.elapsed().as_nanos();

    if action == "emit" || action == "emit_xor" {
        let begin = Instant::now();
        let system = System512::build(&base.subspace_basis, x_of(&target), &kc.curve.b, 5, &table)
            .expect("exact five-summand system");
        let mut circuit =
            construct_circuit(&base.subspace_basis, x_of(&target), &kc.curve.b, &table);
        let planted = input_assignment(&base.points, &PLANTED, &kc);
        compare_expanded(&circuit, &system, planted);
        let pin = if action == "emit_xor" {
            args.get(5).expect("pin mode").as_str()
        } else if mode == "planted" {
            "source"
        } else {
            "none"
        };
        assert!(matches!(pin, "all" | "source" | "none"));
        assert!(mode == "planted" || pin == "none");
        if mode == "planted" {
            assert!(system.all_vanish(&planted));
            assert_eq!(
                target,
                PLANTED.iter().fold(BinaryPoint::Infinity, |sum, &i| kc
                    .add(&sum, &base.points[i]))
            );
            let fixed_bits = match pin {
                "all" => INPUTS,
                "source" => 90,
                "none" => 0,
                _ => unreachable!(),
            };
            for i in 0..fixed_bits {
                let value = planted.0[i / 64] >> (i % 64) & 1 != 0;
                circuit.units.push(if value {
                    Circuit::input(i)
                } else {
                    Circuit::input(i).not()
                });
            }
        }
        let (clauses, unit_clauses, literals) = if action == "emit_xor" {
            circuit.write_xcnf(path).expect("write extended DIMACS")
        } else {
            circuit.write_cnf(path).expect("write DIMACS")
        };
        println!(
            "{}",
            json!({
                "phase":"emitted", "mode":mode, "offset":offset,
                "encoding":if action == "emit_xor" {"native_xor"} else {"cnf"},
                "pin":pin,
                "curve_id":kc.label(), "target":point_json(&target),
                "source_points":base.points.len(), "input_variables":INPUTS,
                "expanded_equations":system.equations.len(),
                "expanded_term_occurrences":system.monomial_count(),
                "cross_checks":9, "xor_gates":circuit.xors,
                "and_gates":circuit.ands, "cnf_variables":circuit.next_var - 1,
                "cnf_clauses":clauses, "cnf_unit_clauses":unit_clauses,
                "cnf_literals":literals, "cnf_bytes":fs::metadata(path).unwrap().len(),
                "setup_ns":setup_ns, "encode_ns":begin.elapsed().as_nanos(),
                "claim_scope":"decomposition_feasibility_only"
            })
        );
        return;
    }

    let contents = fs::read_to_string(path).expect("solver model");
    let sat = contents.lines().any(|line| line.trim() == "s SATISFIABLE");
    let unsat = contents
        .lines()
        .any(|line| line.trim() == "s UNSATISFIABLE");
    if !sat {
        println!(
            "{}",
            json!({"phase":"verified","mode":mode,"offset":offset,"status":if unsat {"unsat"} else {"no_model"}})
        );
        return;
    }
    let values = parse_model(path).unwrap();
    let mut assignment = Mono512::default();
    for i in 0..INPUTS {
        if values[i + 1] {
            assignment.0[i / 64] |= 1u64 << (i % 64);
        }
    }
    let system = System512::build(&base.subspace_basis, x_of(&target), &kc.curve.b, 5, &table)
        .expect("exact five-summand system");
    let equations_zero = system.all_vanish(&assignment);
    let mut source_codes = Vec::new();
    let mut choices: Vec<Vec<(usize, BinaryPoint)>> = Vec::new();
    for summand in 0..5 {
        let mut x = F2mElement::zero(83);
        let mut code = 0u32;
        for bit in 0..18 {
            let index = summand * 18 + bit;
            if assignment.0[index / 64] >> (index % 64) & 1 != 0 {
                x = x.add(&base.subspace_basis[bit]);
                code |= 1 << bit;
            }
        }
        source_codes.push(code);
        choices.push(
            base.points
                .iter()
                .enumerate()
                .filter(|(_, point)| x_of(point) == &x)
                .map(|(index, point)| (index, point.clone()))
                .collect(),
        );
    }
    let mut group_witness = None;
    let mut usable_witness = None;
    if equations_zero && choices.iter().all(|points| !points.is_empty()) {
        for mask in 0..32 {
            let selected: Vec<_> = choices
                .iter()
                .enumerate()
                .map(|(i, points)| &points[((mask >> i) & 1) % points.len()])
                .collect();
            let sum = selected
                .iter()
                .fold(BinaryPoint::Infinity, |acc, (_, point)| kc.add(&acc, point));
            if sum == target {
                let indices = selected.iter().map(|(index, _)| *index).collect::<Vec<_>>();
                if group_witness.is_none() {
                    group_witness = Some(indices.clone());
                }
                if selected
                    .iter()
                    .all(|(_, point)| kc.mul(point, &BigUint::from(4u32)) != BinaryPoint::Infinity)
                {
                    usable_witness = Some(indices);
                    break;
                }
            }
        }
    }
    let status = if usable_witness.is_some() {
        "verified_relation"
    } else if group_witness.is_some() && mode == "planted" {
        "verified_planted"
    } else if group_witness.is_some() {
        "group_witness_unusable"
    } else {
        "algebraic_model_only"
    };
    println!(
        "{}",
        json!({
            "phase":"verified", "mode":mode, "offset":offset,
            "status":status,
            "equations_zero":equations_zero, "source_codes":source_codes,
            "source_choice_counts":choices.iter().map(Vec::len).collect::<Vec<_>>(),
            "group_witness_indices":group_witness,
            "usable_witness_indices":usable_witness,
            "target":point_json(&target),
            "setup_ns":setup_ns, "claim_scope":"decomposition_feasibility_only"
        })
    );
}
