//! Native complete-point Boolean gate for a rotated high-arity Koblitz descent.
//!
//! The concrete n=83 instance is the exact public `E0: y^2 + xy = x^3 + 1`
//! target independently solved by the strong signed-Frobenius rho campaign in
//! the sibling cryptanalysis repository.  Coordinates use its type-II optimal
//! normal basis.  In that basis Frobenius is a permutation, so ten rotated
//! subspaces with dimensions `[9,9,9,8,8,8,8,8,8,8]` partition all 83 field
//! coordinates.  The circuit joins one finite curve point from each subspace
//! with nine complete affine group-law relations.
//!
//! This module is an experimental relation gate, not an ECDLP speed claim.  A
//! useful attack would still need natural-target yield, relation rank, factor
//! logarithms, target recovery, and fully charged comparison with rho.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use serde::Serialize;

pub const N83: usize = 83;
pub const N83_SUBGROUP_ORDER: u128 = 2_417_851_639_230_796_216_685_689;
pub const N83_RHO_SCALAR: u128 = 467_066_815_623_456_506_232_910;
pub const N83_RHO_WALK_ITERATIONS: u128 = 201_733_439_488;
pub const N83_SLOT_DIMENSIONS: [usize; 10] = [9, 9, 9, 8, 8, 8, 8, 8, 8, 8];

pub type NodeId = u32;
const FALSE: NodeId = 0;
const TRUE: NodeId = 1;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Node {
    Constant(bool),
    Input(usize),
    Xor(NodeId, NodeId),
    And(NodeId, NodeId),
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct DagCounts {
    pub primary_inputs: usize,
    pub xor_gates: usize,
    pub and_gates: usize,
    pub total_nodes: usize,
}

/// Hash-consed XOR/AND circuit.  Node zero is false and node one is true.
#[derive(Clone, Debug)]
pub struct BoolDag {
    nodes: Vec<Node>,
    cache: HashMap<(u8, NodeId, NodeId), NodeId>,
    input_names: Vec<String>,
    input_nodes: Vec<NodeId>,
    node_cap: usize,
}

impl BoolDag {
    pub fn new(node_cap: usize) -> Self {
        assert!(node_cap >= 2, "the DAG cap must include both constants");
        Self {
            nodes: vec![Node::Constant(false), Node::Constant(true)],
            cache: HashMap::new(),
            input_names: Vec::new(),
            input_nodes: Vec::new(),
            node_cap,
        }
    }

    fn append(&mut self, node: Node) -> NodeId {
        assert!(
            self.nodes.len() < self.node_cap,
            "Boolean DAG node cap {} exceeded",
            self.node_cap
        );
        let id = self.nodes.len() as NodeId;
        self.nodes.push(node);
        id
    }

    pub fn input(&mut self, name: impl Into<String>) -> NodeId {
        let name = name.into();
        assert!(
            !self.input_names.iter().any(|old| old == &name),
            "duplicate Boolean input {name}"
        );
        let index = self.input_names.len();
        let id = self.append(Node::Input(index));
        self.input_names.push(name);
        self.input_nodes.push(id);
        id
    }

    pub fn xor(&mut self, mut a: NodeId, mut b: NodeId) -> NodeId {
        if a == b {
            return FALSE;
        }
        if a == FALSE {
            return b;
        }
        if b == FALSE {
            return a;
        }
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        let key = (0, a, b);
        if let Some(&id) = self.cache.get(&key) {
            return id;
        }
        let id = self.append(Node::Xor(a, b));
        self.cache.insert(key, id);
        id
    }

    pub fn and(&mut self, mut a: NodeId, mut b: NodeId) -> NodeId {
        if a == FALSE || b == FALSE {
            return FALSE;
        }
        if a == TRUE {
            return b;
        }
        if b == TRUE || a == b {
            return a;
        }
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        let key = (1, a, b);
        if let Some(&id) = self.cache.get(&key) {
            return id;
        }
        let id = self.append(Node::And(a, b));
        self.cache.insert(key, id);
        id
    }

    pub fn not(&mut self, a: NodeId) -> NodeId {
        self.xor(a, TRUE)
    }

    pub fn or(&mut self, a: NodeId, b: NodeId) -> NodeId {
        let na = self.not(a);
        let nb = self.not(b);
        let neither = self.and(na, nb);
        self.not(neither)
    }

    pub fn all(&mut self, bits: impl IntoIterator<Item = NodeId>) -> NodeId {
        let mut out = TRUE;
        for bit in bits {
            out = self.and(out, bit);
        }
        out
    }

    pub fn any(&mut self, bits: impl IntoIterator<Item = NodeId>) -> NodeId {
        let mut out = FALSE;
        for bit in bits {
            out = self.or(out, bit);
        }
        out
    }

    pub fn implies(&mut self, guard: NodeId, consequence: NodeId) -> NodeId {
        let not_guard = self.not(guard);
        self.or(not_guard, consequence)
    }

    pub fn counts(&self) -> DagCounts {
        DagCounts {
            primary_inputs: self.input_names.len(),
            xor_gates: self
                .nodes
                .iter()
                .filter(|node| matches!(node, Node::Xor(_, _)))
                .count(),
            and_gates: self
                .nodes
                .iter()
                .filter(|node| matches!(node, Node::And(_, _)))
                .count(),
            total_nodes: self.nodes.len(),
        }
    }

    pub fn prefix_blake3(&self) -> String {
        let mut hash = blake3::Hasher::new();
        for node in &self.nodes {
            match *node {
                Node::Constant(value) => {
                    hash.update(if value { b"one\n" } else { b"zero\n" });
                }
                Node::Input(index) => {
                    hash.update(format!("var,{index}\n").as_bytes());
                }
                Node::Xor(a, b) => {
                    hash.update(format!("xor,{a},{b}\n").as_bytes());
                }
                Node::And(a, b) => {
                    hash.update(format!("and,{a},{b}\n").as_bytes());
                }
            }
        }
        hash.finalize().to_hex().to_string()
    }

    fn evaluate(&self, inputs: &[bool]) -> Result<Vec<bool>, String> {
        if inputs.len() != self.input_names.len() {
            return Err(format!(
                "got {} primary inputs, expected {}",
                inputs.len(),
                self.input_names.len()
            ));
        }
        let mut values = vec![false; self.nodes.len()];
        values[1] = true;
        for (id, node) in self.nodes.iter().enumerate().skip(2) {
            values[id] = match *node {
                Node::Constant(value) => value,
                Node::Input(index) => inputs[index],
                Node::Xor(a, b) => values[a as usize] ^ values[b as usize],
                Node::And(a, b) => values[a as usize] & values[b as usize],
            };
        }
        Ok(values)
    }
}

/// Type-II optimal normal-basis arithmetic using coordinates of
/// `gamma_i = zeta^i + zeta^-i`, `1 <= i <= m`, with `zeta^(2m+1)=1`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Type2Onb {
    m: usize,
    ring: usize,
    mask: u128,
}

impl Type2Onb {
    pub fn new(m: usize) -> Result<Self, String> {
        if m == 0 || m >= 128 {
            return Err("type-II ONB implementation requires 1 <= m < 128".into());
        }
        let ring = 2 * m + 1;
        if !is_prime(ring) {
            return Err(format!("2m+1={ring} is not prime"));
        }
        let order = multiplicative_order(2, ring);
        if order != m && order != 2 * m {
            return Err(format!("ord_{ring}(2)={order}, so no type-II ONB"));
        }
        Ok(Self {
            m,
            ring,
            mask: (1u128 << m) - 1,
        })
    }

    pub fn degree(self) -> usize {
        self.m
    }

    pub fn zero(self) -> u128 {
        0
    }

    /// In a type-II ONB, `1 = gamma_1 + ... + gamma_m`.
    pub fn one(self) -> u128 {
        self.mask
    }

    fn fold(self, index: usize) -> usize {
        let index = index % self.ring;
        index.min(self.ring - index)
    }

    pub fn add(self, a: u128, b: u128) -> u128 {
        (a ^ b) & self.mask
    }

    pub fn mul(self, a: u128, b: u128) -> u128 {
        debug_assert_eq!(a & !self.mask, 0);
        debug_assert_eq!(b & !self.mask, 0);
        let mut out = 0u128;
        let mut aa = a;
        while aa != 0 {
            let i = aa.trailing_zeros() as usize + 1;
            aa &= aa - 1;
            let mut bb = b;
            while bb != 0 {
                let j = bb.trailing_zeros() as usize + 1;
                bb &= bb - 1;
                let sum = self.fold(i + j);
                if sum != 0 {
                    out ^= 1u128 << (sum - 1);
                }
                let difference = self.fold(i.abs_diff(j));
                if difference != 0 {
                    out ^= 1u128 << (difference - 1);
                }
            }
        }
        out & self.mask
    }

    pub fn square(self, a: u128) -> u128 {
        self.frobenius(a, 1)
    }

    pub fn frobenius(self, mut a: u128, k: usize) -> u128 {
        for _ in 0..(k % self.m) {
            let mut next = 0u128;
            while a != 0 {
                let i = a.trailing_zeros() as usize + 1;
                a &= a - 1;
                let target = self.fold(2 * i);
                next ^= 1u128 << (target - 1);
            }
            a = next;
        }
        a
    }

    pub fn pow(self, mut base: u128, mut exponent: u128) -> u128 {
        let mut out = self.one();
        while exponent != 0 {
            if exponent & 1 == 1 {
                out = self.mul(out, base);
            }
            base = self.square(base);
            exponent >>= 1;
        }
        out
    }

    pub fn inverse(self, a: u128) -> Option<u128> {
        (a != 0).then(|| self.pow(a, (1u128 << self.m) - 2))
    }

    pub fn trace(self, a: u128) -> bool {
        a.count_ones() & 1 == 1
    }

    pub fn half_trace(self, a: u128) -> Result<u128, String> {
        if self.m.is_multiple_of(2) {
            return Err("half-trace requires odd extension degree".into());
        }
        let mut out = 0u128;
        let mut term = a;
        for _ in 0..=self.m / 2 {
            out ^= term;
            term = self.frobenius(term, 2);
        }
        Ok(out)
    }

    pub fn conjugate_coordinate_order(self) -> Vec<usize> {
        let mut order = Vec::with_capacity(self.m);
        let mut coordinate = 1usize;
        for _ in 0..self.m {
            order.push(coordinate - 1);
            coordinate = self.fold(2 * coordinate);
        }
        assert_eq!(order.len(), self.m);
        let mut unique = order.clone();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique, (0..self.m).collect::<Vec<_>>());
        order
    }
}

fn is_prime(n: usize) -> bool {
    if n < 2 {
        return false;
    }
    let mut divisor = 2usize;
    while divisor * divisor <= n {
        if n.is_multiple_of(divisor) {
            return false;
        }
        divisor += 1;
    }
    true
}

fn multiplicative_order(a: usize, modulus: usize) -> usize {
    let mut value = a % modulus;
    let mut order = 1usize;
    while value != 1 {
        value = value * a % modulus;
        order += 1;
        if order > modulus {
            return 0;
        }
    }
    order
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct OnbPoint {
    pub infinity: bool,
    #[serde(serialize_with = "serialize_u128_hex")]
    pub x: u128,
    #[serde(serialize_with = "serialize_u128_hex")]
    pub y: u128,
}

impl OnbPoint {
    pub const INFINITY: Self = Self {
        infinity: true,
        x: 0,
        y: 0,
    };

    pub const fn affine(x: u128, y: u128) -> Self {
        Self {
            infinity: false,
            x,
            y,
        }
    }
}

fn serialize_u128_hex<S>(value: &u128, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    serializer.serialize_str(&format!("0x{value:x}"))
}

#[derive(Clone, Copy, Debug)]
pub struct OnbCurve {
    pub field: Type2Onb,
}

impl OnbCurve {
    pub fn new(field: Type2Onb) -> Self {
        Self { field }
    }

    pub fn is_on_curve(self, point: OnbPoint) -> bool {
        if point.infinity {
            return point.x == 0 && point.y == 0;
        }
        let f = self.field;
        let lhs = f.add(f.square(point.y), f.mul(point.x, point.y));
        let x2 = f.square(point.x);
        let rhs = f.add(f.mul(x2, point.x), f.one());
        lhs == rhs
    }

    pub fn neg(self, point: OnbPoint) -> OnbPoint {
        if point.infinity {
            point
        } else {
            OnbPoint::affine(point.x, point.y ^ point.x)
        }
    }

    pub fn add(self, p: OnbPoint, q: OnbPoint) -> OnbPoint {
        if p.infinity {
            return q;
        }
        if q.infinity {
            return p;
        }
        let f = self.field;
        if p.x == q.x {
            if p.y ^ q.y == p.x {
                return OnbPoint::INFINITY;
            }
            if p.x == 0 {
                return OnbPoint::INFINITY;
            }
            let lambda = p.x ^ f.mul(p.y, f.inverse(p.x).expect("nonzero x"));
            let x3 = f.add(f.square(lambda), lambda);
            let y3 = f.add(f.square(p.x), f.mul(f.add(lambda, f.one()), x3));
            return OnbPoint::affine(x3, y3);
        }
        let denominator = p.x ^ q.x;
        let lambda = f.mul(p.y ^ q.y, f.inverse(denominator).expect("distinct x"));
        let x3 = f.add(f.add(f.square(lambda), lambda), denominator);
        let y3 = f.add(f.add(f.mul(lambda, p.x ^ x3), x3), p.y);
        OnbPoint::affine(x3, y3)
    }

    pub fn scalar_mul(self, point: OnbPoint, scalar: u128) -> OnbPoint {
        let mut out = OnbPoint::INFINITY;
        for bit in (0..128 - scalar.leading_zeros()).rev() {
            out = self.add(out, out);
            if scalar >> bit & 1 == 1 {
                out = self.add(out, point);
            }
        }
        out
    }

    pub fn slope(self, p: OnbPoint, q: OnbPoint) -> u128 {
        if p.infinity || q.infinity || (p.x == q.x && (p.y ^ q.y == p.x)) {
            return 0;
        }
        let f = self.field;
        if p.x == q.x {
            p.x ^ f.mul(p.y, f.inverse(p.x).expect("nonzero double x"))
        } else {
            f.mul(p.y ^ q.y, f.inverse(p.x ^ q.x).expect("distinct x"))
        }
    }

    pub fn points_with_x(self, x: u128) -> Vec<OnbPoint> {
        let f = self.field;
        if x == 0 {
            return vec![OnbPoint::affine(0, f.one())];
        }
        // Set y=xz.  Then z^2+z=x+1/x^2.
        let inverse = f.inverse(x).expect("nonzero x");
        let c = x ^ f.square(inverse);
        if f.trace(c) {
            return Vec::new();
        }
        let z = f.half_trace(c).expect("n83 is odd");
        debug_assert_eq!(f.square(z) ^ z, c);
        let y = f.mul(x, z);
        vec![OnbPoint::affine(x, y), OnbPoint::affine(x, y ^ x)]
    }
}

pub fn n83_public_generator() -> OnbPoint {
    OnbPoint::affine(
        (0x10dbcu128 << 64) | 0xd28e_dfff_9d7c_a5a0,
        (0x6c82au128 << 64) | 0xd9b0_b7b8_bc36_d66c,
    )
}

pub fn n83_public_target() -> OnbPoint {
    OnbPoint::affine(
        (0x78f3du128 << 64) | 0xaa4b_f524_e3e8_e9c7,
        (0x2e08u128 << 64) | 0xb1e3_a8d4_3e5b_ecd4,
    )
}

#[derive(Clone, Debug)]
struct CircuitField {
    degree: usize,
    ring: usize,
}

impl CircuitField {
    fn new(degree: usize) -> Result<Self, String> {
        Type2Onb::new(degree)?;
        Ok(Self {
            degree,
            ring: 2 * degree + 1,
        })
    }

    fn fold(&self, index: usize) -> usize {
        let index = index % self.ring;
        index.min(self.ring - index)
    }

    fn constant(&self, value: u128) -> Vec<NodeId> {
        (0..self.degree)
            .map(|bit| if value >> bit & 1 == 1 { TRUE } else { FALSE })
            .collect()
    }

    fn input(&self, dag: &mut BoolDag, prefix: &str) -> Vec<NodeId> {
        (0..self.degree)
            .map(|bit| dag.input(format!("{prefix}_{bit}")))
            .collect()
    }

    fn add(&self, dag: &mut BoolDag, a: &[NodeId], b: &[NodeId]) -> Vec<NodeId> {
        debug_assert_eq!(a.len(), self.degree);
        debug_assert_eq!(b.len(), self.degree);
        a.iter().zip(b).map(|(&x, &y)| dag.xor(x, y)).collect()
    }

    fn mul(&self, dag: &mut BoolDag, a: &[NodeId], b: &[NodeId]) -> Vec<NodeId> {
        debug_assert_eq!(a.len(), self.degree);
        debug_assert_eq!(b.len(), self.degree);
        let mut out = vec![FALSE; self.degree];
        for (i0, &a_bit) in a.iter().enumerate() {
            if a_bit == FALSE {
                continue;
            }
            let i = i0 + 1;
            for (j0, &b_bit) in b.iter().enumerate() {
                if b_bit == FALSE {
                    continue;
                }
                let term = dag.and(a_bit, b_bit);
                if term == FALSE {
                    continue;
                }
                let sum = self.fold(i + j0 + 1);
                if sum != 0 {
                    out[sum - 1] = dag.xor(out[sum - 1], term);
                }
                let difference = self.fold(i.abs_diff(j0 + 1));
                if difference != 0 {
                    out[difference - 1] = dag.xor(out[difference - 1], term);
                }
            }
        }
        out
    }

    fn square(&self, a: &[NodeId]) -> Vec<NodeId> {
        debug_assert_eq!(a.len(), self.degree);
        let mut out = vec![FALSE; self.degree];
        for (i0, &bit) in a.iter().enumerate() {
            let target = self.fold(2 * (i0 + 1));
            out[target - 1] = bit;
        }
        out
    }

    fn zero(&self, dag: &mut BoolDag, a: &[NodeId]) -> NodeId {
        let negated: Vec<_> = a.iter().map(|&bit| dag.not(bit)).collect();
        dag.all(negated)
    }

    fn equal(&self, dag: &mut BoolDag, a: &[NodeId], b: &[NodeId]) -> NodeId {
        let difference = self.add(dag, a, b);
        self.zero(dag, &difference)
    }
}

#[derive(Clone, Debug)]
struct CircuitPoint {
    infinity: NodeId,
    x: Vec<NodeId>,
    y: Vec<NodeId>,
}

fn point_input(dag: &mut BoolDag, field: &CircuitField, label: &str) -> CircuitPoint {
    CircuitPoint {
        infinity: dag.input(format!("{label}_o")),
        x: field.input(dag, &format!("{label}_x")),
        y: field.input(dag, &format!("{label}_y")),
    }
}

fn point_constant(field: &CircuitField, point: OnbPoint) -> CircuitPoint {
    CircuitPoint {
        infinity: if point.infinity { TRUE } else { FALSE },
        x: field.constant(point.x),
        y: field.constant(point.y),
    }
}

fn point_valid(dag: &mut BoolDag, field: &CircuitField, point: &CircuitPoint) -> NodeId {
    let zero_x = field.zero(dag, &point.x);
    let zero_y = field.zero(dag, &point.y);
    let canonical_infinity = dag.and(zero_x, zero_y);
    let infinity_valid = dag.implies(point.infinity, canonical_infinity);

    let y2 = field.square(&point.y);
    let xy = field.mul(dag, &point.x, &point.y);
    let lhs = field.add(dag, &y2, &xy);
    let x2 = field.square(&point.x);
    let x3 = field.mul(dag, &x2, &point.x);
    let rhs = field.add(dag, &x3, &field.constant((1u128 << field.degree) - 1));
    let equation = field.equal(dag, &lhs, &rhs);
    let finite = dag.not(point.infinity);
    let affine_valid = dag.implies(finite, equation);
    dag.and(infinity_valid, affine_valid)
}

fn equal_points(
    dag: &mut BoolDag,
    field: &CircuitField,
    a: &CircuitPoint,
    b: &CircuitPoint,
) -> NodeId {
    let infinity_difference = dag.xor(a.infinity, b.infinity);
    let same_infinity = dag.not(infinity_difference);
    let same_x = field.equal(dag, &a.x, &b.x);
    let same_y = field.equal(dag, &a.y, &b.y);
    dag.all([same_infinity, same_x, same_y])
}

/// Complete affine addition constraint for `p + q = r` on E0.
fn addition_constraint(
    dag: &mut BoolDag,
    field: &CircuitField,
    p: &CircuitPoint,
    q: &CircuitPoint,
    r: &CircuitPoint,
    lambda: &[NodeId],
    check_point_validity: bool,
) -> NodeId {
    let validity = if check_point_validity {
        let p_valid = point_valid(dag, field, p);
        let q_valid = point_valid(dag, field, q);
        let r_valid = point_valid(dag, field, r);
        dag.all([p_valid, q_valid, r_valid])
    } else {
        TRUE
    };

    let not_po = dag.not(p.infinity);
    let not_qo = dag.not(q.infinity);
    let finite = dag.and(not_po, not_qo);
    let same_x = field.equal(dag, &p.x, &q.x);
    let same_y = field.equal(dag, &p.y, &q.y);
    let py_plus_qy = field.add(dag, &p.y, &q.y);
    let inverse = field.equal(dag, &py_plus_qy, &p.x);
    let nonzero_x = {
        let zero_x = field.zero(dag, &p.x);
        dag.not(zero_x)
    };

    let copy_q_branch = p.infinity;
    let copy_p_branch = dag.and(not_po, q.infinity);
    let inverse_branch = dag.all([finite, same_x, inverse]);
    let not_inverse = dag.not(inverse);
    let double_branch = dag.all([finite, same_x, not_inverse, same_y, nonzero_x]);
    let not_same_x = dag.not(same_x);
    let generic_branch = dag.and(finite, not_same_x);
    let cover = dag.any([
        copy_q_branch,
        copy_p_branch,
        inverse_branch,
        double_branch,
        generic_branch,
    ]);

    let copy_q_equal = equal_points(dag, field, r, q);
    let copy_q = dag.implies(copy_q_branch, copy_q_equal);
    let copy_p_equal = equal_points(dag, field, r, p);
    let copy_p = dag.implies(copy_p_branch, copy_p_equal);
    let zero_rx = field.zero(dag, &r.x);
    let zero_ry = field.zero(dag, &r.y);
    let inverse_value = dag.all([r.infinity, zero_rx, zero_ry]);
    let inverse_result = dag.implies(inverse_branch, inverse_value);

    let lambda_times_px = field.mul(dag, lambda, &p.x);
    let px2_plus_py = field.add(dag, &field.square(&p.x), &p.y);
    let double_slope = field.equal(dag, &lambda_times_px, &px2_plus_py);
    let lambda2_plus_lambda = field.add(dag, &field.square(lambda), lambda);
    let double_x = field.equal(dag, &r.x, &lambda2_plus_lambda);
    let lambda_plus_one = field.add(dag, lambda, &field.constant((1u128 << field.degree) - 1));
    let lambda_one_times_rx = field.mul(dag, &lambda_plus_one, &r.x);
    let double_y_value = field.add(dag, &field.square(&p.x), &lambda_one_times_rx);
    let double_y = field.equal(dag, &r.y, &double_y_value);
    let not_ro = dag.not(r.infinity);
    let double_value = dag.all([not_ro, double_slope, double_x, double_y]);
    let double_result = dag.implies(double_branch, double_value);

    let px_plus_qx = field.add(dag, &p.x, &q.x);
    let generic_slope_lhs = field.mul(dag, lambda, &px_plus_qx);
    let generic_slope_rhs = field.add(dag, &p.y, &q.y);
    let generic_slope = field.equal(dag, &generic_slope_lhs, &generic_slope_rhs);
    let generic_x_value = field.add(dag, &lambda2_plus_lambda, &px_plus_qx);
    let generic_x = field.equal(dag, &r.x, &generic_x_value);
    let px_plus_rx = field.add(dag, &p.x, &r.x);
    let lambda_times_px_rx = field.mul(dag, lambda, &px_plus_rx);
    let generic_y_prefix = field.add(dag, &lambda_times_px_rx, &r.x);
    let generic_y_value = field.add(dag, &generic_y_prefix, &p.y);
    let generic_y = field.equal(dag, &r.y, &generic_y_value);
    let generic_value = dag.all([not_ro, generic_slope, generic_x, generic_y]);
    let generic_result = dag.implies(generic_branch, generic_value);

    dag.all([
        validity,
        cover,
        copy_q,
        copy_p,
        inverse_result,
        double_result,
        generic_result,
    ])
}

#[derive(Clone, Debug, Serialize)]
pub struct RotatedSlot {
    pub dimension: usize,
    /// Zero-based coordinate positions in `gamma_1,...,gamma_83` order.
    pub coordinate_positions: Vec<usize>,
}

pub fn n83_rotated_slots() -> Vec<RotatedSlot> {
    let field = Type2Onb::new(N83).expect("n83 type-II ONB");
    let conjugates = field.conjugate_coordinate_order();
    let mut slots = Vec::with_capacity(10);
    let mut all = Vec::with_capacity(N83);
    for (slot, &dimension) in N83_SLOT_DIMENSIONS.iter().enumerate() {
        let positions: Vec<_> = (0..dimension)
            .map(|row| conjugates[10 * row + slot])
            .collect();
        all.extend_from_slice(&positions);
        slots.push(RotatedSlot {
            dimension,
            coordinate_positions: positions,
        });
    }
    let mut sorted = all;
    sorted.sort_unstable();
    assert_eq!(sorted, (0..N83).collect::<Vec<_>>());
    slots
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum N83TargetKind {
    PublicRhoTarget,
    DeterministicPlanted,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum N83CircuitVariant {
    /// Every addition edge independently checks both inputs and its result.
    CompletePerEdgeValidity,
    /// Check the ten factors once, then derive intermediate validity from the
    /// complete group law.  Exhaustive small-field tests pin the implication.
    InductiveFactorValidity,
}

#[derive(Clone, Debug)]
pub struct RotatedChain {
    pub dag: BoolDag,
    pub output: NodeId,
    pub target: OnbPoint,
    pub target_kind: N83TargetKind,
    pub circuit_variant: N83CircuitVariant,
    pub slots: Vec<RotatedSlot>,
    factor_selectors: Vec<Vec<NodeId>>,
    factor_y: Vec<Vec<NodeId>>,
    intermediate_points: Vec<CircuitPoint>,
    slopes: Vec<Vec<NodeId>>,
    planted_factors: Option<Vec<(u128, OnbPoint)>>,
}

impl RotatedChain {
    pub fn counts(&self) -> DagCounts {
        self.dag.counts()
    }

    pub fn dag_prefix_blake3(&self) -> String {
        self.dag.prefix_blake3()
    }

    pub fn planted_factors(&self) -> Option<Vec<OnbPoint>> {
        self.planted_factors
            .as_ref()
            .map(|rows| rows.iter().map(|row| row.1).collect())
    }

    pub fn known_planted_model(&self) -> Result<Vec<bool>, String> {
        let planted = self
            .planted_factors
            .as_ref()
            .ok_or("the public-target chain has no known decomposition")?;
        let mut inputs = vec![false; self.dag.input_names.len()];
        let input_index: HashMap<NodeId, usize> = self
            .dag
            .input_nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(index, node)| (node, index))
            .collect();
        let mut set_node = |node: NodeId, value: bool| -> Result<(), String> {
            let index = *input_index
                .get(&node)
                .ok_or_else(|| format!("node {node} is not a primary input"))?;
            inputs[index] = value;
            Ok(())
        };

        for (slot, &(mask, point)) in planted.iter().enumerate() {
            for (bit, &node) in self.factor_selectors[slot].iter().enumerate() {
                set_node(node, mask >> bit & 1 == 1)?;
            }
            for (bit, &node) in self.factor_y[slot].iter().enumerate() {
                set_node(node, point.y >> bit & 1 == 1)?;
            }
        }

        let curve = OnbCurve::new(Type2Onb::new(N83)?);
        let mut accumulator = planted[0].1;
        for edge in 0..9 {
            let right = planted[edge + 1].1;
            let slope = curve.slope(accumulator, right);
            for (bit, &node) in self.slopes[edge].iter().enumerate() {
                set_node(node, slope >> bit & 1 == 1)?;
            }
            accumulator = curve.add(accumulator, right);
            if edge < 8 {
                let circuit = &self.intermediate_points[edge];
                set_node(circuit.infinity, accumulator.infinity)?;
                for (bit, &node) in circuit.x.iter().enumerate() {
                    set_node(node, accumulator.x >> bit & 1 == 1)?;
                }
                for (bit, &node) in circuit.y.iter().enumerate() {
                    set_node(node, accumulator.y >> bit & 1 == 1)?;
                }
            }
        }
        if accumulator != self.target {
            return Err("planted factors do not sum to the circuit target".into());
        }
        let values = self.dag.evaluate(&inputs)?;
        if !values[self.output as usize] {
            return Err("known planted model does not satisfy the circuit".into());
        }
        Ok(values)
    }

    pub fn decode_dimacs_model(&self, path: &Path) -> Result<DecodedRelation, String> {
        let values = parse_dimacs_model(path, self.dag.nodes.len())?;
        for (id, node) in self.dag.nodes.iter().enumerate() {
            let expected = match *node {
                Node::Constant(value) => value,
                Node::Input(_) => continue,
                Node::Xor(a, b) => values[a as usize] ^ values[b as usize],
                Node::And(a, b) => values[a as usize] & values[b as usize],
            };
            if values[id] != expected {
                return Err(format!("model violates DAG node {id}"));
            }
        }
        if !values[self.output as usize] {
            return Err("model sets the required chain output false".into());
        }

        let curve = OnbCurve::new(Type2Onb::new(N83)?);
        let mut factors = Vec::with_capacity(10);
        for slot in 0..10 {
            let mut x = 0u128;
            for (&position, &node) in self.slots[slot]
                .coordinate_positions
                .iter()
                .zip(&self.factor_selectors[slot])
            {
                if values[node as usize] {
                    x |= 1u128 << position;
                }
            }
            let mut y = 0u128;
            for (bit, &node) in self.factor_y[slot].iter().enumerate() {
                if values[node as usize] {
                    y |= 1u128 << bit;
                }
            }
            let point = OnbPoint::affine(x, y);
            if !curve.is_on_curve(point) {
                return Err(format!("decoded factor {slot} is off curve"));
            }
            factors.push(point);
        }
        let sum = factors
            .iter()
            .copied()
            .fold(OnbPoint::INFINITY, |acc, point| curve.add(acc, point));
        if sum != self.target {
            return Err("decoded factors do not sum to the fixed target".into());
        }
        Ok(DecodedRelation {
            target: self.target,
            factors,
            sum,
            dag_checked: true,
            curve_checked: true,
            group_sum_checked: true,
        })
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct DecodedRelation {
    pub target: OnbPoint,
    pub factors: Vec<OnbPoint>,
    pub sum: OnbPoint,
    pub dag_checked: bool,
    pub curve_checked: bool,
    pub group_sum_checked: bool,
}

pub fn build_n83_rotated_chain(
    target_kind: N83TargetKind,
    node_cap: usize,
) -> Result<RotatedChain, String> {
    build_n83_rotated_chain_variant(
        target_kind,
        N83CircuitVariant::CompletePerEdgeValidity,
        node_cap,
    )
}

pub fn build_n83_rotated_chain_variant(
    target_kind: N83TargetKind,
    circuit_variant: N83CircuitVariant,
    node_cap: usize,
) -> Result<RotatedChain, String> {
    let numeric_field = Type2Onb::new(N83)?;
    let curve = OnbCurve::new(numeric_field);
    let slots = n83_rotated_slots();
    let planted = match target_kind {
        N83TargetKind::PublicRhoTarget => None,
        N83TargetKind::DeterministicPlanted => Some(deterministic_planted_factors(curve, &slots)?),
    };
    let target = if let Some(rows) = &planted {
        rows.iter()
            .map(|row| row.1)
            .fold(OnbPoint::INFINITY, |acc, point| curve.add(acc, point))
    } else {
        n83_public_target()
    };
    if !curve.is_on_curve(target) || target.infinity {
        return Err("fixed target is not a finite n83 curve point".into());
    }

    let field = CircuitField::new(N83)?;
    let mut dag = BoolDag::new(node_cap);
    let mut factors = Vec::with_capacity(10);
    let mut factor_selectors = Vec::with_capacity(10);
    let mut factor_y = Vec::with_capacity(10);
    for (slot_index, slot) in slots.iter().enumerate() {
        let selectors: Vec<_> = (0..slot.dimension)
            .map(|bit| dag.input(format!("F{slot_index}_a_{bit}")))
            .collect();
        let mut x = vec![FALSE; N83];
        for (&position, &selector) in slot.coordinate_positions.iter().zip(&selectors) {
            x[position] = selector;
        }
        let y = field.input(&mut dag, &format!("F{slot_index}_y"));
        factors.push(CircuitPoint {
            infinity: FALSE,
            x,
            y: y.clone(),
        });
        factor_selectors.push(selectors);
        factor_y.push(y);
    }

    let intermediate_points: Vec<_> = (2..10)
        .map(|index| point_input(&mut dag, &field, &format!("S{index}")))
        .collect();
    let slopes: Vec<_> = (0..9)
        .map(|index| field.input(&mut dag, &format!("L{index}")))
        .collect();
    let target_circuit = point_constant(&field, target);

    let mut constraints = Vec::with_capacity(19);
    if circuit_variant == N83CircuitVariant::InductiveFactorValidity {
        constraints.extend(
            factors
                .iter()
                .map(|factor| point_valid(&mut dag, &field, factor)),
        );
    }
    for edge in 0..9 {
        let left = if edge == 0 {
            &factors[0]
        } else {
            &intermediate_points[edge - 1]
        };
        let right = &factors[edge + 1];
        let output = if edge == 8 {
            &target_circuit
        } else {
            &intermediate_points[edge]
        };
        constraints.push(addition_constraint(
            &mut dag,
            &field,
            left,
            right,
            output,
            &slopes[edge],
            circuit_variant == N83CircuitVariant::CompletePerEdgeValidity,
        ));
    }
    let output = dag.all(constraints);
    let expected_inputs = N83 + 10 * N83 + 8 * (2 * N83 + 1) + 9 * N83;
    if dag.input_names.len() != expected_inputs {
        return Err(format!(
            "primary-input accounting drift: got {}, expected {expected_inputs}",
            dag.input_names.len()
        ));
    }
    Ok(RotatedChain {
        dag,
        output,
        target,
        target_kind,
        circuit_variant,
        slots,
        factor_selectors,
        factor_y,
        intermediate_points,
        slopes,
        planted_factors: planted,
    })
}

fn deterministic_planted_factors(
    curve: OnbCurve,
    slots: &[RotatedSlot],
) -> Result<Vec<(u128, OnbPoint)>, String> {
    let mut factors = Vec::with_capacity(slots.len());
    for (slot_index, slot) in slots.iter().enumerate() {
        let limit = 1u128 << slot.dimension;
        let start = 1 + (slot_index as u128 * 37) % (limit - 1);
        let mut found = None;
        for offset in 0..limit - 1 {
            let mask = 1 + (start - 1 + offset) % (limit - 1);
            let x = slot
                .coordinate_positions
                .iter()
                .enumerate()
                .fold(0u128, |acc, (bit, &position)| {
                    acc | (((mask >> bit) & 1) << position)
                });
            if let Some(point) = curve.points_with_x(x).into_iter().next() {
                found = Some((mask, point));
                break;
            }
        }
        factors.push(found.ok_or_else(|| format!("slot {slot_index} has no lift"))?);
    }
    Ok(factors)
}

#[derive(Clone, Debug, Serialize)]
pub struct CnfReceipt {
    pub variables: usize,
    pub clauses: usize,
    pub bytes: u64,
    pub blake3: String,
    pub relation_output_literal: u32,
    pub dag: DagCounts,
}

struct HashingWriter<W> {
    inner: W,
    hasher: blake3::Hasher,
    bytes: u64,
}

impl<W: Write> HashingWriter<W> {
    fn new(inner: W) -> Self {
        Self {
            inner,
            hasher: blake3::Hasher::new(),
            bytes: 0,
        }
    }
}

impl<W: Write> Write for HashingWriter<W> {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buffer)?;
        self.hasher.update(&buffer[..written]);
        self.bytes += written as u64;
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

pub fn write_dimacs(chain: &RotatedChain, path: &Path) -> Result<CnfReceipt, String> {
    let counts = chain.dag.counts();
    let clauses = 3 + 4 * counts.xor_gates + 3 * counts.and_gates;
    let file = File::create(path).map_err(|error| format!("create {}: {error}", path.display()))?;
    let buffered = BufWriter::with_capacity(1 << 20, file);
    let mut writer = HashingWriter::new(buffered);
    writeln!(writer, "p cnf {} {clauses}", counts.total_nodes)
        .map_err(|error| error.to_string())?;
    writeln!(writer, "-1 0").map_err(|error| error.to_string())?;
    writeln!(writer, "2 0").map_err(|error| error.to_string())?;
    for (id, node) in chain.dag.nodes.iter().enumerate().skip(2) {
        let z = id as i64 + 1;
        let clauses: Vec<[i64; 3]> = match *node {
            Node::Constant(_) | Node::Input(_) => continue,
            Node::Xor(a, b) => {
                let a = a as i64 + 1;
                let b = b as i64 + 1;
                vec![[-a, -b, -z], [a, b, -z], [a, -b, z], [-a, b, z]]
            }
            Node::And(a, b) => {
                let a = a as i64 + 1;
                let b = b as i64 + 1;
                vec![[-a, -b, z], [a, -z, 0], [b, -z, 0]]
            }
        };
        for clause in clauses {
            if clause[2] == 0 {
                writeln!(writer, "{} {} 0", clause[0], clause[1])
                    .map_err(|error| error.to_string())?;
            } else {
                writeln!(writer, "{} {} {} 0", clause[0], clause[1], clause[2])
                    .map_err(|error| error.to_string())?;
            }
        }
    }
    writeln!(writer, "{} 0", chain.output + 1).map_err(|error| error.to_string())?;
    writer.flush().map_err(|error| error.to_string())?;
    let digest = writer.hasher.finalize().to_hex().to_string();
    Ok(CnfReceipt {
        variables: counts.total_nodes,
        clauses,
        bytes: writer.bytes,
        blake3: digest,
        relation_output_literal: chain.output + 1,
        dag: counts,
    })
}

pub fn write_dimacs_model(values: &[bool], path: &Path) -> Result<(), String> {
    let file = File::create(path).map_err(|error| format!("create {}: {error}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "s SATISFIABLE").map_err(|error| error.to_string())?;
    for (chunk_index, chunk) in values.chunks(16).enumerate() {
        write!(writer, "v").map_err(|error| error.to_string())?;
        for (offset, &value) in chunk.iter().enumerate() {
            let variable = chunk_index * 16 + offset + 1;
            let literal = if value {
                variable as i64
            } else {
                -(variable as i64)
            };
            write!(writer, " {literal}").map_err(|error| error.to_string())?;
        }
        writeln!(writer, " 0").map_err(|error| error.to_string())?;
    }
    writer.flush().map_err(|error| error.to_string())
}

fn parse_dimacs_model(path: &Path, variables: usize) -> Result<Vec<bool>, String> {
    let file = File::open(path).map_err(|error| format!("open {}: {error}", path.display()))?;
    let mut assignments = vec![None; variables];
    let mut satisfiable = false;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|error| error.to_string())?;
        if line.starts_with('s') && line.contains("SATISFIABLE") && !line.contains("UNSAT") {
            satisfiable = true;
        }
        if !line.starts_with('v') {
            continue;
        }
        for token in line[1..].split_whitespace() {
            let literal: i64 = token
                .parse()
                .map_err(|_| format!("invalid model literal {token:?}"))?;
            if literal == 0 {
                continue;
            }
            let variable = literal.unsigned_abs() as usize;
            if variable == 0 || variable > variables {
                return Err(format!("model variable {variable} outside 1..={variables}"));
            }
            let value = literal > 0;
            if assignments[variable - 1]
                .replace(value)
                .is_some_and(|old| old != value)
            {
                return Err(format!("model assigns variable {variable} twice"));
            }
        }
    }
    if !satisfiable {
        return Err("model file does not declare SATISFIABLE".into());
    }
    assignments
        .into_iter()
        .enumerate()
        .map(|(index, value)| value.ok_or_else(|| format!("model omits variable {}", index + 1)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set_point_inputs(
        inputs: &mut [bool],
        by_node: &HashMap<NodeId, usize>,
        circuit: &CircuitPoint,
        point: OnbPoint,
    ) {
        inputs[by_node[&circuit.infinity]] = point.infinity;
        for (bit, &node) in circuit.x.iter().enumerate() {
            inputs[by_node[&node]] = point.x >> bit & 1 == 1;
        }
        for (bit, &node) in circuit.y.iter().enumerate() {
            inputs[by_node[&node]] = point.y >> bit & 1 == 1;
        }
    }

    #[test]
    fn type2_onb_field_laws_are_exhaustive_at_n3() {
        let field = Type2Onb::new(3).expect("type-II ONB at n3");
        assert_eq!(field.one(), 7);
        for a in 0..8 {
            assert_eq!(field.mul(a, field.one()), a);
            assert_eq!(field.square(a), field.mul(a, a));
            assert_eq!(field.frobenius(a, 3), a);
            if a != 0 {
                assert_eq!(
                    field.mul(a, field.inverse(a).expect("inverse")),
                    field.one()
                );
            }
            for b in 0..8 {
                assert_eq!(field.mul(a, b), field.mul(b, a));
                for c in 0..8 {
                    assert_eq!(field.mul(a, b ^ c), field.mul(a, b) ^ field.mul(a, c));
                    assert_eq!(field.mul(field.mul(a, b), c), field.mul(a, field.mul(b, c)));
                }
            }
        }
    }

    #[test]
    fn exact_public_n83_target_matches_rho_scalar() {
        let field = Type2Onb::new(N83).expect("n83 ONB");
        let curve = OnbCurve::new(field);
        let generator = n83_public_generator();
        let target = n83_public_target();
        assert!(curve.is_on_curve(generator));
        assert!(curve.is_on_curve(target));
        assert_eq!(
            curve.scalar_mul(generator, N83_SUBGROUP_ORDER),
            OnbPoint::INFINITY
        );
        assert_eq!(curve.scalar_mul(generator, N83_RHO_SCALAR), target);
    }

    #[test]
    fn n83_rotated_slots_partition_every_coordinate() {
        let slots = n83_rotated_slots();
        assert_eq!(slots.len(), 10);
        assert_eq!(slots.iter().map(|slot| slot.dimension).sum::<usize>(), N83);
        assert_eq!(
            slots.iter().map(|slot| slot.dimension).collect::<Vec<_>>(),
            N83_SLOT_DIMENSIONS
        );
    }

    #[test]
    fn validity_free_addition_is_sound_and_complete_for_valid_n3_points() {
        let numeric = Type2Onb::new(3).expect("n3 ONB");
        let curve = OnbCurve::new(numeric);
        let mut points = vec![OnbPoint::INFINITY];
        for x in 0..8 {
            for y in 0..8 {
                let point = OnbPoint::affine(x, y);
                if curve.is_on_curve(point) {
                    points.push(point);
                }
            }
        }

        let field = CircuitField::new(3).expect("circuit field");
        let mut dag = BoolDag::new(100_000);
        let p = point_input(&mut dag, &field, "p");
        let q = point_input(&mut dag, &field, "q");
        let r = point_input(&mut dag, &field, "r");
        let lambda = field.input(&mut dag, "lambda");
        let output = addition_constraint(&mut dag, &field, &p, &q, &r, &lambda, false);
        let by_node: HashMap<_, _> = dag
            .input_nodes
            .iter()
            .copied()
            .enumerate()
            .map(|(index, node)| (node, index))
            .collect();

        for &left in &points {
            for &right in &points {
                let expected = curve.add(left, right);
                let mut expected_has_witness = false;
                for &candidate in &points {
                    for slope in 0..8u128 {
                        let mut inputs = vec![false; dag.input_names.len()];
                        set_point_inputs(&mut inputs, &by_node, &p, left);
                        set_point_inputs(&mut inputs, &by_node, &q, right);
                        set_point_inputs(&mut inputs, &by_node, &r, candidate);
                        for (bit, &node) in lambda.iter().enumerate() {
                            inputs[by_node[&node]] = slope >> bit & 1 == 1;
                        }
                        if dag.evaluate(&inputs).expect("evaluate")[output as usize] {
                            assert_eq!(candidate, expected);
                            expected_has_witness = true;
                        }
                    }
                }
                assert!(
                    expected_has_witness,
                    "missing witness for {left:?}+{right:?}"
                );
            }
        }
    }

    #[test]
    fn planted_chain_has_a_replayable_complete_model() {
        // The n83 chain is intentionally the real circuit.  This catches
        // representation or complete-addition drift before a solver run.
        let chain = build_n83_rotated_chain(N83TargetKind::DeterministicPlanted, 8_000_000)
            .expect("build planted chain");
        let values = chain.known_planted_model().expect("known model");
        assert!(values[chain.output as usize]);
        assert!(chain.counts().total_nodes < 8_000_000);
    }

    #[test]
    fn inductive_planted_chain_replays_and_is_strictly_smaller() {
        let complete = build_n83_rotated_chain(N83TargetKind::DeterministicPlanted, 8_000_000)
            .expect("complete chain");
        let inductive = build_n83_rotated_chain_variant(
            N83TargetKind::DeterministicPlanted,
            N83CircuitVariant::InductiveFactorValidity,
            8_000_000,
        )
        .expect("inductive chain");
        let values = inductive.known_planted_model().expect("known model");
        assert!(values[inductive.output as usize]);
        assert!(inductive.counts().total_nodes < complete.counts().total_nodes);
        assert!(inductive.counts().and_gates < complete.counts().and_gates);
    }
}
