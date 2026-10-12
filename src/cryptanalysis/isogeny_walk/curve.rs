//! Short Weierstrass curves `y² = x³ + ax + b` over [`Field`]: points, the
//! `j`-invariant, isomorphism tests, the canonical model a walk records,
//! and the order audit.

use num_bigint::BigUint;

use super::field::{Fe, Field};

/// A model `y² = x³ + ax + b`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Model {
    pub a: Fe,
    pub b: Fe,
}

/// The name of the rule [`canonical_model`] applies, recorded in every
/// walk so a reader can rebuild a model from its `j`-invariant and order.
pub const CANONICAL_RULE: &str =
    "icwalk-canon/v1: the registered model when the curve is the walk's \
root; else an a=-3 model with the least b when one exists over F_p; else \
y^2=x^3+3j(1728-j)c^2*x+2j(1728-j)^2c^3 with c the least of {1, least \
non-residue} that gives an F_p-isomorphic model";

impl Model {
    pub fn rhs(&self, f: &Field, x: &Fe) -> Fe {
        let x2 = f.sqr(x);
        f.add(&f.mul(&f.add(&x2, &self.a), x), &self.b)
    }

    /// `4a³ + 27b²`.
    pub fn discriminant(&self, f: &Field) -> Fe {
        let a3 = f.mul(&f.sqr(&self.a), &self.a);
        f.add(
            &f.mul(&f.from_u64(4), &a3),
            &f.mul(&f.from_u64(27), &f.sqr(&self.b)),
        )
    }

    /// `j = 1728·4a³ / (4a³ + 27b²)`; `None` when singular.
    pub fn j(&self, f: &Field) -> Option<Fe> {
        let a3 = f.mul(&f.sqr(&self.a), &self.a);
        let num = f.mul(&f.from_u64(1728 * 4), &a3);
        f.div(&num, &self.discriminant(f))
    }

    pub fn on_curve(&self, f: &Field, x: &Fe, y: &Fe) -> bool {
        f.sqr(y) == self.rhs(f, x)
    }
}

/// `s = u²` with `(s²a, s³b) = (a2, b2)` and `s` a square, i.e. the
/// isomorphism `(x, y) ↦ (u²x, u³y)` from `m` to `m2` over `F_p`; `None`
/// when the models are not `F_p`-isomorphic.  Needs `ab ≠ 0`.
pub fn isomorphism(f: &Field, m: &Model, m2: &Model) -> Option<Fe> {
    let s = f.div(&f.mul(&m2.b, &m.a), &f.mul(&m2.a, &m.b))?;
    let s2 = f.sqr(&s);
    let ok =
        f.mul(&s2, &m.a) == m2.a && f.mul(&f.mul(&s2, &s), &m.b) == m2.b && f.legendre(&s) == 1;
    ok.then_some(s)
}

/// The `a = −3` models of `m` over `F_p`, as `(b, s = u²)`, least `b`
/// first.
pub fn a_minus_3_models(f: &Field, m: &Model) -> Vec<(Fe, Fe)> {
    let w = match f.div(&f.from_i64(-3), &m.a) {
        Some(w) => w,
        None => return Vec::new(),
    };
    let Some(r) = f.sqrt(&w) else {
        return Vec::new();
    };
    let mut out: Vec<(Fe, Fe)> = [r, f.neg(&r)]
        .into_iter()
        .filter(|s| f.legendre(s) == 1)
        .map(|s| (f.mul(&f.mul(&f.sqr(&s), &s), &m.b), s))
        .collect();
    out.sort_by_key(|(b, _)| f.to_big(b));
    out.dedup_by_key(|(b, _)| *b);
    out
}

/// The canonical model of `m`'s `F_p`-isomorphism class, and the `s = u²`
/// mapping `m` onto it (see [`CANONICAL_RULE`]).  `None` for `j ∈ {0,
/// 1728}`, which no curve of an ordinary class with a large discriminant
/// reaches.
pub fn canonical_model(f: &Field, m: &Model) -> Option<(Model, Fe)> {
    if f.is_zero(&m.a) || f.is_zero(&m.b) {
        return None;
    }
    if let Some((b, s)) = a_minus_3_models(f, m).into_iter().next() {
        return Some((
            Model {
                a: f.from_i64(-3),
                b,
            },
            s,
        ));
    }
    let j = m.j(f)?;
    let k = f.sub(&f.from_u64(1728), &j);
    let a0 = f.mul(&f.from_u64(3), &f.mul(&j, &k));
    let b0 = f.mul(&f.from_u64(2), &f.mul(&j, &f.sqr(&k)));
    for c in [f.one(), f.from_u64(f.least_nonresidue())] {
        let cand = Model {
            a: f.mul(&a0, &f.sqr(&c)),
            b: f.mul(&b0, &f.mul(&f.sqr(&c), &c)),
        };
        if let Some(s) = isomorphism(f, m, &cand) {
            return Some((cand, s));
        }
    }
    None
}

/// A point in Jacobian coordinates; `z = 0` is the identity.
#[derive(Clone, Copy, Debug)]
struct Jac {
    x: Fe,
    y: Fe,
    z: Fe,
}

fn jac_double(f: &Field, m: &Model, p: &Jac) -> Jac {
    if f.is_zero(&p.z) || f.is_zero(&p.y) {
        return Jac {
            x: f.one(),
            y: f.one(),
            z: f.zero(),
        };
    }
    let xx = f.sqr(&p.x);
    let yy = f.sqr(&p.y);
    let yyyy = f.sqr(&yy);
    let zz = f.sqr(&p.z);
    let s = f.mul(&f.from_u64(4), &f.mul(&p.x, &yy));
    let mm = f.add(&f.mul(&f.from_u64(3), &xx), &f.mul(&m.a, &f.sqr(&zz)));
    let x3 = f.sub(&f.sqr(&mm), &f.add(&s, &s));
    let y3 = f.sub(&f.mul(&mm, &f.sub(&s, &x3)), &f.mul(&f.from_u64(8), &yyyy));
    let z3 = f.mul(&f.from_u64(2), &f.mul(&p.y, &p.z));
    Jac {
        x: x3,
        y: y3,
        z: z3,
    }
}

fn jac_add_affine(f: &Field, m: &Model, p: &Jac, qx: &Fe, qy: &Fe) -> Jac {
    if f.is_zero(&p.z) {
        return Jac {
            x: *qx,
            y: *qy,
            z: f.one(),
        };
    }
    let z1z1 = f.sqr(&p.z);
    let u2 = f.mul(qx, &z1z1);
    let s2 = f.mul(qy, &f.mul(&p.z, &z1z1));
    let h = f.sub(&u2, &p.x);
    let r = f.sub(&s2, &p.y);
    if f.is_zero(&h) {
        return if f.is_zero(&r) {
            jac_double(f, m, p)
        } else {
            Jac {
                x: f.one(),
                y: f.one(),
                z: f.zero(),
            }
        };
    }
    let hh = f.sqr(&h);
    let hhh = f.mul(&h, &hh);
    let v = f.mul(&p.x, &hh);
    let x3 = f.sub(&f.sub(&f.sqr(&r), &hhh), &f.add(&v, &v));
    let y3 = f.sub(&f.mul(&r, &f.sub(&v, &x3)), &f.mul(&p.y, &hhh));
    let z3 = f.mul(&p.z, &h);
    Jac {
        x: x3,
        y: y3,
        z: z3,
    }
}

/// `[k](x, y)` in affine form, `None` for the identity.
pub fn scalar_mul(f: &Field, m: &Model, x: &Fe, y: &Fe, k: &BigUint) -> Option<(Fe, Fe)> {
    let mut acc = Jac {
        x: f.one(),
        y: f.one(),
        z: f.zero(),
    };
    for i in (0..k.bits()).rev() {
        acc = jac_double(f, m, &acc);
        if k.bit(i) {
            acc = jac_add_affine(f, m, &acc, x, y);
        }
    }
    if f.is_zero(&acc.z) {
        return None;
    }
    let zi = f.inv(&acc.z)?;
    let zi2 = f.sqr(&zi);
    Some((f.mul(&acc.x, &zi2), f.mul(&acc.y, &f.mul(&zi2, &zi))))
}

/// The point with the least `x ≥ start` whose `y` is defined and nonzero,
/// with the smaller `y` as an integer.
pub fn least_point(f: &Field, m: &Model, start: u64) -> (Fe, Fe) {
    let mut x = start;
    loop {
        let fx = f.from_u64(x);
        let r = m.rhs(f, &fx);
        if f.legendre(&r) == 1 {
            let y = f.sqrt(&r).expect("residue");
            let ny = f.neg(&y);
            let y = if f.to_big(&y) <= f.to_big(&ny) { y } else { ny };
            return (fx, y);
        }
        x += 1;
    }
}

/// The generator a walk records for a curve it found: `[cofactor]` of the
/// least point that survives it.
pub const GENERATOR_RULE: &str =
    "icwalk-gen/v1: [cofactor]P for the point P with the least x >= 0 \
whose y is nonzero (the smaller y as an integer) and whose multiple is not \
the identity";

pub fn deterministic_generator(f: &Field, m: &Model, cofactor: &BigUint) -> (Fe, Fe) {
    let mut start = 0;
    loop {
        let (x, y) = least_point(f, m, start);
        if let Some(g) = scalar_mul(f, m, &x, &y, cofactor) {
            return g;
        }
        start = f.to_big(&x).try_into().unwrap_or(u64::MAX - 1) + 1;
    }
}

/// How the order of a model was established.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrderAudit {
    /// `[n]P = O` for a point `P ≠ O` and `n` prime above `(√p + 1)²`'s
    /// half width: `n | #E` and the Hasse interval holds one multiple.
    ProvedPrime,
    /// `[n]P = O` for every audited point; `n` composite, so this is
    /// evidence, not proof.
    Consistent,
    /// Some point has `[n]P ≠ O`: the order is not `n`.
    Refuted,
}

/// Audit `#E = n` with `points` deterministic points from `seed_x`.
pub fn audit_order(
    f: &Field,
    m: &Model,
    n: &BigUint,
    n_is_prime: bool,
    points: usize,
    seed_x: u64,
) -> OrderAudit {
    let mut start = seed_x;
    for _ in 0..points.max(1) {
        let (x, y) = least_point(f, m, start);
        if scalar_mul(f, m, &x, &y, n).is_some() {
            return OrderAudit::Refuted;
        }
        start = f.to_big(&x).try_into().unwrap_or(u64::MAX - 1) + 1;
    }
    // n prime and n > 4√p + 1 > the Hasse interval's width: n is the
    // only multiple of n in [p + 1 − 2√p, p + 1 + 2√p].
    let p = f.modulus();
    if n_is_prime && n * n > p * 16u8 {
        OrderAudit::ProvedPrime
    } else {
        OrderAudit::Consistent
    }
}

/// `#{x ∈ [0, 64) : x³ + ax + b is a nonzero square}` — the prefix-residue
/// statistic PR #1330's screening prototype sized.  Model dependent.
pub fn qr_prefix_64(f: &Field, m: &Model) -> u32 {
    (0..64u64)
        .filter(|x| f.legendre(&m.rhs(f, &f.from_u64(*x))) == 1)
        .count() as u32
}
