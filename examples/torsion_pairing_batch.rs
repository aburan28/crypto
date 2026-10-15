//! Native reusable-input Tate-pairing batch on frozen public targets.
//! The solver input deliberately excludes every answer scalar.
#![recursion_limit = "256"]

use std::cell::Cell;
use std::collections::BTreeMap;
use std::fs;
use std::time::Instant;

use crypto_lib::cryptanalysis::ext_curve::{
    jacobi, take_field_counters, ExtCurve, ExtDiagonalAutomorphism, ExtField, ExtPoint, Fp2,
};
use crypto_lib::cryptanalysis::glv_invariant_base::{
    rho_reference_folded, Endomorphism, Negation,
};
use crypto_lib::cryptanalysis::ic_boundary::{rho_reference_negation, CountedGroup, GroupOps};
use crypto_lib::cryptanalysis::residual_walk::{is_prime_u64, pow_mod};
use crypto_lib::hash::sha256::sha256;
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::{Rng, SeedableRng};
use serde::Serialize;
use serde_json::{json, Value};

type F2 = [u64; 2];
type Elt = Vec<F2>;
type SourcePoint = ExtPoint<F2>;

#[derive(Clone, Copy)]
struct Input {
    p: u64,
    nu: u64,
    r: u64,
    degree: usize,
    trace_fp: i64,
    g: SourcePoint,
    q: SourcePoint,
    seed: u64,
}

fn source_point(v: &Value) -> SourcePoint {
    let c = |name: &str| -> F2 {
        let xs = v[name].as_array().expect("coordinate pair");
        [xs[0].as_u64().unwrap(), xs[1].as_u64().unwrap()]
    };
    assert!(!v["infinity"].as_bool().unwrap());
    SourcePoint::affine(c("x"), c("y"))
}

fn read_input(path: &str, cell: usize) -> Input {
    let raw = fs::read_to_string(path).expect("public pairing input");
    let v: Value = serde_json::from_str(raw.lines().nth(cell).expect("cell")).unwrap();
    assert!(
        v.get("answer_scalar").is_none(),
        "solver input contains answer"
    );
    let p = v["p"].as_u64().unwrap();
    let r = v["r"].as_u64().unwrap();
    assert!((p == 61 && r == 19 && cell == 0)
        || (p == 271 && r == 61 && cell == 1)
        || (p == 401 && r == 67 && cell == 2));
    Input {
        p,
        nu: v["nu"].as_u64().unwrap(),
        r,
        degree: if cell == 0 { 9 } else if cell == 1 { 5 } else { 1 },
        trace_fp: v["trace_fp"].as_i64().unwrap(),
        g: source_point(&v["generator"]),
        q: source_point(&v["target"]),
        seed: if cell == 0 {
            0x5041_4952_0000_003d
        } else if cell == 1 {
            0x5041_4952_0000_010f
        } else {
            0x5041_4952_0000_0191
        },
    }
}

struct Field {
    f: Fp2,
    n: usize,
    a: F2,
    order: BigUint,
    muls: Cell<u64>,
    squares: Cell<u64>,
    inversions: Cell<u64>,
    budget_cpu_start_ns: u64,
}

impl Field {
    fn new(input: Input) -> Self {
        let f = Fp2::new(input.p).unwrap();
        assert_eq!(f.nu, input.nu);
        let q = f.order();
        let a = if input.degree == 1 {
            f.zero()
        } else {
            let prime_divisor = if input.degree == 9 { 3 } else { 5 };
            (1..q)
                .map(|i| f.from_coords(&[i % input.p, i / input.p]))
                .find(|&x| f.pow(x, (q - 1) / prime_divisor) != f.one())
                .expect("non-power coefficient")
        };
        Self {
            f,
            n: input.degree,
            a,
            order: BigUint::from(q).pow(input.degree as u32),
            muls: Cell::new(0),
            squares: Cell::new(0),
            inversions: Cell::new(0),
            budget_cpu_start_ns: thread_cpu_ns(),
        }
    }

    fn zero(&self) -> Elt {
        vec![self.f.zero(); self.n]
    }
    fn one(&self) -> Elt {
        self.embed(self.f.one())
    }
    fn embed(&self, x: F2) -> Elt {
        let mut z = self.zero();
        z[0] = x;
        z
    }
    fn add(&self, x: &Elt, y: &Elt) -> Elt {
        x.iter().zip(y).map(|(&a, &b)| self.f.add(a, b)).collect()
    }
    fn sub(&self, x: &Elt, y: &Elt) -> Elt {
        x.iter().zip(y).map(|(&a, &b)| self.f.sub(a, b)).collect()
    }
    fn mul(&self, x: &Elt, y: &Elt) -> Elt {
        let count = self.muls.get() + 1;
        assert!(
            count <= 10_000_000,
            "extension multiplication cap exhausted"
        );
        if count % 2048 == 0 {
            assert!(
                thread_cpu_ns() - self.budget_cpu_start_ns <= 900_000_000_000,
                "15-minute current-thread CPU cap exhausted"
            );
            assert!(
                peak_rss_bytes() <= 512 * 1024 * 1024,
                "512 MiB RSS cap exhausted"
            );
        }
        self.muls.set(count);
        let mut z = self.zero();
        for (i, &xi) in x.iter().enumerate() {
            if xi == self.f.zero() {
                continue;
            }
            for (j, &yj) in y.iter().enumerate() {
                if yj == self.f.zero() {
                    continue;
                }
                let mut product = self.f.mul(xi, yj);
                let mut d = i + j;
                if d >= self.n {
                    product = self.f.mul(product, self.a);
                    d -= self.n;
                }
                z[d] = self.f.add(z[d], product);
            }
        }
        z
    }
    fn square(&self, x: &Elt) -> Elt {
        self.squares.set(self.squares.get() + 1);
        self.mul(x, x)
    }
    fn pow(&self, x: &Elt, exp: &BigUint) -> Elt {
        let mut e = exp.clone();
        let mut base = x.clone();
        let mut out = self.one();
        while !e.is_zero() {
            if (&e & BigUint::one()) == BigUint::one() {
                out = self.mul(&out, &base);
            }
            e >>= 1;
            if !e.is_zero() {
                base = self.square(&base);
            }
        }
        out
    }
    fn inv(&self, x: &Elt) -> Option<Elt> {
        if *x == self.zero() {
            return None;
        }
        self.inversions.set(self.inversions.get() + 1);
        Some(self.pow(x, &(&self.order - BigUint::from(2u32))))
    }
    fn sqrt(&self, x: &Elt) -> Option<Elt> {
        if *x == self.zero() {
            return Some(self.zero());
        }
        let one = self.one();
        let qm1 = &self.order - BigUint::one();
        if self.pow(x, &(&qm1 >> 1)) != one {
            return None;
        }
        let mut odd = qm1;
        let mut s = 0u32;
        while (&odd & BigUint::one()).is_zero() {
            odd >>= 1;
            s += 1;
        }
        let nonsquare = self.embed(self.f.nonresidue());
        assert_ne!(
            self.pow(&nonsquare, &((&self.order - BigUint::one()) >> 1)),
            one
        );
        let mut c = self.pow(&nonsquare, &odd);
        let mut t = self.pow(x, &odd);
        let mut root = self.pow(x, &((&odd + BigUint::one()) >> 1));
        let mut m = s;
        while t != one {
            let mut i = 0u32;
            let mut tt = t.clone();
            while tt != one {
                tt = self.square(&tt);
                i += 1;
                if i == m {
                    return None;
                }
            }
            let mut b = c;
            for _ in 0..(m - i - 1) {
                b = self.square(&b);
            }
            let b2 = self.square(&b);
            root = self.mul(&root, &b);
            t = self.mul(&t, &b2);
            c = b2;
            m = i;
        }
        assert_eq!(self.square(&root), *x);
        Some(root)
    }
    fn counts(&self) -> Value {
        json!({"extension_muls":self.muls.get(),"extension_squares_in_muls":self.squares.get(),"extension_inversions":self.inversions.get()})
    }

    fn irreducibility_receipt(&self) -> Value {
        if self.n == 1 {
            assert_eq!(jacobi(self.f.nu, self.f.p), -1);
            return json!({"criterion":"quadratic_nonresidue",
                "p":self.f.p,"nu":self.f.nu,"degree_over_fp":2,
                "legendre_nu":-1});
        }
        // Rabin's criterion on the monic polynomial z^n-a, evaluated in
        // the quotient ring before any inverse in that ring is trusted.
        let theta = {
            let mut x = self.zero();
            x[1] = self.f.one();
            x
        };
        let q = BigUint::from(self.f.order());
        let top = self.pow(&theta, &q.pow(self.n as u32));
        assert_eq!(top, theta, "Frobenius closure failed");
        let prime = if self.n == 9 { 3 } else { 5 };
        let lower = self.pow(&theta, &q.pow((self.n / prime) as u32));
        let mut modulus = vec![self.f.zero(); self.n + 1];
        modulus[0] = self.f.neg(self.a);
        modulus[self.n] = self.f.one();
        let difference = self.sub(&lower, &theta);
        let gcd = poly_gcd(&self.f, modulus, difference);
        assert_eq!(gcd, vec![self.f.one()], "reducible tower polynomial");
        json!({"criterion":"rabin_frobenius_gcd","q":q.to_string(),"n":self.n,
            "prime_divisor":prime,"z_to_qn_equals_z":true,"gcd_coefficients":gcd})
    }
}

fn poly_trim(f: &Fp2, p: &mut Vec<F2>) {
    while p.last() == Some(&f.zero()) {
        p.pop();
    }
}

fn poly_rem(f: &Fp2, mut a: Vec<F2>, b: &[F2]) -> Vec<F2> {
    assert!(!b.is_empty());
    poly_trim(f, &mut a);
    let inv_lead = f
        .inv(*b.last().unwrap())
        .expect("polynomial leading coefficient");
    while a.len() >= b.len() {
        let offset = a.len() - b.len();
        let factor = f.mul(*a.last().unwrap(), inv_lead);
        for (j, &bj) in b.iter().enumerate() {
            a[offset + j] = f.sub(a[offset + j], f.mul(factor, bj));
        }
        poly_trim(f, &mut a);
    }
    a
}

fn poly_gcd(f: &Fp2, mut a: Vec<F2>, mut b: Vec<F2>) -> Vec<F2> {
    poly_trim(f, &mut a);
    poly_trim(f, &mut b);
    while !b.is_empty() {
        let rem = poly_rem(f, a, &b);
        a = b;
        b = rem;
    }
    let inv = f.inv(*a.last().expect("nonzero gcd")).unwrap();
    for c in &mut a {
        *c = f.mul(*c, inv);
    }
    a
}

#[derive(Clone, PartialEq, Eq)]
struct Point {
    x: Elt,
    y: Elt,
    infinity: bool,
}

impl Point {
    fn infinity(k: &Field) -> Self {
        Self {
            x: k.zero(),
            y: k.zero(),
            infinity: true,
        }
    }
    fn affine(x: Elt, y: Elt) -> Self {
        Self {
            x,
            y,
            infinity: false,
        }
    }
}

fn embed_point(k: &Field, p: SourcePoint) -> Point {
    if p.infinity {
        Point::infinity(k)
    } else {
        Point::affine(k.embed(p.x), k.embed(p.y))
    }
}

fn point_add(k: &Field, p: &Point, q: &Point) -> Point {
    if p.infinity {
        return q.clone();
    }
    if q.infinity {
        return p.clone();
    }
    let zero = k.zero();
    if p.x == q.x && k.add(&p.y, &q.y) == zero {
        return Point::infinity(k);
    }
    let lambda = if p.x == q.x {
        let three = k.embed(k.f.from_fp(3));
        let two = k.embed(k.f.from_fp(2));
        let num = k.mul(&three, &k.square(&p.x));
        let den = k.mul(&two, &p.y);
        k.mul(&num, &k.inv(&den).expect("nonzero tangent denominator"))
    } else {
        let num = k.sub(&q.y, &p.y);
        let den = k.sub(&q.x, &p.x);
        k.mul(&num, &k.inv(&den).expect("nonzero chord denominator"))
    };
    let x3 = k.sub(&k.sub(&k.square(&lambda), &p.x), &q.x);
    let y3 = k.sub(&k.mul(&lambda, &k.sub(&p.x, &x3)), &p.y);
    Point::affine(x3, y3)
}

fn point_mul(k: &Field, p: &Point, scalar: &BigUint) -> Point {
    let mut e = scalar.clone();
    let mut base = p.clone();
    let mut out = Point::infinity(k);
    while !e.is_zero() {
        if (&e & BigUint::one()) == BigUint::one() {
            out = point_add(k, &out, &base);
        }
        e >>= 1;
        if !e.is_zero() {
            base = point_add(k, &base, &base);
        }
    }
    out
}

fn field_curve_order(input: Input) -> BigUint {
    let p = BigInt::from(input.p);
    let q = &p * &p;
    let t = BigInt::from(input.trace_fp).pow(2) - BigInt::from(2) * &p;
    let mut prev = BigInt::from(2);
    let mut curr = t.clone();
    for _ in 2..=input.degree {
        let next = &t * &curr - &q * &prev;
        prev = curr;
        curr = next;
    }
    (q.pow(input.degree as u32) + BigInt::one() - curr)
        .to_biguint()
        .unwrap()
}

fn line_step(
    k: &Field,
    t: SourcePoint,
    u: SourcePoint,
    h: &Point,
) -> Option<(Elt, Elt, SourcePoint)> {
    let f = &k.f;
    assert!(!t.infinity && !u.infinity);
    if t.x == u.x && f.add(t.y, u.y) == f.zero() {
        let vertical = k.sub(&h.x, &k.embed(t.x));
        if vertical == k.zero() {
            return None;
        }
        return Some((vertical, k.one(), SourcePoint::infinity(f.zero())));
    }
    let slope = if t.x == u.x {
        let num = f.scale(f.sqr(t.x), 3);
        let den = f.scale(t.y, 2);
        f.mul(num, f.inv(den)?)
    } else {
        f.mul(f.sub(u.y, t.y), f.inv(f.sub(u.x, t.x))?)
    };
    let x3 = f.sub(f.sub(f.sqr(slope), t.x), u.x);
    let y3 = f.sub(f.mul(slope, f.sub(t.x, x3)), t.y);
    let line = k.sub(
        &k.sub(&h.y, &k.embed(t.y)),
        &k.mul(&k.embed(slope), &k.sub(&h.x, &k.embed(t.x))),
    );
    let vertical = k.sub(&h.x, &k.embed(x3));
    if line == k.zero() || vertical == k.zero() {
        return None;
    }
    Some((line, vertical, SourcePoint::affine(x3, y3)))
}

fn reduced_tate(k: &Field, p: SourcePoint, h: &Point, r: u64) -> Option<Elt> {
    let mut numerator = k.one();
    let mut denominator = k.one();
    let mut t = p;
    let bits = 64 - r.leading_zeros();
    for i in (0..bits - 1).rev() {
        let (line, vertical, next) = line_step(k, t, t, h)?;
        numerator = k.mul(&k.square(&numerator), &line);
        denominator = k.mul(&k.square(&denominator), &vertical);
        t = next;
        if (r >> i) & 1 == 1 {
            let (line, vertical, next) = line_step(k, t, p, h)?;
            numerator = k.mul(&numerator, &line);
            denominator = k.mul(&denominator, &vertical);
            t = next;
        }
    }
    assert!(t.infinity, "Miller terminal point is not identity");
    let raw = k.mul(&numerator, &k.inv(&denominator)?);
    let exp = (&k.order - BigUint::one()) / BigUint::from(r);
    let value = k.pow(&raw, &exp);
    if value == k.zero() {
        None
    } else {
        Some(value)
    }
}

#[derive(Clone, Copy)]
struct PhaseStart {
    cpu_ns: u64,
    k_muls: u64,
    k_squares: u64,
    k_invs: u64,
}

#[derive(Clone, Copy, Default, Serialize)]
struct PhaseCharge {
    cpu_ns: u64,
    fp_muls: u64,
    fp_invs: u64,
    extension_muls: u64,
    extension_squares_in_muls: u64,
    extension_invs: u64,
}

impl PhaseStart {
    fn now(k: &Field) -> Self {
        Self {
            cpu_ns: thread_cpu_ns(),
            k_muls: k.muls.get(),
            k_squares: k.squares.get(),
            k_invs: k.inversions.get(),
        }
    }
    fn finish(self, k: &Field) -> PhaseCharge {
        let fp = take_field_counters();
        PhaseCharge {
            cpu_ns: thread_cpu_ns() - self.cpu_ns,
            fp_muls: fp.0,
            fp_invs: fp.1,
            extension_muls: k.muls.get() - self.k_muls,
            extension_squares_in_muls: k.squares.get() - self.k_squares,
            extension_invs: k.inversions.get() - self.k_invs,
        }
    }
}

impl PhaseCharge {
    fn add(&mut self, other: Self) {
        self.cpu_ns += other.cpu_ns;
        self.fp_muls += other.fp_muls;
        self.fp_invs += other.fp_invs;
        self.extension_muls += other.extension_muls;
        self.extension_squares_in_muls += other.extension_squares_in_muls;
        self.extension_invs += other.extension_invs;
    }
    fn difference(self, earlier: Self) -> Self {
        Self {
            cpu_ns: self.cpu_ns - earlier.cpu_ns,
            fp_muls: self.fp_muls - earlier.fp_muls,
            fp_invs: self.fp_invs - earlier.fp_invs,
            extension_muls: self.extension_muls - earlier.extension_muls,
            extension_squares_in_muls: self.extension_squares_in_muls
                - earlier.extension_squares_in_muls,
            extension_invs: self.extension_invs - earlier.extension_invs,
        }
    }
}

#[derive(Default, Serialize)]
struct SelectionCharges {
    candidate_search: PhaseCharge,
    miller: PhaseCharge,
    final_exponentiation: PhaseCharge,
    pairing_checks: PhaseCharge,
}

impl SelectionCharges {
    fn total(&self) -> PhaseCharge {
        let mut out = self.candidate_search;
        out.add(self.miller);
        out.add(self.final_exponentiation);
        out.add(self.pairing_checks);
        out
    }
}

// The arithmetic is the same as reduced_tate; the two snapshots divide its
// Miller fraction (including denominator inversion) from final exponentiation.
fn reduced_tate_profiled(
    k: &Field,
    p: SourcePoint,
    h: &Point,
    r: u64,
    charges: &mut SelectionCharges,
) -> Option<Elt> {
    let miller_start = PhaseStart::now(k);
    let raw = (|| {
        let mut numerator = k.one();
        let mut denominator = k.one();
        let mut t = p;
        let bits = 64 - r.leading_zeros();
        for i in (0..bits - 1).rev() {
            let (line, vertical, next) = line_step(k, t, t, h)?;
            numerator = k.mul(&k.square(&numerator), &line);
            denominator = k.mul(&k.square(&denominator), &vertical);
            t = next;
            if (r >> i) & 1 == 1 {
                let (line, vertical, next) = line_step(k, t, p, h)?;
                numerator = k.mul(&numerator, &line);
                denominator = k.mul(&denominator, &vertical);
                t = next;
            }
        }
        assert!(t.infinity, "Miller terminal point is not identity");
        Some(k.mul(&numerator, &k.inv(&denominator)?))
    })();
    charges.miller.add(miller_start.finish(k));
    let raw = raw?;
    let final_start = PhaseStart::now(k);
    let exp = (&k.order - BigUint::one()) / BigUint::from(r);
    let value = k.pow(&raw, &exp);
    charges.final_exponentiation.add(final_start.finish(k));
    if value == k.zero() {
        None
    } else {
        Some(value)
    }
}

fn multiplicative_bsgs(k: &Field, alpha: &Elt, beta: &Elt, r: u64) -> Option<u64> {
    let mut width = 1u64;
    while width * width < r {
        width += 1;
    }
    let mut baby = BTreeMap::new();
    let mut value = k.one();
    for j in 0..width {
        baby.entry(value.clone()).or_insert(j);
        value = k.mul(&value, alpha);
    }
    let giant_factor = k.pow(alpha, &BigUint::from(r - width));
    let mut gamma = beta.clone();
    for i in 0..=width {
        if let Some(&j) = baby.get(&gamma) {
            let candidate = (i * width + j) % r;
            if k.pow(alpha, &BigUint::from(candidate)) == *beta {
                return Some(candidate);
            }
        }
        gamma = k.mul(&gamma, &giant_factor);
    }
    None
}

struct ReusableBsgs {
    baby: BTreeMap<Elt, u64>,
    width: u64,
    giant_factor: Elt,
    alpha: Elt,
}

impl ReusableBsgs {
    fn new(k: &Field, alpha: &Elt, r: u64) -> Self {
        let mut width = 1u64;
        while width * width < r {
            width += 1;
        }
        let mut baby = BTreeMap::new();
        let mut value = k.one();
        for j in 0..width {
            baby.entry(value.clone()).or_insert(j);
            value = k.mul(&value, alpha);
        }
        let giant_factor = k.pow(alpha, &BigUint::from(r - width));
        Self {
            baby,
            width,
            giant_factor,
            alpha: alpha.clone(),
        }
    }
    fn solve(&self, k: &Field, beta: &Elt, r: u64) -> Option<u64> {
        let mut gamma = beta.clone();
        for i in 0..=self.width {
            if let Some(&j) = self.baby.get(&gamma) {
                let candidate = (i * self.width + j) % r;
                if k.pow(&self.alpha, &BigUint::from(candidate)) == *beta {
                    return Some(candidate);
                }
            }
            gamma = k.mul(&gamma, &self.giant_factor);
        }
        None
    }
}

fn select_point(
    k: &Field,
    input: Input,
) -> (Point, Value, Elt, Elt, SelectionCharges, PhaseCharge) {
    let mut charges = SelectionCharges::default();
    let mut search_start = PhaseStart::now(k);
    let order = field_curve_order(input);
    let r = BigUint::from(input.r);
    let mut cofactor = order.clone();
    let mut valuation = 0u32;
    while (&cofactor % &r).is_zero() {
        cofactor /= &r;
        valuation += 1;
    }
    assert!(valuation >= 2, "full rational r-torsion expected");
    let mut rng = StdRng::seed_from_u64(input.seed);
    let g = embed_point(k, input.g);
    let mut order_candidates = 0u32;
    let mut pairing_poles = 0u32;
    let mut degenerate = 0u32;
    for attempt in 1..=4096u32 {
        let x: Elt = (0..k.n)
            .map(|_| [rng.gen_range(0..input.p), rng.gen_range(0..input.p)])
            .collect();
        let rhs = k.add(&k.mul(&k.square(&x), &x), &k.one());
        let Some(y) = k.sqrt(&rhs) else {
            continue;
        };
        let raw = Point::affine(x, y);
        let mut h = point_mul(k, &raw, &cofactor);
        for _ in 0..valuation {
            if h.infinity || point_mul(k, &h, &r).infinity {
                break;
            }
            h = point_mul(k, &h, &r);
        }
        if h.infinity || !point_mul(k, &h, &r).infinity {
            continue;
        }
        order_candidates += 1;
        if order_candidates > 512 {
            panic!("512 acceptable-order pairing inputs exhausted");
        }
        if h.x.iter().skip(1).all(|&v| v == k.f.zero())
            && h.y.iter().skip(1).all(|&v| v == k.f.zero())
        {
            let mut multiple = Point::infinity(k);
            let mut dependent = false;
            for _ in 0..input.r {
                if multiple == h {
                    dependent = true;
                    break;
                }
                multiple = point_add(k, &multiple, &g);
            }
            if dependent {
                continue;
            }
        }
        charges.candidate_search.add(search_start.finish(k));
        let Some(alpha) = reduced_tate_profiled(k, input.g, &h, input.r, &mut charges) else {
            pairing_poles += 1;
            search_start = PhaseStart::now(k);
            continue;
        };
        let check_start = PhaseStart::now(k);
        let valid_alpha = alpha != k.one() && k.pow(&alpha, &r) == k.one();
        charges.pairing_checks.add(check_start.finish(k));
        if !valid_alpha {
            degenerate += 1;
            search_start = PhaseStart::now(k);
            continue;
        }
        let before_beta = charges.total();
        let Some(beta) = reduced_tate_profiled(k, input.q, &h, input.r, &mut charges) else {
            pairing_poles += 1;
            search_start = PhaseStart::now(k);
            continue;
        };
        let check_start = PhaseStart::now(k);
        let valid_beta = k.pow(&beta, &r) == k.one();
        charges.pairing_checks.add(check_start.finish(k));
        if !valid_beta {
            degenerate += 1;
            search_start = PhaseStart::now(k);
            continue;
        }
        let receipt = json!({"attempt":attempt,"curve_order":order.to_string(),"r_valuation":valuation,
            "cofactor":cofactor.to_string(),"order_candidates":order_candidates,
            "pairing_poles":pairing_poles,"degenerate":degenerate,
            "h":{"x":h.x,"y":h.y},"alpha":alpha,"beta":beta});
        let first_beta_charge = charges.total().difference(before_beta);
        return (h, receipt, alpha, beta, charges, first_beta_charge);
    }
    panic!("point search reached 4096 candidate x coordinates");
}

fn thread_cpu_ns() -> u64 {
    let mut ts = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    let rc = unsafe { libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut ts) };
    assert_eq!(rc, 0);
    ts.tv_sec as u64 * 1_000_000_000 + ts.tv_nsec as u64
}

fn peak_rss_raw() -> (u64, &'static str) {
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    assert_eq!(unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) }, 0);
    #[cfg(target_os = "macos")]
    {
        (usage.ru_maxrss as u64, "bytes")
    }
    #[cfg(not(target_os = "macos"))]
    {
        (usage.ru_maxrss as u64, "kibibytes")
    }
}

fn peak_rss_bytes() -> u64 {
    let (value, unit) = peak_rss_raw();
    if unit == "bytes" {
        value
    } else {
        value * 1024
    }
}

fn source_curve(input: Input) -> ExtCurve<Fp2> {
    let f = Fp2::new(input.p).unwrap();
    assert_eq!(f.nu, input.nu);
    let curve = ExtCurve {
        f,
        a: f.zero(),
        b: f.one(),
    };
    assert!(curve.is_on_curve(input.g) && curve.is_on_curve(input.q));
    curve
}

fn rho_step_cap(r: u64) -> u64 {
    let root = (r as f64).sqrt().ceil() as u64;
    64.max(16 * root)
}

fn six_class_endomorphism(
    curve: &ExtCurve<Fp2>,
    g: SourcePoint,
    r: u64,
) -> (ExtDiagonalAutomorphism<F2>, GroupOps) {
    let p = curve.f.p;
    let z = (2..p)
        .map(|a| {
            let mut x = 1u64;
            for _ in 0..(p - 1) / 3 {
                x = ((x as u128 * a as u128) % p as u128) as u64;
            }
            x
        })
        .find(|&x| x != 1)
        .expect("nontrivial cube root of unity");
    let image = SourcePoint::affine(curve.f.mul(curve.f.from_fp(z), g.x), g.y);
    let mut ops = GroupOps::default();
    let lambda = (1..r)
        .filter(|&a| (a * a + a + 1) % r == 0)
        .find(|&a| curve.mul(&mut ops, g, a) == image)
        .expect("j=0 subgroup eigenvalue");
    (
        ExtDiagonalAutomorphism {
            cx: curve.f.from_fp(z),
            cy: curve.f.one(),
            eigenvalue: lambda,
            order: 3,
            label: "j0-zeta",
        },
        ops,
    )
}

fn calibrate_source_add(curve: &ExtCurve<Fp2>, input: Input) -> f64 {
    let points = [input.g, input.q, curve.neg(input.g), curve.neg(input.q)];
    let mut ops = GroupOps::default();
    let mut acc = input.g;
    let start = thread_cpu_ns();
    for i in 0..100_000usize {
        acc = curve.add(&mut ops, acc, points[(i * 7 + 1) % points.len()]);
        if acc.infinity {
            acc = points[(i + 3) % points.len()];
        }
    }
    std::hint::black_box(acc);
    (thread_cpu_ns() - start) as f64 / 100_000.0
}

fn run_prepare(input: Input, cell: usize) {
    take_field_counters();
    let wall = Instant::now();
    let cpu = thread_cpu_ns();
    let k = Field::new(input);
    let cert = k.irreducibility_receipt();
    let setup_fp = take_field_counters();
    let (h, selection, alpha, beta, charges, _) = select_point(&k, input);
    let tail_fp = take_field_counters();
    let charged = charges.total();
    let (rss, rss_unit) = peak_rss_raw();
    println!(
        "{}",
        json!({"schema":"torsion_pairing_prepare/v2","cell":cell,"p":input.p,
        "r":input.r,"tower_degree":k.n,"tower_a":k.a,"irreducibility":cert,
        "selection":selection,"selected_h":{"x":h.x,"y":h.y},"alpha":alpha,"beta":beta,
        "field_counts":k.counts(),"fp_counts":{"muls":setup_fp.0+charged.fp_muls+tail_fp.0,
            "invs":setup_fp.1+charged.fp_invs+tail_fp.1},
        "selection_phase_charges":charges,"selection_unclassified_fp":{"muls":tail_fp.0,"invs":tail_fp.1},
        "cpu_ns":thread_cpu_ns()-cpu,"wall_ns":wall.elapsed().as_nanos() as u64,
        "peak_rss_raw":rss,"peak_rss_unit":rss_unit})
    );
}

fn run_pairing(input: Input, cell: usize, frozen_path: &str) {
    let frozen: Value = serde_json::from_str(&fs::read_to_string(frozen_path).unwrap()).unwrap();
    assert_eq!(frozen["cell"].as_u64(), Some(cell as u64));
    take_field_counters();
    let wall = Instant::now();
    let cpu = thread_cpu_ns();
    let k = Field::new(input);
    let cert = k.irreducibility_receipt();
    let setup_cpu = thread_cpu_ns() - cpu;
    let setup_wall = wall.elapsed().as_nanos() as u64;
    let setup_fp = take_field_counters();
    let setup_k = k.counts();
    assert_eq!(frozen["tower_a"], json!(k.a));
    assert_eq!(frozen["irreducibility"], cert);
    let select_wall = Instant::now();
    let select_cpu = thread_cpu_ns();
    let (h, selection, alpha, beta, charges, _) = select_point(&k, input);
    let select_tail_fp = take_field_counters();
    let charged = charges.total();
    let select_fp = (
        charged.fp_muls + select_tail_fp.0,
        charged.fp_invs + select_tail_fp.1,
    );
    let select_k = k.counts();
    let select_cpu_ns = thread_cpu_ns() - select_cpu;
    let select_wall_ns = select_wall.elapsed().as_nanos() as u64;
    let select_unclassified_cpu_ns = select_cpu_ns
        .checked_sub(charged.cpu_ns)
        .expect("phase CPU exceeds selection interval");
    assert_eq!(frozen["selected_h"], json!({"x":h.x,"y":h.y}));
    assert_eq!(frozen["selection"], selection);
    assert_eq!(frozen["alpha"], json!(alpha));
    assert_eq!(frozen["beta"], json!(beta));
    let dlog_wall = Instant::now();
    let dlog_cpu = thread_cpu_ns();
    let recovered = multiplicative_bsgs(&k, &alpha, &beta, input.r);
    let scalar = recovered.expect("multiplicative DLP cap exhausted");
    assert_eq!(k.pow(&alpha, &BigUint::from(scalar)), beta);
    let dlog_fp = take_field_counters();
    let dlog_k = k.counts();
    let dlog_cpu_ns = thread_cpu_ns() - dlog_cpu;
    let dlog_wall_ns = dlog_wall.elapsed().as_nanos() as u64;
    let cold_cpu_ns = thread_cpu_ns() - cpu;
    let cold_wall_ns = wall.elapsed().as_nanos() as u64;
    let (rss, rss_unit) = peak_rss_raw();
    let curve = source_curve(input);
    let ns_per_source_add = calibrate_source_add(&curve, input);
    println!(
        "{}",
        json!({"schema":"torsion_pairing_run/v2","cell":cell,"p":input.p,
        "r":input.r,"recovered":scalar,"selected_h":{"x":h.x,"y":h.y},
        "setup":{"cpu_ns":setup_cpu,"wall_ns":setup_wall,"fp_muls":setup_fp.0,"fp_invs":setup_fp.1,"k_cumulative":setup_k},
        "selection_and_pairings":{"cpu_ns":select_cpu_ns,"wall_ns":select_wall_ns,
            "fp_muls":select_fp.0,"fp_invs":select_fp.1,"k_cumulative":select_k,
            "phase_charges":charges,"unclassified_cpu_ns":select_unclassified_cpu_ns,
            "unclassified_fp":{"muls":select_tail_fp.0,"invs":select_tail_fp.1}},
        "multiplicative_dlp":{"cpu_ns":dlog_cpu_ns,"wall_ns":dlog_wall_ns,"fp_muls":dlog_fp.0,"fp_invs":dlog_fp.1,"k_cumulative":dlog_k},
        "cold_cpu_ns":cold_cpu_ns,"cold_wall_ns":cold_wall_ns,
        "ns_per_source_add_calibration":ns_per_source_add,
        "peak_rss_raw":rss,"peak_rss_unit":rss_unit})
    );
}

fn run_rho(input: Input, cell: usize, round: u64) {
    assert!(round < 3);
    let seed = 0x5041_4952_0000_1000u64 + 256 * cell as u64 + round;
    take_field_counters();
    let wall = Instant::now();
    let cpu = thread_cpu_ns();
    let curve = source_curve(input);
    let result = rho_reference_negation(&curve, input.g, input.q, input.r, seed, 1_000_000);
    let fp = take_field_counters();
    let cpu_ns = thread_cpu_ns() - cpu;
    let wall_ns = wall.elapsed().as_nanos() as u64;
    let (rss, rss_unit) = peak_rss_raw();
    let ns_per_source_add = calibrate_source_add(&curve, input);
    println!(
        "{}",
        json!({"schema":"torsion_pairing_rho/v1","cell":cell,"p":input.p,
        "r":input.r,"round":round,"seed":seed,"rho":result,"cpu_ns":cpu_ns,
        "wall_ns":wall_ns,"fp_muls":fp.0,"fp_invs":fp.1,
        "ns_per_source_add_calibration":ns_per_source_add,
        "peak_rss_raw":rss,"peak_rss_unit":rss_unit})
    );
}

fn freeze_401(previous_inputs: &str, extended_inputs: &str, fixture_path: &str) {
    let p = 401u64;
    let r = 67u64;
    assert!(is_prime_u64(p) && is_prime_u64(r));
    let f = Fp2::new(p).unwrap();
    let curve = ExtCurve { f, a: f.zero(), b: f.one() };
    let mut n1 = 1u64;
    for x in 0..p {
        let rhs = (pow_mod(x,3,p) + 1) % p;
        n1 += if rhs == 0 { 1 } else if jacobi(rhs,p) == 1 { 2 } else { 0 };
    }
    assert_eq!(n1,402);
    let trace = p as i64 + 1 - n1 as i64;
    assert_eq!(trace,0);
    let n2 = (p+1)*(p+1) - (trace*trace) as u64;
    assert_eq!(n2,161604);
    assert_eq!(n2 % (r*r),0);
    let mut ops = GroupOps::default();
    let mut selected = None;
    for x in 0..p {
        let mut lifts = curve.lift_x([x,0]);
        lifts.retain(|pt| pt.y[1] == 0);
        lifts.sort_by_key(|pt| pt.y[0]);
        for point in lifts {
            let g = curve.mul(&mut ops, point,6);
            if !g.infinity {
                selected = Some((point,g));
                break;
            }
        }
        if selected.is_some() { break; }
    }
    let (selected_point,g) = selected.expect("order-67 source point");
    assert!(!g.infinity && curve.mul(&mut ops,g,r).infinity);
    assert_ne!(g.x[0],0);
    assert_eq!(g.x[1],0);
    let target_seed = "pairing-batch/2026-10-15/p=401/r=67";
    let digest = sha256(target_seed.as_bytes());
    let scalar = 1 + u64::from_be_bytes(digest[0..8].try_into().unwrap()) % (r-1);
    let q = curve.mul(&mut ops,g,scalar);
    assert!(!q.infinity && curve.mul(&mut ops,q,r).infinity);
    let public = json!({"p":p,"nu":f.nu,"r":r,"trace_fp":trace,
        "generator":g,"target":q});
    let previous = fs::read_to_string(previous_inputs).unwrap();
    assert_eq!(previous.lines().count(),2);
    let extended = format!("{}\n{}\n",previous.trim_end(),public);
    fs::write(extended_inputs,extended).unwrap();
    fs::write(fixture_path,serde_json::to_vec_pretty(&json!({
        "schema":"torsion_pairing_401_fixture/v1","p":p,"nu":f.nu,
        "r":r,"n1":n1,"n2":n2,"trace_fp":trace,"selected_point":selected_point,
        "generator":g,"target":q,"answer_scalar":scalar,
        "target_seed":target_seed,"source_group_ops":ops})).unwrap()).unwrap();
    println!("{}",json!({"schema":"torsion_pairing_401_freeze/v1",
        "p":p,"r":r,"n1":n1,"n2":n2,"nu":f.nu,"generator":g,
        "target":q,"extended_inputs":extended_inputs,"fixture":fixture_path}));
}

fn freeze_batch(
    input: Input,
    cell: usize,
    fixture_path: &str,
    public_path: &str,
    answer_path: &str,
) {
    let fixture = fs::read_to_string(fixture_path).unwrap();
    let row: Value = if cell == 2 {
        serde_json::from_str(&fixture).unwrap()
    } else {
        serde_json::from_str(fixture.lines().nth(cell).unwrap()).unwrap()
    };
    let first = row["answer_scalar"]
        .as_u64()
        .expect("post-solve fixture scalar");
    let curve = source_curve(input);
    let mut ops = GroupOps::default();
    assert_eq!(curve.mul(&mut ops, input.g, first), input.q);
    let mut remaining: Vec<u64> = (1..input.r).filter(|&x| x != first).collect();
    let seed = 0x4241_5443_4800_0000u64 ^ input.p;
    remaining.shuffle(&mut StdRng::seed_from_u64(seed));
    let mut scalars = vec![first];
    scalars.extend(remaining);
    let targets: Vec<_> = scalars
        .iter()
        .map(|&x| curve.mul(&mut ops, input.g, x))
        .collect();
    assert_eq!(targets.len(), (input.r - 1) as usize);
    assert_eq!(targets[0], input.q);
    let mut unique = std::collections::BTreeSet::new();
    for q in &targets {
        assert!(!q.infinity && unique.insert(curve.key(q)));
    }
    fs::write(
        public_path,
        serde_json::to_vec_pretty(&json!({
            "schema":"torsion_pairing_batch_public/v1","cell":cell,"p":input.p,
            "r":input.r,"seed":seed,"targets":targets
        }))
        .unwrap(),
    )
    .unwrap();
    fs::write(
        answer_path,
        serde_json::to_vec_pretty(&json!({
            "schema":"torsion_pairing_batch_answers/v1","cell":cell,"p":input.p,
            "r":input.r,"scalars":scalars
        }))
        .unwrap(),
    )
    .unwrap();
    println!(
        "{}",
        json!({"schema":"torsion_pairing_batch_freeze/v1","cell":cell,
        "targets":targets.len(),"seed":seed,"group_ops":ops,
        "public_path":public_path,"answer_path":answer_path})
    );
}

fn read_batch(path: &str, input: Input, cell: usize) -> Vec<SourcePoint> {
    let raw = fs::read_to_string(path).unwrap();
    let v: Value = serde_json::from_str(&raw).unwrap();
    assert_eq!(v["schema"], "torsion_pairing_batch_public/v1");
    assert_eq!(v["cell"].as_u64(), Some(cell as u64));
    assert_eq!(v["p"].as_u64(), Some(input.p));
    assert_eq!(v["r"].as_u64(), Some(input.r));
    assert!(v.get("scalars").is_none() && v.get("answer_scalar").is_none());
    let targets: Vec<_> = v["targets"]
        .as_array()
        .unwrap()
        .iter()
        .map(source_point)
        .collect();
    assert_eq!(targets.len(), (input.r - 1) as usize);
    assert_eq!(targets[0], input.q);
    let curve = source_curve(input);
    let mut ops = GroupOps::default();
    let mut seen = std::collections::BTreeSet::new();
    for &q in &targets {
        assert!(curve.is_on_curve(q) && !q.infinity && curve.mul(&mut ops, q, input.r).infinity);
        assert!(seen.insert(curve.key(&q)));
    }
    targets
}

fn run_batch(input: Input, cell: usize, round: u64, batch_path: &str, frozen_path: &str) {
    assert!(round < 3);
    let global_wall = Instant::now();
    let global_cpu = thread_cpu_ns();
    let targets = read_batch(batch_path, input, cell);
    let frozen: Value = serde_json::from_str(&fs::read_to_string(frozen_path).unwrap()).unwrap();
    let curve = source_curve(input);
    take_field_counters();
    let whole_wall = Instant::now();
    let whole_cpu = thread_cpu_ns();
    let k = Field::new(input);
    let cert = k.irreducibility_receipt();
    assert_eq!(frozen["tower_a"], json!(k.a));
    assert_eq!(frozen["irreducibility"], cert);
    let field_cpu = thread_cpu_ns() - whole_cpu;
    let field_wall = whole_wall.elapsed().as_nanos() as u64;
    let field_fp = take_field_counters();
    let select_wall = Instant::now();
    let select_cpu = thread_cpu_ns();
    let (h, selection, alpha, first_beta, selection_charges, first_beta_charge) =
        select_point(&k, input);
    let selection_cpu = thread_cpu_ns() - select_cpu;
    let selection_wall = select_wall.elapsed().as_nanos() as u64;
    let selection_tail_fp = take_field_counters();
    assert_eq!(frozen["selected_h"], json!({"x":h.x,"y":h.y}));
    assert_eq!(frozen["selection"], selection);
    assert_eq!(frozen["alpha"], json!(alpha));
    assert_eq!(frozen["beta"], json!(first_beta));
    let table_wall = Instant::now();
    let table_cpu = thread_cpu_ns();
    let table_start = PhaseStart::now(&k);
    let table = ReusableBsgs::new(&k, &alpha, input.r);
    let table_charge = table_start.finish(&k);
    let table_cpu_ns = thread_cpu_ns() - table_cpu;
    let table_wall_ns = table_wall.elapsed().as_nanos() as u64;
    let setup_cpu = field_cpu + selection_cpu - first_beta_charge.cpu_ns + table_cpu_ns;
    let setup_wall =
        field_wall + selection_wall.saturating_sub(first_beta_charge.cpu_ns) + table_wall_ns;
    let selection_total = selection_charges.total();
    let setup_charge = selection_total.difference(first_beta_charge);
    let (rss, unit) = peak_rss_raw();
    assert!(peak_rss_bytes() <= 512 * 1024 * 1024);
    println!(
        "{}",
        json!({"schema":"torsion_pairing_batch_setup/v1","cell":cell,
        "round":round,"p":input.p,"r":input.r,"field_degree_over_fp":2*k.n,
        "tower_a":k.a,"irreducibility":cert,"selection":selection,
        "setup_cpu_ns":setup_cpu,"setup_wall_ns_estimate":setup_wall,
        "field_cpu_ns":field_cpu,"field_wall_ns":field_wall,"field_fp_muls":field_fp.0,
        "field_fp_invs":field_fp.1,"selection_cpu_ns":selection_cpu,
        "selection_wall_ns":selection_wall,"selection_charges":selection_charges,
        "selection_tail_fp":{"muls":selection_tail_fp.0,"invs":selection_tail_fp.1},
        "first_beta_charge":first_beta_charge,"setup_selection_charge":setup_charge,
        "reusable_bsgs_table":{"entries":table.baby.len(),"charge":table_charge,
            "cpu_ns":table_cpu_ns,"wall_ns":table_wall_ns},
        "peak_rss_raw":rss,"peak_rss_unit":unit})
    );
    let mut six_setup_cpu_ns = 0u64;
    let mut six_setup_wall_ns = 0u64;
    let six_aut = if cell == 2 {
        assert_eq!(input.p % 3,2);
        assert_ne!(input.g.x[0],0);
        assert_eq!(input.g.x[1],0);
        println!("{}",json!({"schema":"torsion_pairing_batch_six_setup/v4",
            "cell":cell,"round":round,"applicable":false,
            "reason":"order-three coordinate automorphism has zeta outside Fp; nonzero Fp x-coordinate leaves the Fp-rational source subgroup"}));
        None
    } else {
        let six_setup_wall = Instant::now();
        let six_setup_cpu = thread_cpu_ns();
        take_field_counters();
        let (aut, eigen_ops) = six_class_endomorphism(&curve, input.g, input.r);
        let six_setup_fp = take_field_counters();
        six_setup_cpu_ns = thread_cpu_ns() - six_setup_cpu;
        six_setup_wall_ns = six_setup_wall.elapsed().as_nanos() as u64;
        println!("{}",json!({"schema":"torsion_pairing_batch_six_setup/v3",
            "cell":cell,"round":round,"applicable":true,
            "eigenvalue":aut.eigenvalue,"eigen_ops":eigen_ops,
            "cpu_ns":six_setup_cpu_ns,"wall_ns":six_setup_wall_ns,
            "fp_muls":six_setup_fp.0,"fp_invs":six_setup_fp.1}));
        Some(aut)
    };
    let mut cumulative_pairing_cpu = setup_cpu;
    let mut cumulative_rho_cpu = 0u64;
    let mut cumulative_six_cpu = six_setup_cpu_ns;
    let mut cumulative_pairing_wall = setup_wall;
    let mut cumulative_rho_wall = 0u64;
    let mut cumulative_six_wall = six_setup_wall_ns;
    let rho_cap = rho_step_cap(input.r);
    for (index, &q) in targets.iter().enumerate() {
        let beta_start_wall = Instant::now();
        let beta_start_cpu = thread_cpu_ns();
        let (beta, beta_charge, beta_cpu, beta_wall) = if index == 0 {
            (
                first_beta.clone(),
                first_beta_charge,
                first_beta_charge.cpu_ns,
                first_beta_charge.cpu_ns,
            )
        } else {
            let mut charge = SelectionCharges::default();
            let beta =
                reduced_tate_profiled(&k, q, &h, input.r, &mut charge).expect("batch pairing pole");
            let check_start = PhaseStart::now(&k);
            assert_eq!(k.pow(&beta, &BigUint::from(input.r)), k.one());
            charge.pairing_checks.add(check_start.finish(&k));
            (
                beta,
                charge.total(),
                thread_cpu_ns() - beta_start_cpu,
                beta_start_wall.elapsed().as_nanos() as u64,
            )
        };
        let dlog_wall = Instant::now();
        let dlog_cpu = thread_cpu_ns();
        let dlog_start = PhaseStart::now(&k);
        let scalar = table
            .solve(&k, &beta, input.r)
            .expect("batch multiplicative DLP exhausted");
        assert_eq!(k.pow(&alpha, &BigUint::from(scalar)), beta);
        let dlog_charge = dlog_start.finish(&k);
        let dlog_cpu_ns = thread_cpu_ns() - dlog_cpu;
        let dlog_wall_ns = dlog_wall.elapsed().as_nanos() as u64;
        let verify_wall = Instant::now();
        let verify_cpu = thread_cpu_ns();
        let mut verify_ops = GroupOps::default();
        assert_eq!(curve.mul(&mut verify_ops, input.g, scalar), q);
        let verify_cpu_ns = thread_cpu_ns() - verify_cpu;
        let verify_wall_ns = verify_wall.elapsed().as_nanos() as u64;
        let target_cpu = beta_cpu + dlog_cpu_ns + verify_cpu_ns;
        let target_wall = beta_wall + dlog_wall_ns + verify_wall_ns;
        cumulative_pairing_cpu += target_cpu;
        cumulative_pairing_wall += target_wall;
        let seed_base = 0x4241_5443_4800_1000u64 + 256 * cell as u64 + 64 * round + index as u64;
        let mut rho_attempts = Vec::new();
        let mut rho_cpu_ns = 0u64;
        let mut rho_wall_ns = 0u64;
        let mut rho_fp_muls = 0u64;
        let mut rho_fp_invs = 0u64;
        let mut final_seed = seed_base;
        let mut final_rho = None;
        for retry in 0..32u64 {
            let seed = seed_base.wrapping_add(retry * 0x01_00_00_00);
            let rho_wall = Instant::now();
            let rho_cpu = thread_cpu_ns();
            take_field_counters();
            let rho = rho_reference_negation(&curve, input.g, q, input.r, seed, rho_cap);
            let rho_fp = take_field_counters();
            let attempt_cpu = thread_cpu_ns() - rho_cpu;
            let attempt_wall = rho_wall.elapsed().as_nanos() as u64;
            rho_cpu_ns += attempt_cpu;
            rho_wall_ns += attempt_wall;
            rho_fp_muls += rho_fp.0;
            rho_fp_invs += rho_fp.1;
            rho_attempts.push(json!({"retry":retry,"seed":seed,"rho":rho,
                "cpu_ns":attempt_cpu,"wall_ns":attempt_wall,"fp_muls":rho_fp.0,
                "fp_invs":rho_fp.1}));
            if rho.verified && rho.recovered == Some(scalar) {
                final_seed = seed;
                final_rho = Some(rho);
                break;
            }
        }
        if final_rho.is_none() {
            println!(
                "{}",
                json!({"schema":"torsion_pairing_batch_rho_exhausted/v3",
                "cell":cell,"round":round,"index":index,"target":q,"recovered":scalar,
                "rho_attempts":rho_attempts,"rho_cpu_ns":rho_cpu_ns})
            );
            panic!("32 charged rho attempts exhausted");
        }
        let rho = final_rho.unwrap();
        cumulative_rho_cpu += rho_cpu_ns;
        cumulative_rho_wall += rho_wall_ns;
        let mut six_attempts = Vec::new();
        let mut six_cpu_ns = 0u64;
        let mut six_wall_ns = 0u64;
        let mut six_fp_muls = 0u64;
        let mut six_fp_invs = 0u64;
        let mut final_six = None;
        if let Some(aut) = &six_aut {
            let neg = Negation { r: input.r };
            let gens: [&dyn Endomorphism<ExtCurve<Fp2>>; 2] = [&neg, aut];
            for retry in 0..32u64 {
                let seed = seed_base.wrapping_add(retry * 0x01_00_00_00);
                let six_wall = Instant::now();
                let six_cpu = thread_cpu_ns();
                take_field_counters();
                let six = rho_reference_folded(&curve, input.g, q, input.r, seed, rho_cap, &gens)
                    .expect("valid six-class subgroup action");
                let six_fp = take_field_counters();
                let attempt_cpu = thread_cpu_ns() - six_cpu;
                let attempt_wall = six_wall.elapsed().as_nanos() as u64;
                six_cpu_ns += attempt_cpu;
                six_wall_ns += attempt_wall;
                six_fp_muls += six_fp.0;
                six_fp_invs += six_fp.1;
                six_attempts.push(json!({"retry":retry,"seed":seed,"rho":six,
                    "cpu_ns":attempt_cpu,"wall_ns":attempt_wall,
                    "fp_muls":six_fp.0,"fp_invs":six_fp.1}));
                if six.verified && six.recovered == Some(scalar) {
                    final_six = Some(six);
                    break;
                }
            }
            if final_six.is_none() {
                println!("{}",json!({"schema":"torsion_pairing_batch_six_exhausted/v3",
                    "cell":cell,"round":round,"index":index,"target":q,
                    "recovered":scalar,"six_attempts":six_attempts,
                    "six_cpu_ns":six_cpu_ns}));
                panic!("32 charged six-class rho attempts exhausted");
            }
        }
        let six = final_six;
        cumulative_six_cpu += six_cpu_ns;
        cumulative_six_wall += six_wall_ns;
        let all_in_pairing_cpu = thread_cpu_ns() - global_cpu - cumulative_rho_cpu - cumulative_six_cpu;
        let all_in_pairing_wall = global_wall.elapsed().as_nanos() as u64
            - cumulative_rho_wall - cumulative_six_wall;
        let (rss, unit) = peak_rss_raw();
        assert!(peak_rss_bytes() <= 512 * 1024 * 1024);
        println!(
            "{}",
            json!({"schema":if cell==2 {"torsion_pairing_batch_target/v4"} else {"torsion_pairing_batch_target/v3"},
            "cell":cell,"round":round,"index":index,"n":index+1,"p":input.p,"r":input.r,
            "target":q,"recovered":scalar,"beta_charge":beta_charge,
            "dlog_charge":dlog_charge,"source_verify_ops":verify_ops,
            "beta_cpu_ns":beta_cpu,"dlog_cpu_ns":dlog_cpu_ns,"verify_cpu_ns":verify_cpu_ns,
            "target_cpu_ns":target_cpu,"target_wall_ns":target_wall,
            "cumulative_pairing_cpu_ns":cumulative_pairing_cpu,
            "cumulative_pairing_wall_ns_estimate":cumulative_pairing_wall,
            "all_in_pairing_cpu_ns":all_in_pairing_cpu,
            "all_in_pairing_wall_ns_l0":all_in_pairing_wall,
            "rho":rho,"rho_seed":final_seed,"rho_attempts":rho_attempts,
            "rho_cpu_ns":rho_cpu_ns,"rho_wall_ns":rho_wall_ns,
            "rho_fp_muls":rho_fp_muls,"rho_fp_invs":rho_fp_invs,
            "cumulative_rho_cpu_ns":cumulative_rho_cpu,
            "cumulative_rho_wall_ns":cumulative_rho_wall,
            "rho_step_cap":rho_cap,"six_rho":six,"six_attempts":six_attempts,
            "six_cpu_ns":six_cpu_ns,"six_wall_ns":six_wall_ns,
            "six_fp_muls":six_fp_muls,"six_fp_invs":six_fp_invs,
            "cumulative_six_cpu_ns":cumulative_six_cpu,
            "cumulative_six_wall_ns":cumulative_six_wall,
            "extension_cumulative":k.counts(),"peak_rss_raw":rss,"peak_rss_unit":unit})
        );
    }
    println!(
        "{}",
        json!({"schema":if cell==2 {"torsion_pairing_batch_end/v4"} else {"torsion_pairing_batch_end/v3"},"cell":cell,
        "round":round,"targets":targets.len(),"pairing_cpu_ns":cumulative_pairing_cpu,
        "rho_cpu_ns":cumulative_rho_cpu,"six_cpu_ns":cumulative_six_cpu,
        "process_cpu_ns":thread_cpu_ns()-whole_cpu,
        "process_wall_ns":whole_wall.elapsed().as_nanos() as u64,
        "extension_totals":k.counts(),"peak_rss_bytes":peak_rss_bytes()})
    );
}

fn verify_batch(
    input: Input,
    cell: usize,
    batch_path: &str,
    answers_path: &str,
    records_path: &str,
) {
    let targets = read_batch(batch_path, input, cell);
    let answers: Value = serde_json::from_str(&fs::read_to_string(answers_path).unwrap()).unwrap();
    assert_eq!(answers["schema"], "torsion_pairing_batch_answers/v1");
    assert_eq!(answers["cell"].as_u64(), Some(cell as u64));
    let scalars: Vec<u64> = answers["scalars"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_u64().unwrap())
        .collect();
    assert_eq!(scalars.len(), targets.len());
    let field = ReplayField::new(input);
    assert!(field.on_curve(input.g) && field.scalar_mul(input.g, input.r).infinity);
    let mut checked = 0u64;
    let mut rounds = [0u64; 3];
    for line in fs::read_to_string(records_path).unwrap().lines() {
        let v: Value = serde_json::from_str(line).unwrap();
        if v["schema"] != "torsion_pairing_batch_target/v3"
            && v["schema"] != "torsion_pairing_batch_target/v4" {
            continue;
        }
        assert_eq!(v["cell"].as_u64(), Some(cell as u64));
        let round = v["round"].as_u64().unwrap() as usize;
        assert!(round < 3);
        let index = v["index"].as_u64().unwrap() as usize;
        let scalar = v["recovered"].as_u64().unwrap();
        assert!(index < targets.len() && scalar < input.r);
        assert_eq!(v["target"], json!(targets[index]));
        assert_eq!(
            field.scalar_mul(input.g, scalar),
            targets[index],
            "independent source-curve replay failed"
        );
        assert_eq!(scalar, scalars[index], "post-solve fixture mismatch");
        assert_eq!(v["rho"]["recovered"].as_u64(), Some(scalar));
        assert_eq!(v["rho"]["verified"], true);
        if cell != 2 {
            assert_eq!(v["six_rho"]["recovered"].as_u64(), Some(scalar));
            assert_eq!(v["six_rho"]["verified"], true);
        } else {
            assert!(v["six_rho"].is_null());
        }
        rounds[round] += 1;
        checked += 1;
    }
    assert!(rounds.iter().all(|&n| n == targets.len() as u64));
    println!(
        "{}",
        json!({"schema":"torsion_pairing_batch_verify/v1","cell":cell,
        "checked":checked,"rounds":rounds,"targets_per_round":targets.len(),
        "independent_fp_muls":field.muls.get(),
        "independent_fp_invs":field.inversions.get(),
        "independent_group_adds":field.group_add_calls.get()})
    );
}

#[derive(Serialize)]
struct CompactPrefix {
    n: usize,
    cpu_ns: u64,
    wall_ns: u64,
    extension_muls: u64,
    extension_inversions: u64,
    rho_attempts: u64,
}

fn compact_pairing(
    input: Input,
    targets: &[SourcePoint],
    frozen: &Value,
    prefixes: &[usize],
) -> (Vec<CompactPrefix>, Vec<u64>, Value) {
    take_field_counters();
    let start_wall = Instant::now();
    let start_cpu = thread_cpu_ns();
    let curve = source_curve(input);
    let k = Field::new(input);
    assert_eq!(frozen["tower_a"], json!(k.a));
    assert_eq!(frozen["irreducibility"], k.irreducibility_receipt());
    let (h, selection, alpha, first_beta, _, _) = select_point(&k, input);
    assert_eq!(frozen["selected_h"], json!({"x":h.x,"y":h.y}));
    assert_eq!(frozen["selection"], selection);
    assert_eq!(frozen["alpha"], json!(alpha));
    assert_eq!(frozen["beta"], json!(first_beta));
    let table = ReusableBsgs::new(&k, &alpha, input.r);
    let mut output = Vec::with_capacity(targets.len());
    let mut samples = Vec::new();
    for (index, &q) in targets.iter().enumerate() {
        let beta = if index == 0 {
            first_beta.clone()
        } else {
            let mut charges = SelectionCharges::default();
            let value = reduced_tate_profiled(&k, q, &h, input.r, &mut charges)
                .expect("compact beta pairing pole");
            assert_eq!(k.pow(&value, &BigUint::from(input.r)), k.one());
            value
        };
        let scalar = table.solve(&k, &beta, input.r).expect("compact FFDLP exhausted");
        assert_eq!(k.pow(&alpha, &BigUint::from(scalar)), beta);
        let mut ops = GroupOps::default();
        assert_eq!(curve.mul(&mut ops, input.g, scalar), q);
        output.push(scalar);
        if prefixes.contains(&(index + 1)) {
            samples.push(CompactPrefix {
                n: index + 1,
                cpu_ns: thread_cpu_ns() - start_cpu,
                wall_ns: start_wall.elapsed().as_nanos() as u64,
                extension_muls: k.muls.get(),
                extension_inversions: k.inversions.get(),
                rho_attempts: 0,
            });
        }
    }
    let fp = take_field_counters();
    (samples, output, json!({"field_counts":k.counts(),
        "fp_muls":fp.0,"fp_invs":fp.1,
        "peak_rss_bytes":peak_rss_bytes()}))
}

fn compact_rho(
    input: Input,
    cell: usize,
    round: u64,
    targets: &[SourcePoint],
    prefixes: &[usize],
) -> (Vec<CompactPrefix>, Vec<u64>, Value) {
    take_field_counters();
    let start_wall = Instant::now();
    let start_cpu = thread_cpu_ns();
    let curve = source_curve(input);
    let cap = rho_step_cap(input.r);
    let mut attempts = 0u64;
    let mut output = Vec::with_capacity(targets.len());
    let mut samples = Vec::new();
    for (index, &q) in targets.iter().enumerate() {
        let seed_base = 0x4241_5443_4800_1000u64 + 256*cell as u64
            + 64*round + index as u64;
        let mut found = None;
        for retry in 0..32u64 {
            attempts += 1;
            let seed = seed_base.wrapping_add(retry * 0x01_00_00_00);
            let result = rho_reference_negation(&curve,input.g,q,input.r,seed,cap);
            if result.verified {
                found = result.recovered;
                break;
            }
        }
        let scalar = found.expect("compact charged rho exhausted");
        output.push(scalar);
        if prefixes.contains(&(index + 1)) {
            samples.push(CompactPrefix {
                n: index + 1,
                cpu_ns: thread_cpu_ns() - start_cpu,
                wall_ns: start_wall.elapsed().as_nanos() as u64,
                extension_muls: 0,
                extension_inversions: 0,
                rho_attempts: attempts,
            });
        }
    }
    let fp = take_field_counters();
    (samples, output, json!({"rho_attempts":attempts,
        "rho_step_cap":cap,"fp_muls":fp.0,"fp_invs":fp.1,
        "peak_rss_bytes":peak_rss_bytes()}))
}

fn run_compact(input: Input, cell: usize, round: u64, batch_path: &str, frozen_path: &str) {
    assert!(round < 3);
    let common_wall = Instant::now();
    let common_cpu = thread_cpu_ns();
    let targets = read_batch(batch_path,input,cell);
    let frozen: Value = serde_json::from_str(&fs::read_to_string(frozen_path).unwrap()).unwrap();
    let input_cpu_ns = thread_cpu_ns() - common_cpu;
    let input_wall_ns = common_wall.elapsed().as_nanos() as u64;
    let prefixes: Vec<usize> = if cell == 0 { vec![1,4,16] }
        else if cell == 1 { vec![1,4,16,60] }
        else { vec![1,4,16,64,66] };
    let (pair_prefix, pair_scalars, pair_counts, rho_prefix, rho_scalars, rho_counts) =
        if round % 2 == 0 {
            let (a,b,c) = compact_pairing(input,&targets,&frozen,&prefixes);
            let (d,e,f) = compact_rho(input,cell,round,&targets,&prefixes);
            (a,b,c,d,e,f)
        } else {
            let (d,e,f) = compact_rho(input,cell,round,&targets,&prefixes);
            let (a,b,c) = compact_pairing(input,&targets,&frozen,&prefixes);
            (a,b,c,d,e,f)
        };
    assert_eq!(pair_scalars,rho_scalars);
    assert_eq!(pair_prefix.len(),prefixes.len());
    assert_eq!(rho_prefix.len(),prefixes.len());
    println!("{}",json!({"schema":"torsion_pairing_compact/v5",
        "cell":cell,"p":input.p,"r":input.r,"round":round,
        "order":if round%2==0 {"pairing,rho"} else {"rho,pairing"},
        "targets":targets.len(),"prefixes":prefixes,
        "input_cpu_ns_common":input_cpu_ns,"input_wall_ns_common":input_wall_ns,
        "pairing_prefix":pair_prefix,"rho_prefix":rho_prefix,
        "pairing_scalars":pair_scalars,"rho_scalars":rho_scalars,
        "pairing_counts":pair_counts,"rho_counts":rho_counts,
        "peak_rss_bytes":peak_rss_bytes()}));
}

fn verify_compact(input: Input, cell: usize, batch_path: &str,
                  answers_path: &str, records_path: &str) {
    let targets = read_batch(batch_path,input,cell);
    let answers: Value = serde_json::from_str(&fs::read_to_string(answers_path).unwrap()).unwrap();
    let expected: Vec<u64> = answers["scalars"].as_array().unwrap().iter()
        .map(|v| v.as_u64().unwrap()).collect();
    assert_eq!(targets.len(),expected.len());
    let checker = ReplayField::new(input);
    let mut rounds = [0u64;3];
    for line in fs::read_to_string(records_path).unwrap().lines() {
        let row: Value = serde_json::from_str(line).unwrap();
        assert_eq!(row["schema"],"torsion_pairing_compact/v5");
        assert_eq!(row["cell"].as_u64(),Some(cell as u64));
        let round = row["round"].as_u64().unwrap() as usize;
        assert!(round < 3);
        for name in ["pairing_scalars","rho_scalars"] {
            let recovered = row[name].as_array().unwrap();
            assert_eq!(recovered.len(),targets.len());
            for (i,v) in recovered.iter().enumerate() {
                let scalar = v.as_u64().unwrap();
                assert_eq!(scalar,expected[i]);
                assert_eq!(checker.scalar_mul(input.g,scalar),targets[i]);
            }
        }
        rounds[round] += 1;
    }
    assert_eq!(rounds,[1,1,1]);
    println!("{}",json!({"schema":"torsion_pairing_compact_replay/v5",
        "cell":cell,"rounds":rounds,"checked_scalars":6*targets.len(),
        "independent_fp_muls":checker.muls.get(),
        "independent_fp_invs":checker.inversions.get(),
        "independent_group_adds":checker.group_add_calls.get()}));
}

// The replay uses separate F_p and F_(p²) arithmetic, as well as separate
// curve formulas. It shares only the public coordinates with the solver.
struct ReplayField {
    p: u64,
    nu: u64,
    muls: Cell<u64>,
    inversions: Cell<u64>,
    group_add_calls: Cell<u64>,
}

impl ReplayField {
    fn new(input: Input) -> Self {
        Self {
            p: input.p,
            nu: input.nu,
            muls: Cell::new(0),
            inversions: Cell::new(0),
            group_add_calls: Cell::new(0),
        }
    }
    fn m(&self, a: u64, b: u64) -> u64 {
        self.muls.set(self.muls.get() + 1);
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    fn addm(&self, a: u64, b: u64) -> u64 {
        (a + b) % self.p
    }
    fn subm(&self, a: u64, b: u64) -> u64 {
        (a + self.p - b) % self.p
    }
    fn add(&self, a: F2, b: F2) -> F2 {
        [self.addm(a[0], b[0]), self.addm(a[1], b[1])]
    }
    fn sub(&self, a: F2, b: F2) -> F2 {
        [self.subm(a[0], b[0]), self.subm(a[1], b[1])]
    }
    fn mul(&self, a: F2, b: F2) -> F2 {
        let real = self.addm(self.m(a[0], b[0]), self.m(self.nu, self.m(a[1], b[1])));
        let imag = self.addm(self.m(a[0], b[1]), self.m(a[1], b[0]));
        [real, imag]
    }
    fn base_pow(&self, mut a: u64, mut n: u64) -> u64 {
        let mut result = 1;
        while n > 0 {
            if n & 1 == 1 {
                result = self.m(result, a);
            }
            n >>= 1;
            if n > 0 {
                a = self.m(a, a);
            }
        }
        result
    }
    fn inv(&self, a: F2) -> F2 {
        let norm = self.subm(self.m(a[0], a[0]), self.m(self.nu, self.m(a[1], a[1])));
        assert_ne!(norm, 0);
        self.inversions.set(self.inversions.get() + 1);
        let ni = self.base_pow(norm, self.p - 2);
        let result = [self.m(a[0], ni), self.m(self.subm(0, a[1]), ni)];
        assert_eq!(self.mul(a, result), [1, 0]);
        result
    }
    fn on_curve(&self, p: SourcePoint) -> bool {
        p.infinity || self.mul(p.y, p.y) == self.add(self.mul(self.mul(p.x, p.x), p.x), [1, 0])
    }
    fn add_points(&self, a: SourcePoint, b: SourcePoint) -> SourcePoint {
        self.group_add_calls.set(self.group_add_calls.get() + 1);
        if a.infinity {
            return b;
        }
        if b.infinity {
            return a;
        }
        if a.x == b.x && self.add(a.y, b.y) == [0, 0] {
            return SourcePoint::infinity([0, 0]);
        }
        let slope = if a.x == b.x {
            self.mul(
                self.mul([3, 0], self.mul(a.x, a.x)),
                self.inv(self.mul([2, 0], a.y)),
            )
        } else {
            self.mul(self.sub(b.y, a.y), self.inv(self.sub(b.x, a.x)))
        };
        let x = self.sub(self.sub(self.mul(slope, slope), a.x), b.x);
        let y = self.sub(self.mul(slope, self.sub(a.x, x)), a.y);
        SourcePoint::affine(x, y)
    }
    fn scalar_mul(&self, p: SourcePoint, mut k: u64) -> SourcePoint {
        let mut out = SourcePoint::infinity([0, 0]);
        let mut base = p;
        while k > 0 {
            if k & 1 == 1 {
                out = self.add_points(out, base);
            }
            k >>= 1;
            if k > 0 {
                base = self.add_points(base, base);
            }
        }
        out
    }
}

fn verify_result(input: Input, cell: usize, path: &str) {
    let v: Value = serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap();
    assert_eq!(v["cell"].as_u64(), Some(cell as u64));
    let scalar = v["recovered"].as_u64().expect("candidate scalar");
    assert!(scalar < input.r);
    let wall = Instant::now();
    let cpu = thread_cpu_ns();
    let field = ReplayField::new(input);
    assert!(field.on_curve(input.g) && field.on_curve(input.q));
    assert!(field.scalar_mul(input.g, input.r).infinity);
    let answer = field.scalar_mul(input.g, scalar);
    assert_eq!(answer, input.q, "independent source-curve replay failed");
    let cpu_ns = thread_cpu_ns() - cpu;
    let wall_ns = wall.elapsed().as_nanos() as u64;
    let (rss, rss_unit) = peak_rss_raw();
    println!(
        "{}",
        json!({"schema":"torsion_pairing_verify/v1","cell":cell,"p":input.p,
        "r":input.r,"recovered":scalar,"verified":true,
        "source_group_add_calls":field.group_add_calls.get(),
        "cpu_ns":cpu_ns,"wall_ns":wall_ns,
        "fp_muls":field.muls.get(),"fp_invs":field.inversions.get(),
        "peak_rss_raw":rss,"peak_rss_unit":rss_unit})
    );
}

fn check_bilinearity(input: Input, cell: usize, frozen_path: &str) {
    let frozen: Value = serde_json::from_str(&fs::read_to_string(frozen_path).unwrap()).unwrap();
    assert_eq!(frozen["cell"].as_u64(), Some(cell as u64));
    take_field_counters();
    let wall = Instant::now();
    let cpu = thread_cpu_ns();
    let k = Field::new(input);
    assert_eq!(frozen["tower_a"], json!(k.a));
    assert_eq!(frozen["irreducibility"], k.irreducibility_receipt());
    let parse_elt = |v: &Value| -> Elt {
        v.as_array()
            .unwrap()
            .iter()
            .map(|c| {
                let pair = c.as_array().unwrap();
                [pair[0].as_u64().unwrap(), pair[1].as_u64().unwrap()]
            })
            .collect()
    };
    let h = Point::affine(
        parse_elt(&frozen["selected_h"]["x"]),
        parse_elt(&frozen["selected_h"]["y"]),
    );
    assert!(point_mul(&k, &h, &BigUint::from(input.r)).infinity);
    let alpha = parse_elt(&frozen["alpha"]);
    let curve = source_curve(input);
    let mut source_ops = GroupOps::default();
    let mut checks = Vec::new();
    for scalar in [1u64, 2, 3, 5] {
        let multiple = curve.mul(&mut source_ops, input.g, scalar);
        let value = reduced_tate(&k, multiple, &h, input.r)
            .expect("bilinearity check encountered a pairing pole");
        assert_eq!(
            value,
            k.pow(&alpha, &BigUint::from(scalar)),
            "pairing bilinearity check failed"
        );
        checks.push(json!({"scalar":scalar,"passed":true}));
    }
    let fp = take_field_counters();
    let (rss, rss_unit) = peak_rss_raw();
    println!(
        "{}",
        json!({"schema":"torsion_pairing_bilinearity/v1",
        "cell":cell,"p":input.p,"r":input.r,"checks":checks,
        "source_group_ops":source_ops,"field_counts":k.counts(),
        "fp_muls":fp.0,"fp_invs":fp.1,
        "cpu_ns":thread_cpu_ns()-cpu,"wall_ns":wall.elapsed().as_nanos() as u64,
        "peak_rss_raw":rss,"peak_rss_unit":rss_unit})
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("freeze-401") {
        freeze_401(&args[2], &args[3], &args[4]);
        return;
    }
    assert!(
        args.len() >= 4,
        "usage: torsion_pairing_transfer prepare|run INPUTS.jsonl CELL"
    );
    let cell: usize = args[3].parse().unwrap();
    let input = read_input(&args[2], cell);
    match args[1].as_str() {
        "freeze-401" => unreachable!(),
        "freeze-batch" => freeze_batch(
            input,
            cell,
            args.get(4).expect("fixture path"),
            args.get(5).expect("public output path"),
            args.get(6).expect("answer output path"),
        ),
        "batch-run" => run_batch(
            input,
            cell,
            args.get(4).expect("round").parse().unwrap(),
            args.get(5).expect("public batch path"),
            args.get(6).expect("selection path"),
        ),
        "compact-run" => run_compact(
            input, cell, args.get(4).expect("round").parse().unwrap(),
            args.get(5).expect("public batch path"),
            args.get(6).expect("selection path"),
        ),
        "compact-verify" => verify_compact(
            input, cell, args.get(4).expect("public batch path"),
            args.get(5).expect("answer path"),
            args.get(6).expect("records path"),
        ),
        "batch-verify" => verify_batch(
            input,
            cell,
            args.get(4).expect("public batch path"),
            args.get(5).expect("answer path"),
            args.get(6).expect("records path"),
        ),
        "prepare" => run_prepare(input, cell),
        "run" => run_pairing(input, cell, args.get(4).expect("frozen selection path")),
        "rho" => run_rho(input, cell, args.get(4).expect("round").parse().unwrap()),
        "verify" => verify_result(input, cell, args.get(4).expect("run record path")),
        "bilinear" => check_bilinearity(input, cell, args.get(4).expect("frozen selection path")),
        _ => panic!("unknown command"),
    }
}
