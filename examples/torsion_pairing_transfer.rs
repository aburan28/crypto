//! Native, cold Tate-pairing transfer on the frozen quadratic-extension targets.
//! The input file deliberately excludes the known logarithm. See the study
//! protocol for field, point-selection, accounting, and stop conditions.

use std::cell::Cell;
use std::collections::BTreeMap;
use std::fs;
use std::time::Instant;

use crypto_lib::cryptanalysis::ext_curve::{
    take_field_counters, ExtCurve, ExtField, ExtPoint, Fp2,
};
use crypto_lib::cryptanalysis::ic_boundary::{rho_reference_negation, CountedGroup, GroupOps};
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use rand::rngs::StdRng;
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
    assert!((p == 61 && r == 19 && cell == 0) || (p == 271 && r == 61 && cell == 1));
    Input {
        p,
        nu: v["nu"].as_u64().unwrap(),
        r,
        degree: if cell == 0 { 9 } else { 5 },
        trace_fp: v["trace_fp"].as_i64().unwrap(),
        g: source_point(&v["generator"]),
        q: source_point(&v["target"]),
        seed: if cell == 0 {
            0x5041_4952_0000_003d
        } else {
            0x5041_4952_0000_010f
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
        let prime_divisor = if input.degree == 9 { 3 } else { 5 };
        let a = (1..q)
            .map(|i| f.from_coords(&[i % input.p, i / input.p]))
            .find(|&x| f.pow(x, (q - 1) / prime_divisor) != f.one())
            .expect("non-power coefficient");
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

fn select_point(k: &Field, input: Input) -> (Point, Value, Elt, Elt, SelectionCharges) {
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
        return (h, receipt, alpha, beta, charges);
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
    let (h, selection, alpha, beta, charges) = select_point(&k, input);
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
    let (h, selection, alpha, beta, charges) = select_point(&k, input);
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
    for scalar in [2u64, 3, 5] {
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
    assert!(
        args.len() >= 4,
        "usage: torsion_pairing_transfer prepare|run INPUTS.jsonl CELL"
    );
    let cell: usize = args[3].parse().unwrap();
    let input = read_input(&args[2], cell);
    match args[1].as_str() {
        "prepare" => run_prepare(input, cell),
        "run" => run_pairing(input, cell, args.get(4).expect("frozen selection path")),
        "rho" => run_rho(input, cell, args.get(4).expect("round").parse().unwrap()),
        "verify" => verify_result(input, cell, args.get(4).expect("run record path")),
        "bilinear" => check_bilinearity(input, cell, args.get(4).expect("frozen selection path")),
        _ => panic!("unknown command"),
    }
}
