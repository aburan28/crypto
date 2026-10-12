//! Library entry points over prime fields given as integers: the operations the `isogeny-algos`
//! command line exposes, usable directly from Rust. Every result carries the checks that were
//! run on it; nothing here is reported without a verification status.
//!
//! * [`count_points`]: #E(F_p) by SEA (Elkies/Atkin primes, isogeny cycles), by complex
//!   multiplication for j = 0 and 1728, or by baby-step giant-step for p < 2^40, followed by
//!   [`certify_order`].
//! * [`certify_order`]: a claimed #E checked against random points, and *certified* when the
//!   orders of the sampled points force it: if D divides the order of a point and D > 4 sqrt(p),
//!   only one multiple of D lies in the Hasse interval. D comes from a partial factorisation of
//!   the claimed order (trial division, Pollard–Brent, Miller–Rabin), so the certificate is as
//!   strong as those probable-prime tests.
//! * [`isogenies`]: the F_p-rational l-isogenies of E with kernel polynomials, codomains and
//!   rational maps, each map checked on random points (lands on the codomain, respects sums).
//! * [`modular_polynomial`] and [`modular_polynomial_integer`]: Phi_l mod p and over Z.
//! * [`isogeny_from_kernel`]: Kohel/Velu from a kernel polynomial or a kernel x-coordinate.
//!
//! Fields are chosen by the size of p: u64 arithmetic below 2^62, and Montgomery with
//! exactly the required number of limbs, up to 640 bits.
use crate::bigint::Big;
use crate::curve::{
    jinv, on_curve, padd, pmul_big, random_point_f, Curve, Isogeny, Pt, RatIsogeny,
};
use crate::field::{Field, Rng, Zp};
use crate::find::divpoly::kernel_polys;
use crate::find::modpoly::{integer_coeffs, Phi};
use crate::fpm::FpM;
use crate::int::Int;
use crate::kernel::kohel::kohel;
use crate::poly;
use std::collections::HashMap;

/// y^2 = x^3 + a x + b over GF(p), p an odd prime above 3, coefficients reduced mod p.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PrimeCurve {
    pub p: Int,
    pub a: Int,
    pub b: Int,
}

impl PrimeCurve {
    /// Checks that p is a (probable) prime above 3 and the curve is non-singular.
    pub fn new(p: Int, a: Int, b: Int) -> Result<Self, String> {
        if p <= Int::from(3i64) || !p.is_probable_prime() {
            return Err(format!("p = {p} is not a prime above 3"));
        }
        if p.bits() > 640 {
            return Err("p above 640 bits is not supported".into());
        }
        let (a, b) = (a.modulo(&p), b.modulo(&p));
        let c = PrimeCurve { p, a, b };
        if c.discriminant().is_zero() {
            return Err("singular curve (4a^3 + 27b^2 = 0 mod p)".into());
        }
        Ok(c)
    }
    /// Standard curves by name: p192, p224, p256 (NIST), secp256k1.
    pub fn preset(name: &str) -> Option<Self> {
        let (p, a, b) = match name.to_ascii_lowercase().as_str() {
            "p192" | "p-192" | "secp192r1" => (
                "6277101735386680763835789423207666416083908700390324961279",
                "6277101735386680763835789423207666416083908700390324961276",
                "2455155546008943817740293915197451784769108058161191238065",
            ),
            "p224" | "p-224" | "secp224r1" => (
                "26959946667150639794667015087019630673557916260026308143510066298881",
                "26959946667150639794667015087019630673557916260026308143510066298878",
                "18958286285566608000408668544493926415504680968679321075787234672564",
            ),
            "p256" | "p-256" | "secp256r1" => (
                "115792089210356248762697446949407573530086143415290314195533631308867097853951",
                "115792089210356248762697446949407573530086143415290314195533631308867097853948",
                "41058363725152142129326129780047268409114441015993725554835256314039467401291",
            ),
            "secp256k1" => (
                "115792089237316195423570985008687907853269984665640564039457584007908834671663",
                "0",
                "7",
            ),
            _ => return None,
        };
        let d = |s: &str| Int::from_big(&Big::from_dec(s));
        PrimeCurve::new(d(p), d(a), d(b)).ok()
    }
    fn discriminant(&self) -> Int {
        let p = &self.p;
        let a3 = self.a.pow_mod(&Int::from(3i64), p);
        (&(&Int::from(4i64) * &a3) + &(&Int::from(27i64) * &(&self.b * &self.b))).modulo(p)
    }
    /// j = 1728 * 4a^3 / (4a^3 + 27b^2) mod p.
    pub fn j(&self) -> Int {
        let p = &self.p;
        let a3 = self.a.pow_mod(&Int::from(3i64), p);
        let inv = self.discriminant().inv_mod(p).expect("non-singular");
        (&(&Int::from(6912i64) * &a3) * &inv).modulo(p)
    }
    /// The quadratic twist y^2 = x^3 + a d^2 x + b d^3 by the least non-residue d.
    pub fn twist(&self) -> PrimeCurve {
        let p = &self.p;
        let half = (p - &Int::one()).div_floor(&Int::from(2i64));
        let mut d = Int::from(2i64);
        while d.pow_mod(&half, p) == Int::one() {
            d = &d + &Int::one();
        }
        let d2 = (&d * &d).modulo(p);
        let d3 = (&d2 * &d).modulo(p);
        PrimeCurve {
            p: p.clone(),
            a: (&self.a * &d2).modulo(p),
            b: (&self.b * &d3).modulo(p),
        }
    }
}

/// A prime field whose elements convert from and to integers.
pub trait PrimeFieldInt: Field {
    fn elem(&self, v: &Int) -> Self::E;
    fn int(&self, e: Self::E) -> Int;
}
impl PrimeFieldInt for Zp {
    fn elem(&self, v: &Int) -> u64 {
        v.modulo(&Int::from(self.p)).to_i128().unwrap() as u64
    }
    fn int(&self, e: u64) -> Int {
        Int::from(e)
    }
}
impl<const N: usize> PrimeFieldInt for FpM<N> {
    fn elem(&self, v: &Int) -> [u64; N] {
        let m = Int::from_big(self.modulus());
        self.from_big(v.modulo(&m).mag())
    }
    fn int(&self, e: [u64; N]) -> Int {
        Int::from_big(&self.to_big(&e))
    }
}

/// Work generic over the field type, run by [`with_field`].
pub trait FieldTask {
    type Out;
    fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out;
}

/// Run `task` over GF(p) with the smallest field implementation that holds p.
pub fn with_field<T: FieldTask>(p: &Int, task: T) -> T::Out {
    let bits = p.bits();
    if bits <= 62 {
        task.run(&Zp::new(p.to_i128().unwrap() as u64))
    } else if bits <= 64 {
        task.run(&FpM::<1>::from_dec(&p.to_string()))
    } else if bits <= 128 {
        task.run(&FpM::<2>::from_dec(&p.to_string()))
    } else if bits <= 192 {
        task.run(&FpM::<3>::from_dec(&p.to_string()))
    } else if bits <= 256 {
        task.run(&FpM::<4>::from_dec(&p.to_string()))
    } else if bits <= 320 {
        task.run(&FpM::<5>::from_dec(&p.to_string()))
    } else if bits <= 384 {
        task.run(&FpM::<6>::from_dec(&p.to_string()))
    } else if bits <= 448 {
        task.run(&FpM::<7>::from_dec(&p.to_string()))
    } else if bits <= 512 {
        task.run(&FpM::<8>::from_dec(&p.to_string()))
    } else if bits <= 576 {
        task.run(&FpM::<9>::from_dec(&p.to_string()))
    } else {
        assert!(bits <= 640, "prime fields above 640 bits are unsupported");
        task.run(&FpM::<10>::from_dec(&p.to_string()))
    }
}

fn curve_in<F: PrimeFieldInt>(f: &F, c: &PrimeCurve) -> Curve<F::E> {
    Curve::new(f.elem(&c.a), f.elem(&c.b))
}
fn point_out<F: PrimeFieldInt>(f: &F, p: &Pt<F::E>) -> Option<(Int, Int)> {
    match *p {
        Pt::Inf => None,
        Pt::Aff(x, y) => Some((f.int(x), f.int(y))),
    }
}

// ------------------------------------------------------------------ checks

/// The status vocabulary of the curve-audit records in the toolset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Pass,
    Fail,
    NotRun,
    Indeterminate,
}
impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Status::Pass => "PASS",
            Status::Fail => "FAIL",
            Status::NotRun => "NOT_RUN",
            Status::Indeterminate => "INDETERMINATE",
        }
    }
}

/// One named check with its outcome and a short explanation.
#[derive(Clone, Debug)]
pub struct Check {
    pub name: String,
    pub status: Status,
    pub detail: String,
}
impl Check {
    fn new(name: &str, status: Status, detail: String) -> Check {
        Check {
            name: name.into(),
            status,
            detail,
        }
    }
}

// ------------------------------------------------------------------ integers

fn rand_below(rng: &mut Rng, n: &Int) -> Int {
    let limbs: Vec<u64> = (0..n.bits() / 64 + 2).map(|_| rng.next()).collect();
    Int::from_big(&Big::from_limbs(&limbs)).modulo(n)
}

/// One non-trivial factor of the odd composite n by Pollard–Brent, or None within `budget`
/// iterations.
fn pollard_brent(n: &Int, budget: u64, rng: &mut Rng) -> Option<Int> {
    let one = Int::one();
    for _ in 0..4 {
        let c = &rand_below(rng, &(n - &one)) + &one;
        let f = |v: &Int| (&(v * v) + &c).modulo(n);
        let mut y = rand_below(rng, n);
        let (mut g, mut r, mut q) = (one.clone(), 1u64, one.clone());
        let (mut x, mut ys) = (y.clone(), y.clone());
        let mut used = 0u64;
        while g == one {
            x = y.clone();
            for _ in 0..r {
                y = f(&y);
            }
            let mut k = 0;
            while k < r && g == one {
                ys = y.clone();
                for _ in 0..128.min(r - k) {
                    y = f(&y);
                    q = (&q * &(&x - &y).abs()).modulo(n);
                }
                g = Int::gcd(&q, n);
                k += 128;
                used += 128;
                if used > budget {
                    return None;
                }
            }
            r *= 2;
        }
        if g == *n {
            loop {
                ys = f(&ys);
                g = Int::gcd(&(&x - &ys).abs(), n);
                if g > one {
                    break;
                }
            }
        }
        if g != *n {
            return Some(g);
        }
    }
    None
}

/// Prime powers of n found by trial division, Pollard–Brent and Miller–Rabin, and the
/// remaining composite part (1 if n is fully factored).
pub fn factor_partial(n: &Int, budget: u64, rng: &mut Rng) -> (Vec<(Int, u32)>, Int) {
    let mut found: Vec<Int> = vec![];
    let mut m = n.abs();
    let mut d = 2u64;
    while d < 1 << 16 && m > Int::one() {
        let dd = Int::from(d);
        while m.mod_u64(d) == 0 {
            found.push(dd.clone());
            m = m.div_floor(&dd);
        }
        d += if d == 2 { 1 } else { 2 };
    }
    let mut stack = vec![m];
    let mut rest = Int::one();
    while let Some(c) = stack.pop() {
        if c == Int::one() {
            continue;
        }
        if c.is_probable_prime() {
            found.push(c);
            continue;
        }
        match pollard_brent(&c, budget, rng) {
            Some(g) => {
                let other = c.div_floor(&g);
                stack.push(g);
                stack.push(other);
            }
            None => rest = &rest * &c,
        }
    }
    found.sort();
    let mut out: Vec<(Int, u32)> = vec![];
    for f in found {
        match out.last_mut() {
            Some((q, e)) if *q == f => *e += 1,
            _ => out.push((f, 1)),
        }
    }
    (out, rest)
}

// ------------------------------------------------------------------ point counts

/// How a group order was obtained.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CountMethod {
    /// baby-step giant-step in the Hasse interval (p < 2^40)
    Bsgs,
    /// complex multiplication, j = 0 or 1728 (Cornacchia, trace chosen by points)
    Cm,
    /// Schoof–Elkies–Atkin; the primes used and their congruences
    Sea {
        elkies: Vec<(u64, u64)>,
        atkin: Vec<(u64, usize, usize)>,
        candidates: u128,
    },
    /// supplied by the caller (only checked)
    Given,
}

/// A group order with the evidence for it.
#[derive(Clone, Debug)]
pub struct OrderCertificate {
    pub order: Int,
    /// the order is the only multiple of `divisor` in the Hasse interval, or (`via_twist`) the
    /// twist's order is certified and #E + #E^t = 2p + 2 fixes this one
    pub certified: bool,
    /// certified through the twist rather than by this curve's own point orders
    pub via_twist: bool,
    /// lcm of divisors of point orders established by the sampled points
    pub divisor: Int,
    /// prime factorisation of the claimed order found (probable primes), and its unfactored part
    pub factors: Vec<(Int, u32)>,
    pub unfactored: Int,
    /// sample points (x, y) used
    pub points: Vec<(Int, Int)>,
    pub checks: Vec<Check>,
}

/// Check a claimed #E(F_p) and certify it when the sampled point orders force it.
pub fn certify_order(curve: &PrimeCurve, order: &Int, seed: u64) -> OrderCertificate {
    struct T<'a>(&'a PrimeCurve, &'a Int, u64);
    impl FieldTask for T<'_> {
        type Out = OrderCertificate;
        fn run<F: PrimeFieldInt>(self, f: &F) -> OrderCertificate {
            certify_in(f, self.0, self.1, self.2)
        }
    }
    with_field(&curve.p, T(curve, order, seed))
}

fn certify_in<F: PrimeFieldInt>(f: &F, c: &PrimeCurve, order: &Int, seed: u64) -> OrderCertificate {
    let p = &c.p;
    let e = curve_in(f, c);
    let mut rng = Rng::new(seed ^ 0xCE27);
    let mut checks = vec![];
    let t = &(p + &Int::one()) - order;
    let hasse = &t * &t <= &Int::from(4i64) * p;
    checks.push(Check::new(
        "hasse_interval",
        if hasse { Status::Pass } else { Status::Fail },
        format!(
            "|p + 1 - N| = |{t}| {} 2 sqrt(p)",
            if hasse { "<=" } else { ">" }
        ),
    ));
    let (factors, unfactored) = factor_partial(order, 1 << 18, &mut rng);
    let mut divisor = Int::one();
    let mut points = vec![];
    let mut annihilated = hasse;
    let bound_sq = &Int::from(16i64) * p; // D > 4 sqrt(p)  <=>  D^2 > 16 p
    let ord_big = order.mag().clone();
    for _ in 0..12 {
        if !annihilated || &divisor * &divisor > bound_sq {
            break;
        }
        let pt = random_point_f(f, &e, &mut rng);
        if pmul_big(f, &e, &pt, &ord_big) != Pt::Inf {
            annihilated = false;
            break;
        }
        if let Some(xy) = point_out(f, &pt) {
            points.push(xy);
        }
        // r-part of ord(P) for every prime r found in N
        let mut d = Int::one();
        for (r, k) in &factors {
            let rk = r.pow(*k);
            let mut q = pmul_big(f, &e, &pt, order.div_floor(&rk).mag());
            let mut v = 0;
            while q != Pt::Inf && v < *k {
                q = pmul_big(f, &e, &q, r.mag());
                v += 1;
            }
            d = &d * &r.pow(v);
        }
        divisor = lcm(&divisor, &d);
    }
    checks.push(Check::new(
        "points_annihilated",
        if annihilated {
            Status::Pass
        } else {
            Status::Fail
        },
        format!("[N]P = O for {} random points", points.len()),
    ));
    let certified = annihilated && &divisor * &divisor > bound_sq;
    checks.push(Check::new(
        "order_unique_in_hasse_interval",
        if certified {
            Status::Pass
        } else if annihilated {
            Status::Indeterminate
        } else {
            Status::Fail
        },
        if certified {
            format!("a divisor D = {divisor} of the point orders exceeds 4 sqrt(p); prime factors by Miller-Rabin")
        } else {
            format!("divisor of the point orders established: {divisor}; unfactored part of N: {unfactored}")
        },
    ));
    OrderCertificate {
        order: order.clone(),
        certified,
        via_twist: false,
        divisor,
        factors,
        unfactored,
        points,
        checks,
    }
}

fn lcm(a: &Int, b: &Int) -> Int {
    (a * b).div_floor(&Int::gcd(a, b))
}

/// #E(F_p), its trace and twist, with the method and the certificate.
#[derive(Clone, Debug)]
pub struct PointCount {
    pub order: Int,
    pub trace: Int,
    pub twist_order: Int,
    pub method: CountMethod,
    pub certificate: OrderCertificate,
    pub twist_certificate: OrderCertificate,
}

/// Count points (method by size and j), then certify the order and the twist's order.
pub fn count_points(curve: &PrimeCurve, seed: u64) -> Result<PointCount, String> {
    struct T<'a>(&'a PrimeCurve, u64);
    impl FieldTask for T<'_> {
        type Out = Result<(Int, CountMethod), String>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            count_in(f, self.0, self.1)
        }
    }
    let (order, method) = with_field(&curve.p, T(curve, seed))?;
    Ok(finish_count(curve, order, method, seed))
}

/// Certify a supplied order (and the twist's) without counting.
pub fn check_order(curve: &PrimeCurve, order: &Int, seed: u64) -> PointCount {
    finish_count(curve, order.clone(), CountMethod::Given, seed)
}

fn finish_count(curve: &PrimeCurve, order: Int, method: CountMethod, seed: u64) -> PointCount {
    let p = &curve.p;
    let trace = &(p + &Int::one()) - &order;
    let twist_order = &(p + &Int::one()) + &trace;
    let mut certificate = certify_order(curve, &order, seed);
    let mut twist_certificate = certify_order(&curve.twist(), &twist_order, seed ^ 0x7715);
    // #E + #E^t = 2p + 2: a certificate for either order fixes the other
    let (ce, ct) = (certificate.certified, twist_certificate.certified);
    let points_ok = |c: &OrderCertificate| c.checks.iter().all(|k| k.status != Status::Fail);
    for (cert, other_certified, other) in [
        (&mut certificate, ct, "twist's"),
        (&mut twist_certificate, ce, "curve's"),
    ] {
        if !cert.certified && other_certified && points_ok(cert) {
            cert.certified = true;
            cert.via_twist = true;
            for k in cert
                .checks
                .iter_mut()
                .filter(|k| k.name == "order_unique_in_hasse_interval")
            {
                k.status = Status::Pass;
                k.detail = format!(
                    "follows from the certified {other} order: #E + #E^t = 2p + 2 ({})",
                    k.detail
                );
            }
        }
    }
    PointCount {
        order,
        trace,
        twist_order,
        method,
        certificate,
        twist_certificate,
    }
}

fn count_in<F: PrimeFieldInt>(
    f: &F,
    c: &PrimeCurve,
    seed: u64,
) -> Result<(Int, CountMethod), String> {
    let mut rng = Rng::new(seed);
    let p = &c.p;
    if p.bits() <= 40 {
        let fp = Zp::new(p.to_i128().unwrap() as u64);
        let e = Curve::new(fp.elem(&c.a), fp.elem(&c.b));
        return Ok((
            Int::from(crate::curve::order(&fp, &e, &mut rng)),
            CountMethod::Bsgs,
        ));
    }
    if c.a.is_zero() || c.b.is_zero() {
        return cm_count(f, c, &mut rng).map(|n| (n, CountMethod::Cm));
    }
    let e = curve_in(f, c);
    let mut phis: HashMap<usize, Phi<F>> = HashMap::new();
    let (n, st) =
        crate::find::sea::sea(f, &e, 1000, &mut phis, &mut rng).ok_or("SEA found no order")?;
    Ok((
        Int::from_big(&n),
        CountMethod::Sea {
            elkies: st.elkies,
            atkin: st.atkin,
            candidates: st.candidates,
        },
    ))
}

/// x^2 + d y^2 = p (d = 1, 3) by Cornacchia.
fn cornacchia(p: &Int, d: i64) -> Option<(Int, Int)> {
    let r0 = Int::sqrt_mod_prime(&Int::from(-d).modulo(p), p)?;
    let (mut a, mut b) = (p.clone(), r0);
    let lim = p.isqrt();
    while b > lim {
        let r = a.modulo(&b);
        a = b;
        b = r;
    }
    let rest = p - &(&b * &b);
    if rest.modulo(&Int::from(d)) != Int::zero() {
        return None;
    }
    let y2 = rest.div_floor(&Int::from(d));
    let y = y2.isqrt();
    (&y * &y == y2).then_some((b, y))
}

/// j = 0 (a = 0) or j = 1728 (b = 0): the trace is one of the norms' associates; the points
/// choose it.
fn cm_count<F: PrimeFieldInt>(f: &F, c: &PrimeCurve, rng: &mut Rng) -> Result<Int, String> {
    let p = &c.p;
    let p1 = p + &Int::one();
    let j0 = c.a.is_zero();
    let (supersingular, d) = if j0 {
        (p.mod_u64(3) == 2, 3)
    } else {
        (p.mod_u64(4) == 3, 1)
    };
    if supersingular {
        return Ok(p1);
    }
    let (x, y) = cornacchia(p, d).ok_or("Cornacchia failed")?;
    let mut traces: Vec<Int> = if j0 {
        let three_y = &Int::from(3i64) * &y;
        vec![&x + &x, &x + &three_y, &x - &three_y]
    } else {
        vec![&x + &x, &y + &y]
    };
    traces.extend(traces.clone().into_iter().map(|t| -t));
    let e = curve_in(f, c);
    let mut alive: Vec<Int> = traces.into_iter().map(|t| &p1 - &t).collect();
    alive.sort();
    alive.dedup();
    for _ in 0..40 {
        if alive.len() <= 1 {
            break;
        }
        let pt = random_point_f(f, &e, rng);
        alive.retain(|n| pmul_big(f, &e, &pt, n.mag()) == Pt::Inf);
    }
    match alive.len() {
        1 => Ok(alive.pop().unwrap()),
        0 => Err("no CM candidate annihilates the sampled points".into()),
        _ => Err("several CM candidates remain".into()),
    }
}

// ------------------------------------------------------------------ isogenies

/// One rational l-isogeny E -> E'.
#[derive(Clone, Debug)]
pub struct IsogenyRecord {
    pub degree: u64,
    /// codomain y^2 = x^3 + a' x + b'
    pub codomain: (Int, Int),
    pub j_codomain: Int,
    /// monic kernel polynomial, coefficients low to high
    pub kernel: Vec<Int>,
    /// x-map num(x) / den(x), coefficients low to high (den = kernel^2)
    pub map_num: Vec<Int>,
    pub map_den: Vec<Int>,
    pub checks: Vec<Check>,
}

impl IsogenyRecord {
    pub fn verified(&self) -> bool {
        self.checks.iter().all(|c| c.status == Status::Pass)
    }
}

fn record_from<F: PrimeFieldInt>(
    f: &F,
    e: &Curve<F::E>,
    iso: &RatIsogeny<F>,
    phi: Option<&Phi<F>>,
    rng: &mut Rng,
) -> IsogenyRecord {
    let cod = iso.cod;
    let mut checks = vec![];
    let mut on = 0;
    let mut hom = 0;
    let trials = 8;
    for _ in 0..trials {
        let (p1, p2) = (random_point_f(f, e, rng), random_point_f(f, e, rng));
        let (i1, i2) = (iso.eval(f, &p1), iso.eval(f, &p2));
        if on_curve(f, &cod, &i1) && on_curve(f, &cod, &i2) {
            on += 1;
        }
        if iso.eval(f, &padd(f, e, &p1, &p2)) == padd(f, &cod, &i1, &i2) {
            hom += 1;
        }
    }
    checks.push(Check::new(
        "image_on_codomain",
        if on == trials {
            Status::Pass
        } else {
            Status::Fail
        },
        format!("{on}/{trials} pairs of random points map onto E'"),
    ));
    checks.push(Check::new(
        "homomorphism",
        if hom == trials {
            Status::Pass
        } else {
            Status::Fail
        },
        format!("phi(P + Q) = phi(P) + phi(Q) for {hom}/{trials} random pairs"),
    ));
    let kdeg = iso.ker.len().saturating_sub(1) as u64;
    let want = if iso.deg == 2 { 1 } else { (iso.deg - 1) / 2 };
    checks.push(Check::new(
        "kernel_degree",
        if kdeg == want {
            Status::Pass
        } else {
            Status::Fail
        },
        format!(
            "kernel polynomial of degree {kdeg} for an isogeny of degree {}",
            iso.deg
        ),
    ));
    let jt = jinv(f, &cod);
    if let Some(phi) = phi {
        let ok = f.is_zero(poly::eval(f, &phi.y_poly(f, jinv(f, e)), jt));
        checks.push(Check::new(
            "modular_polynomial",
            if ok { Status::Pass } else { Status::Fail },
            "Phi_l(j(E), j(E')) = 0".into(),
        ));
    }
    let ints = |v: &[F::E]| v.iter().map(|&x| f.int(x)).collect::<Vec<Int>>();
    IsogenyRecord {
        degree: iso.deg,
        codomain: (f.int(cod.a), f.int(cod.b)),
        j_codomain: f.int(jt),
        kernel: ints(&iso.ker),
        map_num: ints(&iso.num),
        map_den: ints(&iso.den),
        checks,
    }
}

/// The F_p-rational isogenies of prime degree l from E: l = 2 from the rational 2-torsion; odd
/// l from the roots of Phi_l(j(E), Y) (Elkies codomain, BMSS kernel) or, for j = 0, 1728 and
/// for small fields where Phi_l does not apply, from the factors of the division polynomial.
pub fn isogenies(curve: &PrimeCurve, ell: u64, seed: u64) -> Result<Vec<IsogenyRecord>, String> {
    if ell < 2 || !crate::field::is_prime(ell) {
        return Err(format!("l = {ell} is not prime"));
    }
    struct T<'a>(&'a PrimeCurve, u64, u64);
    impl FieldTask for T<'_> {
        type Out = Result<Vec<IsogenyRecord>, String>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            isogenies_in(f, self.0, self.1, self.2)
        }
    }
    with_field(&curve.p, T(curve, ell, seed))
}

fn isogenies_in<F: PrimeFieldInt>(
    f: &F,
    c: &PrimeCurve,
    ell: u64,
    seed: u64,
) -> Result<Vec<IsogenyRecord>, String> {
    let mut rng = Rng::new(seed);
    let e = curve_in(f, c);
    if ell == 2 {
        let cubic = vec![e.b, e.a, f.zero(), f.one()];
        let roots = poly::roots(f, &cubic, &mut rng);
        return Ok(roots
            .into_iter()
            .map(|x0| {
                let iso = kohel(f, &e, &vec![f.neg(x0), f.one()], 2);
                record_from(f, &e, &iso, None, &mut rng)
            })
            .collect());
    }
    let special = c.a.is_zero() || c.b.is_zero();
    let phi_ok = c.p > Int::from(ell as i64 + 1);
    if !special && phi_ok {
        let phi = Phi::compute(f, ell as usize);
        let j = jinv(f, &e);
        let mut out = vec![];
        let mut failed = false;
        for jt in phi.neighbors(f, j, &mut rng) {
            let iso = crate::find::elkies::elkies_codomain(f, &phi, &e, jt).and_then(|et| {
                crate::find::bmss::isogeny(
                    f,
                    crate::find::bmss::Method::FastElkiesPrime,
                    &e,
                    &et,
                    ell as usize,
                    None,
                )
            });
            match iso {
                Some(iso) => out.push(record_from(f, &e, &iso, Some(&phi), &mut rng)),
                None => failed = true,
            }
        }
        if !failed {
            return Ok(out);
        }
        // a repeated root of Phi_l(j, Y) defeats Elkies' formulas: fall through to the
        // division polynomial, which sees every rational kernel
    }
    if ell > 31 {
        return Err("this curve needs the division-polynomial route (j = 0, 1728, a repeated root of Phi_l(j, Y), or p <= l + 1), limited to l <= 31".into());
    }
    Ok(kernel_polys(f, &e, ell, &mut rng)
        .into_iter()
        .map(|h| {
            let iso = kohel(f, &e, &h, ell);
            record_from(f, &e, &iso, None, &mut rng)
        })
        .collect())
}

/// The isogeny with the given monic kernel polynomial (coefficients low to high; degree
/// (l-1)/2 for odd l, 1 for l = 2), by Kohel's formulas.
pub fn isogeny_from_kernel(
    curve: &PrimeCurve,
    kernel: &[Int],
    ell: u64,
    seed: u64,
) -> Result<IsogenyRecord, String> {
    struct T<'a>(&'a PrimeCurve, &'a [Int], u64, u64);
    impl FieldTask for T<'_> {
        type Out = Result<IsogenyRecord, String>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            let (c, ker, ell) = (self.0, self.1, self.2);
            let e = curve_in(f, c);
            let h: Vec<F::E> = ker.iter().map(|v| f.elem(v)).collect();
            if h.last().is_none_or(|&l| l != f.one()) {
                return Err("kernel polynomial must be monic".into());
            }
            let mut rng = Rng::new(self.3);
            Ok(record_from(f, &e, &kohel(f, &e, &h, ell), None, &mut rng))
        }
    }
    with_field(&curve.p, T(curve, kernel, ell, seed))
}

/// The kernel polynomial prod (x - x([k]P)), k = 1..(l-1)/2, of the subgroup generated by a
/// point with x-coordinate `x0` (odd l); None if x0 is not the x-coordinate of a point of
/// order l (over F_p or its quadratic extension).
pub fn kernel_from_x(curve: &PrimeCurve, x0: &Int, ell: u64) -> Option<Vec<Int>> {
    struct T<'a>(&'a PrimeCurve, &'a Int, u64);
    impl FieldTask for T<'_> {
        type Out = Option<Vec<Int>>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            let e = curve_in(f, self.0);
            let n = ((self.2 - 1) / 2) as usize;
            let xs = crate::kernel::xonly::multiples_x_fast(f, &e, f.elem(self.1), n + 1)?;
            // x([n+1]P) = x([n]P) exactly when P has order 2n + 1 = l
            if xs[n] != xs[n - 1] {
                return None;
            }
            let h = poly::from_roots(f, &xs[..n]);
            Some(h.iter().map(|&c| f.int(c)).collect())
        }
    }
    if ell < 3 || ell.is_multiple_of(2) {
        return None;
    }
    with_field(&curve.p, T(curve, x0, ell))
}

// ------------------------------------------------------------------ modular polynomials

/// Phi_l mod p as the (l+2) x (l+2) matrix c[i][j] of X^i Y^j (Hecke/Newton; char > l + 1).
pub fn modular_polynomial(p: &Int, ell: usize) -> Result<Vec<Vec<Int>>, String> {
    if *p <= Int::from(ell as i64 + 1) {
        return Err("Phi_l mod p by q-expansions needs p > l + 1".into());
    }
    struct T(usize);
    impl FieldTask for T {
        type Out = Vec<Vec<Int>>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            Phi::compute(f, self.0)
                .c
                .iter()
                .map(|row| row.iter().map(|&x| f.int(x)).collect())
                .collect()
        }
    }
    Ok(with_field(p, T(ell)))
}

/// Phi_l over Z (CRT of Hecke/Newton over 62-bit primes) and a check: its reduction modulo a
/// further prime outside the CRT set equals Phi_l computed there.
pub fn modular_polynomial_integer(ell: usize) -> (Vec<Vec<Int>>, Check) {
    let c = integer_coeffs(ell);
    let q = (1u64 << 61) - 1; // a Mersenne prime; the CRT set uses primes just below 2^62
    let fq = Zp::new(q);
    let direct = Phi::compute(&fq, ell).c;
    let qi = Int::from(q);
    let same = c
        .iter()
        .zip(&direct)
        .all(|(r, d)| r.iter().zip(d).all(|(x, &y)| x.modulo(&qi) == Int::from(y)));
    let check = Check::new(
        "reduction_mod_independent_prime",
        if same { Status::Pass } else { Status::Fail },
        format!(
            "Phi_l over Z reduced mod 2^61 - 1 equals Phi_l computed mod 2^61 - 1 ({})",
            if same { "equal" } else { "differs" }
        ),
    );
    (c, check)
}

/// Coefficients of Phi_l(X, j) mod p (low to high) and its roots in F_p (the j-invariants of the
/// F_p-rational l-isogenous curves).
pub fn modular_polynomial_at(
    p: &Int,
    ell: usize,
    j: &Int,
    seed: u64,
) -> Result<(Vec<Int>, Vec<Int>), String> {
    if *p <= Int::from(ell as i64 + 1) {
        return Err("Phi_l mod p by q-expansions needs p > l + 1".into());
    }
    struct T<'a>(usize, &'a Int, u64);
    impl FieldTask for T<'_> {
        type Out = (Vec<Int>, Vec<Int>);
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            let phi = Phi::compute(f, self.0);
            let g = phi.y_poly(f, f.elem(self.1));
            let mut rng = Rng::new(self.2);
            let mut roots: Vec<Int> = poly::roots(f, &g, &mut rng)
                .into_iter()
                .map(|r| f.int(r))
                .collect();
            roots.sort();
            (g.iter().map(|&x| f.int(x)).collect(), roots)
        }
    }
    Ok(with_field(p, T(ell, j, seed)))
}
