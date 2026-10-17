//! Entry points for the isogeny families beyond counting and the Elkies/Kohel route: the
//! kernel-to-isogeny formulas (Kohel, Vélu, x-only Vélu, sqrt-Vélu), l^n chains with
//! strategies, every BMSS variant, dual isogenies and the CSIDH class-group action. Each result
//! carries checks that do not reuse the code path that produced it.
use super::*;
use crate::curve::{pmul, rhs};
use crate::find::bmss::{self, sigma_from_phi};
use crate::find::dual::dual_isogeny;
use crate::find::elkies::elkies_codomain;
use crate::kernel::chain::{ell_power_isogeny, Strategy};
use crate::kernel::montgomery::{a24, ladder_xz};
use crate::kernel::sqrt_velu::sqrt_velu_fast;
use crate::kernel::velu::{velu_cyclic, IntoX};
use crate::kernel::xonly::velu_xonly_fast;
use crate::path::csidh::Csidh;

pub use crate::find::bmss::Method as BmssMethod;

/// The random-point checks shared by every map: images land on the codomain and the map
/// respects sums.
fn map_checks<F: PrimeFieldInt, I: Isogeny<F> + ?Sized>(
    f: &F,
    iso: &I,
    rng: &mut Rng,
) -> Vec<Check> {
    let (e, cod) = (*iso.domain(), *iso.codomain());
    let trials = 8;
    let (mut on, mut hom) = (0, 0);
    for _ in 0..trials {
        let (p1, p2) = (random_point_f(f, &e, rng), random_point_f(f, &e, rng));
        let (i1, i2) = (iso.eval(f, &p1), iso.eval(f, &p2));
        if on_curve(f, &cod, &i1) && on_curve(f, &cod, &i2) {
            on += 1;
        }
        if iso.eval(f, &padd(f, &e, &p1, &p2)) == padd(f, &cod, &i1, &i2) {
            hom += 1;
        }
    }
    let st = |ok: bool| if ok { Status::Pass } else { Status::Fail };
    vec![
        Check::new(
            "image_on_codomain",
            st(on == trials),
            format!("{on}/{trials} pairs of random points map onto E'"),
        ),
        Check::new(
            "homomorphism",
            st(hom == trials),
            format!("phi(P + Q) = phi(P) + phi(Q) for {hom}/{trials} random pairs"),
        ),
    ]
}

fn curve_out<F: PrimeFieldInt>(f: &F, c: &Curve<F::E>) -> (Int, Int, Int) {
    (f.int(c.a), f.int(c.b), f.int(jinv(f, c)))
}

// ------------------------------------------------------------------ kernel to isogeny

/// The formulas that turn a kernel point into its isogeny.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VeluMethod {
    /// Kohel (1996): kernel polynomial, explicit rational maps
    Kohel,
    /// Vélu (1971) from the affine multiples of P (P must be F_p-rational)
    Velu,
    /// Vélu from x(P) with projective x-only multiples and one inversion
    XOnly,
    /// Bernstein–De Feo–Leroux–Smith (2020), O~(sqrt l) per codomain and evaluation
    SqrtVelu,
}

impl VeluMethod {
    pub const ALL: [VeluMethod; 4] = [
        VeluMethod::Kohel,
        VeluMethod::Velu,
        VeluMethod::XOnly,
        VeluMethod::SqrtVelu,
    ];
    pub fn name(self) -> &'static str {
        match self {
            VeluMethod::Kohel => "kohel",
            VeluMethod::Velu => "velu",
            VeluMethod::XOnly => "xonly",
            VeluMethod::SqrtVelu => "sqrtvelu",
        }
    }
    pub fn from_name(s: &str) -> Option<VeluMethod> {
        VeluMethod::ALL.into_iter().find(|m| m.name() == s)
    }
}

/// An isogeny given by its codomain (and checked as a map), without explicit rational maps.
#[derive(Clone, Debug)]
pub struct MapRecord {
    pub method: String,
    pub degree: Int,
    pub codomain: (Int, Int),
    pub j_codomain: Int,
    pub checks: Vec<Check>,
}

/// The odd-degree l-isogeny with kernel <P>, x(P) = x0, by `method`. Checked against Kohel's
/// codomain from the kernel polynomial built independently from the x-multiples of P.
pub fn isogeny_from_x(
    curve: &PrimeCurve,
    x0: &Int,
    ell: u64,
    method: VeluMethod,
    seed: u64,
) -> Result<MapRecord, String> {
    if ell < 3 || !crate::field::is_prime(ell) {
        return Err(format!("l = {ell} must be an odd prime"));
    }
    let ker = kernel_from_x(curve, x0, ell)
        .ok_or("x is not the x-coordinate of a point of order l (over F_p or F_p^2)")?;
    struct T<'a>(&'a PrimeCurve, &'a Int, u64, VeluMethod, u64, Vec<Int>);
    impl FieldTask for T<'_> {
        type Out = Result<MapRecord, String>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            let T(c, x0, ell, method, seed, ker) = self;
            let mut rng = Rng::new(seed);
            let e = curve_in(f, c);
            let x = f.elem(x0);
            let h: Vec<F::E> = ker.iter().map(|v| f.elem(v)).collect();
            let reference = kohel(f, &e, &h, ell);
            let (cod, mut checks): (Curve<F::E>, Vec<Check>) = match method {
                VeluMethod::Kohel => (reference.cod, map_checks(f, &reference, &mut rng)),
                VeluMethod::Velu => {
                    let y = f.sqrt(rhs(f, &e, x)).ok_or(
                        "the kernel point is not F_p-rational (y in F_p^2); use --method xonly, sqrtvelu or kohel",
                    )?;
                    let iso = velu_cyclic(f, &e, &Pt::Aff(x, y), ell);
                    (iso.cod, map_checks(f, &iso, &mut rng))
                }
                VeluMethod::XOnly => {
                    let iso = velu_xonly_fast(f, &e, x, ell).ok_or("x-only Vélu failed")?;
                    (iso.cod, map_checks(f, &iso, &mut rng))
                }
                VeluMethod::SqrtVelu => {
                    let iso = sqrt_velu_fast(f, &e, x, ell).ok_or("sqrt-Vélu failed")?;
                    (iso.cod, map_checks(f, &iso, &mut rng))
                }
            };
            if method != VeluMethod::Kohel {
                let same = cod == reference.cod;
                checks.push(Check::new(
                    "codomain_matches_kohel",
                    if same { Status::Pass } else { Status::Fail },
                    "the codomain equals Kohel's from the kernel polynomial".into(),
                ));
            }
            let (a, b, j) = curve_out(f, &cod);
            Ok(MapRecord {
                method: method.name().into(),
                degree: Int::from(ell),
                codomain: (a, b),
                j_codomain: j,
                checks,
            })
        }
    }
    with_field(&curve.p, T(curve, x0, ell, method, seed, ker))
}

// ------------------------------------------------------------------ chains

/// How an l^n chain walks its kernel: multiply-then-push choices of De Feo–Jao–Plût.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChainStrategy {
    Naive,
    Balanced,
    Optimal,
}

impl ChainStrategy {
    pub fn from_name(s: &str) -> Option<ChainStrategy> {
        match s {
            "naive" => Some(ChainStrategy::Naive),
            "balanced" => Some(ChainStrategy::Balanced),
            "optimal" => Some(ChainStrategy::Optimal),
            _ => None,
        }
    }
    fn strategy(self, ell: u64) -> Strategy {
        match self {
            ChainStrategy::Naive => Strategy::Naive,
            ChainStrategy::Balanced => Strategy::Balanced,
            // costs in field multiplications: [l] by a ladder, one Vélu evaluation and one
            // Vélu construction, each linear in l
            ChainStrategy::Optimal => {
                let l = ell as f64;
                Strategy::Optimal {
                    mul: 12.0 * l.log2().ceil(),
                    eval: 4.0 * l,
                    iso: 6.0 * l,
                }
            }
        }
    }
}

/// An l^n-isogeny computed as a chain of n l-isogenies.
#[derive(Clone, Debug)]
pub struct ChainRecord {
    pub ell: u64,
    pub n: usize,
    pub strategy: ChainStrategy,
    /// the kernel generator, of exact order l^n
    pub kernel_point: (Int, Int),
    pub codomain: (Int, Int),
    pub j_codomain: Int,
    /// j-invariants after each step (the path in the l-isogeny graph)
    pub path: Vec<Int>,
    /// multiplications by l, point evaluations and isogeny constructions
    pub stats: (usize, usize, usize),
    pub checks: Vec<Check>,
}

/// The l^n-isogeny with kernel <P>. P is given, or found from the group order (counted by SEA
/// when not given): a random point times the prime-to-l part of #E, then times l until its
/// order is l^n.
pub fn ell_power_chain(
    curve: &PrimeCurve,
    ell: u64,
    n: usize,
    kernel: Option<(Int, Int)>,
    order: Option<Int>,
    strategy: ChainStrategy,
    seed: u64,
) -> Result<ChainRecord, String> {
    if ell < 2 || !crate::field::is_prime(ell) || n == 0 {
        return Err("l must be prime and n at least 1".into());
    }
    if (ell as f64).powi(n as i32) > 2f64.powi(63) {
        return Err("l^n must stay below 2^63".into());
    }
    let order = match (&kernel, order) {
        (Some(_), o) => o,
        (None, Some(o)) => Some(o),
        (None, None) => Some(count_points(curve, seed)?.order),
    };
    struct T<'a>(
        &'a PrimeCurve,
        u64,
        usize,
        Option<(Int, Int)>,
        Option<Int>,
        ChainStrategy,
        u64,
    );
    impl FieldTask for T<'_> {
        type Out = Result<ChainRecord, String>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            let T(c, ell, n, kernel, order, strategy, seed) = self;
            let mut rng = Rng::new(seed);
            let e = curve_in(f, c);
            let ln = (ell as u128).pow(n as u32);
            let p = match kernel {
                Some((x, y)) => {
                    let p = Pt::Aff(f.elem(&x), f.elem(&y));
                    if !on_curve(f, &e, &p) {
                        return Err("the kernel point is not on the curve".into());
                    }
                    p
                }
                None => find_ell_power_point(f, &e, ell, n, &order.unwrap(), &mut rng)?,
            };
            let mut checks = vec![];
            let exact =
                pmul(f, &e, &p, ln) == Pt::Inf && pmul(f, &e, &p, ln / ell as u128) != Pt::Inf;
            checks.push(Check::new(
                "kernel_order",
                if exact { Status::Pass } else { Status::Fail },
                format!("[l^n]P = O and [l^(n-1)]P != O for l^n = {ln}"),
            ));
            if !exact {
                return Err(format!("the kernel point does not have order exactly {ln}"));
            }
            let (chain, pushed, st) =
                ell_power_isogeny(f, &e, &p, ell, n, &[p], strategy.strategy(ell));
            checks.extend(map_checks(f, &chain, &mut rng));
            checks.push(Check::new(
                "kernel_maps_to_zero",
                if pushed[0] == Pt::Inf {
                    Status::Pass
                } else {
                    Status::Fail
                },
                "phi(P) = O for the kernel generator".into(),
            ));
            // the naive strategy recomputes every kernel point from scratch: an independent
            // route to the same codomain
            if strategy != ChainStrategy::Naive {
                let (naive, _, _) = ell_power_isogeny(f, &e, &p, ell, n, &[], Strategy::Naive);
                checks.push(Check::new(
                    "matches_naive_strategy",
                    if naive.cod == chain.cod {
                        Status::Pass
                    } else {
                        Status::Fail
                    },
                    "the codomain equals the naive strategy's".into(),
                ));
            }
            let mut path = vec![f.int(jinv(f, &e))];
            path.extend(chain.steps.iter().map(|s| f.int(jinv(f, s.codomain()))));
            let (a, b, j) = curve_out(f, &chain.cod);
            let (kx, ky) = point_out(f, &p).unwrap();
            Ok(ChainRecord {
                ell,
                n,
                strategy,
                kernel_point: (kx, ky),
                codomain: (a, b),
                j_codomain: j,
                path,
                stats: (st.l_mults, st.evals, st.builds),
                checks,
            })
        }
    }
    with_field(&curve.p, T(curve, ell, n, kernel, order, strategy, seed))
}

fn find_ell_power_point<F: PrimeFieldInt>(
    f: &F,
    e: &Curve<F::E>,
    ell: u64,
    n: usize,
    order: &Int,
    rng: &mut Rng,
) -> Result<Pt<F::E>, String> {
    let l = Int::from(ell);
    let (mut cof, mut v) = (order.clone(), 0usize);
    loop {
        let (q, r) = cof.divrem_trunc(&l);
        if !r.is_zero() {
            break;
        }
        cof = q;
        v += 1;
    }
    if v < n {
        return Err(format!(
            "l^n does not divide #E = {order} (l-adic valuation {v})"
        ));
    }
    for _ in 0..64 {
        let mut q = pmul_big(f, e, &random_point_f(f, e, rng), cof.mag());
        // the order of q is l^m for some m <= v; reduce it to l^n if m >= n
        let mut m = 0;
        let mut t = q;
        while t != Pt::Inf {
            t = pmul(f, e, &t, ell as u128);
            m += 1;
        }
        if m >= n {
            for _ in n..m {
                q = pmul(f, e, &q, ell as u128);
            }
            return Ok(q);
        }
    }
    Err(format!(
        "no point of order l^n found in 64 tries: E[l^infinity](F_p) has exponent below l^{n}"
    ))
}

// ------------------------------------------------------------------ BMSS and duals

/// The rational l-isogenies of E (odd prime l > 3 where Phi_l applies, j(E) not 0 or 1728) by
/// one BMSS variant: Elkies gives each normalised codomain, `method` the kernel.
pub fn isogenies_by(
    curve: &PrimeCurve,
    ell: u64,
    method: BmssMethod,
    seed: u64,
) -> Result<Vec<IsogenyRecord>, String> {
    struct T<'a>(&'a PrimeCurve, u64, BmssMethod, u64);
    impl FieldTask for T<'_> {
        type Out = Result<Vec<IsogenyRecord>, String>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            let T(c, ell, method, seed) = self;
            let mut rng = Rng::new(seed);
            let e = curve_in(f, c);
            let (phi, isos) = elkies_isogenies(f, c, &e, ell, method, &mut rng)?;
            Ok(isos
                .iter()
                .map(|iso| record_from(f, &e, iso, Some(&phi), &mut rng))
                .collect())
        }
    }
    with_field(&curve.p, T(curve, ell, method, seed))
}

fn elkies_isogenies<F: PrimeFieldInt>(
    f: &F,
    c: &PrimeCurve,
    e: &Curve<F::E>,
    ell: u64,
    method: BmssMethod,
    rng: &mut Rng,
) -> Result<(Phi<F>, Vec<RatIsogeny<F>>), String> {
    if ell < 3 || !crate::field::is_prime(ell) {
        return Err(format!("l = {ell} must be an odd prime"));
    }
    if c.a.is_zero() || c.b.is_zero() {
        return Err("BMSS needs j(E) other than 0 and 1728".into());
    }
    if c.p <= Int::from(2 * ell + 1) {
        return Err("BMSS needs p > 2l + 1".into());
    }
    let phi = Phi::compute(f, ell as usize);
    let mut out = vec![];
    for jt in phi.neighbors(f, jinv(f, e), rng) {
        let et = elkies_codomain(f, &phi, e, jt)
            .ok_or("Elkies' codomain is undefined here (repeated root or j' in {0, 1728})")?;
        let sigma = if method.needs_sigma() {
            Some(sigma_from_phi(f, &phi, e, &et).ok_or("sigma from Phi_l is undefined here")?)
        } else {
            None
        };
        out.push(
            bmss::isogeny(f, method, e, &et, ell as usize, sigma)
                .ok_or(format!("{} found no kernel", method.name()))?,
        );
    }
    Ok((phi, out))
}

/// An l-isogeny together with its dual.
#[derive(Clone, Debug)]
pub struct DualRecord {
    pub isogeny: IsogenyRecord,
    /// kernel polynomial of the dual (on the codomain E'), low to high
    pub dual_kernel: Vec<Int>,
    /// x-map of the dual E' -> E, num / den, low to high
    pub dual_map_num: Vec<Int>,
    pub dual_map_den: Vec<Int>,
    pub checks: Vec<Check>,
}

/// Every rational l-isogeny phi of E (as in [`isogenies_by`] with fast_elkies_prime) and its
/// dual from Phi_l alone, checked by dual(phi(P)) = [l]P.
pub fn duals(curve: &PrimeCurve, ell: u64, seed: u64) -> Result<Vec<DualRecord>, String> {
    struct T<'a>(&'a PrimeCurve, u64, u64);
    impl FieldTask for T<'_> {
        type Out = Result<Vec<DualRecord>, String>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            let T(c, ell, seed) = self;
            let mut rng = Rng::new(seed);
            let e = curve_in(f, c);
            let (phi, isos) =
                elkies_isogenies(f, c, &e, ell, BmssMethod::FastElkiesPrime, &mut rng)?;
            let mut out = vec![];
            for iso in &isos {
                let rec = record_from(f, &e, iso, Some(&phi), &mut rng);
                let d = dual_isogeny(f, &phi, iso).ok_or("the dual is undefined here")?;
                let mut checks = map_checks(f, &d, &mut rng);
                checks.push(Check::new(
                    "codomain_is_domain",
                    if d.cod == e {
                        Status::Pass
                    } else {
                        Status::Fail
                    },
                    "the dual maps E' back onto E".into(),
                ));
                // u = +-l fixes the dual up to [-1]: compare x-coordinates
                let trials = 8;
                let ok = (0..trials)
                    .filter(|_| {
                        let p = random_point_f(f, &e, &mut rng);
                        d.eval(f, &iso.eval(f, &p)).into_x()
                            == pmul(f, &e, &p, ell as u128).into_x()
                    })
                    .count();
                checks.push(Check::new(
                    "dual_composition",
                    if ok == trials {
                        Status::Pass
                    } else {
                        Status::Fail
                    },
                    format!("x(dual(phi(P))) = x([l]P) for {ok}/{trials} random points"),
                ));
                let u2 = f.mul(d.u, d.u);
                let ints = |v: &[F::E]| v.iter().map(|&x| f.int(x)).collect::<Vec<Int>>();
                let den: Vec<F::E> = d.inner.den.iter().map(|&x| f.mul(x, u2)).collect();
                out.push(DualRecord {
                    isogeny: rec,
                    dual_kernel: ints(&d.inner.ker),
                    dual_map_num: ints(&d.inner.num),
                    dual_map_den: ints(&den),
                    checks,
                });
            }
            Ok(out)
        }
    }
    with_field(&curve.p, T(curve, ell, seed))
}

// ------------------------------------------------------------------ CSIDH

/// The result of the CSIDH class-group action on a Montgomery curve y^2 = x^3 + A x^2 + x.
#[derive(Clone, Debug)]
pub struct CsidhRecord {
    pub p: Int,
    pub primes: Vec<u64>,
    pub a_in: Int,
    pub a_out: Int,
    pub checks: Vec<Check>,
}

/// [prod l_i^{e_i}] E_A over p = 4 prod(l_i) - 1, by the tree-strategy action, checked against
/// the one-isogeny-at-a-time action and for supersingularity of the result.
pub fn csidh_action(
    primes: &[u64],
    a: &Int,
    exps: &[i32],
    seed: u64,
) -> Result<CsidhRecord, String> {
    if exps.len() != primes.len() {
        return Err(format!(
            "{} exponents for {} primes",
            exps.len(),
            primes.len()
        ));
    }
    if primes.iter().any(|&l| l < 3 || !crate::field::is_prime(l)) {
        return Err("every l_i must be an odd prime".into());
    }
    let mut p = Int::from(4u64);
    for &l in primes {
        p = &p * &Int::from(l);
    }
    let p = &p - &Int::one();
    if !p.is_probable_prime() {
        return Err(format!("4 prod(l_i) - 1 = {p} is not prime"));
    }
    if p.bits() > 640 {
        return Err("p above 640 bits is not supported".into());
    }
    struct T<'a>(&'a [u64], &'a Int, &'a [i32], u64, Int);
    impl FieldTask for T<'_> {
        type Out = Result<CsidhRecord, String>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            let T(primes, a, exps, seed, p) = self;
            let cs = Csidh {
                fp: f.clone(),
                primes: primes.to_vec(),
                p1: (&p + &Int::one()).mag().clone(),
                sqrt_velu_min: crate::path::csidh::SQRT_VELU_MIN,
            };
            let ai = f.elem(a);
            let mut rng = Rng::new(seed);
            let mut checks = vec![supersingular(f, &cs, ai, &mut rng, "input_supersingular")];
            if checks[0].status != Status::Pass {
                return Err("the input curve is not supersingular".into());
            }
            let out = cs.action_fast(ai, exps, &mut rng);
            let slow = cs.action(ai, exps, &mut rng);
            checks.push(Check::new(
                "matches_one_at_a_time",
                if out == slow {
                    Status::Pass
                } else {
                    Status::Fail
                },
                "the tree-strategy action equals the one-isogeny-at-a-time action".into(),
            ));
            checks.push(supersingular(f, &cs, out, &mut rng, "output_supersingular"));
            Ok(CsidhRecord {
                p: p.clone(),
                primes: primes.to_vec(),
                a_in: a.modulo(&p),
                a_out: f.int(out),
                checks,
            })
        }
    }
    let p2 = p.clone();
    with_field(&p, T(primes, a, exps, seed, p2))
}

/// CSIDH-512's primes: the first 73 odd primes and 587.
pub fn csidh512_primes() -> Vec<u64> {
    let mut primes = vec![];
    let mut l = 3u64;
    while primes.len() < 73 {
        if crate::field::is_prime(l) {
            primes.push(l);
        }
        l += 2;
    }
    primes.push(587);
    primes
}

/// [p + 1] x = O for random x on the curve or its twist (both have p + 1 points when E_A is
/// supersingular), and some x of order above 4 sqrt(p), which rules out every ordinary curve
/// (CSIDH's validation).
fn supersingular<F: PrimeFieldInt>(
    f: &F,
    cs: &Csidh<F>,
    a: F::E,
    rng: &mut Rng,
    name: &str,
) -> Check {
    let e24 = a24(f, a);
    let two = f.from_u64(2);
    if a == two || a == f.neg(two) {
        return Check::new(name, Status::Fail, "A = +-2 is singular".into());
    }
    let p1 = &cs.p1;
    let bound = Int::from_big(p1).bits() / 2 + 2; // log2(4 sqrt(p)) <= this
    for _ in 0..32 {
        let x = f.random(rng);
        let pt = (x, f.one());
        if !f.is_zero(ladder_xz(f, e24, pt, p1).1) {
            return Check::new(
                name,
                Status::Fail,
                "a point whose order does not divide p + 1".into(),
            );
        }
        // order of x: divide out each l_i while [(p+1)/l_i] x = O
        let mut ord = Int::one();
        for &l in &cs.primes {
            let k = p1.divrem_small(l).0;
            if !f.is_zero(ladder_xz(f, e24, pt, &k).1) {
                ord = &ord * &Int::from(l);
            }
        }
        if ord.bits() > bound {
            return Check::new(
                name,
                Status::Pass,
                format!(
                    "a point of order dividing p + 1 and above 4 sqrt(p) ({} bits)",
                    ord.bits()
                ),
            );
        }
    }
    Check::new(
        name,
        Status::Indeterminate,
        "every sampled point had order dividing p + 1 but none above 4 sqrt(p)".into(),
    )
}

// ------------------------------------------------------------------ path finding

/// Path finding in the l-isogeny graphs of ordinary curves on the j-line.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathMethod {
    /// Galbraith (1999): bidirectional breadth-first search
    Galbraith,
    /// Galbraith–Hess–Smart (2002): two random walks until they meet
    Ghs,
    /// Galbraith–Stolbunov (2013): GHS walks with small degrees favoured
    GalbraithStolbunov,
}

impl PathMethod {
    pub fn from_name(s: &str) -> Option<PathMethod> {
        match s {
            "galbraith" => Some(PathMethod::Galbraith),
            "ghs" => Some(PathMethod::Ghs),
            "gs" | "galbraith_stolbunov" => Some(PathMethod::GalbraithStolbunov),
            _ => None,
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            PathMethod::Galbraith => "galbraith",
            PathMethod::Ghs => "ghs",
            PathMethod::GalbraithStolbunov => "galbraith_stolbunov",
        }
    }
}

/// A path j1 -> j2 of isogenies of the given degrees.
#[derive(Clone, Debug)]
pub struct PathRecord {
    pub method: PathMethod,
    pub js: Vec<Int>,
    pub ells: Vec<usize>,
    /// nodes expanded (BFS) or walk steps (GHS, Galbraith–Stolbunov)
    pub work: usize,
    pub checks: Vec<Check>,
}

/// A curve with the given j-invariant: y^2 = x^3 + 3k x + 2k, k = j / (1728 - j).
fn curve_with_j<F: PrimeFieldInt>(f: &F, j: F::E) -> Curve<F::E> {
    let c1728 = f.from_u64(1728);
    if f.is_zero(j) {
        return Curve::new(f.zero(), f.one());
    }
    if j == c1728 {
        return Curve::new(f.one(), f.zero());
    }
    let k = f.div(j, f.sub(c1728, j));
    Curve::new(f.mul(f.from_u64(3), k), f.mul(f.from_u64(2), k))
}

/// A path of l-isogenies (l in `ells`) from j1 to j2 over F_p, or an error when the budget runs
/// out. Checked edge by edge on Phi_l, and by building the isogenies explicitly (Elkies + BMSS)
/// from a curve with j-invariant j1 and arriving at j2.
pub fn isogeny_path(
    p: &Int,
    j1: &Int,
    j2: &Int,
    ells: &[usize],
    method: PathMethod,
    max_steps: usize,
    seed: u64,
) -> Result<PathRecord, String> {
    if *p <= Int::from(3i64) || !p.is_probable_prime() || p.bits() > 640 {
        return Err("p must be a prime above 3 (at most 640 bits)".into());
    }
    if ells.is_empty() || ells.iter().any(|&l| !crate::field::is_prime(l as u64)) {
        return Err("--ells must list primes".into());
    }
    struct T<'a>(&'a Int, &'a Int, &'a [usize], PathMethod, usize, u64);
    impl FieldTask for T<'_> {
        type Out = Result<PathRecord, String>;
        fn run<F: PrimeFieldInt>(self, f: &F) -> Self::Out {
            use crate::path::graph::{verify_path, PhiCache};
            let T(j1, j2, ells, method, max_steps, seed) = self;
            let mut rng = Rng::new(seed);
            let cache = PhiCache::new(f, ells);
            let (a, b) = (f.elem(j1), f.elem(j2));
            let (path, work) = match method {
                PathMethod::Galbraith => {
                    let (p, st) = crate::path::galbraith::galbraith(
                        f, &cache, ells, a, b, max_steps, &mut rng,
                    );
                    (p, st.nodes_expanded)
                }
                PathMethod::Ghs => {
                    let (p, st) = crate::path::ghs::ghs(f, &cache, ells, a, b, max_steps, &mut rng);
                    (p, st.steps)
                }
                PathMethod::GalbraithStolbunov => {
                    let w: Vec<u32> = ells
                        .iter()
                        .map(|&l| (1_000_000 / (l * l)).max(1) as u32)
                        .collect();
                    let (p, st) = crate::path::ghs::galbraith_stolbunov(
                        f, &cache, ells, &w, a, b, max_steps, &mut rng,
                    );
                    (p, st.steps)
                }
            };
            let path = without_cycles(path.ok_or(format!(
                "no path found within {max_steps} steps (the curves may not be isogenous)"
            ))?);
            let mut checks = vec![];
            let ends = path.js.first() == Some(&a) && path.js.last() == Some(&b);
            checks.push(Check::new(
                "endpoints",
                if ends { Status::Pass } else { Status::Fail },
                "the path starts at j1 and ends at j2".into(),
            ));
            checks.push(Check::new(
                "modular_polynomial_edges",
                if verify_path(f, &cache, &path) {
                    Status::Pass
                } else {
                    Status::Fail
                },
                format!("Phi_l(j_i, j_(i+1)) = 0 on all {} edges", path.ells.len()),
            ));
            let e = curve_with_j(f, a);
            checks.push(match explicit_path(f, &cache, &e, &path, &mut rng) {
                Some(end) => Check::new(
                    "explicit_isogenies",
                    if end == b { Status::Pass } else { Status::Fail },
                    "every edge built as an isogeny of its degree (Elkies + BMSS, or Kohel on a division-polynomial kernel), ending at j2".into(),
                ),
                None => Check::new(
                    "explicit_isogenies",
                    Status::Indeterminate,
                    "an edge through j = 0 or 1728 of degree above 31 could not be built".into(),
                ),
            });
            Ok(PathRecord {
                method,
                js: path.js.iter().map(|&j| f.int(j)).collect(),
                ells: path.ells.clone(),
                work,
                checks,
            })
        }
    }
    with_field(p, T(j1, j2, ells, method, max_steps, seed))
}

/// The path with every closed sub-walk cut out (random walks revisit vertices).
fn without_cycles<E: Copy + Eq + std::hash::Hash>(
    path: crate::path::graph::Path<E>,
) -> crate::path::graph::Path<E> {
    let (mut js, mut ells) = (vec![path.js[0]], vec![]);
    let mut at: HashMap<E, usize> = HashMap::from([(path.js[0], 0)]);
    for (i, &j) in path.js[1..].iter().enumerate() {
        if let Some(&k) = at.get(&j) {
            for gone in js.drain(k + 1..) {
                at.remove(&gone);
            }
            ells.truncate(k);
        } else {
            at.insert(j, js.len());
            js.push(j);
            ells.push(path.ells[i]);
        }
    }
    crate::path::graph::Path { js, ells }
}

/// Walk `path` with explicit isogenies from `e`; the j-invariant reached, or None when an edge
/// cannot be built. Each edge uses Elkies + BMSS, or Kohel on a factor of the division polynomial
/// where Elkies' codomain is undefined (j in {0, 1728}, repeated roots of Phi_l).
fn explicit_path<F: PrimeFieldInt>(
    f: &F,
    cache: &crate::path::graph::PhiCache<F>,
    e: &Curve<F::E>,
    path: &crate::path::graph::Path<F::E>,
    rng: &mut Rng,
) -> Option<F::E> {
    let mut cur = *e;
    for (i, &l) in path.ells.iter().enumerate() {
        let target = path.js[i + 1];
        let fast = elkies_codomain(f, cache.get(l), &cur, target)
            .and_then(|et| crate::find::elkies::normalised_isogeny(f, &cur, &et, l));
        let iso = match fast {
            Some(iso) if iso.deg == l as u64 && jinv(f, &iso.cod) == target => iso,
            _ => {
                if l > 31 {
                    return None;
                }
                let kernels = if l == 2 {
                    let cubic = vec![cur.b, cur.a, f.zero(), f.one()];
                    poly::roots(f, &cubic, rng)
                        .into_iter()
                        .map(|x0| vec![f.neg(x0), f.one()])
                        .collect()
                } else {
                    kernel_polys(f, &cur, l as u64, rng)
                };
                kernels
                    .into_iter()
                    .map(|h| kohel(f, &cur, &h, l as u64))
                    .find(|iso| jinv(f, &iso.cod) == target)?
            }
        };
        cur = iso.cod;
    }
    Some(jinv(f, &cur))
}

// ------------------------------------------------------------------ KLPT

/// A quaternion (c0 + c1 i + c2 j + c3 k) / d in B_{p,inf} (i^2 = -1, j^2 = -p).
pub type QuatInts = ([Int; 4], Int);

/// An ideal I of norm N and its KLPT equivalent J = I xi of norm l^e.
#[derive(Clone, Debug)]
pub struct KlptRecord {
    pub p: Int,
    pub norm_in: Int,
    /// Z-basis of the input left O_0-ideal I
    pub ideal_in: Vec<QuatInts>,
    /// Z-basis of the output ideal J
    pub ideal_out: Vec<QuatInts>,
    pub e: u32,
    pub xi: QuatInts,
    pub attempts: usize,
    pub checks: Vec<Check>,
}

fn quat_ints(q: &crate::quat::Quat) -> QuatInts {
    (q.c.clone(), q.d.clone())
}

/// KLPT (Kohel–Lauter–Petit–Tignol 2014) on a random left O_0-ideal of norm `n` in the
/// quaternion algebra ramified at p and infinity (p = 3 mod 4, O_0 = <1, i, (i + j)/2, (1 + k)/2>):
/// an equivalent ideal of l-power norm. Checked by recomputing the norm, inclusion in O_0, the
/// left order, J = I xi and primitivity.
pub fn klpt_random(p: &Int, n: &Int, ell: u64, seed: u64) -> Result<KlptRecord, String> {
    use crate::quat::klpt::{ideal_from, klpt};
    use crate::quat::{Alg, Quat};
    if !p.is_probable_prime() || p.mod_u64(4) != 3 {
        return Err("p must be a prime with p = 3 mod 4".into());
    }
    if *n <= Int::one() || ell < 2 || !crate::field::is_prime(ell) {
        return Err("--norm must exceed 1 and --ell must be prime".into());
    }
    if n.modulo(&Int::from(ell)).is_zero() {
        return Err("--norm must be coprime to l".into());
    }
    let alg = Alg::new(p);
    let o0 = alg.o0();
    let mut rng = Rng::new(seed);
    // a random left O_0-ideal of norm n: O_0 n + O_0 alpha with n | Nrd(alpha)
    let mut ideal = None;
    for _ in 0..100_000 {
        let x = if n.is_probable_prime() {
            let (b, c, d) = (
                Int::from(rng.next() >> 2),
                Int::from(rng.next() >> 2),
                Int::from(rng.next() >> 2),
            );
            let t = -&(&(&b * &b) + &(&alg.p * &(&(&c * &c) + &(&d * &d))));
            let Some(a) = Int::sqrt_mod_prime(&t, n) else {
                continue;
            };
            Quat::new([a, b, c, d], Int::one())
        } else {
            let c = [0; 4].map(|_| Int::from((rng.next() % 2001) as i64 - 1000));
            let x = Quat::new(c, Int::one());
            if x.is_zero() || !alg.nrd(&x).0.modulo(n).is_zero() {
                continue;
            }
            x
        };
        let i = ideal_from(&alg, &o0, n, &x);
        if i.norm(&o0) == (n.clone(), Int::one()) {
            ideal = Some(i);
            break;
        }
    }
    let i = ideal.ok_or("no ideal of this norm found by sampling (try a prime norm)")?;
    let res = klpt(&alg, &i, ell, &mut rng).ok_or("KLPT found no equivalent ideal")?;
    let l = Int::from(ell);
    let st = |ok: bool| if ok { Status::Pass } else { Status::Fail };
    let checks = vec![
        Check::new(
            "norm",
            st(res.j.norm(&o0) == (l.pow(res.e), Int::one())),
            format!("N(J) = {ell}^{}", res.e),
        ),
        Check::new(
            "integral",
            st(o0.contains_lattice(&res.j)),
            "J is contained in O_0".into(),
        ),
        Check::new(
            "left_order",
            st(res.j.left_order(&alg, &o0) == o0),
            "J is a left O_0-ideal".into(),
        ),
        Check::new(
            "equivalent",
            st(i.rmul(&alg, &res.xi) == res.j),
            "J = I xi".into(),
        ),
        Check::new(
            "primitive",
            st(!o0.scale(&l, &Int::one()).contains_lattice(&res.j)),
            "J is not contained in l O_0".into(),
        ),
    ];
    Ok(KlptRecord {
        p: p.clone(),
        norm_in: n.clone(),
        ideal_in: i.basis().iter().map(quat_ints).collect(),
        ideal_out: res.j.basis().iter().map(quat_ints).collect(),
        e: res.e,
        xi: quat_ints(&res.xi),
        attempts: res.attempts,
        checks,
    })
}

// ------------------------------------------------------------------ Richelot

/// The codomain of one (2,2)-isogeny out of the Jacobian of y^2 = lc prod (x - r_i).
#[derive(Clone, Debug)]
pub struct RichelotRecord {
    /// the pairing of the six roots into the quadratic factors G1, G2, G3 (indices into roots)
    pub matching: [(usize, usize); 3],
    /// Some(sextic, low to high) for a Jacobian codomain y^2 = Delta H1 H2 H3
    pub jacobian: Option<Vec<Int>>,
    /// Some((cubic1, cubic2)) when the codomain splits as E1 x E2 (Delta = 0)
    pub product: Option<(Vec<Int>, Vec<Int>)>,
    pub checks: Vec<Check>,
}

/// Every Richelot (2,2)-isogeny out of J(C), C: y^2 = lc prod_{i<6} (x - r_i) over F_p
/// (p < 2^62), one per matching of the roots. Checked by L-polynomials for p < 512 (isogenous
/// surfaces have equal L-polynomials), and for every p by mapping random points of C through the
/// Richelot correspondence onto the codomain.
pub fn richelot_isogenies(
    p: &Int,
    roots: &[Int],
    lc: &Int,
    seed: u64,
) -> Result<Vec<RichelotRecord>, String> {
    use crate::field::Zp2;
    use crate::genus2::{elliptic_trace, lpoly, matchings6, richelot, split, LPoly};
    if p.bits() > 62 || *p <= Int::from(3i64) || !p.is_probable_prime() {
        return Err("Richelot needs a prime 3 < p < 2^62".into());
    }
    if roots.len() != 6 {
        return Err("--roots needs six roots".into());
    }
    let pu = p.to_i128().unwrap() as u64;
    let fp = Zp::new(pu);
    let r: Vec<u64> = roots.iter().map(|v| fp.elem(v)).collect();
    for i in 0..6 {
        if r[..i].contains(&r[i]) {
            return Err("the six roots must be distinct mod p".into());
        }
    }
    let lcv = fp.elem(lc);
    if lcv == 0 {
        return Err("--lc must be nonzero mod p".into());
    }
    let f = poly::scale(&fp, &poly::from_roots(&fp, &r), lcv);
    let f2 = Zp2::new(pu);
    let small = pu < 512;
    let l0 = small.then(|| lpoly(&fp, &f));
    let mut rng = Rng::new(seed);
    let ints = |v: &[u64]| v.iter().map(|&c| Int::from(c)).collect::<Vec<Int>>();
    let mut out = vec![];
    for (a, b, c) in matchings6() {
        let g = [
            poly::scale(&fp, &poly::from_roots(&fp, &[r[a.0], r[a.1]]), lcv),
            poly::from_roots(&fp, &[r[b.0], r[b.1]]),
            poly::from_roots(&fp, &[r[c.0], r[c.1]]),
        ];
        let rl = richelot(&fp, &g);
        let mut checks = vec![];
        let (jac, prod) = if rl.delta != 0 {
            let fc = rl.codomain(&fp);
            if let Some(l0) = l0 {
                checks.push(Check::new(
                    "l_polynomial",
                    if lpoly(&fp, &fc) == l0 {
                        Status::Pass
                    } else {
                        Status::Fail
                    },
                    "L(J(C')) = L(J(C)) by point counts over F_p and F_p^2".into(),
                ));
            }
            let cod: Vec<(u64, u64)> = fc.iter().map(|&c| (c, 0)).collect();
            let (mut on, mut tried) = (0, 0);
            while tried < 16 {
                let x = fp.random(&mut rng);
                let Some(y) = fp.sqrt(poly::eval(&fp, &f, x)) else {
                    continue;
                };
                if y == 0 {
                    continue;
                }
                tried += 1;
                let img = rl.image_point(&fp, &g, &f2, |c| (c, 0), x, y, &mut rng);
                if !img.is_empty()
                    && img
                        .iter()
                        .all(|&(xp, yp)| f2.mul(yp, yp) == poly::eval(&f2, &cod, xp))
                {
                    on += 1;
                }
            }
            checks.push(Check::new(
                "image_on_codomain",
                if on == tried {
                    Status::Pass
                } else {
                    Status::Fail
                },
                format!("{on}/{tried} random points of C map onto C' through the correspondence"),
            ));
            (Some(ints(&fc)), None)
        } else {
            match split(&fp, &f, &g) {
                Some((e1, e2)) => {
                    checks.push(match l0 {
                        Some(l0) => {
                            let lp = LPoly::product(
                                elliptic_trace(&fp, &e1),
                                elliptic_trace(&fp, &e2),
                                pu as i64,
                            );
                            Check::new(
                                "l_polynomial",
                                if lp == l0 { Status::Pass } else { Status::Fail },
                                "L(E1) L(E2) = L(J(C)) by point counts".into(),
                            )
                        }
                        None => Check::new(
                            "l_polynomial",
                            Status::Indeterminate,
                            "split codomain: the L-polynomial check runs only for p < 512".into(),
                        ),
                    });
                    (None, Some((ints(&e1), ints(&e2))))
                }
                None => {
                    checks.push(Check::new(
                        "split_rational",
                        Status::Indeterminate,
                        "Delta = 0 but the splitting involution is not defined over F_p".into(),
                    ));
                    (None, None)
                }
            }
        };
        out.push(RichelotRecord {
            matching: [a, b, c],
            jacobian: jac,
            product: prod,
            checks,
        });
    }
    Ok(out)
}
