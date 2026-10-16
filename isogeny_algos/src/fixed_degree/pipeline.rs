//! The end-to-end fixed-degree algorithm of the slides ("Our strategy", steps 1–6):
//!
//! 1. End(E_1), End(E_2) as left O_0-ideals I_1, I_2 (`endring::end_ring`);
//! 2. a connecting ideal I between O_1 = O_R(I_1) and O_2 = O_R(I_2) (Kirschmer–Voight,
//!    `connecting::connecting_ideal`);
//! 3. the normalised norm form Q of I on an LLL-reduced basis (`connecting::norm_form`);
//! 4. x with Q(x) = 2d (Q = 2 Nrd / N(I)) by Cornacchia with guessing (`cornacchia`), Coron's
//!    bivariate or the trivariate integer Coppersmith method (`coppersmith`);
//! 5. alpha = sum x_i sigma_i in I with Nrd(alpha) = d N(I) and J = I conj(alpha) / N(I), an ideal
//!    of norm d with left order O_1 and right order conjugate to O_2;
//! 6. the isogeny of J (`represent::ideal_to_isogeny`: kernel / Vélu chain for smooth d, torsion
//!    images certified by Kani for any odd d < 2^a).
use super::connecting::{combine, connecting_ideal, connects, norm_form, NormForm};
use super::coppersmith::{solve_norm_eq_coppersmith, NormEqMode};
use super::cornacchia::{solve_norm_eq_cornacchia, CornacchiaOutcome, CornacchiaStats};
use super::endring::{end_ring, EndRing, JGraph, J};
use super::represent::{ideal_to_isogeny, IsogenyRep, ReprOptions, E4};
use crate::curve::Curve;
use crate::field::Rng;
use crate::int::Int;
use crate::quat::deuring::Deuring;
use crate::quat::{Alg, Lattice, Quat};
use std::time::Instant;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Solver {
    Cornacchia,
    CoppersmithBivariate,
    CoppersmithTrivariate,
}

#[derive(Clone, Debug)]
pub struct DegreeSolution {
    pub x: [Int; 4],
    pub alpha: Quat,
    /// J = I conj(alpha) / N(I), norm d
    pub ideal: Lattice,
    pub cornacchia_stats: Option<CornacchiaStats>,
}

/// A box bound for the coordinates of solutions of Q(x) = t on an LLL-reduced form:
/// |x_i| <= 2 sqrt(t / min_i G_ii) + 1 (Q(x) >= c min G_ii |x|_inf^2 on reduced bases).
pub fn default_bound(nf: &NormForm, t: &Int) -> Int {
    let gmin = (0..4).map(|i| nf.gram[i][i].clone()).min().unwrap();
    &(&Int::from(2i64) * &t.div_floor(&gmin).isqrt()) + &Int::one()
}

/// Steps 3–5 on a given connecting ideal: an equivalent ideal of norm d.
pub fn solve_degree(
    alg: &Alg,
    o1: &Lattice,
    conn: &Lattice,
    d: &Int,
    solver: Solver,
    bound: Option<Int>,
) -> Option<DegreeSolution> {
    let nf = norm_form(alg, o1, conn);
    let t = d * &Int::from(2i64);
    let bound = bound.unwrap_or_else(|| default_bound(&nf, &t));
    let (x, stats) = match solver {
        Solver::Cornacchia => match solve_norm_eq_cornacchia(&nf.gram, &t, &bound) {
            CornacchiaOutcome::Found(x, s) => (x, Some(s)),
            CornacchiaOutcome::NotFound(_) => return None,
        },
        Solver::CoppersmithBivariate => (
            solve_norm_eq_coppersmith(&nf.gram, &t, &bound, NormEqMode::Bivariate)?,
            None,
        ),
        Solver::CoppersmithTrivariate => (
            solve_norm_eq_coppersmith(&nf.gram, &t, &bound, NormEqMode::Trivariate)?,
            None,
        ),
    };
    let alpha = combine(&nf.basis, &x);
    let want = d * &nf.norm;
    assert!(alg.nrd(&alpha) == (want, Int::one()), "norm equation solution is wrong");
    let ideal = conn.rmul(alg, &alpha.conj().scale(&Int::one(), &nf.norm));
    assert!(ideal.norm(o1) == (d.clone(), Int::one()));
    Some(DegreeSolution {
        x,
        alpha,
        ideal,
        cornacchia_stats: stats,
    })
}

#[derive(Clone, Debug)]
pub struct Instance {
    pub j1: J,
    pub j2: J,
    /// optional curve models (over F_{p^2}, embedded in F_{p^4}) to express the answer on
    pub e1: Option<Curve<E4>>,
    pub e2: Option<Curve<E4>>,
    pub d: Int,
}

#[derive(Clone, Debug, Default)]
pub struct StepTimes {
    pub end_ring_ms: f64,
    pub connecting_ms: f64,
    pub solve_ms: f64,
    pub isogeny_ms: f64,
}

#[derive(Clone, Debug)]
pub struct Solution {
    pub end1: EndRing,
    pub end2: EndRing,
    pub connecting: Lattice,
    pub connecting_norm: Int,
    pub degree: DegreeSolution,
    pub isogeny: IsogenyRep,
    pub times: StepTimes,
}

fn ms(t: Instant) -> f64 {
    t.elapsed().as_secs_f64() * 1e3
}

/// Steps 1–6.
pub fn find_isogeny_of_degree(
    dg: &Deuring,
    graph: &mut JGraph,
    inst: &Instance,
    solver: Solver,
    opts: &ReprOptions,
    rng: &mut Rng,
) -> Result<Solution, String> {
    let alg = &dg.alg;
    let o0 = &dg.o0;
    let mut times = StepTimes::default();
    let t = Instant::now();
    let end1 = end_ring(dg, graph, inst.j1, rng).ok_or("End(E1) failed")?;
    let end2 = end_ring(dg, graph, inst.j2, rng).ok_or("End(E2) failed")?;
    times.end_ring_ms = ms(t);
    let t = Instant::now();
    let o1 = end1.ideal.right_order(alg, o0);
    let o2 = end2.ideal.right_order(alg, o0);
    let conn = connecting_ideal(alg, &o1, &o2);
    if !connects(alg, &conn, &o1, &o2) {
        return Err("connecting ideal has the wrong orders".into());
    }
    let connecting_norm = conn.norm_int(&o1);
    times.connecting_ms = ms(t);
    let t = Instant::now();
    let degree = solve_degree(alg, &o1, &conn, &inst.d, solver, None)
        .ok_or("no solution of the norm equation in the search box")?;
    times.solve_ms = ms(t);
    let t = Instant::now();
    let targets = match (&inst.e1, &inst.e2) {
        (Some(a), Some(b)) => Some((a, b)),
        _ => None,
    };
    let isogeny = ideal_to_isogeny(dg, &end1.ideal, &degree.ideal, targets, opts, rng)
        .ok_or("ideal to isogeny failed")?;
    times.isogeny_ms = ms(t);
    Ok(Solution {
        end1,
        end2,
        connecting: conn,
        connecting_norm,
        degree,
        isogeny,
        times,
    })
}

// ------------------------------------------------------------------ instance generation

/// sqrt(a) mod an odd prime q.
fn sqrt_mod(a: &Int, q: &Int) -> Option<Int> {
    if a.modulo(q).is_zero() {
        return Some(Int::zero());
    }
    Int::sqrt_mod_prime(&a.modulo(q), q)
}

fn rand_int(rng: &mut Rng, q: &Int) -> Int {
    let mut v = Int::zero();
    for _ in 0..q.bits().div_ceil(62) + 1 {
        v = &(&v * &Int::from(1i64 << 62)) + &Int::from((rng.next() >> 2) as i64);
    }
    v.modulo(q)
}

/// A random left ideal of the maximal order `o` of prime norm q (q odd): o q + o beta with
/// beta = sum c_t b_t, (c0, c1, c2) random and c3 solving Nrd(beta) = 0 mod q.
pub fn random_prime_norm_ideal(alg: &Alg, o: &Lattice, q: &Int, rng: &mut Rng) -> Lattice {
    let b = o.basis();
    // Nrd(sum c_t b_t) * den^2 is a quadratic form in c with Gram bil(b_r, b_s) (integer rows)
    let rows = o.m.clone();
    let g = |r: usize, s: usize| alg.bil(&rows[r], &rows[s]);
    loop {
        let c: Vec<Int> = (0..3).map(|_| rand_int(rng, q)).collect();
        // f(t) = A t^2 + B t + C with t = c3
        let a2 = g(3, 3);
        let mut bb = Int::zero();
        let mut cc = Int::zero();
        for r in 0..3 {
            bb = &bb + &(&(&Int::from(2i64) * &c[r]) * &g(r, 3));
            for s in 0..3 {
                cc = &cc + &(&(&c[r] * &c[s]) * &g(r, s));
            }
        }
        let a2m = a2.modulo(q);
        if a2m.is_zero() {
            continue;
        }
        let disc = (&(&bb * &bb) - &(&(&Int::from(4i64) * &a2) * &cc)).modulo(q);
        let Some(s) = sqrt_mod(&disc, q) else {
            continue;
        };
        let inv2a = (&Int::from(2i64) * &a2m).inv_mod(q).unwrap();
        let c3 = (&(&(-&bb) + &s) * &inv2a).modulo(q);
        let mut beta = Quat::zero();
        for (k, cv) in c.iter().chain(std::iter::once(&c3)).enumerate() {
            beta = beta.add(&b[k].scale(cv, &Int::one()));
        }
        let mut gens: Vec<Quat> = b.iter().map(|x| x.scale(q, &Int::one())).collect();
        gens.extend(o.rmul(alg, &beta).basis());
        let i = Lattice::from_gens(&gens);
        if i.norm(o) == (q.clone(), Int::one()) {
            return i;
        }
    }
}

/// A random left O_R(I_1)-ideal of norm d = prod q_i (product of prime-norm ideals, each in the
/// right order of the previous), i.e. a random degree-d isogeny out of E_1 on the quaternion
/// side; the planted E_2 is its Deuring image.
pub fn plant_ideal(
    alg: &Alg,
    o0: &Lattice,
    i1: &Lattice,
    primes: &[Int],
    rng: &mut Rng,
) -> Lattice {
    let o1 = i1.right_order(alg, o0);
    let mut jd = o1.clone();
    for q in primes {
        let or = jd.right_order(alg, &o1);
        let step = random_prime_norm_ideal(alg, &or, q, rng);
        jd = jd.mul(alg, &step);
    }
    jd
}
