//! The Cornacchia step of the fixed-degree algorithm (Benčina–Kutas–Merz–Petit–Stopar–Weitkämper
//! 2023; slides "Solving Step 4: Cornacchia's algorithm I/II").
//!
//! The norm equation is Q(x) = d for a positive-definite quaternary form
//! Q(x) = sum_{i,j} G[i][j] x_i x_j (G symmetric, LLL-reduced in practice). Guess (x3, x4), which
//! leaves a binary quadratic equation f(x1, x2) = Q(x1, x2, x3, x4) - d = 0 with positive-definite
//! quadratic part. Completing the square reduces it to u^2 + |D| v^2 = N, solved by Cornacchia
//! (all roots of -|D| mod N/g^2 for every g^2 | N, so N must be factored); the guess is abandoned
//! when N has too many prime factors (the slides' heuristic omega(N) > B log log N, B = 11) or when
//! factoring N exceeds a work bound.
//!
//! Exactness: every returned solution is verified by evaluation. [`cornacchia_all`] is checked
//! against brute force in the tests; it relies on the Cornacchia–Smith reduction being complete
//! for composite moduli (each primitive solution (x, y) with y > 0 corresponds to the root
//! r = x / y of -d mod m, and Euclid on (m, r) stops at x; Basilla 2004 / Nitaj 1995). Both r and
//! m - r are tried for every root, so the convention for the representative does not matter.
//! The search over guesses is a heuristic (as in the paper): `NotFound` only means no solution
//! with |x3|, |x4| <= bound was found among the guesses not abandoned.
use super::factor::{factor_bounded, product, sqrt_mod_composite, Factors};
use crate::int::Int;
use std::collections::BTreeSet;

fn int(v: i64) -> Int {
    Int::from(v)
}

/// All divisors g >= 1 with g^2 | m, as (g, factorisation of m / g^2).
fn square_divisors(factors: &[(Int, u32)]) -> Vec<(Int, Factors)> {
    let mut out = vec![(Int::one(), factors.to_vec())];
    for (i, (p, k)) in factors.iter().enumerate() {
        let mut next = vec![];
        for (g, fm) in &out {
            let mut gg = g.clone();
            for t in 0..=(k / 2) {
                let mut f2 = fm.clone();
                f2[i].1 = k - 2 * t;
                next.push((gg.clone(), f2));
                gg = &gg * p;
            }
        }
        out = next;
    }
    out.into_iter()
        .map(|(g, f)| (g, f.into_iter().filter(|(_, e)| *e > 0).collect()))
        .collect()
}

/// Primitive solutions (gcd(x, y) = 1), x, y >= 0, of x^2 + d y^2 = m (d > 0, m >= 1).
fn cornacchia_primitive(d: &Int, m: &Int, factors: &[(Int, u32)]) -> BTreeSet<(Int, Int)> {
    let mut out = BTreeSet::new();
    if *m == Int::one() {
        out.insert((Int::one(), Int::zero()));
    }
    if m == d {
        out.insert((Int::zero(), Int::one()));
    }
    if d > m {
        return out; // only y = 0 is possible
    }
    let roots = sqrt_mod_composite(&(-d), m, factors);
    let s = m.isqrt();
    for r in roots {
        if r.is_zero() {
            continue;
        }
        for r0 in [r.clone(), m - &r] {
            let (mut a, mut b) = (m.clone(), r0);
            while b > s {
                let t = a.modulo(&b);
                a = b;
                b = t;
            }
            let rest = m - &(&b * &b);
            if rest.is_neg() || !rest.modulo(d).is_zero() {
                continue;
            }
            let y2 = rest.div_floor(d);
            if y2.is_square() {
                let y = y2.isqrt();
                if b.gcd(&y) == Int::one() {
                    // d = 1: Z[i] has the extra units +-i, so each root class also
                    // carries the swapped solution (y, x), which Euclid never outputs.
                    if *d == Int::one() {
                        out.insert((y.clone(), b.clone()));
                    }
                    out.insert((b, y));
                }
            }
        }
    }
    out
}

/// All solutions (x, y), x, y >= 0, of x^2 + d y^2 = m for d > 0, given the factorisation of
/// m >= 1 (m = 0 gives [(0, 0)], m < 0 gives []): for every g with g^2 | m, the primitive
/// solutions of x^2 + d y^2 = m / g^2 scaled by g. Sorted, without duplicates.
pub fn cornacchia_all(d: &Int, m: &Int, factors: &[(Int, u32)]) -> Vec<(Int, Int)> {
    assert!(!d.is_neg() && !d.is_zero(), "cornacchia_all needs d > 0");
    if m.is_neg() {
        return vec![];
    }
    if m.is_zero() {
        return vec![(Int::zero(), Int::zero())];
    }
    debug_assert_eq!(product(factors), *m);
    let mut out = BTreeSet::new();
    for (g, fm) in square_divisors(factors) {
        let mg = product(&fm);
        for (x, y) in cornacchia_primitive(d, &mg, &fm) {
            out.insert((&x * &g, &y * &g));
        }
    }
    out.into_iter().collect()
}

/// Limits for the binary-quadratic solver.
#[derive(Clone, Debug)]
pub struct BqParams {
    /// Abandon when omega(N) > b * max(1, ln ln N). The slides use B = 11.
    pub b: f64,
    /// Pollard-rho iteration budget for factoring N.
    pub rho_budget: u64,
}
impl Default for BqParams {
    fn default() -> Self {
        BqParams {
            b: 11.0,
            rho_budget: 1 << 20,
        }
    }
}

/// Why a binary quadratic equation was not solved.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AbandonReason {
    /// omega(N) exceeded B log log N (number of square roots 2^omega would be too large).
    TooManyPrimeFactors { omega: usize },
    /// Factoring N exceeded the rho budget.
    FactoringBudget,
}

/// Result of [`solve_binary_quadratic`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BqOutcome {
    /// All integer solutions (possibly none). `factored` says whether N had to be factored.
    Solved {
        sols: Vec<(Int, Int)>,
        factored: bool,
    },
    Abandoned(AbandonReason),
}

fn eval_bq(c: &[Int; 6], x: &Int, y: &Int) -> Int {
    let [a, b, cc, d, e, f] = c;
    &(&(&(&(a * &(x * x)) + &(b * &(x * y))) + &(cc * &(y * y))) + &(&(d * x) + &(e * y))) + f
}

/// The reduction of a x^2 + b xy + c y^2 + d x + e y + f = 0 (D = b^2 - 4ac < 0, a > 0) to
/// u^2 + |D| v^2 = N. Returns (|D|, E, N) where E = 2ae - bd and
///   v = 2a x + b y + d,       u = D y - E,       N = E^2 + D (4af - d^2).
/// Derivation: 4a F = v^2 - D y^2 + 2E y + (4af - d^2); multiplying by -D (> 0) and writing
/// D^2 y^2 - 2DE y = (D y - E)^2 - E^2 gives -4aD F = u^2 + |D| v^2 - N. (Up to the sign of u this
/// is the slides' "x^2 - D y^2 = N" with D < 0.)
pub fn bq_reduction(c: &[Int; 6]) -> (Int, Int, Int) {
    let [a, b, cc, d, e, f] = c;
    let disc = &(b * b) - &(&int(4) * &(a * cc));
    let ee = &(&int(2) * &(a * e)) - &(b * d);
    let n = &(&ee * &ee) + &(&disc * &(&(&int(4) * &(a * f)) - &(d * d)));
    (disc.abs(), ee, n)
}

/// All integer solutions of a x^2 + b xy + c y^2 + d x + e y + f = 0 with b^2 - 4ac < 0 and a > 0,
/// with the default limits (B = 11).
pub fn solve_binary_quadratic(a: &Int, b: &Int, c: &Int, d: &Int, e: &Int, f: &Int) -> BqOutcome {
    solve_binary_quadratic_with(&[a, b, c, d, e, f].map(|v| v.clone()), &BqParams::default())
}

/// As [`solve_binary_quadratic`], coefficients `[a, b, c, d, e, f]`, explicit limits.
pub fn solve_binary_quadratic_with(coef: &[Int; 6], prm: &BqParams) -> BqOutcome {
    let [a, b, _, d, _, _] = coef;
    assert!(!a.is_neg() && !a.is_zero(), "need a > 0");
    let (absd, ee, n) = bq_reduction(coef);
    assert!(!absd.is_zero(), "need b^2 - 4ac < 0");
    let disc = -&absd;
    if n.is_neg() {
        return BqOutcome::Solved {
            sols: vec![],
            factored: false,
        };
    }
    let uv: Vec<(Int, Int)> = if n.is_zero() {
        vec![(Int::zero(), Int::zero())]
    } else {
        let Some(fac) = factor_bounded(&n, prm.rho_budget) else {
            return BqOutcome::Abandoned(AbandonReason::FactoringBudget);
        };
        let lnln = (n.to_f64().ln()).ln().max(1.0);
        if fac.len() as f64 > prm.b * lnln {
            return BqOutcome::Abandoned(AbandonReason::TooManyPrimeFactors { omega: fac.len() });
        }
        cornacchia_all(&absd, &n, &fac)
    };
    let two_a = &int(2) * a;
    let mut sols = BTreeSet::new();
    for (u0, v0) in uv {
        for su in [1i64, -1] {
            for sv in [1i64, -1] {
                let u = &u0 * &int(su);
                let v = &v0 * &int(sv);
                // y = (u + E) / D
                let num = &u + &ee;
                if !num.modulo(&disc).is_zero() {
                    continue;
                }
                let y = num.div_floor(&disc);
                // x = (v - b y - d) / (2a)
                let num = &(&v - &(b * &y)) - d;
                if !num.modulo(&two_a).is_zero() {
                    continue;
                }
                let x = num.div_floor(&two_a);
                if eval_bq(coef, &x, &y).is_zero() {
                    sols.insert((x, y));
                }
            }
        }
    }
    BqOutcome::Solved {
        sols: sols.into_iter().collect(),
        factored: !n.is_zero(),
    }
}

/// Q(x) = sum_{i,j} G[i][j] x_i x_j.
pub fn eval_form(g: &[[Int; 4]; 4], x: &[Int; 4]) -> Int {
    let mut s = Int::zero();
    for i in 0..4 {
        for j in 0..4 {
            s = &s + &(&g[i][j] * &(&x[i] * &x[j]));
        }
    }
    s
}

/// Coefficients [a, b, c, d, e, f] of f(x1, x2) = Q(x1, x2, x3, x4) - target:
///   a = g11, b = 2 g12, c = g22, d = 2(g13 x3 + g14 x4), e = 2(g23 x3 + g24 x4),
///   f = g33 x3^2 + 2 g34 x3 x4 + g44 x4^2 - target
/// (0-based indices in the code). G is assumed symmetric.
pub fn guess_quadratic(g: &[[Int; 4]; 4], target: &Int, x3: &Int, x4: &Int) -> [Int; 6] {
    let two = int(2);
    [
        g[0][0].clone(),
        &two * &g[0][1],
        g[1][1].clone(),
        &two * &(&(&g[0][2] * x3) + &(&g[0][3] * x4)),
        &two * &(&(&g[1][2] * x3) + &(&g[1][3] * x4)),
        &(&(&(&g[2][2] * &(x3 * x3)) + &(&two * &(&g[2][3] * &(x3 * x4))))
            + &(&g[3][3] * &(x4 * x4)))
            - target,
    ]
}

/// Counters for [`solve_norm_eq_cornacchia`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CornacchiaStats {
    /// (x3, x4) guesses tried.
    pub guesses: usize,
    /// Guesses with N < 0 (no solution, no factoring needed).
    pub negative_n: usize,
    /// Guesses for which N was factored.
    pub factorisations: usize,
    /// Guesses abandoned (too many prime factors or factoring budget).
    pub abandoned: usize,
    /// Of which by the omega(N) > B log log N rule.
    pub abandoned_omega: usize,
}

/// Outcome of [`solve_norm_eq_cornacchia`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CornacchiaOutcome {
    Found([Int; 4], CornacchiaStats),
    NotFound(CornacchiaStats),
}

/// The guesses (x3, x4) with max(|x3|, |x4|) <= bound, one of each pair +-(x3, x4) (Q is even,
/// so x is a solution iff -x is), ordered by increasing max(|x3|, |x4|) then |x3| + |x4|.
pub fn guess_order(bound: i64) -> Vec<(i64, i64)> {
    let mut v = vec![];
    for x3 in -bound..=bound {
        for x4 in -bound..=bound {
            if x3 > 0 || (x3 == 0 && x4 >= 0) {
                v.push((x3, x4));
            }
        }
    }
    v.sort_by_key(|&(a, b)| (a.abs().max(b.abs()), a.abs() + b.abs(), a, b));
    v
}

/// Solve Q(x) = d (Q(x) = sum G[i][j] x_i x_j, G symmetric positive definite) by guessing
/// (x3, x4) with |x_i| <= bound (see [`guess_order`]) and solving the binary quadratic equation in
/// (x1, x2) with Cornacchia; returns the first verified solution. Default limits (B = 11).
pub fn solve_norm_eq_cornacchia(gram: &[[Int; 4]; 4], d: &Int, bound: &Int) -> CornacchiaOutcome {
    solve_norm_eq_cornacchia_with(gram, d, bound, &BqParams::default())
}

pub fn solve_norm_eq_cornacchia_with(
    gram: &[[Int; 4]; 4],
    d: &Int,
    bound: &Int,
    prm: &BqParams,
) -> CornacchiaOutcome {
    let bound = bound.to_i128().expect("bound must be small") as i64;
    let mut st = CornacchiaStats::default();
    for (x3, x4) in guess_order(bound) {
        st.guesses += 1;
        let (x3, x4) = (int(x3), int(x4));
        let coef = guess_quadratic(gram, d, &x3, &x4);
        match solve_binary_quadratic_with(&coef, prm) {
            BqOutcome::Abandoned(r) => {
                st.factorisations +=
                    matches!(r, AbandonReason::TooManyPrimeFactors { .. }) as usize;
                st.abandoned += 1;
                st.abandoned_omega +=
                    matches!(r, AbandonReason::TooManyPrimeFactors { .. }) as usize;
            }
            BqOutcome::Solved { sols, factored } => {
                if factored {
                    st.factorisations += 1;
                } else if sols.is_empty() {
                    st.negative_n += 1;
                }
                if let Some((x1, x2)) = sols.into_iter().next() {
                    let x = [x1, x2, x3, x4];
                    assert_eq!(eval_form(gram, &x), *d, "unverified solution");
                    return CornacchiaOutcome::Found(x, st);
                }
            }
        }
    }
    CornacchiaOutcome::NotFound(st)
}
