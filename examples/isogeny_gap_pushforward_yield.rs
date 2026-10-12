//! # Lead E6 of research/isogeny_conductor_gap_round5_20261012
//!
//! The yield of a pushed-forward factor base, against the direct one, at a
//! matched number of columns.
//!
//! For a rational subgroup `H ⊆ E(F_p)[h]` and its Vélu isogeny
//! `ψ: E → E/H`:
//! - an arbitrary base `V` of `c` abscissae on `E` gives the **direct**
//!   system `S₃(x₁, x₂, x_T)` over `V`;
//! - the **quotient** system is `S₃'(X₁, X₂, X(x_T))` over `W = X(V)`.
//!
//! Both have `c` columns and degree `c`.  A quotient solution is a
//! decomposition `±P₁ ± P₂ = ±T + K` with `K ∈ H`, which is a valid row
//! after the `[h]` projection.  The prediction is `|H|` times the direct
//! yield at the same F4 cost.
//!
//! Per curve, base size and arm, the runner records:
//! - over 20,000 shared targets `[h]·(random point)`, the targets with a
//!   decomposition and the total decompositions, by exact enumeration of
//!   base pairs;
//! - every hit lifted back to points of `E` and checked;
//! - one F4 system's solving degree, widest matrix and multiplications.
//!
//! The field, curve, torsion and Vélu code is round 3's
//! (`examples/isogeny_gap_two_torsion_f4.rs`), copied unchanged.
//!
//! ```bash
//! cargo run --release --example isogeny_gap_pushforward_yield -- \
//!   --pairs research/isogeny_conductor_gap_ic_20261010/pairs/prime_pairs.json \
//!   --out research/isogeny_conductor_gap_round5_20261012/e6/records.jsonl
//! ```

use std::collections::HashSet;
use std::io::Write as _;
use std::time::Duration;

use crypto_lib::cryptanalysis::f4_fp::{self, F4Options, Ordering, Poly};
use serde_json::{json, Value};

// ── F_p ────────────────────────────────────────────────────────────

fn mul(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
fn add(a: u64, b: u64, p: u64) -> u64 {
    (a + b) % p
}
fn sub(a: u64, b: u64, p: u64) -> u64 {
    (a + p - b % p) % p
}
fn pow(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1;
    a %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mul(r, a, p);
        }
        a = mul(a, a, p);
        e >>= 1;
    }
    r
}
fn inv(a: u64, p: u64) -> u64 {
    assert!(!a.is_multiple_of(p), "inverse of zero");
    pow(a, p - 2, p)
}
fn is_square(a: u64, p: u64) -> bool {
    a.is_multiple_of(p) || pow(a, (p - 1) / 2, p) == 1
}
fn sqrt(a: u64, p: u64) -> Option<u64> {
    let a = a % p;
    if a == 0 {
        return Some(0);
    }
    if !is_square(a, p) {
        return None;
    }
    // Tonelli–Shanks.
    let (mut q, mut s) = (p - 1, 0u32);
    while q.is_multiple_of(2) {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while is_square(z, p) {
        z += 1;
    }
    let (mut m, mut c, mut t, mut r) = (s, pow(z, q, p), pow(a, q, p), pow(a, q.div_ceil(2), p));
    while t != 1 {
        let mut i = 0;
        let mut tt = t;
        while tt != 1 {
            tt = mul(tt, tt, p);
            i += 1;
        }
        let b = pow(c, 1 << (m - i - 1), p);
        m = i;
        c = mul(b, b, p);
        t = mul(t, c, p);
        r = mul(r, b, p);
    }
    Some(r)
}

// ── the curve y² = x³ + ax + b ─────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
enum Pt {
    O,
    A(u64, u64),
}

#[derive(Clone, Copy)]
struct Curve {
    p: u64,
    a: u64,
    b: u64,
}

impl Curve {
    fn rhs(&self, x: u64) -> u64 {
        let p = self.p;
        add(
            add(mul(mul(x, x, p), x, p), mul(self.a, x, p), p),
            self.b,
            p,
        )
    }
    fn neg(&self, q: Pt) -> Pt {
        match q {
            Pt::O => Pt::O,
            Pt::A(x, y) => Pt::A(x, (self.p - y) % self.p),
        }
    }
    fn add(&self, s: Pt, t: Pt) -> Pt {
        let p = self.p;
        match (s, t) {
            (Pt::O, q) | (q, Pt::O) => q,
            (Pt::A(x1, y1), Pt::A(x2, y2)) => {
                let lam = if x1 == x2 {
                    if add(y1, y2, p) == 0 {
                        return Pt::O;
                    }
                    mul(
                        add(mul(3, mul(x1, x1, p), p), self.a, p),
                        inv(mul(2, y1, p), p),
                        p,
                    )
                } else {
                    mul(sub(y2, y1, p), inv(sub(x2, x1, p), p), p)
                };
                let x3 = sub(sub(mul(lam, lam, p), x1, p), x2, p);
                let y3 = sub(mul(lam, sub(x1, x3, p), p), y1, p);
                Pt::A(x3, y3)
            }
        }
    }
    fn mul(&self, mut k: u64, q: Pt) -> Pt {
        let mut acc = Pt::O;
        let mut base = q;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(acc, base);
            }
            base = self.add(base, base);
            k >>= 1;
        }
        acc
    }
    /// A point with abscissa `x`, if there is one.
    fn lift(&self, x: u64) -> Option<Pt> {
        sqrt(self.rhs(x), self.p).map(|y| Pt::A(x % self.p, y))
    }
    /// `S₃(x₁, x₂, x₃)`, Semaev's third summation polynomial.
    fn s3(&self, x1: u64, x2: u64, x3: u64) -> u64 {
        let p = self.p;
        let (a, b) = (self.a, self.b);
        let d = sub(x1, x2, p);
        let t1 = mul(mul(d, d, p), mul(x3, x3, p), p);
        let inner = add(
            mul(add(x1, x2, p), add(mul(x1, x2, p), a, p), p),
            mul(2, b, p),
            p,
        );
        let t2 = mul(mul(2, inner, p), x3, p);
        let e = sub(mul(x1, x2, p), a, p);
        let t3 = sub(mul(e, e, p), mul(mul(4, b, p), add(x1, x2, p), p), p);
        add(sub(t1, t2, p), t3, p)
    }
    /// `S₃(x₁, x₂, c)` as a polynomial in `(x₁, x₂)`.
    fn s3_poly(&self, c: u64) -> Vec<(Vec<u32>, u64)> {
        let p = self.p;
        let (a, b) = (self.a, self.b);
        let neg = |v: u64| (p - v % p) % p;
        let c2 = mul(c, c, p);
        vec![
            (vec![2, 2], 1),
            (vec![2, 1], neg(mul(2, c, p))),
            (vec![1, 2], neg(mul(2, c, p))),
            (vec![2, 0], c2),
            (vec![0, 2], c2),
            (vec![1, 1], neg(add(mul(2, c2, p), mul(2, a, p), p))),
            (
                vec![1, 0],
                neg(add(mul(mul(2, c, p), a, p), mul(4, b, p), p)),
            ),
            (
                vec![0, 1],
                neg(add(mul(mul(2, c, p), a, p), mul(4, b, p), p)),
            ),
            (vec![0, 0], sub(mul(a, a, p), mul(mul(4, c, p), b, p), p)),
        ]
    }
}

// ── the 2-primary torsion and its subgroups ────────────────────────

/// A deterministic stream of points: abscissae from a 64-bit LCG.
struct Points {
    state: u64,
}
impl Points {
    fn next(&mut self, e: &Curve) -> Pt {
        loop {
            self.state = self
                .state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            if let Some(q) = e.lift((self.state >> 11) % e.p) {
                return q;
            }
        }
    }
}

/// `E(F_p)[2^∞]`, generated from `[N_odd]P` over sampled points.
fn two_primary(e: &Curve, order: u64, pts: &mut Points) -> Vec<Pt> {
    let v = order.trailing_zeros();
    let n_odd = order >> v;
    let full = 1usize << v;
    let mut group = vec![Pt::O];
    let mut keys: HashSet<Pt> = HashSet::from([Pt::O]);
    let mut tries = 0;
    while group.len() < full && tries < 256 {
        tries += 1;
        let k = e.mul(n_odd, pts.next(e));
        if keys.contains(&k) {
            continue;
        }
        let old = group.clone();
        let mut m = k;
        while !keys.contains(&m) {
            for &o in &old {
                let q = e.add(o, m);
                if keys.insert(q) {
                    group.push(q);
                }
            }
            m = e.add(m, k);
        }
    }
    assert_eq!(group.len(), full, "the 2-primary part has 2^{v} points");
    group
}

fn order_of(e: &Curve, q: Pt) -> u64 {
    let mut k = 1;
    let mut m = q;
    while m != Pt::O {
        m = e.add(m, q);
        k += 1;
    }
    k
}

/// The subgroups tested on a curve: `C2`, `C4`, `E[2]` (rank 2 only).
fn subgroups(e: &Curve, two: &[Pt]) -> Vec<(&'static str, Vec<Pt>)> {
    let mut out = Vec::new();
    let t2: Vec<Pt> = two
        .iter()
        .copied()
        .filter(|&q| q != Pt::O && order_of(e, q) == 2)
        .collect();
    if let Some(&t) = t2.first() {
        out.push(("C2", vec![Pt::O, t]));
    }
    if let Some(&q) = two.iter().find(|&&q| order_of(e, q) == 4) {
        let q2 = e.add(q, q);
        out.push(("C4", vec![Pt::O, q, q2, e.add(q2, q)]));
    }
    if t2.len() == 3 {
        let mut v4 = vec![Pt::O];
        v4.extend(t2);
        out.push(("E[2]", v4));
    }
    out
}

// ── Vélu ───────────────────────────────────────────────────────────

/// The quotient `E/H` and its abscissa map (Vélu; Washington, Thm 12.16).
struct Velu {
    target: Curve,
    /// `(x_Q, v_Q, u_Q)` over one representative of each `±Q ≠ O` in `H`.
    terms: Vec<(u64, u64, u64)>,
    kernel_x: HashSet<u64>,
}

impl Velu {
    fn new(e: &Curve, h: &[Pt]) -> Self {
        let p = e.p;
        let mut terms = Vec::new();
        let mut seen: HashSet<Pt> = HashSet::new();
        let (mut t, mut w) = (0u64, 0u64);
        let mut kernel_x = HashSet::new();
        for &q in h {
            let Pt::A(x, y) = q else { continue };
            kernel_x.insert(x);
            if seen.contains(&e.neg(q)) {
                continue;
            }
            seen.insert(q);
            let gx = add(mul(3, mul(x, x, p), p), e.a, p);
            let two_torsion = y == 0;
            let v = if two_torsion { gx } else { mul(2, gx, p) };
            let u = mul(4, mul(y, y, p), p);
            t = add(t, v, p);
            w = add(w, add(u, mul(x, v, p), p), p);
            terms.push((x, v, u));
        }
        let target = Curve {
            p,
            a: sub(e.a, mul(5, t, p), p),
            b: sub(e.b, mul(7, w, p), p),
        };
        Velu {
            target,
            terms,
            kernel_x,
        }
    }
    fn x(&self, x: u64) -> u64 {
        let p = self.target.p;
        let mut acc = x % p;
        for &(xq, v, u) in &self.terms {
            let d = inv(sub(x, xq, p), p);
            acc = add(acc, add(mul(v, d, p), mul(u, mul(d, d, p), p), p), p);
        }
        acc
    }
}

// ── systems ────────────────────────────────────────────────────────

/// `Π (X − w)` in variable `var` of two.
fn root_poly(roots: &[u64], var: usize, p: u64) -> Vec<(Vec<u32>, u64)> {
    let mut c: Vec<u64> = vec![1];
    for &w in roots {
        let mut next = vec![0u64; c.len() + 1];
        for (k, &ck) in c.iter().enumerate() {
            next[k + 1] = add(next[k + 1], ck, p);
            next[k] = sub(next[k], mul(ck, w, p), p);
        }
        c = next;
    }
    c.iter()
        .enumerate()
        .filter(|(_, &ck)| ck != 0)
        .map(|(k, &ck)| {
            let mut e = vec![0u32; 2];
            e[var] = k as u32;
            (e, ck)
        })
        .collect()
}

struct Outcome {
    rec: Value,
}

fn run_system(
    curve: &Curve,
    roots: &[u64],
    xr: u64,
    planted: Option<(u64, u64)>,
    budget: Duration,
) -> Outcome {
    let p = curve.p;
    let input: Vec<Poly> = vec![
        f4_fp::normalise(&curve.s3_poly(xr), p, Ordering::Grevlex),
        f4_fp::normalise(&root_poly(roots, 0, p), p, Ordering::Grevlex),
        f4_fp::normalise(&root_poly(roots, 1, p), p, Ordering::Grevlex),
    ];
    let max_degree = 2 * roots.len() as u32 + 8;
    let opts = F4Options::new(Ordering::Grevlex, max_degree).with_budget(budget);
    let rep = f4_fp::f4(&input, 2, p, &opts);
    // Solutions by enumeration over the base, independently of F4.
    let mut solutions = 0u64;
    for (i, &w1) in roots.iter().enumerate() {
        for &w2 in &roots[i..] {
            if curve.s3(w1, w2, xr) == 0 {
                solutions += 1;
            }
        }
    }
    let planted_ok = planted.map(|(w1, w2)| {
        !rep.inconsistent && rep.basis.iter().all(|f| f4_fp::eval(f, &[w1, w2], p) == 0)
    });
    let consistent_with_enum = rep.timed_out || (rep.inconsistent == (solutions == 0));
    Outcome {
        rec: json!({
            "base_degree": roots.len(),
            "solving_degree_max": rep.solving_degree_max,
            "solving_degree": rep.solving_degree,
            "degree_reached": rep.degree_reached,
            "max_cols_to_solution": rep.max_cols_to_solution,
            "max_rows": rep.max_rows, "max_cols": rep.max_cols,
            "steps_to_solution": rep.steps_to_solution,
            "field_ops": rep.field_ops,
            "ms": rep.ms,
            "timed_out": rep.timed_out,
            "pairs_above_bound": rep.pairs_above_bound,
            "inconsistent": rep.inconsistent,
            "solutions_by_enumeration": solutions,
            "planted_annihilates_basis": planted_ok,
            "f4_agrees_with_enumeration": consistent_with_enum,
        }),
    }
}

// ── E6 driver ──────────────────────────────────────────────────────

/// `E(F_p)[h]`, generated from `[r]P` over sampled points (coset closure).
fn cofactor_torsion(e: &Curve, r: u64, h: u64, pts: &mut Points) -> Vec<Pt> {
    let mut group = vec![Pt::O];
    let mut keys: HashSet<Pt> = HashSet::from([Pt::O]);
    let mut tries = 0;
    while (group.len() as u64) < h && tries < 512 {
        tries += 1;
        let k = e.mul(r, pts.next(e));
        if keys.contains(&k) {
            continue;
        }
        let old = group.clone();
        let mut m = k;
        while !keys.contains(&m) {
            for &o in &old {
                let q = e.add(o, m);
                if keys.insert(q) {
                    group.push(q);
                }
            }
            m = e.add(m, k);
        }
    }
    assert_eq!(group.len() as u64, h, "E(F_p)[h] has h points");
    group
}

fn arg(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .map(|i| args[i + 1].clone())
}

/// Decompositions of every target over a base, with each one lifted back.
/// `xs` are the system's base abscissae (on the curve the system is
/// written on), `pts_of` maps one to the base points of `E` behind it,
/// `target_x` gives the system's target abscissa, and `in_h` tells whether
/// a point of `E` is in `H` (`{O}` for the direct system).
struct Yield {
    targets_hit: u64,
    decompositions: u64,
    lift_failures: u64,
}

#[allow(clippy::too_many_arguments)]
fn yield_of(
    sys: &Curve,
    e: &Curve,
    xs: &[u64],
    pts_of: &std::collections::HashMap<u64, Vec<Pt>>,
    targets: &[Pt],
    target_x: &dyn Fn(u64) -> u64,
    in_h: &dyn Fn(Pt) -> bool,
) -> Yield {
    let mut y = Yield {
        targets_hit: 0,
        decompositions: 0,
        lift_failures: 0,
    };
    for &t in targets {
        let Pt::A(xt, _) = t else { continue };
        let xr = target_x(xt);
        let mut hit = false;
        for (i, &w1) in xs.iter().enumerate() {
            for &w2 in &xs[i..] {
                if sys.s3(w1, w2, xr) != 0 {
                    continue;
                }
                y.decompositions += 1;
                hit = true;
                // Lift: some base points and signs with ±P₁ ± P₂ − T ∈ H.
                let mut ok = false;
                'lift: for &p1 in &pts_of[&w1] {
                    for &p2 in &pts_of[&w2] {
                        for s1 in [false, true] {
                            for s2 in [false, true] {
                                let a = if s1 { e.neg(p1) } else { p1 };
                                let b = if s2 { e.neg(p2) } else { p2 };
                                let d = e.add(e.add(a, b), e.neg(t));
                                if in_h(d) || in_h(e.add(e.add(a, b), t)) {
                                    ok = true;
                                    break 'lift;
                                }
                            }
                        }
                    }
                }
                if !ok {
                    y.lift_failures += 1;
                }
            }
        }
        y.targets_hit += u64::from(hit);
    }
    y
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let pairs: Value = serde_json::from_str(
        &std::fs::read_to_string(arg(&args, "--pairs").expect("--pairs")).unwrap(),
    )
    .unwrap();
    let out_path = arg(&args, "--out").expect("--out");
    let n_targets: usize = arg(&args, "--targets").map_or(20000, |s| s.parse().unwrap());
    let selected = [(-7i64, 3u64), (-7, 5), (-3, 2), (-8, 3), (-4, 5)];
    let mut out = std::fs::File::create(&out_path).expect("create out");
    for pr in pairs["pairs"].as_array().unwrap() {
        let key = (pr["d_k"].as_i64().unwrap(), pr["ell"].as_u64().unwrap());
        if pr["bits"].as_u64() != Some(22) || !selected.contains(&key) {
            continue;
        }
        let pair = format!(
            "D={} {}^{} {}b h={}",
            pr["d_k"], pr["ell"], pr["h"], pr["bits"], pr["cofactor"]
        );
        let h = pr["cofactor"].as_u64().unwrap();
        let mut roles: Vec<(String, Value)> = pr["curves"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|c| c["role"] == "crater" || c["role"] == "floor")
            .map(|c| (c["role"].as_str().unwrap().to_string(), c["spec"].clone()))
            .collect();
        roles.push((
            "floor-iso".into(),
            pr["floor_isomorphic_model"]["spec"].clone(),
        ));
        for (role, spec) in roles {
            let e = Curve {
                p: spec["p"].as_u64().unwrap(),
                a: spec["a"].as_u64().unwrap(),
                b: spec["b"].as_u64().unwrap(),
            };
            let order = spec["group_order"].as_u64().unwrap();
            let r = spec["r"].as_u64().unwrap();
            let mut pts = Points {
                state: e.p ^ e.a.rotate_left(13) ^ e.b.rotate_left(37) ^ 0xe6,
            };
            let two = two_primary(&e, order, &mut pts);
            let mut groups: Vec<(&str, Vec<Pt>)> = subgroups(&e, &two);
            groups.push(("E[h]", cofactor_torsion(&e, r, h, &mut pts)));
            let targets: Vec<Pt> = (0..n_targets)
                .map(|_| loop {
                    let t = e.mul(h, pts.next(&e));
                    if t != Pt::O {
                        break t;
                    }
                })
                .collect();
            for &c in &[32usize, 64] {
                // The direct base: the c smallest abscissae outside E[h].
                let full_h: HashSet<Pt> = groups.last().unwrap().1.iter().copied().collect();
                let mut base: Vec<Pt> = Vec::new();
                let mut x0 = 1u64;
                while base.len() < c {
                    if let Some(q) = e.lift(x0) {
                        if !full_h.contains(&q) {
                            base.push(q);
                        }
                    }
                    x0 += 1;
                }
                let v: Vec<u64> = base
                    .iter()
                    .map(|q| match q {
                        Pt::A(x, _) => *x,
                        Pt::O => unreachable!(),
                    })
                    .collect();
                let mut direct_pts: std::collections::HashMap<u64, Vec<Pt>> =
                    std::collections::HashMap::new();
                for (&x, &q) in v.iter().zip(&base) {
                    direct_pts.entry(x).or_default().push(q);
                }
                let is_o = |q: Pt| q == Pt::O;
                let d = yield_of(&e, &e, &v, &direct_pts, &targets, &|x| x, &is_o);
                let f4_direct = run_system(
                    &e,
                    &v,
                    match targets[0] {
                        Pt::A(x, _) => x,
                        Pt::O => unreachable!(),
                    },
                    None,
                    Duration::from_secs(20),
                )
                .rec;
                let mut write = |arm: &str, h_order: usize, cols: usize, y: &Yield, f4: &Value| {
                    let line = json!({
                        "kind": "arm", "pair": pair, "role": role, "p": e.p, "group_order": order,
                        "r": r, "cofactor": h, "base_size": c, "arm": arm, "subgroup_order": h_order,
                        "columns": cols, "targets": targets.len(),
                        "targets_hit": y.targets_hit, "decompositions": y.decompositions,
                        "lift_failures": y.lift_failures, "f4": f4,
                    });
                    writeln!(out, "{line}").unwrap();
                };
                write("direct", 1, v.len(), &d, &f4_direct);
                for (hname, hset) in &groups {
                    let velu = Velu::new(&e, hset);
                    let hs: HashSet<Pt> = hset.iter().copied().collect();
                    let mut w_pts: std::collections::HashMap<u64, Vec<Pt>> =
                        std::collections::HashMap::new();
                    for (&x, &q) in v.iter().zip(&base) {
                        assert!(!velu.kernel_x.contains(&x), "base points lie outside H");
                        w_pts.entry(velu.x(x)).or_default().push(q);
                    }
                    let mut w: Vec<u64> = w_pts.keys().copied().collect();
                    w.sort();
                    let in_h = |q: Pt| hs.contains(&q);
                    let qy = yield_of(
                        &velu.target,
                        &e,
                        &w,
                        &w_pts,
                        &targets,
                        &|x| velu.x(x),
                        &in_h,
                    );
                    let xt0 = match targets[0] {
                        Pt::A(x, _) => velu.x(x),
                        Pt::O => unreachable!(),
                    };
                    let f4 = run_system(&velu.target, &w, xt0, None, Duration::from_secs(20)).rec;
                    write(&format!("quotient {hname}"), hset.len(), w.len(), &qy, &f4);
                }
                out.flush().unwrap();
            }
            eprintln!("{pair} {role}: done");
        }
    }
}
