//! # Lead E4 of research/isogeny_conductor_gap_round3_20261011
//!
//! Two-torsion-symmetrised summation-polynomial systems over `F_p`, crater
//! against floor, on the certified pairs of round 1.
//!
//! For a rational subgroup `H ⊆ E(F_p)[2^∞]` the Vélu isogeny
//! `ψ: E → E' = E/H` turns a decomposition `P₁ + P₂ = R` over an
//! `H`-stable factor base `V` into `ψ(P₁) + ψ(P₂) = ψ(R)` over
//! `W = X(V)`, where `X` is the Vélu abscissa map.  The quotient is up to
//! `H`, so the rows it yields are rows of the cofactor-projected relation
//! matrix (`H ⊆ E(F_p)[h]`).
//!
//! Translation by `T ∈ H` preserves `S₃` on `E` in the sense that
//! `(P₁ + T) + (P₂ − T) = R`.  So writing the system on `E'` *is* the
//! symmetrisation: the factor-base polynomial falls from degree `|V|` to
//! `|V| / |H|`, and `S₃'` is the same shape.
//!
//! The run compares, at a matched factor-base size `B = |V|` on `E`:
//!
//! * **direct**: `{S₃(x₁, x₂, x_R), F(x₁), F(x₂)}` on `E`, with
//!   `F = Π_{v∈V}(x − v)`;
//! * **quotient by `H`**: `{S₃'(X₁, X₂, X(x_R)), G(X₁), G(X₂)}` on
//!   `E/H`, with `G = Π_{w∈W}(X − w)`, for every `H` the curve has
//!   among `C2`, `C4` and `E[2]` (the last only where `E[2]` is rational,
//!   which a ℓ = 2 volcano gives the crater and not the floor).
//!
//! It also covers crater, floor and the isomorphic floor model of every
//! even-cofactor prime pair.  Each system is solved by the repository's
//! degree-bounded F4 over `F_p` (`cryptanalysis::f4_fp`).  The run records
//! the solving degree, the widest matrix up to the solution, the `F_p`
//! multiplications of row reduction, and wall time.
//!
//! Each curve gets two random targets and one planted target
//! `R = P_i + P_j` with `P_i, P_j` in the base.  The planted solution must
//! annihilate the basis; that checks the Vélu map and the system end to
//! end.  Before any system is built, the Vélu map is checked on 32 random
//! sums: `S₃'(X(P₁), X(P₂), X(P₁ + P₂)) = 0`.
//!
//! ```bash
//! cargo run --release --example isogeny_gap_two_torsion_f4 -- \
//!   --pairs research/isogeny_conductor_gap_ic_20261010/pairs/prime_pairs.json \
//!   --out research/isogeny_conductor_gap_round3_20261011/e4/records.jsonl \
//!   --sizes 32,64 --budget-ms 20000
//! ```
//!
//! Records append one JSON line per system; a restart skips every curve
//! whose records are complete.

use std::collections::{BTreeSet, HashSet};
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

// ── driver ─────────────────────────────────────────────────────────

fn arg(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .map(|i| args[i + 1].clone())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let pairs_path = arg(&args, "--pairs").expect("--pairs");
    let out_path = arg(&args, "--out").expect("--out");
    let sizes: Vec<usize> = arg(&args, "--sizes")
        .unwrap_or("32,64".into())
        .split(',')
        .map(|s| s.parse().expect("size"))
        .collect();
    let budget = Duration::from_millis(
        arg(&args, "--budget-ms").map_or(20000, |s| s.parse().expect("budget")),
    );
    let pairs: Value =
        serde_json::from_str(&std::fs::read_to_string(&pairs_path).unwrap()).unwrap();

    // Curves already complete in the output.
    let mut done: BTreeSet<String> = BTreeSet::new();
    if let Ok(text) = std::fs::read_to_string(&out_path) {
        for line in text.lines() {
            let v: Value = serde_json::from_str(line).expect("record");
            if v["kind"] == "curve_done" {
                done.insert(v["curve_key"].as_str().unwrap().to_string());
            }
        }
    }
    let mut out = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&out_path)
        .expect("open out");

    for pr in pairs["pairs"].as_array().unwrap() {
        if !pr["cofactor"].as_u64().unwrap().is_multiple_of(2) {
            continue;
        }
        let pair = format!(
            "D={} {}^{} {}b h={}",
            pr["d_k"], pr["ell"], pr["h"], pr["bits"], pr["cofactor"]
        );
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
            let key = format!("{pair}|{role}");
            if done.contains(&key) {
                continue;
            }
            let e = Curve {
                p: spec["p"].as_u64().unwrap(),
                a: spec["a"].as_u64().unwrap(),
                b: spec["b"].as_u64().unwrap(),
            };
            let order = spec["group_order"].as_u64().unwrap();
            let mut pts = Points {
                state: e.p ^ e.a.rotate_left(17) ^ e.b.rotate_left(41),
            };
            let two = two_primary(&e, order, &mut pts);
            let rank2 = two
                .iter()
                .filter(|&&q| q != Pt::O && order_of(&e, q) == 2)
                .count();
            let two_rank = if rank2 == 3 { 2 } else { 1 };
            let mut groups: Vec<(&str, Vec<Pt>)> = vec![("1", vec![Pt::O])];
            groups.extend(subgroups(&e, &two));
            for (hname, h) in &groups {
                let velu = Velu::new(&e, h);
                // The Vélu map, checked on random sums.
                for _ in 0..32 {
                    let (p1, p2) = (pts.next(&e), pts.next(&e));
                    let p3 = e.add(p1, p2);
                    let xs: Vec<u64> = [p1, p2, p3]
                        .iter()
                        .filter_map(|q| match q {
                            Pt::A(x, _) if !velu.kernel_x.contains(x) => Some(*x),
                            _ => None,
                        })
                        .collect();
                    if xs.len() == 3 {
                        let (a1, a2, a3) = (velu.x(xs[0]), velu.x(xs[1]), velu.x(xs[2]));
                        assert_eq!(
                            velu.target.s3(a1, a2, a3),
                            0,
                            "Vélu map fails S₃' on {key} {hname}"
                        );
                        assert!(
                            is_square(velu.target.rhs(a1), e.p),
                            "image off E' on {key} {hname}"
                        );
                    }
                }
                for &size in &sizes {
                    // An H-stable base of `size` abscissae: the smallest
                    // abscissae on the curve, each with its H-orbit.
                    let mut base: Vec<Pt> = Vec::new();
                    let mut xs: BTreeSet<u64> = BTreeSet::new();
                    let mut x0 = 1u64;
                    while xs.len() < size {
                        if let Some(q) = e.lift(x0) {
                            if !velu.kernel_x.contains(&x0) && !xs.contains(&x0) {
                                let orbit: Vec<Pt> = h.iter().map(|&k| e.add(q, k)).collect();
                                // Only orbits of |H| distinct abscissae, so the
                                // base reaches exactly `size` (an orbit with
                                // `P + K = −(P + K')` repeats an abscissa).
                                let orbit_xs: BTreeSet<u64> = orbit
                                    .iter()
                                    .filter_map(|o| match o {
                                        Pt::A(x, _) => Some(*x),
                                        Pt::O => None,
                                    })
                                    .collect();
                                if orbit_xs.len() == h.len()
                                    && orbit_xs.iter().all(|x| !xs.contains(x))
                                    && xs.len() + orbit.len() <= size
                                {
                                    for o in orbit {
                                        if let Pt::A(x, _) = o {
                                            xs.insert(x);
                                            base.push(o);
                                        }
                                    }
                                }
                            }
                        }
                        x0 += 1;
                    }
                    let v: Vec<u64> = xs.iter().copied().collect();
                    let w: Vec<u64> = v
                        .iter()
                        .map(|&x| velu.x(x))
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect();
                    // Targets: two random, one planted.
                    let mut targets: Vec<(&str, Pt, Option<(Pt, Pt)>)> = vec![
                        ("random", pts.next(&e), None),
                        ("random", pts.next(&e), None),
                    ];
                    let (pi, pj) = (base[1 % base.len()], base[base.len() / 2 + 1]);
                    targets.push(("planted", e.add(pi, pj), Some((pi, pj))));
                    for (ti, (tkind, r, planted)) in targets.into_iter().enumerate() {
                        let Pt::A(xr, _) = r else { continue };
                        if velu.kernel_x.contains(&xr) {
                            continue;
                        }
                        let mut write = |method: &str, rec: Value| {
                            let mut line = json!({
                                "kind": "system", "pair": pair, "role": role, "curve_key": key,
                                "p": e.p, "a": e.a, "b": e.b, "group_order": order,
                                "cofactor": pr["cofactor"], "ell": pr["ell"], "d_k": pr["d_k"], "bits": pr["bits"],
                                "two_rank": two_rank, "two_part": two.len(),
                                "subgroup": hname, "subgroup_order": h.len(),
                                "base_size": size, "method": method,
                                "target_index": ti, "target_kind": tkind,
                            });
                            line.as_object_mut()
                                .unwrap()
                                .extend(rec.as_object().unwrap().clone());
                            writeln!(out, "{line}").unwrap();
                        };
                        let pl = |q: Pt, f: &dyn Fn(u64) -> u64| match q {
                            Pt::A(x, _) => f(x),
                            Pt::O => unreachable!(),
                        };
                        if *hname == "1" {
                            let planted_x = planted.map(|(a, b)| (pl(a, &|x| x), pl(b, &|x| x)));
                            write("direct", run_system(&e, &v, xr, planted_x, budget).rec);
                        } else {
                            let planted_x = planted
                                .map(|(a, b)| (pl(a, &|x| velu.x(x)), pl(b, &|x| velu.x(x))));
                            write(
                                "quotient",
                                run_system(&velu.target, &w, velu.x(xr), planted_x, budget).rec,
                            );
                        }
                    }
                    out.flush().unwrap();
                }
            }
            writeln!(
                out,
                "{}",
                json!({"kind": "curve_done", "curve_key": key, "two_rank": two_rank,
                       "two_part": two.len(), "subgroups": groups.iter().map(|(n, _)| *n).collect::<Vec<_>>()})
            )
            .unwrap();
            out.flush().unwrap();
            eprintln!(
                "{key}: rank {two_rank}, 2-part {}, subgroups {:?}",
                two.len(),
                groups.iter().map(|(n, _)| *n).collect::<Vec<_>>()
            );
        }
    }
}
