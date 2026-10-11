//! Round-4 lead 16.7: sanity checks for the isogeny-oracle AGM note.
//!
//! Toy curve `E: y² = x³ + a x² + b x` over `F_p` (`p ≡ 3 mod 4`, about
//! 2^20) with one rational 2-torsion point `T = (0,0)` and order `N = 4r`,
//! `r` prime, together with the explicit 2-isogeny
//! `φ: E → E'`, `E': Y² = X³ − 2a X² + (a² − 4b) X`,
//! `φ(x,y) = (y²/x², y(b − x²)/x²)` and its dual
//! `φ̂(X,Y) = (Y²/(4X²), Y(a² − 4b − X²)/(8X²))` (Silverman, AEC III.4.5).
//!
//! Checks (every count is exact; nothing is sampled except the AGM query
//! script and the random points of check 4):
//!
//! 1. `φ`, `φ̂` land on the right curves, `ker φ = {O, T}`, and
//!    `φ̂∘φ = [2]`, `φ∘φ̂ = [2]` on random points; the scalar law
//!    `λ_φ̂·λ_φ = deg φ = 2` on the order-`r` subgroup.
//! 2. AGM oracle simulation (Theorem 1 of the note): a random algebraic
//!    query script against the real oracles `φ`, `φ̂` and against the
//!    bookkeeping simulator that evaluates `φ` once (on the challenge) and
//!    never again.  Every answer must agree.
//! 3. Pushforward of decompositions (Theorem 5): for `m ∈ {2, 3}` and three
//!    factor bases on `E` (plain small-x, φ-injective, T-saturated), the
//!    fiber-sum identity `Σ_{φ(R)=Q'} d_F(R) = (w_F)^{*m}(Q')` at every
//!    `Q'`, the injective form, and the saturated form
//!    `d_F(R) = 2^{m−1} d_{φ(F)}(φ(R))` at every `R`.
//! 4. The cofactor split used for unrestricted queries:
//!    `φ(Z) = φ(Z_G) + φ(Z_h)` with the idempotent decomposition.
//! 5. A class-factor census (distinct m-sums landing in the order-r
//!    subgroup) for the plain and image bases on `E'` and the saturated
//!    base on `E`, to set against the heuristic image/saturation laws.
//!
//! ```sh
//! cargo run --release --example r4_16_7_agm_pushforward -- \
//!     --pmin 1048576 --queries 5000 --seed 7 > out.json
//! ```

use crypto_lib::cryptanalysis::curve_id;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::json;
use std::collections::HashMap;

type Pt = Option<(u64, u64)>;

fn mulm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
fn addm(a: u64, b: u64, p: u64) -> u64 {
    (a + b) % p
}
fn subm(a: u64, b: u64, p: u64) -> u64 {
    (a + p - b) % p
}
fn powm(mut b: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mulm(r, b, p);
        }
        b = mulm(b, b, p);
        e >>= 1;
    }
    r
}
fn invm(a: u64, p: u64) -> u64 {
    powm(a, p - 2, p)
}
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    let mut d = 2u64;
    while d * d <= n {
        if n.is_multiple_of(d) {
            return false;
        }
        d += 1;
    }
    true
}
fn chi(v: u64, p: u64) -> i64 {
    if v == 0 {
        0
    } else if powm(v, (p - 1) / 2, p) == 1 {
        1
    } else {
        -1
    }
}

/// `y² = x³ + a x² + b x` over `F_p`.
#[derive(Clone, Copy)]
struct Curve {
    p: u64,
    a: u64,
    b: u64,
}

impl Curve {
    fn rhs(&self, x: u64) -> u64 {
        let p = self.p;
        let x2 = mulm(x, x, p);
        addm(
            addm(mulm(x2, x, p), mulm(self.a, x2, p), p),
            mulm(self.b, x, p),
            p,
        )
    }
    fn on(&self, q: Pt) -> bool {
        match q {
            None => true,
            Some((x, y)) => mulm(y, y, self.p) == self.rhs(x),
        }
    }
    fn neg(&self, q: Pt) -> Pt {
        q.map(|(x, y)| (x, (self.p - y) % self.p))
    }
    fn add(&self, s: Pt, t: Pt) -> Pt {
        let p = self.p;
        let (x1, y1) = match s {
            None => return t,
            Some(v) => v,
        };
        let (x2, y2) = match t {
            None => return s,
            Some(v) => v,
        };
        let lam = if x1 == x2 {
            if addm(y1, y2, p) == 0 {
                return None;
            }
            // (3x² + 2a x + b) / 2y
            let num = addm(
                addm(
                    mulm(3, mulm(x1, x1, p), p),
                    mulm(2, mulm(self.a, x1, p), p),
                    p,
                ),
                self.b,
                p,
            );
            mulm(num, invm(mulm(2, y1, p), p), p)
        } else {
            mulm(subm(y2, y1, p), invm(subm(x2, x1, p), p), p)
        };
        let x3 = subm(subm(subm(mulm(lam, lam, p), self.a, p), x1, p), x2, p);
        let y3 = subm(mulm(lam, subm(x1, x3, p), p), y1, p);
        Some((x3, y3))
    }
    fn mul(&self, q: Pt, k: u64) -> Pt {
        let mut r = None;
        let mut b = q;
        let mut k = k;
        while k > 0 {
            if k & 1 == 1 {
                r = self.add(r, b);
            }
            b = self.add(b, b);
            k >>= 1;
        }
        r
    }
    fn count(&self) -> u64 {
        let p = self.p;
        let s: i64 = (0..p).map(|x| chi(self.rhs(x), p)).sum();
        (p as i64 + 1 + s) as u64
    }
    fn lift(&self, x: u64) -> Option<(u64, u64)> {
        let v = self.rhs(x);
        if chi(v, self.p) < 0 {
            return None;
        }
        let y = powm(v, (self.p + 1) / 4, self.p);
        Some((x, y))
    }
    fn random_point(&self, rng: &mut StdRng) -> Pt {
        loop {
            let x = rng.gen_range(1..self.p);
            if let Some((x, y)) = self.lift(x) {
                if y != 0 {
                    return Some(if rng.gen::<bool>() {
                        (x, y)
                    } else {
                        (x, self.p - y)
                    });
                }
            }
        }
    }
    /// Short Weierstrass (A, B) after `x = u − a/3`.
    fn short(&self) -> (u64, u64) {
        let p = self.p;
        let a3 = mulm(self.a, invm(3, p), p);
        let big_a = subm(self.b, mulm(self.a, a3, p), p);
        // B = 2a³/27 − ab/3
        let t = mulm(mulm(a3, a3, p), a3, p);
        let big_b = subm(mulm(2, t, p), mulm(self.b, a3, p), p);
        (big_a, big_b)
    }
}

/// `φ: E → E'` and `φ̂: E' → E` for `E: y² = x³ + a x² + b x`.
struct TwoIsogeny {
    e: Curve,
    e2: Curve,
}

impl TwoIsogeny {
    fn new(e: Curve) -> Self {
        let p = e.p;
        let a2 = (p - mulm(2, e.a, p)) % p;
        let b2 = subm(mulm(e.a, e.a, p), mulm(4, e.b, p), p);
        TwoIsogeny {
            e,
            e2: Curve { p, a: a2, b: b2 },
        }
    }
    fn phi(&self, q: Pt) -> Pt {
        let p = self.e.p;
        match q {
            None => None,
            Some((0, _)) => None,
            Some((x, y)) => {
                let ix2 = invm(mulm(x, x, p), p);
                let xx = mulm(mulm(y, y, p), ix2, p);
                let yy = mulm(mulm(y, subm(self.e.b, mulm(x, x, p), p), p), ix2, p);
                Some((xx, yy))
            }
        }
    }
    fn dual(&self, q: Pt) -> Pt {
        let p = self.e.p;
        match q {
            None => None,
            Some((0, _)) => None,
            Some((x, y)) => {
                let ix2 = invm(mulm(x, x, p), p);
                let xx = mulm(mulm(mulm(y, y, p), ix2, p), invm(4, p), p);
                let c = subm(
                    subm(mulm(self.e.a, self.e.a, p), mulm(4, self.e.b, p), p),
                    mulm(x, x, p),
                    p,
                );
                let yy = mulm(mulm(mulm(y, c, p), ix2, p), invm(8, p), p);
                Some((xx, yy))
            }
        }
    }
}

/// Counts of ordered m-tuples of `base` by their sum.
fn sums(c: &Curve, base: &[Pt], m: usize) -> HashMap<Pt, u64> {
    let mut cur: HashMap<Pt, u64> = HashMap::new();
    cur.insert(None, 1);
    for _ in 0..m {
        let mut nxt: HashMap<Pt, u64> = HashMap::new();
        for (&s, &k) in &cur {
            for &b in base {
                *nxt.entry(c.add(s, b)).or_insert(0) += k;
            }
        }
        cur = nxt;
    }
    cur
}

/// m-fold convolution of a weight function on a curve.
fn convolve(c: &Curve, w: &HashMap<Pt, u64>, m: usize) -> HashMap<Pt, u64> {
    let mut cur: HashMap<Pt, u64> = HashMap::new();
    cur.insert(None, 1);
    for _ in 0..m {
        let mut nxt: HashMap<Pt, u64> = HashMap::new();
        for (&s, &k) in &cur {
            for (&b, &wb) in w {
                *nxt.entry(c.add(s, b)).or_insert(0) += k * wb;
            }
        }
        cur = nxt;
    }
    cur
}

fn small_x_base(c: &Curve, xmax: u64) -> Vec<Pt> {
    let mut v = Vec::new();
    for x in 1..=xmax {
        if let Some((x, y)) = c.lift(x) {
            if y == 0 {
                continue;
            }
            v.push(Some((x, y)));
            v.push(Some((x, c.p - y)));
        }
    }
    v
}

/// `c = D_m(G)·m!/|F|^m`: distinct m-sums in the order-r subgroup per
/// unordered m-multiset of the base.
fn class_factor(dg: usize, size: usize, m: usize) -> f64 {
    let fact: f64 = (1..=m).map(|k| k as f64).product();
    dg as f64 * fact / (size as f64).powi(m as i32)
}

fn arg(args: &[String], name: &str, default: u64) -> u64 {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

/// First `(p, a, b)` with `p > pmin`, `p ≡ 3 mod 4`, `a² − 4b` a
/// non-square (so `E(F_p)[2] = {O, T}`) and `N = 4r`, `r` prime.
fn find_curve(pmin: u64) -> (Curve, u64) {
    let mut p = pmin;
    loop {
        p += 1;
        if p % 4 != 3 || !is_prime(p) {
            continue;
        }
        for a in 1..20u64 {
            for b in 1..20u64 {
                let disc = subm(mulm(a, a, p), mulm(4, b, p), p);
                if chi(disc, p) != -1 {
                    continue;
                }
                let e = Curve { p, a, b };
                let n = e.count();
                if n.is_multiple_of(4) && !n.is_multiple_of(8) && is_prime(n / 4) {
                    return (e, n);
                }
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let pmin = arg(&args, "--pmin", 1 << 20);
    let queries = arg(&args, "--queries", 5000);
    let seed = arg(&args, "--seed", 7);
    let mut rng = StdRng::seed_from_u64(seed);

    // ---- curve search: first (p, a, b) with p ≡ 3 mod 4, a² − 4b a
    // non-square (so E(F_p)[2] = {O, T}), N = 4r with r prime.
    let (e, n) = find_curve(pmin);
    let h = 4u64;
    let r = n / h;
    let iso = TwoIsogeny::new(e);
    let e2 = iso.e2;
    let n2 = e2.count();
    let t: Pt = Some((0, 0));
    let slug = |c: &Curve, order: u64| {
        let (sa, sb) = c.short();
        curve_id::prime(
            &BigUint::from(c.p),
            &BigUint::from(sa),
            &BigUint::from(sb),
            &BigUint::from(order),
        )
        .map(|id| (id.slug, id.icv1))
        .unwrap_or_default()
    };
    let (slug_e, icv_e) = slug(&e, n);
    let (slug_e2, icv_e2) = slug(&e2, n2);

    // ---- check 1: maps, kernel, duals, scalar law.
    let samples = 2000;
    let mut ok_on = 0;
    let mut ok_hom = 0;
    let mut ok_dual = 0;
    let mut ok_dual2 = 0;
    for _ in 0..samples {
        let s = e.random_point(&mut rng);
        let u = e.random_point(&mut rng);
        let fs = iso.phi(s);
        if e2.on(fs) && e.on(iso.dual(fs)) {
            ok_on += 1;
        }
        if iso.phi(e.add(s, u)) == e2.add(fs, iso.phi(u)) {
            ok_hom += 1;
        }
        if iso.dual(fs) == e.mul(s, 2) {
            ok_dual += 1;
        }
        let s2 = e2.random_point(&mut rng);
        if iso.phi(iso.dual(s2)) == e2.mul(s2, 2) {
            ok_dual2 += 1;
        }
    }
    let kernel_ok = iso.phi(t).is_none() && iso.phi(None).is_none();
    // generator of G = E(F_p)[r]
    let gen = loop {
        let q = e.mul(e.random_point(&mut rng), h);
        if q.is_some() {
            break q;
        }
    };
    assert!(e.mul(gen, r).is_none());
    let gen2 = iso.phi(gen); // P' := φ(P): λ_φ = 1 by normalisation
    let lam_dual_ok = iso.dual(gen2) == e.mul(gen, 2);
    let n_order_ok = n2 == n;

    // ---- check 2: AGM oracle simulation.
    let x_secret = rng.gen_range(1..r);
    let chal = e.mul(gen, x_secret);
    let mut evals_sim = 0u64;
    let chal2 = iso.phi(chal); // the only per-instance evaluation
    evals_sim += 1;
    // lists: (point, u, v) with point = u·P + v·X on E, u·P' + v·X' on E'.
    let mut l1: Vec<(Pt, u64, u64)> = vec![(gen, 1, 0), (chal, 0, 1)];
    let mut l2: Vec<(Pt, u64, u64)> = vec![(gen2, 1, 0)];
    let mut agree = 0u64;
    let mut q_phi = 0u64;
    let mut q_dual = 0u64;
    let mut sim_group_ops_scalar_mults = 0u64;
    for _ in 0..queries {
        let use_phi = rng.gen::<bool>();
        let (src, src_curve) = if use_phi { (&l1, &e) } else { (&l2, &e2) };
        let k = src.len().min(4);
        let mut z: Pt = None;
        let (mut u, mut v) = (0u64, 0u64);
        for _ in 0..k {
            let idx = rng.gen_range(0..src.len());
            let c = rng.gen_range(0..r);
            let (pt, pu, pv) = src[idx];
            z = src_curve.add(z, src_curve.mul(pt, c));
            u = (u + mulm(c, pu, r)) % r;
            v = (v + mulm(c, pv, r)) % r;
        }
        if use_phi {
            q_phi += 1;
            let truth = iso.phi(z);
            // λ_φ = 1: answer u·P' + v·X'
            let sim = e2.add(e2.mul(gen2, u), e2.mul(chal2, v));
            sim_group_ops_scalar_mults += 2;
            if sim == truth {
                agree += 1;
            }
            l2.push((truth, u, v));
        } else {
            q_dual += 1;
            let truth = iso.dual(z);
            // λ_φ̂ = deg φ / λ_φ = 2: answer 2u·P + 2v·X
            let (u2, v2) = (mulm(2, u, r), mulm(2, v, r));
            let sim = e.add(e.mul(gen, u2), e.mul(chal, v2));
            sim_group_ops_scalar_mults += 2;
            if sim == truth {
                agree += 1;
            }
            l1.push((truth, u2, v2));
        }
    }

    // ---- check 4: cofactor split for unrestricted queries.
    let s_inv = invm(h % r, r); // h^{-1} mod r
    let e_g = (s_inv as u128 * h as u128 % n as u128) as u64; // idempotent onto G
    let mut tors: Vec<Pt> = Vec::new(); // E(F_p)[4]
    while tors.len() < 4 {
        let q = e.mul(e.random_point(&mut rng), r);
        let mut c = q;
        for _ in 0..4 {
            if !tors.contains(&c) {
                tors.push(c);
            }
            c = e.add(c, q);
        }
    }
    let phi_tab: HashMap<Pt, Pt> = tors.iter().map(|&q| (q, iso.phi(q))).collect();
    let mut split_ok = 0;
    let split_n = 2000;
    for _ in 0..split_n {
        let z = e.random_point(&mut rng);
        let zg = e.mul(z, e_g);
        let zh = e.add(z, e.neg(zg));
        if e.mul(zg, r).is_none()
            && e.mul(zh, h).is_none()
            && iso.phi(z) == e2.add(iso.phi(zg), phi_tab[&zh])
        {
            split_ok += 1;
        }
    }

    // ---- check 3 and 5: pushforward identities and census.
    let in_g = |c: &Curve, q: Pt| c.mul(q, r).is_none();
    let mut push_rows = Vec::new();
    let mut census_rows = Vec::new();
    for &(m, xmax) in &[(2usize, 300u64), (3usize, 40u64)] {
        let plain = small_x_base(&e, xmax);
        // φ-injective: drop P when P + T is already in.
        let mut inj: Vec<Pt> = Vec::new();
        for &q in &plain {
            if !inj.contains(&e.add(q, t)) {
                inj.push(q);
            }
        }
        // T-saturated from half the x-range, so sizes are comparable.
        let half = small_x_base(&e, xmax / 2);
        let mut sat: Vec<Pt> = Vec::new();
        for &q in &half {
            for c in [q, e.add(q, t)] {
                if !sat.contains(&c) {
                    sat.push(c);
                }
            }
        }
        for (name, base) in [("plain", &plain), ("injective", &inj), ("saturated", &sat)] {
            let d = sums(&e, base, m);
            // pushforward of d
            let mut push: HashMap<Pt, u64> = HashMap::new();
            for (&rr, &k) in &d {
                *push.entry(iso.phi(rr)).or_insert(0) += k;
            }
            // w_F and its m-fold convolution on E'
            let mut w: HashMap<Pt, u64> = HashMap::new();
            for &q in base.iter() {
                *w.entry(iso.phi(q)).or_insert(0) += 1;
            }
            let conv = convolve(&e2, &w, m);
            let mut keys: Vec<Pt> = push.keys().chain(conv.keys()).copied().collect();
            keys.sort();
            keys.dedup();
            let fiber_sum_mismatch = keys
                .iter()
                .filter(|q| push.get(q).copied().unwrap_or(0) != conv.get(q).copied().unwrap_or(0))
                .count();
            // image base φ(F) (as a set) and its sums
            let mut img: Vec<Pt> = w.keys().copied().collect();
            img.sort();
            let dimg = sums(&e2, &img, m);
            // pointwise laws over the support of d (and of φ*d_img)
            let mut rset: Vec<Pt> = d.keys().copied().collect();
            // add the full fibers (R + T) of every support point
            let extra: Vec<Pt> = rset.iter().map(|&q| e.add(q, t)).collect();
            rset.extend(extra);
            rset.sort();
            rset.dedup();
            let mut sat_law_ok = 0u64;
            let mut inj_law_ok = 0u64;
            for &rr in &rset {
                let lhs = d.get(&rr).copied().unwrap_or(0);
                let rhs = dimg.get(&iso.phi(rr)).copied().unwrap_or(0);
                if lhs == (1u64 << (m - 1)) * rhs {
                    sat_law_ok += 1;
                }
                let fsum = lhs + d.get(&e.add(rr, t)).copied().unwrap_or(0);
                if fsum == rhs {
                    inj_law_ok += 1;
                }
            }
            // m-sums of an image base never leave φ(E(F_p)).  Descent
            // criterion: Q' ∈ φ(E(F_p)) iff Q' = O or x(Q') is a non-zero
            // square ((0,0) is not in the image because a² − 4b is a
            // non-square).
            let in_image = |q: &Pt| match q {
                None => true,
                Some((x, _)) => *x != 0 && chi(*x, e2.p) == 1,
            };
            let base_outside = img.iter().filter(|q| !in_image(q)).count();
            let img_outside = dimg.keys().filter(|q| !in_image(q)).count();
            push_rows.push(json!({
                "m": m, "base": name, "size": base.len(), "image_size": img.len(),
                "phi_injective_on_base": img.len() == base.len(),
                "support_points_checked": keys.len(),
                "fiber_sum_identity_mismatches": fiber_sum_mismatch,
                "points_R_checked": rset.len(),
                "saturated_law_holds_at": sat_law_ok,
                "injective_law_holds_at": inj_law_ok,
                "image_base_points_outside_image": base_outside,
                "image_base_sums_outside_image": img_outside,
            }));
            // census: distinct m-sums landing in G (or G')
            let dg = d.keys().filter(|q| in_g(&e, **q)).count();
            let dgi = dimg.keys().filter(|q| in_g(&e2, **q)).count();
            census_rows.push(json!({
                "m": m, "base": name, "curve": "E", "size": base.len(),
                "distinct_sums": d.len(), "distinct_sums_in_G": dg,
                "class_factor": class_factor(dg, base.len(), m),
            }));
            census_rows.push(json!({
                "m": m, "base": format!("phi({name})"), "curve": "E'", "size": img.len(),
                "distinct_sums": dimg.len(), "distinct_sums_in_G": dgi,
                "class_factor": class_factor(dgi, img.len(), m),
            }));
        }
        // plain base on E' of size close to |φ(sat)| for the image law
        let target = sat.len() / 2;
        let mut xm = 1u64;
        let mut plain2 = small_x_base(&e2, xm);
        while plain2.len() < target {
            xm += 1;
            plain2 = small_x_base(&e2, xm);
        }
        let d2 = sums(&e2, &plain2, m);
        let dg2 = d2.keys().filter(|q| in_g(&e2, **q)).count();
        census_rows.push(json!({
            "m": m, "base": "plain'", "curve": "E'", "size": plain2.len(),
            "distinct_sums": d2.len(), "distinct_sums_in_G": dg2,
            "class_factor": class_factor(dg2, plain2.len(), m),
        }));
    }
    // Derived ratios of class factors, against the heuristic predictions
    // (image law: ×i with i = 2; saturation law: ÷2^{m−1}) and the exact
    // prediction D_m(G)[saturated on E] = D_m(G')[φ(saturated) on E'].
    let cf = |m: usize, base: &str, curve: &str| -> f64 {
        census_rows
            .iter()
            .find(|v| v["m"] == m && v["base"] == base && v["curve"] == curve)
            .and_then(|v| v["class_factor"].as_f64())
            .unwrap_or(f64::NAN)
    };
    let dgv = |m: usize, base: &str, curve: &str| -> u64 {
        census_rows
            .iter()
            .find(|v| v["m"] == m && v["base"] == base && v["curve"] == curve)
            .and_then(|v| v["distinct_sums_in_G"].as_u64())
            .unwrap_or(0)
    };
    let mut ratios = Vec::new();
    for m in [2usize, 3] {
        ratios.push(json!({
            "m": m,
            "image_law_phi_plain_over_plain_E": cf(m, "phi(plain)", "E'") / cf(m, "plain", "E"),
            "image_law_phi_sat_over_plain_E2": cf(m, "phi(saturated)", "E'") / cf(m, "plain'", "E'"),
            "image_law_predicted": 2.0,
            "saturation_law_sat_over_plain_E": cf(m, "saturated", "E") / cf(m, "plain", "E"),
            "saturation_law_predicted": 1.0 / f64::from(1u32 << (m - 1)),
            "exact_DG_saturated_E": dgv(m, "saturated", "E"),
            "exact_DG_phi_saturated_E2": dgv(m, "phi(saturated)", "E'"),
        }));
    }

    let out = json!({
        "schema": "r4_16_7_agm_pushforward/v1",
        "seed": seed,
        "curve_E": {"p": e.p, "a": e.a, "b": e.b, "model": "y^2 = x^3 + a x^2 + b x",
                    "order": n, "h": h, "r": r, "icv1": icv_e, "slug": slug_e},
        "curve_E2": {"p": e2.p, "a": e2.a, "b": e2.b, "model": "y^2 = x^3 + a x^2 + b x",
                     "order": n2, "icv1": icv_e2, "slug": slug_e2},
        "check1_maps": {
            "samples": samples, "on_curve": ok_on, "homomorphism": ok_hom,
            "dual_after_phi_is_2": ok_dual, "phi_after_dual_is_2": ok_dual2,
            "kernel_is_O_T": kernel_ok, "orders_equal": n_order_ok,
            "lambda_dual_equals_deg_over_lambda_phi": lam_dual_ok,
        },
        "check2_agm_simulation": {
            "queries": queries, "phi_queries": q_phi, "dual_queries": q_dual,
            "answers_agree": agree,
            "simulator_isogeny_evaluations_per_instance": evals_sim,
            "instance_independent_evaluations": 1,
            "naive_simulator_isogeny_evaluations": queries,
            "simulator_scalar_multiplications": sim_group_ops_scalar_mults,
        },
        "check3_pushforward": push_rows,
        "check4_cofactor_split": {"points": split_n, "holds": split_ok, "e_G": e_g},
        "check5_census": census_rows,
        "check5_ratios_derived": ratios,
    });
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
}
