//! Counted Cantor addition and doubling on genus-`g` Jacobians over `F_p`, `p`
//! odd, against the cost model's `10 · g²` field multiplications per group
//! operation (the former `jacobian_group_op_mults_per_genus_squared` in
//! `docs/curves/jacobian-transfers.json`, since replaced by this run's counts
//! as `jacobian_group_op_mults_measured`).
//!
//! ```text
//! cargo run --release --example jacobian_op_count -- run \
//!     --out research/jacobian_op_count_20261014/rows.csv --seed 20261014 \
//!     --trials 200 --genera 2,3,4,5,6,8,10,12,14,16,20,24,32 --primes 1009,65521
//! cargo run --release --example jacobian_op_count -- summary \
//!     --rows research/jacobian_op_count_20261014/rows.csv \
//!     --jacobian docs/curves/jacobian-transfers.json \
//!     --out research/jacobian_op_count_20261014/summary.csv
//! ```
//!
//! The counted implementation is textbook Cantor (extended Euclid for the two
//! gcds, composition, then reduction) over `u64` residues, with one counter per
//! field multiplication, addition and inversion.  Every result is compared with
//! the repository's `MumfordDivisorP::add`, so a wrong count cannot hide behind
//! a wrong group law.  Protocol: `research/jacobian_op_count_20261014/PROTOCOL.md`.

use crypto_lib::prime_hyperelliptic::{FpPoly, HyperellipticCurveP, MumfordDivisorP};
use num_bigint::BigUint;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::{fs, io::Write};

/// Cost convention for an inversion, in multiplications (extended Euclid on a
/// residue of this size).  Reported counts keep inversions separate.
const INV_WEIGHT: f64 = 8.0;

#[derive(Default, Clone, Copy, Debug)]
struct Cnt {
    mul: u64,
    add: u64,
    inv: u64,
    reduce: u64,
}

struct Fp {
    p: u64,
}

type Poly = Vec<u64>; // little-endian, no trailing zeros; zero is empty

impl Fp {
    fn mul(&self, c: &mut Cnt, a: u64, b: u64) -> u64 {
        c.mul += 1;
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    fn add(&self, c: &mut Cnt, a: u64, b: u64) -> u64 {
        c.add += 1;
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    fn sub(&self, c: &mut Cnt, a: u64, b: u64) -> u64 {
        c.add += 1;
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    fn inv(&self, c: &mut Cnt, a: u64) -> u64 {
        c.inv += 1;
        let (mut base, mut e, mut acc) = (a % self.p, self.p - 2, 1u64);
        while e > 0 {
            if e & 1 == 1 {
                acc = ((acc as u128 * base as u128) % self.p as u128) as u64;
            }
            base = ((base as u128 * base as u128) % self.p as u128) as u64;
            e >>= 1;
        }
        acc
    }
}

fn trim(v: &mut Poly) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

fn deg(v: &Poly) -> Option<usize> {
    v.len().checked_sub(1)
}

fn padd(f: &Fp, c: &mut Cnt, a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        out.push(match (a.get(i), b.get(i)) {
            (Some(&x), Some(&y)) => f.add(c, x, y),
            (Some(&x), None) | (None, Some(&x)) => x,
            (None, None) => 0,
        });
    }
    trim(&mut out);
    out
}

fn psub(f: &Fp, c: &mut Cnt, a: &Poly, b: &Poly) -> Poly {
    let n = a.len().max(b.len());
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        out.push(match (a.get(i), b.get(i)) {
            (Some(&x), Some(&y)) => f.sub(c, x, y),
            (Some(&x), None) => x,
            (None, Some(&y)) => f.p - y,
            (None, None) => 0,
        });
    }
    trim(&mut out);
    out
}

fn pneg(f: &Fp, a: &Poly) -> Poly {
    a.iter()
        .map(|&x| if x == 0 { 0 } else { f.p - x })
        .collect()
}

fn pmul(f: &Fp, c: &mut Cnt, a: &Poly, b: &Poly) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            let prod = f.mul(c, x, y);
            out[i + j] = f.add(c, out[i + j], prod);
        }
    }
    trim(&mut out);
    out
}

fn pscale(f: &Fp, c: &mut Cnt, a: &Poly, s: u64) -> Poly {
    let mut out: Poly = a.iter().map(|&x| f.mul(c, x, s)).collect();
    trim(&mut out);
    out
}

fn pdivrem(f: &Fp, c: &mut Cnt, a: &Poly, b: &Poly) -> (Poly, Poly) {
    let db = deg(b).expect("division by zero polynomial");
    let lead = b[db];
    let inv_lead = if lead == 1 { 1 } else { f.inv(c, lead) };
    let mut r = a.clone();
    let mut q = vec![0u64; a.len().saturating_sub(db).max(1)];
    while let Some(dr) = deg(&r) {
        if dr < db {
            break;
        }
        let shift = dr - db;
        let coef = if lead == 1 {
            r[dr]
        } else {
            f.mul(c, r[dr], inv_lead)
        };
        q[shift] = coef;
        for i in 0..=db {
            let t = f.mul(c, coef, b[i]);
            r[i + shift] = f.sub(c, r[i + shift], t);
        }
        r[dr] = 0;
        trim(&mut r);
    }
    trim(&mut q);
    (q, r)
}

fn pmonic(f: &Fp, c: &mut Cnt, a: &Poly) -> Poly {
    match a.last() {
        Some(&l) if l != 1 => {
            let inv = f.inv(c, l);
            pscale(f, c, a, inv)
        }
        _ => a.clone(),
    }
}

/// `(g, s, t)` with `g = s a + t b`, `g` monic.
fn pxgcd(f: &Fp, c: &mut Cnt, a: &Poly, b: &Poly) -> (Poly, Poly, Poly) {
    let (mut r0, mut r1) = (a.clone(), b.clone());
    let (mut s0, mut s1): (Poly, Poly) = (vec![1], Vec::new());
    let (mut t0, mut t1): (Poly, Poly) = (Vec::new(), vec![1]);
    while !r1.is_empty() {
        let (q, r) = pdivrem(f, c, &r0, &r1);
        let qs = pmul(f, c, &q, &s1);
        let s2 = psub(f, c, &s0, &qs);
        let qt = pmul(f, c, &q, &t1);
        let t2 = psub(f, c, &t0, &qt);
        r0 = r1;
        r1 = r;
        s0 = s1;
        s1 = s2;
        t0 = t1;
        t1 = t2;
    }
    match r0.last() {
        Some(&l) if l != 1 => {
            let inv = f.inv(c, l);
            (
                pscale(f, c, &r0, inv),
                pscale(f, c, &s0, inv),
                pscale(f, c, &t0, inv),
            )
        }
        _ => (r0, s0, t0),
    }
}

/// Cantor composition then reduction of two reduced divisors, counted.
fn cantor_add(
    f: &Fp,
    c: &mut Cnt,
    curve_f: &Poly,
    genus: usize,
    d1: &(Poly, Poly),
    d2: &(Poly, Poly),
) -> (Poly, Poly) {
    let (u1, v1) = d1;
    let (u2, v2) = d2;
    let (g1, e1, e2) = pxgcd(f, c, u1, u2);
    let w = padd(f, c, v1, v2);
    let (d, c1, c2) = pxgcd(f, c, &g1, &w);
    let s1 = pmul(f, c, &c1, &e1);
    let s2 = pmul(f, c, &c1, &e2);
    let s3 = c2;
    let mut u = pmul(f, c, u1, u2);
    let mut sum = {
        let s1u1 = pmul(f, c, &s1, u1);
        let t1 = pmul(f, c, &s1u1, v2);
        let s2u2 = pmul(f, c, &s2, u2);
        let t2 = pmul(f, c, &s2u2, v1);
        let v1v2 = pmul(f, c, v1, v2);
        let v1v2f = padd(f, c, &v1v2, curve_f);
        let t3 = pmul(f, c, &s3, &v1v2f);
        let t12 = padd(f, c, &t1, &t2);
        padd(f, c, &t12, &t3)
    };
    if d.len() > 1 {
        let u_over_d = pdivrem(f, c, &u, &d).0;
        u = pdivrem(f, c, &u_over_d, &d).0;
        sum = pdivrem(f, c, &sum, &d).0;
    }
    let mut u = pmonic(f, c, &u);
    let mut v = pdivrem(f, c, &sum, &u).1;
    while deg(&u).is_some_and(|d| d > genus) {
        c.reduce += 1;
        let vv = pmul(f, c, &v, &v);
        let tmp = psub(f, c, curve_f, &vv);
        let quotient = pdivrem(f, c, &tmp, &u).0;
        let uq = pmonic(f, c, &quotient);
        v = pdivrem(f, c, &pneg(f, &v), &uq).1;
        u = uq;
    }
    (u, v)
}

fn to_poly(p: &FpPoly) -> Poly {
    match p.degree() {
        None => Vec::new(),
        Some(d) => (0..=d)
            .map(|i| p.coeff(i).to_u64_digits().first().copied().unwrap_or(0))
            .collect(),
    }
}

fn big(x: u64) -> BigUint {
    BigUint::from(x)
}

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter()
        .position(|a| a == k)
        .map(|i| args[i + 1].clone())
}

fn list(s: &str) -> Vec<u64> {
    s.split(',')
        .filter(|x| !x.is_empty())
        .map(|x| x.parse().expect("integer list"))
        .collect()
}

const ROW_HEADER: &str = "g,p,op,trial,mul,add,inv,reduce_steps,matches_repo";

/// One random genus-`g` curve over `F_p` (smooth, monic of degree `2g+1`).
fn random_curve(rng: &mut StdRng, p: u64, g: usize) -> HyperellipticCurveP {
    loop {
        let mut coeffs: Vec<BigUint> = (0..=2 * g).map(|_| big(rng.gen_range(0..p))).collect();
        coeffs.push(big(1));
        let f = FpPoly::from_coeffs(coeffs, big(p));
        if let Ok(c) = HyperellipticCurveP::try_new(big(p), f, g as u32) {
            return c;
        }
    }
}

/// A divisor of degree `g` (sum of `g` random points), generic with high
/// probability; resampled until `deg u = g`.
fn random_divisor(
    rng: &mut StdRng,
    curve: &HyperellipticCurveP,
    sqrt: &[Option<u64>],
    p: u64,
    g: usize,
) -> MumfordDivisorP {
    loop {
        let mut d = MumfordDivisorP::identity(big(p));
        let mut count = 0;
        while count < g {
            let x = rng.gen_range(0..p);
            let rhs = curve.f.eval(&big(x));
            let rhs_u = rhs.to_u64_digits().first().copied().unwrap_or(0);
            if let Some(y) = sqrt[rhs_u as usize] {
                if let Some(pt) = MumfordDivisorP::from_point(curve, &big(x), &big(y)) {
                    d = d.add(&pt, curve);
                    count += 1;
                }
            }
        }
        if d.u.degree() == Some(g) {
            return d;
        }
    }
}

fn run(args: &[String]) {
    let out = arg(args, "--out").expect("--out");
    let seed: u64 = arg(args, "--seed").map_or(20261014, |s| s.parse().unwrap());
    let trials: usize = arg(args, "--trials").map_or(200, |s| s.parse().unwrap());
    let genera = list(&arg(args, "--genera").expect("--genera"));
    let primes = list(&arg(args, "--primes").expect("--primes"));
    let mut rng = StdRng::seed_from_u64(seed);
    let mut rows = vec![ROW_HEADER.to_string()];
    let mut mismatches = 0usize;
    for &p in &primes {
        let mut sqrt = vec![None; p as usize];
        for y in 0..p {
            sqrt[((y as u128 * y as u128) % p as u128) as usize] = Some(y);
        }
        let fp = Fp { p };
        for &g in &genera {
            let g = g as usize;
            let curve = random_curve(&mut rng, p, g);
            let curve_f = to_poly(&curve.f);
            for trial in 0..trials {
                let a = random_divisor(&mut rng, &curve, &sqrt, p, g);
                let b = random_divisor(&mut rng, &curve, &sqrt, p, g);
                for (op, x, y) in [("add", &a, &b), ("double", &a, &a)] {
                    let mut c = Cnt::default();
                    let mine = cantor_add(
                        &fp,
                        &mut c,
                        &curve_f,
                        g,
                        &(to_poly(&x.u), to_poly(&x.v)),
                        &(to_poly(&y.u), to_poly(&y.v)),
                    );
                    let repo = x.add(y, &curve);
                    let ok = mine.0 == to_poly(&repo.u) && mine.1 == to_poly(&repo.v);
                    if !ok {
                        mismatches += 1;
                    }
                    rows.push(format!(
                        "{g},{p},{op},{trial},{},{},{},{},{ok}",
                        c.mul, c.add, c.inv, c.reduce
                    ));
                }
            }
            eprintln!("p = {p}, g = {g}: {trials} trials done");
        }
    }
    let mut f = fs::File::create(&out).expect("out");
    for r in &rows {
        writeln!(f, "{r}").unwrap();
    }
    println!(
        "wrote {} rows to {out}; mismatches against the repository group law: {mismatches}",
        rows.len() - 1
    );
    if mismatches > 0 {
        std::process::exit(1);
    }
}

// ── Summary ─────────────────────────────────────────────────────────────

struct Row {
    g: u64,
    p: u64,
    op: String,
    mul: f64,
    inv: f64,
}

fn read_rows(path: &str) -> Vec<Row> {
    fs::read_to_string(path)
        .expect("rows")
        .lines()
        .skip(1)
        .map(|l| {
            let t: Vec<&str> = l.split(',').collect();
            Row {
                g: t[0].parse().unwrap(),
                p: t[1].parse().unwrap(),
                op: t[2].to_string(),
                mul: t[4].parse().unwrap(),
                inv: t[6].parse().unwrap(),
            }
        })
        .collect()
}

/// Least-squares line `y = a + e x`; returns `(a, e)`.
fn fit(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let (mx, my) = (xs.iter().sum::<f64>() / n, ys.iter().sum::<f64>() / n);
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let e = sxy / sxx;
    (my - e * mx, e)
}

fn summary(args: &[String]) {
    let rows = read_rows(&arg(args, "--rows").expect("--rows"));
    let out = arg(args, "--out").expect("--out");
    let json: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(arg(args, "--jacobian").expect("--jacobian")).unwrap(),
    )
    .unwrap();
    let mut lines = vec!["kind,op,p,g,mean_mul,mean_inv,mean_total,model_10g2,log2_ratio_to_model,fit_log2_c,fit_exponent,exponent_ci_lo,exponent_ci_hi,local_exponent".to_string()];
    let mut rng = StdRng::seed_from_u64(7);
    let mut fits: Vec<(String, u64, f64, f64)> = Vec::new(); // (op, p, c, e) for the sensitivity step
    let mut ops: Vec<(String, u64)> = rows.iter().map(|r| (r.op.clone(), r.p)).collect();
    ops.sort();
    ops.dedup();
    for (op, p) in &ops {
        let sel: Vec<&Row> = rows.iter().filter(|r| &r.op == op && r.p == *p).collect();
        let mut gs: Vec<u64> = sel.iter().map(|r| r.g).collect();
        gs.sort();
        gs.dedup();
        let by_g: Vec<Vec<f64>> = gs
            .iter()
            .map(|g| {
                sel.iter()
                    .filter(|r| r.g == *g)
                    .map(|r| r.mul + INV_WEIGHT * r.inv)
                    .collect()
            })
            .collect();
        let means: Vec<f64> = by_g
            .iter()
            .map(|v| v.iter().sum::<f64>() / v.len() as f64)
            .collect();
        let (xs, ys): (Vec<f64>, Vec<f64>) = gs
            .iter()
            .zip(&means)
            .map(|(g, m)| ((*g as f64).log2(), m.log2()))
            .unzip();
        let (a, e) = fit(&xs, &ys);
        // bootstrap the exponent: resample trials within each genus
        let mut es: Vec<f64> = Vec::new();
        for _ in 0..2000 {
            let bm: Vec<f64> = by_g
                .iter()
                .map(|v| {
                    (0..v.len())
                        .map(|_| v[rng.gen_range(0..v.len())])
                        .sum::<f64>()
                        / v.len() as f64
                })
                .collect();
            let by: Vec<f64> = bm.iter().map(|m| m.log2()).collect();
            es.push(fit(&xs, &by).1);
        }
        es.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let (lo, hi) = (
            es[(es.len() as f64 * 0.025) as usize],
            es[(es.len() as f64 * 0.975) as usize],
        );
        for (i, g) in gs.iter().enumerate() {
            let mean_mul: f64 =
                sel.iter().filter(|r| r.g == *g).map(|r| r.mul).sum::<f64>() / by_g[i].len() as f64;
            let mean_inv: f64 =
                sel.iter().filter(|r| r.g == *g).map(|r| r.inv).sum::<f64>() / by_g[i].len() as f64;
            let model = 10.0 * (*g as f64) * (*g as f64);
            // exponent between this genus and the previous one
            let local = if i == 0 {
                String::new()
            } else {
                let num = means[i].log2() - means[i - 1].log2();
                let den = (*g as f64).log2() - (gs[i - 1] as f64).log2();
                format!("{:.3}", num / den)
            };
            lines.push(format!(
                "genus,{op},{p},{g},{mean_mul:.1},{mean_inv:.2},{:.1},{model:.0},{:.3},{a:.3},{e:.3},{lo:.3},{hi:.3},{local}",
                means[i],
                (means[i] / model).log2()
            ));
        }
        fits.push((op.clone(), *p, a, e));
    }
    // Sensitivity: how many advantage verdicts sit within the op-cost error?
    let curves = json["curves"].as_array().expect("curves");
    let mut slacks: Vec<(f64, f64)> = Vec::new(); // (slack bits over the margin, log2 genus)
    let mut marginal = 0usize;
    for c in curves {
        let best = &c["routes"][c["best_route"].as_u64().unwrap() as usize];
        match c["verdict"].as_str() {
            Some("formula_advantage") => {
                let slack = best["cost"]["ratio_log2"].as_f64().unwrap() - 8.0;
                slacks.push((
                    slack,
                    best["destination"]["log2_genus"].as_f64().unwrap_or(0.0),
                ));
            }
            Some("formula_marginal") => marginal += 1,
            _ => {}
        }
    }
    let add_fit = fits
        .iter()
        .find(|f| f.0 == "add" && f.1 == ops.iter().map(|o| o.1).max().unwrap())
        .expect("add fit");
    let mut flipped = 0usize;
    for (slack, lg) in &slacks {
        let g = (2f64).powf(*lg);
        let measured = add_fit.2 + add_fit.3 * g.log2();
        let delta = measured - (10.0 * g * g).log2();
        if delta > 0.0 && *slack <= delta {
            flipped += 1;
        }
    }
    let within = |b: f64| slacks.iter().filter(|s| s.0 <= b).count();
    let detail = format!(
        "advantage_curves={};marginal_curves={};slack_le_1bit={};slack_le_2bits={};slack_le_3bits={};would_drop_below_margin_if_measured_cost_applied_to_the_whole_route={}",
        slacks.len(),
        marginal,
        within(1.0),
        within(2.0),
        within(3.0),
        flipped
    );
    lines.push(format!(
        "sensitivity,,,,,,,,,{:.3},{:.3},,{detail}",
        add_fit.2, add_fit.3
    ));
    fs::write(&out, lines.join("\n") + "\n").expect("summary");
    println!("wrote {out}");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("run") => run(&args),
        Some("summary") => summary(&args),
        _ => {
            eprintln!("usage: jacobian_op_count run --out F --seed S --trials N --genera 2,3,.. --primes P,.. | summary --rows F --jacobian F --out F");
            std::process::exit(2);
        }
    }
}
