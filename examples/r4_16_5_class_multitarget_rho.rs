//! Round-4 lead 16.5: class-wide multi-target rho, measured.
//!
//! Builds a toy `F_p`-isogeny class (prime-field short Weierstrass curves of
//! order `h·r`, `r` prime, `6 | h`) by a breadth-first walk over rational
//! 2- and 3-isogenies (kernels from the rational roots of the 2- and
//! 3-division polynomials, Vélu's formulas), with an explicit dual for every
//! edge.  Generators are transported consistently: `P_i = ψ_i(P_0)` along
//! the forward path, so the dual path sends `P_i` to `[d_i]P_0` with
//! `d_i = deg ψ_i`.
//!
//! Three arms, all on the same distinguished-point rho with a
//! target-independent 32-adding walk `X ← X + α_j P`:
//!
//! * **A** — `L` independent single-target runs, target `i` on curve `E_i`
//!   (`i = 1..L`), a fresh DP table per run.
//! * **B** — the `L` targets transported to `E_0` by the dual isogeny paths
//!   and solved by one multi-target (Kuhn–Struik) rho that keeps its DP table
//!   across targets; each recovered log is mapped back by `d_i^{-1}` and
//!   verified on its original curve.
//! * **C** — the single-curve control: `L` targets planted directly on `E_0`
//!   and solved by the same multi-target rho.
//!
//! Counted unit: group operations (affine additions/doublings) in the order-`r`
//! subgroup.  Transport is counted separately in field multiplications and
//! inversions.  Construction uses `crypto_lib`'s isogeny-walk field, root
//! finder and odd-degree Vélu codomain (as an independent cross-check of the
//! 3-isogeny codomains), and `curve_id` for the base curve's ICV1 slug.
//!
//!     cargo run --release --example r4_16_5_class_multitarget_rho -- \
//!         --rbits 32 --trials 20 --dpbits 8 --seed 1 --out out.json

use crypto_lib::cryptanalysis::curve_id;
use crypto_lib::cryptanalysis::isogeny_walk::curve::Model;
use crypto_lib::cryptanalysis::isogeny_walk::field::{is_probable_prime, Field};
use crypto_lib::cryptanalysis::isogeny_walk::kernel::velu_codomain;
use crypto_lib::cryptanalysis::isogeny_walk::poly;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::collections::{HashMap, HashSet};
use std::time::Instant;

// ---------------------------------------------------------------- field ----

#[derive(Clone, Copy)]
struct Fp {
    p: u64,
}

impl Fp {
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    fn mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    fn inv(&self, a: u64) -> u64 {
        inv_mod(a, self.p)
    }
    fn pow(&self, a: u64, e: u64) -> u64 {
        pow_mod(a, e, self.p)
    }
}

fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn pow_mod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1 % m;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mul_mod(r, a, m);
        }
        a = mul_mod(a, a, m);
        e >>= 1;
    }
    r
}

fn inv_mod(a: u64, m: u64) -> u64 {
    let (mut t, mut nt) = (0i128, 1i128);
    let (mut rr, mut nr) = (m as i128, (a % m) as i128);
    while nr != 0 {
        let q = rr / nr;
        (t, nt) = (nt, t - q * nt);
        (rr, nr) = (nr, rr - q * nr);
    }
    assert!(rr == 1, "not invertible");
    if t < 0 {
        t += m as i128;
    }
    t as u64
}

// ---------------------------------------------------------------- curve ----

type Pt = Option<(u64, u64)>;

#[derive(Clone, Copy, Debug)]
struct Curve {
    a: u64,
    b: u64,
}

impl Curve {
    fn rhs(&self, f: &Fp, x: u64) -> u64 {
        f.add(f.mul(f.add(f.mul(x, x), self.a), x), self.b)
    }
    fn on(&self, f: &Fp, pt: Pt) -> bool {
        match pt {
            None => true,
            Some((x, y)) => f.mul(y, y) == self.rhs(f, x),
        }
    }
    fn j(&self, f: &Fp) -> u64 {
        let a3 = f.mul(f.mul(self.a, self.a), self.a);
        let d = f.add(f.mul(4, a3), f.mul(27, f.mul(self.b, self.b)));
        f.mul(f.mul(1728 * 4, a3), f.inv(d))
    }
    fn neg(&self, f: &Fp, p: Pt) -> Pt {
        p.map(|(x, y)| (x, f.neg(y)))
    }
    /// One group operation (counted by the caller).
    fn add(&self, f: &Fp, p: Pt, q: Pt) -> Pt {
        let (x1, y1) = match p {
            None => return q,
            Some(v) => v,
        };
        let (x2, y2) = match q {
            None => return p,
            Some(v) => v,
        };
        let lam = if x1 == x2 {
            if y1 != y2 || y1 == 0 {
                return None;
            }
            let num = f.add(f.mul(3, f.mul(x1, x1)), self.a);
            f.mul(num, f.inv(f.add(y1, y1)))
        } else {
            f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1)))
        };
        let x3 = f.sub(f.sub(f.mul(lam, lam), x1), x2);
        let y3 = f.sub(f.mul(lam, f.sub(x1, x3)), y1);
        Some((x3, y3))
    }
    /// Double-and-add; adds the group operations used to `ops`.
    fn smul(&self, f: &Fp, p: Pt, k: u64, ops: &mut u64) -> Pt {
        let mut acc: Pt = None;
        if k == 0 {
            return None;
        }
        let top = 63 - k.leading_zeros();
        for i in (0..=top).rev() {
            if acc.is_some() {
                acc = self.add(f, acc, acc);
                *ops += 1;
            }
            if (k >> i) & 1 == 1 {
                if acc.is_some() {
                    *ops += 1;
                }
                acc = self.add(f, acc, p);
            }
        }
        acc
    }
    fn smul_big(&self, f: &Fp, p: Pt, k: u128) -> Pt {
        let mut acc: Pt = None;
        let mut base = p;
        let mut k = k;
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(f, acc, base);
            }
            base = self.add(f, base, base);
            k >>= 1;
        }
        acc
    }
}

fn sqrt_lib(field: &Field, v: u64) -> Option<u64> {
    let s = field.sqrt(&field.from_u64(v))?;
    field.to_big(&s).to_u64()
}

fn random_point(f: &Fp, field: &Field, e: &Curve, rng: &mut StdRng) -> (u64, u64) {
    loop {
        let x = rng.gen_range(0..f.p);
        let r = e.rhs(f, x);
        if r == 0 {
            continue;
        }
        if f.pow(r, (f.p - 1) / 2) != 1 {
            continue;
        }
        let y = sqrt_lib(field, r).expect("sqrt");
        return (x, y);
    }
}

/// `#E(F_p)` by baby-step giant-step over the Hasse interval (Mestre style),
/// confirmed on extra random points.  `None` when ambiguous.
fn count_points(f: &Fp, field: &Field, e: &Curve, rng: &mut StdRng) -> Option<u64> {
    let w = (2.0 * (f.p as f64).sqrt()).ceil() as u64 + 1;
    let low = f.p + 1 - w;
    for _ in 0..4 {
        let p = Some(random_point(f, field, e, rng));
        let q = e.smul_big(f, p, low as u128);
        let m = ((2 * w + 1) as f64).sqrt().ceil() as u64 + 1;
        let mut baby: HashMap<(u64, u64), u64> = HashMap::new();
        let mut cur: Pt = None;
        for i in 0..m {
            let key = cur.unwrap_or((u64::MAX, u64::MAX));
            baby.entry(key).or_insert(i);
            cur = e.add(f, cur, p);
        }
        let step = e.neg(f, cur); // -m P
        let mut g = e.neg(f, q);
        let mut sols = Vec::new();
        for k in 0..=m {
            let key = g.unwrap_or((u64::MAX, u64::MAX));
            if let Some(&i) = baby.get(&key) {
                let jj = k * m + i;
                if jj <= 2 * w {
                    sols.push(jj);
                }
            }
            g = e.add(f, g, step);
        }
        sols.sort();
        sols.dedup();
        if sols.len() != 1 {
            continue;
        }
        let n = low + sols[0];
        let mut ok = true;
        for _ in 0..3 {
            let t = Some(random_point(f, field, e, rng));
            if e.smul_big(f, t, n as u128).is_some() {
                ok = false;
            }
        }
        if ok {
            return Some(n);
        }
    }
    None
}

// ------------------------------------------------------------- isogenies ----

/// A Vélu isogeny with a single kernel `x`-coordinate (ℓ = 2: the 2-torsion
/// point; ℓ = 3: the pair ±Q), followed by `(X, Y) ↦ (s X, t Y)`.
#[derive(Clone, Copy, Debug)]
struct Iso {
    ell: u64,
    x0: u64,
    v: u64,
    u: u64,
    cod: Curve,
    s: u64,
    t: u64,
}

#[derive(Default, Clone, Copy)]
struct FieldOps {
    mul: u64,
    inv: u64,
}

fn velu(f: &Fp, e: &Curve, ell: u64, x0: u64) -> Iso {
    let gx = f.add(f.mul(3, f.mul(x0, x0)), e.a);
    let (v, u) = if ell == 2 {
        (gx, 0)
    } else {
        (f.add(gx, gx), f.mul(4, e.rhs(f, x0)))
    };
    let w = f.add(u, f.mul(x0, v));
    let cod = Curve {
        a: f.sub(e.a, f.mul(5, v)),
        b: f.sub(e.b, f.mul(7, w)),
    };
    Iso {
        ell,
        x0,
        v,
        u,
        cod,
        s: 1,
        t: 1,
    }
}

impl Iso {
    /// `X = x + v z + u z²`, `Y = y (1 − v z² − 2u z³)` with `z = 1/(x − x0)`,
    /// then the isomorphism.  Counts field operations.
    fn eval(&self, f: &Fp, p: Pt, c: &mut FieldOps) -> Pt {
        let (x, y) = p?;
        if x == self.x0 {
            return None;
        }
        let z = f.inv(f.sub(x, self.x0));
        let z2 = f.mul(z, z);
        let z3 = f.mul(z2, z);
        let xx = f.add(f.add(x, f.mul(self.v, z)), f.mul(self.u, z2));
        let d = f.sub(
            f.sub(1, f.mul(self.v, z2)),
            f.mul(f.add(self.u, self.u), z3),
        );
        let yy = f.mul(y, d);
        c.inv += 1;
        c.mul += 7;
        if self.s == 1 && self.t == 1 {
            return Some((xx, yy));
        }
        c.mul += 2;
        Some((f.mul(self.s, xx), f.mul(self.t, yy)))
    }
}

/// Rational kernel x-coordinates of ℓ-isogenies (ℓ ∈ {2, 3}) via the library
/// root finder.
fn kernel_xs(field: &Field, f: &Fp, e: &Curve, ell: u64) -> Vec<u64> {
    let c = |v: u64| field.from_u64(v);
    let pol = if ell == 2 {
        vec![c(e.b), c(e.a), c(0), c(1)]
    } else {
        // ψ3 = 3x⁴ + 6a x² + 12 b x − a².
        vec![
            c(f.neg(f.mul(e.a, e.a))),
            c(f.mul(12, e.b)),
            c(f.mul(6, e.a)),
            c(0),
            c(3),
        ]
    };
    poly::roots(field, &pol, 7)
        .iter()
        .map(|r| field.to_big(r).to_u64().unwrap())
        .collect()
}

struct Node {
    e: Curve,
    j: u64,
    parent: usize,
    ell: u64,
    depth: usize,
    deg: u128,
    gen: (u64, u64),
    /// Isogeny from this node's curve back to the parent's curve (the dual).
    dual: Option<Iso>,
}

// -------------------------------------------------------------------- rho ----

fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

const PARTS: usize = 32;

#[derive(Default, Clone)]
struct RhoOut {
    steps: u64,
    setup_ops: u64,
    restarts: u64,
    useless: u64,
    abandoned: u64,
    verify_ops: u64,
    per_target_steps: Vec<u64>,
    logs: Vec<u64>,
    all_verified: bool,
}

/// Distinguished-point rho in `⟨P⟩` (order `r`) for a list of targets, solved
/// in order with one shared DP table (Kuhn–Struik).  With one target this is
/// the ordinary single-target DP rho.
#[allow(clippy::too_many_arguments)]
fn multi_rho(
    f: &Fp,
    e: &Curve,
    r: u64,
    p: Pt,
    targets: &[Pt],
    dp_bits: u32,
    rng: &mut StdRng,
) -> RhoOut {
    let salt: u64 = rng.gen();
    let dp_mask = (1u64 << dp_bits) - 1;
    let max_walk = 20u64 << dp_bits;
    let mut out = RhoOut {
        all_verified: true,
        ..Default::default()
    };
    let mut alpha = [0u64; PARTS];
    let mut rstep: [Pt; PARTS] = [None; PARTS];
    for j in 0..PARTS {
        alpha[j] = rng.gen_range(1..r);
        rstep[j] = e.smul(f, p, alpha[j], &mut out.setup_ops);
    }
    // x → (y, target, a, b): log = a + b·k_target.
    let mut table: HashMap<u64, (u64, u32, u64, u64)> = HashMap::new();
    let mut solved: Vec<Option<u64>> = vec![None; targets.len()];
    for (t, &q) in targets.iter().enumerate() {
        let before = out.steps;
        let a0 = rng.gen_range(0..r);
        let b0 = rng.gen_range(1..r);
        let sa = rng.gen_range(0..r);
        let sb = rng.gen_range(1..r);
        let s1 = e.smul(f, p, a0, &mut out.setup_ops);
        let s2 = e.smul(f, q, b0, &mut out.setup_ops);
        let mut start = e.add(f, s1, s2);
        out.setup_ops += 1;
        let t1 = e.smul(f, p, sa, &mut out.setup_ops);
        let t2 = e.smul(f, q, sb, &mut out.setup_ops);
        let shift = e.add(f, t1, t2);
        out.setup_ops += 1;
        let (mut sa_cur, mut sb_cur) = (a0, b0);
        'walks: loop {
            let (mut x, mut a, b) = (start, sa_cur, sb_cur);
            let mut len = 0u64;
            loop {
                let (px, py) = match x {
                    None => {
                        out.abandoned += 1;
                        break;
                    }
                    Some(v) => v,
                };
                let h = splitmix(px ^ salt);
                if h & dp_mask == 0 {
                    match table.get(&px).copied() {
                        None => {
                            table.insert(px, (py, t as u32, a, b));
                        }
                        Some((ey, et, ea, eb)) => {
                            let sgn_neg = ey != py;
                            let k = if (et as usize) != t {
                                let ke = solved[et as usize].expect("older targets solved");
                                let mut le = (ea as u128 + mul_mod(eb, ke, r) as u128) as u64 % r;
                                if sgn_neg {
                                    le = (r - le) % r;
                                }
                                Some(mul_mod((le + r - a) % r, inv_mod(b, r), r))
                            } else {
                                let (ea2, eb2) = if sgn_neg {
                                    ((r - ea) % r, (r - eb) % r)
                                } else {
                                    (ea, eb)
                                };
                                let den = (b + r - eb2) % r;
                                if den == 0 {
                                    None
                                } else {
                                    let num = (ea2 + r - a) % r;
                                    Some(mul_mod(num, inv_mod(den, r), r))
                                }
                            };
                            match k {
                                Some(k) => {
                                    let chk = e.smul(f, p, k, &mut out.verify_ops);
                                    if chk == q {
                                        solved[t] = Some(k);
                                        break 'walks;
                                    }
                                    out.all_verified = false;
                                    out.useless += 1;
                                }
                                None => out.useless += 1,
                            }
                        }
                    }
                    break;
                }
                if len >= max_walk {
                    out.abandoned += 1;
                    break;
                }
                let jdx = (h >> 59) as usize;
                x = e.add(f, x, rstep[jdx]);
                a = (a + alpha[jdx]) % r;
                out.steps += 1;
                len += 1;
            }
            start = e.add(f, start, shift);
            out.setup_ops += 1;
            out.restarts += 1;
            sa_cur = (sa_cur + sa) % r;
            sb_cur = (sb_cur + sb) % r;
            if sb_cur == 0 {
                sb_cur = 1;
                start = e.add(f, start, q);
                out.setup_ops += 1;
            }
        }
        out.per_target_steps.push(out.steps - before);
    }
    out.logs = solved.into_iter().map(|k| k.unwrap()).collect();
    out
}

// ------------------------------------------------------------------- main ----

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn mean_sd(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (n - 1.0).max(1.0);
    (m, var.sqrt())
}

fn ks_reference(l: usize) -> f64 {
    // Σ_{i<L} binom(2i,i)/4^i · √(π/2): expected rho steps / √r for L targets.
    let mut c = 1.0f64;
    let mut s = 0.0;
    for i in 0..l {
        s += c;
        c *= (2 * i + 1) as f64 / (2 * i + 2) as f64;
    }
    s * (std::f64::consts::PI / 2.0).sqrt()
}

/// Least-squares slope and intercept of `ln y` on `ln x`.
fn loglog_fit(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let lx: Vec<f64> = xs.iter().map(|v| v.ln()).collect();
    let ly: Vec<f64> = ys.iter().map(|v| v.ln()).collect();
    let mx = lx.iter().sum::<f64>() / n;
    let my = ly.iter().sum::<f64>() / n;
    let sxy: f64 = lx.iter().zip(&ly).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = lx.iter().map(|a| (a - mx) * (a - mx)).sum();
    let beta = sxy / sxx;
    (beta, (my - beta * mx).exp())
}

/// `--fit a.json,b.json --out fit.json`: fit recorded totals against L.
/// Per arm: free power law `S = c·L^β` (with a 95% trial-bootstrap interval
/// for β) and the fixed-exponent prefactors `S/L` (arm A) and `S/√L`.
fn fit_mode(files: &str, out_path: &str) {
    let mut res = Vec::new();
    let mut rng = StdRng::seed_from_u64(0x165);
    for path in files.split(',') {
        let v: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        let r = v["base_curve"]["r"].as_f64().unwrap();
        let sr = r.sqrt();
        let ls: Vec<usize> = v["params"]["ls"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_u64().unwrap() as usize)
            .collect();
        // per arm, per L: per-trial S.
        let mut per: Vec<Vec<Vec<f64>>> = vec![vec![Vec::new(); ls.len()]; 3];
        for row in v["trials"].as_array().unwrap() {
            let l = row["L"].as_u64().unwrap() as usize;
            let li = ls.iter().position(|&x| x == l).unwrap();
            let g = |arm: &str, k: &str| row[arm][k].as_f64().unwrap_or(0.0);
            per[0][li].push((g("A", "steps") + g("A", "setup_ops")) / sr);
            per[1][li].push((g("B", "steps") + g("B", "setup_ops") + g("B", "transport_goe")) / sr);
            per[2][li].push((g("C", "steps") + g("C", "setup_ops")) / sr);
        }
        let lx: Vec<f64> = ls.iter().map(|&l| l as f64).collect();
        let mut arms = Vec::new();
        for (ai, name) in ["A", "B", "C"].iter().enumerate() {
            let means: Vec<f64> = per[ai].iter().map(|t| mean_sd(t).0).collect();
            let (beta, c) = loglog_fit(&lx, &means);
            let mut boots = Vec::new();
            for _ in 0..2000 {
                let bm: Vec<f64> = per[ai]
                    .iter()
                    .map(|t| {
                        (0..t.len())
                            .map(|_| t[rng.gen_range(0..t.len())])
                            .sum::<f64>()
                            / t.len() as f64
                    })
                    .collect();
                boots.push(loglog_fit(&lx, &bm).0);
            }
            boots.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let lo = boots[50];
            let hi = boots[1949];
            let pref_lin: Vec<f64> = means.iter().zip(&lx).map(|(m, l)| m / l).collect();
            let pref_sqrt: Vec<f64> = means.iter().zip(&lx).map(|(m, l)| m / l.sqrt()).collect();
            eprintln!(
                "{path} r=2^{:.2} arm {name}: beta={beta:.3} [{lo:.3},{hi:.3}] c={c:.3} S/L={:?} S/sqrtL={:?}",
                r.log2(),
                pref_lin.iter().map(|x| (x * 1000.0).round() / 1000.0).collect::<Vec<_>>(),
                pref_sqrt.iter().map(|x| (x * 1000.0).round() / 1000.0).collect::<Vec<_>>()
            );
            arms.push(json!({"arm": name, "S_means": means, "beta": beta,
                             "beta_ci95_bootstrap": [lo, hi], "c": c,
                             "S_over_L": pref_lin, "S_over_sqrtL": pref_sqrt}));
        }
        res.push(json!({"file": path, "r": r, "log2_r": r.log2(), "ls": ls, "arms": arms}));
    }
    // Derived Kuhn–Struik expected rho steps / √r for the i-th target (i < 64).
    let ks_per_target: Vec<f64> = (0..64)
        .map(|i| ks_reference(i + 1) - ks_reference(i))
        .collect();
    let out = json!({"ks_per_target_S_derived": ks_per_target,
                     "fit": "ln S_mean = ln c + beta ln L, least squares over the recorded L; \
                      bootstrap resamples trials within each L (2000 draws, seed 0x165)",
                     "results": res});
    std::fs::write(out_path, serde_json::to_string_pretty(&out).unwrap()).unwrap();
    eprintln!("wrote {out_path}");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let fit_files: String = arg(&args, "--fit", String::new());
    if !fit_files.is_empty() {
        let out: String = arg(&args, "--out", String::from("fit.json"));
        fit_mode(&fit_files, &out);
        return;
    }
    let rbits: u32 = arg(&args, "--rbits", 32);
    let trials: usize = arg(&args, "--trials", 20);
    let dp_bits: u32 = arg(&args, "--dpbits", 8);
    let seed: u64 = arg(&args, "--seed", 1);
    let out_path: String = arg(&args, "--out", String::from("r4_16_5.json"));
    let ls: Vec<usize> = arg(&args, "--ls", String::from("1,4,16,64"))
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let lmax = *ls.iter().max().unwrap();
    let t0 = Instant::now();
    let mut rng = StdRng::seed_from_u64(seed);

    // ---- 1. Base curve: #E = h·r, r prime, 6 | h, class with > lmax curves.
    let (f, field, e0, n, r, h, nodes, curves_tried) = {
        let mut tried = 0u64;
        loop {
            tried += 1;
            let p = loop {
                let c = rng.gen_range((1u64 << (rbits + 3))..(1u64 << (rbits + 4))) | 1;
                if c % 4 == 3 && is_probable_prime(&BigUint::from(c)) {
                    break c;
                }
            };
            let f = Fp { p };
            let field = Field::new(&BigUint::from(p)).expect("field");
            let e = Curve {
                a: rng.gen_range(1..p),
                b: rng.gen_range(1..p),
            };
            let disc = f.add(
                f.mul(4, f.mul(f.mul(e.a, e.a), e.a)),
                f.mul(27, f.mul(e.b, e.b)),
            );
            if disc == 0 {
                continue;
            }
            let n = match count_points(&f, &field, &e, &mut rng) {
                Some(n) => n,
                None => continue,
            };
            let mut found = None;
            for h in (6..=48u64).step_by(6) {
                if n % h == 0 {
                    let r = n / h;
                    if r >= (1u64 << (rbits - 1)) && is_probable_prime(&BigUint::from(r)) {
                        found = Some((h, r));
                        break;
                    }
                }
            }
            let (h, r) = match found {
                Some(v) => v,
                None => continue,
            };
            // Generator of the order-r subgroup.
            let g0 = loop {
                let q = Some(random_point(&f, &field, &e, &mut rng));
                let g = e.smul_big(&f, q, h as u128);
                if let Some(gv) = g {
                    assert!(e.smul_big(&f, g, r as u128).is_none());
                    break gv;
                }
            };
            // ---- 2. BFS over rational 2- and 3-isogenies.
            let mut nodes = vec![Node {
                e,
                j: e.j(&f),
                parent: 0,
                ell: 1,
                depth: 0,
                deg: 1,
                gen: g0,
                dual: None,
            }];
            let mut seen: HashSet<u64> = HashSet::new();
            seen.insert(nodes[0].j);
            let mut qi = 0;
            let mut bad = false;
            while qi < nodes.len() && nodes.len() <= lmax {
                let cur = nodes[qi].e;
                let cur_gen = nodes[qi].gen;
                for ell in [2u64, 3] {
                    for x0 in kernel_xs(&field, &f, &cur, ell) {
                        if nodes.len() > lmax {
                            break;
                        }
                        let fwd = velu(&f, &cur, ell, x0);
                        let cod = fwd.cod;
                        if cod.a == 0 || cod.b == 0 {
                            continue;
                        }
                        let jc = cod.j(&f);
                        if seen.contains(&jc) || jc == 1728 {
                            continue;
                        }
                        // Cross-check odd-degree codomains against the library.
                        if ell == 3 {
                            let lib = velu_codomain(
                                &field,
                                &Model {
                                    a: field.from_u64(cur.a),
                                    b: field.from_u64(cur.b),
                                },
                                &vec![field.from_u64(f.neg(x0)), field.one()],
                            );
                            let la = field.to_big(&lib.a).to_u64().unwrap();
                            let lb = field.to_big(&lib.b).to_u64().unwrap();
                            assert_eq!((la, lb), (cod.a, cod.b), "Vélu codomain mismatch");
                            let lj = Model { a: lib.a, b: lib.b }.j(&field).unwrap();
                            assert_eq!(field.to_big(&lj).to_u64().unwrap(), jc);
                        }
                        let mut dummy = FieldOps::default();
                        let pg = fwd.eval(&f, Some(cur_gen), &mut dummy);
                        assert!(cod.on(&f, pg));
                        // Homomorphism check on a random full-group point.
                        let tpt = Some(random_point(&f, &field, &cur, &mut rng));
                        let lhs = fwd.eval(&f, cur.add(&f, tpt, Some(cur_gen)), &mut dummy);
                        let rhs = cod.add(
                            &f,
                            fwd.eval(&f, tpt, &mut dummy),
                            fwd.eval(&f, Some(cur_gen), &mut dummy),
                        );
                        if lhs != rhs {
                            bad = true;
                            break;
                        }
                        // Dual: a kernel on the codomain back to cur's j, then
                        // the isomorphism, certified by dual∘fwd = [ℓ].
                        let target = cur.smul_big(&f, Some(cur_gen), ell as u128);
                        let mut dual = None;
                        for y0 in kernel_xs(&field, &f, &cod, ell) {
                            let mut d = velu(&f, &cod, ell, y0);
                            if d.cod.a == 0 || d.cod.b == 0 || d.cod.j(&f) != nodes[qi].j {
                                continue;
                            }
                            // cur: a = s² a'', b = s³ b''  ⇒  s = b a''/(b'' a).
                            let s = f.mul(f.mul(cur.b, d.cod.a), f.inv(f.mul(d.cod.b, cur.a)));
                            if f.mul(f.mul(s, s), d.cod.a) != cur.a
                                || f.mul(f.mul(f.mul(s, s), s), d.cod.b) != cur.b
                            {
                                continue;
                            }
                            let u = match sqrt_lib(&field, s) {
                                Some(u) if f.mul(u, u) == s => u,
                                _ => continue,
                            };
                            d.s = s;
                            d.t = f.mul(s, u);
                            let img = d.eval(&f, pg, &mut dummy);
                            if img == target {
                                dual = Some(d);
                                break;
                            }
                            d.t = f.neg(d.t);
                            if d.eval(&f, pg, &mut dummy) == target {
                                dual = Some(d);
                                break;
                            }
                        }
                        let dual = match dual {
                            Some(d) => d,
                            None => {
                                bad = true;
                                break;
                            }
                        };
                        seen.insert(jc);
                        let (depth, deg) = (nodes[qi].depth + 1, nodes[qi].deg * ell as u128);
                        nodes.push(Node {
                            e: cod,
                            j: jc,
                            parent: qi,
                            ell,
                            depth,
                            deg,
                            gen: pg.unwrap(),
                            dual: Some(dual),
                        });
                    }
                }
                if bad {
                    break;
                }
                qi += 1;
            }
            if bad {
                panic!("isogeny verification failed");
            }
            if nodes.len() > lmax {
                break (f, field, e, n, r, h, nodes, tried);
            }
        }
    };
    let _ = &field;
    let p = f.p;
    let sqrt_r = (r as f64).sqrt();
    eprintln!(
        "p={p} a={} b={} N={n} h={h} r={r} (2^{:.2}); class nodes={} after {curves_tried} curves; {:.1}s",
        e0.a,
        e0.b,
        (r as f64).log2(),
        nodes.len(),
        t0.elapsed().as_secs_f64()
    );
    let slug = curve_id::prime(
        &BigUint::from(p),
        &BigUint::from(e0.a),
        &BigUint::from(e0.b),
        &BigUint::from(n),
    )
    .map(|c| c.slug)
    .unwrap_or_default();

    // Transport along the dual path; counts field ops.
    let transport = |i: usize, q: Pt, c: &mut FieldOps| -> Pt {
        let mut cur = i;
        let mut pt = q;
        while cur != 0 {
            pt = nodes[cur].dual.as_ref().unwrap().eval(&f, pt, c);
            cur = nodes[cur].parent;
        }
        pt
    };
    // Construction check: transport(P_i) = [d_i] P_0 for every node.
    let g0: Pt = Some(nodes[0].gen);
    for i in 1..nodes.len() {
        let mut c = FieldOps::default();
        let img = transport(i, Some(nodes[i].gen), &mut c);
        let want = e0.smul_big(&f, g0, nodes[i].deg % r as u128);
        assert_eq!(img, want, "transport of generator {i}");
    }
    let curves_json: Vec<Value> = nodes
        .iter()
        .enumerate()
        .map(|(i, nd)| {
            json!({"i": i, "a": nd.e.a, "b": nd.e.b, "j": nd.j, "parent": nd.parent,
                   "ell": nd.ell, "depth": nd.depth, "deg": nd.deg.to_string(),
                   "gen": [nd.gen.0, nd.gen.1],
                   "dual_ell": nd.dual.map(|d| d.ell),
                   "dual_kernel_x": nd.dual.map(|d| d.x0)})
        })
        .collect();

    // ---- 3. Trials.
    let mut rows = Vec::new();
    let mut summary = Vec::new();
    for &l in &ls {
        let mut tot = [Vec::new(), Vec::new(), Vec::new()];
        let mut steps_only = [Vec::new(), Vec::new(), Vec::new()];
        let mut tr_mul = Vec::new();
        let mut tr_inv = Vec::new();
        let mut per_target_b = vec![0f64; l];
        let mut per_target_c = vec![0f64; l];
        let mut all_ok = true;
        for trial in 0..trials {
            let tseed = splitmix(seed ^ ((l as u64) << 32) ^ trial as u64);
            let mut trng = StdRng::seed_from_u64(tseed);
            // Planted keys on E_1..E_L.
            let keys: Vec<u64> = (0..l).map(|_| trng.gen_range(1..r)).collect();
            let mut dummy = 0u64;
            let qs: Vec<Pt> = (0..l)
                .map(|i| {
                    let nd = &nodes[i + 1];
                    nd.e.smul(&f, Some(nd.gen), keys[i], &mut dummy)
                })
                .collect();
            // Arm A: independent runs on each curve.
            let mut a_steps = 0u64;
            let mut a_setup = 0u64;
            let mut a_ok = true;
            for i in 0..l {
                let nd = &nodes[i + 1];
                let o = multi_rho(&f, &nd.e, r, Some(nd.gen), &[qs[i]], dp_bits, &mut trng);
                a_ok &= o.all_verified && o.logs[0] == keys[i];
                a_steps += o.steps;
                a_setup += o.setup_ops;
            }
            // Arm B: transport to E_0, one multi-target run, map back, verify.
            let mut c = FieldOps::default();
            let tq: Vec<Pt> = (0..l).map(|i| transport(i + 1, qs[i], &mut c)).collect();
            let ob = multi_rho(&f, &e0, r, g0, &tq, dp_bits, &mut trng);
            let mut b_ok = ob.all_verified;
            let mut back_ops = 0u64;
            for i in 0..l {
                let nd = &nodes[i + 1];
                let dinv = inv_mod((nd.deg % r as u128) as u64, r);
                let k = mul_mod(ob.logs[i], dinv, r);
                let chk = nd.e.smul(&f, Some(nd.gen), k, &mut back_ops);
                b_ok &= chk == qs[i] && k == keys[i];
            }
            // Arm C: L targets planted directly on E_0.
            let keys_c: Vec<u64> = (0..l).map(|_| trng.gen_range(1..r)).collect();
            let qc: Vec<Pt> = keys_c
                .iter()
                .map(|&k| e0.smul(&f, g0, k, &mut dummy))
                .collect();
            let oc = multi_rho(&f, &e0, r, g0, &qc, dp_bits, &mut trng);
            let c_ok = oc.all_verified && oc.logs == keys_c;
            all_ok &= a_ok && b_ok && c_ok;
            for (k, v) in ob.per_target_steps.iter().enumerate() {
                per_target_b[k] += *v as f64 / trials as f64;
            }
            for (k, v) in oc.per_target_steps.iter().enumerate() {
                per_target_c[k] += *v as f64 / trials as f64;
            }
            // Transport in group-op equivalents: one affine addition costs one
            // inversion plus three multiplications in this code.
            let tr_goe = c.inv as f64 + (c.mul as f64 - 3.0 * c.inv as f64).max(0.0) / 3.0;
            let totals = [
                (a_steps + a_setup) as f64,
                (ob.steps + ob.setup_ops) as f64 + tr_goe,
                (oc.steps + oc.setup_ops) as f64,
            ];
            for k in 0..3 {
                tot[k].push(totals[k]);
            }
            steps_only[0].push(a_steps as f64);
            steps_only[1].push(ob.steps as f64);
            steps_only[2].push(oc.steps as f64);
            tr_mul.push(c.mul as f64);
            tr_inv.push(c.inv as f64);
            rows.push(json!({
                "L": l, "trial": trial, "seed": tseed,
                "A": {"steps": a_steps, "setup_ops": a_setup, "verified": a_ok},
                "B": {"steps": ob.steps, "setup_ops": ob.setup_ops, "restarts": ob.restarts,
                      "useless": ob.useless, "abandoned": ob.abandoned,
                      "transport_mul": c.mul, "transport_inv": c.inv,
                      "transport_goe": tr_goe, "back_verify_ops": back_ops,
                      "per_target_steps": ob.per_target_steps, "verified": b_ok},
                "C": {"steps": oc.steps, "setup_ops": oc.setup_ops, "restarts": oc.restarts,
                      "useless": oc.useless, "abandoned": oc.abandoned,
                      "per_target_steps": oc.per_target_steps, "verified": c_ok},
            }));
        }
        let ms: Vec<(f64, f64)> = tot.iter().map(|v| mean_sd(v)).collect();
        let st: Vec<(f64, f64)> = steps_only.iter().map(|v| mean_sd(v)).collect();
        let nt = trials as f64;
        // Ratio B/C of means with a delta-method 95% interval (independent arms).
        let (mb, sb) = ms[1];
        let (mc, sc) = ms[2];
        let ratio = mb / mc;
        let se = ratio * ((sb / mb).powi(2) / nt + (sc / mc).powi(2) / nt).sqrt();
        let (mtm, _) = mean_sd(&tr_mul);
        let (mti, _) = mean_sd(&tr_inv);
        eprintln!(
            "L={l:3}: A {:.0} (S={:.3})  B {:.0} (S={:.3})  C {:.0} (S={:.3})  B/C={ratio:.3}±{:.3}  A/B={:.2}  KS={:.3} ok={all_ok} [{:.0}s]",
            ms[0].0,
            ms[0].0 / sqrt_r,
            mb,
            mb / sqrt_r,
            mc,
            mc / sqrt_r,
            1.96 * se,
            ms[0].0 / mb,
            ks_reference(l),
            t0.elapsed().as_secs_f64()
        );
        summary.push(json!({
            "L": l, "trials": trials, "all_verified": all_ok,
            "A_total_mean": ms[0].0, "A_total_sd": ms[0].1,
            "B_total_mean": mb, "B_total_sd": sb,
            "C_total_mean": mc, "C_total_sd": sc,
            "A_steps_mean": st[0].0, "B_steps_mean": st[1].0, "C_steps_mean": st[2].0,
            "A_S": ms[0].0 / sqrt_r, "B_S": mb / sqrt_r, "C_S": mc / sqrt_r,
            "A_S_sd": ms[0].1 / sqrt_r, "B_S_sd": sb / sqrt_r, "C_S_sd": sc / sqrt_r,
            "B_over_C": ratio, "B_over_C_ci95": [ratio - 1.96 * se, ratio + 1.96 * se],
            "A_over_B": ms[0].0 / mb,
            "per_target_A_S": ms[0].0 / sqrt_r / l as f64,
            "per_target_B_S": mb / sqrt_r / l as f64,
            "ks_reference_S_derived": ks_reference(l),
            "independent_reference_S_derived": l as f64 * ks_reference(1),
            "transport_mul_mean": mtm, "transport_inv_mean": mti,
            "per_target_steps_B_mean": per_target_b, "per_target_steps_C_mean": per_target_c,
        }));
    }
    let depths: Vec<usize> = nodes.iter().map(|nd| nd.depth).collect();
    let report = json!({
        "experiment": "r4_16_5_class_multitarget_rho",
        "unit": "group operations in the order-r subgroup (affine add/double); transport in field mul/inv",
        "params": {"rbits": rbits, "trials": trials, "dp_bits": dp_bits, "seed": seed,
                   "ls": ls, "partitions": PARTS},
        "base_curve": {"p": p, "a": e0.a, "b": e0.b, "order": n, "cofactor": h, "r": r,
                       "log2_r": (r as f64).log2(), "icv1_slug": slug,
                       "curves_tried": curves_tried},
        "class": {"nodes": nodes.len(), "max_depth": depths.iter().max(), "curves": curves_json},
        "summary": summary,
        "trials": rows,
        "contended_wall_s": t0.elapsed().as_secs_f64(),
    });
    std::fs::write(&out_path, serde_json::to_string_pretty(&report).unwrap()).unwrap();
    eprintln!("wrote {out_path}");
}
