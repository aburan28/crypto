//! # Curves and specs for research/fiber_ic_scaling_20261011
//!
//! Finds random prime-field curves `y² = x³ + ax + b` whose group order is
//! `h·r` with `r` prime and a **fixed** cofactor `h`, at several field
//! sizes.  It writes one `ecbench` spec per cofactor holding every size,
//! with rho and the fiber-closed and fiber-blind MITM index calculus over a
//! sweep of factor-base sizes.
//!
//! `#E` comes from baby-step giant-step on the Hasse interval: the order
//! of a random point.  It is accepted only when exactly one multiple lies
//! in the interval, and confirmed on a second point.  `ecbench`'s
//! `prime_explicit` constructor then re-checks `[r]G = O` on load.
//!
//! ```bash
//! cargo run --release --example fiber_scaling_curves -- \
//!   --cofactors 4,72 --bits 18,22,26,30 --per-size 2 --seed 20261011 \
//!   --out research/fiber_ic_scaling_20261011
//! ```

use std::collections::HashMap;

use serde_json::{json, Value};

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
    pow(a, p - 2, p)
}
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(q) {
            return n == q;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = pow(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mul(x, x, n);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}
fn sqrt(a: u64, p: u64) -> Option<u64> {
    let a = a % p;
    if a == 0 {
        return Some(0);
    }
    if pow(a, (p - 1) / 2, p) != 1 {
        return None;
    }
    let (mut q, mut s) = (p - 1, 0u32);
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while pow(z, (p - 1) / 2, p) == 1 {
        z += 1;
    }
    let (mut m, mut c, mut t, mut r) = (s, pow(z, q, p), pow(a, q, p), pow(a, q.div_ceil(2), p));
    while t != 1 {
        let (mut i, mut tt) = (0, t);
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

type Pt = Option<(u64, u64)>;

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
    fn add(&self, s: Pt, t: Pt) -> Pt {
        let p = self.p;
        let (Some((x1, y1)), Some((x2, y2))) = (s, t) else {
            return s.or(t);
        };
        let lam = if x1 == x2 {
            if add(y1, y2, p) == 0 {
                return None;
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
        Some((x3, sub(mul(lam, sub(x1, x3, p), p), y1, p)))
    }
    fn neg(&self, s: Pt) -> Pt {
        s.map(|(x, y)| (x, (self.p - y) % self.p))
    }
    fn mul(&self, mut k: u64, q: Pt) -> Pt {
        let (mut acc, mut base) = (None, q);
        while k > 0 {
            if k & 1 == 1 {
                acc = self.add(acc, base);
            }
            base = self.add(base, base);
            k >>= 1;
        }
        acc
    }
    /// Every `m` in the Hasse interval with `[m]P = O` (baby-step giant-step).
    fn hasse_multiples(&self, pt: Pt) -> Vec<u64> {
        let p = self.p;
        let w = 2 * ((p as f64).sqrt() as u64 + 1);
        let lo = (p + 1).saturating_sub(w);
        let hi = p + 1 + w;
        let m = ((hi - lo) as f64).sqrt() as u64 + 1;
        // Baby steps: j·P for 0 ≤ j < m.
        let mut baby: HashMap<Pt, Vec<u64>> = HashMap::new();
        let mut cur: Pt = None;
        for j in 0..m {
            baby.entry(cur).or_default().push(j);
            cur = self.add(cur, pt);
        }
        // Giant steps: is −(lo + i·m)P a baby step?
        let step = self.neg(self.mul(m, pt));
        let mut g = self.neg(self.mul(lo, pt));
        let mut out = Vec::new();
        let mut i = 0;
        while lo + i * m <= hi {
            if let Some(js) = baby.get(&g) {
                for &j in js {
                    let n = lo + i * m + j;
                    if (lo..=hi).contains(&n) {
                        out.push(n);
                    }
                }
            }
            g = self.add(g, step);
            i += 1;
        }
        out.sort();
        out.dedup();
        out
    }
}

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
}

fn find(bits: u32, h: u64, rng: &mut Lcg) -> Value {
    loop {
        let p = (rng.next() & ((1 << bits) - 1)) | (1 << (bits - 1)) | 1;
        if !is_prime(p) {
            continue;
        }
        for _ in 0..200 {
            let e = Curve {
                p,
                a: rng.next() % p,
                b: 1 + rng.next() % (p - 1),
            };
            let disc = add(
                mul(4, mul(mul(e.a, e.a, p), e.a, p), p),
                mul(27, mul(e.b, e.b, p), p),
                p,
            );
            if disc == 0 {
                continue;
            }
            let point = |rng: &mut Lcg| loop {
                let x = rng.next() % p;
                if let Some(y) = sqrt(e.rhs(x), p) {
                    return Some((x, y));
                }
            };
            let p1 = point(rng);
            let ms = e.hasse_multiples(p1);
            if ms.len() != 1 {
                continue;
            }
            let n = ms[0];
            if !n.is_multiple_of(h) || !is_prime(n / h) {
                continue;
            }
            let r = n / h;
            // Confirm #E on a second point, and take a generator of E[r].
            if e.mul(n, point(rng)).is_some() {
                continue;
            }
            let g = loop {
                let q = e.mul(h, point(rng));
                if q.is_some() {
                    break q;
                }
            };
            assert!(e.mul(r, g).is_none());
            let (gx, gy) = g.unwrap();
            return json!({"kind": "prime_explicit", "p": p, "a": e.a, "b": e.b,
                          "group_order": n, "r": r, "gx": gx, "gy": gy});
        }
    }
}

fn arg(args: &[String], name: &str) -> String {
    let i = args
        .iter()
        .position(|a| a == name)
        .unwrap_or_else(|| panic!("missing {name}"));
    args[i + 1].clone()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let list = |s: String| -> Vec<u64> { s.split(',').map(|v| v.parse().unwrap()).collect() };
    let cofactors = list(arg(&args, "--cofactors"));
    let bits = list(arg(&args, "--bits"));
    let per: usize = arg(&args, "--per-size").parse().unwrap();
    let seed: u64 = arg(&args, "--seed").parse().unwrap();
    let out = arg(&args, "--out");
    std::fs::create_dir_all(format!("{out}/specs")).unwrap();
    let rho = json!({"name": "rho-neg", "role": "reference", "method": {"id": "rho.negation"}});
    let ic = |name: String, size: u64, fiber: &str, role: &str| {
        let oracle = if fiber == "closed" {
            "mitm-fiber:fiber=closed,negation_folded=1,fiber_gen=coset".to_string()
        } else {
            "mitm-fiber:fiber=blind,negation_folded=1".to_string()
        };
        json!({"name": name, "role": role, "method": {"id": "ic.pipeline", "params": {
            "factor_base": format!("prime-abscissa:size={size}"), "oracle": oracle}}})
    };
    for &h in &cofactors {
        let mut rng = Lcg(seed ^ h.rotate_left(29));
        let mut curves = Vec::new();
        for &b in &bits {
            for _ in 0..per {
                let c = find(b as u32, h, &mut rng);
                eprintln!(
                    "h={h} bits={b}: p={} #E={} r={}",
                    c["p"], c["group_order"], c["r"]
                );
                curves.push(c);
            }
        }
        let mut arms = vec![rho.clone()];
        for size in [16u64, 32, 64, 128, 256] {
            arms.push(ic(format!("closed-b{size}"), size, "closed", "candidate"));
        }
        arms.push(ic("blind-b128".into(), 128, "blind", "control"));
        let spec = json!({
            "schema": "ecbench.spec/v1",
            "label": format!("fiber-closed IC against rho across field sizes, cofactor {h}"),
            "description": format!("research/fiber_ic_scaling_20261011/PROTOCOL.md: random prime curves with #E = {h}·r at {bits:?} bits ({per} per size); rho against the coset fiber-closed MITM at factor-base sizes 16 to 256 and the blind MITM at 128."),
            "arms": arms,
            "measurement": {"isolation_required": "L2", "order": "alternate", "rounds": 2, "seed": 1,
                            "timeout_seconds": 900, "warmup": 0},
            "workloads": {"curves": curves, "target_seed": 20261011, "targets_per_curve": 3},
        });
        std::fs::write(
            format!("{out}/specs/scaling_h{h}.json"),
            serde_json::to_string_pretty(&spec).unwrap() + "\n",
        )
        .unwrap();
    }
}
