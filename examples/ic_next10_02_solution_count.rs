//! IC next-10, item 2: combined-system solution-count census for algebraic factor bases.
//!
//! For a prime-order curve E: y^2 = x^3 + a x + b over F_p, a factor base V (|V| = B = 2^t)
//! and a target R, count the distinct ordered x-tuples (x_1..x_m) in V^m with
//! S_{m+1}(x_1, .., x_m, x_R) = 0, i.e. the F_p-points (= F̄_p points, since V ⊂ F_p) of the
//! combined system {S_{m+1} = 0, L(x_i) = 0}. Enumeration is through the group law:
//! m = 2: x(R - P) ∈ V for every signed base point P;
//! m = 3: a sorted table of x(P_a + P_b) over ordered signed pairs, matched to x(R - P_c).
//!
//! Arms: kummer (g·μ_B), subgroup (μ_B), dickson (split), dickson-torus (needs B | p+1),
//! smallx ({1..B}), random (null), onlycurve and dloginterval (positive controls).
//!
//! Output: one CSV row per (arm, m, t) cell on stdout. With `--dump-bases DIR`, also writes each
//! base and the targets for the msolve validation arm.
//!
//! Usage: ic_next10_02_solution_count --p P --a A --b B --n N --tmin 2 --tmax 10
//!        --targets2 2000 --targets3 500 --seed S [--dump-bases DIR --dump-targets K]

use std::collections::{HashMap, HashSet};
use std::io::Write;

#[derive(Clone, Copy)]
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        // splitmix64
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, m: u64) -> u64 {
        self.next() % m
    }
}

#[derive(Clone, Copy)]
struct Fp {
    p: u64,
}
impl Fp {
    fn add(&self, x: u64, y: u64) -> u64 {
        let s = x + y;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    fn sub(&self, x: u64, y: u64) -> u64 {
        if x >= y {
            x - y
        } else {
            x + self.p - y
        }
    }
    fn mul(&self, x: u64, y: u64) -> u64 {
        ((x as u128 * y as u128) % self.p as u128) as u64
    }
    fn pow(&self, mut x: u64, mut e: u64) -> u64 {
        let mut r = 1;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, x);
            }
            x = self.mul(x, x);
            e >>= 1;
        }
        r
    }
    fn inv(&self, x: u64) -> u64 {
        assert!(x != 0);
        self.pow(x, self.p - 2)
    }
    fn is_sq(&self, x: u64) -> bool {
        x == 0 || self.pow(x, (self.p - 1) / 2) == 1
    }
    fn sqrt(&self, x: u64) -> u64 {
        // Tonelli-Shanks
        if x == 0 {
            return 0;
        }
        let p = self.p;
        let mut q = p - 1;
        let mut s = 0;
        while q.is_multiple_of(2) {
            q /= 2;
            s += 1;
        }
        let mut z = 2;
        while self.is_sq(z) {
            z += 1;
        }
        let mut m = s;
        let mut c = self.pow(z, q);
        let mut t = self.pow(x, q);
        let mut r = self.pow(x, q.div_ceil(2));
        while t != 1 {
            let mut i = 0;
            let mut tt = t;
            while tt != 1 {
                tt = self.mul(tt, tt);
                i += 1;
            }
            let b = self.pow(c, 1 << (m - i - 1));
            m = i;
            c = self.mul(b, b);
            t = self.mul(t, c);
            r = self.mul(r, b);
        }
        r
    }
}

// F_{p^2} = F_p[i]/(i^2 - nr)
#[derive(Clone, Copy, PartialEq)]
struct F2(u64, u64);
fn f2_mul(f: &Fp, nr: u64, x: F2, y: F2) -> F2 {
    F2(
        f.add(f.mul(x.0, y.0), f.mul(nr, f.mul(x.1, y.1))),
        f.add(f.mul(x.0, y.1), f.mul(x.1, y.0)),
    )
}
fn f2_pow(f: &Fp, nr: u64, mut x: F2, mut e: u64) -> F2 {
    let mut r = F2(1, 0);
    while e > 0 {
        if e & 1 == 1 {
            r = f2_mul(f, nr, r, x);
        }
        x = f2_mul(f, nr, x, x);
        e >>= 1;
    }
    r
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Pt {
    O,
    A(u64, u64),
}

struct Curve {
    f: Fp,
    a: u64,
    b: u64,
}
impl Curve {
    fn rhs(&self, x: u64) -> u64 {
        let f = &self.f;
        f.add(f.add(f.mul(f.mul(x, x), x), f.mul(self.a, x)), self.b)
    }
    fn neg(&self, p: Pt) -> Pt {
        match p {
            Pt::O => Pt::O,
            Pt::A(x, y) => Pt::A(x, if y == 0 { 0 } else { self.f.p - y }),
        }
    }
    fn add(&self, p: Pt, q: Pt) -> Pt {
        let f = &self.f;
        match (p, q) {
            (Pt::O, _) => q,
            (_, Pt::O) => p,
            (Pt::A(x1, y1), Pt::A(x2, y2)) => {
                let l = if x1 == x2 {
                    if f.add(y1, y2) == 0 {
                        return Pt::O;
                    }
                    let num = f.add(f.mul(3, f.mul(x1, x1)), self.a);
                    f.mul(num, f.inv(f.mul(2, y1)))
                } else {
                    f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1)))
                };
                let x3 = f.sub(f.sub(f.mul(l, l), x1), x2);
                let y3 = f.sub(f.mul(l, f.sub(x1, x3)), y1);
                Pt::A(x3, y3)
            }
        }
    }
    fn mul(&self, mut k: u64, mut p: Pt) -> Pt {
        let mut r = Pt::O;
        while k > 0 {
            if k & 1 == 1 {
                r = self.add(r, p);
            }
            p = self.add(p, p);
            k >>= 1;
        }
        r
    }
    fn random_point(&self, rng: &mut Rng) -> Pt {
        loop {
            let x = rng.below(self.f.p);
            let r = self.rhs(x);
            if r != 0 && self.f.is_sq(r) {
                let y = self.f.sqrt(r);
                return if rng.next() & 1 == 0 {
                    Pt::A(x, y)
                } else {
                    Pt::A(x, self.f.p - y)
                };
            }
        }
    }
}

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter()
        .position(|s| s == k)
        .map(|i| args[i + 1].clone())
}

/// Element of exact order B (power of two) in F_p^*, B | p-1.
fn root_of_unity(f: &Fp, bsz: u64) -> u64 {
    let mut h = 2;
    while f.is_sq(h) {
        h += 1;
    }
    f.pow(h, (f.p - 1) / bsz)
}

fn build_base(arm: &str, c: &Curve, n: u64, bsz: u64, rng: &mut Rng) -> Option<Vec<u64>> {
    let f = &c.f;
    let p = f.p;
    let v: Vec<u64> = match arm {
        "kummer" | "subgroup" => {
            if !(p - 1).is_multiple_of(bsz) {
                return None;
            }
            let z = root_of_unity(f, bsz);
            let g = if arm == "subgroup" {
                1
            } else {
                loop {
                    let g = 2 + rng.below(p - 2);
                    if f.pow(g, bsz) != 1 {
                        break g;
                    }
                }
            };
            let mut w = g;
            (0..bsz)
                .map(|_| {
                    let x = w;
                    w = f.mul(w, z);
                    x
                })
                .collect()
        }
        "dickson" => {
            if !(p - 1).is_multiple_of(bsz) {
                return None;
            }
            let z = root_of_unity(f, bsz);
            let u = loop {
                let u = 2 + rng.below(p - 2);
                if f.pow(u, 2 * bsz) != 1 {
                    break u;
                }
            };
            let mut w = u;
            (0..bsz)
                .map(|_| {
                    let x = f.add(w, f.inv(w));
                    w = f.mul(w, z);
                    x
                })
                .collect()
        }
        "dickson-torus" => {
            if !(p + 1).is_multiple_of(bsz) {
                return None;
            }
            let mut nr = 2;
            while f.is_sq(nr) {
                nr += 1;
            }
            // torus element: z^(p-1) for random z in F_{p^2}^*
            let torus = |rng: &mut Rng| loop {
                let z = F2(rng.below(p), 1 + rng.below(p - 1));
                let w = f2_pow(f, nr, z, p - 1);
                if w != F2(1, 0) {
                    break w;
                }
            };
            let zeta = loop {
                let w = torus(rng);
                let ze = f2_pow(f, nr, w, (p + 1) / bsz);
                if bsz == 1 || f2_pow(f, nr, ze, bsz / 2) != F2(1, 0) {
                    break ze;
                }
            };
            let u = loop {
                let u = torus(rng);
                if f2_pow(f, nr, u, 2 * bsz) != F2(1, 0) {
                    break u;
                }
            };
            let mut w = u;
            (0..bsz)
                .map(|_| {
                    // w + w^{-1} = w + conj(w) = 2 Re(w) on the norm-1 torus
                    let x = f.add(w.0, w.0);
                    w = f2_mul(f, nr, w, zeta);
                    x
                })
                .collect()
        }
        "smallx" => (1..=bsz).collect(),
        "random" => {
            let mut s = HashSet::new();
            let mut v = Vec::new();
            while v.len() < bsz as usize {
                let x = rng.below(p);
                if s.insert(x) {
                    v.push(x);
                }
            }
            v
        }
        "onlycurve" => {
            let mut s = HashSet::new();
            let mut v = Vec::new();
            while v.len() < bsz as usize {
                let x = rng.below(p);
                let r = c.rhs(x);
                if r != 0 && f.is_sq(r) && s.insert(x) {
                    v.push(x);
                }
            }
            v
        }
        "dloginterval" => {
            let g = c.random_point(rng);
            assert!(c.mul(n, g) == Pt::O);
            let mut q = g;
            (0..bsz)
                .map(|_| {
                    let x = match q {
                        Pt::A(x, _) => x,
                        Pt::O => unreachable!(),
                    };
                    q = c.add(q, g);
                    x
                })
                .collect()
        }
        _ => panic!("unknown arm {arm}"),
    };
    let s: HashSet<u64> = v.iter().copied().collect();
    assert_eq!(s.len(), v.len(), "base {arm} B={bsz} not distinct");
    Some(v)
}

fn signed_points(c: &Curve, v: &[u64]) -> Vec<Pt> {
    let mut out = Vec::new();
    for &x in v {
        let r = c.rhs(x);
        if r != 0 && c.f.is_sq(r) {
            let y = c.f.sqrt(r);
            out.push(Pt::A(x, y));
            out.push(Pt::A(x, c.f.p - y));
        }
    }
    out
}

fn xof(p: Pt) -> u64 {
    match p {
        Pt::A(x, _) => x,
        Pt::O => u64::MAX,
    }
}

fn count_m2(c: &Curve, vs: &HashSet<u64>, sp: &[Pt], r: Pt) -> usize {
    let mut set = HashSet::new();
    for &pt in sp {
        let q = c.add(r, c.neg(pt));
        if q != Pt::O && vs.contains(&xof(q)) {
            set.insert((xof(pt), xof(q)));
        }
    }
    set.len()
}

fn pair_table(c: &Curve, sp: &[Pt]) -> Vec<(u64, u32, u32)> {
    let mut t = Vec::with_capacity(sp.len() * sp.len());
    for (i, &pa) in sp.iter().enumerate() {
        for (j, &pb) in sp.iter().enumerate() {
            let s = c.add(pa, pb);
            if s != Pt::O {
                t.push((xof(s), i as u32, j as u32));
            }
        }
    }
    t.sort_unstable();
    t
}

fn count_m3(c: &Curve, sp: &[Pt], tab: &[(u64, u32, u32)], r: Pt) -> usize {
    let mut set = HashSet::new();
    for &pc in sp {
        let q = c.add(r, c.neg(pc));
        if q == Pt::O {
            continue;
        }
        let x = xof(q);
        let lo = tab.partition_point(|e| e.0 < x);
        let mut k = lo;
        while k < tab.len() && tab[k].0 == x {
            let (_, a, b) = tab[k];
            set.insert((xof(sp[a as usize]), xof(sp[b as usize]), xof(pc)));
            k += 1;
        }
    }
    set.len()
}

fn stats(counts: &[usize], rng: &mut Rng) -> (f64, f64, f64, f64, f64) {
    let k = counts.len() as f64;
    let mean = counts.iter().sum::<usize>() as f64 / k;
    let var = counts
        .iter()
        .map(|&c| (c as f64 - mean).powi(2))
        .sum::<f64>()
        / (k - 1.0);
    let frac = counts.iter().filter(|&&c| c >= 1).count() as f64 / k;
    let mut boots: Vec<f64> = (0..2000)
        .map(|_| {
            let mut s = 0usize;
            for _ in 0..counts.len() {
                s += counts[rng.below(counts.len() as u64) as usize];
            }
            s as f64 / k
        })
        .collect();
    boots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (mean, boots[49], boots[1949], var, frac)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let g = |k: &str| {
        arg(&args, k)
            .unwrap_or_else(|| panic!("missing {k}"))
            .parse::<u64>()
            .unwrap()
    };
    let (p, a, b, n) = (g("--p"), g("--a"), g("--b"), g("--n"));
    let tmin = g("--tmin");
    let tmax = g("--tmax");
    let targets2 = g("--targets2") as usize;
    let targets3 = g("--targets3") as usize;
    let seed = g("--seed");
    let dump = arg(&args, "--dump-bases");
    let dump_targets = arg(&args, "--dump-targets")
        .map(|s| s.parse::<usize>().unwrap())
        .unwrap_or(5);
    let ms: Vec<usize> = arg(&args, "--ms")
        .unwrap_or("2,3".into())
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let c = Curve { f: Fp { p }, a, b };

    // shared targets
    let mut trng = Rng(seed ^ 0x7a7a);
    let tmax_targets = targets2.max(targets3);
    let targets: Vec<Pt> = (0..tmax_targets)
        .map(|_| c.random_point(&mut trng))
        .collect();
    // every target must lie in the prime-order group (cofactor 1, so automatic) -- spot check
    assert!(c.mul(n, targets[0]) == Pt::O);

    if let Some(dir) = &dump {
        std::fs::create_dir_all(dir).unwrap();
        let mut fh = std::fs::File::create(format!("{dir}/targets.csv")).unwrap();
        writeln!(fh, "idx,xR,yR").unwrap();
        for (i, t) in targets.iter().take(dump_targets).enumerate() {
            if let Pt::A(x, y) = t {
                writeln!(fh, "{i},{x},{y}").unwrap();
            }
        }
    }

    let arms = [
        "kummer",
        "subgroup",
        "dickson",
        "dickson-torus",
        "smallx",
        "random",
        "onlycurve",
        "dloginterval",
    ];
    println!("p,a,b,n,arm,m,t,B,B_E,targets,skipped,mean,ci_lo,ci_hi,var,frac_ge1,pred_exact,pred_product,Bm1,hist");
    for t in tmin..=tmax {
        let bsz = 1u64 << t;
        for (ai, arm) in arms.iter().enumerate() {
            let mut brng = Rng(seed.wrapping_mul(1000003) ^ (t << 8) ^ ai as u64);
            let v = match build_base(arm, &c, n, bsz, &mut brng) {
                Some(v) => v,
                None => continue,
            };
            if let Some(dir) = &dump {
                let mut fh = std::fs::File::create(format!("{dir}/base_{arm}_t{t}.txt")).unwrap();
                for x in &v {
                    writeln!(fh, "{x}").unwrap();
                }
            }
            let vs: HashSet<u64> = v.iter().copied().collect();
            let sp = signed_points(&c, &v);
            let be = sp.len() / 2;
            for &m in &ms {
                let ntar = if m == 2 { targets2 } else { targets3 };
                let tab = if m == 3 {
                    pair_table(&c, &sp)
                } else {
                    Vec::new()
                };
                let mut counts = Vec::with_capacity(ntar);
                let mut skipped = 0;
                for &r in targets.iter().take(ntar) {
                    if vs.contains(&xof(r)) {
                        skipped += 1;
                        continue;
                    }
                    counts.push(if m == 2 {
                        count_m2(&c, &vs, &sp, r)
                    } else {
                        count_m3(&c, &sp, &tab, r)
                    });
                }
                let mut srng = Rng(seed ^ 0xb00757 ^ (t << 16) ^ ((ai as u64) << 4) ^ m as u64);
                let (mean, lo, hi, var, frac) = stats(&counts, &mut srng);
                let pred_exact = (2.0 * be as f64).powi(m as i32) / n as f64;
                let pred_product = (bsz as f64).powi(m as i32) / p as f64;
                let bm1 = (bsz as f64).powi(m as i32 - 1);
                let mut h: HashMap<usize, usize> = HashMap::new();
                for &k in &counts {
                    *h.entry(k).or_default() += 1;
                }
                let mut hv: Vec<_> = h.into_iter().collect();
                hv.sort();
                let hist: Vec<String> = hv.iter().map(|(k, v)| format!("{k}:{v}")).collect();
                println!(
                    "{p},{a},{b},{n},{arm},{m},{t},{bsz},{be},{},{skipped},{mean:.6},{lo:.6},{hi:.6},{var:.6},{frac:.6},{pred_exact:.6},{pred_product:.6},{bm1},{}",
                    counts.len(),
                    hist.join(" ")
                );
                std::io::stdout().flush().unwrap();
            }
        }
    }
}
