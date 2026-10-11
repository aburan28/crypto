//! Multiplier closure in the 4-sum pair-table index calculus on toy Koblitz
//! curves `y² + xy = x³ + 1` over `F_{2^n}`.
//!
//! Three multiplier families generate the summand set `S` from `K` hidden
//! columns `P_i`:
//!
//!   frob   M = {±τ^k}                       (the landed compact-orbit design)
//!   tau    M = {±α τ^k : α ∈ A_c}           A_c = the c smallest-norm
//!                                            τ-indivisible elements of Z[τ]
//!   int    M = {±a τ^k : a ∈ {1, …, c}}     plain integer multiples
//!
//! Every arm of a cell has the same `K·c`, hence the same `|S| = 2n·K·c` and
//! the same pair table.  A relation is a 4-sum `ρG = s_a + s_b + s_c + s_d`
//! found by probing random pairs `(c, d)` against the table of `x(s_a + s_b)`.
//! Columns are solved mod `r` and checked against the hidden scalars; then one
//! hidden target is decomposed and its scalar verified.
//!
//! Standalone: `rustc -O multipliers.rs && ./multipliers <n> <seed> <Kc> <hits> <out.json>`
//! Pure std; no crate.  See PROTOCOL.md for the frozen predictions.

use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

// ---------------------------------------------------------------- field

#[derive(Clone, Copy)]
struct Field {
    n: u32,
    poly: u64,
}

impl Field {
    fn new(n: u32) -> Field {
        // Find the lowest-weight irreducible x^n + x^k + 1 (or pentanomial).
        for k in 1..n {
            let f = Field { n, poly: (1u64 << n) | (1u64 << k) | 1 };
            if f.is_irreducible() {
                return f;
            }
        }
        for k1 in 1..n {
            for k2 in 1..k1 {
                for k3 in 1..k2 {
                    let f = Field { n, poly: (1u64 << n) | (1u64 << k1) | (1u64 << k2) | (1u64 << k3) | 1 };
                    if f.is_irreducible() {
                        return f;
                    }
                }
            }
        }
        panic!("no irreducible found for n={}", n);
    }

    /// For prime n: f irreducible iff x^(2^n) ≡ x (mod f) and f has no root.
    fn is_irreducible(&self) -> bool {
        let mut x = 2u64;
        for _ in 0..self.n {
            x = self.mul(x, x);
        }
        x == 2 && (self.poly.count_ones() % 2 == 1)
    }

    #[inline]
    fn mul(&self, mut a: u64, mut b: u64) -> u64 {
        let n = self.n;
        let mut r = 0u64;
        while b != 0 {
            if b & 1 == 1 {
                r ^= a;
            }
            b >>= 1;
            a <<= 1;
            if (a >> n) & 1 == 1 {
                a ^= self.poly;
            }
        }
        r
    }

    #[inline]
    fn sq(&self, a: u64) -> u64 {
        self.mul(a, a)
    }

    fn inv(&self, a: u64) -> u64 {
        assert!(a != 0);
        // a^(2^n - 2)
        let mut r = 1u64;
        let mut base = a;
        let mut e = (1u64 << self.n) - 2;
        while e != 0 {
            if e & 1 == 1 {
                r = self.mul(r, base);
            }
            base = self.mul(base, base);
            e >>= 1;
        }
        r
    }

    fn tr(&self, mut a: u64) -> u64 {
        let mut s = 0;
        for _ in 0..self.n {
            s ^= a;
            a = self.sq(a);
        }
        s
    }

    /// n odd: z with z² + z = a when Tr(a) = 0.
    fn htr(&self, mut a: u64) -> u64 {
        let mut s = 0;
        for _ in 0..(self.n + 1) / 2 {
            s ^= a;
            a = self.sq(self.sq(a));
        }
        s
    }

    fn canon(&self, mut x: u64) -> u64 {
        let mut m = x;
        for _ in 1..self.n {
            x = self.sq(x);
            if x < m {
                m = x;
            }
        }
        m
    }
}

// ---------------------------------------------------------------- curve

type Pt = Option<(u64, u64)>;

fn neg(p: Pt) -> Pt {
    p.map(|(x, y)| (x, x ^ y))
}

fn add(f: &Field, p: Pt, q: Pt) -> Pt {
    match (p, q) {
        (None, q) => q,
        (p, None) => p,
        (Some((x1, y1)), Some((x2, y2))) => {
            if x1 == x2 {
                if y1 ^ y2 == x1 || x1 == 0 {
                    return None;
                }
                let lam = x1 ^ f.mul(y1, f.inv(x1));
                let x3 = f.sq(lam) ^ lam;
                Some((x3, f.sq(x1) ^ f.mul(lam ^ 1, x3)))
            } else {
                let lam = f.mul(y1 ^ y2, f.inv(x1 ^ x2));
                let x3 = f.sq(lam) ^ lam ^ x1 ^ x2;
                Some((x3, f.mul(lam, x1 ^ x3) ^ x3 ^ y1))
            }
        }
    }
}

fn sub(f: &Field, p: Pt, q: Pt) -> Pt {
    add(f, p, neg(q))
}

fn frob(f: &Field, p: Pt) -> Pt {
    p.map(|(x, y)| (f.sq(x), f.sq(y)))
}

fn smul(f: &Field, mut k: u64, p: Pt) -> Pt {
    let mut r: Pt = None;
    let mut b = p;
    while k != 0 {
        if k & 1 == 1 {
            r = add(f, r, b);
        }
        b = add(f, b, b);
        k >>= 1;
    }
    r
}

fn on_curve(f: &Field, p: Pt) -> bool {
    match p {
        None => true,
        Some((x, y)) => f.sq(y) ^ f.mul(x, y) == f.mul(f.sq(x), x) ^ 1,
    }
}

fn lift(f: &Field, x: u64) -> Pt {
    if x == 0 {
        return Some((0, 1));
    }
    let k = x ^ f.inv(f.sq(x));
    if f.tr(k) != 0 {
        return None;
    }
    Some((x, f.mul(x, f.htr(k))))
}

// ---------------------------------------------------------------- Z/r

fn mulmod(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 * b as u128) % r as u128) as u64
}

fn powmod(mut a: u64, mut e: u64, r: u64) -> u64 {
    let mut x = 1u64;
    a %= r;
    while e != 0 {
        if e & 1 == 1 {
            x = mulmod(x, a, r);
        }
        a = mulmod(a, a, r);
        e >>= 1;
    }
    x
}

fn invmod(a: u64, r: u64) -> u64 {
    powmod(a, r - 2, r)
}

fn negmod(a: u64, r: u64) -> u64 {
    if a == 0 { 0 } else { r - a }
}

/// Tonelli–Shanks square root mod prime r; None if non-residue.
fn sqrtmod(a: u64, r: u64) -> Option<u64> {
    let a = a % r;
    if a == 0 {
        return Some(0);
    }
    if powmod(a, (r - 1) / 2, r) != 1 {
        return None;
    }
    if r % 4 == 3 {
        return Some(powmod(a, (r + 1) / 4, r));
    }
    let mut q = r - 1;
    let mut s = 0;
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while powmod(z, (r - 1) / 2, r) != r - 1 {
        z += 1;
    }
    let mut m = s;
    let mut c = powmod(z, q, r);
    let mut t = powmod(a, q, r);
    let mut x = powmod(a, (q + 1) / 2, r);
    while t != 1 {
        let mut i = 0;
        let mut tt = t;
        while tt != 1 {
            tt = mulmod(tt, tt, r);
            i += 1;
        }
        let mut b = c;
        for _ in 0..(m - i - 1) {
            b = mulmod(b, b, r);
        }
        m = i;
        c = mulmod(b, b, r);
        t = mulmod(t, c, r);
        x = mulmod(x, b, r);
    }
    Some(x)
}

// ---------------------------------------------------------------- rng

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

// ---------------------------------------------------------------- group order

fn koblitz_order(n: u32) -> u64 {
    // tau + taubar = -1, tau*taubar = 2; V_k = tau^k + taubar^k.
    let (mut v0, mut v1): (i128, i128) = (2, -1);
    for _ in 2..=n {
        let v2 = -v1 - 2 * v0;
        v0 = v1;
        v1 = v2;
    }
    ((1i128 << n) + 1 - v1) as u64
}

fn factor(mut v: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut d = 2u64;
    while d * d <= v {
        let mut e = 0;
        while v % d == 0 {
            v /= d;
            e += 1;
        }
        if e > 0 {
            out.push((d, e));
        }
        d += 1;
    }
    if v > 1 {
        out.push((v, 1));
    }
    out
}

// ---------------------------------------------------------------- multipliers

#[derive(Clone, Copy, PartialEq, Debug)]
enum Kind {
    Frob,
    Tau,
    Int,
}

/// Elements a + bτ of Z[τ], τ² + τ + 2 = 0, norm a² − ab + 2b², with a odd
/// (τ-indivisible) and a > 0 (one of ±α), the c smallest by norm.
fn tau_set(c: usize) -> Vec<(i64, i64)> {
    let mut v = Vec::new();
    let bound = 4 * (c as i64) + 8;
    for a in (1..=bound).step_by(2) {
        for b in -bound..=bound {
            let norm = a * a - a * b + 2 * b * b;
            if norm <= bound {
                v.push((norm, b.abs(), b, a));
            }
        }
    }
    v.sort();
    v.iter().take(c).map(|&(_, _, b, a)| (a, b)).collect()
}

struct Arm {
    kind: Kind,
    k: usize,
    c: usize,
}

// ---------------------------------------------------------------- main

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 6 {
        eprintln!("usage: multipliers <n> <seed> <Kc> <hits> <out.json>");
        std::process::exit(2);
    }
    let n: u32 = args[1].parse().unwrap();
    let seed: u64 = args[2].parse().unwrap();
    let kc: usize = args[3].parse().unwrap();
    let hits_wanted: usize = args[4].parse().unwrap();
    let out_path = &args[5];

    let f = Field::new(n);
    let order = koblitz_order(n);
    let fac = factor(order);
    let r = fac.last().unwrap().0;
    let cof = order / r;
    assert!(fac.last().unwrap().1 == 1 && r > 1 << 16, "largest prime factor must be simple and large");
    eprintln!("n={n} poly={:#x} #E={order} = {:?}  r={r} cofactor={cof}", f.poly, fac);

    let mut rng = Rng(seed ^ 0xA5A5_0000_0000 ^ n as u64);

    // Generator of the order-r subgroup.
    let g = loop {
        let x = rng.below(1u64 << n);
        if let Some(p) = lift(&f, x) {
            let q = smul(&f, cof, Some(p));
            if q.is_some() {
                assert!(smul(&f, r, q).is_none());
                break q;
            }
        }
    };
    assert!(on_curve(&f, g));

    // Frobenius eigenvalue mu on <G>: mu^2 + mu + 2 = 0 mod r, tau(G) = [mu]G.
    let s = sqrtmod(r - 7 % r, r).expect("-7 must be a QR mod r");
    let inv2 = invmod(2, r);
    let cands = [mulmod((s + r - 1) % r, inv2, r), mulmod((negmod(s, r) + r - 1) % r, inv2, r)];
    let tg = frob(&f, g);
    let mu = *cands.iter().find(|&&m| smul(&f, m, g) == tg).expect("no eigenvalue matched");
    assert_eq!(mulmod(mu, mu, r) + mu + 2 - r * ((mulmod(mu, mu, r) + mu + 2) / r), 0);
    assert_eq!(powmod(mu, n as u64, r), 1);

    // Arms: all divisors K of Kc with K*c = Kc, kinds tau and int, plus frob at K = Kc.
    let mut arms = vec![Arm { kind: Kind::Frob, k: kc, c: 1 }];
    let mut ks: Vec<usize> = (1..=kc).filter(|k| kc % k == 0 && *k < kc).collect();
    ks.sort_by(|a, b| b.cmp(a));
    // keep at most 4 K values below Kc: largest, two middle, and 1
    if ks.len() > 4 {
        let l = ks.len();
        ks = vec![ks[0], ks[l / 3], ks[2 * l / 3], ks[l - 1]];
    }
    for &k in &ks {
        arms.push(Arm { kind: Kind::Tau, k, c: kc / k });
        arms.push(Arm { kind: Kind::Int, k, c: kc / k });
    }

    let mut results: Vec<String> = Vec::new();
    for arm in &arms {
        let t0 = Instant::now();
        let kind_s = format!("{:?}", arm.kind).to_lowercase();
        eprintln!("--- arm {kind_s} K={} c={}", arm.k, arm.c);

        // Multiplier set as (scalar mod r, description) built faithfully from the point.
        let alphas: Vec<(i64, i64)> = match arm.kind {
            Kind::Frob => vec![(1, 0)],
            Kind::Tau => tau_set(arm.c),
            Kind::Int => (1..=arm.c as i64).map(|a| (a, 0)).collect(),
        };
        // hidden column scalars
        let cols: Vec<u64> = (0..arm.k).map(|_| 1 + rng.below(r - 1)).collect();
        let col_pts: Vec<Pt> = cols.iter().map(|&c| smul(&f, c, g)).collect();

        // S: summands (point, column, scalar multiplier mod r)
        let mut summands: Vec<(Pt, usize, u64)> = Vec::new();
        let mut seen: HashMap<u64, ()> = HashMap::new();
        let mut dup_scalars = 0usize;
        let mut mismatches = 0usize;
        for (i, &p) in col_pts.iter().enumerate() {
            for &(a, b) in &alphas {
                // alpha(P) = a P + b tau(P), computed from the point alone
                let ap = smul(&f, a as u64, p);
                let bp = if b >= 0 { smul(&f, b as u64, frob(&f, p)) } else { neg(smul(&f, (-b) as u64, frob(&f, p))) };
                let mut q = add(&f, ap, bp);
                let alpha_scalar = ((a as i128 + b as i128 * mu as i128).rem_euclid(r as i128)) as u64;
                // self-check: alpha(P) == [alpha(mu)] P
                if smul(&f, alpha_scalar, p) != q {
                    mismatches += 1;
                }
                let mut m = alpha_scalar;
                for _ in 0..n {
                    for sign in [false, true] {
                        let (pt, sc) = if sign { (neg(q), negmod(m, r)) } else { (q, m) };
                        let key = mulmod(sc, cols[i], r);
                        if seen.insert(key, ()).is_some() {
                            dup_scalars += 1;
                            continue;
                        }
                        summands.push((pt, i, sc));
                    }
                    q = frob(&f, q);
                    m = mulmod(m, mu, r);
                }
            }
        }
        assert_eq!(mismatches, 0, "alpha(mu) eigenvalue self-check failed");
        let l = summands.len();
        eprintln!("|S|={l} (2nKc={}) dup={dup_scalars}", 2 * n as usize * arm.k * arm.c);

        // Pair table: (x(s_a+s_b), a, b) for a < b, batched inversion.
        let t1 = Instant::now();
        let mut table: Vec<(u64, u32, u32)> = Vec::with_capacity(l * (l - 1) / 2);
        let mut skipped = 0usize;
        // batch over a: for fixed a, all b > a; denominators x_a ^ x_b
        for a in 0..l {
            let (xa, ya) = summands[a].0.unwrap();
            let mut dens: Vec<u64> = Vec::with_capacity(l - a);
            let mut idx: Vec<usize> = Vec::with_capacity(l - a);
            for b in (a + 1)..l {
                let (xb, yb) = summands[b].0.unwrap();
                if xa == xb {
                    if ya ^ yb == xa {
                        skipped += 1; // s_a + s_b = O
                    } else {
                        // doubling (same point twice: cannot happen, scalars deduped)
                        let q = add(&f, summands[a].0, summands[b].0);
                        table.push((q.unwrap().0, a as u32, b as u32));
                    }
                    continue;
                }
                dens.push(xa ^ xb);
                idx.push(b);
            }
            // Montgomery batch inversion
            let m = dens.len();
            if m == 0 {
                continue;
            }
            let mut pref = vec![0u64; m];
            let mut acc = 1u64;
            for i in 0..m {
                acc = f.mul(acc, dens[i]);
                pref[i] = acc;
            }
            let mut inv_acc = f.inv(acc);
            for i in (0..m).rev() {
                let inv_i = if i == 0 { inv_acc } else { f.mul(inv_acc, pref[i - 1]) };
                inv_acc = f.mul(inv_acc, dens[i]);
                let b = idx[i];
                let (xb, yb) = summands[b].0.unwrap();
                let lam = f.mul(ya ^ yb, inv_i);
                let x3 = f.sq(lam) ^ lam ^ xa ^ xb;
                table.push((x3, a as u32, b as u32));
            }
        }
        table.sort_unstable();
        let pairs = table.len();
        let table_secs = t1.elapsed().as_secs_f64();

        // memory metric: distinct canonical x under Frobenius
        let t2 = Instant::now();
        let mut canon: Vec<u64> = table.iter().map(|e| f.canon(e.0)).collect();
        canon.sort_unstable();
        canon.dedup();
        let canon_count = canon.len();
        drop(canon);
        let mut distinct_x = 1usize;
        for i in 1..table.len() {
            if table[i].0 != table[i - 1].0 {
                distinct_x += 1;
            }
        }
        let canon_secs = t2.elapsed().as_secs_f64();
        eprintln!("pairs={pairs} distinct_x={distinct_x} canonical={canon_count} (L^2/4n={:.0}) table {:.1}s canon {:.1}s",
            (l * l) as f64 / (4.0 * n as f64), table_secs, canon_secs);

        // relation collection
        let probe = |rng: &mut Rng, target: Pt| -> (u64, Vec<(usize, u64)>) {
            // returns (probes, relation as list of (column, coefficient)) ; target = sum coeff*P_col
            let mut probes = 0u64;
            loop {
                let c = rng.below(l as u64) as usize;
                let d = rng.below(l as u64) as usize;
                if c >= d {
                    continue;
                }
                probes += 1;
                let t = sub(&f, sub(&f, target, summands[c].0), summands[d].0);
                let tx = match t {
                    None => continue,
                    Some((x, _)) => x,
                };
                // binary search
                let lo = table.partition_point(|e| e.0 < tx);
                let mut i = lo;
                while i < table.len() && table[i].0 == tx {
                    let (a, b) = (table[i].1 as usize, table[i].2 as usize);
                    let sab = add(&f, summands[a].0, summands[b].0);
                    let sign_ab = if sab == t { 1u64 } else if sab == neg(t) { r - 1 } else { i += 1; continue };
                    let mut rel: Vec<(usize, u64)> = vec![
                        (summands[a].1, mulmod(sign_ab, summands[a].2, r)),
                        (summands[b].1, mulmod(sign_ab, summands[b].2, r)),
                        (summands[c].1, summands[c].2),
                        (summands[d].1, summands[d].2),
                    ];
                    rel.sort();
                    return (probes, rel);
                }
            }
        };

        let t3 = Instant::now();
        let k = arm.k;
        // echelon over Z/r: rows of length k+1 (last = rhs)
        let mut ech: Vec<Vec<u64>> = Vec::new();
        let mut pivots: Vec<usize> = Vec::new();
        let mut rank = 0usize;
        let mut rows_to_full_rank: Option<usize> = None;
        let mut first_hit_probes: Vec<u64> = Vec::new();
        let mut degenerate = 0usize;
        let mut dependent = 0usize;
        let mut rows_seen = 0usize;
        let mut total_probes_to_rank = 0u64;
        let mut touches = vec![0u64; k];
        let mut rank_trace: Vec<usize> = Vec::new();
        for h in 0..hits_wanted {
            let rho = 1 + rng.below(r - 1);
            let target = smul(&f, rho, g);
            let (probes, rel) = probe(&mut rng, target);
            first_hit_probes.push(probes);
            rows_seen += 1;
            // row vector
            let mut row = vec![0u64; k + 1];
            for &(col, coef) in &rel {
                row[col] = (row[col] + coef) % r;
                touches[col] += 1;
            }
            row[k] = rho;
            // sanity against hidden scalars
            let mut lhs = 0u64;
            for i in 0..k {
                lhs = (lhs + mulmod(row[i], cols[i], r)) % r;
            }
            assert_eq!(lhs, rho, "relation does not hold against hidden scalars");
            if row[..k].iter().all(|&v| v == 0) {
                degenerate += 1;
                assert_eq!(rho, 0);
                continue;
            }
            // reduce
            let mut v = row.clone();
            for (pi, prow) in pivots.iter().zip(ech.iter()) {
                if v[*pi] != 0 {
                    let fct = v[*pi];
                    for j in 0..=k {
                        v[j] = (v[j] + r - mulmod(fct, prow[j], r)) % r;
                    }
                }
            }
            if let Some(p) = (0..k).find(|&j| v[j] != 0) {
                let iv = invmod(v[p], r);
                for j in 0..=k {
                    v[j] = mulmod(v[j], iv, r);
                }
                // eliminate from earlier rows
                for (pi, prow) in pivots.iter().zip(ech.iter_mut()) {
                    let _ = pi;
                    if prow[p] != 0 {
                        let fct = prow[p];
                        for j in 0..=k {
                            prow[j] = (prow[j] + r - mulmod(fct, v[j], r)) % r;
                        }
                    }
                }
                pivots.push(p);
                ech.push(v);
                rank += 1;
                if rank == k && rows_to_full_rank.is_none() {
                    rows_to_full_rank = Some(rows_seen);
                    total_probes_to_rank = first_hit_probes.iter().sum();
                }
            } else {
                dependent += 1;
            }
            rank_trace.push(rank);
            let _ = h;
        }
        let untouched = touches.iter().filter(|&&t| t == 0).count();
        eprintln!("column touches min={} max={} untouched={untouched}", touches.iter().min().unwrap(), touches.iter().max().unwrap());
        let collect_secs = t3.elapsed().as_secs_f64();

        // solve and verify columns
        let mut solved = vec![0u64; k];
        let mut logs_ok = false;
        if rank == k {
            for (pi, prow) in pivots.iter().zip(ech.iter()) {
                solved[*pi] = prow[k];
            }
            logs_ok = solved == cols;
        }

        // online: hidden target
        let d_hidden = 1 + rng.below(r - 1);
        let q = smul(&f, d_hidden, g);
        let (online_probes, rel) = probe(&mut rng, q);
        let mut d_rec = 0u64;
        for &(col, coef) in &rel {
            d_rec = (d_rec + mulmod(coef, solved[col], r)) % r;
        }
        let target_ok = logs_ok && d_rec == d_hidden;

        let mean_probes = first_hit_probes.iter().sum::<u64>() as f64 / first_hit_probes.len() as f64;
        let mut sorted = first_hit_probes.clone();
        sorted.sort();
        let median = sorted[sorted.len() / 2];
        let model = 6.0 * r as f64 / (l as f64 * l as f64);
        eprintln!("hits={} mean_probes={mean_probes:.1} median={median} model 6r/L^2={model:.1} rank={rank}/{k} rows_to_rank={:?} degenerate={degenerate} dependent={dependent} logs_ok={logs_ok} target_ok={target_ok} online_probes={online_probes} {:.1}s",
            first_hit_probes.len(), rows_to_full_rank, t0.elapsed().as_secs_f64());

        let alphas_s: Vec<String> = alphas.iter().map(|(a, b)| format!("[{a},{b}]")).collect();
        results.push(format!(
            "{{\"kind\":\"{kind_s}\",\"K\":{k},\"c\":{},\"alphas\":[{}],\"L\":{l},\"L_nominal\":{},\"dup_scalars\":{dup_scalars},\"pairs\":{pairs},\"pairs_skipped\":{skipped},\"distinct_x\":{distinct_x},\"canonical\":{canon_count},\"canonical_model\":{:.1},\"hits\":{},\"first_hit_probes\":[{}],\"mean_first_hit_probes\":{mean_probes:.3},\"median_first_hit_probes\":{median},\"model_6r_over_L2\":{model:.3},\"rank\":{rank},\"rows_to_full_rank\":{},\"total_probes_to_full_rank\":{total_probes_to_rank},\"degenerate_rows\":{degenerate},\"dependent_rows\":{dependent},\"column_touches\":[{}],\"rank_trace\":[{}],\"logs_ok\":{logs_ok},\"target_ok\":{target_ok},\"online_probes\":{online_probes},\"secs\":{{\"table\":{table_secs:.2},\"canon\":{canon_secs:.2},\"collect\":{collect_secs:.2},\"arm\":{:.2}}}}}",
            arm.c,
            alphas_s.join(","),
            2 * n as usize * arm.k * arm.c,
            (l * l) as f64 / (4.0 * n as f64),
            first_hit_probes.len(),
            first_hit_probes.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(","),
            rows_to_full_rank.map(|v| v.to_string()).unwrap_or("null".into()),
            touches.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(","),
            rank_trace.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(","),
            t0.elapsed().as_secs_f64(),
        ));
    }

    let fac_s: Vec<String> = fac.iter().map(|(p, e)| format!("[{p},{e}]")).collect();
    let json = format!(
        "{{\"n\":{n},\"poly\":{},\"order\":{order},\"factorization\":[{}],\"r\":{r},\"cofactor\":{cof},\"mu\":{mu},\"seed\":{seed},\"Kc\":{kc},\"hits_wanted\":{hits_wanted},\"arms\":[\n{}\n]}}\n",
        f.poly,
        fac_s.join(","),
        results.join(",\n")
    );
    let mut fh = std::fs::File::create(out_path).unwrap();
    fh.write_all(json.as_bytes()).unwrap();
    eprintln!("wrote {out_path}");
}
