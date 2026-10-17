// 01-endo-chain: factor bases closed under the degree-5 CM endomorphism w=(1+sqrt(-19))/2.
// Standalone instrument (see PROTOCOL.md). All arithmetic mod p < 2^63 via u128.
use std::collections::HashMap;
use std::time::Instant;

// ---------------- RNG ----------------
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        ((self.next() as u128 * n as u128) >> 64) as u64
    }
}

// ---------------- modular arithmetic ----------------
#[inline]
fn mulm(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}
#[inline]
fn addm(a: u64, b: u64, m: u64) -> u64 {
    let s = a as u128 + b as u128;
    (if s >= m as u128 { s - m as u128 } else { s }) as u64
}
#[inline]
fn subm(a: u64, b: u64, m: u64) -> u64 {
    if a >= b { a - b } else { a + (m - b) }
}
fn powm(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1 % m;
    b %= m;
    while e > 0 {
        if e & 1 == 1 { r = mulm(r, b, m); }
        b = mulm(b, b, m);
        e >>= 1;
    }
    r
}
fn invm(a: u64, m: u64) -> u64 {
    let (mut t, mut nt) = (0i128, 1i128);
    let (mut r, mut nr) = (m as i128, (a % m) as i128);
    while nr != 0 {
        let q = r / nr;
        let tmp = t - q * nt; t = nt; nt = tmp;
        let tmp = r - q * nr; r = nr; nr = tmp;
    }
    assert!(r == 1, "not invertible");
    if t < 0 { t += m as i128; }
    t as u64
}
fn is_prime(n: u64) -> bool {
    if n < 2 { return false; }
    for &q in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 { return n == q; }
    }
    let mut d = n - 1; let mut s = 0;
    while d % 2 == 0 { d /= 2; s += 1; }
    'w: for &a in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powm(a, d, n);
        if x == 1 || x == n - 1 { continue; }
        for _ in 1..s { x = mulm(x, x, n); if x == n - 1 { continue 'w; } }
        return false;
    }
    true
}
fn legendre(a: u64, p: u64) -> i32 {
    if a % p == 0 { return 0; }
    if powm(a, (p - 1) / 2, p) == 1 { 1 } else { -1 }
}
fn sqrtm(a: u64, p: u64) -> Option<u64> {
    let a = a % p;
    if a == 0 { return Some(0); }
    if legendre(a, p) != 1 { return None; }
    if p % 4 == 3 { return Some(powm(a, (p + 1) / 4, p)); }
    let (mut q, mut s) = (p - 1, 0u32);
    while q % 2 == 0 { q /= 2; s += 1; }
    let mut z = 2; while legendre(z, p) != -1 { z += 1; }
    let (mut m, mut c, mut t, mut r) = (s, powm(z, q, p), powm(a, q, p), powm(a, (q + 1) / 2, p));
    while t != 1 {
        let mut i = 0; let mut tt = t;
        while tt != 1 { tt = mulm(tt, tt, p); i += 1; }
        let b = powm(c, 1u64 << (m - i - 1), p);
        m = i; c = mulm(b, b, p); t = mulm(t, c, p); r = mulm(r, b, p);
    }
    Some(r)
}
fn factor(mut n: u64) -> Vec<u64> {
    let mut f = vec![]; let mut q = 2u64;
    while q * q <= n { if n % q == 0 { f.push(q); while n % q == 0 { n /= q; } } q += if q == 2 { 1 } else { 2 }; }
    if n > 1 { f.push(n); }
    f
}
fn mult_order(a: u64, n: u64) -> u64 {
    // n prime; order of a in (Z/n)^*
    let mut o = n - 1;
    for q in factor(n - 1) { while o % q == 0 && powm(a, o / q, n) == 1 { o /= q; } }
    o
}

// ---------------- curve ----------------
#[derive(Clone, Copy)]
struct Curve { p: u64, a: u64, b: u64 }
type Pt = Option<(u64, u64)>;
#[derive(Default, Clone, Copy, Debug)]
struct Ctr { a: u64, e: u64, l: u64 }

impl Curve {
    fn rhs(&self, x: u64) -> u64 { let p = self.p; addm(addm(mulm(mulm(x, x, p), x, p), mulm(self.a, x, p), p), self.b, p) }
    fn neg(&self, q: Pt) -> Pt { q.map(|(x, y)| (x, if y == 0 { 0 } else { self.p - y })) }
    fn add(&self, p1: Pt, p2: Pt, c: &mut Ctr) -> Pt {
        c.a += 1;
        let p = self.p;
        match (p1, p2) {
            (None, q) | (q, None) => q,
            (Some((x1, y1)), Some((x2, y2))) => {
                let l = if x1 == x2 {
                    if addm(y1, y2, p) == 0 { return None; }
                    mulm(addm(mulm(3, mulm(x1, x1, p), p), self.a, p), invm(mulm(2, y1, p), p), p)
                } else {
                    mulm(subm(y2, y1, p), invm(subm(x2, x1, p), p), p)
                };
                let x3 = subm(subm(mulm(l, l, p), x1, p), x2, p);
                let y3 = subm(mulm(l, subm(x1, x3, p), p), y1, p);
                Some((x3, y3))
            }
        }
    }
    fn mul(&self, k: u64, q: Pt, c: &mut Ctr) -> Pt {
        let mut r: Pt = None;
        if k == 0 { return None; }
        let bits = 64 - k.leading_zeros();
        for i in (0..bits).rev() {
            if r.is_some() { r = self.add(r, r, c); }
            if (k >> i) & 1 == 1 { r = if r.is_none() { q } else { self.add(r, q, c) }; }
        }
        r
    }
    fn rand_point(&self, rng: &mut Rng) -> Pt {
        loop {
            let x = rng.below(self.p);
            if let Some(y) = sqrtm(self.rhs(x), self.p) { if y != 0 { return Some((x, if rng.next() & 1 == 0 { y } else { self.p - y })); } }
        }
    }
}

// x-map of the endomorphism: R(x) = Num(x)/Den(x), Num deg 5, Den monic deg 4
#[derive(Clone)]
struct XMap { num: [u64; 6], den: [u64; 5], p: u64 }
impl XMap {
    #[inline]
    fn eval(&self, x: u64, c: &mut Ctr) -> Option<u64> {
        c.e += 1;
        let p = self.p;
        let mut n = self.num[5];
        for i in (0..5).rev() { n = addm(mulm(n, x, p), self.num[i], p); }
        let mut d = self.den[4];
        for i in (0..4).rev() { d = addm(mulm(d, x, p), self.den[i], p); }
        if d == 0 { return None; }
        Some(mulm(n, invm(d, p), p))
    }
}

fn solve_linear(mut m: Vec<Vec<u64>>, p: u64) -> Option<Vec<u64>> {
    let n = m.len();
    for col in 0..n {
        let piv = (col..n).find(|&r| m[r][col] != 0)?;
        m.swap(col, piv);
        let iv = invm(m[col][col], p);
        for j in col..=n { m[col][j] = mulm(m[col][j], iv, p); }
        for r in 0..n {
            if r != col && m[r][col] != 0 {
                let f = m[r][col];
                for j in col..=n { let t = mulm(f, m[col][j], p); m[r][j] = subm(m[r][j], t, p); }
            }
        }
    }
    Some(m.iter().map(|row| row[n]).collect())
}

// interpolate x(lambda P) = Num/Den, verify on `verify` fresh points
fn interpolate_xmap(e: &Curve, n_ord: u64, lam: u64, rng: &mut Rng, verify: usize) -> Option<XMap> {
    let p = e.p; let mut c = Ctr::default();
    let mut rows = vec![];
    while rows.len() < 10 {
        let pt = e.rand_point(rng);
        let q = e.mul(lam, pt, &mut c);
        let (x, _) = pt.unwrap(); let z = q?.0;
        // sum n_j x^j - z sum_{j<4} d_j x^j = z x^4
        let mut row = vec![0u64; 11]; let mut xp = 1u64;
        for j in 0..6 { row[j] = xp; if j < 4 { row[6 + j] = subm(0, mulm(z, xp, p), p); } xp = mulm(xp, x, p); }
        row[10] = mulm(z, powm(x, 4, p), p);
        rows.push(row);
    }
    let sol = solve_linear(rows, p)?;
    let mut num = [0u64; 6]; let mut den = [0u64; 5];
    num.copy_from_slice(&sol[0..6]); den[..4].copy_from_slice(&sol[6..10]); den[4] = 1;
    let xm = XMap { num, den, p };
    let _ = n_ord;
    for _ in 0..verify {
        let pt = e.rand_point(rng);
        let q = e.mul(lam, pt, &mut c);
        let want = q.map(|v| v.0);
        if xm.eval(pt.unwrap().0, &mut c) != want { return None; }
    }
    if num[5] == 0 { return None; }
    Some(xm)
}

// ---------------- toy curve construction ----------------
struct Toy { bits: u32, t: i64, e: Curve, n: u64, lam: u64, xm: XMap, ord_lam: u64, twist: bool }

fn build_toy(bits: u32, rng: &mut Rng) -> Toy {
    let mut t: u64 = 2;
    loop {
        let p = 5 * t * t + 19 * t + 19;
        if 64 - p.leading_zeros() > bits { panic!("no curve at {bits}"); }
        if 64 - p.leading_zeros() == bits && is_prime(p) {
            let j = p - (884736 % p);
            let k = mulm(j, invm(subm(1728 % p, j, p), p), p);
            let e0 = Curve { p, a: mulm(3, k, p), b: mulm(2, k, p) };
            let mut nr = 2; while legendre(nr, p) != -1 { nr += 1; }
            let e1 = Curve { p, a: mulm(e0.a, mulm(nr, nr, p), p), b: mulm(e0.b, powm(nr, 3, p), p) };
            for (tw, e) in [(false, e0), (true, e1)] {
                for n in [p + 1 - t, p + 1 + t] {
                    if !is_prime(n) { continue; }
                    let mut c = Ctr::default();
                    let pt = e.rand_point(rng);
                    if e.mul(n, pt, &mut c).is_some() { continue; }
                    // lambda: roots of X^2 - X + 5 mod n
                    let s = match sqrtm(n - 19 % n, n) { Some(s) => s, None => continue };
                    let i2 = invm(2, n);
                    for lam in [mulm(addm(1, s, n), i2, n), mulm(subm(1, s, n), i2, n)] {
                        if let Some(xm) = interpolate_xmap(&e, n, lam, rng, 200) {
                            let ord_lam = mult_order(lam, n);
                            return Toy { bits, t: t as i64, e, n, lam, xm, ord_lam, twist: tw };
                        }
                    }
                }
            }
        }
        t += 1;
    }
}

// ---------------- timing (phase 0) ----------------
fn median(v: &mut Vec<f64>) -> f64 { v.sort_by(|a, b| a.partial_cmp(b).unwrap()); v[v.len() / 2] }

fn calibrate(e: &Curve, xm: &XMap, rng: &mut Rng, reps: usize, ops: usize) -> (f64, f64, f64, f64, f64) {
    let base = e.rand_point(rng);
    let (mut ta, mut te, mut tl) = (vec![], vec![], vec![]);
    let mut sink = 0u64;
    for _ in 0..reps {
        let mut c = Ctr::default();
        let mut r = e.rand_point(rng);
        let t0 = Instant::now();
        for _ in 0..ops { r = e.add(r, base, &mut c); }
        ta.push(t0.elapsed().as_nanos() as f64 / ops as f64);
        sink ^= r.map(|v| v.0).unwrap_or(0);
        let mut x = e.rand_point(rng).unwrap().0;
        let t0 = Instant::now();
        for _ in 0..ops { x = xm.eval(x, &mut c).unwrap_or(1); }
        te.push(t0.elapsed().as_nanos() as f64 / ops as f64);
        sink ^= x;
        let mut x = rng.below(e.p); let mut acc = 0i64;
        let t0 = Instant::now();
        for _ in 0..ops / 10 { acc += legendre(e.rhs(x), e.p) as i64; x = addm(x, 1, e.p); }
        tl.push(t0.elapsed().as_nanos() as f64 / (ops / 10) as f64);
        sink ^= acc as u64;
    }
    if sink == 42 { eprintln!("."); }
    let (ma, me, ml) = (median(&mut ta), median(&mut te), median(&mut tl));
    (ma, me, ml, me / ma, ml / ma)
}

// ---------------- arms ----------------
struct Out { a: u64, e: u64, l: u64, ok: bool }

fn verify_d(e: &Curve, d: u64, p0: Pt, q: Pt, c: &mut Ctr) -> bool { e.mul(d, p0, c) == q }

fn rho(e: &Curve, n: u64, p0: Pt, q: Pt, rng: &mut Rng, negation: bool, cap: u64) -> (Out, u64) {
    let mut c = Ctr::default();
    let r = 20usize;
    let mut ms = vec![];
    for _ in 0..r { let (u, v) = (rng.below(n), rng.below(n)); let m = e.add(e.mul(u, p0, &mut c), e.mul(v, q, &mut c), &mut c); ms.push((m, u, v)); }
    let canon = |pt: Pt, a: u64, b: u64| -> (Pt, u64, u64) {
        if !negation { return (pt, a, b); }
        match pt { Some((x, y)) if y > (e.p - 1) / 2 => (Some((x, e.p - y)), (n - a) % n, (n - b) % n), _ => (pt, a, b) }
    };
    let mut table: HashMap<(u64, u64), (u64, u64)> = HashMap::new();
    let start = |rng: &mut Rng, c: &mut Ctr| { let (a, b) = (rng.below(n), 1 + rng.below(n - 1)); (e.add(e.mul(a, p0, c), e.mul(b, q, c), c), a, b) };
    let (mut pt, mut a, mut b) = start(rng, &mut c);
    let (pt2, a2, b2) = canon(pt, a, b); pt = pt2; a = a2; b = b2;
    let mut steps = 0u64;
    loop {
        if c.a > cap { return (Out { a: c.a, e: 0, l: 0, ok: false }, steps); }
        let key = match pt { Some(v) => v, None => { let s = start(rng, &mut c); let s2 = canon(s.0, s.1, s.2); pt = s2.0; a = s2.1; b = s2.2; continue; } };
        if let Some(&(a1, b1)) = table.get(&key) {
            if b1 != b {
                let d = mulm(subm(a, a1, n), invm(subm(b1, b, n), n), n);
                let ok = verify_d(e, d, p0, q, &mut c);
                return (Out { a: c.a, e: 0, l: 0, ok }, steps);
            }
            // useless (fruitless cycle or same rep): restart
            let s = start(rng, &mut c); let s2 = canon(s.0, s.1, s.2); pt = s2.0; a = s2.1; b = s2.2;
            continue;
        }
        table.insert(key, (a, b));
        let h = (key.0 as usize).wrapping_mul(0x9E37) % r;
        let (m, u, v) = ms[h];
        let np = e.add(pt, m, &mut c);
        let (np2, na, nb) = canon(np, addm(a, u, n), addm(b, v, n));
        pt = np2; a = na; b = nb; steps += 1;
    }
}

fn bsgs_pm(e: &Curve, n: u64, p0: Pt, q: Pt, cap: u64) -> Out {
    let mut c = Ctr::default();
    let m = ((n as f64).sqrt() / 2.0).ceil() as u64;
    let mut tab: HashMap<u64, u64> = HashMap::new();
    let mut r = p0;
    for j in 1..=m { tab.insert(r.unwrap().0, j); if j < m { r = e.add(r, p0, &mut c); } }
    let stride = e.mul(2 * m + 1, p0, &mut c);
    let mstride = e.neg(stride);
    let mut g = q; let mut i = 0u64;
    loop {
        if c.a > cap { return Out { a: c.a, e: 0, l: 0, ok: false }; }
        if g.is_none() { let d = mulm(i % n, (2 * m + 1) % n, n); let ok = verify_d(e, d, p0, q, &mut c); return Out { a: c.a, e: 0, l: 0, ok }; }
        if let Some(&j) = tab.get(&g.unwrap().0) {
            let base = mulm(i % n, (2 * m + 1) % n, n);
            for d in [addm(base, j, n), subm(base, j, n)] { if verify_d(e, d, p0, q, &mut c) { return Out { a: c.a, e: 0, l: 0, ok: true }; } }
        }
        g = e.add(g, mstride, &mut c); i += 1;
    }
}

// table of chains: x -> (chain id, exponent i). heads at smallest on-curve x.
fn chain_table(e: &Curve, xm: &XMap, total: u64, k: u64, c: &mut Ctr) -> (HashMap<u64, (u32, u32)>, u64) {
    let mut tab: HashMap<u64, (u32, u32)> = HashMap::new();
    let mut x0 = 0u64; let mut chains = 0u32; let mut dup = 0u64;
    while (tab.len() as u64) < total {
        // next head
        loop { c.l += 1; let r = e.rhs(x0); x0 += 1; if legendre(r, e.p) == 1 { break; } }
        let head = x0 - 1;
        if tab.contains_key(&head) { dup += 1; continue; }
        let mut x = head; let cid = chains; chains += 1;
        for i in 0..k {
            if (tab.len() as u64) >= total { break; }
            if tab.contains_key(&x) { dup += 1; break; }
            tab.insert(x, (cid, i as u32));
            if i + 1 < k { match xm.eval(x, c) { Some(nx) => x = nx, None => break } }
        }
    }
    (tab, dup)
}

// walk R_j = Q + (a0+j)P; decoder: fold=true pairs hits on the same chain, fold=false same element only
fn chain_solve(e: &Curve, n: u64, lam: u64, p0: Pt, q: Pt, tab: &HashMap<u64, (u32, u32)>, fold: bool, rng: &mut Rng, c: &mut Ctr, cap: u64) -> bool {
    let a0 = rng.below(n);
    let mut r = e.add(q, e.mul(a0, p0, c), c);
    let mut a = a0;
    let mut hits: HashMap<(u32, u32), Vec<(u64, u32)>> = HashMap::new();
    // cache lambda powers
    let mut lp: HashMap<u32, u64> = HashMap::new();
    let mut lpow = |i: u32| -> u64 { *lp.entry(i).or_insert_with(|| powm(lam, i as u64, n)) };
    loop {
        if c.a > cap { return false; }
        if let Some((x, _)) = r {
            if let Some(&(cid, i)) = tab.get(&x) {
                let key = if fold { (cid, 0) } else { (cid, i) };
                let prev = hits.entry(key).or_default().clone();
                for (a1, i1) in prev {
                    // (a1+d) lam^i = sigma (a+d) lam^i1
                    let (l1, l2) = (lpow(i1), lpow(i));
                    for sigma in [1u64, n - 1] {
                        let coef = subm(l2, mulm(sigma, l1, n), n);
                        if coef == 0 { continue; }
                        let rhs = subm(mulm(sigma, mulm(a, l1, n), n), mulm(a1, l2, n), n);
                        let d = mulm(rhs, invm(coef, n), n);
                        if verify_d(e, d, p0, q, c) { return true; }
                    }
                }
                hits.get_mut(&key).unwrap().push((a, i));
            }
        }
        r = e.add(r, p0, c); a = addm(a, 1, n);
    }
}

fn known_chain_solve(e: &Curve, n: u64, lam: u64, xm: &XMap, p0: Pt, q: Pt, k: u64, rng: &mut Rng, cap: u64) -> Out {
    let mut c = Ctr::default();
    let mut tab: HashMap<u64, u32> = HashMap::new();
    let mut x = p0.unwrap().0;
    for i in 0..k { if tab.contains_key(&x) { break; } tab.insert(x, i as u32); if i + 1 < k { match xm.eval(x, &mut c) { Some(v) => x = v, None => break } } }
    let a0 = rng.below(n);
    let mut r = e.add(q, e.mul(a0, p0, &mut c), &mut c); let mut a = a0;
    loop {
        if c.a > cap { return Out { a: c.a, e: c.e, l: c.l, ok: false }; }
        if let Some((xx, _)) = r {
            if let Some(&i) = tab.get(&xx) {
                let li = powm(lam, i as u64, n);
                for d in [subm(li, a, n), subm(n - li, a, n)] { if verify_d(e, d, p0, q, &mut c) { return Out { a: c.a, e: c.e, l: c.l, ok: true }; } }
            }
        }
        r = e.add(r, p0, &mut c); a = addm(a, 1, n);
    }
}

// ---------------- bootstrap ----------------
fn boot_ci(v: &[f64], seed: u64) -> (f64, f64, f64) {
    let m = v.iter().sum::<f64>() / v.len() as f64;
    let mut rng = Rng(seed); let mut ms = vec![];
    for _ in 0..2000 { let mut s = 0.0; for _ in 0..v.len() { s += v[rng.below(v.len() as u64) as usize]; } ms.push(s / v.len() as f64); }
    ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (m, ms[49], ms[1949])
}
fn boot_ratio(a: &[f64], b: &[f64], seed: u64) -> (f64, f64, f64) {
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let r0 = mean(a) / mean(b);
    let mut rng = Rng(seed); let mut rs = vec![];
    for _ in 0..2000 {
        let sa: f64 = (0..a.len()).map(|_| a[rng.below(a.len() as u64) as usize]).sum::<f64>() / a.len() as f64;
        let sb: f64 = (0..b.len()).map(|_| b[rng.below(b.len() as u64) as usize]).sum::<f64>() / b.len() as f64;
        rs.push(sa / sb);
    }
    rs.sort_by(|x, y| x.partial_cmp(y).unwrap());
    (r0, rs[49], rs[1949])
}
fn wilson(k: u64, n: u64) -> (f64, f64, f64) {
    let (k, n) = (k as f64, n as f64); let z = 1.96; let ph = k / n;
    let den = 1.0 + z * z / n; let cen = (ph + z * z / (2.0 * n)) / den;
    let half = z * (ph * (1.0 - ph) / n + z * z / (4.0 * n * n)).sqrt() / den;
    (ph, cen - half, cen + half)
}

fn jarr(v: &[f64]) -> String { format!("[{}]", v.iter().map(|x| format!("{:.6}", x)).collect::<Vec<_>>().join(",")) }


// POST-HOC (not preregistered): random-relation walk. Segments of length seg from fresh random
// aP+bQ, r=20 adding walk with multipliers c_iP+e_iQ; hits decoded with general (a,b).
fn chain_solve_rw(e: &Curve, n: u64, lam: u64, p0: Pt, q: Pt, tab: &HashMap<u64, (u32, u32)>, rng: &mut Rng, c: &mut Ctr, cap: u64, seg: u64, xdet: bool) -> bool {
    let mut ms = vec![];
    for _ in 0..20 { let (u, v) = (rng.below(n), rng.below(n)); let m = e.add(e.mul(u, p0, c), e.mul(v, q, c), c); ms.push((m, u, v)); }
    let mut hits: HashMap<u32, Vec<(u64, u64, u32)>> = HashMap::new();
    loop {
        let (mut a, mut b) = (rng.below(n), 1 + rng.below(n - 1));
        let mut r = e.add(e.mul(a, p0, c), e.mul(b, q, c), c);
        for _ in 0..seg {
            if c.a > cap { return false; }
            let x = match r { Some((x, _)) => x, None => break };
            if let Some(&(cid, i)) = tab.get(&x) {
                let prev = hits.entry(cid).or_default().clone();
                for (a1, b1, i1) in prev {
                    let (l1, l2) = (powm(lam, i1 as u64, n), powm(lam, i as u64, n));
                    // (a1+b1 d) l2 = sigma (a+b d) l1
                    for sigma in [1u64, n - 1] {
                        let coef = subm(mulm(b1, l2, n), mulm(sigma, mulm(b, l1, n), n), n);
                        if coef == 0 { continue; }
                        let rhs = subm(mulm(sigma, mulm(a, l1, n), n), mulm(a1, l2, n), n);
                        let d = mulm(rhs, invm(coef, n), n);
                        if verify_d(e, d, p0, q, c) { return true; }
                    }
                }
                hits.get_mut(&cid).unwrap().push((a, b, i));
            }
            let h = if xdet { (x as usize).wrapping_mul(0x9E37) % 20 } else { rng.below(20) as usize };
            let (m, u, v) = ms[h];
            r = e.add(r, m, c); a = addm(a, u, n); b = addm(b, v, n);
        }
    }
}

fn posthoc(targets: usize, xdet: bool) {
    let mut json = format!("{{\"schema\":\"pair46.endo_chain.posthoc/v1\",\"note\":\"POST-HOC supplementary, not preregistered: random-relation walk replaces the fixed-b walk for the FB arms\",\"x_determined_walk\":{xdet},\"toys\":[\n");
    let mut first = true;
    for &bits in &[16u32, 18, 20, 22] {
        let mut rng = Rng(1000 + bits as u64);
        let toy = build_toy(bits, &mut rng);
        let (e, n, lam) = (toy.e, toy.n, toy.lam);
        let cal = calibrate(&e, &toy.xm, &mut rng, 10, 200000);
        let (c_e, c_l) = (cal.3, cal.4);
        let sq = (n as f64).sqrt(); let f = sq.ceil() as u64; let cap = 40 * n;
        let seg = (sq / 4.0).ceil() as u64;
        let mut trng = Rng(20261015 + bits as u64);
        let mut g = e.rand_point(&mut trng); while g.is_none() { g = e.rand_point(&mut trng); }
        let ds: Vec<u64> = (0..targets).map(|_| 1 + trng.below(n - 1)).collect();
        let mut cc = Ctr::default();
        let qs: Vec<Pt> = ds.iter().map(|&d| e.mul(d, g, &mut cc)).collect();
        let mut arms = vec![];
        let mut arm_rng = Rng(888 + bits as u64);
        for (name, k) in [("abscissa-k1-rw", 1u64), ("endo-chain-k4-rw", 4), ("endo-chain-k16-rw", 16), ("endo-chain-k64-rw", 64), ("endo-chain-kF-rw", f)] {
            let mut tc = Ctr::default();
            let (t, _) = chain_table(&e, &toy.xm, f, k, &mut tc);
            let mut sv = vec![]; let mut nok = 0;
            for (ti, &q) in qs.iter().enumerate() {
                let mut rr = Rng(arm_rng.next() ^ ti as u64);
                let mut c2 = tc;
                if chain_solve_rw(&e, n, lam, g, q, &t, &mut rr, &mut c2, cap, if xdet { seg } else { u64::MAX }, xdet) { nok += 1; }
                sv.push((c2.a as f64 + c_e * c2.e as f64 + c_l * c2.l as f64) / sq);
            }
            let (m, lo, hi) = boot_ci(&sv, 7);
            println!("posthoc b={bits} {name:20} ok={nok}/{targets} S={m:.3} [{lo:.3},{hi:.3}]");
            arms.push(format!("{{\"arm\":\"{name}\",\"k\":{k},\"verified\":{nok},\"n\":{targets},\"S_mean\":{m:.5},\"S_ci95\":[{lo:.5},{hi:.5}],\"S\":{}}}", jarr(&sv)));
        }
        if !first { json += ",\n"; } first = false;
        json += &format!("{{\"bits\":{bits},\"N\":{n},\"F\":{f},\"segment\":{seg},\"c_E\":{c_e:.4},\"arms\":[{}]}}", arms.join(",\n"));
    }
    json += "\n]}\n";
    std::fs::write(if xdet { "posthoc_v1_xdetermined.json" } else { "posthoc_v2_randomwalk.json" }, json).unwrap();
}

// derived ratios from results.json (recorded per-target S arrays), bootstrap seed 19
fn arrays_by_bits(txt: &str) -> Vec<(u32, Vec<(String, Vec<f64>)>)> {
    let mut out = vec![];
    for chunk in txt.split("{\"bits\":").skip(1) {
        let bits: u32 = chunk.split(',').next().unwrap().trim().parse().unwrap();
        let mut arms = vec![];
        for a in chunk.split("{\"arm\":\"").skip(1) {
            let name = a.split('"').next().unwrap().to_string();
            let sidx = a.find("\"S\":[").unwrap() + 5;
            let end = a[sidx..].find(']').unwrap() + sidx;
            let v: Vec<f64> = a[sidx..end].split(',').map(|x| x.trim().parse().unwrap()).collect();
            arms.push((name, v));
        }
        if !arms.is_empty() { out.push((bits, arms)); }
    }
    out
}
fn ratios_mode() {
    let txt = std::fs::read_to_string("results.json").unwrap();
    let ph = std::fs::read_to_string("posthoc_v2_randomwalk.json").unwrap();
    let mut out = String::from("{\"schema\":\"pair46.endo_chain.derived_ratios/v1\",\"method\":\"ratio of means, percentile bootstrap 95% (2000 resamples, seed 19), from recorded per-target S\",\"rows\":[\n");
    let mut rows = vec![];
    let get = |arms: &Vec<(String, Vec<f64>)>, n: &str| arms.iter().find(|x| x.0 == n).map(|x| x.1.clone());
    let main = arrays_by_bits(&txt); let post = arrays_by_bits(&ph);
    for (bits, arms) in &main {
        for (num, den) in [("endo-chain-known", "bsgs-pm"), ("endo-chain-kF", "bsgs-pm"), ("endo-chain-k16", "bsgs-pm"), ("endo-chain-known", "rho")] {
            let (r, lo, hi) = boot_ratio(&get(arms, num).unwrap(), &get(arms, den).unwrap(), 19);
            rows.push(format!("{{\"bits\":{bits},\"num\":\"{num}\",\"den\":\"{den}\",\"ratio\":[{r:.4},{lo:.4},{hi:.4}]}}"));
            println!("b={bits} {num}/{den} = {r:.3} [{lo:.3},{hi:.3}]");
        }
    }
    for (bits, arms) in &post {
        for num in ["endo-chain-k4-rw", "endo-chain-k16-rw", "endo-chain-k64-rw", "endo-chain-kF-rw"] {
            let (r, lo, hi) = boot_ratio(&get(arms, num).unwrap(), &get(arms, "abscissa-k1-rw").unwrap(), 19);
            rows.push(format!("{{\"bits\":{bits},\"num\":\"{num}\",\"den\":\"abscissa-k1-rw\",\"posthoc\":true,\"ratio\":[{r:.4},{lo:.4},{hi:.4}]}}"));
            println!("posthoc b={bits} {num}/abscissa-k1-rw = {r:.3} [{lo:.3},{hi:.3}]");
        }
    }
    out += &rows.join(",\n"); out += "\n]}\n";
    std::fs::write("derived_ratios.json", out).unwrap();
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(|s| s.as_str()) == Some("ratios") { ratios_mode(); return; }
    if args.get(1).map(|s| s.as_str()) == Some("posthoc") { posthoc(40, args.get(2).map(|s| s.as_str()) == Some("xdet")); return; }
    let targets: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(40);
    let budget_s: f64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1800.0);
    let t_start = Instant::now();
    let mut json = String::from("{\n\"schema\":\"pair46.endo_chain.results/v1\",\n\"protocol_sha256\":\"0ce98a09f63265cb8fda67f1803209c2061476f7fc8116f57239fdb195edfcf2\",\n");

    // ---------- pair 46 ----------
    let mut rng = Rng(46);
    let p46 = 5947698030651343u64; let n46 = 5947697996161643u64; let lam46 = 2973849084305079u64;
    let e46 = Curve { p: p46, a: 5757116383569871, b: 3560173328417442 };
    let mut c = Ctr::default();
    let pt = e46.rand_point(&mut rng);
    let order_ok = e46.mul(n46, pt, &mut c).is_none() && is_prime(n46) && is_prime(p46);
    let lam_root = (mulm(lam46, lam46, n46) + 5) % n46 == lam46 % n46;
    let t46 = 34489701u64;
    let fam = 5 * t46 * t46 + 19 * t46 + 19 == p46 && p46 + 1 - t46 == n46;
    let xm46 = interpolate_xmap(&e46, n46, lam46, &mut rng, 1000);
    let ord46 = mult_order(lam46, n46);
    println!("pair46: order_ok={order_ok} lam_root={lam_root} family={fam} xmap_deg5_verified_1000={} ord(lambda)={ord46}", xm46.is_some());
    let xm46 = xm46.expect("pair46 xmap");
    let cal46 = calibrate(&e46, &xm46, &mut rng, 10, 200000);
    println!("pair46 timing ns: add={:.1} psi_x={:.1} legendre={:.1} c_E={:.3} c_L={:.3}", cal46.0, cal46.1, cal46.2, cal46.3, cal46.4);
    json += &format!("\"pair46\":{{\"p\":{p46},\"N\":{n46},\"lambda\":{lam46},\"order_verified\":{order_ok},\"lambda_root_of_x2_minus_x_plus_5\":{lam_root},\"family_p_5t2_19t_19\":{fam},\"xmap_degree5_verified_points\":1000,\"ord_lambda\":{ord46},\"xmap_num\":{:?},\"xmap_den\":{:?},\"ns_add\":{:.3},\"ns_psi_x\":{:.3},\"ns_legendre\":{:.3},\"c_E\":{:.4},\"c_L\":{:.4},\"field_ops_add_affine\":\"1I+2M+1S\",\"field_ops_psi_x\":\"1I+10M (Horner 5M num, 4M den, 1M final)\"}},\n",
        xm46.num, xm46.den, cal46.0, cal46.1, cal46.2, cal46.3, cal46.4);

    // ---------- toys ----------
    json += "\"toys\":[\n";
    let mut first = true;
    let mut cut: Vec<String> = vec![];
    let mut toys_done = vec![];
    for &bits in &[16u32, 18, 20, 22] {
        if t_start.elapsed().as_secs_f64() > budget_s { cut.push(format!("phase1 size {bits}")); continue; }
        let mut rng = Rng(1000 + bits as u64);
        let toy = build_toy(bits, &mut rng);
        let (e, n, lam) = (toy.e, toy.n, toy.lam);
        let cal = calibrate(&e, &toy.xm, &mut rng, 10, 200000);
        let (c_e, c_l) = (cal.3, cal.4);
        println!("toy b={bits} t={} p={} N={} twist={} lam={} ord_lam={} c_E={:.3} c_L={:.3}", toy.t, e.p, n, toy.twist, lam, toy.ord_lam, c_e, c_l);
        let sq = (n as f64).sqrt();
        let f = sq.ceil() as u64;
        let kstar = ((n as f64) / (2.0 * c_e)).sqrt().ceil() as u64;
        let cap = 40 * n.max(1);
        let mut trng = Rng(20261015 + bits as u64);
        let mut g = e.rand_point(&mut trng); while g.is_none() { g = e.rand_point(&mut trng); }
        let ds: Vec<u64> = (0..targets).map(|_| 1 + trng.below(n - 1)).collect();
        let mut cc = Ctr::default();
        let qs: Vec<Pt> = ds.iter().map(|&d| e.mul(d, g, &mut cc)).collect();
        // prebuild chain tables (build cost charged per solve: table is rebuilt cost-wise for each target)
        let arm_names = ["rho", "rho-neg", "bsgs-pm", "abscissa-k1", "endo-chain-k4", "endo-chain-k16", "endo-chain-k64", "endo-chain-kF", "endo-chain-known", "endo-chain-k16-nofold"];
        let mut tabs: HashMap<&str, (HashMap<u64, (u32, u32)>, Ctr, u64)> = HashMap::new();
        for (name, k) in [("abscissa-k1", 1u64), ("endo-chain-k4", 4), ("endo-chain-k16", 16), ("endo-chain-k64", 64), ("endo-chain-kF", f)] {
            let mut tc = Ctr::default();
            let (t, dup) = chain_table(&e, &toy.xm, f, k, &mut tc);
            tabs.insert(name, (t, tc, dup));
        }
        let mut arm_json = vec![];
        let mut s_by_arm: HashMap<&str, (Vec<f64>, Vec<f64>)> = HashMap::new();
        let mut arm_rng = Rng(777 + bits as u64);
        for &arm in &arm_names {
            let (mut sv, mut sv1, mut av, mut ev, mut lv) = (vec![], vec![], vec![], vec![], vec![]);
            let mut stepv: Vec<f64> = vec![];
            let mut nok = 0;
            for (ti, &q) in qs.iter().enumerate() {
                let mut rr = Rng(arm_rng.next() ^ ti as u64);
                let out = match arm {
                    "rho" | "rho-neg" => { let (o, st) = rho(&e, n, g, q, &mut rr, arm == "rho-neg", cap); stepv.push(st as f64 / sq); o }
                    "bsgs-pm" => bsgs_pm(&e, n, g, q, cap),
                    "endo-chain-known" => known_chain_solve(&e, n, lam, &toy.xm, g, q, kstar, &mut rr, cap),
                    _ => {
                        let src = if arm == "endo-chain-k16-nofold" { "endo-chain-k16" } else { arm };
                        let (t, tc, _) = &tabs[src];
                        let mut c2 = *tc;
                        let ok = chain_solve(&e, n, lam, g, q, t, arm != "endo-chain-k16-nofold", &mut rr, &mut c2, cap);
                        Out { a: c2.a, e: c2.e, l: c2.l, ok }
                    }
                };
                if out.ok { nok += 1; }
                let gops = out.a as f64 + c_e * out.e as f64 + c_l * out.l as f64;
                let g1 = out.a as f64 + out.e as f64 + c_l * out.l as f64;
                sv.push(gops / sq); sv1.push(g1 / sq); av.push(out.a as f64); ev.push(out.e as f64); lv.push(out.l as f64);
            }
            let (m, lo, hi) = boot_ci(&sv, 7); let (m1, lo1, hi1) = boot_ci(&sv1, 7);
            let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
            println!("  {arm:24} ok={nok}/{targets} S={m:.3} [{lo:.3},{hi:.3}] S(c_E=1)={m1:.3} A={:.0} E={:.0} L={:.0}", mean(&av), mean(&ev), mean(&lv));
            let dup = tabs.get(arm).map(|x| x.2).unwrap_or(0);
            arm_json.push(format!("{{\"arm\":\"{arm}\",\"verified\":{nok},\"n\":{targets},\"S_mean\":{m:.5},\"S_ci95\":[{lo:.5},{hi:.5}],\"S_cE1_mean\":{m1:.5},\"S_cE1_ci95\":[{lo1:.5},{hi1:.5}],\"mean_A\":{:.1},\"mean_E\":{:.1},\"mean_L\":{:.1},\"table_duplicates\":{dup},\"walk_steps_over_sqrtN\":{},\"S\":{}}}", mean(&av), mean(&ev), mean(&lv), jarr(&stepv), jarr(&sv)));
            s_by_arm.insert(arm, (sv, sv1));
        }
        let mut ratios = vec![];
        for &arm in &arm_names {
            let (r, lo, hi) = boot_ratio(&s_by_arm[arm].0, &s_by_arm["rho-neg"].0, 11);
            ratios.push(format!("\"{arm}\":[{r:.4},{lo:.4},{hi:.4}]"));
        }
        let (rk, rklo, rkhi) = boot_ratio(&s_by_arm["endo-chain-k16"].0, &s_by_arm["abscissa-k1"].0, 13);
        let (rn, rnlo, rnhi) = boot_ratio(&s_by_arm["endo-chain-k16-nofold"].0, &s_by_arm["endo-chain-k16"].0, 17);
        println!("  ratio k16/abscissa = {rk:.3} [{rklo:.3},{rkhi:.3}]  nofold/k16 = {rn:.3} [{rnlo:.3},{rnhi:.3}]");
        if !first { json += ",\n"; } first = false;
        json += &format!("{{\"bits\":{bits},\"t\":{},\"p\":{},\"a\":{},\"b\":{},\"N\":{n},\"twist_of_j_model\":{},\"lambda\":{lam},\"ord_lambda\":{},\"index_of_lambda_subgroup\":{},\"xmap_num\":{:?},\"xmap_den\":{:?},\"ns_add\":{:.3},\"ns_psi_x\":{:.3},\"c_E\":{c_e:.4},\"c_L\":{c_l:.4},\"F\":{f},\"k_star\":{kstar},\"arms\":[{}],\"ratio_to_rho_neg\":{{{}}},\"ratio_k16_over_abscissa\":[{rk:.4},{rklo:.4},{rkhi:.4}],\"ratio_nofold_over_k16\":[{rn:.4},{rnlo:.4},{rnhi:.4}]}}",
            toy.t, e.p, e.a, e.b, toy.twist, toy.ord_lam, (n - 1) / toy.ord_lam, toy.xm.num, toy.xm.den, cal.0, cal.1, arm_json.join(",\n"), ratios.join(","));
        toys_done.push(toy);
    }
    json += "\n],\n";

    // ---------- phase 2: m=2 decomposition rate at 20 bits ----------
    let mut p2 = String::from("null");
    if t_start.elapsed().as_secs_f64() > budget_s { cut.push("phase2".into()); }
    else if let Some(toy) = toys_done.iter().find(|t| t.bits == 20) {
        let e = toy.e; let n = toy.n;
        let mut rng = Rng(2020);
        let f2 = ((n as f64).sqrt() / 4.0).ceil() as u64;
        let mut c = Ctr::default();
        let g = e.rand_point(&mut rng);
        let lift = |x: u64, e: &Curve| -> Pt { sqrtm(e.rhs(x), e.p).map(|y| (x, y)) };
        let mut sets: Vec<(&str, Vec<Pt>)> = vec![];
        let mut rs = vec![]; while rs.len() < f2 as usize { rs.push(e.rand_point(&mut rng)); } sets.push(("random", rs));
        let (t1, _) = chain_table(&e, &toy.xm, f2, 1, &mut c); sets.push(("small-x", t1.keys().map(|&x| lift(x, &e)).collect()));
        let (t16, _) = chain_table(&e, &toy.xm, f2, 16, &mut c); sets.push(("endo-chain-k16", t16.keys().map(|&x| lift(x, &e)).collect()));
        let (tf, _) = chain_table(&e, &toy.xm, f2, f2, &mut c); sets.push(("endo-chain-kF2", tf.keys().map(|&x| lift(x, &e)).collect()));
        let mut ap = vec![]; let mut r = g; for _ in 0..f2 { ap.push(r); r = e.add(r, g, &mut c); } sets.push(("arith-progression (positive control)", ap));
        let trials = 3000u64;
        let mut rows = vec![];
        let mut rrng = Rng(31337);
        let rpts: Vec<Pt> = (0..trials).map(|_| e.rand_point(&mut rrng)).collect();
        for (name, s) in &sets {
            let xs: std::collections::HashSet<u64> = s.iter().map(|q| q.unwrap().0).collect();
            let mut hit = 0u64;
            for &rp in &rpts {
                let mut found = false;
                for &pi in s {
                    for cand in [e.add(rp, e.neg(pi), &mut c), e.add(rp, pi, &mut c)] {
                        if let Some((x, _)) = cand { if xs.contains(&x) { found = true; break; } }
                    }
                    if found { break; }
                }
                if found { hit += 1; }
            }
            let (ph, lo, hi) = wilson(hit, trials);
            let distinct = xs.len();
            // heuristic for a random set: number of distinct +-Pi+-Pj ~ 2 F^2 (classes up to sign counted both)
            let pred = 1.0 - (-(2.0 * (distinct as f64).powi(2)) / n as f64).exp();
            println!("phase2 {name:38} F={distinct} rate={ph:.4} [{lo:.4},{hi:.4}] random-set heuristic={pred:.4}");
            rows.push(format!("{{\"fb\":\"{name}\",\"F\":{distinct},\"trials\":{trials},\"hits\":{hit},\"rate\":{ph:.5},\"wilson95\":[{lo:.5},{hi:.5}],\"random_set_heuristic\":{pred:.5}}}"));
        }
        p2 = format!("{{\"bits\":20,\"N\":{n},\"F2\":{f2},\"rows\":[{}]}}", rows.join(",\n"));
    } else { cut.push("phase2 (no 20-bit toy)".into()); }
    json += &format!("\"phase2_m2_rate\":{p2},\n\"cut\":{:?},\n\"wall_seconds\":{:.1}\n}}\n", cut, t_start.elapsed().as_secs_f64());
    std::fs::write("results.json", json).unwrap();
    println!("cut={:?} wall={:.1}s", cut, t_start.elapsed().as_secs_f64());
}
