// 10-bounded-degree: exact enumeration of bounded-degree algebraic factor bases
// over prime fields, with the extension-field positive control and a pair-46
// Monte Carlo arm. See PROTOCOL.md (frozen; sha256 in PROTOCOL.sha256).
use std::fmt::Write as _;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

static MULS: AtomicU64 = AtomicU64::new(0);

// ---------- scalar arithmetic ----------
#[inline]
fn mm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
fn pw(mut b: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64 % p;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mm(r, b, p);
        }
        b = mm(b, b, p);
        e >>= 1;
    }
    r
}
fn inv(a: u64, p: u64) -> u64 {
    pw(a, p - 2, p)
}
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for &q in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % q == 0 {
            return n == q;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'w: for &a in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = pw(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mm(x, x, n);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}
fn gcd_u(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd_u(b, a % b)
    }
}
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, p: u64) -> u64 {
        self.next() % p
    }
}
fn sm(x: u64) -> u64 {
    Rng(x).next()
}

// ---------- polynomials over F_p (low -> high) ----------
type Poly = Vec<u64>;
fn trim(mut a: Poly) -> Poly {
    while a.len() > 1 && *a.last().unwrap() == 0 {
        a.pop();
    }
    if a.is_empty() {
        a.push(0);
    }
    a
}
fn deg(a: &Poly) -> i64 {
    let a = trim(a.clone());
    if a.len() == 1 && a[0] == 0 {
        -1
    } else {
        a.len() as i64 - 1
    }
}
fn padd(a: &Poly, b: &Poly, p: u64) -> Poly {
    let n = a.len().max(b.len());
    trim((0..n)
        .map(|i| (a.get(i).copied().unwrap_or(0) + b.get(i).copied().unwrap_or(0)) % p)
        .collect())
}
fn psub(a: &Poly, b: &Poly, p: u64) -> Poly {
    let n = a.len().max(b.len());
    trim((0..n)
        .map(|i| (a.get(i).copied().unwrap_or(0) + p - b.get(i).copied().unwrap_or(0)) % p)
        .collect())
}
fn pmul(a: &Poly, b: &Poly, p: u64) -> Poly {
    let mut r = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = (r[i + j] + mm(x, y, p)) % p;
        }
    }
    MULS.fetch_add((a.len() * b.len()) as u64, Ordering::Relaxed);
    trim(r)
}
fn pscale(a: &Poly, c: u64, p: u64) -> Poly {
    trim(a.iter().map(|&x| mm(x, c, p)).collect())
}
fn pdivrem(a: &Poly, b: &Poly, p: u64) -> (Poly, Poly) {
    let a = trim(a.clone());
    let b = trim(b.clone());
    let db = deg(&b);
    assert!(db >= 0);
    let mut r = a.clone();
    let da = deg(&a);
    if da < db {
        return (vec![0], r);
    }
    let mut q = vec![0u64; (da - db + 1) as usize];
    let il = inv(b[db as usize], p);
    let mut cnt = 0u64;
    for k in (0..=(da - db) as usize).rev() {
        let c = mm(r[k + db as usize], il, p);
        q[k] = c;
        if c != 0 {
            for j in 0..=db as usize {
                r[k + j] = (r[k + j] + p - mm(c, b[j], p)) % p;
            }
            cnt += db as u64 + 2;
        }
    }
    MULS.fetch_add(cnt, Ordering::Relaxed);
    (trim(q), trim(r))
}
fn pmonic(a: &Poly, p: u64) -> Poly {
    let a = trim(a.clone());
    let l = *a.last().unwrap();
    if l == 0 {
        return a;
    }
    pscale(&a, inv(l, p), p)
}
fn pgcd(a: &Poly, b: &Poly, p: u64) -> Poly {
    let mut x = trim(a.clone());
    let mut y = trim(b.clone());
    while deg(&y) >= 0 {
        let (_, r) = pdivrem(&x, &y, p);
        x = y;
        y = r;
    }
    pmonic(&x, p)
}
fn pderiv(a: &Poly, p: u64) -> Poly {
    if a.len() <= 1 {
        return vec![0];
    }
    trim((1..a.len()).map(|i| mm(a[i], i as u64 % p, p)).collect())
}
fn peval(a: &Poly, x: u64, p: u64) -> u64 {
    let mut r = 0u64;
    for &c in a.iter().rev() {
        r = (mm(r, x, p) + c) % p;
    }
    r
}
/// Yun squarefree decomposition (char p > deg). Returns degree of the
/// odd-multiplicity part and whether a is squarefree.
fn odd_part_degree(a: &Poly, p: u64) -> (i64, bool) {
    let a = pmonic(a, p);
    let da = pderiv(&a, p);
    let a0 = pgcd(&a, &da, p);
    let sqfree = deg(&a0) == 0;
    let mut b = pdivrem(&a, &a0, p).0;
    let c = pdivrem(&da, &a0, p).0;
    let mut dd = psub(&c, &pderiv(&b, p), p);
    let mut i = 1i64;
    let mut odd = 0i64;
    while deg(&b) > 0 {
        let ai = pgcd(&b, &dd, p);
        if i % 2 == 1 {
            odd += deg(&ai);
        }
        let nb = pdivrem(&b, &ai, p).0;
        let nc = pdivrem(&dd, &ai, p).0;
        dd = psub(&nc, &pderiv(&nb, p), p);
        b = nb;
        i += 1;
        if i > 1000 {
            break;
        }
    }
    (odd, sqfree)
}
/// genus of y^2 = h (normalisation), None if h is const * square (reducible)
fn hyper_genus(h: &Poly, p: u64) -> (Option<i64>, bool) {
    let (r, sf) = odd_part_degree(h, p);
    if r == 0 {
        return (None, sf);
    }
    let br = r + (r % 2);
    (Some(br / 2 - 1), sf)
}

// ---------- elliptic curve (affine, for #E checks) ----------
#[derive(Clone, Copy, PartialEq)]
enum Pt {
    O,
    A(u64, u64),
}
fn ec_add(P: Pt, Q: Pt, a: u64, p: u64) -> Pt {
    match (P, Q) {
        (Pt::O, _) => Q,
        (_, Pt::O) => P,
        (Pt::A(x1, y1), Pt::A(x2, y2)) => {
            let l;
            if x1 == x2 {
                if (y1 + y2) % p == 0 {
                    return Pt::O;
                }
                l = mm((3 * mm(x1, x1, p) % p + a) % p, inv(2 * y1 % p, p), p);
            } else {
                l = mm((y2 + p - y1) % p, inv((x2 + p - x1) % p, p), p);
            }
            let x3 = (mm(l, l, p) + 2 * p - x1 - x2) % p;
            let y3 = (mm(l, (x1 + p - x3) % p, p) + p - y1) % p;
            Pt::A(x3, y3)
        }
    }
}
fn ec_mul(mut k: u64, P: Pt, a: u64, p: u64) -> Pt {
    let mut r = Pt::O;
    let mut b = P;
    while k > 0 {
        if k & 1 == 1 {
            r = ec_add(r, b, a, p);
        }
        b = ec_add(b, b, a, p);
        k >>= 1;
    }
    r
}

fn next_prime_3mod4(k: u32) -> u64 {
    let mut n = (1u64 << k) + 1;
    loop {
        if n % 4 == 3 && is_prime(n) {
            return n;
        }
        n += 1;
    }
}
fn ln(x: f64) -> f64 {
    x.ln()
}
fn cd(d: usize) -> f64 {
    let mut s = 0.0;
    let mut f = 1.0;
    for k in 1..=d {
        f *= k as f64;
        s += if k % 2 == 1 { 1.0 } else { -1.0 } / f;
    }
    s
}

struct Field {
    p: u64,
    a: u64,
    b: u64,
    chi: Vec<i8>,   // chi(g(x)) for x in F_p
    root: Vec<u32>, // a square root of g(x) when chi >= 0
    ne: u64,        // #E(F_p) including O
}

fn setup_field(k: u32, js: &mut String) -> Field {
    let p = next_prime_3mod4(k);
    let mut rng = Rng(1825564914 + 1000 * k as u64);
    let (a, b) = loop {
        let a = rng.below(p);
        let b = rng.below(p);
        let disc = (4 * mm(mm(a, a, p), a, p) + 27 * mm(b, b, p)) % p;
        if disc != 0 {
            break (a, b);
        }
    };
    // square table
    let mut sqrt_of: Vec<u32> = vec![u32::MAX; p as usize];
    for y in 0..p {
        let s = mm(y, y, p) as usize;
        if sqrt_of[s] == u32::MAX {
            sqrt_of[s] = y as u32;
        }
    }
    let mut chi = vec![0i8; p as usize];
    let mut root = vec![0u32; p as usize];
    let mut aff: u64 = 0;
    for x in 0..p {
        let gx = (mm(mm(x, x, p), x, p) + mm(a, x, p) + b) % p;
        let c: i8 = if gx == 0 {
            0
        } else if sqrt_of[gx as usize] != u32::MAX {
            1
        } else {
            -1
        };
        chi[x as usize] = c;
        if c >= 0 {
            root[x as usize] = sqrt_of[gx as usize];
        }
        aff += (1 + c as i64) as u64;
    }
    let ne = aff + 1;
    // checks: Hasse + [ne]R = O for 4 random points
    let hasse = ((ne as f64) - (p as f64) - 1.0).abs() <= 2.0 * (p as f64).sqrt();
    let mut okmul = true;
    let mut cnt = 0;
    while cnt < 4 {
        let x = rng.below(p);
        if chi[x as usize] == 1 {
            let P = Pt::A(x, root[x as usize] as u64);
            if ec_mul(ne, P, a, p) != Pt::O {
                okmul = false;
            }
            cnt += 1;
        }
    }
    let mr = is_prime(p);
    let _ = write!(
        js,
        "{{\"k\":{},\"p\":{},\"bits\":{},\"p_mod4\":{},\"prime_mr12\":{},\"a\":{},\"b\":{},\"nE\":{},\"hasse_ok\":{},\"order_check_4pts\":{}}}",
        k, p, 64 - p.leading_zeros(), p % 4, mr, a, b, ne, hasse, okmul
    );
    eprintln!("field k={} p={} a={} b={} #E={} hasse={} order={}", k, p, a, b, ne, hasse, okmul);
    Field { p, a, b, chi, root, ne }
}

/// count image and |F| from an image indicator
fn count_f(img: &[bool], f: &Field) -> (u64, u64) {
    let mut im = 0u64;
    let mut fb = 0u64;
    for x in 0..f.p as usize {
        if img[x] {
            im += 1;
            fb += (1 + f.chi[x] as i64) as u64;
        }
    }
    (im, fb)
}

fn rand_poly(d: usize, rng: &mut Rng, p: u64, monic: bool) -> Poly {
    let mut c: Poly = (0..=d).map(|_| rng.below(p)).collect();
    if monic {
        c[d] = 1;
    } else {
        while c[d] == 0 {
            c[d] = rng.below(p);
        }
    }
    c
}

fn row_json(
    js: &mut String,
    fam: &str,
    k: u32,
    p: u64,
    d: usize,
    seed: i64,
    im: u64,
    fb: u64,
    genus: Option<i64>,
    sqfree: Option<bool>,
    lb: Option<f64>,
    size_bits: Option<u64>,
    extra: &str,
) {
    let rho = fb as f64 / p as f64;
    let e = if fb > 0 { ln(fb as f64) / ln(p as f64) } else { 0.0 };
    let g = genus.map(|g| g.to_string()).unwrap_or("null".into());
    let sf = sqfree.map(|s| s.to_string()).unwrap_or("null".into());
    let l = lb.map(|v| format!("{:.3}", v)).unwrap_or("null".into());
    let lbok = lb.map(|v| ((fb as f64) >= v).to_string()).unwrap_or("null".into());
    let sb = size_bits.map(|v| v.to_string()).unwrap_or("null".into());
    let _ = write!(
        js,
        "{{\"family\":\"{}\",\"k\":{},\"p\":{},\"d\":{},\"seed\":{},\"image\":{},\"F\":{},\"rho\":{:.6},\"exponent\":{:.6},\"c_d\":{:.6},\"genus\":{},\"h_squarefree\":{},\"LB\":{},\"LB_ok\":{},\"size_bits\":{}{}}}",
        fam, k, p, d, seed, im, fb, rho, e, cd(d), g, sf, l, lbok, sb, extra
    );
}

fn lb_of(p: u64, g: Option<i64>, sf: bool, d: usize, degq: usize) -> Option<f64> {
    match g {
        Some(g) if sf => {
            let v = (p as f64 + 1.0 - 2.0 * g as f64 * (p as f64).sqrt() - 3.0 - 2.0 * degq as f64)
                / d as f64;
            Some(v.max(0.0))
        }
        _ => None,
    }
}

fn run_prime_field(k: u32, rows: &mut Vec<String>, fields: &mut Vec<String>, skip_rr: bool) {
    let mut fj = String::new();
    let f = setup_field(k, &mut fj);
    fields.push(fj);
    let p = f.p;
    let bl = 64 - p.leading_zeros() as u64;
    let gpoly: Poly = vec![f.b, f.a, 0, 1];
    let mut img = vec![false; p as usize];
    // negative control (i): identity
    {
        img.iter_mut().for_each(|v| *v = true);
        let (im, fb) = count_f(&img, &f);
        let mut s = String::new();
        let ok = fb == f.ne - 1;
        row_json(&mut s, "ID", k, p, 1, -1, im, fb, Some(1), Some(true), None, Some(2 * bl),
            &format!(",\"control_ok\":{}", ok));
        rows.push(s);
    }
    for d in 2..=8usize {
        let t0 = Instant::now();
        // RP
        for si in 0..3u64 {
            let seed = 1825564914 + si;
            let mut rng = Rng(seed ^ (k as u64) << 40 ^ (d as u64) << 48);
            let fp = rand_poly(d, &mut rng, p, true);
            img.iter_mut().for_each(|v| *v = false);
            for t in 0..p {
                let mut r = 0u64;
                for &c in fp.iter().rev() {
                    r = (r * t + c) % p;
                }
                img[r as usize] = true;
            }
            let (im, fb) = count_f(&img, &f);
            // h = g(f)
            let f2 = pmul(&fp, &fp, p);
            let h = padd(&padd(&pmul(&f2, &fp, p), &pscale(&fp, f.a, p), p), &vec![f.b], p);
            let (g, sf) = hyper_genus(&h, p);
            let mut s = String::new();
            row_json(&mut s, "RP", k, p, d, seed as i64, im, fb, g, Some(sf), lb_of(p, g, sf, d, 0),
                Some((d as u64 + 1) * bl), "");
            rows.push(s);
        }
        // RR
        if !(skip_rr) {
            for si in 0..3u64 {
                let seed = 1825564914 + si;
                let mut rng = Rng(seed ^ (k as u64) << 40 ^ (d as u64) << 48 ^ 0xABCDEF);
                let (P, Q) = loop {
                    let P = rand_poly(d, &mut rng, p, false);
                    let Q = rand_poly(d - 1, &mut rng, p, true);
                    if deg(&pgcd(&P, &Q, p)) == 0 {
                        break (P, Q);
                    }
                };
                // evaluate Q at all t, batch invert
                let qv: Vec<u64> = (0..p)
                    .map(|t| {
                        let mut r = 0u64;
                        for &c in Q.iter().rev() {
                            r = (r * t + c) % p;
                        }
                        r
                    })
                    .collect();
                let mut pref = vec![1u64; p as usize + 1];
                for i in 0..p as usize {
                    let q = if qv[i] == 0 { 1 } else { qv[i] };
                    pref[i + 1] = pref[i] * q % p;
                }
                let mut acc = inv(pref[p as usize], p);
                img.iter_mut().for_each(|v| *v = false);
                for i in (0..p as usize).rev() {
                    let q = if qv[i] == 0 { 1 } else { qv[i] };
                    let qi = acc * pref[i] % p;
                    acc = acc * q % p;
                    if qv[i] == 0 {
                        continue; // pole -> infinity (excluded)
                    }
                    let t = i as u64;
                    let mut r = 0u64;
                    for &c in P.iter().rev() {
                        r = (r * t + c) % p;
                    }
                    img[(r * qi % p) as usize] = true;
                }
                // f(infinity) = infinity since deg P > deg Q: excluded
                let (im, fb) = count_f(&img, &f);
                let Q2 = pmul(&Q, &Q, p);
                let P3 = pmul(&pmul(&P, &P, p), &P, p);
                let inner = padd(&padd(&P3, &pscale(&pmul(&P, &Q2, p), f.a, p), p),
                    &pscale(&pmul(&Q2, &Q, p), f.b, p), p);
                let h = pmul(&Q, &inner, p);
                let (g, sf) = hyper_genus(&h, p);
                let mut s = String::new();
                row_json(&mut s, "RR", k, p, d, seed as i64, im, fb, g, Some(sf),
                    lb_of(p, g, sf, d, d - 1), Some((2 * d as u64 + 1) * bl), "");
                rows.push(s);
            }
        }
        // PW
        {
            img.iter_mut().for_each(|v| *v = false);
            for t in 0..p {
                img[pw(t, d as u64, p) as usize] = true;
            }
            let (im, fb) = count_f(&img, &f);
            let expect = (p - 1) / gcd_u(d as u64, p - 1) + 1;
            let mut fp = vec![0u64; d + 1];
            fp[d] = 1;
            let f2 = pmul(&fp, &fp, p);
            let h = padd(&padd(&pmul(&f2, &fp, p), &pscale(&fp, f.a, p), p), &vec![f.b], p);
            let (g, sf) = hyper_genus(&h, p);
            let mut s = String::new();
            row_json(&mut s, "PW", k, p, d, -1, im, fb, g, Some(sf), lb_of(p, g, sf, d, 0),
                Some((d as u64 + 1) * bl),
                &format!(",\"image_expected\":{},\"control_ok\":{}", expect, expect == im));
            rows.push(s);
        }
        // DK: D_d(t,1)
        {
            let mut dk: Vec<Poly> = vec![vec![2], vec![0, 1]];
            for n in 2..=d {
                let a = pmul(&vec![0, 1], &dk[n - 1], p);
                dk.push(psub(&a, &dk[n - 2], p));
            }
            let fp = dk[d].clone();
            img.iter_mut().for_each(|v| *v = false);
            for t in 0..p {
                let (mut x0, mut x1) = (2u64, t);
                for _ in 1..d {
                    let x2 = (t * x1 % p + p - x0) % p;
                    x0 = x1;
                    x1 = x2;
                }
                img[x1 as usize] = true;
            }
            let (im, fb) = count_f(&img, &f);
            let perm = gcd_u(d as u64, p - 1) == 1 && gcd_u(d as u64, p + 1) == 1;
            let f2 = pmul(&fp, &fp, p);
            let h = padd(&padd(&pmul(&f2, &fp, p), &pscale(&fp, f.a, p), p), &vec![f.b], p);
            let (g, sf) = hyper_genus(&h, p);
            let ok = if perm { im == p } else { true };
            let mut s = String::new();
            row_json(&mut s, "DK", k, p, d, -1, im, fb, g, Some(sf), lb_of(p, g, sf, d, 0),
                Some((d as u64 + 1) * bl),
                &format!(",\"perm_predicted\":{},\"control_ok\":{}", perm, ok));
            rows.push(s);
        }
        // ZD
        for si in 0..3u64 {
            let seed = 1825564914 + si;
            let mut rng = Rng(seed ^ (k as u64) << 40 ^ (d as u64) << 48 ^ 0x5A5A);
            let A = rand_poly(d, &mut rng, p, false);
            let B = if d >= 2 { rand_poly(d - 2, &mut rng, p, false) } else { vec![0] };
            let mut fb = 0u64;
            for x in 0..p {
                let c = f.chi[x as usize];
                if c < 0 {
                    continue;
                }
                let y = f.root[x as usize] as u64;
                let ax = peval(&A, x, p);
                let bx = peval(&B, x, p);
                let byv = mm(bx, y, p);
                if (ax + byv) % p == 0 {
                    fb += 1;
                }
                if c == 1 && (ax + p - byv) % p == 0 {
                    fb += 1;
                }
            }
            let mut s = String::new();
            row_json(&mut s, "ZD", k, p, d, seed as i64, fb, fb, None, None, None,
                Some((2 * d as u64) * bl),
                &format!(",\"UB\":{},\"UB_ok\":{}", 2 * d, fb <= 2 * d as u64));
            rows.push(s);
        }
        // RF: random function (non-algebraic) -- independent of d, run once per d slot for d==2 only
        if d == 2 {
            for si in 0..3u64 {
                let seed = 1825564914 + si;
                img.iter_mut().for_each(|v| *v = false);
                for t in 0..p {
                    img[(sm(seed.wrapping_mul(0x1000193) ^ t ^ (k as u64) << 50) % p) as usize] = true;
                }
                let (im, fb) = count_f(&img, &f);
                let mut s = String::new();
                row_json(&mut s, "RF", k, p, 0, seed as i64, im, fb, None, None, None, None, "");
                rows.push(s);
            }
        }
        eprintln!("k={} d={} done in {:.1}s", k, d, t0.elapsed().as_secs_f64());
    }
    drop(f);
    let _ = (gpoly,);
}

// ---------- extension control over F_{q^2} = F_q[i]/(i^2+1) ----------
fn run_ext(k: u32, rows: &mut Vec<String>) {
    let q = next_prime_3mod4(k);
    let irred = pw(q - 1, (q - 1) / 2, q) == q - 1;
    let mut rng = Rng(1825564914 + 7777 * k as u64);
    let nz = |r: &mut Rng| loop {
        let v = r.below(q);
        if v != 0 {
            break v;
        }
    };
    let (ar, ai, br, bi) = (rng.below(q), nz(&mut rng), rng.below(q), nz(&mut rng));
    // legendre table over F_q
    let mut leg = vec![-1i8; q as usize];
    leg[0] = 0;
    for y in 1..q {
        leg[(y * y % q) as usize] = 1;
    }
    let mul = |a: (u64, u64), b: (u64, u64)| {
        ((a.0 * b.0 % q + q * q - a.1 * b.1 % q) % q, (a.0 * b.1 + a.1 * b.0) % q)
    };
    let gfun = |x: (u64, u64)| {
        let x2 = mul(x, x);
        let x3 = mul(x2, x);
        let ax = mul((ar, ai), x);
        ((x3.0 + ax.0 + br) % q, (x3.1 + ax.1 + bi) % q)
    };
    let chi2 = |z: (u64, u64)| -> i64 {
        let n = (z.0 * z.0 + z.1 * z.1) % q;
        leg[n as usize] as i64
    };
    // nonsingular check: disc = 4A^3 + 27B^2 != 0
    let a3 = mul(mul((ar, ai), (ar, ai)), (ar, ai));
    let b2 = mul((br, bi), (br, bi));
    let disc = ((4 * a3.0 + 27 * b2.0) % q, (4 * a3.1 + 27 * b2.1) % q);
    let nonsing = disc != (0, 0);
    let mut aff: i64 = 0;
    for u in 0..q {
        for v in 0..q {
            aff += 1 + chi2(gfun((u, v)));
        }
    }
    let ne = aff as u64 + 1;
    let qq = q * q;
    let hasse = ((ne as f64) - (qq as f64) - 1.0).abs() <= 2.0 * q as f64;
    let mut fb: i64 = 0;
    for u in 0..q {
        fb += 1 + chi2(gfun((u, 0)));
    }
    let e = (fb as f64).ln() / (qq as f64).ln();
    let valid = (0.40..=0.60).contains(&e) && irred && nonsing && hasse;
    let mut s = String::new();
    let _ = write!(
        s,
        "{{\"family\":\"EXT\",\"k\":{},\"q\":{},\"Q\":{},\"i2p1_irreducible\":{},\"A\":[{},{}],\"B\":[{},{}],\"nonsingular\":{},\"nE\":{},\"hasse_ok\":{},\"F\":{},\"F_over_q\":{:.6},\"rho\":{:.8},\"exponent\":{:.6},\"control_valid\":{}}}",
        k, q, qq, irred, ar, ai, br, bi, nonsing, ne, hasse, fb, fb as f64 / q as f64,
        fb as f64 / qq as f64, e, valid
    );
    eprintln!("ext k={} q={} F={} e={:.4} valid={}", k, q, fb, e, valid);
    rows.push(s);
}

// ---------- pair 46 Monte Carlo ----------
fn p46(rows: &mut Vec<String>, nsamp: u64) {
    let p: u64 = 5947698030651343;
    let a: u64 = 5757116383569871;
    let b: u64 = 3560173328417442;
    assert!(is_prime(p));
    for d in 2..=8usize {
        let t0 = Instant::now();
        let mut rng = Rng(1825564914);
        let mut rng = Rng(rng.next() ^ (d as u64) << 48 ^ 46);
        let fp = rand_poly(d, &mut rng, p, true);
        let mut srng = Rng(1825564915 + d as u64);
        let mut sum = 0f64;
        let mut sum2 = 0f64;
        let mut in_img = 0u64;
        let mut muls_member = 0u64;
        for _ in 0..nsamp {
            let x = srng.below(p);
            let gx = (mm(mm(x, x, p), x, p) + mm(a, x, p) + b) % p;
            let c: i64 = if gx == 0 { 0 } else if pw(gx, (p - 1) / 2, p) == 1 { 1 } else { -1 };
            // membership: gcd(t^p - t, f(t) - x) != 1
            let before = MULS.load(Ordering::Relaxed);
            let mut h = fp.clone();
            h[0] = (h[0] + p - x) % p;
            // t^p mod h by square-and-multiply
            let mut r: Poly = vec![1];
            let mut base: Poly = vec![0, 1];
            let mut e = p;
            while e > 0 {
                if e & 1 == 1 {
                    r = pdivrem(&pmul(&r, &base, p), &h, p).1;
                }
                base = pdivrem(&pmul(&base, &base, p), &h, p).1;
                e >>= 1;
            }
            let rt = psub(&r, &vec![0, 1], p);
            let gg = pgcd(&h, &rt, p);
            let member = deg(&gg) >= 1;
            muls_member += MULS.load(Ordering::Relaxed) - before;
            let v = if member { (1 + c) as f64 } else { 0.0 };
            if member {
                in_img += 1;
            }
            sum += v;
            sum2 += v * v;
        }
        let n = nsamp as f64;
        let mean = sum / n;
        let var = (sum2 / n - mean * mean) * n / (n - 1.0);
        let se = (var / n).sqrt();
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"family\":\"P46\",\"p\":{},\"d\":{},\"n\":{},\"rho_hat\":{:.6},\"ci95\":[{:.6},{:.6}],\"image_frac_hat\":{:.6},\"c_d\":{:.6},\"muls_per_membership\":{:.1},\"size_bits\":{},\"coeffs\":\"{:?}\"}}",
            p, d, nsamp, mean, mean - 1.96 * se, mean + 1.96 * se, in_img as f64 / n, cd(d),
            muls_member as f64 / n, (d as u64 + 1) * 53, fp
        );
        eprintln!("p46 d={} rho_hat={:.4} ({:.1}s)", d, mean, t0.elapsed().as_secs_f64());
        rows.push(s);
    }
}

fn main() {
    let t0 = Instant::now();
    let args: Vec<String> = std::env::args().collect();
    let nsamp: u64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(20000);
    let mut rows = Vec::new();
    let mut fields = Vec::new();
    for &k in &[16u32, 20, 24] {
        run_prime_field(k, &mut rows, &mut fields, false);
    }
    for &k in &[8u32, 10, 12] {
        run_ext(k, &mut rows);
    }
    p46(&mut rows, nsamp);
    let mut out = String::new();
    let _ = write!(out, "{{\"protocol_sha256\":\"{}\",\"wall_s\":{:.1},\"fields\":[{}],\"rows\":[\n{}\n]}}\n",
        std::fs::read_to_string("PROTOCOL.sha256").unwrap_or_default().split_whitespace().next().unwrap_or(""),
        t0.elapsed().as_secs_f64(), fields.join(","), rows.join(",\n"));
    std::fs::write("raw.json", out).unwrap();
    eprintln!("total {:.1}s", t0.elapsed().as_secs_f64());
}
