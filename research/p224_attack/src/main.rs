//! Biehl-Meyer-Mueller invalid-curve attack on NIST P-224, in Rust.
//!
//! Mirrors the Python demo: scans invalid curves E_{b'}, finds small prime
//! factors q in #E_{b'}, queries a vulnerable ECDH oracle, solves DLPs via
//! BSGS, CRTs the residues.
//!
//! Point counting (#E_{b'}) is delegated to PARI/GP via subprocess — no
//! mature SEA implementation exists in pure Rust, and reimplementing
//! Schoof is out of scope.

use num_bigint::{BigInt, RandBigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::collections::HashMap;
use std::io::Write;
use std::process::{Command, Stdio};
use std::time::Instant;

// ---- Arithmetic helpers ----------------------------------------------------

fn modp(x: &BigInt, p: &BigInt) -> BigInt {
    let r = x % p;
    if r.sign() == Sign::Minus { r + p } else { r }
}

fn modinv(a: &BigInt, p: &BigInt) -> BigInt {
    let eg = modp(a, p).extended_gcd(p);
    assert!(eg.gcd.is_one(), "no inverse: gcd != 1");
    modp(&eg.x, p)
}

fn modpow(base: &BigInt, exp: &BigInt, m: &BigInt) -> BigInt {
    base.modpow(exp, m)
}

fn legendre(a: &BigInt, p: &BigInt) -> i32 {
    let r = modpow(&modp(a, p), &((p - 1) / 2), p);
    if r.is_zero() { 0 } else if r.is_one() { 1 } else { -1 }
}

/// Tonelli-Shanks modular square root.
fn sqrt_mod(n: &BigInt, p: &BigInt) -> Option<BigInt> {
    let n = modp(n, p);
    if n.is_zero() { return Some(BigInt::zero()); }
    if legendre(&n, p) != 1 { return None; }
    let mut q: BigInt = p - 1;
    let mut s: u32 = 0;
    while q.is_even() { q /= 2; s += 1; }
    if s == 1 {
        return Some(modpow(&n, &((p + BigInt::one()) / 4), p));
    }
    let mut z = BigInt::from(2);
    while legendre(&z, p) != -1 { z += 1; }
    let mut m = s;
    let mut c = modpow(&z, &q, p);
    let mut t = modpow(&n, &q, p);
    let mut r = modpow(&n, &((&q + 1) / 2), p);
    loop {
        if t.is_one() { return Some(r); }
        let mut i: u32 = 0;
        let mut tmp = t.clone();
        while !tmp.is_one() {
            tmp = (&tmp * &tmp) % p;
            i += 1;
            if i == m { return None; }
        }
        let mut b = c.clone();
        for _ in 0..(m - i - 1) { b = (&b * &b) % p; }
        m = i;
        c = (&b * &b) % p;
        t = (&t * &c) % p;
        r = (&r * &b) % p;
    }
}

// ---- Elliptic curve arithmetic --------------------------------------------
// Short-Weierstrass y^2 = x^3 + a*x + b. Formulas reference a but never b.

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Point(Option<(BigInt, BigInt)>);

impl Point {
    fn identity() -> Self { Point(None) }
    fn affine(x: BigInt, y: BigInt) -> Self { Point(Some((x, y))) }
    fn is_identity(&self) -> bool { self.0.is_none() }
    fn neg(&self, p: &BigInt) -> Self {
        match &self.0 {
            None => Point::identity(),
            Some((x, y)) => Point::affine(x.clone(), modp(&-y, p)),
        }
    }
}

fn ec_add(p1: &Point, p2: &Point, a: &BigInt, p: &BigInt) -> Point {
    let (px, qy) = match (&p1.0, &p2.0) {
        (None, _) => return p2.clone(),
        (_, None) => return p1.clone(),
        (Some(u), Some(v)) => (u, v),
    };
    let (x1, y1) = px;
    let (x2, y2) = qy;
    let m: BigInt = if x1 == x2 {
        if modp(&(y1 + y2), p).is_zero() {
            return Point::identity();
        }
        let num = modp(&(BigInt::from(3) * x1 * x1 + a), p);
        let den = modinv(&modp(&(BigInt::from(2) * y1), p), p);
        modp(&(num * den), p)
    } else {
        let num = modp(&(y2 - y1), p);
        let den = modinv(&modp(&(x2 - x1), p), p);
        modp(&(num * den), p)
    };
    let x3 = modp(&(&m * &m - x1 - x2), p);
    let y3 = modp(&(&m * (x1 - &x3) - y1), p);
    Point::affine(x3, y3)
}

fn ec_mul(k: &BigInt, base: &Point, a: &BigInt, p: &BigInt) -> Point {
    let mut r = Point::identity();
    let mut q = base.clone();
    let mut k = k.clone();
    while k.sign() == Sign::Plus {
        if k.is_odd() {
            r = ec_add(&r, &q, a, p);
        }
        q = ec_add(&q, &q, a, p);
        k >>= 1;
    }
    r
}

// ---- PARI/GP subprocess for #E(F_p) ---------------------------------------

fn ellcard(a: &BigInt, b: &BigInt, p: &BigInt) -> BigInt {
    let script = format!("print(ellcard(ellinit([{}, {}], {})));", a, b, p);
    let mut child = Command::new("gp")
        .args(["-q", "-s", "1000000000"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("failed to spawn gp");
    child.stdin.as_mut().unwrap().write_all(script.as_bytes()).unwrap();
    let out = child.wait_with_output().unwrap();
    let s = String::from_utf8(out.stdout).unwrap();
    let s: String = s.chars().filter(|c| !c.is_whitespace() && *c != '\\').collect();
    BigInt::parse_bytes(s.as_bytes(), 10)
        .unwrap_or_else(|| panic!("gp output not numeric: {s:?}"))
}

// ---- Trial factoring + ECM + BSGS + CRT -----------------------------------

/// Trial-divide n by primes up to `limit`. Returns the list of prime factors
/// found AND the remaining cofactor (which may still be composite).
fn trial_split(mut n: BigInt, limit: u64) -> (Vec<u64>, BigInt) {
    let mut out = Vec::new();
    let mut d: u64 = 2;
    while d <= limit && n > BigInt::one() {
        let dd = BigInt::from(d);
        if (&n % &dd).is_zero() {
            out.push(d);
            while (&n % &dd).is_zero() { n /= &dd; }
        }
        d = if d == 2 { 3 } else { d + 2 };
    }
    (out, n)
}

fn primes_up_to(n: u64) -> Vec<u64> {
    if n < 2 { return vec![]; }
    let n = n as usize;
    let mut sieve = vec![true; n + 1];
    sieve[0] = false; sieve[1] = false;
    for i in 2..=n {
        if sieve[i] {
            let mut j = i.saturating_mul(i);
            while j <= n { sieve[j] = false; j += i; }
        }
    }
    (2..=n).filter(|&i| sieve[i]).map(|i| i as u64).collect()
}

/// Miller-Rabin primality test. Deterministic for n < 3.3 · 10^24 with these
/// witnesses (covers everything we'd care about post-ECM split).
fn miller_rabin(n: &BigInt) -> bool {
    let one = BigInt::one();
    let two = BigInt::from(2);
    if n <= &one { return false; }
    if n == &two || n == &BigInt::from(3) { return true; }
    if n.is_even() { return false; }
    let n_minus_1: BigInt = n - 1;
    let mut d = n_minus_1.clone();
    let mut r: u32 = 0;
    while d.is_even() { d /= 2; r += 1; }
    for &w in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let wb = BigInt::from(w);
        if &wb >= n { continue; }
        let mut x = modpow(&wb, &d, n);
        if x.is_one() || x == n_minus_1 { continue; }
        let mut composite = true;
        for _ in 0..r.saturating_sub(1) {
            x = (&x * &x) % n;
            if x == n_minus_1 { composite = false; break; }
        }
        if composite { return false; }
    }
    true
}

/// Modular inverse that *reports* a factor of n if gcd(a, n) ∉ {1, n}.
/// This is the engine of ECM: when EC addition over Z/nZ needs to invert a
/// number that happens to share a factor with n, we recover the factor.
fn modinv_or_split(a: &BigInt, n: &BigInt) -> Result<BigInt, BigInt> {
    let am = modp(a, n);
    let eg = am.extended_gcd(n);
    if eg.gcd.is_one() {
        Ok(modp(&eg.x, n))
    } else if &eg.gcd != n {
        Err(eg.gcd.clone())
    } else {
        Err(n.clone())
    }
}

/// EC addition over Z/nZ — propagates the factor-discovery error.
fn ec_add_ecm(p1: &Point, p2: &Point, a: &BigInt, n: &BigInt) -> Result<Point, BigInt> {
    let (px, qy) = match (&p1.0, &p2.0) {
        (None, _) => return Ok(p2.clone()),
        (_, None) => return Ok(p1.clone()),
        (Some(u), Some(v)) => (u, v),
    };
    let (x1, y1) = px;
    let (x2, y2) = qy;
    let m: BigInt = if x1 == x2 {
        if modp(&(y1 + y2), n).is_zero() {
            return Ok(Point::identity());
        }
        let num = modp(&(BigInt::from(3) * x1 * x1 + a), n);
        let den = modp(&(BigInt::from(2) * y1), n);
        let den_inv = modinv_or_split(&den, n)?;
        modp(&(num * den_inv), n)
    } else {
        let num = modp(&(y2 - y1), n);
        let den = modp(&(x2 - x1), n);
        let den_inv = modinv_or_split(&den, n)?;
        modp(&(num * den_inv), n)
    };
    let x3 = modp(&(&m * &m - x1 - x2), n);
    let y3 = modp(&(&m * (x1 - &x3) - y1), n);
    Ok(Point::affine(x3, y3))
}

fn ec_mul_ecm(k: u64, base: &Point, a: &BigInt, n: &BigInt) -> Result<Point, BigInt> {
    let mut r = Point::identity();
    let mut q = base.clone();
    let mut k = k;
    while k > 0 {
        if k & 1 == 1 {
            r = ec_add_ecm(&r, &q, a, n)?;
        }
        q = ec_add_ecm(&q, &q, a, n)?;
        k >>= 1;
    }
    Ok(r)
}

enum Stage1 {
    Factor(BigInt),
    Continue { a: BigInt, q: Point },
    Dead,
}

/// Lenstra ECM stage 1. Pick a random curve E_a,b over Z/nZ through a
/// random point P=(x,y) (b is implicit, formulas don't reference it), then
/// compute k·P for k = ∏ p^⌊log_p B1⌋. If at any step a denominator shares
/// a factor with n, we return that factor; otherwise we hand the point off
/// to stage 2.
fn ecm_stage1(
    n: &BigInt,
    b1: u64,
    primes: &[u64],
    rng: &mut ChaCha8Rng,
) -> Stage1 {
    let a = rng.gen_bigint_range(&BigInt::one(), n);
    let x = rng.gen_bigint_range(&BigInt::one(), n);
    let y = rng.gen_bigint_range(&BigInt::one(), n);
    let mut p = Point::affine(x, y);
    for &prime in primes {
        if prime > b1 { break; }
        let mut pe: u64 = prime;
        while let Some(next) = pe.checked_mul(prime) {
            if next > b1 { break; }
            pe = next;
        }
        match ec_mul_ecm(pe, &p, &a, n) {
            Ok(np) => p = np,
            Err(g) => {
                if g > BigInt::one() && &g < n { return Stage1::Factor(g); }
                return Stage1::Dead;
            }
        }
    }
    Stage1::Continue { a, q: p }
}

fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 { let t = b; b = a % b; a = t; }
    a
}

/// Montgomery-style ECM stage 2. After stage 1's k·P, we want to catch the
/// case where the residual order on E mod p is *exactly one prime* q in
/// (B1, B2]. Use primorial D=210 to enumerate primes as q = i·D ± j with
/// j coprime to D, j ≤ D/2.
///
/// Baby steps: x(jQ) for j ∈ S = {j ≤ D/2 : gcd(j,D)=1}.
/// Giant steps: walk i·D·Q from ⌈B1/D⌉ to ⌊B2/D⌋ by adding D·Q each step.
/// Accumulate ∏ (x(iDQ) - x(jQ)) mod n. gcd with n recovers the factor.
fn ecm_stage2(
    n: &BigInt,
    a: &BigInt,
    q_point: &Point,
    b1: u64,
    b2: u64,
) -> Option<BigInt> {
    if q_point.is_identity() { return None; }
    let d: u64 = 210;
    let s: Vec<u64> = (1..=d/2).filter(|j| gcd_u64(*j, d) == 1).collect();

    // Baby steps: collect x-coordinates of j·Q.
    let mut baby_x: Vec<BigInt> = Vec::with_capacity(s.len());
    for &j in &s {
        match ec_mul_ecm(j, q_point, a, n) {
            Ok(pt) => {
                if let Some((x, _)) = pt.0 { baby_x.push(x); }
                else { return None; }
            }
            Err(g) => return if g > BigInt::one() && &g < n { Some(g) } else { None },
        }
    }

    // Pre-compute step D·Q and starting point ⌈B1/D⌉ · D · Q.
    let dq = match ec_mul_ecm(d, q_point, a, n) {
        Ok(p) => p,
        Err(g) => return if g > BigInt::one() && &g < n { Some(g) } else { None },
    };
    let i_start = (b1 + d - 1) / d;
    let i_end = b2 / d;
    let mut idq = match ec_mul_ecm(i_start * d, q_point, a, n) {
        Ok(p) => p,
        Err(g) => return if g > BigInt::one() && &g < n { Some(g) } else { None },
    };

    let mut acc = BigInt::one();
    let mut since_check: u32 = 0;
    for _ in i_start..=i_end {
        if let Some((x_idq, _)) = &idq.0 {
            for x_jq in &baby_x {
                let diff = modp(&(x_idq - x_jq), n);
                if !diff.is_zero() {
                    acc = (&acc * &diff) % n;
                }
            }
        }
        since_check += 1;
        if since_check >= 64 {
            since_check = 0;
            let g = acc.gcd(n);
            if g > BigInt::one() && &g < n { return Some(g); }
            if &g == n {
                // acc divisible by n: every factor was already in some prior step
                acc = BigInt::one();
            }
        }
        idq = match ec_add_ecm(&idq, &dq, a, n) {
            Ok(p) => p,
            Err(g) => return if g > BigInt::one() && &g < n { Some(g) } else { None },
        };
    }
    let g = acc.gcd(n);
    if g > BigInt::one() && &g < n { Some(g) } else { None }
}

#[derive(Debug, Clone, Copy)]
enum Stage { S1, S2 }

/// Combined ECM: stage 1, then stage 2 on the surviving point.
fn ecm(
    n: &BigInt,
    b1: u64,
    b2: u64,
    primes: &[u64],
    rng: &mut ChaCha8Rng,
) -> Option<(BigInt, Stage)> {
    match ecm_stage1(n, b1, primes, rng) {
        Stage1::Factor(g) => Some((g, Stage::S1)),
        Stage1::Continue { a, q } => ecm_stage2(n, &a, &q, b1, b2).map(|g| (g, Stage::S2)),
        Stage1::Dead => None,
    }
}

/// Find prime factors of n that are ≤ accept_bound. Combines trial division
/// (cheap up to ~10⁴) with ECM (finds factors up to ~2^40 with B1=50000).
/// Returns (factors_with_source). Source is "td" (trial-div), "s1", or "s2".
fn find_small_factors(
    n_orig: &BigInt,
    accept_bound: u64,
    b1: u64,
    b2: u64,
    primes: &[u64],
    rng: &mut ChaCha8Rng,
) -> Vec<(u64, &'static str)> {
    let mut out: Vec<(u64, &'static str)> = Vec::new();
    let (small, rest) = trial_split(n_orig.clone(), 10_000);
    for q in small { out.push((q, "td")); }
    let mut stack: Vec<(BigInt, &'static str)> = vec![(rest, "?")];
    let accept_big = BigInt::from(accept_bound);
    'outer: while let Some((m, src)) = stack.pop() {
        if m <= BigInt::one() { continue; }
        if miller_rabin(&m) {
            if m <= accept_big {
                if let Ok(v) = u64::try_from(m) { out.push((v, src)); }
            }
            continue;
        }
        for _ in 0..30 {
            if let Some((g, stage)) = ecm(&m, b1, b2, primes, rng) {
                if g > BigInt::one() && g < m {
                    let h = &m / &g;
                    let tag = match stage { Stage::S1 => "s1", Stage::S2 => "s2" };
                    stack.push((g, tag));
                    stack.push((h, tag));
                    continue 'outer;
                }
            }
        }
    }
    out.sort_by_key(|&(q, _)| q);
    out.dedup_by_key(|&mut (q, _)| q);
    out
}

fn bsgs(target: &Point, base: &Point, ord: u64, a: &BigInt, p: &BigInt) -> u64 {
    if target.is_identity() { return 0; }
    let m = ((ord as f64).sqrt() as u64) + 1;
    let mut baby: HashMap<Point, u64> = HashMap::new();
    let mut cur = Point::identity();
    baby.insert(cur.clone(), 0);
    for j in 1..m {
        cur = ec_add(&cur, base, a, p);
        baby.insert(cur.clone(), j);
    }
    let m_big = BigInt::from(m);
    let giant = ec_mul(&m_big, base, a, p).neg(p);
    let mut cur = target.clone();
    for i in 0..=m {
        if let Some(&j) = baby.get(&cur) {
            return (i * m + j) % ord;
        }
        cur = ec_add(&cur, &giant, a, p);
    }
    panic!("BSGS failed (ord={ord})");
}

fn crt(rems: &[BigInt], mods: &[BigInt]) -> (BigInt, BigInt) {
    let mut big_m = BigInt::one();
    for m in mods { big_m *= m; }
    let mut x = BigInt::zero();
    for (r, m) in rems.iter().zip(mods.iter()) {
        let mi = &big_m / m;
        let mi_inv = modinv(&mi, m);
        x += r * &mi * mi_inv;
    }
    (modp(&x, &big_m), big_m)
}

// ---- Main -----------------------------------------------------------------

fn main() {
    // NIST P-224
    let p = BigInt::parse_bytes(
        b"26959946667150639794667015087019630673557916260026308143510066298881", 10,
    ).unwrap();
    let a = modp(&BigInt::from(-3), &p);
    let b = BigInt::parse_bytes(
        b"B4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4", 16,
    ).unwrap();
    let n = BigInt::parse_bytes(
        b"FFFFFFFFFFFFFFFFFFFFFFFFFFFF16A2E0B8F03E13DD29455C5C2A3D", 16,
    ).unwrap();

    let mut rng = ChaCha8Rng::seed_from_u64(0xC0FFEE);
    let d = rng.gen_bigint_range(&BigInt::from(2), &n);
    println!("[*] P-224 server private key:");
    println!("    d = 0x{}", d.to_str_radix(16));

    let oracle = |q: &Point| ec_mul(&d, q, &a, &p);

    const TARGET_BITS: u64 = 80;
    println!("\n[*] Scanning invalid curves; target = recover {TARGET_BITS} bits of d.\n");

    let mut factors: HashMap<u64, (i64, Point)> = HashMap::new();
    let mut big_m = BigInt::one();
    let start = Instant::now();
    let mut bp: i64 = 0;
    let mut scanned = 0;

    // ECM machinery: precompute primes up to B1, share across calls.
    const ECM_B1: u64 = 500;            // intentionally tiny: forces stage 2 to do work
    const ECM_B2: u64 = 200_000;        // stage 2 ceiling — catches one extra prime in (B1, B2]
    const ACCEPT_BOUND: u64 = 1 << 40;
    let ecm_primes = primes_up_to(ECM_B1);

    let disc_nonzero = |bp_big: &BigInt| -> bool {
        let v = modp(
            &(BigInt::from(4) * &a * &a * &a + BigInt::from(27) * bp_big * bp_big),
            &p,
        );
        !v.is_zero()
    };

    while big_m.bits() < TARGET_BITS {
        bp += 1;
        let bp_big = BigInt::from(bp);
        if bp_big == b { continue; }
        if !disc_nonzero(&bp_big) { continue; }

        let n_prime = ellcard(&a, &bp_big, &p);
        scanned += 1;

        let small_qs = find_small_factors(&n_prime, ACCEPT_BOUND, ECM_B1, ECM_B2, &ecm_primes, &mut rng);
        let new_qs: Vec<(u64, &'static str)> = small_qs.iter().copied()
            .filter(|(q, _)| *q > 2 && !factors.contains_key(q))
            .collect();
        if new_qs.is_empty() { continue; }

        for (q, src) in new_qs {
            let q_big = BigInt::from(q);
            let cof = &n_prime / &q_big;
            let mut found: Option<Point> = None;
            for _ in 0..200 {
                let x = rng.gen_bigint_range(&BigInt::zero(), &p);
                let rhs = modp(&(&x * &x * &x + &a * &x + &bp_big), &p);
                let y = match sqrt_mod(&rhs, &p) { Some(y) => y, None => continue };
                let pt = Point::affine(x, y);
                let t = ec_mul(&cof, &pt, &a, &p);
                if t.is_identity() { continue; }
                if ec_mul(&q_big, &t, &a, &p).is_identity() {
                    found = Some(t);
                    break;
                }
            }
            let Some(pq) = found else { continue; };
            factors.insert(q, (bp, pq));
            big_m *= q;
            let qbits = 64 - q.leading_zeros();
            println!(
                "    b'={bp:3}  q={q:>10} ({qbits:>2}b, {src})  M = {:>3}b   [{:>5.1}s]",
                big_m.bits(),
                start.elapsed().as_secs_f64()
            );
            if big_m.bits() >= TARGET_BITS { break; }
        }
    }

    println!(
        "\n[*] {scanned} curves counted, {} subgroups, M = {}b  ({:.1}s)\n",
        factors.len(),
        big_m.bits(),
        start.elapsed().as_secs_f64()
    );

    println!("[*] Querying oracle and solving DLPs (BSGS) per subgroup...");
    let mut mods_vec = Vec::new();
    let mut rems_vec = Vec::new();
    let mut qs: Vec<u64> = factors.keys().copied().collect();
    qs.sort();
    for q in qs {
        let (bp_used, pq) = &factors[&q];
        let r = oracle(pq);
        let k = bsgs(&r, pq, q, &a, &p);
        rems_vec.push(BigInt::from(k));
        mods_vec.push(BigInt::from(q));
        println!("    via b'={bp_used:3}:  d mod {q:>6} = {k}");
    }

    let (d_rec, m_total) = crt(&rems_vec, &mods_vec);
    let d_mod_m = modp(&d, &m_total);
    println!(
        "\n[*] CRT result:    d ≡ 0x{}  (mod M, M is {}b)",
        d_rec.to_str_radix(16),
        m_total.bits()
    );
    println!("[*] True d mod M = 0x{}", d_mod_m.to_str_radix(16));
    println!("[*] MATCH: {}", d_rec == d_mod_m);
    println!("\n[*] Bits recovered: {} / {}", m_total.bits(), n.bits());
}
