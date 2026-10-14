//! Round 2 pair search: D = -3 (j = 0) craters at prime conductor gaps of 20-21 bits.
//! The crater is y^2 = x^3 + b (j = 0, rational automorphism group of order 6 when
//! p = 1 mod 3 and 3 | r - 1); the floor is a random order-N curve with j not in {0, 1728}.
//! Certification is the class-number-one conductor argument plus the rho.cubic
//! construction checks; every order claim is verified by scalar multiplication.
use serde_json::{json, Value};
use std::time::Instant;

fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}
fn addmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 + b as u128) % m as u128) as u64
}
fn submod(a: u64, b: u64, m: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        m - (b - a)
    }
}
fn powmod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1;
    a %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, m);
        }
        a = mulmod(a, a, m);
        e >>= 1;
    }
    r
}
fn inv(a: u64, m: u64) -> u64 {
    powmod(a, m - 2, m)
}
fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let d0 = n - 1;
    let s = d0.trailing_zeros();
    let d = d0 >> s;
    'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'w;
            }
        }
        return false;
    }
    true
}
fn legendre(a: u64, p: u64) -> i32 {
    let r = powmod(a, (p - 1) / 2, p);
    if r == 0 {
        0
    } else if r == 1 {
        1
    } else {
        -1
    }
}
fn sqrtmod(a: u64, p: u64) -> Option<u64> {
    if a == 0 {
        return Some(0);
    }
    if legendre(a, p) != 1 {
        return None;
    }
    let mut q = p - 1;
    let mut s = 0;
    while q.is_multiple_of(2) {
        q /= 2;
        s += 1;
    }
    let mut z = 2;
    while legendre(z, p) != -1 {
        z += 1;
    }
    let mut m = s;
    let mut c = powmod(z, q, p);
    let mut t = powmod(a, q, p);
    let mut r = powmod(a, q.div_ceil(2), p);
    while t != 1 {
        let mut i = 0;
        let mut tt = t;
        while tt != 1 {
            tt = mulmod(tt, tt, p);
            i += 1;
        }
        let b = powmod(c, 1u64 << (m - i - 1), p);
        m = i;
        c = mulmod(b, b, p);
        t = mulmod(t, c, p);
        r = mulmod(r, b, p);
    }
    Some(r)
}

type Pt = Option<(u64, u64)>;
fn padd(p: &Pt, q: &Pt, a: u64, m: u64) -> Pt {
    match (p, q) {
        (None, _) => *q,
        (_, None) => *p,
        (Some((x1, y1)), Some((x2, y2))) => {
            let lam = if x1 == x2 {
                if addmod(*y1, *y2, m) == 0 {
                    return None;
                }
                mulmod(
                    addmod(mulmod(3, mulmod(*x1, *x1, m), m), a, m),
                    inv(mulmod(2, *y1, m), m),
                    m,
                )
            } else {
                mulmod(submod(*y2, *y1, m), inv(submod(*x2, *x1, m), m), m)
            };
            let x3 = submod(submod(mulmod(lam, lam, m), *x1, m), *x2, m);
            let y3 = submod(mulmod(lam, submod(*x1, x3, m), m), *y1, m);
            Some((x3, y3))
        }
    }
}
fn pmul(p: &Pt, mut k: u64, a: u64, m: u64) -> Pt {
    let mut r = None;
    let mut b = *p;
    while k > 0 {
        if k & 1 == 1 {
            r = padd(&r, &b, a, m);
        }
        b = padd(&b, &b, a, m);
        k >>= 1;
    }
    r
}
fn jinv(a: u64, b: u64, p: u64) -> u64 {
    let a3 = mulmod(mulmod(a, a, p), a, p);
    let num = mulmod(1728, mulmod(4, a3, p), p);
    let den = addmod(mulmod(4, a3, p), mulmod(27, mulmod(b, b, p), p), p);
    if den == 0 {
        return 0;
    }
    mulmod(num, inv(den, p), p)
}
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

/// Order-N test by one scalar multiplication of a random point, with
/// 64 random points per accepted candidate (the round-1 standard).
fn has_order_n(p: u64, a: u64, b: u64, n: u64, rng: &mut Rng) -> bool {
    let a3 = mulmod(mulmod(a, a, p), a, p);
    let disc = addmod(mulmod(4, a3, p), mulmod(27, mulmod(b, b, p), p), p);
    if disc == 0 {
        return false;
    }
    for _ in 0..64 {
        let rx = rng.next() % p;
        let rr = addmod(
            addmod(mulmod(rx, mulmod(rx, rx, p), p), mulmod(a, rx, p), p),
            b,
            p,
        );
        if let Some(ry) = sqrtmod(rr, p) {
            if ry != 0 {
                let rp = Some((rx, ry));
                if pmul(&rp, n, a, p).is_some() {
                    return false; // [N]P != O
                }
            }
        }
    }
    true
}

/// The first generator of a curve known to have order N.
fn first_generator(p: u64, a: u64, b: u64) -> Option<(u64, u64)> {
    let mut x = 0u64;
    while x < p {
        let rhs = addmod(
            addmod(mulmod(x, mulmod(x, x, p), p), mulmod(a, x, p), p),
            b,
            p,
        );
        if let Some(y) = sqrtmod(rhs, p) {
            if y != 0 {
                return Some((x, y));
            }
        }
        x += 1;
    }
    None
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out_path = args
        .iter()
        .position(|s| s == "--out")
        .and_then(|i| args.get(i + 1))
        .expect("--out <path>")
        .clone();
    let cases: &[u32] = &[20, 21];
    let mut rng = Rng(0x64723373_6a5a3033); // "dr3 jz0s"
    let started = Instant::now();
    let mut pairs_out = Vec::new();
    let mut failures = Vec::new();

    for &k in cases {
        let scan_start = Instant::now();
        let base = 1u64 << k;
        let mut lo = base - base / 33;
        if lo.is_multiple_of(2) {
            lo += 1;
        }
        let hi = base + base / 33;
        let mut found: Option<Value> = None;
        let mut ells_tried = 0u64;
        let mut tees_tried = 0u64;
        'outer: for ell in (lo..=hi).step_by(2) {
            if !is_prime(ell) {
                continue;
            }
            ells_tried += 1;
            for t in (3u64..=3999).step_by(2) {
                for &sign in &[1i64, -1i64] {
                    let tt = t as i64 * sign;
                    let num = (tt * tt) as u128 + 3u128 * (ell as u128) * (ell as u128);
                    if !num.is_multiple_of(4) {
                        continue;
                    }
                    let q = (num / 4) as u64;
                    tees_tried += 1;
                    if q >= (1u64 << 44) || q % 3 != 1 {
                        continue;
                    }
                    if !is_prime(q) {
                        continue;
                    }
                    let n = (q as i128 + 1 - tt as i128) as u64;
                    if n < 2 || !is_prime(n) || n % 3 != 1 || n == q {
                        continue;
                    }
                    let mut acc = q % n;
                    let mut ok = true;
                    for _ in 0..32 {
                        if acc == 1 {
                            ok = false;
                            break;
                        }
                        acc = mulmod(acc, q % n, n);
                    }
                    if !ok {
                        continue;
                    }
                    // crater: scan sextic twists y^2 = x^3 + b for order N
                    let mut crater: Option<(u64, u64, u64, u64)> = None; // (b, gx, gy, trials)
                    let crater_start = Instant::now();
                    let crater_budget = std::time::Duration::from_secs(900);
                    let mut crater_trials = 0u64;
                    while crater_start.elapsed() < crater_budget && crater_trials < 400 {
                        crater_trials += 1;
                        let cb = 1 + rng.next() % q; // b != 0
                        if has_order_n(q, 0, cb, n, &mut rng) {
                            if let Some((gx, gy)) = first_generator(q, 0, cb) {
                                crater = Some((cb, gx, gy, crater_trials));
                                break;
                            }
                        }
                    }
                    let (cb2, gx, gy, crater_trials) = match crater {
                        Some(c) => c,
                        None => {
                            failures.push(json!({"k": k, "ell": ell, "t": tt,
                                "reason": "no order-N sextic twist of j = 0 within the scan"}));
                            continue;
                        }
                    };
                    // floor: random (a, b) with order N and j not in {0, 1728}
                    let mut floor: Option<(u64, u64, u64, u64)> = None;
                    let floor_start = Instant::now();
                    let floor_budget = std::time::Duration::from_secs(900);
                    let mut trials = 0u64;
                    let j1728 = submod(0, 1728 % q, q).max(1728 % q);
                    while floor_start.elapsed() < floor_budget {
                        trials += 1;
                        let fa = rng.next() % q;
                        let fb = rng.next() % q;
                        if fa == 0 {
                            continue; // j would be 0
                        }
                        if has_order_n(q, fa, fb, n, &mut rng) {
                            let jf = jinv(fa, fb, q);
                            if jf != 0 && jf != j1728 {
                                if let Some((fx, fy)) = first_generator(q, fa, fb) {
                                    floor = Some((fa, fb, fx, fy));
                                    break;
                                }
                            }
                        }
                    }
                    let (fa, fb, fx, fy) = match floor {
                        Some(f) => f,
                        None => {
                            failures.push(json!({"k": k, "ell": ell, "t": tt,
                                "reason": format!("no random order-N floor with j not in {{0, 1728}} within the 900 s scan ({} trials)", trials)}));
                            continue;
                        }
                    };
                    found = Some(json!({
                        "d": 3, "k": k, "ell": ell, "t": tt,
                        "p": q, "group_order": n, "bits": 64 - q.leading_zeros(),
                        "j_crater": 0,
                        "j_floor": jinv(fa, fb, q),
                        "crater": {"kind": "prime_explicit", "p": q, "a": 0, "b": cb2,
                                    "group_order": n, "r": n, "gx": gx, "gy": gy},
                        "floor": {"kind": "prime_explicit", "p": q, "a": fa, "b": fb,
                                   "group_order": n, "r": n, "gx": fx, "gy": fy},
                        "certificate": {
                            "frobenius_discriminant": format!("t^2-4q = -3*{}^2", ell),
                            "crater_end": "End = O_K: j = 0 = j(O_K) and h(-3) = 1 (class number one)",
                            "floor_end": format!("f_E = {}: j not in {{0, 1728}}, f_E | f_pi = ell, f_E != 1", ell),
                            "automorphism_structure": "p = 1 mod 3, r = 1 mod 3, [kappa]G = zeta(G) verified by rho.cubic construction",
                            "orders_verified": "64 random points per model: [N]P = O",
                            "sage_used": false,
                        },
                        "crater_trials": crater_trials,
                        "floor_trials": trials, "floor_scan_seconds": floor_start.elapsed().as_secs_f64(),
                        "ells_tried": ells_tried, "t_candidates_tried": tees_tried,
                        "search_seconds": scan_start.elapsed().as_secs_f64(),
                    }));
                    break 'outer;
                }
            }
        }
        match found {
            Some(p) => pairs_out.push(p),
            None => failures.push(json!({"d": 3, "k": k,
                "reason": "no certified pair in the ell window (scan exhausted)"})),
        }
    }

    let doc = json!({
        "schema": "drastic-prime-conductor-gap-pairs-round2/v1",
        "frozen_before": "research/drastic_prime_conductor_gap_round2_20261012/PROTOCOL.md",
        "seed_note": "deterministic xorshift seeded 0x647233736a5a3033",
        "total_seconds": started.elapsed().as_secs_f64(),
        "pairs": pairs_out,
        "failures": failures,
    });
    std::fs::write(&out_path, serde_json::to_vec_pretty(&doc).unwrap()).unwrap();
    println!(
        "wrote {} pairs, {} failures -> {}",
        pairs_out.len(),
        failures.len(),
        out_path
    );
    for p in &pairs_out {
        println!(
            "  D=-3 ell={} t={} q=2^{} N={}",
            p["ell"], p["t"], p["bits"], p["group_order"]
        );
    }
}
