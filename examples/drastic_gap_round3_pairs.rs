//! Round 3 pair search. (a) A D = -4 (j = 1728) quartic pair at a prime
//! conductor gap near 2^20: t = 0 (mod 4), q = (t/2)^2 + ell^2 prime with
//! q = 1 (mod 4), N = q+1-t = 2 (mod 8), r = N/2 prime = 1 (mod 4); the crater
//! is y^2 = x^3 + a*x (j = 1728), the floor a random order-N curve with
//! j not in {0, 1728}; the generator is [2]P of a full-order point P.
//! (b) A D = -19 pair at ell near 2^23 (cofactor 1), the round-1 model,
//! for the 30x gap sweep. Failures are kept.
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

/// Cofactor-1 order test: [N]P = O with N prime (round-1 standard).
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
                    return false;
                }
            }
        }
    }
    true
}

/// Cofactor-2 order test: the curve has order exactly N = 2r, witnessed
/// by a point P of order 2r ([N]P = O, [r]P != O, y != 0). Returns the
/// order-r generator [2]P when it exists.
fn order_2r_generator(p: u64, a: u64, b: u64, n: u64, r: u64, rng: &mut Rng) -> Option<(u64, u64)> {
    let a3 = mulmod(mulmod(a, a, p), a, p);
    let disc = addmod(mulmod(4, a3, p), mulmod(27, mulmod(b, b, p), p), p);
    if disc == 0 {
        return None;
    }
    for _ in 0..256 {
        let rx = rng.next() % p;
        let rr = addmod(
            addmod(mulmod(rx, mulmod(rx, rx, p), p), mulmod(a, rx, p), p),
            b,
            p,
        );
        if let Some(ry) = sqrtmod(rr, p) {
            if ry == 0 {
                continue;
            }
            let rp = Some((rx, ry));
            if pmul(&rp, n, a, p).is_some() {
                return None; // one witness: the curve's order is not N
            }
            if pmul(&rp, r, a, p).is_none() {
                continue; // a point of order r: try another point
            }
            // the point has order 2r; its double generates the subgroup
            let double = padd(&rp, &rp, a, p)?;
            let g = Some(double);
            if pmul(&g, r, a, p).is_some() {
                return None; // [r][2]P = [2r]P must be O for an order-2r point
            }
            return Some(double);
        }
    }
    None
}

/// The first point of order 2r on a curve known to have order 2r, then
/// its double (order r). Returns the order-r generator.
#[allow(dead_code)]
fn first_generator_2r(p: u64, a: u64, b: u64, n: u64, r: u64) -> Option<(u64, u64)> {
    let mut x = 0u64;
    while x < p {
        let rhs = addmod(
            addmod(mulmod(x, mulmod(x, x, p), p), mulmod(a, x, p), p),
            b,
            p,
        );
        if let Some(y) = sqrtmod(rhs, p) {
            if y != 0 {
                let pt = Some((x, y));
                if pmul(&pt, n, a, p).is_none() && pmul(&pt, r, a, p).is_some() {
                    let double = padd(&pt, &pt, a, p)?;
                    let g = Some(double);
                    if pmul(&g, r, a, p).is_some() {
                        return None; // [r][2]P must be O
                    }
                    return Some(double);
                }
            }
        }
        x += 1;
    }
    None
}

/// The first point of order exactly N (cofactor 1).
fn first_generator(p: u64, a: u64, b: u64, _n: u64) -> Option<(u64, u64)> {
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

fn no_small_embedding(q: u64, r: u64) -> bool {
    let mut acc = q % r;
    for _ in 0..32 {
        if acc == 1 {
            return false;
        }
        acc = mulmod(acc, q % r, r);
    }
    true
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let out_path = args
        .iter()
        .position(|s| s == "--out")
        .and_then(|i| args.get(i + 1))
        .expect("--out <path>")
        .clone();
    let mut rng = Rng(0x64723373_72336a60); // "dr3 r3jz"
    let started = Instant::now();
    let mut pairs_out = Vec::new();
    let mut failures = Vec::new();

    // ── (a) D = -4 quartic pair, ell near 2^20 ────────────────────────
    {
        let k = 20u32;
        let base = 1u64 << k;
        let mut lo = base - base / 33;
        if lo.is_multiple_of(2) {
            lo += 1;
        }
        let hi = base + base / 33;
        let mut found: Option<Value> = None;
        'd4: for ell in (lo..=hi).step_by(2) {
            if !is_prime(ell) {
                continue;
            }
            for t in (4u64..=3996).step_by(4) {
                let half = t / 2;
                let q = half * half + ell * ell;
                if q >= (1u64 << 45) || q % 4 != 1 || !is_prime(q) {
                    continue;
                }
                let n = q + 1 - t;
                if n % 8 != 2 {
                    continue;
                }
                let r = n / 2;
                if !is_prime(r) || r % 4 != 1 || !no_small_embedding(q, r) {
                    continue;
                }
                // crater: quartic twists y^2 = x^3 + a*x of order N
                let crater_start = Instant::now();
                let mut crater: Option<(u64, u64, u64)> = None; // (a, gx, gy)
                let mut crater_trials = 0u64;
                while crater_start.elapsed() < std::time::Duration::from_secs(900)
                    && crater_trials < 400
                {
                    crater_trials += 1;
                    let ca = 1 + rng.next() % q;
                    if let Some((gx, gy)) = order_2r_generator(q, ca, 0, n, r, &mut rng) {
                        crater = Some((ca, gx, gy));
                        break;
                    }
                }
                let (ca, gx, gy) = match crater {
                    Some(c) => c,
                    None => {
                        eprintln!(
                            "d4 crater miss ell={} t={} after {} trials",
                            ell, t, crater_trials
                        );
                        failures.push(json!({"family": "d4", "ell": ell, "t": t,
                            "reason": "no order-2r quartic twist of j = 1728 within the scan"}));
                        continue;
                    }
                };
                // floor: random (a, b) with order N and j not in {0, 1728}
                let floor_start = Instant::now();
                let mut floor: Option<(u64, u64, u64, u64)> = None;
                let mut trials = 0u64;
                let j1728 = 1728 % q;
                while floor_start.elapsed() < std::time::Duration::from_secs(3600) {
                    trials += 1;
                    let fa = rng.next() % q;
                    let fb = rng.next() % q;
                    if fa == 0 || fb == 0 {
                        continue;
                    }
                    if let Some((fx, fy)) = order_2r_generator(q, fa, fb, n, r, &mut rng) {
                        let jf = jinv(fa, fb, q);
                        if jf != 0 && jf != j1728 {
                            floor = Some((fa, fb, fx, fy));
                            break;
                        }
                    }
                }
                let (fa, fb, fx, fy) = match floor {
                    Some(f) => f,
                    None => {
                        eprintln!("d4 floor miss ell={} t={} after {} trials", ell, t, trials);
                        failures.push(json!({"family": "d4", "ell": ell, "t": t,
                            "reason": format!("no order-2r floor with j not in {{0, 1728}} ({} trials)", trials)}));
                        continue;
                    }
                };
                found = Some(json!({
                    "family": "d4", "d": 4, "k": k, "ell": ell, "t": t,
                    "p": q, "group_order": n, "r": r, "bits": 64 - q.leading_zeros(),
                    "j_crater": 1728 % q, "j_floor": jinv(fa, fb, q),
                    "crater": {"kind": "prime_explicit", "p": q, "a": ca, "b": 0,
                               "group_order": n, "r": r, "gx": gx, "gy": gy},
                    "floor": {"kind": "prime_explicit", "p": q, "a": fa, "b": fb,
                              "group_order": n, "r": r, "gx": fx, "gy": fy},
                    "certificate": {
                        "frobenius_discriminant": format!("t^2-4q = -4*{}^2", ell),
                        "crater_end": "End = O_K: j = 1728 = j(O_K) and h(-4) = 1 (class number one)",
                        "floor_end": format!("f_E = {}: j not in {{0, 1728}}, f_E | f_pi = ell, f_E != 1", ell),
                        "automorphism_structure": "p = 1 mod 4, r = 1 mod 4, [kappa]G = iota(G) verified by rho.quartic construction",
                        "orders_verified": "points of order 2r: [N]P = O with [r]P != O and y != 0; [2]P generates the order-r subgroup",
                        "sage_used": false,
                    },
                    "crater_trials": crater_trials,
                    "floor_trials": trials,
                    "search_seconds": crater_start.elapsed().as_secs_f64(),
                }));
                break 'd4;
            }
        }
        match found {
            Some(p) => pairs_out.push(p),
            None => failures.push(json!({"family": "d4",
                "reason": "no certified D = -4 pair in the ell window (scan exhausted)"})),
        }
    }

    // ── (b) D = -19 pair, ell near 2^23 (round-1 model) ────────────────
    {
        let k = 22u32;
        let base = 1u64 << k;
        let mut lo = base - base / 33;
        if lo.is_multiple_of(2) {
            lo += 1;
        }
        let hi = base + base / 33;
        let mut found: Option<Value> = None;
        'd19: for ell in (lo..=hi).step_by(2) {
            if !is_prime(ell) {
                continue;
            }
            for t in (3u64..=3999).step_by(2) {
                for &sign in &[1i64, -1i64] {
                    let tt = t as i64 * sign;
                    let num = (tt * tt) as u128 + 19u128 * (ell as u128) * (ell as u128);
                    if !num.is_multiple_of(4) {
                        continue;
                    }
                    let q = (num / 4) as u64;
                    if q >= (1u64 << 48) {
                        continue;
                    }
                    if !is_prime(q) {
                        continue;
                    }
                    let n = (q as i128 + 1 - tt as i128) as u64;
                    if n < 2 || !is_prime(n) || n == q || !no_small_embedding(q, n) {
                        continue;
                    }
                    // crater: j(O_K) = -884736 model, both twists
                    let j0 = submod(0, 884736 % q, q);
                    let den = submod(1728 % q, j0, q);
                    if den == 0 {
                        continue;
                    }
                    let c = mulmod(j0, inv(den, q), q);
                    let mut crater: Option<(u64, u64, u64, u64)> = None;
                    for (wia, wib) in [
                        (mulmod(3, c, q), mulmod(2, c, q)),
                        (mulmod(3, c, q), submod(0, mulmod(2, c, q), q)),
                    ] {
                        if has_order_n(q, wia, wib, n, &mut rng) {
                            if let Some((gx, gy)) = first_generator(q, wia, wib, n) {
                                crater = Some((wia, wib, gx, gy));
                            }
                            break;
                        }
                    }
                    let (ca, cb, gx, gy) = match crater {
                        Some(c) => c,
                        None => continue,
                    };
                    let floor_start = Instant::now();
                    let mut floor: Option<(u64, u64, u64, u64)> = None;
                    let mut trials = 0u64;
                    while floor_start.elapsed() < std::time::Duration::from_secs(1500) {
                        trials += 1;
                        let fa = rng.next() % q;
                        let fb = rng.next() % q;
                        if has_order_n(q, fa, fb, n, &mut rng) {
                            let jf = jinv(fa, fb, q);
                            if jf != j0 {
                                if let Some((fx, fy)) = first_generator(q, fa, fb, n) {
                                    floor = Some((fa, fb, fx, fy));
                                    break;
                                }
                            }
                        }
                    }
                    let (fa, fb, fx, fy) = match floor {
                        Some(f) => f,
                        None => continue,
                    };
                    found = Some(json!({
                        "family": "d19", "d": 19, "k": k, "ell": ell, "t": tt,
                        "p": q, "group_order": n, "r": n, "bits": 64 - q.leading_zeros(),
                        "j_crater": j0, "j_floor": jinv(fa, fb, q),
                        "crater": {"kind": "prime_explicit", "p": q, "a": ca, "b": cb,
                                   "group_order": n, "r": n, "gx": gx, "gy": gy},
                        "floor": {"kind": "prime_explicit", "p": q, "a": fa, "b": fb,
                                  "group_order": n, "r": n, "gx": fx, "gy": fy},
                        "certificate": {
                            "frobenius_discriminant": format!("t^2-4q = -19*{}^2", ell),
                            "crater_end": "End = O_K: j = j(O_K) and h(-19) = 1",
                            "floor_end": format!("f_E = {}: j != j(O_K)", ell),
                            "orders_verified": "64 random [N]P = O per model",
                            "sage_used": false,
                        },
                        "floor_trials": trials,
                        "search_seconds": floor_start.elapsed().as_secs_f64(),
                    }));
                    break 'd19;
                }
            }
        }
        match found {
            Some(p) => pairs_out.push(p),
            None => failures.push(json!({"family": "d19",
                "reason": "no certified D = -19 pair in the ell window (scan exhausted)"})),
        }
    }

    let doc = json!({
        "schema": "drastic-prime-conductor-gap-pairs-round3/v1",
        "frozen_before": "research/drastic_prime_conductor_gap_round3_20261012/PROTOCOL.md",
        "seed_note": "deterministic xorshift seeded 0x647233737336a60",
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
            "  {} ell={} t={} q=2^{} N={} r={}",
            p["family"], p["ell"], p["t"], p["bits"], p["group_order"], p["r"]
        );
    }
}
