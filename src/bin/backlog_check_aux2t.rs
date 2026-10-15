//! `backlog_check_aux2t` — the independent check of a claimed auxiliary curve
//! for `backlog/XP-20261015-p192-isogeny-aux-search-t48`.
//!
//! ```text
//! backlog_check_aux2t ARTIFACT.json --expect-t N [--allow-test-prime]
//! ```
//!
//! The artifact is `{"p", "a", "b", "Qx", "Qy", "t"}`: integers as decimal
//! (or `0x` hex) strings, `t` a JSON integer. It is accepted when
//!
//! - `p` is the NIST P-192 prime `2^192 − 2^64 − 1` (with
//!   `--allow-test-prime`, any odd prime `p > 3`, for small-field tests);
//! - `0 ≤ a, b, Qx, Qy < p` and `4a³ + 27b² ≠ 0 (mod p)`;
//! - `Q = (Qx, Qy)` is on `y² = x³ + ax + b`;
//! - `t` equals `--expect-t` (the card pins 48: a point of order 2 proves nothing);
//! - `[2^(t−1)]Q ≠ O` and `[2^t]Q = O`, so `Q` has order exactly `2^t`.
//!
//! Deterministic, no randomness, no network. One JSON line on stdout; exit 0
//! accept, 1 reject, 2 bad input (unreadable or malformed artifact). It shares
//! no code with any solver: affine arithmetic over `num_bigint::BigUint`.

use std::process::ExitCode;

use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::Value;

/// `2^192 − 2^64 − 1`.
fn p192() -> BigUint {
    (BigUint::one() << 192u32) - (BigUint::one() << 64u32) - BigUint::one()
}

#[derive(Debug, PartialEq, Eq)]
enum Verdict {
    Accept(String),
    Reject(String),
    BadInput(String),
}

fn parse_int(v: &Value, name: &str) -> Result<BigUint, String> {
    let s = v
        .get(name)
        .and_then(Value::as_str)
        .ok_or_else(|| format!("{name} must be a string holding an integer"))?;
    // Strict: decimal digits, or `0x` and hex digits. No sign, no
    // underscores, no whitespace (`parse_bytes` alone would take `_`).
    let parsed = if let Some(hex) = s.strip_prefix("0x") {
        (!hex.is_empty() && hex.bytes().all(|b| b.is_ascii_hexdigit()))
            .then(|| BigUint::parse_bytes(hex.as_bytes(), 16))
            .flatten()
    } else if !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()) {
        BigUint::parse_bytes(s.as_bytes(), 10)
    } else {
        None
    };
    parsed.ok_or_else(|| format!("{name} = {s:?} is not a decimal or 0x-hex integer"))
}

/// Deterministic Miller–Rabin over the first twelve primes as bases: exact
/// below 3.3·10²⁴, a strong probable-prime test above. Only used for test
/// primes; the production prime is compared for equality.
fn is_prime(n: &BigUint) -> bool {
    let two = BigUint::from(2u32);
    if *n < two {
        return false;
    }
    let small = [2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for &q in &small {
        let q = BigUint::from(q);
        if *n == q {
            return true;
        }
        if (n % &q).is_zero() {
            return false;
        }
    }
    let n1 = n - BigUint::one();
    let mut d = n1.clone();
    let mut s = 0u32;
    while (&d % &two).is_zero() {
        d >>= 1u32;
        s += 1;
    }
    'bases: for &q in &small {
        let mut x = BigUint::from(q).modpow(&d, n);
        if x.is_one() || x == n1 {
            continue;
        }
        for _ in 1..s {
            x = (&x * &x) % n;
            if x == n1 {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

/// An affine point, `None` the point at infinity.
type Pt = Option<(BigUint, BigUint)>;

struct Curve {
    p: BigUint,
    a: BigUint,
}

impl Curve {
    fn sub(&self, x: &BigUint, y: &BigUint) -> BigUint {
        ((x + &self.p) - (y % &self.p)) % &self.p
    }
    fn inv(&self, x: &BigUint) -> BigUint {
        // p is prime: x^(p−2).
        x.modpow(&(&self.p - BigUint::from(2u32)), &self.p)
    }
    fn double(&self, pt: &Pt) -> Pt {
        let (x, y) = pt.as_ref()?;
        if y.is_zero() {
            return None;
        }
        let p = &self.p;
        let num = (BigUint::from(3u32) * x * x + &self.a) % p;
        let lambda = (num * self.inv(&((BigUint::from(2u32) * y) % p))) % p;
        let x3 = self.sub(&((&lambda * &lambda) % p), &((BigUint::from(2u32) * x) % p));
        let y3 = self.sub(&((&lambda * self.sub(x, &x3)) % p), y);
        Some((x3, y3))
    }
}

fn check(artifact: &Value, allow_test_prime: bool, expect_t: u64) -> Verdict {
    let ints = (|| -> Result<_, String> {
        Ok((
            parse_int(artifact, "p")?,
            parse_int(artifact, "a")?,
            parse_int(artifact, "b")?,
            parse_int(artifact, "Qx")?,
            parse_int(artifact, "Qy")?,
        ))
    })();
    let (p, a, b, x, y) = match ints {
        Ok(v) => v,
        Err(e) => return Verdict::BadInput(e),
    };
    let t = match artifact.get("t").and_then(Value::as_u64) {
        Some(t) if (1..=1024).contains(&t) => t,
        _ => return Verdict::BadInput("t must be an integer in 1..=1024".into()),
    };
    // The claim is about the card's t, not any t: a point of order 2 is not
    // a point of order 2^48.
    if t != expect_t {
        return Verdict::Reject(format!("t = {t}, but this check is for t = {expect_t}"));
    }
    if p != p192() {
        if !allow_test_prime {
            return Verdict::Reject("p is not the NIST P-192 prime 2^192 - 2^64 - 1".into());
        }
        if p <= BigUint::from(3u32) || !is_prime(&p) {
            return Verdict::Reject("p is not a prime > 3".into());
        }
    }
    for (name, v) in [("a", &a), ("b", &b), ("Qx", &x), ("Qy", &y)] {
        if *v >= p {
            return Verdict::Reject(format!("{name} is not reduced mod p"));
        }
    }
    let disc = (BigUint::from(4u32) * &a * &a * &a + BigUint::from(27u32) * &b * &b) % &p;
    if disc.is_zero() {
        return Verdict::Reject("singular curve: 4a^3 + 27b^2 = 0 mod p".into());
    }
    let lhs = (&y * &y) % &p;
    let rhs = (&x * &x * &x + &a * &x + &b) % &p;
    if lhs != rhs {
        return Verdict::Reject("Q is not on y^2 = x^3 + ax + b".into());
    }
    let curve = Curve { p, a };
    let mut q: Pt = Some((x, y));
    for i in 1..t {
        q = curve.double(&q);
        if q.is_none() {
            return Verdict::Reject(format!("[2^{i}]Q = O: the order of Q is below 2^{t}"));
        }
    }
    // q = [2^(t−1)]Q ≠ O.
    if curve.double(&q).is_some() {
        return Verdict::Reject(format!("[2^{t}]Q != O: the order of Q is not 2^{t}"));
    }
    Verdict::Accept(format!("Q has order exactly 2^{t}"))
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut allow = false;
    let mut expect_t: Option<u64> = None;
    let mut files: Vec<&String> = Vec::new();
    let mut bad = false;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--allow-test-prime" => allow = true,
            "--expect-t" => match it.next().map(|v| (v, v.parse::<u64>())) {
                Some((v, Ok(t)))
                    if v.bytes().all(|b| b.is_ascii_digit()) && (1..=1024).contains(&t) =>
                {
                    expect_t = Some(t)
                }
                _ => bad = true,
            },
            s if s.starts_with("--") => bad = true,
            _ => files.push(a),
        }
    }
    let verdict = if files.len() != 1 || bad || expect_t.is_none() {
        Verdict::BadInput(
            "usage: backlog_check_aux2t ARTIFACT.json --expect-t N [--allow-test-prime]".into(),
        )
    } else {
        let expect_t = expect_t.unwrap_or(0);
        match std::fs::read_to_string(files[0]) {
            Err(e) => Verdict::BadInput(format!("{}: {e}", files[0])),
            Ok(text) => match serde_json::from_str::<Value>(&text) {
                Err(e) => Verdict::BadInput(format!("not JSON: {e}")),
                Ok(v) if !v.is_object() => {
                    Verdict::BadInput("the artifact is not a JSON object".into())
                }
                Ok(v) => check(&v, allow, expect_t),
            },
        }
    };
    let (word, reason, code) = match &verdict {
        Verdict::Accept(r) => ("accept", r, 0),
        Verdict::Reject(r) => ("reject", r, 1),
        Verdict::BadInput(r) => ("bad_input", r, 2),
    };
    println!(
        "{}",
        serde_json::json!({"checker": "backlog_check_aux2t", "verdict": word, "reason": reason, "test_prime": allow, "expect_t": expect_t})
    );
    ExitCode::from(code)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// Order of an affine point by repeated *addition* in u64 arithmetic,
    /// independent of the doubling chain the checker uses.
    fn order_by_addition(p: u64, a: u64, x: u64, y: u64) -> u64 {
        let inv = |v: u64| -> u64 {
            let (mut r, mut b, mut e) = (1u64, v % p, p - 2);
            while e > 0 {
                if e & 1 == 1 {
                    r = r * b % p;
                }
                b = b * b % p;
                e >>= 1;
            }
            r
        };
        let add = |p1: Option<(u64, u64)>, p2: (u64, u64)| -> Option<(u64, u64)> {
            let (x1, y1) = p1?;
            let (x2, y2) = p2;
            let l = if x1 == x2 {
                if (y1 + y2).is_multiple_of(p) {
                    return None;
                }
                (3 * x1 % p * x1 % p + a) % p * inv(2 * y1 % p) % p
            } else {
                (y2 + p - y1) % p * inv((x2 + p - x1) % p) % p
            };
            let x3 = (l * l % p + 2 * p - x1 - x2) % p;
            let y3 = (l * ((x1 + p - x3) % p) % p + p - y1) % p;
            Some((x3, y3))
        };
        let mut acc = Some((x, y));
        let mut n = 1;
        while let Some(pt) = acc {
            acc = add(Some(pt), (x, y));
            n += 1;
            assert!(n < 4 * p, "runaway");
        }
        n
    }

    /// The first point (by a, b, x) over F_p of order exactly 2^t.
    fn find(p: u64, t: u32) -> (u64, u64, u64, u64) {
        let mut sqrt = vec![None; p as usize];
        for y in 0..p {
            sqrt[(y * y % p) as usize].get_or_insert(y);
        }
        for a in 0..p {
            for b in 1..p {
                if (4 * a * a % p * a + 27 * b % p * b).is_multiple_of(p) {
                    continue;
                }
                for x in 0..p {
                    let r = (x * x % p * x + a * x + b) % p;
                    if let Some(y) = sqrt[r as usize] {
                        if order_by_addition(p, a, x, y) == 1 << t {
                            return (a, b, x, y);
                        }
                    }
                }
            }
        }
        panic!("no point of order 2^{t} over F_{p}");
    }

    fn art(p: u64, a: u64, b: u64, x: u64, y: u64, t: u32) -> Value {
        json!({"p": p.to_string(), "a": a.to_string(), "b": b.to_string(), "Qx": x.to_string(), "Qy": y.to_string(), "t": t})
    }

    #[test]
    fn accepts_points_of_order_exactly_two_to_the_t_on_small_primes() {
        for (p, t) in [(103u64, 3u32), (241, 4), (521, 5), (1031, 6)] {
            let (a, b, x, y) = find(p, t);
            assert!(
                matches!(
                    check(&art(p, a, b, x, y, t), true, u64::from(t)),
                    Verdict::Accept(_)
                ),
                "p={p} t={t}"
            );
            // Wrong order claimed, both directions.
            assert!(
                matches!(check(&art(p, a, b, x, y, t + 1), true, u64::from(t + 1)), Verdict::Reject(r) if r.contains("below 2^"))
            );
            assert!(
                matches!(check(&art(p, a, b, x, y, t - 1), true, u64::from(t - 1)), Verdict::Reject(r) if r.contains("not 2^"))
            );
            // Off the curve.
            assert!(
                matches!(check(&art(p, a, b, x, (y + 1) % p, t), true, u64::from(t)), Verdict::Reject(r) if r.contains("not on"))
            );
            // Without --allow-test-prime a small prime is the wrong prime.
            assert!(
                matches!(check(&art(p, a, b, x, y, t), false, u64::from(t)), Verdict::Reject(r) if r.contains("P-192"))
            );
        }
    }

    #[test]
    fn a_point_of_smaller_order_is_rejected_early() {
        let (p, t) = (241u64, 4u32);
        let (a, b, x, y) = find(p, 2);
        let v = check(&art(p, a, b, x, y, t), true, u64::from(t));
        assert!(
            matches!(v, Verdict::Reject(ref r) if r.contains("[2^2]Q = O")),
            "{v:?}"
        );
    }

    #[test]
    fn singular_composite_unreduced_and_malformed_inputs() {
        // y^2 = x^3 over F_103 is singular; (1, 1) lies on it.
        assert!(
            matches!(check(&art(103, 0, 0, 1, 1, 3), true, 3), Verdict::Reject(r) if r.contains("singular"))
        );
        assert!(
            matches!(check(&art(105, 1, 1, 0, 1, 3), true, 3), Verdict::Reject(r) if r.contains("prime"))
        );
        assert!(
            matches!(check(&art(103, 103, 1, 0, 1, 3), true, 3), Verdict::Reject(r) if r.contains("reduced"))
        );
        assert!(matches!(
            check(&json!({"p": "103"}), true, 3),
            Verdict::BadInput(_)
        ));
        let mut bad = art(103, 1, 1, 0, 1, 3);
        bad["a"] = json!(1);
        assert!(matches!(check(&bad, true, 3), Verdict::BadInput(_)));
        bad = art(103, 1, 1, 0, 1, 0);
        assert!(matches!(check(&bad, true, 3), Verdict::BadInput(_)));
    }

    /// Review regression: (0, 0) on y² = x³ + x has order 2 over the P-192
    /// prime. It is a true point of order 2^1 and must not pass a t = 48 check.
    #[test]
    fn t_is_pinned_by_the_caller_not_the_artifact() {
        let a = json!({"p": p192().to_string(), "a": "1", "b": "0", "Qx": "0", "Qy": "0", "t": 1});
        assert!(matches!(check(&a, false, 48), Verdict::Reject(r) if r.contains("t = 1")));
        assert!(matches!(check(&a, false, 1), Verdict::Accept(_)));
    }

    #[test]
    fn integers_parse_strictly() {
        for bad in [
            "+5", "-5", "1_0", "0x", "0x_1", "0x1g", " 5", "5 ", "", "0X10",
        ] {
            let a = json!({"p": "103", "a": bad, "b": "1", "Qx": "0", "Qy": "1", "t": 3});
            assert!(
                matches!(check(&a, true, 3), Verdict::BadInput(_)),
                "{bad:?}"
            );
        }
        let ok = json!({"p": "0x67", "a": "0x0", "b": "1", "Qx": "0", "Qy": "1", "t": 3});
        assert!(!matches!(check(&ok, true, 3), Verdict::BadInput(_)));
    }

    #[test]
    fn the_p192_prime_is_the_standard_one_and_prime() {
        assert_eq!(
            p192().to_string(),
            "6277101735386680763835789423207666416083908700390324961279"
        );
        assert!(is_prime(&p192()));
        // A Carmichael number and a strong pseudoprime to bases 2, 3, 5, 7.
        assert!(!is_prime(&BigUint::from(561u32)));
        assert!(!is_prime(&BigUint::from(3_215_031_751u64)));
        assert!(is_prime(&BigUint::from(29_686_813_949_953u64)));
        // The P-192 generator has odd prime order, so it never satisfies t >= 1.
        let g = json!({
            "p": p192().to_string(), "a": (p192() - BigUint::from(3u32)).to_string(),
            "b": "0x64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1",
            "Qx": "0x188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012",
            "Qy": "0x07192b95ffc8da78631011ed6b24cdd573f977a11e794811", "t": 1
        });
        assert!(matches!(check(&g, false, 1), Verdict::Reject(r) if r.contains("not 2^1")));
    }
}
