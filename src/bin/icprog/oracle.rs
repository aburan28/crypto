//! An independent checker for binary Koblitz curves in a polynomial basis:
//! the tournament's `oracle.py`, natively, with the same checks and the same
//! refusals.
//!
//! Its arithmetic is its own: a field element is a `u64` and a product is
//! shift-and-add, reduced bit by bit, as `oracle.py` did with Python
//! integers.  Nothing here calls `crypto_lib`'s fields or curves, so a
//! certificate `ic` writes is checked by code `ic` does not share.

use super::json::J;

pub type Point = Option<(u64, u64)>;

fn require(cond: bool, message: &str) -> Result<(), String> {
    if cond {
        Ok(())
    } else {
        Err(message.to_string())
    }
}

/// The bit length of `a`, as Python's `int.bit_length`.
fn bits(a: u128) -> u32 {
    128 - a.leading_zeros()
}

/// GF(2)[x] division's remainder in the integer encoding.
fn poly_rem(mut a: u128, b: u128) -> u128 {
    while bits(a) >= bits(b) {
        a ^= b << (bits(a) - bits(b));
    }
    a
}

/// HAC Algorithm 4.69 for p = 2: no irreducible factor of degree at most
/// `n / 2`, by `gcd(f, x^(2^i) − x)`.
fn irreducible(modulus: u128) -> bool {
    let degree = bits(modulus) as i64 - 1;
    if degree < 1 {
        return false;
    }
    let mut power: u128 = 2;
    for _ in 0..degree / 2 {
        // `power` squared in GF(2)[x]: each bit moves to twice its place.
        let mut square = 0u128;
        for i in 0..bits(power) {
            if power >> i & 1 == 1 {
                square |= 1u128 << (2 * i);
            }
        }
        power = poly_rem(square, modulus);
        let (mut a, mut b) = (modulus, power ^ 2);
        while b != 0 {
            let r = poly_rem(a, b);
            a = b;
            b = r;
        }
        if a != 1 {
            return false;
        }
    }
    true
}

fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

fn pow_mod(mut a: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64 % m;
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

/// Whether `n` is prime: Miller–Rabin with the first twelve primes as
/// bases, which is exact for every 64-bit integer.  (`oracle.py` divided by
/// every integer up to the square root; the answer is the same.)
fn prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for p in BASES {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let (mut d, mut s) = (n - 1, 0);
    while d.is_multiple_of(2) {
        d /= 2;
        s += 1;
    }
    'bases: for a in BASES {
        let mut x = pow_mod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mul_mod(x, x, n);
            if x == n - 1 {
                continue 'bases;
            }
        }
        return false;
    }
    true
}

/// A fixture's integer field, which a report may write as a number or a
/// decimal string.
fn int(v: &J, name: &str) -> Result<i128, String> {
    match v {
        J::Int(i) => Ok(*i),
        J::Str(s) => s
            .trim()
            .parse()
            .map_err(|_| format!("`{name}` is not an integer")),
        _ => Err(format!("`{name}` is not an integer")),
    }
}

fn field<'a>(fixture: &'a J, key: &str) -> Result<&'a J, String> {
    fixture
        .get(key)
        .ok_or_else(|| format!("fixture has no `{key}`"))
}

#[derive(Clone, Debug)]
pub struct Curve {
    pub n: u32,
    pub a: u64,
    pub r: u64,
    pub h: u64,
    pub modulus: u64,
    pub g: (u64, u64),
    pub lam: u64,
}

impl Curve {
    /// The curve a fixture names, after every check `oracle.Curve` makes.
    pub fn new(fixture: &J) -> Result<Curve, String> {
        let n = int(field(fixture, "degree")?, "degree")?;
        let a = int(field(fixture, "curve_a")?, "curve_a")?;
        let r = int(field(fixture, "subgroup_order")?, "subgroup_order")?;
        require((5..=61).contains(&n) && n % 2 == 1, "unsupported degree")?;
        require(a == 0 || a == 1, "unsupported coefficient")?;
        let irr = field(fixture, "irreducible")?;
        // The exponents are JSON integers, and the degree is compared as
        // `oracle.py` compares it, by value: "41" is not 41.
        let terms: Vec<i128> = irr
            .get("low_terms")
            .and_then(J::as_arr)
            .ok_or("bad modulus")?
            .iter()
            .map(|t| t.as_i128().ok_or_else(|| "bad modulus".to_string()))
            .collect::<Result<_, _>>()?;
        require(
            irr.get("degree")
                .is_some_and(|d| super::json::py_eq(d, &J::Int(n))),
            "field degree mismatch",
        )?;
        let mut sorted = terms.clone();
        sorted.sort_unstable();
        sorted.dedup();
        require(
            sorted.len() == terms.len() && terms.iter().all(|&t| 0 <= t && t < n),
            "bad modulus",
        )?;
        let n = n as u32;
        let modulus = terms.iter().fold(1u64 << n, |m, &t| m | 1u64 << t);
        require(irreducible(modulus as u128), "reducible field modulus")?;
        require(r > 2 && r < 1 << 64 && prime(r as u64), "nonprime subgroup")?;
        let t: i128 = if a == 0 { -1 } else { 1 };
        let (mut s0, mut s1): (i128, i128) = (2, t);
        for _ in 1..n {
            (s0, s1) = (s1, t * s1 - 2 * s0);
        }
        let order = (1i128 << n) + 1 - s1;
        require(
            order == int(field(fixture, "group_order")?, "group_order")?,
            "wrong group order",
        )?;
        let h = int(field(fixture, "cofactor")?, "cofactor")?;
        require(order == h * r, "wrong cofactor")?;
        let mut c = Curve {
            n,
            a: a as u64,
            r: r as u64,
            h: h as u64,
            modulus,
            g: (0, 0),
            lam: 0,
        };
        let g = c.decode(field(fixture, "generator")?)?;
        require(
            g.is_some() && c.mul(g, r as u128).is_none(),
            "bad generator",
        )?;
        c.g = g.expect("checked");
        let lam = int(field(fixture, "lambda")?, "lambda")?;
        require(lam >= 0, "negative scalar")?;
        c.lam = lam as u64;
        require(
            c.mul(g, lam as u128) == c.frob(g),
            "wrong Frobenius eigenvalue",
        )?;
        Ok(c)
    }

    /// A product in the field: shift and add, reduced bit by bit.
    pub fn fm(&self, mut a: u64, mut b: u64) -> u64 {
        let mut z = 0u64;
        while b != 0 {
            if b & 1 == 1 {
                z ^= a;
            }
            b >>= 1;
            a <<= 1;
            if a >> self.n != 0 {
                a ^= self.modulus;
            }
        }
        z
    }

    /// An inverse, by the extended Euclidean algorithm over GF(2)[x].
    pub fn inv(&self, x: u64) -> Result<u64, String> {
        require(x != 0, "zero denominator")?;
        let (mut u, mut v, mut a, mut b) = (x as u128, self.modulus as u128, 1u128, 0u128);
        while u != 1 {
            require(u != 0, "reducible field modulus")?;
            let mut shift = bits(u) as i64 - bits(v) as i64;
            if shift < 0 {
                std::mem::swap(&mut u, &mut v);
                std::mem::swap(&mut a, &mut b);
                shift = -shift;
            }
            u ^= v << shift;
            a ^= b << shift;
        }
        while bits(a) > self.n {
            a ^= (self.modulus as u128) << (bits(a) - self.n - 1);
        }
        Ok(a as u64)
    }

    /// A point from its JSON pair, checked to lie on the curve.
    pub fn decode(&self, p: &J) -> Result<Point, String> {
        if matches!(p, J::Null) {
            return Ok(None);
        }
        let pair = p
            .as_arr()
            .filter(|v| v.len() == 2)
            .ok_or("malformed point")?;
        let (x, y) = (int(&pair[0], "x")?, int(&pair[1], "y")?);
        let top = 1i128 << self.n;
        require(
            (0..top).contains(&x) && (0..top).contains(&y),
            "point outside field",
        )?;
        let (x, y) = (x as u64, y as u64);
        let xx = self.fm(x, x);
        require(
            self.fm(y, y) ^ self.fm(x, y) == self.fm(xx, x) ^ self.fm(self.a, xx) ^ 1,
            "point does not lift to curve",
        )?;
        Ok(Some((x, y)))
    }

    pub fn neg(&self, p: Point) -> Point {
        p.map(|(x, y)| (x, x ^ y))
    }

    pub fn frob(&self, p: Point) -> Point {
        p.map(|(x, y)| (self.fm(x, x), self.fm(y, y)))
    }

    pub fn add(&self, p: Point, q: Point) -> Point {
        let (Some((x, y)), Some((u, v))) = (p, q) else {
            return p.or(q);
        };
        if x == u {
            if y != v || x == 0 {
                return None;
            }
            let lam = x ^ self.fm(y, self.inv(x).expect("x is nonzero"));
            let z = self.fm(lam, lam) ^ lam ^ self.a;
            return Some((z, self.fm(x, x) ^ self.fm(lam ^ 1, z)));
        }
        let lam = self.fm(y ^ v, self.inv(x ^ u).expect("x and u differ"));
        let z = self.fm(lam, lam) ^ lam ^ x ^ u ^ self.a;
        Some((z, self.fm(lam, x ^ z) ^ z ^ y))
    }

    pub fn mul(&self, mut p: Point, mut k: u128) -> Point {
        let mut q = None;
        while k != 0 {
            if k & 1 == 1 {
                q = self.add(q, p);
            }
            p = self.add(p, p);
            k >>= 1;
        }
        q
    }
}

#[cfg(test)]
mod tests {
    use super::super::json::parse;
    use super::*;

    #[test]
    fn primality_and_irreducibility_agree_with_small_cases() {
        let primes: Vec<u64> = (0..60).filter(|&n| prime(n)).collect();
        assert_eq!(
            primes,
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]
        );
        assert!(prime(2_305_843_009_213_693_951)); // 2^61 − 1
        assert!(!prime(3_215_031_751)); // a strong pseudoprime to 2, 3, 5, 7
        assert!(irreducible(0b1011)); // x^3 + x + 1
        assert!(!irreducible(0b101)); // x^2 + 1 = (x + 1)^2
        assert!(irreducible((1 << 41) | (1 << 3) | 1)); // x^41 + x^3 + 1
    }

    #[test]
    fn a_small_curves_group_law_closes() {
        // K_1 over GF(2^5), modulus x^5 + x^2 + 1: 22 = 2 · 11 points.
        let modulus = 0b100101u64;
        let c = Curve {
            n: 5,
            a: 1,
            r: 11,
            h: 2,
            modulus,
            g: (0, 0),
            lam: 0,
        };
        let mut points = vec![];
        for x in 0..32u64 {
            for y in 0..32u64 {
                let p = parse(&format!("[{x}, {y}]")).unwrap();
                if let Ok(Some(q)) = c.decode(&p) {
                    points.push(Some(q));
                }
            }
        }
        assert_eq!(points.len() + 1, 22);
        for &p in &points {
            assert_eq!(c.mul(p, 22), None);
            assert_eq!(c.add(p, c.neg(p)), None);
            assert_eq!(c.add(c.add(p, p), p), c.mul(p, 3));
        }
    }
}
