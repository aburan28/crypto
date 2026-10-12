//! Independent polynomial-basis and group arithmetic for frozen IC receipts.
//! This module intentionally does not call the producer's Koblitz operations.

use crypto_lib::hash::sha256::Sha256;
use serde_json::Value;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub type Check<T> = Result<T, String>;
pub type Point = Option<(u64, u64)>;

pub fn ensure(ok: bool, message: impl Into<String>) -> Check<()> {
    if ok {
        Ok(())
    } else {
        Err(message.into())
    }
}

pub fn integer(v: &Value, field: &str) -> Check<u64> {
    v[field]
        .as_u64()
        .ok_or_else(|| format!("missing integer {field}"))
}

pub fn array<'a>(v: &'a Value, field: &str) -> Check<&'a [Value]> {
    v[field]
        .as_array()
        .map(Vec::as_slice)
        .ok_or_else(|| format!("missing array {field}"))
}

pub fn point(v: &Value) -> Check<Point> {
    if v.is_null() {
        return Ok(None);
    }
    let a = v.as_array().ok_or("point is not a pair")?;
    ensure(a.len() == 2, "point is not a pair")?;
    Ok(Some((
        a[0].as_u64().ok_or("invalid x")?,
        a[1].as_u64().ok_or("invalid y")?,
    )))
}

pub fn pair(v: &Value) -> Check<(u64, u64)> {
    point(v)?.ok_or_else(|| "unexpected point at infinity".to_string())
}

#[cfg(test)]
pub fn rows(path: &Path) -> Check<Vec<Value>> {
    let input = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    BufReader::new(input)
        .lines()
        .enumerate()
        .filter_map(|(i, line)| match line {
            Ok(s) if s.trim().is_empty() => None,
            other => Some((i, other)),
        })
        .map(|(i, line)| {
            let line = line.map_err(|e| format!("{}:{}: {e}", path.display(), i + 1))?;
            serde_json::from_str(&line).map_err(|e| format!("{}:{}: {e}", path.display(), i + 1))
        })
        .collect()
}

pub fn one(path: &Path, kind: &str) -> Check<Value> {
    let input = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut found = None;
    for (i, line) in BufReader::new(input).lines().enumerate() {
        let line = line.map_err(|e| format!("{}:{}: {e}", path.display(), i + 1))?;
        if line.trim().is_empty() {
            continue;
        }
        let row: Value = serde_json::from_str(&line)
            .map_err(|e| format!("{}:{}: {e}", path.display(), i + 1))?;
        if row["kind"] == kind {
            ensure(
                found.is_none(),
                format!("duplicate {kind} in {}", path.display()),
            )?;
            found = Some(row);
        }
    }
    found.ok_or_else(|| format!("expected one {kind} in {}", path.display()))
}

pub fn sha256_file(path: &Path) -> Check<String> {
    let mut file = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut hash = Sha256::new();
    let mut buf = [0u8; 64 * 1024];
    loop {
        let count = file.read(&mut buf).map_err(|e| e.to_string())?;
        if count == 0 {
            break;
        }
        hash.update(&buf[..count]);
    }
    Ok(hex::encode(hash.finalize()))
}

#[derive(Clone, Copy)]
pub struct BinaryCurve {
    pub n: u32,
    pub a: u64,
    modulus: u128,
    mask: u64,
}

impl BinaryCurve {
    pub fn new(n: u32, a: u64, terms: &[Value]) -> Check<Self> {
        ensure((2..=63).contains(&n), "unsupported polynomial degree")?;
        ensure(a <= 1, "unsupported binary curve a")?;
        let mut modulus = 1u128 << n;
        for t in terms {
            let bit = t.as_u64().ok_or("invalid modulus term")?;
            ensure(bit < u64::from(n), "modulus term out of range")?;
            modulus ^= 1u128 << bit;
        }
        Ok(Self {
            n,
            a,
            modulus,
            mask: (1u64 << n) - 1,
        })
    }

    pub fn mul(self, a: u64, mut b: u64) -> u64 {
        let mut a = u128::from(a);
        let mut out = 0u128;
        while b != 0 {
            if b & 1 != 0 {
                out ^= a;
            }
            b >>= 1;
            a <<= 1;
            if a & (1u128 << self.n) != 0 {
                a ^= self.modulus;
            }
        }
        (out as u64) & self.mask
    }

    fn inv(self, a: u64) -> Check<u64> {
        ensure(a != 0, "division by zero")?;
        let (mut u, mut v, mut left, mut right) = (u128::from(a), self.modulus, 1u128, 0u128);
        while u != 1 {
            ensure(u != 0, "reducible field modulus")?;
            let mut shift = u.ilog2() as i32 - v.ilog2() as i32;
            if shift < 0 {
                std::mem::swap(&mut u, &mut v);
                std::mem::swap(&mut left, &mut right);
                shift = -shift;
            }
            u ^= v << shift;
            left ^= right << shift;
        }
        while left.ilog2() >= self.n {
            left ^= self.modulus << (left.ilog2() - self.n);
        }
        Ok(left as u64)
    }

    pub fn on_curve(self, p: Point) -> bool {
        let Some((x, y)) = p else {
            return true;
        };
        x <= self.mask
            && y <= self.mask
            && (self.mul(y, y) ^ self.mul(x, y))
                == (self.mul(self.mul(x, x), x)
                    ^ (if self.a == 1 { self.mul(x, x) } else { 0 })
                    ^ 1)
    }

    pub fn neg(self, p: Point) -> Point {
        p.map(|(x, y)| (x, x ^ y))
    }
    pub fn frob(self, p: Point) -> Point {
        p.map(|(x, y)| (self.mul(x, x), self.mul(y, y)))
    }

    pub fn add(self, p: Point, q: Point) -> Check<Point> {
        let (Some((x1, y1)), Some((x2, y2))) = (p, q) else {
            return Ok(p.or(q));
        };
        if x1 == x2 && y1 ^ y2 == x1 {
            return Ok(None);
        }
        let (slope, x3, y3) = if x1 == x2 {
            ensure(y1 == y2 && x1 != 0, "invalid doubling pair")?;
            let slope = x1 ^ self.mul(y1, self.inv(x1)?);
            let x3 = self.mul(slope, slope) ^ slope ^ self.a;
            let y3 = self.mul(x1, x1) ^ self.mul(slope ^ 1, x3);
            (slope, x3, y3)
        } else {
            let slope = self.mul(y1 ^ y2, self.inv(x1 ^ x2)?);
            let x3 = self.mul(slope, slope) ^ slope ^ x1 ^ x2 ^ self.a;
            let y3 = self.mul(slope, x1 ^ x3) ^ x3 ^ y1;
            (slope, x3, y3)
        };
        let _ = slope;
        Ok(Some((x3, y3)))
    }

    pub fn scale(self, mut k: u64, mut p: Point) -> Check<Point> {
        let mut out = None;
        while k != 0 {
            if k & 1 != 0 {
                out = self.add(out, p)?;
            }
            p = self.add(p, p)?;
            k >>= 1;
        }
        Ok(out)
    }

    pub fn sum(self, points: impl IntoIterator<Item = Point>) -> Check<Point> {
        points.into_iter().try_fold(None, |acc, p| self.add(acc, p))
    }
}

pub fn mod_add(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 + b as u128) % r as u128) as u64
}
pub fn mod_sub(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 + r as u128 - b as u128) % r as u128) as u64
}
pub fn mod_mul(a: u64, b: u64, r: u64) -> u64 {
    ((a as u128 * b as u128) % r as u128) as u64
}
pub fn mod_pow(mut a: u64, mut e: u64, r: u64) -> u64 {
    let mut out = 1;
    while e != 0 {
        if e & 1 != 0 {
            out = mod_mul(out, a, r);
        }
        a = mod_mul(a, a, r);
        e >>= 1;
    }
    out
}
