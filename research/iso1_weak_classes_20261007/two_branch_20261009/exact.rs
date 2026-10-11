//! Independent complete affine-x counts in the saved PARI polynomial representation.
use serde_json::Value;
use std::{collections::BTreeMap, fs};
const Q: usize = 117649;
#[derive(Clone)]
struct Field {
    modulus: [u32; 6],
}
fn decode(mut code: u32) -> [u32; 6] {
    let mut a = [0; 6];
    for v in &mut a {
        *v = code % 7;
        code /= 7;
    }
    a
}
fn encode(a: [u32; 6]) -> u32 {
    a.iter().rev().fold(0, |s, &v| 7 * s + v)
}
impl Field {
    fn add(&self, a: u32, b: u32) -> u32 {
        let a = decode(a);
        let b = decode(b);
        encode(std::array::from_fn(|i| (a[i] + b[i]) % 7))
    }
    fn mul(&self, a: u32, b: u32) -> u32 {
        let a = decode(a);
        let b = decode(b);
        let mut c = [0u32; 11];
        for (i, &x) in a.iter().enumerate() {
            for (j, &y) in b.iter().enumerate() {
                c[i + j] += x * y;
            }
        }
        for k in (6..=10).rev() {
            let v = c[k] % 7;
            for i in 0..6 {
                c[k - 6 + i] += v * (7 - self.modulus[i]);
            }
        }
        encode(std::array::from_fn(|i| c[i] % 7))
    }
}
fn kv(s: &str) -> BTreeMap<&str, &str> {
    s.split('|')
        .skip(1)
        .filter_map(|v| v.split_once('='))
        .collect()
}
fn arr(s: &str) -> [u32; 6] {
    let v: Value = serde_json::from_str(s).unwrap();
    std::array::from_fn(|i| v[i].as_str().unwrap().parse::<u32>().unwrap())
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 1);
    let raw = fs::read_to_string(&a[0]).unwrap();
    let meta = kv(raw.lines().next().unwrap());
    assert_eq!(meta["p"], "7");
    let f = Field {
        modulus: arr(meta["modulus"]),
    };
    let mut chi = vec![-1i64; Q];
    chi[0] = 0;
    let mut cubes = vec![0u32; Q];
    for x in 0..Q as u32 {
        let sq = f.mul(x, x);
        chi[sq as usize] = 1;
        cubes[x as usize] = f.mul(sq, x);
    }
    chi[0] = 0;
    assert_eq!(chi.iter().filter(|&&c| c == 1).count(), (Q - 1) / 2);
    let mut count = 0;
    let mut minus = 0;
    for line in raw
        .lines()
        .filter(|l| l.starts_with("PLUS|") || l.starts_with("MINUS|"))
    {
        let m = kv(line);
        let ca = encode(arr(m["short_a"]));
        let cb = encode(arr(m["short_b"]));
        let t: i64 = m["trace"].parse().unwrap();
        let mut sum = 0i64;
        for (x, &cube) in cubes.iter().enumerate() {
            let value = f.add(f.add(cube, f.mul(ca, x as u32)), cb);
            sum += chi[value as usize];
        }
        assert_eq!(
            -sum,
            t,
            "independent all-x mismatch at {}",
            line.split('|').next().unwrap()
        );
        count += 1;
        minus += usize::from(line.starts_with("MINUS|"));
    }
    assert_eq!((count, minus), (409, 196));
    println!("EXACT_COMPLETE|p=7|models={count}|quadratic_models={minus}|affine_x_evaluations={}|all_traces_match=1",count*Q);
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn multiplication_and_character() {
        let f = Field {
            modulus: [5, 6, 5, 3, 5, 5],
        };
        for x in [0, 1, 7, 202, 117648] {
            assert_eq!(f.mul(x, 1), x);
            assert_eq!(f.mul(x, 0), 0);
            assert_eq!(f.mul(x, 7), f.mul(7, x));
        }
        let z6 = (0..6).fold(1, |v, _| f.mul(v, 7));
        assert_eq!(z6, encode([2, 1, 2, 4, 2, 2]));
    }
}
