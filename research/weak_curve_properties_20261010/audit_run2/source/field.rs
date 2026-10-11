//! Polynomial-basis arithmetic for the hash-bound odd-characteristic fixtures.
use num_bigint::BigUint;
use num_traits::Zero;
pub type Element = Vec<u64>;
#[derive(Clone)]
pub struct Field {
    pub p: u64,
    pub n: usize,
    pub modulus: Element,
    pub q: BigUint,
    frobenius: Vec<Element>,
}
impl Field {
    pub fn new(p: u64, modulus: Element) -> Self {
        let n = modulus.len();
        assert!(p > 3 && p < (1 << 52) && (1..=16).contains(&n));
        assert!(modulus.iter().all(|c| *c < p));
        let q = BigUint::from(p).pow(n as u32);
        let mut f = Self {
            p,
            n,
            modulus,
            q,
            frobenius: vec![],
        };
        let x = f.generator();
        let xp = f.pow64(&x, p);
        let mut column = f.one();
        for _ in 0..n {
            f.frobenius.push(column.clone());
            column = f.mul(&column, &xp);
        }
        assert!(f.irreducible(), "recorded modulus is reducible");
        f
    }
    pub fn zero(&self) -> Element {
        vec![0; self.n]
    }
    pub fn one(&self) -> Element {
        self.constant(1)
    }
    pub fn constant(&self, c: u64) -> Element {
        let mut v = self.zero();
        v[0] = c % self.p;
        v
    }
    pub fn generator(&self) -> Element {
        if self.n == 1 {
            return self.constant((self.p - self.modulus[0]) % self.p);
        }
        let mut v = self.zero();
        v[1] = 1;
        v
    }
    pub fn add(&self, a: &Element, b: &Element) -> Element {
        a.iter().zip(b).map(|(a, b)| (a + b) % self.p).collect()
    }
    pub fn neg(&self, a: &Element) -> Element {
        a.iter().map(|a| (self.p - a) % self.p).collect()
    }
    pub fn sub(&self, a: &Element, b: &Element) -> Element {
        self.add(a, &self.neg(b))
    }
    pub fn mul(&self, a: &Element, b: &Element) -> Element {
        assert_eq!(a.len(), self.n);
        assert_eq!(b.len(), self.n);
        let mut sums = vec![0u128; 2 * self.n - 1];
        for i in 0..self.n {
            for j in 0..self.n {
                sums[i + j] += a[i] as u128 * b[j] as u128;
            }
        }
        let mut v: Vec<u64> = sums.iter().map(|x| (x % self.p as u128) as u64).collect();
        for k in (self.n..v.len()).rev() {
            let c = v[k];
            for i in 0..self.n {
                let sub = (c as u128 * self.modulus[i] as u128 % self.p as u128) as u64;
                v[k - self.n + i] = (v[k - self.n + i] + self.p - sub) % self.p;
            }
        }
        v.truncate(self.n);
        v
    }
    pub fn pow64(&self, a: &Element, mut e: u64) -> Element {
        let mut b = a.clone();
        let mut r = self.one();
        while e != 0 {
            if e & 1 != 0 {
                r = self.mul(&r, &b);
            }
            e >>= 1;
            if e != 0 {
                b = self.mul(&b, &b);
            }
        }
        r
    }
    pub fn pow(&self, a: &Element, e: &BigUint) -> Element {
        let mut e = e.clone();
        let mut b = a.clone();
        let mut r = self.one();
        while !e.is_zero() {
            if e.bit(0) {
                r = self.mul(&r, &b);
            }
            e >>= 1usize;
            if !e.is_zero() {
                b = self.mul(&b, &b);
            }
        }
        r
    }
    pub fn inv(&self, a: &Element) -> Element {
        assert_ne!(*a, self.zero());
        self.pow(a, &(&self.q - 2u32))
    }
    pub fn square(&self, a: &Element) -> bool {
        *a == self.zero() || self.pow(a, &((&self.q - 1u32) >> 1usize)) == self.one()
    }
    pub fn frob(&self, a: &Element) -> Element {
        let mut sums = vec![0u128; self.n];
        for (i, a) in a.iter().enumerate() {
            for (j, col) in self.frobenius[i].iter().enumerate() {
                sums[j] += *a as u128 * *col as u128;
            }
        }
        sums.into_iter()
            .map(|x| (x % self.p as u128) as u64)
            .collect()
    }
    pub fn degree(&self, values: &[Element]) -> usize {
        let mut v = values.to_vec();
        for d in 1..=self.n {
            for a in &mut v {
                *a = self.frob(a);
            }
            if v == values {
                return d;
            }
        }
        panic!("Frobenius period mismatch")
    }
    fn irreducible(&self) -> bool {
        let x = self.generator();
        let mut v = x.clone();
        let mut powers = vec![v.clone()];
        for _ in 0..self.n {
            v = self.frob(&v);
            powers.push(v.clone());
        }
        if v != x {
            return false;
        }
        let mut modulus = self.modulus.clone();
        modulus.push(1);
        let mut remaining = self.n;
        let mut primes = vec![];
        for r in 2..=self.n {
            if remaining % r == 0 {
                primes.push(r);
                while remaining % r == 0 {
                    remaining /= r;
                }
            }
        }
        for r in primes {
            let g = self.sub(&powers[self.n / r], &x);
            if gcd_degree(self.p, modulus.clone(), g) != 0 {
                return false;
            }
        }
        true
    }
}
fn trim(v: &mut Element) {
    while v.last() == Some(&0) {
        v.pop();
    }
}
fn scalar_pow(p: u64, mut a: u64, mut e: u64) -> u64 {
    let mut r = 1;
    while e != 0 {
        if e & 1 != 0 {
            r = (r as u128 * a as u128 % p as u128) as u64;
        }
        e >>= 1;
        a = (a as u128 * a as u128 % p as u128) as u64;
    }
    r
}
fn gcd_degree(p: u64, mut a: Element, mut b: Element) -> usize {
    trim(&mut a);
    trim(&mut b);
    while !b.is_empty() {
        let inverse = scalar_pow(p, *b.last().unwrap(), p - 2);
        while a.len() >= b.len() {
            let at = a.len() - b.len();
            let c = (*a.last().unwrap() as u128 * inverse as u128 % p as u128) as u64;
            for (i, b) in b.iter().enumerate() {
                let d = (c as u128 * *b as u128 % p as u128) as u64;
                a[at + i] = (a[at + i] + p - d) % p;
            }
            trim(&mut a);
        }
        std::mem::swap(&mut a, &mut b);
    }
    a.len() - 1
}
pub fn divisors(n: usize) -> Vec<usize> {
    (1..=n).filter(|d| n % d == 0).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn prime_and_extension_arithmetic() {
        let f = Field::new(7, vec![1, 0]);
        for a in 0..7 {
            for b in 0..7 {
                let v = vec![a, b];
                assert_eq!(f.frob(&f.frob(&v)), v);
                if v != f.zero() {
                    assert_eq!(f.mul(&v, &f.inv(&v)), f.one());
                }
            }
        }
        assert_eq!(f.degree(&[f.generator()]), 2);
    }
    #[test]
    #[should_panic(expected = "reducible")]
    fn reject_reducible_modulus() {
        Field::new(7, vec![6, 0]);
    }
}
