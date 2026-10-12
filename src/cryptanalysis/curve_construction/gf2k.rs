//! Small binary fields `F_{2^k}` (`1 ≤ k ≤ 16`) with log/antilog tables.
//!
//! Elements are `u16` bit vectors in a polynomial basis whose modulus is the
//! first primitive polynomial of degree `k` (searched, not tabulated), so `z`
//! generates the multiplicative group and multiplication is two table reads.

/// `F_{2^k}` with table arithmetic.
#[derive(Clone, Debug)]
pub struct SmallField {
    pub k: u32,
    pub q: usize,
    /// Primitive modulus, bit `i` = coefficient of `z^i`, leading bit included.
    pub poly: u32,
    exp: Vec<u16>,
    log: Vec<u16>,
    trace: Vec<u8>,
    /// `half[c]` is a root of `z^2 + z = c`, or `u16::MAX` if there is none.
    half: Vec<u16>,
}

impl SmallField {
    pub fn new(k: u32) -> Self {
        assert!((1..=16).contains(&k), "k must be in 1..=16");
        let q = 1usize << k;
        let order = q - 1;
        let mut found = None;
        for poly in (1u32 << k)..(1u32 << (k + 1)) {
            if poly & 1 == 0 {
                continue;
            }
            let mut exp = vec![0u16; 2 * order];
            let mut x = 1u32;
            let mut primitive = true;
            for (i, slot) in exp.iter_mut().take(order).enumerate() {
                if i > 0 && x == 1 {
                    primitive = false;
                    break;
                }
                *slot = x as u16;
                x <<= 1;
                if (x >> k) & 1 == 1 {
                    x ^= poly;
                }
            }
            if primitive && x == 1 {
                found = Some((poly, exp));
                break;
            }
        }
        let (poly, mut exp) = found.expect("a primitive polynomial exists for every degree");
        for i in order..2 * order {
            exp[i] = exp[i - order];
        }
        let mut log = vec![0u16; q];
        for (i, &e) in exp.iter().take(order).enumerate() {
            log[e as usize] = i as u16;
        }
        let mut f = SmallField {
            k,
            q,
            poly,
            exp,
            log,
            trace: vec![0; q],
            half: vec![u16::MAX; q],
        };
        for a in 0..q {
            let mut acc = 0u16;
            let mut c = a as u16;
            for _ in 0..k {
                acc ^= c;
                c = f.mul(c, c);
            }
            assert!(acc <= 1, "absolute trace must lie in F_2");
            f.trace[a] = acc as u8;
        }
        for z in 0..q {
            let z = z as u16;
            let c = f.mul(z, z) ^ z;
            if f.half[c as usize] == u16::MAX {
                f.half[c as usize] = z;
            }
        }
        f
    }

    #[inline]
    pub fn mul(&self, a: u16, b: u16) -> u16 {
        if a == 0 || b == 0 {
            0
        } else {
            self.exp[self.log[a as usize] as usize + self.log[b as usize] as usize]
        }
    }

    #[inline]
    pub fn inv(&self, a: u16) -> u16 {
        assert!(a != 0, "zero has no inverse");
        let order = self.q - 1;
        self.exp[(order - self.log[a as usize] as usize) % order]
    }

    #[inline]
    pub fn sqr(&self, a: u16) -> u16 {
        self.mul(a, a)
    }

    pub fn pow(&self, a: u16, e: u64) -> u16 {
        if e == 0 {
            return 1;
        }
        if a == 0 {
            return 0;
        }
        let order = (self.q - 1) as u64;
        self.exp[((self.log[a as usize] as u64 * (e % order)) % order) as usize]
    }

    #[inline]
    pub fn trace(&self, a: u16) -> u8 {
        self.trace[a as usize]
    }

    /// A root `z` of `z^2 + z = c`, if one exists in this field.
    #[inline]
    pub fn solve_as(&self, c: u16) -> Option<u16> {
        let z = self.half[c as usize];
        (z != u16::MAX).then_some(z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_axioms_small_degrees() {
        for k in 1..=9 {
            let f = SmallField::new(k);
            for a in 1..f.q as u16 {
                assert_eq!(f.mul(a, f.inv(a)), 1, "k={k} a={a}");
            }
            // the trace is F_2-linear and onto for k >= 1
            let ones = (0..f.q as u16).filter(|&a| f.trace(a) == 1).count();
            assert_eq!(ones, f.q / 2);
            // z^2 + z = c is solvable exactly when Tr(c) = 0
            for c in 0..f.q as u16 {
                assert_eq!(f.solve_as(c).is_some(), f.trace(c) == 0);
            }
        }
    }
}
