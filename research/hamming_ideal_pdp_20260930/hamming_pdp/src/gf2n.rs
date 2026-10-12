//! `F_{2^n}` on `u64`, `1 ≤ n ≤ 63`, polynomial basis, plus a normal basis.
//!
//! `n = 83` does not fit. The irreducible needs bit `n` set in a `u64`, and
//! the normal-basis inverse is a `u128` row of width `2n` (`n + (n − 1) < 128`).
//! Callers that scan `1 << n` are gated on their own; this type stops at 63.

#[derive(Clone, Debug)]
pub struct Field {
    pub n: u32,
    /// Reduction polynomial with its top bit set.
    pub irr: u64,
    pub mask: u64,
    /// Normal element `α`: `nb[i] = α^{2^i}` in polynomial-basis bits.
    pub nb: Vec<u64>,
    /// `to_nb[i]`: row `i` of the inverse change-of-basis; the `i`-th normal
    /// coordinate of `x` is `parity(x & to_nb[i])`.
    pub to_nb: Vec<u64>,
}

fn clmul_mod(mut a: u64, mut b: u64, irr: u64, n: u32) -> u64 {
    let mut r = 0u64;
    let top = 1u64 << n;
    while b != 0 {
        if b & 1 == 1 {
            r ^= a;
        }
        b >>= 1;
        a <<= 1;
        if a & top != 0 {
            a ^= irr;
        }
    }
    r
}

pub fn is_irreducible(f: u64, n: u32) -> bool {
    // Rabin: x^(2^n) = x mod f and gcd(x^(2^(n/p)) - x, f) = 1 for primes p | n.
    let powx = |k: u32| -> u64 {
        let mut h = 2u64;
        for _ in 0..k {
            h = clmul_mod(h, h, f, n);
        }
        h
    };
    if powx(n) != 2 {
        return false;
    }
    let mut primes = Vec::new();
    let mut m = n;
    let mut p = 2;
    while p * p <= m {
        if m % p == 0 {
            primes.push(p);
            while m % p == 0 {
                m /= p;
            }
        }
        p += 1;
    }
    if m > 1 {
        primes.push(m);
    }
    for p in primes {
        let g = poly_gcd(powx(n / p) ^ 2, f);
        if g != 1 {
            return false;
        }
    }
    true
}

fn deg(a: u64) -> i32 {
    63 - a.leading_zeros() as i32
}

fn poly_gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        while a != 0 && deg(a) >= deg(b) {
            a ^= b << (deg(a) - deg(b));
        }
        std::mem::swap(&mut a, &mut b);
    }
    a
}

/// The `find_irreducible` convention of `scripts/ecc2k130_point_decomposition.py`.
pub fn find_irreducible(n: u32) -> u64 {
    for extra in 1u64.. {
        let f = (1u64 << n) | (extra << 1) | 1;
        if deg(f) as u32 == n && is_irreducible(f, n) {
            return f;
        }
    }
    unreachable!()
}

/// Rank of a set of bit-vectors over `F_2`.
pub fn rank_f2(rows: &[u64]) -> usize {
    let mut rows = rows.to_vec();
    let mut r = 0;
    for bit in (0..64).rev() {
        if let Some(p) = (r..rows.len()).find(|&i| (rows[i] >> bit) & 1 == 1) {
            rows.swap(r, p);
            let piv = rows[r];
            for (i, row) in rows.iter_mut().enumerate() {
                if i != r && (*row >> bit) & 1 == 1 {
                    *row ^= piv;
                }
            }
            r += 1;
        }
    }
    r
}

impl Field {
    pub fn new(n: u32) -> Field {
        assert!(
            (1..64).contains(&n),
            "F_2^n is a u64 polynomial basis; n = 83 needs a multi-limb field"
        );
        let irr = find_irreducible(n);
        let mut f = Field {
            n,
            irr,
            mask: (1u64 << n) - 1,
            nb: vec![],
            to_nb: vec![],
        };
        // Smallest normal element in integer order.
        for alpha in 2u64..(1u64 << n) {
            let mut conj = Vec::with_capacity(n as usize);
            let mut a = alpha;
            for _ in 0..n {
                conj.push(a);
                a = f.sqr(a);
            }
            if rank_f2(&conj) == n as usize {
                f.nb = conj;
                break;
            }
        }
        assert!(!f.nb.is_empty(), "no normal element found");
        // Invert the change of basis: columns of B are nb[i]; solve B u = x.
        // Build the matrix M with M[j][i] = bit j of nb[i]; invert it.
        let nn = n as usize;
        let mut aug: Vec<u128> = (0..nn)
            .map(|j| {
                let mut row = 0u128;
                for i in 0..nn {
                    if (f.nb[i] >> j) & 1 == 1 {
                        row |= 1u128 << i;
                    }
                }
                row | (1u128 << (nn + j))
            })
            .collect();
        for col in 0..nn {
            let p = (col..nn)
                .find(|&i| (aug[i] >> col) & 1 == 1)
                .expect("singular");
            aug.swap(col, p);
            let piv = aug[col];
            for i in 0..nn {
                if i != col && (aug[i] >> col) & 1 == 1 {
                    aug[i] ^= piv;
                }
            }
        }
        // Now aug[i] = e_i | (row i of M^{-1}) << nn.  u_i = Σ_j Minv[i][j] x_j.
        f.to_nb = (0..nn).map(|i| ((aug[i] >> nn) as u64) & f.mask).collect();
        // Self-check: Frobenius is a cyclic shift of normal coordinates.
        for x in [1u64, 3, 5, 7, 11, 13] {
            let x = x & f.mask;
            let u = f.to_normal(x);
            let u2 = f.to_normal(f.sqr(x));
            let shifted = ((u << 1) | (u >> (n - 1))) & f.mask;
            assert_eq!(
                u2, shifted,
                "normal basis: Frobenius must be a cyclic shift"
            );
            assert_eq!(f.from_normal(u), x);
        }
        f
    }

    #[inline]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        clmul_mod(a, b, self.irr, self.n)
    }
    #[inline]
    pub fn sqr(&self, a: u64) -> u64 {
        self.mul(a, a)
    }
    pub fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.sqr(a);
            e >>= 1;
        }
        r
    }
    pub fn inv(&self, a: u64) -> u64 {
        assert!(a != 0);
        self.pow(a, (1u64 << self.n) - 2)
    }
    pub fn trace(&self, a: u64) -> u64 {
        let mut s = 0;
        let mut v = a;
        for _ in 0..self.n {
            s ^= v;
            v = self.sqr(v);
        }
        s
    }
    /// For odd `n`: `z` with `z² + z = c` when `Tr(c) = 0`.
    pub fn half_trace(&self, c: u64) -> u64 {
        assert!(self.n % 2 == 1);
        let mut s = 0;
        let mut v = c;
        for _ in 0..=(self.n - 1) / 2 {
            s ^= v;
            v = self.sqr(self.sqr(v));
        }
        s
    }
    pub fn to_normal(&self, x: u64) -> u64 {
        let mut u = 0u64;
        for (i, &row) in self.to_nb.iter().enumerate() {
            if (x & row).count_ones() & 1 == 1 {
                u |= 1 << i;
            }
        }
        u
    }
    pub fn from_normal(&self, u: u64) -> u64 {
        let mut x = 0;
        for (i, &c) in self.nb.iter().enumerate() {
            if (u >> i) & 1 == 1 {
                x ^= c;
            }
        }
        x
    }
}
