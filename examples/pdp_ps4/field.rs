//! Type-II optimal normal basis arithmetic with counted operations.
//!
//! An element of `F_{2^n}` is the `n`-bit word of its coordinates in the
//! basis `γ_i = ζ^i + ζ^{-i}`, `1 ≤ i ≤ n`, `ζ^{2n+1} = 1`, bit `i − 1` for
//! `γ_i`: exactly the representation of `Type2Onb` in the library, so every
//! word here is comparable to the library's.  Since `γ_i γ_j = γ_{i+j} +
//! γ_{i−j}` (indices folded into `1..=n`, `γ_0 = 0`), the product is the cyclic
//! convolution of the two symmetric `(2n+1)`-bit words `E(a) = Σ a_i (z^i +
//! z^{-i})`.  Squaring is the coordinate permutation `i ↦ fold(2i)`; it costs
//! no gate and is kept as a table lookup only for speed.
//!
//! The counters mirror the library's `semaev_decomp::opcount` (`mul`, `sqr`,
//! `sqr_k`, `inv`, `batch_*`) so that the control kernel and this field are
//! charged in the same vocabulary, plus `add` (a field addition, one XOR of
//! `n`-bit words), `popcnt` (one weight test) and `lookup` (one table probe).

use crypto_lib::cryptanalysis::koblitz_rotated_chain::Type2Onb;
use std::cell::Cell;

#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct Counts {
    pub mul: u64,
    pub sqr: u64,
    pub sqr_k: u64,
    pub inv: u64,
    pub batch_calls: u64,
    pub batch_elems: u64,
    pub batch_mul: u64,
    pub add: u64,
    pub popcnt: u64,
    pub lookup: u64,
}

impl Counts {
    pub fn sub(&self, o: &Counts) -> Counts {
        Counts {
            mul: self.mul - o.mul,
            sqr: self.sqr - o.sqr,
            sqr_k: self.sqr_k - o.sqr_k,
            inv: self.inv - o.inv,
            batch_calls: self.batch_calls - o.batch_calls,
            batch_elems: self.batch_elems - o.batch_elems,
            batch_mul: self.batch_mul - o.batch_mul,
            add: self.add - o.add,
            popcnt: self.popcnt - o.popcnt,
            lookup: self.lookup - o.lookup,
        }
    }

    pub fn plus(&mut self, o: &Counts) {
        self.mul += o.mul;
        self.sqr += o.sqr;
        self.sqr_k += o.sqr_k;
        self.inv += o.inv;
        self.batch_calls += o.batch_calls;
        self.batch_elems += o.batch_elems;
        self.batch_mul += o.batch_mul;
        self.add += o.add;
        self.popcnt += o.popcnt;
        self.lookup += o.lookup;
    }

    /// Field multiplications in the library's accounting: the products of
    /// `mul` (which include the ones inside `inv`) and the `3` per element
    /// that Montgomery's trick spends.
    pub fn mulq(&self) -> u64 {
        self.mul + self.batch_mul
    }

    pub fn json(&self) -> String {
        format!(
            "{{\"mul\":{},\"sqr\":{},\"sqr_k\":{},\"inv\":{},\"batch_calls\":{},\"batch_elems\":{},\"batch_mul\":{},\"add\":{},\"popcnt\":{},\"lookup\":{},\"mulq\":{}}}",
            self.mul, self.sqr, self.sqr_k, self.inv, self.batch_calls, self.batch_elems,
            self.batch_mul, self.add, self.popcnt, self.lookup, self.mulq()
        )
    }
}

/// What the curve code needs of a field.
pub trait Fld {
    fn n(&self) -> usize;
    fn one(&self) -> u64;
    fn mul(&self, a: u64, b: u64) -> u64;
    fn sqr(&self, a: u64) -> u64;
    fn inv(&self, a: u64) -> u64;
}

pub struct Onb {
    pub n: usize,
    pub mask: u64,
    nn: u32,
    maskn: u128,
    sq: Vec<[u64; 256]>,
    pub lib: Type2Onb,
    c_mul: Cell<u64>,
    c_sqr: Cell<u64>,
    c_sqr_k: Cell<u64>,
    c_inv: Cell<u64>,
    c_batch_calls: Cell<u64>,
    c_batch_elems: Cell<u64>,
    c_batch_mul: Cell<u64>,
    c_add: Cell<u64>,
    c_popcnt: Cell<u64>,
    c_lookup: Cell<u64>,
}

#[inline(always)]
fn bump(c: &Cell<u64>, k: u64) {
    c.set(c.get() + k);
}

impl Onb {
    pub fn new(n: usize) -> Result<Onb, String> {
        if !(3..=60).contains(&n) {
            return Err(format!("n = {n} outside 3..=60"));
        }
        let lib = Type2Onb::new(n)?;
        let mut perm = vec![0usize; n];
        for (t, p) in perm.iter_mut().enumerate() {
            *p = lib.frobenius(1u128 << t, 1).trailing_zeros() as usize;
        }
        let bytes = n.div_ceil(8);
        let mut sq = vec![[0u64; 256]; bytes];
        for (j, table) in sq.iter_mut().enumerate() {
            for (v, slot) in table.iter_mut().enumerate() {
                let mut out = 0u64;
                for k in 0..8 {
                    let t = 8 * j + k;
                    if t < n && (v >> k) & 1 == 1 {
                        out |= 1u64 << perm[t];
                    }
                }
                *slot = out;
            }
        }
        let nn = (2 * n + 1) as u32;
        Ok(Onb {
            n,
            mask: (1u64 << n) - 1,
            nn,
            maskn: (1u128 << nn) - 1,
            sq,
            lib,
            c_mul: Cell::new(0),
            c_sqr: Cell::new(0),
            c_sqr_k: Cell::new(0),
            c_inv: Cell::new(0),
            c_batch_calls: Cell::new(0),
            c_batch_elems: Cell::new(0),
            c_batch_mul: Cell::new(0),
            c_add: Cell::new(0),
            c_popcnt: Cell::new(0),
            c_lookup: Cell::new(0),
        })
    }

    pub fn snapshot(&self) -> Counts {
        Counts {
            mul: self.c_mul.get(),
            sqr: self.c_sqr.get(),
            sqr_k: self.c_sqr_k.get(),
            inv: self.c_inv.get(),
            batch_calls: self.c_batch_calls.get(),
            batch_elems: self.c_batch_elems.get(),
            batch_mul: self.c_batch_mul.get(),
            add: self.c_add.get(),
            popcnt: self.c_popcnt.get(),
            lookup: self.c_lookup.get(),
        }
    }

    // ---- uncounted primitives -------------------------------------------

    #[inline(always)]
    fn embed(&self, b: u64) -> u128 {
        let u = (b as u128) << 1;
        u | (u << (127 - self.nn)).reverse_bits()
    }

    #[inline(always)]
    fn rotl(&self, x: u128, k: u32) -> u128 {
        ((x << k) | (x >> (self.nn - k))) & self.maskn
    }

    /// The product, not counted.
    #[inline]
    pub fn mul_raw(&self, a: u64, b: u64) -> u64 {
        if a == 0 || b == 0 {
            return 0;
        }
        let (a, b) = if a.count_ones() <= b.count_ones() {
            (a, b)
        } else {
            (b, a)
        };
        let e = self.embed(b);
        let mut p = 0u128;
        let mut aa = a;
        while aa != 0 {
            let i = aa.trailing_zeros() + 1;
            aa &= aa - 1;
            p ^= self.rotl(e, i) ^ self.rotl(e, self.nn - i);
        }
        ((p >> 1) as u64) & self.mask
    }

    /// The square, not counted: the coordinate permutation `i ↦ fold(2i)`.
    #[inline(always)]
    pub fn sqr_raw(&self, a: u64) -> u64 {
        let mut out = 0u64;
        for (j, table) in self.sq.iter().enumerate() {
            out ^= table[((a >> (8 * j)) & 255) as usize];
        }
        out
    }

    /// `a^(2^k)`, not counted.
    pub fn frob_raw(&self, mut a: u64, k: usize) -> u64 {
        for _ in 0..(k % self.n) {
            a = self.sqr_raw(a);
        }
        a
    }

    // ---- counted operations ----------------------------------------------

    #[inline(always)]
    pub fn mul_c(&self, a: u64, b: u64) -> u64 {
        bump(&self.c_mul, 1);
        self.mul_raw(a, b)
    }

    #[inline(always)]
    pub fn sqr_c(&self, a: u64) -> u64 {
        bump(&self.c_sqr, 1);
        self.sqr_raw(a)
    }

    #[inline(always)]
    pub fn sqr_k(&self, a: u64, k: u32) -> u64 {
        bump(&self.c_sqr_k, k as u64);
        self.frob_raw(a, k as usize)
    }

    /// One field addition (a word XOR), counted.
    #[inline(always)]
    pub fn add_c(&self, a: u64, b: u64) -> u64 {
        bump(&self.c_add, 1);
        a ^ b
    }

    /// The weight test `wt(t) ≤ w`, one popcount.
    #[inline(always)]
    pub fn weight_le(&self, t: u64, w: u32) -> bool {
        bump(&self.c_popcnt, 1);
        t.count_ones() <= w
    }

    #[inline(always)]
    pub fn note_lookup(&self) {
        bump(&self.c_lookup, 1);
    }

    /// Itoh–Tsujii, as the library's `Gf2::inv`: `⌊log₂(n−1)⌋ + wt(n−1) − 1`
    /// products and `n − 1` squarings (free here).
    pub fn inv_c(&self, a: u64) -> u64 {
        bump(&self.c_inv, 1);
        if a == 0 {
            return 0;
        }
        let e = (self.n - 1) as u32;
        let mut beta = a;
        let mut len = 1u32;
        for bit in (0..(31 - e.leading_zeros())).rev() {
            beta = self.mul_c(self.sqr_k(beta, len), beta);
            len *= 2;
            if (e >> bit) & 1 == 1 {
                beta = self.mul_c(self.sqr_c(beta), a);
                len += 1;
            }
        }
        debug_assert_eq!(len, e);
        self.sqr_c(beta)
    }

    /// Montgomery's trick, charged as the library charges it: `3` products
    /// per non-zero element in `batch_mul`, one `inv`.
    pub fn batch_inv(&self, xs: &mut [u64], scratch: &mut Vec<u64>) {
        let nz = xs.iter().filter(|&&x| x != 0).count() as u64;
        bump(&self.c_batch_calls, 1);
        bump(&self.c_batch_elems, nz);
        bump(&self.c_batch_mul, 3 * nz);
        if nz == 0 {
            return;
        }
        scratch.clear();
        let mut acc = self.mask; // the field's 1
        for &x in xs.iter() {
            if x != 0 {
                acc = self.mul_raw(acc, x);
            }
            scratch.push(acc);
        }
        let mut inv_acc = self.inv_c(acc);
        for i in (0..xs.len()).rev() {
            if xs[i] != 0 {
                let prev = if i > 0 { scratch[i - 1] } else { self.mask };
                let inv_i = self.mul_raw(inv_acc, prev);
                inv_acc = self.mul_raw(inv_acc, xs[i]);
                xs[i] = inv_i;
            }
        }
    }
}

impl Fld for Onb {
    fn n(&self) -> usize {
        self.n
    }
    fn one(&self) -> u64 {
        self.mask
    }
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.mul_c(a, b)
    }
    #[inline]
    fn sqr(&self, a: u64) -> u64 {
        self.sqr_c(a)
    }
    fn inv(&self, a: u64) -> u64 {
        self.inv_c(a)
    }
}

/// A small field with full tables, for the exhaustive checks (`n ≤ 11`).
pub struct Small {
    pub f: Onb,
    pub size: usize,
    mulq: Vec<u16>,
    sqt: Vec<u16>,
    invt: Vec<u16>,
}

impl Small {
    pub fn new(n: usize) -> Result<Small, String> {
        if n > 11 {
            return Err("Small supports n <= 11".into());
        }
        let f = Onb::new(n)?;
        let size = 1usize << n;
        let mut mulq = vec![0u16; size * size];
        for a in 0..size {
            for b in a..size {
                let v = f.mul_raw(a as u64, b as u64) as u16;
                mulq[a * size + b] = v;
                mulq[b * size + a] = v;
            }
        }
        let sqt: Vec<u16> = (0..size).map(|a| f.sqr_raw(a as u64) as u16).collect();
        let mut invt = vec![0u16; size];
        for a in 1..size {
            for b in a..size {
                if mulq[a * size + b] as u64 == f.mask {
                    invt[a] = b as u16;
                    invt[b] = a as u16;
                }
            }
        }
        Ok(Small {
            f,
            size,
            mulq,
            sqt,
            invt,
        })
    }
}

impl Fld for Small {
    fn n(&self) -> usize {
        self.f.n
    }
    fn one(&self) -> u64 {
        self.f.mask
    }
    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.mulq[a as usize * self.size + b as usize] as u64
    }
    #[inline(always)]
    fn sqr(&self, a: u64) -> u64 {
        self.sqt[a as usize] as u64
    }
    #[inline(always)]
    fn inv(&self, a: u64) -> u64 {
        self.invt[a as usize] as u64
    }
}

/// `⌊log₂(n−1)⌋ + popcount(n−1) − 1`: the products of one Itoh–Tsujii inversion.
pub fn itoh_tsujii_products(n: usize) -> u32 {
    let e = (n - 1) as u32;
    (31 - e.leading_zeros()) + e.count_ones() - 1
}
