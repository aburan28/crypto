//! PS-8: the meet-in-the-middle pair table (`protocols/PS-8.md`).
//!
//! The table holds `x(Pᵢ + Pⱼ)` for `i ≤ j` and `x(Pᵢ − Pⱼ)` for `i < j` over the rational base
//! points `B` of `V`, one `(i, j)` payload each.  A target `R` is probed once per signed base
//! point: `R = ±Pᵢ ± Pⱼ ± P_k` exactly when `x(R ∓ P_k)` equals a key `x(Pᵢ ± Pⱼ)`, so the two
//! abscissae `x(R + P_k)`, `x(R − P_k)` are looked up for every `k`.  The keys do not depend on
//! `R`; the group-law kernels of PS-2/PS-3 recompute the equivalent quantities per target.
//! Same exactness conditions as the kernels: the caller sends `degenerate_target`s to the
//! baseline, and every hit is confirmed with the exact reference [`thirds`].

use crate::group::*;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;

const EMPTY: u64 = u64::MAX;

pub struct MitmTable<'a> {
    pub gf: &'a Gf2,
    pub l: u32,
    /// Rational base points `(x, y)`, `x` ascending (the same `B` as the group kernels).
    pub base: Vec<(u64, u64)>,
    keys: Vec<u64>,
    /// `(i << 16) | j` of the pair that produced the key.
    vals: Vec<u32>,
    mask: usize,
    shift: u32,
    /// Keys stored.
    pub entries: u64,
}

#[inline(always)]
fn slot(key: u64, shift: u32) -> usize {
    (key.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> shift) as usize
}

impl<'a> MitmTable<'a> {
    pub fn new(gf: &'a Gf2, solver: &ArtinSchreier, l: u32) -> Self {
        let base = base_points(gf, solver, l);
        let m = base.len();
        assert!(m < (1 << 16), "pair index does not fit in 16 bits");
        // At most m² keys; load factor at most 0.7.
        let want = (m * m * 10 / 7).max(16);
        let bits = usize::BITS - (want - 1).leading_zeros();
        let size = 1usize << bits;
        let mut t = Self {
            gf,
            l,
            base,
            keys: vec![EMPTY; size],
            vals: vec![0; size],
            mask: size - 1,
            shift: 64 - bits,
            entries: 0,
        };
        t.build();
        t
    }

    /// Bytes of the table (keys and payloads).
    pub fn bytes(&self) -> u64 {
        (self.keys.len() * 8 + self.vals.len() * 4) as u64
    }

    /// Probes per refuted target: `x(R + P_k)` and `x(R − P_k)` for every `k`.
    pub fn probes_per_target(&self) -> u64 {
        2 * self.base.len() as u64
    }

    #[inline(always)]
    fn insert(&mut self, key: u64, i: usize, j: usize) {
        let mut p = slot(key, self.shift);
        while self.keys[p] != EMPTY {
            p = (p + 1) & self.mask;
        }
        self.keys[p] = key;
        self.vals[p] = ((i << 16) | j) as u32;
        self.entries += 1;
    }

    fn build(&mut self) {
        #[cfg(target_arch = "x86_64")]
        if self.gf.folds() && std::arch::is_x86_feature_detected!("pclmulqdq") {
            // SAFETY: pclmulqdq was just detected.
            return unsafe { self.build_pclmul() };
        }
        self.build_body();
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "pclmulqdq")]
    unsafe fn build_pclmul(&mut self) {
        self.build_body()
    }

    #[inline(always)]
    fn build_body(&mut self) {
        let gf = self.gf;
        let m = self.base.len();
        let mut den = vec![0u64; m];
        let mut scratch: Vec<u64> = Vec::with_capacity(m);
        for i in 0..m {
            let (xi, yi) = self.base[i];
            let cnt = m - i;
            for k in 0..cnt {
                den[k] = xi ^ self.base[i + k].0;
            }
            gf.batch_inv(&mut den[..cnt], &mut scratch);
            // j = i: the double (P − P is the point at infinity: no key).
            if let Pt::Aff(x3, _) = add(gf, Pt::Aff(xi, yi), Pt::Aff(xi, yi)) {
                self.insert(x3, i, i);
            }
            for k in 1..cnt {
                let (xj, yj) = self.base[i + k];
                let inv = den[k];
                let lam = gf.mul(yi ^ yj, inv);
                let sum = gf.sqr(lam) ^ lam ^ xi ^ xj ^ 1;
                let lam2 = lam ^ gf.mul(xj, inv);
                let dif = gf.sqr(lam2) ^ lam2 ^ xi ^ xj ^ 1;
                self.insert(sum, i, i + k);
                self.insert(dif, i, i + k);
            }
        }
    }

    /// Whether `R` decomposes over the base: the first witness, or every sorted root when `all`.
    /// Returns `(found, probes)`.
    pub fn run(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, u64) {
        #[cfg(target_arch = "x86_64")]
        if self.gf.folds() && std::arch::is_x86_feature_detected!("pclmulqdq") {
            // SAFETY: pclmulqdq was just detected.
            return unsafe { self.run_pclmul(r, all, out) };
        }
        self.run_body(r, all, out)
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "pclmulqdq")]
    unsafe fn run_pclmul(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, u64) {
        self.run_body(r, all, out)
    }

    #[inline(always)]
    fn run_body(&self, r: (u64, u64), all: bool, out: &mut Vec<[u64; 3]>) -> (bool, u64) {
        let gf = self.gf;
        let m = self.base.len();
        let (xr, yr) = r;
        let mut den = vec![0u64; m];
        let mut scratch: Vec<u64> = Vec::with_capacity(m);
        for k in 0..m {
            den[k] = xr ^ self.base[k].0;
        }
        // x_R ∈ V is excluded by the caller, so no denominator is zero.
        gf.batch_inv(&mut den, &mut scratch);
        let mut found = false;
        let mut probes = 0u64;
        for k in 0..m {
            let (xk, yk) = self.base[k];
            let inv = den[k];
            if inv == 0 {
                continue;
            }
            let lam = gf.mul(yr ^ yk, inv);
            let plus = gf.sqr(lam) ^ lam ^ xr ^ xk ^ 1;
            let lam2 = lam ^ gf.mul(xk, inv);
            let minus = gf.sqr(lam2) ^ lam2 ^ xr ^ xk ^ 1;
            for key in [plus, minus] {
                probes += 1;
                let mut p = slot(key, self.shift);
                while self.keys[p] != EMPTY {
                    if self.keys[p] == key {
                        let v = self.vals[p] as usize;
                        if self.confirm(v >> 16, v & 0xFFFF, k, r, out) {
                            found = true;
                            if !all {
                                return (true, probes);
                            }
                        }
                    }
                    p = (p + 1) & self.mask;
                }
            }
        }
        (found, probes)
    }

    /// Exact confirmation of a table hit: `x(P_k)` must be one of the four thirds of `(P_i, P_j)`.
    #[cold]
    fn confirm(&self, i: usize, j: usize, k: usize, r: (u64, u64), out: &mut Vec<[u64; 3]>) -> bool {
        let xk = self.base[k].0;
        let ok = thirds(self.gf, self.base[i], self.base[j], r)
            .into_iter()
            .flatten()
            .any(|t| t == xk);
        assert!(ok, "table hit not confirmed by the exact reference");
        let mut t = [self.base[i].0, self.base[j].0, xk];
        t.sort_unstable();
        out.push(t);
        true
    }
}
