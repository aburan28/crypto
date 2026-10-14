//! Scalar row kernel (PS-2): pairs-and-solve by the group law with the two
//! row points `U± = P₁ ± R` hoisted out of the pair loop and the pair
//! denominators inverted a row at a time (Montgomery's trick).
//!
//! For a base pair `(P₁, P₂)` the four candidate third abscissae are
//! `x(U₊ ± P₂)` and `x(U₋ ± P₂)`; with `λ = (y_U + y₂)/(x_U + x₂)`:
//!
//! ```text
//!   x(U + P₂) = λ² + λ + (x_U + x₂) + 1
//!   x(U − P₂) = λ̃² + λ̃ + (x_U + x₂) + 1,   λ̃ = λ + x₂/(x_U + x₂)
//! ```
//!
//! so a pair costs two inversions (batched), four multiplications and four
//! squarings, and the test `X₃ ∈ V` is `X₃ < 2^ℓ`.

use crate::group::*;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;

pub struct ScalarKernel<'a> {
    pub gf: &'a Gf2,
    pub l: u32,
    /// Rational base points `(x, y)`, `x` ascending.
    pub base: Vec<(u64, u64)>,
}

impl<'a> ScalarKernel<'a> {
    pub fn new(gf: &'a Gf2, solver: &ArtinSchreier, l: u32) -> Self {
        Self {
            gf,
            l,
            base: base_points(gf, solver, l),
        }
    }

    pub fn pairs_total(&self) -> u64 {
        let m = self.base.len() as u64;
        m * (m + 1) / 2
    }

    /// Whether `R` decomposes over the base: the first witness and the number
    /// of pairs visited.  `all` collects every sorted root instead of
    /// stopping.
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
        let span = 1u64 << self.l;
        let base = &self.base;
        let m = base.len();
        let rp = Pt::Aff(r.0, r.1);
        let nr = neg(rp);
        let mut den = vec![0u64; 2 * m];
        let mut scratch: Vec<u64> = Vec::with_capacity(2 * m);
        let mut pairs = 0u64;
        let mut found = false;
        for i1 in 0..m {
            let p1 = Pt::Aff(base[i1].0, base[i1].1);
            let (Pt::Aff(xu, yu), Pt::Aff(xv, yv)) = (add(gf, p1, rp), add(gf, p1, nr)) else {
                // P₁ = ±R: only when x_R ∈ V, which the caller excludes.
                for i2 in i1..m {
                    pairs += 1;
                    self.exact_pair(i1, i2, r, all, out, &mut found);
                    if found && !all {
                        return (true, pairs);
                    }
                }
                continue;
            };
            let cnt = m - i1;
            for k in 0..cnt {
                let x2 = base[i1 + k].0;
                den[k] = xu ^ x2;
                den[cnt + k] = xv ^ x2;
            }
            gf.batch_inv(&mut den[..2 * cnt], &mut scratch);
            for k in 0..cnt {
                let (x2, y2) = base[i1 + k];
                pairs += 1;
                let (j, jj) = (den[k], den[cnt + k]);
                if j == 0 || jj == 0 {
                    self.exact_pair(i1, i1 + k, r, all, out, &mut found);
                    if found && !all {
                        return (true, pairs);
                    }
                    continue;
                }
                let d = xu ^ x2;
                let lam = gf.mul(yu ^ y2, j);
                let t1 = gf.sqr(lam) ^ lam ^ d ^ 1;
                let lam2 = lam ^ gf.mul(x2, j);
                let t2 = gf.sqr(lam2) ^ lam2 ^ d ^ 1;
                let dd = xv ^ x2;
                let mu = gf.mul(yv ^ y2, jj);
                let t3 = gf.sqr(mu) ^ mu ^ dd ^ 1;
                let mu2 = mu ^ gf.mul(x2, jj);
                let t4 = gf.sqr(mu2) ^ mu2 ^ dd ^ 1;
                if (t1 < span) | (t2 < span) | (t3 < span) | (t4 < span) {
                    found = true;
                    if !all {
                        let t = [t1, t2, t3, t4].into_iter().find(|&t| t < span).unwrap();
                        out.push(sorted(base[i1].0, x2, t));
                        return (true, pairs);
                    }
                    for t in [t1, t2, t3, t4] {
                        if t < span {
                            out.push(sorted(base[i1].0, x2, t));
                        }
                    }
                }
            }
        }
        (found, pairs)
    }

    #[cold]
    fn exact_pair(
        &self,
        i1: usize,
        i2: usize,
        r: (u64, u64),
        all: bool,
        out: &mut Vec<[u64; 3]>,
        found: &mut bool,
    ) {
        let span = 1u64 << self.l;
        for t in thirds(self.gf, self.base[i1], self.base[i2], r)
            .into_iter()
            .flatten()
        {
            if t < span {
                *found = true;
                out.push(sorted(self.base[i1].0, self.base[i2].0, t));
                if !all {
                    return;
                }
            }
        }
    }
}

fn sorted(a: u64, b: u64, c: u64) -> [u64; 3] {
    let mut t = [a, b, c];
    t.sort_unstable();
    t
}
