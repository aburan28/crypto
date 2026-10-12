//! The subspace-leaf control: the PS-2 scalar group-law kernel, copied from
//! `examples/pdp_speedup/{group,kernel}.rs`, run on the PS-0 verdict corpus
//! with the library's field-operation counters (`--features pdp-opcount`).
//!
//! The curve of the corpus is `y² + xy = x³ + x² + 1` (`a = b = 1`) in a
//! polynomial basis, the base is `V = ⟨1, z, …, z^{ℓ−1}⟩` and the leaf is the
//! comparison `x₃ < 2^ℓ`.  The weight-base kernel (`kernel.rs`) has the same
//! row structure; what differs is the field, the curve constant `a`, and the
//! leaf.

use crypto_lib::binary_ecc::IrreduciblePoly;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;

/// Solver for `z² + z = c` by Gaussian elimination on the `F₂`-linear map `z ↦ z² + z`.
pub struct ArtinSchreier {
    rows: Vec<Option<(u64, u64)>>,
}

impl ArtinSchreier {
    pub fn new(gf: &Gf2) -> Self {
        let n = gf.n;
        let mut rows: Vec<Option<(u64, u64)>> = vec![None; n as usize];
        for i in 0..n {
            let e = 1u64 << i;
            let (mut v, mut w) = (gf.sqr(e) ^ e, e);
            while v != 0 {
                let b = 63 - v.leading_zeros() as usize;
                match rows[b] {
                    Some((pv, pw)) => {
                        v ^= pv;
                        w ^= pw;
                    }
                    None => {
                        rows[b] = Some((v, w));
                        break;
                    }
                }
            }
        }
        Self { rows }
    }

    pub fn solve(&self, c: u64) -> Option<u64> {
        let (mut v, mut w) = (c, 0u64);
        while v != 0 {
            let b = 63 - v.leading_zeros() as usize;
            let (pv, pw) = self.rows[b]?;
            v ^= pv;
            w ^= pw;
        }
        Some(w)
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CPt {
    Inf,
    Aff(u64, u64),
}

fn cneg(p: CPt) -> CPt {
    match p {
        CPt::Inf => CPt::Inf,
        CPt::Aff(x, y) => CPt::Aff(x, y ^ x),
    }
}

fn cadd(gf: &Gf2, p: CPt, q: CPt) -> CPt {
    match (p, q) {
        (CPt::Inf, r) | (r, CPt::Inf) => r,
        (CPt::Aff(x1, y1), CPt::Aff(x2, y2)) => {
            if x1 == x2 {
                if y1 ^ y2 == x1 {
                    return CPt::Inf;
                }
                let lam = x1 ^ gf.mul(y1, gf.inv(x1));
                let x3 = gf.sqr(lam) ^ lam ^ 1;
                let y3 = gf.sqr(x1) ^ gf.mul(lam ^ 1, x3);
                return CPt::Aff(x3, y3);
            }
            let lam = gf.mul(y1 ^ y2, gf.inv(x1 ^ x2));
            let x3 = gf.sqr(lam) ^ lam ^ x1 ^ x2 ^ 1;
            let y3 = gf.mul(lam, x1 ^ x3) ^ x3 ^ y1;
            CPt::Aff(x3, y3)
        }
    }
}

fn clift(gf: &Gf2, solver: &ArtinSchreier, x: u64) -> Option<u64> {
    if x == 0 {
        return Some(1);
    }
    let c = x ^ 1 ^ gf.sqr(gf.inv(x));
    let z = solver.solve(c)?;
    let y = gf.mul(x, z);
    Some(y.min(y ^ x))
}

fn cthirds(gf: &Gf2, p1: (u64, u64), p2: (u64, u64), r: (u64, u64)) -> [Option<u64>; 4] {
    let (p1, p2) = (CPt::Aff(p1.0, p1.1), CPt::Aff(p2.0, p2.1));
    let r = CPt::Aff(r.0, r.1);
    let mut out = [None; 4];
    for (k, z) in [cadd(gf, p1, p2), cadd(gf, p1, cneg(p2))]
        .into_iter()
        .enumerate()
    {
        for (m, rr) in [r, cneg(r)].into_iter().enumerate() {
            if let CPt::Aff(x, _) = cadd(gf, z, rr) {
                out[2 * k + m] = Some(x);
            }
        }
    }
    out
}

pub struct Control<'a> {
    pub gf: &'a Gf2,
    pub l: u32,
    pub base: Vec<(u64, u64)>,
}

pub struct CtlRun {
    pub pairs: u64,
    pub found: bool,
    pub hits: u64,
}

impl<'a> Control<'a> {
    pub fn new(gf: &'a Gf2, l: u32) -> Self {
        let solver = ArtinSchreier::new(gf);
        let base = (0..(1u64 << l))
            .filter_map(|x| clift(gf, &solver, x).map(|y| (x, y)))
            .collect();
        Self { gf, l, base }
    }

    /// Full scan over the pairs `i₁ ≤ i₂`.
    pub fn run(&self, r: (u64, u64)) -> CtlRun {
        #[cfg(target_arch = "x86_64")]
        if self.gf.folds() && std::arch::is_x86_feature_detected!("pclmulqdq") {
            // SAFETY: pclmulqdq was just detected.
            return unsafe { self.run_pclmul(r) };
        }
        self.run_body(r)
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "pclmulqdq")]
    unsafe fn run_pclmul(&self, r: (u64, u64)) -> CtlRun {
        self.run_body(r)
    }

    #[inline(always)]
    fn run_body(&self, r: (u64, u64)) -> CtlRun {
        let gf = self.gf;
        let span = 1u64 << self.l;
        let base = &self.base;
        let m = base.len();
        let rp = CPt::Aff(r.0, r.1);
        let nr = cneg(rp);
        let mut den = vec![0u64; 2 * m];
        let mut scratch: Vec<u64> = Vec::with_capacity(2 * m);
        let mut res = CtlRun {
            pairs: 0,
            found: false,
            hits: 0,
        };
        for i1 in 0..m {
            let p1 = CPt::Aff(base[i1].0, base[i1].1);
            let (CPt::Aff(xu, yu), CPt::Aff(xv, yv)) = (cadd(gf, p1, rp), cadd(gf, p1, nr)) else {
                for i2 in i1..m {
                    res.pairs += 1;
                    for t in cthirds(gf, base[i1], base[i2], r).into_iter().flatten() {
                        if t < span {
                            res.found = true;
                            res.hits += 1;
                        }
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
                res.pairs += 1;
                let (j, jj) = (den[k], den[cnt + k]);
                if j == 0 || jj == 0 {
                    for t in cthirds(gf, base[i1], base[i1 + k], r).into_iter().flatten() {
                        if t < span {
                            res.found = true;
                            res.hits += 1;
                        }
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
                for t in [t1, t2, t3, t4] {
                    if t < span {
                        res.found = true;
                        res.hits += 1;
                    }
                }
            }
        }
        res
    }
}

pub struct CorpusTarget {
    pub planted: bool,
    pub xr: u64,
    pub yr: u64,
    pub has_rational_root: bool,
}

pub struct CorpusFile {
    pub l: u32,
    pub n: u32,
    pub modulus: u64,
    pub targets: Vec<CorpusTarget>,
}

pub fn read_corpus(text: &str) -> CorpusFile {
    let h = |s: &str| u64::from_str_radix(s.trim_start_matches("0x"), 16).unwrap();
    let mut c = CorpusFile {
        l: 0,
        n: 0,
        modulus: 0,
        targets: Vec::new(),
    };
    for line in text.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f.first().copied() {
            Some("l") => c.l = f[1].parse().unwrap(),
            Some("n") => c.n = f[1].parse().unwrap(),
            Some("modulus") => c.modulus = h(f[1]),
            Some("T") => c.targets.push(CorpusTarget {
                planted: f[2] == "p",
                xr: h(f[3]),
                yr: h(f[4]),
                has_rational_root: false,
            }),
            Some("R") if f[4] == "1" => {
                c.targets.last_mut().unwrap().has_rational_root = true;
            }
            _ => {}
        }
    }
    c
}

pub fn field_from_modulus(n: u32, modulus: u64) -> Gf2 {
    let low_terms: Vec<u32> = (0..n).filter(|&b| (modulus >> b) & 1 == 1).collect();
    Gf2::new(&IrreduciblePoly {
        degree: n,
        low_terms,
    })
}
