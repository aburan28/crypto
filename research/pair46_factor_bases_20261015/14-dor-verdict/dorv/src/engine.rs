//! Port of the m = 2 path of crypto/sage/factor_base_dor.sage.

// ── CPython random.Random (MT19937 + init_by_array + _randbelow) ──────────
pub struct PyRandom {
    mt: [u32; 624],
    idx: usize,
}

impl PyRandom {
    pub fn new(seed: u64) -> Self {
        let mut key: Vec<u32> = Vec::new();
        let mut s = seed;
        if s == 0 { key.push(0); }
        while s > 0 { key.push((s & 0xffff_ffff) as u32); s >>= 32; }
        let mut r = PyRandom { mt: [0; 624], idx: 624 };
        r.init_genrand(19650218);
        let n = 624usize;
        let mut i = 1usize;
        let mut j = 0usize;
        let kl = key.len();
        let mut k = if n > kl { n } else { kl };
        while k > 0 {
            let prev = r.mt[i - 1];
            r.mt[i] = (r.mt[i] ^ ((prev ^ (prev >> 30)).wrapping_mul(1664525)))
                .wrapping_add(key[j]).wrapping_add(j as u32);
            i += 1; j += 1;
            if i >= n { r.mt[0] = r.mt[n - 1]; i = 1; }
            if j >= kl { j = 0; }
            k -= 1;
        }
        k = n - 1;
        while k > 0 {
            let prev = r.mt[i - 1];
            r.mt[i] = (r.mt[i] ^ ((prev ^ (prev >> 30)).wrapping_mul(1566083941)))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= n { r.mt[0] = r.mt[n - 1]; i = 1; }
            k -= 1;
        }
        r.mt[0] = 0x8000_0000;
        r.idx = 624;
        r
    }
    fn init_genrand(&mut self, s: u32) {
        self.mt[0] = s;
        for i in 1..624 {
            let p = self.mt[i - 1];
            self.mt[i] = 1812433253u32.wrapping_mul(p ^ (p >> 30)).wrapping_add(i as u32);
        }
        self.idx = 624;
    }
    fn gen(&mut self) -> u32 {
        if self.idx >= 624 {
            for k in 0..624 {
                let y = (self.mt[k] & 0x8000_0000) | (self.mt[(k + 1) % 624] & 0x7fff_ffff);
                let mut v = self.mt[(k + 397) % 624] ^ (y >> 1);
                if y & 1 != 0 { v ^= 0x9908_b0df; }
                self.mt[k] = v;
            }
            self.idx = 0;
        }
        let mut y = self.mt[self.idx];
        self.idx += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }
    fn getrandbits(&mut self, k: u32) -> u64 {
        assert!(k >= 1 && k <= 32);
        (self.gen() >> (32 - k)) as u64
    }
    pub fn randbelow(&mut self, n: u64) -> u64 {
        let k = 64 - n.leading_zeros();
        let mut r = self.getrandbits(k);
        while r >= n { r = self.getrandbits(k); }
        r
    }
    pub fn randrange(&mut self, a: u64, b: u64) -> u64 { a + self.randbelow(b - a) }
}

// ── GF(2^n) with Sage's default (Conway) modulus ──────────────────────────
pub fn conway(n: usize) -> u64 {
    let terms: &[u32] = match n {
        7 => &[7, 1, 0],
        9 => &[9, 4, 0],
        11 => &[11, 2, 0],
        13 => &[13, 4, 3, 1, 0],
        15 => &[15, 5, 4, 2, 0],
        17 => &[17, 3, 0],
        19 => &[19, 5, 2, 1, 0],
        21 => &[21, 6, 5, 2, 0],
        23 => &[23, 5, 0],
        _ => panic!("no modulus for n={}", n),
    };
    terms.iter().fold(0u64, |a, &e| a | (1u64 << e))
}

#[derive(Clone)]
pub struct Field {
    pub n: usize,
    pub f: u64,
    pub tensor: Vec<Vec<u64>>,
}

impl Field {
    pub fn new(n: usize) -> Self { Self::with_modulus(n, conway(n)) }
    pub fn with_modulus(n: usize, f: u64) -> Self {
        let mut fl = Field { n, f, tensor: vec![] };
        fl.tensor = (0..n).map(|i| (0..n).map(|j| fl.mul(1 << i, 1 << j)).collect()).collect();
        fl
    }
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        let (mut a, mut b, mut r) = (a, b, 0u64);
        while b != 0 {
            if b & 1 != 0 { r ^= a; }
            b >>= 1;
            a <<= 1;
            if (a >> self.n) & 1 != 0 { a ^= self.f; }
        }
        r
    }
    pub fn pow(&self, a: u64, mut e: u64) -> u64 {
        let (mut b, mut r) = (a, 1u64);
        while e > 0 { if e & 1 != 0 { r = self.mul(r, b); } b = self.mul(b, b); e >>= 1; }
        r
    }
    pub fn sq(&self, a: u64) -> u64 { self.mul(a, a) }
    pub fn inv(&self, a: u64) -> u64 { self.pow(a, (1u64 << self.n) - 2) }
    pub fn trace(&self, a: u64) -> u64 {
        let (mut t, mut c) = (0u64, a);
        for _ in 0..self.n { t ^= c; c = self.sq(c); }
        t // 0 or 1
    }
}

// ── GF(2) linear algebra on packed ints (field-size vectors) ──────────────
/// gf2_reduce / subspace_from_generators exactly as the engine: basis keyed
/// by lowest set bit, returned sorted by pivot.
pub fn subspace(gens: &[u64]) -> Vec<u64> {
    let mut basis: std::collections::BTreeMap<u32, u64> = Default::default();
    for &g in gens {
        let mut v = g;
        loop {
            if v == 0 { break; }
            let piv = v.trailing_zeros();
            match basis.get(&piv) {
                Some(&b) => v ^= b,
                None => { basis.insert(piv, v); break; }
            }
        }
    }
    basis.values().cloned().collect()
}

pub fn rank_u64(gens: &[u64]) -> usize { subspace(gens).len() }

pub fn span(basis: &[u64]) -> Vec<u64> {
    let d = basis.len();
    (0..(1u64 << d)).map(|bits| {
        let mut v = 0;
        for i in 0..d { if (bits >> i) & 1 != 0 { v ^= basis[i]; } }
        v
    }).collect()
}

// ── symbolic multilinear polynomials ──────────────────────────────────────
pub type Poly = Vec<u32>; // sorted set of monomial masks; 0 = constant 1

pub fn padd(a: &Poly, b: &Poly) -> Poly {
    let (mut i, mut j) = (0, 0);
    let mut o = Vec::with_capacity(a.len() + b.len());
    while i < a.len() && j < b.len() {
        if a[i] < b[j] { o.push(a[i]); i += 1; }
        else if a[i] > b[j] { o.push(b[j]); j += 1; }
        else { i += 1; j += 1; }
    }
    o.extend_from_slice(&a[i..]);
    o.extend_from_slice(&b[j..]);
    o
}

pub fn pmul(a: &Poly, b: &Poly) -> Poly {
    let mut v: Vec<u32> = Vec::with_capacity(a.len() * b.len());
    for &x in a { for &y in b { v.push(x | y); } }
    v.sort_unstable();
    let mut o = Vec::with_capacity(v.len());
    let mut i = 0;
    while i < v.len() {
        let mut j = i;
        while j < v.len() && v[j] == v[i] { j += 1; }
        if (j - i) % 2 == 1 { o.push(v[i]); }
        i = j;
    }
    o
}

pub type Fe = Vec<Poly>;

pub fn fe_add(u: &Fe, v: &Fe) -> Fe { u.iter().zip(v).map(|(a, b)| padd(a, b)).collect() }

pub fn fe_mul(f: &Field, u: &Fe, v: &Fe) -> Fe {
    let n = f.n;
    let mut w: Fe = vec![vec![]; n];
    for i in 0..n {
        if u[i].is_empty() { continue; }
        for j in 0..n {
            if v[j].is_empty() { continue; }
            let p = pmul(&u[i], &v[j]);
            if p.is_empty() { continue; }
            let mut bits = f.tensor[i][j];
            let mut k = 0;
            while bits != 0 {
                if bits & 1 != 0 { w[k] = padd(&w[k], &p); }
                bits >>= 1;
                k += 1;
            }
        }
    }
    w
}

pub fn fe_const(n: usize, v: u64) -> Fe {
    (0..n).map(|k| if (v >> k) & 1 != 0 { vec![0u32] } else { vec![] }).collect()
}

pub fn fe_linear(n: usize, basis: &[u64], off: usize) -> Fe {
    let mut out: Fe = vec![vec![]; n];
    for (i, &vi) in basis.iter().enumerate() {
        let m = vec![1u32 << (off + i)];
        for k in 0..n { if (vi >> k) & 1 != 0 { out[k] = padd(&out[k], &m); } }
    }
    out
}

pub fn descent_equations(f: &Field, b1: &[u64], b2: &[u64], t: u64, b: u64) -> Vec<Poly> {
    let n = f.n;
    let x1 = fe_linear(n, b1, 0);
    let x2 = fe_linear(n, b2, b1.len());
    let x3 = fe_const(n, t);
    let bv = fe_const(n, b);
    let s12 = fe_add(&x1, &x2);
    let p12 = fe_mul(f, &x1, &x2);
    let s12sq = fe_mul(f, &s12, &s12);
    let x3sq = fe_mul(f, &x3, &x3);
    let t1 = fe_mul(f, &s12sq, &x3sq);
    let t2 = fe_mul(f, &p12, &x3);
    let t3 = fe_add(&fe_mul(f, &p12, &p12), &bv);
    let s3 = fe_add(&fe_add(&t1, &t2), &t3);
    s3.into_iter().filter(|e| !e.is_empty()).collect()
}

pub fn pdeg(p: &Poly) -> u32 { p.iter().map(|m| m.count_ones()).max().unwrap_or(0) }

// ── Macaulay refutation degree ────────────────────────────────────────────
pub struct Basis {
    words: usize,
    piv: Vec<Option<Vec<u64>>>,
    pub rank: usize,
}

impl Basis {
    pub fn new(cols: usize) -> Self { Basis { words: (cols + 63) / 64, piv: vec![None; cols], rank: 0 } }
    /// Reduce by lowest-set-bit pivots; insert if nonzero. Returns true if inserted.
    pub fn insert(&mut self, mut v: Vec<u64>) -> bool {
        let mut w = 0;
        loop {
            while w < self.words && v[w] == 0 { w += 1; }
            if w == self.words { return false; }
            let p = w * 64 + v[w].trailing_zeros() as usize;
            match &self.piv[p] {
                Some(b) => { for k in w..self.words { v[k] ^= b[k]; } }
                None => { self.piv[p] = Some(v); self.rank += 1; return true; }
            }
        }
    }
    pub fn contains(&self, mut v: Vec<u64>) -> bool {
        let mut w = 0;
        loop {
            while w < self.words && v[w] == 0 { w += 1; }
            if w == self.words { return true; }
            let p = w * 64 + v[w].trailing_zeros() as usize;
            match &self.piv[p] {
                Some(b) => { for k in w..self.words { v[k] ^= b[k]; } }
                None => return false,
            }
        }
    }
}

pub struct Refute {
    pub dstar: Option<u32>,
    pub profile: Vec<(u32, usize)>,
    pub shapes: Vec<(u32, usize, usize)>,
}

pub fn monomials_upto(nvars: usize, d: u32) -> Vec<u32> {
    // engine order: degree 0, then combinations by degree (order irrelevant for rank)
    let mut v: Vec<u32> = (0..(1u32 << nvars)).filter(|m| m.count_ones() <= d).collect();
    v.sort_by_key(|m| (m.count_ones(), *m));
    v
}

pub fn refutation_degree(eqs: &[Poly], nvars: usize, dmax: u32) -> Refute {
    let mut out = Refute { dstar: None, profile: vec![], shapes: vec![] };
    let mut index = vec![u32::MAX; 1usize << nvars];
    for d in 2..=dmax {
        let monos = monomials_upto(nvars, d);
        for (i, &m) in monos.iter().enumerate() { index[m as usize] = i as u32; }
        let cols = monos.len();
        let words = (cols + 63) / 64;
        let mut basis = Basis::new(cols);
        let mut nrows = 0usize;
        for eq in eqs {
            let ed = pdeg(eq);
            for &mu in &monos {
                if mu.count_ones() + ed > d { continue; }
                let mut row = vec![0u64; words];
                let mut nz = false;
                for &mono in eq {
                    let c = index[(mu | mono) as usize] as usize;
                    row[c / 64] ^= 1u64 << (c % 64);
                }
                for x in &row { if *x != 0 { nz = true; break; } }
                if nz { nrows += 1; basis.insert(row); }
            }
        }
        let mut e0 = vec![0u64; words];
        e0[0] = 1; // constant monomial has index 0
        out.profile.push((d, basis.rank));
        out.shapes.push((d, nrows, cols));
        if basis.contains(e0) { out.dstar = Some(d); return out; }
        for &m in &monos { index[m as usize] = u32::MAX; }
    }
    out
}

/// Exact E1 invariants of an equation list: (rank of all equations,
/// rank of linear parts, rank of quadratic parts) over GF(2).
pub fn part_ranks(eqs: &[Poly], nvars: usize) -> (usize, usize, usize) {
    let monos = monomials_upto(nvars, 2);
    let mut index = std::collections::HashMap::new();
    for (i, &m) in monos.iter().enumerate() { index.insert(m, i); }
    let cols = monos.len();
    let words = (cols + 63) / 64;
    let mut ball = Basis::new(cols);
    let mut blin = Basis::new(cols);
    let mut bq = Basis::new(cols);
    for eq in eqs {
        let mut r = vec![0u64; words];
        let mut rl = vec![0u64; words];
        let mut rq = vec![0u64; words];
        for &m in eq {
            let c = index[&m];
            r[c / 64] ^= 1 << (c % 64);
            match m.count_ones() { 1 => rl[c / 64] ^= 1 << (c % 64), 2 => rq[c / 64] ^= 1 << (c % 64), _ => {} }
        }
        ball.insert(r); blin.insert(rl); bq.insert(rq);
    }
    (ball.rank, blin.rank, bq.rank)
}

// ── structures ────────────────────────────────────────────────────────────
pub fn normal_conjugates(f: &Field, seed: u64) -> Vec<u64> {
    let n = f.n;
    let mut found: Vec<Vec<u64>> = vec![];
    for v in 1..(1u64 << n) {
        let mut conj = vec![];
        let mut c = v;
        for _ in 0..n { conj.push(c); c = f.sq(c); }
        if rank_u64(&conj) == n {
            found.push(conj);
            if found.len() >= 64 { break; }
        }
    }
    let mut r = PyRandom::new(seed);
    let k = r.randbelow(found.len() as u64) as usize;
    found[k].clone()
}

pub fn structure_bases(f: &Field, name: &str, dim: usize, seed: u64) -> Option<(Vec<u64>, Vec<u64>)> {
    let n = f.n;
    match name {
        "nb_window" | "nb_break" => {
            if 2 * dim > n { return None; }
            let conj = normal_conjugates(f, seed);
            let v1 = subspace(&conj[..dim]);
            if name == "nb_window" { return Some((v1.clone(), v1)); }
            let v2 = subspace(&conj[dim..2 * dim]);
            Some((v1, v2))
        }
        "standard" => {
            let g: Vec<u64> = (0..dim).map(|i| 1u64 << i).collect();
            let b = subspace(&g);
            Some((b.clone(), b))
        }
        "subfield" => {
            if n % dim != 0 { return None; }
            let w = f.pow(2, ((1u64 << n) - 1) / ((1u64 << dim) - 1));
            let g: Vec<u64> = (0..dim as u64).map(|i| f.pow(w, i)).collect();
            let b = subspace(&g);
            Some((b.clone(), b))
        }
        "random" => {
            let mut r = PyRandom::new(seed);
            for _ in 0..256 {
                let g: Vec<u64> = (0..dim).map(|_| r.randrange(1, 1u64 << n)).collect();
                let b = subspace(&g);
                if b.len() == dim { return Some((b.clone(), b)); }
            }
            None
        }
        _ => None,
    }
}

pub fn targets(n: usize, count: usize, seed: u64) -> Vec<u64> {
    let mut r = PyRandom::new(seed);
    let mut seen = std::collections::HashSet::new();
    let mut out = vec![];
    while out.len() < count && seen.len() < (1usize << n) - 1 {
        let t = r.randrange(1, 1u64 << n);
        if seen.insert(t) { out.push(t); }
    }
    out
}

/// True if target has NO decomposition with x1 in span(b1), x2 in span(b2).
pub fn unsat(f: &Field, m1: &[u64], m2: &[u64], t: u64, b: u64) -> bool {
    let t2 = f.sq(t);
    for &x1 in m1 {
        for &x2 in m2 {
            let s = f.sq(x1 ^ x2);
            let p = f.mul(x1, x2);
            let v = f.mul(s, t2) ^ f.mul(p, t) ^ f.sq(p) ^ b;
            if v == 0 { return false; }
        }
    }
    true
}

/// Points on y^2 + xy = x^3 + 1 with x in span(basis), x != 0 (lift_x(all)).
pub fn fb_points(f: &Field, basis: &[u64]) -> usize {
    let mut c = 0;
    for x in span(basis) {
        if x == 0 { continue; }
        let xi2 = f.sq(f.inv(x));
        if f.trace(x ^ xi2) == 0 { c += 2; }
    }
    c
}

/// dim of F2-span of {v*w : v in V1 basis, w in V2 basis}
pub fn product_dim(f: &Field, b1: &[u64], b2: &[u64]) -> (usize, Vec<u64>) {
    let mut g = vec![];
    for &v in b1 { for &w in b2 { g.push(f.mul(v, w)); } }
    let s = subspace(&g);
    (s.len(), s)
}

pub fn in_span(basis: &[u64], x: u64) -> bool {
    let mut g = basis.to_vec();
    let r = rank_u64(&g);
    g.push(x);
    rank_u64(&g) == r
}
