//! The quaternion algebra B_{p,inf} = (-1, -p | Q) for p = 3 mod 4 (i^2 = -1, j^2 = -p, k = ij),
//! its maximal order O_0 = <1, i, (i+j)/2, (1+k)/2> (End(E_0) for E_0: y^2 = x^3 + x, with i the
//! automorphism (x, y) -> (-x, sqrt(-1) y) and j the Frobenius), lattices in Hermite normal form,
//! ideals, norms, left/right orders, LLL on the norm form and short-vector enumeration. Built on
//! `int::Int` (exact arithmetic throughout). Submodules: `klpt` (KLPT for O_0), `brandt` (class
//! sets, Brandt matrices, the supersingular j-graph), `deuring` (ideals to isogenies on E_0).
pub mod brandt;
pub mod deuring;
pub mod klpt;

use crate::int::Int;

fn int(v: i64) -> Int {
    Int::from(v)
}

/// x = (c0 + c1 i + c2 j + c3 k) / d, d > 0, gcd(c, d) = 1.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Quat {
    pub c: [Int; 4],
    pub d: Int,
}

impl Quat {
    pub fn new(c: [Int; 4], d: Int) -> Quat {
        assert!(!d.is_zero());
        let (mut c, mut d) = (c, d);
        if d.is_neg() {
            c = c.map(|x| -x);
            d = -d;
        }
        let mut g = d.clone();
        for x in &c {
            g = g.gcd(x);
        }
        if g != Int::one() {
            c = c.map(|x| x.div_floor(&g));
            d = d.div_floor(&g);
        }
        Quat { c, d }
    }
    pub fn from_ints(c: [i64; 4]) -> Quat {
        Quat::new(c.map(int), Int::one())
    }
    pub fn scalar(v: &Int) -> Quat {
        Quat::new([v.clone(), Int::zero(), Int::zero(), Int::zero()], Int::one())
    }
    pub fn zero() -> Quat {
        Quat::from_ints([0, 0, 0, 0])
    }
    pub fn one() -> Quat {
        Quat::from_ints([1, 0, 0, 0])
    }
    pub fn is_zero(&self) -> bool {
        self.c.iter().all(|x| x.is_zero())
    }
    pub fn conj(&self) -> Quat {
        Quat::new(
            [self.c[0].clone(), -&self.c[1], -&self.c[2], -&self.c[3]],
            self.d.clone(),
        )
    }
    pub fn scale(&self, num: &Int, den: &Int) -> Quat {
        Quat::new(self.c.clone().map(|x| &x * num), &self.d * den)
    }
    pub fn add(&self, o: &Quat) -> Quat {
        Quat::new(
            std::array::from_fn(|t| &(&self.c[t] * &o.d) + &(&o.c[t] * &self.d)),
            &self.d * &o.d,
        )
    }
    pub fn sub(&self, o: &Quat) -> Quat {
        self.add(&o.scale(&int(-1), &Int::one()))
    }
}

/// B_{p,inf} with i^2 = -1, j^2 = -p (p = 3 mod 4).
#[derive(Clone, Debug)]
pub struct Alg {
    pub p: Int,
}

impl Alg {
    pub fn new(p: &Int) -> Alg {
        assert!(p.mod_u64(4) == 3, "p must be 3 mod 4");
        Alg { p: p.clone() }
    }
    /// Product of integer coordinate vectors (no denominators).
    pub fn mul_raw(&self, x: &[Int; 4], y: &[Int; 4]) -> [Int; 4] {
        let p = &self.p;
        let [a1, b1, c1, d1] = x;
        let [a2, b2, c2, d2] = y;
        // q = 1: real a1a2 - b1b2 - p c1c2 - p d1d2; i a1b2 + b1a2 + p(c1d2 - d1c2);
        // j a1c2 + c1a2 + (d1b2 - b1d2); k a1d2 + d1a2 + b1c2 - c1b2
        [
            &(&(a1 * a2) - &(b1 * b2)) - &(p * &(&(c1 * c2) + &(d1 * d2))),
            &(&(a1 * b2) + &(b1 * a2)) + &(p * &(&(c1 * d2) - &(d1 * c2))),
            &(&(a1 * c2) + &(c1 * a2)) + &(&(d1 * b2) - &(b1 * d2)),
            &(&(a1 * d2) + &(d1 * a2)) + &(&(b1 * c2) - &(c1 * b2)),
        ]
    }
    pub fn mul(&self, x: &Quat, y: &Quat) -> Quat {
        Quat::new(self.mul_raw(&x.c, &y.c), &x.d * &y.d)
    }
    /// Reduced norm as a fraction (num, den), den > 0.
    pub fn nrd(&self, x: &Quat) -> (Int, Int) {
        let n = self.qf(&x.c);
        let d = &x.d * &x.d;
        let g = n.gcd(&d);
        (n.div_floor(&g), d.div_floor(&g))
    }
    /// a^2 + b^2 + p c^2 + p d^2 on integer coordinates.
    pub fn qf(&self, c: &[Int; 4]) -> Int {
        &(&(&c[0] * &c[0]) + &(&c[1] * &c[1])) + &(&self.p * &(&(&c[2] * &c[2]) + &(&c[3] * &c[3])))
    }
    /// Bilinear form <u, v> = (qf(u+v) - qf(u) - qf(v))/2 = u0v0 + u1v1 + p(u2v2 + u3v3).
    pub fn bil(&self, u: &[Int; 4], v: &[Int; 4]) -> Int {
        &(&(&u[0] * &v[0]) + &(&u[1] * &v[1])) + &(&self.p * &(&(&u[2] * &v[2]) + &(&u[3] * &v[3])))
    }
    /// The maximal order O_0 = <1, i, (i+j)/2, (1+k)/2>.
    pub fn o0(&self) -> Lattice {
        Lattice::from_gens(&[
            Quat::from_ints([1, 0, 0, 0]),
            Quat::from_ints([0, 1, 0, 0]),
            Quat::new([int(0), int(1), int(1), int(0)], int(2)),
            Quat::new([int(1), int(0), int(0), int(1)], int(2)),
        ])
    }
    pub fn i(&self) -> Quat {
        Quat::from_ints([0, 1, 0, 0])
    }
    pub fn j(&self) -> Quat {
        Quat::from_ints([0, 0, 1, 0])
    }
    pub fn k(&self) -> Quat {
        Quat::from_ints([0, 0, 0, 1])
    }
}

/// Full-rank lattice (1/den) * rowspace(m), m in upper-triangular Hermite normal form
/// (positive pivots m[r][r], entries above each pivot reduced into [0, pivot)).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Lattice {
    pub m: [[Int; 4]; 4],
    pub den: Int,
}

/// Row-style HNF of integer rows (at least four, full rank 4).
fn hnf(rows: Vec<[Int; 4]>) -> [[Int; 4]; 4] {
    let mut rows: Vec<[Int; 4]> = rows.into_iter().filter(|r| r.iter().any(|x| !x.is_zero())).collect();
    let mut out: Vec<[Int; 4]> = vec![];
    for col in 0..4 {
        // combine all rows' entries in this column into one pivot row by gcd steps
        loop {
            let nz: Vec<usize> = (0..rows.len()).filter(|&r| !rows[r][col].is_zero()).collect();
            if nz.len() <= 1 {
                break;
            }
            // pick the row with smallest |entry| and reduce the others by it
            let piv = *nz.iter().min_by(|&&a, &&b| rows[a][col].abs().cmp(&rows[b][col].abs())).unwrap();
            let pv = rows[piv][col].clone();
            for &r in &nz {
                if r == piv {
                    continue;
                }
                let q = rows[r][col].div_floor(&pv);
                let pr = rows[piv].clone();
                for t in 0..4 {
                    rows[r][t] = &rows[r][t] - &(&q * &pr[t]);
                }
            }
            rows.retain(|r| r.iter().any(|x| !x.is_zero()));
        }
        let Some(pos) = rows.iter().position(|r| !r[col].is_zero()) else {
            panic!("HNF: lattice not of full rank");
        };
        let mut pr = rows.swap_remove(pos);
        if pr[col].is_neg() {
            pr = pr.map(|x| -x);
        }
        out.push(pr);
    }
    let mut m: [[Int; 4]; 4] = std::array::from_fn(|r| out[r].clone());
    // reduce above the pivots
    for c in 0..4 {
        let pv = m[c][c].clone();
        for r in 0..c {
            let q = m[r][c].div_floor(&pv);
            if !q.is_zero() {
                let pr = m[c].clone();
                for t in 0..4 {
                    m[r][t] = &m[r][t] - &(&q * &pr[t]);
                }
            }
        }
    }
    m
}

impl Lattice {
    pub fn from_gens(gens: &[Quat]) -> Lattice {
        let mut den = Int::one();
        for g in gens {
            den = (&den * &g.d).div_floor(&den.gcd(&g.d));
        }
        let rows: Vec<[Int; 4]> = gens
            .iter()
            .map(|g| {
                let f = den.div_floor(&g.d);
                g.c.clone().map(|x| &x * &f)
            })
            .collect();
        Lattice::from_int_rows(rows, den)
    }
    fn from_int_rows(rows: Vec<[Int; 4]>, den: Int) -> Lattice {
        let m = hnf(rows);
        // remove common content
        let mut g = den.clone();
        for r in &m {
            for x in r {
                g = g.gcd(x);
            }
        }
        if g != Int::one() {
            let m2 = m.map(|r| r.map(|x| x.div_floor(&g)));
            Lattice { m: m2, den: den.div_floor(&g) }
        } else {
            Lattice { m, den }
        }
    }
    pub fn basis(&self) -> Vec<Quat> {
        self.m.iter().map(|r| Quat::new(r.clone(), self.den.clone())).collect()
    }
    /// Is x in the lattice? (back substitution in the triangular basis)
    pub fn contains(&self, x: &Quat) -> bool {
        // x * den = sum u_r m[r] with u integral; x * den = c * den / x.d
        let mut v: [Int; 4] = x.c.clone().map(|c| &c * &self.den);
        // need v / x.d integral combination
        for col in 0..4 {
            // coefficient of row col: v[col] / (x.d * m[col][col])
            let piv = &self.m[col][col] * &x.d;
            let (u, r) = v[col].divrem_trunc(&piv);
            if !r.is_zero() {
                return false;
            }
            for t in col..4 {
                v[t] = &v[t] - &(&(&u * &self.m[col][t]) * &x.d);
            }
        }
        true
    }
    pub fn contains_lattice(&self, o: &Lattice) -> bool {
        o.basis().iter().all(|b| self.contains(b))
    }
    pub fn mul(&self, alg: &Alg, o: &Lattice) -> Lattice {
        let mut rows = vec![];
        for a in &self.m {
            for b in &o.m {
                rows.push(alg.mul_raw(a, b));
            }
        }
        Lattice::from_int_rows(rows, &self.den * &o.den)
    }
    /// x * L (left multiplication by an element).
    pub fn lmul(&self, alg: &Alg, x: &Quat) -> Lattice {
        let rows = self.m.iter().map(|b| alg.mul_raw(&x.c, b)).collect();
        Lattice::from_int_rows(rows, &self.den * &x.d)
    }
    /// L * x.
    pub fn rmul(&self, alg: &Alg, x: &Quat) -> Lattice {
        let rows = self.m.iter().map(|b| alg.mul_raw(b, &x.c)).collect();
        Lattice::from_int_rows(rows, &self.den * &x.d)
    }
    pub fn conj(&self) -> Lattice {
        Lattice::from_gens(&self.basis().iter().map(|b| b.conj()).collect::<Vec<_>>())
    }
    pub fn scale(&self, num: &Int, den: &Int) -> Lattice {
        Lattice::from_int_rows(self.m.iter().map(|r| r.clone().map(|x| &x * num)).collect(), &self.den * den)
    }
    /// Covolume as a fraction (num, den) (determinant of the basis).
    pub fn det(&self) -> (Int, Int) {
        let mut n = Int::one();
        for r in 0..4 {
            n = &n * &self.m[r][r];
        }
        let d = self.den.pow(4);
        let g = n.gcd(&d);
        (n.div_floor(&g), d.div_floor(&g))
    }
    /// Reduced norm of an ideal of the order `o` (N(I)^2 = det(I)/det(O)), as a fraction.
    pub fn norm(&self, o: &Lattice) -> (Int, Int) {
        let (a, b) = self.det();
        let (c, d) = o.det();
        let (n, m) = (&a * &d, &b * &c);
        let g = n.gcd(&m);
        let (n, m) = (n.div_floor(&g), m.div_floor(&g));
        let (rn, rm) = (n.isqrt(), m.isqrt());
        assert!(&rn * &rn == n && &rm * &rm == m, "not an ideal: index is not a square");
        (rn, rm)
    }
    /// Integral norm for integral ideals (panics if not integral).
    pub fn norm_int(&self, o: &Lattice) -> Int {
        let (n, d) = self.norm(o);
        assert!(d == Int::one());
        n
    }
    /// Left order I Ibar / N(I) and right order Ibar I / N(I) (I invertible).
    pub fn left_order(&self, alg: &Alg, o: &Lattice) -> Lattice {
        let (n, d) = self.norm(o);
        self.mul(alg, &self.conj()).scale(&d, &n)
    }
    pub fn right_order(&self, alg: &Alg, o: &Lattice) -> Lattice {
        let (n, d) = self.norm(o);
        self.conj().mul(alg, self).scale(&d, &n)
    }
    /// Integer Gram matrix of the norm form on the basis, scaled by den^2:
    /// Nrd(sum u_r b_r) = u^T G u / den^2.
    pub fn gram(&self, alg: &Alg) -> [[Int; 4]; 4] {
        std::array::from_fn(|r| std::array::from_fn(|s| alg.bil(&self.m[r], &self.m[s])))
    }
    /// LLL-reduced basis (as integer rows over `den`) w.r.t. the norm form.
    pub fn reduced_basis(&self, alg: &Alg) -> Vec<[Int; 4]> {
        let h = lll_gram(&self.gram(alg));
        (0..4)
            .map(|r| {
                std::array::from_fn(|t| {
                    let mut s = Int::zero();
                    for k in 0..4 {
                        s = &s + &(&h[r][k] * &self.m[k][t]);
                    }
                    s
                })
            })
            .collect()
    }
    /// Elements with Nrd(x) <= bound (bound as a fraction num/den of the lattice scale), by
    /// Fincke-Pohst enumeration on the LLL-reduced basis (f64 bounds, exact recheck).
    /// Returns (x, Nrd(x) * lattice den^2) pairs, one of each +-x, excluding 0; at most `cap`.
    pub fn short_vectors(&self, alg: &Alg, bound_scaled: &Int, cap: usize) -> Vec<([Int; 4], Int)> {
        let rb = self.reduced_basis(alg);
        let g: Vec<Vec<f64>> = (0..4).map(|r| (0..4).map(|s| alg.bil(&rb[r], &rb[s]).to_f64()).collect()).collect();
        let bound = bound_scaled.to_f64() * (1.0 + 1e-9) + 1.0;
        // Cholesky-like decomposition q_ii, q_ij (Fincke-Pohst)
        let n = 4;
        let mut q = vec![vec![0f64; n]; n];
        let mut a = g.clone();
        for i in 0..n {
            q[i][i] = a[i][i];
            for j in i + 1..n {
                q[i][j] = a[i][j] / a[i][i];
            }
            for k in i + 1..n {
                for l in k..n {
                    a[k][l] -= q[i][k] * q[i][l] * q[i][i];
                }
            }
        }
        let mut out = vec![];
        let mut x = vec![0i64; n];
        fn rec(
            i: isize,
            x: &mut Vec<i64>,
            rem: f64,
            q: &Vec<Vec<f64>>,
            rb: &Vec<[Int; 4]>,
            alg: &Alg,
            bound: &Int,
            out: &mut Vec<([Int; 4], Int)>,
            cap: usize,
        ) {
            if out.len() >= cap {
                return;
            }
            if i < 0 {
                if x.iter().all(|&v| v == 0) {
                    return;
                }
                let v: [Int; 4] = std::array::from_fn(|t| {
                    let mut s = Int::zero();
                    for k in 0..4 {
                        s = &s + &(&Int::from(x[k]) * &rb[k][t]);
                    }
                    s
                });
                let nv = alg.qf(&v);
                if nv <= *bound {
                    out.push((v, nv));
                }
                return;
            }
            let iu = i as usize;
            let mut c = 0f64;
            for j in iu + 1..4 {
                c += q[iu][j] * x[j] as f64;
            }
            let r = (rem / q[iu][iu]).max(0.0).sqrt();
            let mut lo = (-c - r).ceil() as i64;
            // one of +-x: while every higher coordinate is zero, this one is >= 0
            if x[iu + 1..].iter().all(|&v| v == 0) {
                lo = lo.max(0);
            }
            let hi = (-c + r).floor() as i64;
            for v in lo..=hi {
                x[iu] = v;
                let t = v as f64 + c;
                let nrem = rem - q[iu][iu] * t * t;
                if nrem < -1e-6 * rem.abs().max(1.0) {
                    continue;
                }
                rec(i - 1, x, nrem, q, rb, alg, bound, out, cap);
                if out.len() >= cap {
                    return;
                }
            }
            x[iu] = 0;
        }
        rec(3, &mut x, bound, &q, &rb, alg, bound_scaled, &mut out, cap);
        out
    }
}

/// Integral LLL (Cohen, Alg. 2.6.7) on a positive-definite integer Gram matrix; returns the
/// unimodular transformation H (rows: new basis in terms of the old), delta = 3/4.
pub fn lll_gram(g: &[[Int; 4]; 4]) -> [[Int; 4]; 4] {
    let n = 4usize;
    let mut b: Vec<Vec<Int>> = g.iter().map(|r| r.to_vec()).collect(); // current Gram
    let mut h: Vec<Vec<Int>> = (0..n).map(|i| (0..n).map(|j| Int::from((i == j) as i64)).collect()).collect();
    let mut d = vec![Int::one(); n + 1]; // d[0] = 1, d[i+1] = det of leading (i+1) block
    let mut lam = vec![vec![Int::zero(); n]; n];
    // incremental Gram-Schmidt in integral form
    let gs = |b: &Vec<Vec<Int>>, d: &mut Vec<Int>, lam: &mut Vec<Vec<Int>>, k: usize| {
        for j in 0..=k {
            let mut u = b[k][j].clone();
            for i in 0..j {
                u = (&(&d[i + 1] * &u) - &(&lam[k][i] * &lam[j][i])).div_floor(&d[i]);
            }
            if j < k {
                lam[k][j] = u;
            } else {
                d[k + 1] = u;
            }
        }
    };
    for k in 0..n {
        gs(&b, &mut d, &mut lam, k);
    }
    let mut k = 1usize;
    let swap = |k: usize, b: &mut Vec<Vec<Int>>, h: &mut Vec<Vec<Int>>, d: &mut Vec<Int>, lam: &mut Vec<Vec<Int>>| {
        h.swap(k, k - 1);
        b.swap(k, k - 1);
        for row in b.iter_mut() {
            row.swap(k, k - 1);
        }
        for j in 0..k - 1 {
            let t = lam[k][j].clone();
            lam[k][j] = lam[k - 1][j].clone();
            lam[k - 1][j] = t;
        }
        let l = lam[k][k - 1].clone();
        let bb = (&(&d[k - 1] * &d[k + 1]) + &(&l * &l)).div_floor(&d[k]);
        for i in k + 1..n {
            let t = lam[i][k].clone();
            lam[i][k] = (&(&d[k + 1] * &lam[i][k - 1]) - &(&l * &t)).div_floor(&d[k]);
            lam[i][k - 1] = (&(&bb * &t) + &(&l * &lam[i][k])).div_floor(&d[k + 1]);
        }
        d[k] = bb;
    };
    let red = |k: usize, l: usize, b: &mut Vec<Vec<Int>>, h: &mut Vec<Vec<Int>>, d: &Vec<Int>, lam: &mut Vec<Vec<Int>>| {
        let two = Int::from(2i64);
        if (&lam[k][l] * &two).abs() > d[l + 1] {
            let q = lam[k][l].div_round(&d[l + 1]);
            // b_k -= q b_l (on H and on the Gram matrix)
            for t in 0..n {
                h[k][t] = &h[k][t] - &(&q * &h[l][t]);
            }
            // Gram: row/col k
            let bkl = b[k][l].clone();
            let bll = b[l][l].clone();
            let bkk = &(&b[k][k] - &(&(&q * &bkl) * &two)) + &(&(&q * &q) * &bll);
            for t in 0..n {
                if t != k {
                    let v = &b[k][t] - &(&q * &b[l][t]);
                    b[k][t] = v.clone();
                    b[t][k] = v;
                }
            }
            b[k][k] = bkk;
            lam[k][l] = &lam[k][l] - &(&q * &d[l + 1]);
            for i in 0..l {
                lam[k][i] = &lam[k][i] - &(&q * &lam[l][i]);
            }
        }
    };
    while k < n {
        red(k, k - 1, &mut b, &mut h, &d, &mut lam);
        // Lovasz: 4 d_{k+1} d_{k-1} < 3 d_k^2 - 4 lam^2  -> swap
        let lhs = &(&int(4) * &d[k + 1]) * &d[k - 1];
        let rhs = &(&int(3) * &(&d[k] * &d[k])) - &(&int(4) * &(&lam[k][k - 1] * &lam[k][k - 1]));
        if lhs < rhs {
            swap(k, &mut b, &mut h, &mut d, &mut lam);
            if k > 1 {
                k -= 1;
            }
        } else {
            for l in (0..k - 1).rev() {
                red(k, l, &mut b, &mut h, &d, &mut lam);
            }
            k += 1;
        }
    }
    std::array::from_fn(|r| std::array::from_fn(|t| h[r][t].clone()))
}
