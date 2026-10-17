//! Exact integer LLL in arbitrary dimension.
//!
//! Rows of the input matrix are the basis (or generating) vectors, entries are [`Int`]s.  The core
//! is the fraction-free *integral LLL* of de Weger / Cohen (Cohen, "A Course in Computational
//! Algebraic Number Theory", Alg. 2.6.7), i.e. the arbitrary-dimension version of
//! [`crate::quat::lll_gram`]: Gram–Schmidt data is kept as the integers
//! `d_i = det(<b_a, b_b>)_{a,b < i}` and `lambda_{k,j} = d_{j+1} mu_{k,j}`, every division is
//! exact, and the Lovász test for a rational `delta = num/den` is
//! `den d_{k+1} d_{k-1} >= num d_k^2 - den lambda_{k,k-1}^2`.
//!
//! Speed-up that keeps exactness: by default the basis is first pre-reduced by a floating-point
//! LLL in the style of Nguyen–Stehlé's L² (exact integer basis and exact integer Gram matrix,
//! `f64` Gram–Schmidt coefficients recomputed from the exact Gram row, lazy size reduction).
//! The floating-point phase only ever applies unimodular integer row operations; whatever it
//! produces is then handed to the exact integral LLL, which certifies (and if necessary
//! finishes) the reduction.  The returned basis therefore always satisfies, exactly,
//! `|mu_{k,j}| <= 1/2` and the Lovász condition for the requested `delta`, independently of any
//! rounding behaviour (see [`is_lll_reduced`]).
//!
//! Linearly dependent rows: the rank is first tested modulo the prime `2^61 - 1` (full rank
//! modulo p proves independence over Q).  If that test fails, an independent basis of the same
//! lattice is computed by unimodular integer row echelon (Euclid on columns) and zero rows are
//! dropped, then LLL runs on that basis.  The output has `rank` rows.
use crate::int::Int;

fn int(v: i64) -> Int {
    Int::from(v)
}

/// Exact inner product.
pub fn dot(a: &[Int], b: &[Int]) -> Int {
    let mut s = Int::zero();
    for (x, y) in a.iter().zip(b) {
        if !x.is_zero() && !y.is_zero() {
            s = &s + &(x * y);
        }
    }
    s
}

/// Parameters for [`lll_with`].
#[derive(Clone, Debug)]
pub struct LllParams {
    /// Lovász constant `delta = delta_num / delta_den`, must satisfy `1/4 < delta < 1`.
    pub delta_num: u64,
    pub delta_den: u64,
    /// Run the floating-point pre-reduction before the exact phase (exactness is unaffected).
    pub fp_prereduce: bool,
    /// Also return the transformation matrix `U` with `U * input = output`.
    pub want_transform: bool,
}

impl Default for LllParams {
    /// `delta = 99/100`, floating-point pre-reduction on, no transform.
    fn default() -> Self {
        LllParams {
            delta_num: 99,
            delta_den: 100,
            fp_prereduce: true,
            want_transform: false,
        }
    }
}

impl LllParams {
    pub fn new(delta_num: u64, delta_den: u64) -> Self {
        LllParams {
            delta_num,
            delta_den,
            ..Default::default()
        }
    }
    /// The textbook `delta = 3/4`.
    pub fn classic() -> Self {
        LllParams::new(3, 4)
    }
    pub fn exact_only(mut self) -> Self {
        self.fp_prereduce = false;
        self
    }
    pub fn with_transform(mut self) -> Self {
        self.want_transform = true;
        self
    }
}

/// Result of [`lll_with`].
#[derive(Clone, Debug)]
pub struct LllOutput {
    /// Reduced basis (`rank` nonzero rows).
    pub basis: Vec<Vec<Int>>,
    /// `U` (rank x #input rows) with `U * input = basis`, when requested.
    pub transform: Option<Vec<Vec<Int>>>,
    pub rank: usize,
    /// Number of swaps done by the exact phase (0 means the fp phase already produced a
    /// reduced basis and the exact phase only certified it).
    pub exact_swaps: usize,
    /// Whether the floating-point phase ran to completion (false: it bailed out and the exact
    /// phase did the remaining work).
    pub fp_completed: bool,
}

/// LLL-reduce the rows with the default parameters (`delta = 99/100`).
pub fn lll(rows: &[Vec<Int>]) -> Vec<Vec<Int>> {
    lll_with(rows, &LllParams::default()).basis
}

/// LLL-reduce the rows (see the module documentation).
pub fn lll_with(rows: &[Vec<Int>], p: &LllParams) -> LllOutput {
    assert!(
        4 * p.delta_num > p.delta_den && p.delta_num < p.delta_den,
        "need 1/4 < delta < 1"
    );
    let r_in = rows.len();
    let m = rows.first().map_or(0, |r| r.len());
    assert!(rows.iter().all(|r| r.len() == m), "ragged matrix");
    let ident = |r: usize| -> Vec<Vec<Int>> {
        (0..r)
            .map(|i| (0..r).map(|j| int((i == j) as i64)).collect())
            .collect()
    };
    // independent basis of the same lattice
    let (mut b, mut u) = if rank_mod_p(rows) == r_in {
        (
            rows.to_vec(),
            if p.want_transform {
                Some(ident(r_in))
            } else {
                None
            },
        )
    } else {
        let (b, u) = echelon(rows);
        (b, if p.want_transform { Some(u) } else { None })
    };
    let mut fp_completed = true;
    if p.fp_prereduce && b.len() > 1 {
        let delta = (p.delta_num as f64 / p.delta_den as f64 + 1.0) / 2.0;
        // aim a bit higher than the target so the exact phase has little left to do
        let delta = delta
            .max(p.delta_num as f64 / p.delta_den as f64)
            .min(0.999);
        // staged: a cheap delta = 3/4 pass removes most of the potential, then the target
        fp_lll(&mut b, u.as_mut(), 0.75);
        fp_completed = fp_lll(&mut b, u.as_mut(), delta);
    }
    let exact_swaps = lll_exact(
        &mut b,
        u.as_mut(),
        &Int::from(p.delta_num),
        &Int::from(p.delta_den),
    );
    let rank = b.len();
    LllOutput {
        basis: b,
        transform: u,
        rank,
        exact_swaps,
        fp_completed,
    }
}

// ---------------------------------------------------------------------------------------------
// rank modulo a prime and integer echelon form

const P61: u64 = (1u64 << 61) - 1;

fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
fn powmod(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, a, p);
        }
        a = mulmod(a, a, p);
        e >>= 1;
    }
    r
}

/// Rank of the rows modulo `2^61 - 1` (a lower bound for the rank over Q).
pub fn rank_mod_p(rows: &[Vec<Int>]) -> usize {
    let p = P61;
    let mut a: Vec<Vec<u64>> = rows
        .iter()
        .map(|r| r.iter().map(|x| x.mod_u64(p)).collect())
        .collect();
    let n = a.len();
    let m = a.first().map_or(0, |r| r.len());
    let mut rank = 0;
    for c in 0..m {
        if rank == n {
            break;
        }
        let Some(piv) = (rank..n).find(|&i| a[i][c] != 0) else {
            continue;
        };
        a.swap(rank, piv);
        let inv = powmod(a[rank][c], p - 2, p);
        for i in rank + 1..n {
            if a[i][c] != 0 {
                let f = mulmod(a[i][c], inv, p);
                for j in c..m {
                    let t = mulmod(f, a[rank][j], p);
                    a[i][j] = (a[i][j] + p - t) % p;
                }
            }
        }
        rank += 1;
    }
    rank
}

fn row_sub(b: &mut [Vec<Int>], k: usize, j: usize, q: &Int) {
    // b_k -= q b_j   (k != j)
    let (bk, bj) = if k > j {
        let (lo, hi) = b.split_at_mut(k);
        (&mut hi[0], &lo[j])
    } else {
        let (lo, hi) = b.split_at_mut(j);
        (&mut lo[k], &hi[0])
    };
    for (x, y) in bk.iter_mut().zip(bj.iter()) {
        if !y.is_zero() {
            *x = &*x - &(q * y);
        }
    }
}

/// Independent basis of the lattice generated by the rows (unimodular integer row echelon,
/// zero rows dropped) together with `U` such that `U * rows = basis`.
pub fn echelon(rows: &[Vec<Int>]) -> (Vec<Vec<Int>>, Vec<Vec<Int>>) {
    let r = rows.len();
    let m = rows.first().map_or(0, |x| x.len());
    let mut a = rows.to_vec();
    let mut u: Vec<Vec<Int>> = (0..r)
        .map(|i| (0..r).map(|j| int((i == j) as i64)).collect())
        .collect();
    let mut pr = 0usize;
    for c in 0..m {
        if pr == r {
            break;
        }
        loop {
            let mut best: Option<usize> = None;
            for i in pr..r {
                if !a[i][c].is_zero() && best.is_none_or(|b| a[i][c].abs() < a[b][c].abs()) {
                    best = Some(i);
                }
            }
            let Some(bi) = best else { break };
            a.swap(pr, bi);
            u.swap(pr, bi);
            let mut done = true;
            for i in pr + 1..r {
                if !a[i][c].is_zero() {
                    let q = a[i][c].div_round(&a[pr][c]);
                    row_sub(&mut a, i, pr, &q);
                    row_sub(&mut u, i, pr, &q);
                    if !a[i][c].is_zero() {
                        done = false;
                    }
                }
            }
            if done {
                pr += 1;
                break;
            }
        }
    }
    a.truncate(pr);
    u.truncate(pr);
    (a, u)
}

// ---------------------------------------------------------------------------------------------
// floating-point pre-reduction (L²-style, exact Gram matrix)

fn ldexp(mut x: f64, mut e: i64) -> f64 {
    while e > 1000 {
        x *= 2f64.powi(1000);
        e -= 1000;
    }
    while e < -1000 {
        x *= 2f64.powi(-1000);
        e += 1000;
    }
    x * 2f64.powi(e as i32)
}

/// `v * 2^-s` as an `f64` (correct exponent for any size of `v`).
fn to_f64_scaled(v: &Int, s: usize) -> f64 {
    let nb = v.bits();
    if nb == 0 {
        return 0.0;
    }
    let (t, e) = if nb > 64 {
        (v.abs().shr(nb - 64), nb as i64 - 64)
    } else {
        (v.abs(), 0)
    };
    let m = t.to_i128().unwrap() as f64;
    let r = ldexp(m, e - s as i64);
    if v.is_neg() {
        -r
    } else {
        r
    }
}

/// Exact conversion of an integral finite `f64`.
fn f64_to_int(x: f64) -> Int {
    if x.abs() < 4.0e18 {
        return int(x as i64);
    }
    let bits = x.to_bits();
    let exp = ((bits >> 52) & 0x7ff) as i64 - 1075;
    let mant = (bits & ((1u64 << 52) - 1)) | (1u64 << 52);
    let v = Int::from(mant).shl(exp as usize);
    if x < 0.0 {
        -v
    } else {
        v
    }
}

/// Returns false if it gave up (non-finite values or no progress); the basis is a valid basis
/// of the same lattice either way.
fn fp_lll(b: &mut [Vec<Int>], mut u: Option<&mut Vec<Vec<Int>>>, delta: f64) -> bool {
    let n = b.len();
    let mut g = vec![vec![Int::zero(); n]; n];
    for i in 0..n {
        for j in 0..=i {
            let v = dot(&b[i], &b[j]);
            g[j][i] = v.clone();
            g[i][j] = v;
        }
    }
    let maxbits = (0..n).map(|i| g[i][i].bits()).max().unwrap_or(0);
    let s = maxbits.saturating_sub(900);
    // f64 mirror of the exact Gram matrix, refreshed whenever an entry changes
    let mut gc: Vec<Vec<f64>> = g
        .iter()
        .map(|row| row.iter().map(|v| to_f64_scaled(v, s)).collect())
        .collect();
    let mut r = vec![vec![0f64; n]; n];
    let mut mu = vec![vec![0f64; n]; n];
    r[0][0] = gc[0][0];
    let mut k = 1usize;
    let mut iters = 0usize;
    let cap = 200_000 + 20 * n * n * (maxbits + 10);
    while k < n {
        iters += 1;
        if iters > cap {
            return false;
        }
        let mut rounds = 0;
        loop {
            for j in 0..k {
                let mut x = gc[k][j];
                for l in 0..j {
                    x -= mu[j][l] * r[k][l];
                }
                r[k][j] = x;
                mu[k][j] = x / r[j][j];
            }
            if (0..k).any(|j| !mu[k][j].is_finite()) {
                return false;
            }
            if (0..k).all(|j| mu[k][j].abs() <= 0.51) {
                break;
            }
            rounds += 1;
            if rounds > 400 {
                return false;
            }
            for j in (0..k).rev() {
                let x = mu[k][j].round();
                if x != 0.0 {
                    let xi = f64_to_int(x);
                    row_sub(b, k, j, &xi);
                    // exact incremental Gram update for b_k -= x b_j
                    let gkk =
                        &(&g[k][k] - &(&(&xi * &g[k][j]) * &int(2))) + &(&(&xi * &xi) * &g[j][j]);
                    for i in 0..n {
                        if i != k {
                            let v = &g[k][i] - &(&xi * &g[j][i]);
                            let vf = to_f64_scaled(&v, s);
                            gc[i][k] = vf;
                            gc[k][i] = vf;
                            g[i][k] = v.clone();
                            g[k][i] = v;
                        }
                    }
                    gc[k][k] = to_f64_scaled(&gkk, s);
                    g[k][k] = gkk;
                    if let Some(u) = u.as_deref_mut() {
                        row_sub(u, k, j, &xi);
                    }
                    for l in 0..j {
                        mu[k][l] -= x * mu[j][l];
                    }
                    mu[k][j] -= x;
                }
            }
        }
        let mut rkk = gc[k][k];
        for j in 0..k {
            rkk -= mu[k][j] * r[k][j];
        }
        r[k][k] = rkk;
        // A non-positive value is cancellation noise for a tiny r_kk: the Lovász test below
        // then fails and the swap is the right move (as in L²); only non-finite values abort.
        if !rkk.is_finite() {
            return false;
        }
        let m1 = mu[k][k - 1];
        if delta * r[k - 1][k - 1] > rkk + m1 * m1 * r[k - 1][k - 1] {
            b.swap(k, k - 1);
            if let Some(u) = u.as_deref_mut() {
                u.swap(k, k - 1);
            }
            g.swap(k, k - 1);
            for row in g.iter_mut() {
                row.swap(k, k - 1);
            }
            gc.swap(k, k - 1);
            for row in gc.iter_mut() {
                row.swap(k, k - 1);
            }
            if k == 1 {
                r[0][0] = gc[0][0];
            } else {
                k -= 1;
            }
        } else {
            k += 1;
        }
    }
    true
}

// ---------------------------------------------------------------------------------------------
// exact integral LLL (Cohen 2.6.7)

struct Exact<'a> {
    b: &'a mut [Vec<Int>],
    u: Option<&'a mut Vec<Vec<Int>>>,
    d: Vec<Int>,
    lam: Vec<Vec<Int>>,
    kmax: usize,
}

impl Exact<'_> {
    fn red(&mut self, k: usize, l: usize) {
        let two_l = &self.lam[k][l] * &int(2);
        if two_l.abs() > self.d[l + 1] {
            let q = self.lam[k][l].div_round(&self.d[l + 1]);
            row_sub(self.b, k, l, &q);
            if let Some(u) = self.u.as_deref_mut() {
                row_sub(u, k, l, &q);
            }
            self.lam[k][l] = &self.lam[k][l] - &(&q * &self.d[l + 1]);
            for i in 0..l {
                let t = &q * &self.lam[l][i];
                self.lam[k][i] = &self.lam[k][i] - &t;
            }
        }
    }
    fn swap(&mut self, k: usize) {
        self.b.swap(k, k - 1);
        if let Some(u) = self.u.as_deref_mut() {
            u.swap(k, k - 1);
        }
        for j in 0..k - 1 {
            let t = std::mem::replace(&mut self.lam[k][j], Int::zero());
            self.lam[k][j] = std::mem::replace(&mut self.lam[k - 1][j], t);
        }
        let l = self.lam[k][k - 1].clone();
        let d = &self.d;
        let bb = (&(&d[k - 1] * &d[k + 1]) + &(&l * &l)).div_floor(&d[k]);
        for i in k + 1..=self.kmax {
            let t = self.lam[i][k].clone();
            let nk = (&(&self.d[k + 1] * &self.lam[i][k - 1]) - &(&l * &t)).div_floor(&self.d[k]);
            let nk1 = (&(&bb * &t) + &(&l * &nk)).div_floor(&self.d[k + 1]);
            self.lam[i][k] = nk;
            self.lam[i][k - 1] = nk1;
        }
        self.d[k] = bb;
    }
    fn gs_row(&mut self, k: usize) {
        for j in 0..=k {
            let mut uu = dot(&self.b[k], &self.b[j]);
            for i in 0..j {
                uu = (&(&self.d[i + 1] * &uu) - &(&self.lam[k][i] * &self.lam[j][i]))
                    .div_floor(&self.d[i]);
            }
            if j < k {
                self.lam[k][j] = uu;
            } else {
                assert!(!uu.is_zero(), "lll_exact: dependent vectors");
                self.d[k + 1] = uu;
            }
        }
    }
}

/// Exact integral LLL on linearly independent rows; returns the number of swaps.
fn lll_exact(b: &mut [Vec<Int>], u: Option<&mut Vec<Vec<Int>>>, num: &Int, den: &Int) -> usize {
    let n = b.len();
    if n <= 1 {
        return 0;
    }
    let mut st = Exact {
        b,
        u,
        d: vec![Int::one(); n + 1],
        lam: vec![vec![Int::zero(); n]; n],
        kmax: 0,
    };
    st.gs_row(0);
    let mut k = 1usize;
    let mut swaps = 0usize;
    while k < n {
        if k > st.kmax {
            st.kmax = k;
            st.gs_row(k);
        }
        st.red(k, k - 1);
        let lhs = &(den * &st.d[k + 1]) * &st.d[k - 1];
        let lk = &st.lam[k][k - 1];
        let rhs = &(num * &(&st.d[k] * &st.d[k])) - &(den * &(lk * lk));
        if lhs < rhs {
            st.swap(k);
            swaps += 1;
            if k > 1 {
                k -= 1;
            }
        } else {
            for l in (0..k - 1).rev() {
                st.red(k, l);
            }
            k += 1;
        }
    }
    swaps
}

// ---------------------------------------------------------------------------------------------
// exact checks

/// Integral Gram–Schmidt data `(d, lambda)` (`d[0] = 1`, `d[i]` = Gram determinant of the first
/// `i` rows, `lambda[k][j] = d[j+1] mu_{k,j}`), or `None` if the rows are dependent.
pub fn integral_gram_schmidt(rows: &[Vec<Int>]) -> Option<(Vec<Int>, Vec<Vec<Int>>)> {
    let n = rows.len();
    let mut d = vec![Int::one(); n + 1];
    let mut lam = vec![vec![Int::zero(); n]; n];
    for k in 0..n {
        for j in 0..=k {
            let mut uu = dot(&rows[k], &rows[j]);
            for i in 0..j {
                uu = (&(&d[i + 1] * &uu) - &(&lam[k][i] * &lam[j][i])).div_floor(&d[i]);
            }
            if j < k {
                lam[k][j] = uu;
            } else {
                if uu.is_zero() {
                    return None;
                }
                d[k + 1] = uu;
            }
        }
    }
    Some((d, lam))
}

/// Exact test: rows independent, `|mu_{k,j}| <= 1/2` and Lovász with `delta = num/den`.
pub fn is_lll_reduced(rows: &[Vec<Int>], num: u64, den: u64) -> bool {
    let Some((d, lam)) = integral_gram_schmidt(rows) else {
        return false;
    };
    let (num, den) = (Int::from(num), Int::from(den));
    let n = rows.len();
    for k in 0..n {
        for j in 0..k {
            if (&lam[k][j] * &int(2)).abs() > d[j + 1] {
                return false;
            }
        }
        if k >= 1 {
            let lhs = &(&den * &d[k + 1]) * &d[k - 1];
            let rhs = &(&num * &(&d[k] * &d[k])) - &(&den * &(&lam[k][k - 1] * &lam[k][k - 1]));
            if lhs < rhs {
                return false;
            }
        }
    }
    true
}

/// Gram determinant `det(B B^T)` of the rows (0 if dependent).
pub fn gram_det(rows: &[Vec<Int>]) -> Int {
    match integral_gram_schmidt(rows) {
        Some((d, _)) => d[rows.len()].clone(),
        None => Int::zero(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Rng;
    use crate::quat::lll_gram;

    fn rnd(rng: &mut Rng, bits: u32) -> Int {
        let v = Int::from(rng.next() >> (64 - bits.min(63)));
        if rng.next() & 1 == 1 {
            -v
        } else {
            v
        }
    }
    fn rnd_big(rng: &mut Rng, bits: usize) -> Int {
        let mut v = Int::zero();
        let mut got = 0;
        while got < bits {
            v = &v.shl(32) + &Int::from(rng.next() >> 32);
            got += 32;
        }
        v.shr(got - bits)
    }

    #[test]
    fn matches_lll_gram_4x4() {
        let mut rng = Rng::new(7);
        for _ in 0..30 {
            let b: Vec<Vec<Int>> = (0..4)
                .map(|_| (0..4).map(|_| rnd(&mut rng, 30)).collect())
                .collect();
            if rank_mod_p(&b) < 4 {
                continue;
            }
            let g: [[Int; 4]; 4] =
                std::array::from_fn(|i| std::array::from_fn(|j| dot(&b[i], &b[j])));
            let h = lll_gram(&g);
            let want: Vec<Vec<Int>> = (0..4)
                .map(|i| {
                    (0..4)
                        .map(|t| {
                            let mut s = Int::zero();
                            for k in 0..4 {
                                s = &s + &(&h[i][k] * &b[k][t]);
                            }
                            s
                        })
                        .collect()
                })
                .collect();
            let out = lll_with(&b, &LllParams::classic().exact_only().with_transform());
            assert_eq!(out.basis, want);
            let u = out.transform.unwrap();
            for i in 0..4 {
                assert_eq!(u[i], h[i].to_vec());
            }
            // fp path: still reduced, same lattice (same Gram determinant)
            let out2 = lll_with(&b, &LllParams::classic());
            assert!(is_lll_reduced(&out2.basis, 3, 4));
            assert_eq!(gram_det(&out2.basis), gram_det(&b));
        }
    }

    #[test]
    fn reduced_and_transform() {
        let mut rng = Rng::new(11);
        for n in [2usize, 5, 9, 14] {
            let b: Vec<Vec<Int>> = (0..n)
                .map(|_| (0..n + 2).map(|_| rnd(&mut rng, 40)).collect())
                .collect();
            for p in [
                LllParams::default().with_transform(),
                LllParams::classic().exact_only().with_transform(),
            ] {
                let out = lll_with(&b, &p);
                assert_eq!(out.rank, n);
                assert!(is_lll_reduced(&out.basis, p.delta_num, p.delta_den));
                let u = out.transform.unwrap();
                for i in 0..n {
                    let row: Vec<Int> = (0..n + 2)
                        .map(|t| {
                            let mut s = Int::zero();
                            for k in 0..n {
                                s = &s + &(&u[i][k] * &b[k][t]);
                            }
                            s
                        })
                        .collect();
                    assert_eq!(row, out.basis[i]);
                }
                assert_eq!(gram_det(&out.basis), gram_det(&b));
            }
        }
    }

    #[test]
    fn dependent_rows() {
        let mut rng = Rng::new(3);
        let n = 6;
        let base: Vec<Vec<Int>> = (0..n)
            .map(|_| (0..8).map(|_| rnd(&mut rng, 20)).collect())
            .collect();
        let mut rows = base.clone();
        // extra dependent rows and a zero row
        for _ in 0..4 {
            let mut v = vec![Int::zero(); 8];
            for b in &base {
                let c = rnd(&mut rng, 4);
                for t in 0..8 {
                    v[t] = &v[t] + &(&c * &b[t]);
                }
            }
            rows.push(v);
        }
        rows.push(vec![Int::zero(); 8]);
        rows.swap(0, 7);
        let out = lll_with(&rows, &LllParams::default().with_transform());
        assert_eq!(out.rank, n);
        assert!(is_lll_reduced(&out.basis, 99, 100));
        assert_eq!(gram_det(&out.basis), gram_det(&base));
        let u = out.transform.unwrap();
        for i in 0..n {
            let row: Vec<Int> = (0..8)
                .map(|t| {
                    let mut s = Int::zero();
                    for k in 0..rows.len() {
                        s = &s + &(&u[i][k] * &rows[k][t]);
                    }
                    s
                })
                .collect();
            assert_eq!(row, out.basis[i]);
        }
        // lattice with a sublattice index: 2 e1, 3 e1 -> e1
        let r2 = vec![
            vec![int(2), int(0)],
            vec![int(3), int(0)],
            vec![int(0), int(5)],
        ];
        let o = lll(&r2);
        assert_eq!(o.len(), 2);
        assert_eq!(gram_det(&o), int(25));
    }

    #[test]
    fn knapsack_planted() {
        let mut rng = Rng::new(99);
        for &(n, bits) in &[(20usize, 60usize), (30, 100), (40, 140)] {
            let a: Vec<Int> = (0..n).map(|_| rnd_big(&mut rng, bits)).collect();
            let s: Vec<bool> = (0..n).map(|_| rng.next() & 1 == 1).collect();
            let mut tot = Int::zero();
            for i in 0..n {
                if s[i] {
                    tot = &tot + &a[i];
                }
            }
            let big = Int::from(1u64 << 20);
            let mut rows = vec![];
            for i in 0..n {
                let mut r = vec![Int::zero(); n + 1];
                r[i] = int(2);
                r[n] = &a[i] * &big;
                rows.push(r);
            }
            let mut r = vec![int(1); n + 1];
            r[n] = &tot * &big;
            rows.push(r);
            let t0 = std::time::Instant::now();
            let out = lll_with(&rows, &LllParams::default());
            assert!(is_lll_reduced(&out.basis, 99, 100));
            let found = out.basis.iter().any(|v| {
                v[n].is_zero()
                    && (0..n).all(|i| {
                        let want = if s[i] { int(-1) } else { int(1) };
                        v[i] == want
                    })
                    || v[n].is_zero()
                        && (0..n).all(|i| {
                            let want = if s[i] { int(1) } else { int(-1) };
                            v[i] == want
                        })
            });
            eprintln!(
                "knapsack n={n} bits={bits}: found={found} exact_swaps={} fp_ok={} {:?}",
                out.exact_swaps,
                out.fp_completed,
                t0.elapsed()
            );
            assert!(found);
        }
    }
}
