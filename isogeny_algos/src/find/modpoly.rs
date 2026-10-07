//! Classical modular polynomial Phi_l(X,Y) over F_p, built from q-expansions of j reduced mod p:
//! solve for the unique (up to scale) symmetric Phi with Phi(j(q), j(q^l)) having no polar part
//! or constant term. Practical up to l ~ 31 with this baseline linear algebra.
use crate::field::Field;
use crate::poly::{self, Poly};
use crate::series;

pub struct Phi<F: Field> {
    pub ell: usize,
    /// c[i][j] coefficient of X^i Y^j, 0 <= i,j <= l+1
    pub c: Vec<Vec<F::E>>,
}

impl<F: Field> Clone for Phi<F> {
    fn clone(&self) -> Self {
        Phi {
            ell: self.ell,
            c: self.c.clone(),
        }
    }
}

fn sigma3(n: usize) -> u64 {
    let mut s = 0u64;
    for d in 1..=n {
        if n % d == 0 {
            s += (d as u64).pow(3);
        }
    }
    s
}

/// S(q) = q*j(q) = E4^3 / prod(1-q^n)^24 as a power series over F, `len` coefficients.
fn s_series<F: Field>(f: &F, len: usize) -> Vec<F::E> {
    let mut e4 = vec![f.zero(); len];
    e4[0] = f.one();
    let c240 = f.from_u64(240);
    for n in 1..len {
        e4[n] = f.mul(c240, f.from_u64(sigma3(n)));
    }
    let e4_3 = series::mul(f, &series::mul(f, &e4, &e4, len), &e4, len);
    // prod (1-q^n) via pentagonal numbers
    let mut eta = vec![f.zero(); len];
    let mut k: i64 = 0;
    loop {
        let mut any = false;
        for (idx, kk) in [k, -k].into_iter().enumerate() {
            if k == 0 && idx == 1 {
                continue;
            }
            let e = kk * (3 * kk - 1) / 2;
            if e >= 0 && (e as usize) < len {
                any = true;
                let sgn = if kk % 2 == 0 { f.one() } else { f.neg(f.one()) };
                eta[e as usize] = f.add(eta[e as usize], sgn);
            }
        }
        if !any && k > 0 {
            break;
        }
        k += 1;
    }
    let e2 = series::mul(f, &eta, &eta, len);
    let e4s = series::mul(f, &e2, &e2, len);
    let e8 = series::mul(f, &e4s, &e4s, len);
    let e16 = series::mul(f, &e8, &e8, len);
    let e24 = series::mul(f, &e16, &e8, len);
    series::mul(f, &e4_3, &series::inv(f, &e24, len), len)
}

/// Integer coefficients of Phi_l (symmetric residues of a CRT over 62-bit primes; the
/// Bröker–Sutherland height bound log|c| <= 6 l log l + 16 l + 14 sqrt(l) log l sets the count).
pub fn integer_coeffs(ell: usize) -> Vec<Vec<crate::int::Int>> {
    use crate::bigint::Big;
    use crate::field::{is_prime, Zp};
    use crate::int::Int;
    let l = ell as f64;
    let nats = 6.0 * l * l.ln() + 16.0 * l + 14.0 * l.sqrt() * l.ln();
    let bits = (nats / std::f64::consts::LN_2) as usize + 16;
    let mut primes = vec![];
    let mut p = (1u64 << 62) - 57;
    while primes.len() * 61 < bits {
        if is_prime(p) {
            primes.push(p);
        }
        p -= 2;
    }
    let phis: Vec<Phi<Zp>> = primes.iter().map(|&p| Phi::compute(&Zp::new(p), ell)).collect();
    let (rows, cols) = (phis[0].c.len(), phis[0].c[0].len());
    let mut m = Big::from_u64(1);
    for &p in &primes {
        m = m.mul_small(p);
    }
    let mi = Int::from_big(&m);
    (0..rows)
        .map(|i| {
            (0..cols)
                .map(|k| {
                    let mut x = Big::zero();
                    let mut mm = Big::from_u64(1);
                    for (t, &p) in primes.iter().enumerate() {
                        let fp = Zp::new(p);
                        let r = phis[t].c[i][k];
                        let delta = fp.mul(fp.sub(r, x.rem_small(p)), fp.inv(mm.rem_small(p)));
                        x = x.add(&mm.mul_small(delta));
                        mm = mm.mul_small(p);
                    }
                    let xi = Int::from_big(&x);
                    if x.add(&x) > m {
                        &xi - &mi
                    } else {
                        xi
                    }
                })
                .collect()
        })
        .collect()
}

impl<F: Field> Phi<F> {
    /// Phi_l over any field, from the integer coefficients reduced by the characteristic
    /// (valid in every characteristic, including small p and 2).
    pub fn via_crt(f: &F, ell: usize) -> Phi<F> {
        let ch = f.char();
        let c = integer_coeffs(ell)
            .iter()
            .map(|row| row.iter().map(|x| f.from_u64(x.mod_u64(ch))).collect())
            .collect();
        Phi { ell, c }
    }

    /// Phi_l over F: the Hecke/Newton construction (`compute_hecke`), valid for char F > l + 1.
    pub fn compute(f: &F, ell: usize) -> Phi<F> {
        Self::compute_hecke(f, ell)
    }

    /// Phi_l from the factorisation Phi_l(j(q), Y) = (Y - j(q^l)) G(Y), where G's roots are the l
    /// conjugates j(zeta^r q^{1/l}). Their power sums are s_m = l U_l(j^m) (keep the coefficients
    /// of q^{l n} of j^m, as a series in q^n), so Newton's identities give the coefficients of G as
    /// Laurent series with at most a simple pole, to precision q^l; then
    /// E_k = e_k + j(q^l) e_{k-1} is a polynomial of degree <= l + 1 in j, read off its polar part
    /// and constant term. Cost: l + 1 products of series of length ~l^2 plus O(l^2) short
    /// products, against the O(l^6) dense solve of `compute_linear_algebra`.
    pub fn compute_hecke(f: &F, ell: usize) -> Phi<F> {
        assert!(f.char() > ell as u64 + 1 || f.char() == 0, "characteristic must exceed l + 1");
        let l = ell;
        let l1 = l + 1;
        let prec = l + 1; // non-negative powers q^0 .. q^l of the e_k(j_r)
        let len = l * prec + l1 + 2; // indices of S^m needed: n + m with n <= l (prec - 1), m <= l + 1
        let s = s_series(f, len);
        // S^m, m = 0..=l+1 (j^m = q^{-m} S^m)
        let mut sp: Vec<Vec<F::E>> = vec![vec![f.zero(); len]];
        sp[0][0] = f.one();
        for m in 1..=l1 {
            let next = series::mul(f, &sp[m - 1], &s, len);
            sp.push(next);
        }
        // Laurent series with offset 1: index i <-> q^{i-1}, i = 0 ..= prec
        let w = prec + 1;
        let lf = f.from_u64(l as u64);
        // s_m(q) = l * sum_{n: l | n} [j^m]_n q^{n/l}; [j^m]_n = [S^m]_{n+m}
        let mut pw: Vec<Vec<F::E>> = vec![vec![f.zero(); w]]; // pw[m], m = 1..=l
        for m in 1..=l {
            let mut v = vec![f.zero(); w];
            for (i, slot) in v.iter_mut().enumerate() {
                let e = i as i64 - 1; // exponent n/l
                let n = e * l as i64;
                let idx = n + m as i64;
                if idx >= 0 && (idx as usize) < len {
                    *slot = f.mul(lf, sp[m][idx as usize]);
                }
            }
            pw.push(v);
        }
        let lmul = |a: &[F::E], b: &[F::E]| -> Vec<F::E> {
            // offset-1 Laurent product truncated to q^{prec-1}: (q^{i-1})(q^{k-1}) = q^{i+k-2}
            let mut r = vec![f.zero(); w];
            for (i, &x) in a.iter().enumerate() {
                if f.is_zero(x) {
                    continue;
                }
                for (k, &y) in b.iter().enumerate() {
                    let t = i + k;
                    if t >= 1 && t - 1 < w {
                        r[t - 1] = f.add(r[t - 1], f.mul(x, y));
                    }
                }
            }
            r
        };
        // Newton: k e_k = sum_{i=1}^k (-1)^{i-1} e_{k-i} s_i
        let mut e: Vec<Vec<F::E>> = vec![{
            let mut one = vec![f.zero(); w];
            one[1] = f.one();
            one
        }];
        for k in 1..=l {
            let mut acc = vec![f.zero(); w];
            for i in 1..=k {
                let t = lmul(&e[k - i], &pw[i]);
                for (a, b) in acc.iter_mut().zip(t) {
                    *a = if i % 2 == 1 { f.add(*a, b) } else { f.sub(*a, b) };
                }
            }
            let ki = f.inv(f.from_u64(k as u64));
            e.push(acc.into_iter().map(|x| f.mul(x, ki)).collect());
        }
        e.push(vec![f.zero(); w]); // e_{l+1} of l conjugates = 0
        // j(q^l) = sum_n [j]_n q^{l n}, n >= -1
        // E_k = e_k + j(q^l) e_{k-1}: polar order <= l + 1; represent with offset l+1, up to q^0
        let wo = l1 + 1;
        let jpow = |d: usize, ex: i64| -> F::E {
            // coefficient of q^ex in j^d
            let idx = ex + d as i64;
            if idx >= 0 && (idx as usize) < len {
                sp[d][idx as usize]
            } else {
                f.zero()
            }
        };
        let mut cmat = vec![vec![f.zero(); l1 + 1]; l1 + 1];
        cmat[0][l1] = f.one();
        for k in 1..=l1 {
            let mut big = vec![f.zero(); wo]; // index i <-> q^{i - (l+1)}
            for i in 0..w {
                let ex = i as i64 - 1;
                if ex <= 0 {
                    let pos = (ex + l1 as i64) as usize;
                    big[pos] = f.add(big[pos], e[k][i]);
                }
            }
            // j(q^l) e_{k-1}: terms [j]_n q^{l n} times q^{ex}
            for n in -1i64..=1 {
                let c = jpow(1, n);
                if f.is_zero(c) {
                    continue;
                }
                for i in 0..w {
                    let ex = i as i64 - 1 + l as i64 * n;
                    if ex <= 0 && ex >= -(l1 as i64) {
                        let pos = (ex + l1 as i64) as usize;
                        big[pos] = f.add(big[pos], f.mul(c, e[k - 1][i]));
                    }
                }
            }
            // polynomial in X: subtract c_d j^d from the top pole down
            let mut coef = vec![f.zero(); l1 + 1];
            for d in (0..=l1).rev() {
                let c = big[l1 - d]; // coefficient of q^{-d}
                coef[d] = c;
                if !f.is_zero(c) {
                    for ex in -(d as i64)..=0 {
                        let pos = (ex + l1 as i64) as usize;
                        big[pos] = f.sub(big[pos], f.mul(c, jpow(d, ex)));
                    }
                }
            }
            let sign_neg = k % 2 == 1;
            for (d, &c) in coef.iter().enumerate() {
                cmat[d][l1 - k] = if sign_neg { f.neg(c) } else { c };
            }
        }
        Phi { ell, c: cmat }
    }

    /// Phi_l by dense linear algebra on q-expansions (the original method; field of size > 10^4
    /// and > 4l).
    pub fn compute_linear_algebra(f: &F, ell: usize) -> Phi<F> {
        assert!(f.q() > crate::bigint::Big::from_u64(10_000.max(4 * ell as u64)));
        let l1 = ell + 1;
        let top = l1 * l1; // max pole order
        let len = top + 2;
        let s = s_series(f, len);
        // A_i = S^i
        let mut a: Vec<Vec<F::E>> = vec![vec![f.zero(); len]; l1 + 1];
        a[0][0] = f.one();
        for i in 1..=l1 {
            a[i] = series::mul(f, &a[i - 1], &s, len);
        }
        let mut pairs = vec![];
        for i in 0..=l1 {
            for j in i..=l1 {
                pairs.push((i, j));
            }
        }
        let norm = pairs.iter().position(|&(i, j)| i == 0 && j == l1).unwrap();
        let ncols = pairs.len();
        let nrows = top + 1;
        // T(i,j)[k] = coefficient of q^{k-(i+l j)} in J1^i J2^j = (A_i * A_j(q^l))[k]
        let term = |i: usize, j: usize| -> Vec<F::E> {
            let kmax = i + ell * j;
            let mut r = vec![f.zero(); kmax + 1];
            for m in 0..=(kmax / ell) {
                let cj = a[j][m];
                if f.is_zero(cj) {
                    continue;
                }
                let off = m * ell;
                for k in off..=kmax {
                    r[k] = f.add(r[k], f.mul(cj, a[i][k - off]));
                }
            }
            r
        };
        let mut mat = vec![vec![f.zero(); ncols + 1]; nrows];
        for (ci, &(i, j)) in pairs.iter().enumerate() {
            let mut add_term = |ii: usize, jj: usize| {
                let t = term(ii, jj);
                let shift = ii + ell * jj;
                for (k, &v) in t.iter().enumerate() {
                    let e = k as i64 - shift as i64;
                    if e <= 0 {
                        let row = (e + top as i64) as usize;
                        mat[row][ci] = f.add(mat[row][ci], v);
                    }
                }
            };
            add_term(i, j);
            if i != j {
                add_term(j, i);
            }
        }
        for row in mat.iter_mut() {
            row[ncols] = f.neg(row[norm]);
        }
        let cols: Vec<usize> = (0..ncols).filter(|&c| c != norm).collect();
        let mut piv_row = 0;
        let mut piv_of_col = vec![usize::MAX; ncols];
        for &col in &cols {
            let pr = (piv_row..nrows).find(|&r| !f.is_zero(mat[r][col]));
            let Some(pr) = pr else {
                panic!("Phi_{ell}: rank-deficient system (field too small?)")
            };
            mat.swap(piv_row, pr);
            let inv = f.inv(mat[piv_row][col]);
            for c in col..=ncols {
                mat[piv_row][c] = f.mul(mat[piv_row][c], inv);
            }
            let pivot = mat[piv_row].clone();
            for (r, row) in mat.iter_mut().enumerate() {
                if r != piv_row && !f.is_zero(row[col]) {
                    let fct = row[col];
                    for c in col..=ncols {
                        row[c] = f.sub(row[c], f.mul(fct, pivot[c]));
                    }
                }
            }
            piv_of_col[col] = piv_row;
            piv_row += 1;
        }
        for row in mat.iter().skip(piv_row) {
            assert!(f.is_zero(row[ncols]), "Phi_{ell}: inconsistent system");
        }
        let mut cmat = vec![vec![f.zero(); l1 + 1]; l1 + 1];
        for (ci, &(i, j)) in pairs.iter().enumerate() {
            let v = if ci == norm {
                f.one()
            } else {
                mat[piv_of_col[ci]][ncols]
            };
            cmat[i][j] = v;
            cmat[j][i] = v;
        }
        Phi { ell, c: cmat }
    }

    /// The same polynomial over another field through an embedding (e.g. F_p into F_{p^2}).
    pub fn lift<G: Field>(&self, emb: impl Fn(F::E) -> G::E) -> Phi<G> {
        Phi {
            ell: self.ell,
            c: self
                .c
                .iter()
                .map(|row| row.iter().map(|&x| emb(x)).collect())
                .collect(),
        }
    }

    /// Phi(x0, Y) as a polynomial in Y over F.
    pub fn y_poly(&self, f: &F, x0: F::E) -> Poly<F> {
        let l1 = self.ell + 1;
        let mut pw = vec![f.one()];
        for _ in 0..l1 {
            pw.push(f.mul(*pw.last().unwrap(), x0));
        }
        let mut out = vec![];
        for j in 0..=l1 {
            let mut s = f.zero();
            for i in 0..=l1 {
                s = f.add(s, f.mul(self.c[i][j], pw[i]));
            }
            out.push(s);
        }
        poly::trim(f, &mut out);
        out
    }

    /// (Phi_X, Phi_Y) at (x0, y0).
    pub fn grad(&self, f: &F, x0: F::E, y0: F::E) -> (F::E, F::E) {
        let l1 = self.ell + 1;
        let mut xp = vec![f.one()];
        let mut yp = vec![f.one()];
        for _ in 0..l1 {
            xp.push(f.mul(*xp.last().unwrap(), x0));
            yp.push(f.mul(*yp.last().unwrap(), y0));
        }
        let (mut px, mut py) = (f.zero(), f.zero());
        for i in 0..=l1 {
            for j in 0..=l1 {
                let c = self.c[i][j];
                if i > 0 {
                    px = f.add(
                        px,
                        f.mul(f.mul(c, f.from_u64(i as u64)), f.mul(xp[i - 1], yp[j])),
                    );
                }
                if j > 0 {
                    py = f.add(
                        py,
                        f.mul(f.mul(c, f.from_u64(j as u64)), f.mul(xp[i], yp[j - 1])),
                    );
                }
            }
        }
        (px, py)
    }

    /// (Phi_XX, Phi_XY, Phi_YY) at (x0, y0).
    pub fn hessian(&self, f: &F, x0: F::E, y0: F::E) -> (F::E, F::E, F::E) {
        let l1 = self.ell + 1;
        let (mut xp, mut yp) = (vec![f.one()], vec![f.one()]);
        for _ in 0..l1 {
            xp.push(f.mul(*xp.last().unwrap(), x0));
            yp.push(f.mul(*yp.last().unwrap(), y0));
        }
        let (mut pxx, mut pxy, mut pyy) = (f.zero(), f.zero(), f.zero());
        for i in 0..=l1 {
            for j in 0..=l1 {
                let c = self.c[i][j];
                if i > 1 {
                    let w = f.mul(f.from_u64((i * (i - 1)) as u64), f.mul(xp[i - 2], yp[j]));
                    pxx = f.add(pxx, f.mul(c, w));
                }
                if i > 0 && j > 0 {
                    let w = f.mul(f.from_u64((i * j) as u64), f.mul(xp[i - 1], yp[j - 1]));
                    pxy = f.add(pxy, f.mul(c, w));
                }
                if j > 1 {
                    let w = f.mul(f.from_u64((j * (j - 1)) as u64), f.mul(xp[i], yp[j - 2]));
                    pyy = f.add(pyy, f.mul(c, w));
                }
            }
        }
        (pxx, pxy, pyy)
    }

    pub fn eval(&self, f: &F, x0: F::E, y0: F::E) -> F::E {
        poly::eval(f, &self.y_poly(f, x0), y0)
    }

    /// Distinct F-rational roots j' of Phi_l(j, Y).
    pub fn neighbors(&self, f: &F, j: F::E, rng: &mut crate::field::Rng) -> Vec<F::E> {
        poly::roots(f, &self.y_poly(f, j), rng)
    }
}
