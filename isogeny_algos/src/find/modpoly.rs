//! Classical modular polynomial Phi_l(X,Y) over F_p, built from q-expansions of j reduced mod p:
//! solve for the unique (up to scale) symmetric Phi with Phi(j(q), j(q^l)) having no polar part
//! or constant term. Practical up to l ~ 31 with this baseline linear algebra.
use crate::field::{Field, Zp};
use crate::poly::{self, Poly};

pub struct Phi {
    pub ell: usize,
    pub p: u64,
    /// c[i][j] coefficient of X^i Y^j, 0 <= i,j <= l+1
    pub c: Vec<Vec<u64>>,
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

fn smul(fp: &Zp, a: &[u64], b: &[u64], len: usize) -> Vec<u64> {
    let mut r = vec![0u64; len];
    for i in 0..a.len().min(len) {
        if a[i] == 0 {
            continue;
        }
        for j in 0..b.len().min(len - i) {
            r[i + j] = fp.add(r[i + j], fp.mul(a[i], b[j]));
        }
    }
    r
}

/// S(q) = q*j(q) = E4^3 / prod(1-q^n)^24 as a power series mod p, `len` coefficients.
fn s_series(fp: &Zp, len: usize) -> Vec<u64> {
    let mut e4 = vec![0u64; len];
    e4[0] = 1;
    for n in 1..len {
        e4[n] = fp.mul(240 % fp.p, sigma3(n) % fp.p);
    }
    let e4_3 = smul(fp, &smul(fp, &e4, &e4, len), &e4, len);
    // prod (1-q^n) via pentagonal numbers
    let mut eta = vec![0u64; len];
    let mut k: i64 = 0;
    loop {
        let mut any = false;
        for (idx, kk) in [k, -k].into_iter().enumerate() {
            if k == 0 && idx == 1 {
                continue;
            }
            let e = (kk * (3 * kk - 1) / 2) as i64;
            if e >= 0 && (e as usize) < len {
                any = true;
                let sgn = if kk % 2 == 0 { 1 } else { fp.p - 1 };
                eta[e as usize] = fp.add(eta[e as usize], sgn);
            }
        }
        if !any && k > 0 {
            break;
        }
        k += 1;
    }
    let e2 = smul(fp, &eta, &eta, len);
    let e4s = smul(fp, &e2, &e2, len);
    let e8 = smul(fp, &e4s, &e4s, len);
    let e16 = smul(fp, &e8, &e8, len);
    let e24 = smul(fp, &e16, &e8, len);
    // invert e24
    let mut inv = vec![0u64; len];
    inv[0] = 1;
    for kx in 1..len {
        let mut s = 0u64;
        for i in 1..=kx {
            s = fp.add(s, fp.mul(e24[i], inv[kx - i]));
        }
        inv[kx] = fp.neg(s);
    }
    smul(fp, &e4_3, &inv, len)
}

impl Phi {
    pub fn compute(fp: &Zp, ell: usize) -> Phi {
        assert!(fp.p > 10_000 && (fp.p as usize) > 4 * ell);
        let l1 = ell + 1;
        let top = l1 * l1; // max pole order
        let len = top + 2;
        let s = s_series(fp, len);
        // A_i = S^i
        let mut a: Vec<Vec<u64>> = vec![vec![0; len]; l1 + 1];
        a[0][0] = 1;
        for i in 1..=l1 {
            a[i] = smul(fp, &a[i - 1], &s, len);
        }
        // pairs
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
        let term = |i: usize, j: usize| -> Vec<u64> {
            let kmax = i + ell * j; // need k <= kmax
            let mut r = vec![0u64; kmax + 1];
            for m in 0..=(kmax / ell) {
                let cj = a[j][m];
                if cj == 0 {
                    continue;
                }
                let off = m * ell;
                for k in off..=kmax {
                    r[k] = fp.add(r[k], fp.mul(cj, a[i][k - off]));
                }
            }
            r
        };
        let mut mat = vec![vec![0u64; ncols + 1]; nrows];
        for (ci, &(i, j)) in pairs.iter().enumerate() {
            let mut add_term = |ii: usize, jj: usize| {
                let t = term(ii, jj);
                let shift = ii + ell * jj; // exponent e = k - shift, row = e + top
                for (k, &v) in t.iter().enumerate() {
                    let e = k as i64 - shift as i64;
                    if e <= 0 {
                        let row = (e + top as i64) as usize;
                        mat[row][ci] = fp.add(mat[row][ci], v);
                    }
                }
            };
            add_term(i, j);
            if i != j {
                add_term(j, i);
            }
        }
        // move normalisation column to RHS: sum c_k col_k = - col_norm
        for r in 0..nrows {
            mat[r][ncols] = fp.neg(mat[r][norm]);
        }
        let cols: Vec<usize> = (0..ncols).filter(|&c| c != norm).collect();
        // Gaussian elimination on columns `cols`
        let mut piv_row = 0;
        let mut piv_of_col = vec![usize::MAX; ncols];
        for &col in &cols {
            let mut pr = None;
            for r in piv_row..nrows {
                if mat[r][col] != 0 {
                    pr = Some(r);
                    break;
                }
            }
            let Some(pr) = pr else { panic!("Phi_{ell}: rank-deficient system (p too small?)") };
            mat.swap(piv_row, pr);
            let inv = fp.inv(mat[piv_row][col]);
            for c in 0..=ncols {
                mat[piv_row][c] = fp.mul(mat[piv_row][c], inv);
            }
            for r in 0..nrows {
                if r != piv_row && mat[r][col] != 0 {
                    let fct = mat[r][col];
                    for c in col..=ncols {
                        let t = fp.mul(fct, mat[piv_row][c]);
                        mat[r][c] = fp.sub(mat[r][c], t);
                    }
                }
            }
            piv_of_col[col] = piv_row;
            piv_row += 1;
        }
        for r in piv_row..nrows {
            assert_eq!(mat[r][ncols], 0, "Phi_{ell}: inconsistent system");
        }
        let mut cmat = vec![vec![0u64; l1 + 1]; l1 + 1];
        for (ci, &(i, j)) in pairs.iter().enumerate() {
            let v = if ci == norm { 1 } else { mat[piv_of_col[ci]][ncols] };
            cmat[i][j] = v;
            cmat[j][i] = v;
        }
        Phi { ell, p: fp.p, c: cmat }
    }

    /// Phi(x0, Y) as a polynomial in Y over F.
    pub fn y_poly<F: Field>(&self, f: &F, x0: F::E) -> Poly<F> {
        let l1 = self.ell + 1;
        let mut pw = vec![f.one()];
        for _ in 0..l1 {
            pw.push(f.mul(*pw.last().unwrap(), x0));
        }
        let mut out = vec![];
        for j in 0..=l1 {
            let mut s = f.zero();
            for i in 0..=l1 {
                s = f.add(s, f.mul(f.from_u64(self.c[i][j]), pw[i]));
            }
            out.push(s);
        }
        poly::trim(f, &mut out);
        out
    }

    /// (Phi_X, Phi_Y) at (x0, y0).
    pub fn grad<F: Field>(&self, f: &F, x0: F::E, y0: F::E) -> (F::E, F::E) {
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
                let c = f.from_u64(self.c[i][j]);
                if i > 0 {
                    px = f.add(px, f.mul(f.mul(c, f.from_u64(i as u64)), f.mul(xp[i - 1], yp[j])));
                }
                if j > 0 {
                    py = f.add(py, f.mul(f.mul(c, f.from_u64(j as u64)), f.mul(xp[i], yp[j - 1])));
                }
            }
        }
        (px, py)
    }

    pub fn eval<F: Field>(&self, f: &F, x0: F::E, y0: F::E) -> F::E {
        poly::eval(f, &self.y_poly(f, x0), y0)
    }

    /// Distinct F-rational roots j' of Phi_l(j, Y).
    pub fn neighbors<F: Field>(&self, f: &F, j: F::E, rng: &mut crate::field::Rng) -> Vec<F::E> {
        poly::roots(f, &self.y_poly(f, j), rng)
    }
}
