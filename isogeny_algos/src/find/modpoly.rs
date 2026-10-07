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

impl<F: Field> Phi<F> {
    /// Phi_l over F (characteristic must be large compared with the coefficients' denominators:
    /// we require a field of size > 10^4 and > 4l).
    pub fn compute(f: &F, ell: usize) -> Phi<F> {
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
