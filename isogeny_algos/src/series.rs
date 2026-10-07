//! Truncated power series over a field (dense, schoolbook): the shared substrate of the
//! Elkies/Atkin/Stark family. All functions return exactly `n` coefficients.
use crate::field::Field;

pub fn mul<F: Field>(f: &F, a: &[F::E], b: &[F::E], n: usize) -> Vec<F::E> {
    let mut r = vec![f.zero(); n];
    for i in 0..a.len().min(n) {
        if f.is_zero(a[i]) {
            continue;
        }
        for j in 0..b.len().min(n - i) {
            r[i + j] = f.add(r[i + j], f.mul(a[i], b[j]));
        }
    }
    r
}

/// 1/a (requires a[0] != 0).
pub fn inv<F: Field>(f: &F, a: &[F::E], n: usize) -> Vec<F::E> {
    let a0i = f.inv(a[0]);
    let mut r = vec![f.zero(); n];
    r[0] = a0i;
    for k in 1..n {
        let mut s = f.zero();
        for i in 1..=k.min(a.len() - 1) {
            s = f.add(s, f.mul(a[i], r[k - i]));
        }
        r[k] = f.neg(f.mul(s, a0i));
    }
    r
}

pub fn pow<F: Field>(f: &F, a: &[F::E], mut e: usize, n: usize) -> Vec<F::E> {
    let mut r = vec![f.zero(); n];
    r[0] = f.one();
    let mut b: Vec<F::E> = (0..n)
        .map(|i| if i < a.len() { a[i] } else { f.zero() })
        .collect();
    while e > 0 {
        if e & 1 == 1 {
            r = mul(f, &r, &b, n);
        }
        e >>= 1;
        if e > 0 {
            b = mul(f, &b, &b, n);
        }
    }
    r
}

/// exp(F) for F with F[0] = 0, via G' = F' G.
pub fn exp<F: Field>(f: &F, fs: &[F::E], n: usize) -> Vec<F::E> {
    let mut g = vec![f.zero(); n];
    g[0] = f.one();
    for m in 1..n {
        let mut s = f.zero();
        for k in 1..=m.min(fs.len() - 1) {
            s = f.add(s, f.mul(f.mul(f.from_u64(k as u64), fs[k]), g[m - k]));
        }
        g[m] = f.div(s, f.from_u64(m as u64));
    }
    g
}

/// a(b(x)) with b[0] = 0.
pub fn compose<F: Field>(f: &F, a: &[F::E], b: &[F::E], n: usize) -> Vec<F::E> {
    let mut res = vec![f.zero(); n];
    for k in (0..a.len()).rev() {
        res = mul(f, &res, b, n);
        res[0] = f.add(res[0], a[k]);
    }
    res
}
