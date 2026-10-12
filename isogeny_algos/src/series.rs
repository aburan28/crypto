//! Truncated power series over a field (dense, schoolbook): the shared substrate of the
//! Elkies/Atkin/Stark family. All functions return exactly `n` coefficients.
use crate::field::Field;

pub fn mul<F: Field>(f: &F, a: &[F::E], b: &[F::E], n: usize) -> Vec<F::E> {
    let (la, lb) = (a.len().min(n), b.len().min(n));
    if la >= 128 && lb >= 128 {
        let mut r = crate::poly::mul_raw(f, &a[..la], &b[..lb]);
        r.resize(n.max(r.len()), f.zero());
        r.truncate(n);
        return r;
    }
    if la == 0 || lb == 0 {
        return vec![f.zero(); n];
    }
    f.conv_trunc(&a[..la], &b[..lb], n)
}

/// 1/a (requires a[0] != 0): Newton iteration g <- g (2 - a g), O(M(n)).
pub fn inv<F: Field>(f: &F, a: &[F::E], n: usize) -> Vec<F::E> {
    if n < 512 {
        return inv_quadratic(f, a, n);
    }
    let mut g = inv_quadratic(f, a, 128);
    let mut k = 128;
    while k < n {
        let k2 = (2 * k).min(n);
        let ag = mul(f, &a[..a.len().min(k2)], &g, k2);
        // 2 - a g
        let mut t: Vec<F::E> = ag.iter().map(|&c| f.neg(c)).collect();
        t[0] = f.add(t[0], f.from_u64(2));
        g = mul(f, &g, &t, k2);
        k = k2;
    }
    g.truncate(n);
    g
}

/// Quadratic-time inverse (reference).
pub fn inv_quadratic<F: Field>(f: &F, a: &[F::E], n: usize) -> Vec<F::E> {
    let a0i = f.inv(a[0]);
    let mut r = Vec::with_capacity(n);
    r.push(a0i);
    for k in 1..n {
        // s = sum_{i=1}^{k} a[i] r[k-i]
        let top = k.min(a.len() - 1);
        let s = f.dot_rev(&a[1..=top], &r[k - top..k]);
        r.push(f.neg(f.mul(s, a0i)));
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
    // k F_k, precomputed
    let kf: Vec<F::E> = (0..fs.len().min(n))
        .map(|k| f.mul(f.from_u64(k as u64), fs[k]))
        .collect();
    let invs = f.batch_inv(&(1..n).map(|m| f.from_u64(m as u64)).collect::<Vec<_>>());
    for m in 1..n {
        let top = m.min(kf.len() - 1);
        let s = f.dot_rev(&kf[1..=top], &g[m - top..m]);
        g[m] = f.mul(s, invs[m - 1]);
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
