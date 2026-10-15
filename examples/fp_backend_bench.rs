//! Timing harness for `crypto::cryptanalysis::fp_backend`.
//!
//! Every arm is checked for correctness before it is timed. Arms are
//! interleaved: each repetition runs every arm once, in an order rotated
//! by the repetition index, so drift on a shared host spreads over all
//! arms. Each measurement runs the arm's inner loop for about
//! `--target-ms` milliseconds (iterations calibrated once, then frozen)
//! and records nanoseconds per operation.
//!
//! ```text
//! cargo run --release --example fp_backend_bench -- \
//!     --reps 7 --target-ms 30 --out DIR [--filter SUBSTR]
//! ```
//!
//! Writes `DIR/raw.csv` (one row per arm and repetition) and
//! `DIR/summary.csv` (median, quartiles, min, max per arm), and prints
//! the summary. The f4_fp_tower kernels below are verbatim copies of the
//! private `Fp::reduce`, `lazy_axpy`, `Kernel::axpy`/`sub` and
//! `Acc::sub_tail` arithmetic in `src/cryptanalysis/f4_fp_tower.rs` at the
//! commit this file was added (the engine does not export them).

use std::hint::black_box;
use std::io::Write;
use std::time::Instant;

use crypto_lib::cryptanalysis::fp_backend::{
    sqrt_tonelli_shanks_plain, Fp192, Fp224, Fp256, Fp32, Fp64, FpP192, FpWide, PrimeField,
};
use crypto_lib::cryptanalysis::prime_fast;
use crypto_lib::ct_bignum::{MontgomeryContext, Uint};
use num_bigint::{BigUint, RandBigInt};
use num_traits::One;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

// ── f4_fp_tower baseline (copied arithmetic) ─────────────────────────

#[derive(Clone, Copy)]
struct TowerFp {
    p: u64,
    m: u64,
}

impl TowerFp {
    fn new(p: u64) -> Self {
        TowerFp {
            p,
            m: (u128::from(u64::MAX) + 1).div_euclid(u128::from(p)) as u64,
        }
    }
    #[inline(always)]
    fn reduce(self, x: u64) -> u64 {
        let q = ((u128::from(x) * u128::from(self.m)) >> 64) as u64;
        let r = x - q * self.p;
        if r >= self.p {
            r - self.p
        } else {
            r
        }
    }
    #[inline(always)]
    fn mul(self, a: u64, b: u64) -> u64 {
        self.reduce(a * b)
    }
    #[inline(always)]
    fn neg(self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
}

#[inline(always)]
fn tower_lazy_axpy(acc: &mut [u64], a: u32, b: &[u32], m: u64) {
    let n = acc.len().min(b.len());
    for (x, &y) in acc[..n].iter_mut().zip(&b[..n]) {
        let s = *x + u64::from(a) * u64::from(y);
        *x = s - ((s >> 63).wrapping_neg() & m);
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn tower_lazy_axpy_avx2(acc: &mut [u64], a: u32, b: &[u32], m: u64) {
    tower_lazy_axpy(acc, a, b, m)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn tower_lazy_axpy_avx512(acc: &mut [u64], a: u32, b: &[u32], m: u64) {
    tower_lazy_axpy(acc, a, b, m)
}

#[derive(Clone, Copy, PartialEq)]
enum TowerLanes {
    Scalar,
    Avx2,
    Avx512,
}

#[derive(Clone, Copy)]
struct TowerKernel {
    fp: TowerFp,
    lazy: Option<u64>,
    lanes: TowerLanes,
}

impl TowerKernel {
    fn new(fp: TowerFp) -> Self {
        let mut lanes = TowerLanes::Scalar;
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("avx512f") {
                lanes = TowerLanes::Avx512;
            } else if std::is_x86_feature_detected!("avx2") {
                lanes = TowerLanes::Avx2;
            }
        }
        TowerKernel {
            fp,
            lazy: (fp.p < 1 << 31).then(|| ((1u64 << 63) / fp.p) * fp.p),
            lanes,
        }
    }
    #[inline]
    fn axpy(self, acc: &mut [u64], a: u32, b: &[u32]) {
        match (self.lazy, self.lanes) {
            #[cfg(target_arch = "x86_64")]
            // SAFETY: detected at construction.
            (Some(m), TowerLanes::Avx512) => unsafe { tower_lazy_axpy_avx512(acc, a, b, m) },
            #[cfg(target_arch = "x86_64")]
            // SAFETY: detected at construction.
            (Some(m), TowerLanes::Avx2) => unsafe { tower_lazy_axpy_avx2(acc, a, b, m) },
            (Some(m), _) => tower_lazy_axpy(acc, a, b, m),
            (None, _) => {
                let fp = self.fp;
                for (x, &y) in acc.iter_mut().zip(b) {
                    *x = fp.reduce(*x + u64::from(a) * u64::from(y));
                }
            }
        }
    }
    #[inline]
    fn sub(self, acc: &mut [u64], v: u64, b: &[u32]) -> u64 {
        self.axpy(acc, self.fp.neg(self.fp.reduce(v)) as u32, b);
        b.len() as u64
    }
}

/// `Acc::sub_tail`'s arithmetic and bitset work on one sparse row.
#[inline]
fn tower_sub_tail(
    fp: TowerFp,
    vals: &mut [u64],
    bits: &mut [u64],
    f: u64,
    cols: &[u32],
    v: &[u32],
) {
    let nf = fp.neg(f);
    for (&c, &x) in cols.iter().zip(v) {
        let ci = c as usize;
        let s = vals[ci] + nf * u64::from(x);
        vals[ci] = fp.reduce(s);
        bits[ci / 64] |= 1u64 << (ci % 64);
    }
}

// ── harness ──────────────────────────────────────────────────────────

struct Arm {
    name: String,
    backend: String,
    prime: String,
    op: String,
    /// Operations performed by one call of `run(1)`.
    ops: u64,
    iters: u64,
    run: Box<dyn FnMut(u64)>,
}

fn arm(
    arms: &mut Vec<Arm>,
    backend: &str,
    prime: &str,
    op: &str,
    ops: u64,
    run: impl FnMut(u64) + 'static,
) {
    arms.push(Arm {
        name: format!("{backend}/{prime}/{op}"),
        backend: backend.into(),
        prime: prime.into(),
        op: op.into(),
        ops,
        iters: 0,
        run: Box::new(run),
    });
}

fn hex(s: &str) -> BigUint {
    BigUint::parse_bytes(s.as_bytes(), 16).unwrap()
}

fn is_probable_prime(n: &BigUint) -> bool {
    let one = BigUint::one();
    let nm1 = n - &one;
    let s = nm1.trailing_zeros().unwrap();
    let d = &nm1 >> s;
    'a: for a in [2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = BigUint::from(a).modpow(&d, n);
        if x == one || x == nm1 {
            continue;
        }
        for _ in 1..s {
            x = (&x * &x) % n;
            if x == nm1 {
                continue 'a;
            }
        }
        return false;
    }
    true
}

const N_ARR: usize = 4096;
const K_ROWS: usize = 16;

/// The generic arms for one backend at one prime.
fn field_arms<F: PrimeField + 'static>(arms: &mut Vec<Arm>, f: F, prime: &str, seed: u64) {
    let backend = f.name().to_string();
    let p = f.modulus().clone();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut rand_vec = |n: usize| -> Vec<F::Elem> {
        (0..n)
            .map(|_| f.from_biguint(&rng.gen_biguint_below(&p)))
            .collect()
    };
    let a = rand_vec(N_ARR);
    let b = rand_vec(N_ARR);
    let nz: Vec<F::Elem> = a
        .iter()
        .copied()
        .filter(|&x| !f.is_zero(x))
        .take(256)
        .collect();
    let squares: Vec<F::Elem> = a.iter().take(64).map(|&x| f.sqr(x)).collect();
    // Correctness spot checks against BigUint before timing.
    for i in 0..64 {
        let (x, y) = (f.to_biguint(a[i]), f.to_biguint(b[i]));
        assert_eq!(f.to_biguint(f.mul(a[i], b[i])), (&x * &y) % &p);
        let r = f.sqrt(squares[i]).unwrap();
        assert_eq!(f.sqr(r), squares[i]);
    }

    {
        let (f, a, b) = (f.clone(), a.clone(), b.clone());
        let mut out = vec![f.zero(); N_ARR];
        arm(arms, &backend, prime, "mul_tp", N_ARR as u64, move |it| {
            for _ in 0..it {
                for i in 0..N_ARR {
                    out[i] = f.mul(black_box(a[i]), b[i]);
                }
                black_box(&out);
            }
        });
    }
    {
        let (f, x0, y) = (f.clone(), a[1], b[2]);
        arm(arms, &backend, prime, "mul_lat", 1024, move |it| {
            let mut x = black_box(x0);
            for _ in 0..it * 1024 {
                x = f.mul(x, y);
            }
            black_box(x);
        });
    }
    {
        let (f, x0) = (f.clone(), a[3]);
        arm(arms, &backend, prime, "sqr_lat", 1024, move |it| {
            let mut x = black_box(x0);
            for _ in 0..it * 1024 {
                x = f.sqr(x);
            }
            black_box(x);
        });
    }
    {
        let (f, a, b) = (f.clone(), a.clone(), b.clone());
        let mut out = vec![f.zero(); N_ARR];
        arm(arms, &backend, prime, "add_tp", N_ARR as u64, move |it| {
            for _ in 0..it {
                for i in 0..N_ARR {
                    out[i] = f.add(black_box(a[i]), b[i]);
                }
                black_box(&out);
            }
        });
    }
    {
        let (f, nz) = (f.clone(), nz.clone());
        let n = nz.len() as u64;
        arm(arms, &backend, prime, "inv", n, move |it| {
            for _ in 0..it {
                for &x in &nz {
                    black_box(f.inv(black_box(x)));
                }
            }
        });
    }
    {
        let (f, a) = (f.clone(), a.clone());
        let mut v = a.clone();
        let mut scratch = Vec::new();
        arm(
            arms,
            &backend,
            prime,
            "batch_inv",
            N_ARR as u64,
            move |it| {
                for _ in 0..it {
                    v.copy_from_slice(&a);
                    f.batch_inv(&mut v, &mut scratch);
                    black_box(&v);
                }
            },
        );
    }
    {
        let (f, sq) = (f.clone(), squares.clone());
        let n = sq.len() as u64;
        arm(arms, &backend, prime, "sqrt", n, move |it| {
            for _ in 0..it {
                for &x in &sq {
                    black_box(f.sqrt(black_box(x)));
                }
            }
        });
    }
    {
        let (f, a) = (f.clone(), a.clone());
        arm(arms, &backend, prime, "legendre", 64, move |it| {
            for _ in 0..it {
                for &x in &a[..64] {
                    black_box(f.legendre(black_box(x)));
                }
            }
        });
    }
    // Row kernels: K_ROWS updates of one row of length L from 4 distinct
    // pivot rows, then one reduction; ns per entry update.
    for &len in &[10_000usize, 100_000] {
        let rows: Vec<Vec<F::Elem>> = (0..4).map(|_| rand_vec(len)).collect();
        let start = rand_vec(len);
        let cs = rand_vec(K_ROWS);
        // Correctness: lazy equals eager.
        {
            let mut acc = vec![f.wide_zero(); len];
            f.row_load(&start, &mut acc);
            let mut eager = start.clone();
            for k in 0..K_ROWS {
                f.row_sub(&mut acc, cs[k], &rows[k % 4]);
                f.row_axpy_eager(&mut eager, f.neg(cs[k]), &rows[k % 4]);
            }
            let mut out = vec![f.zero(); len];
            f.row_reduce(&acc, &mut out);
            assert_eq!(out, eager, "{backend} lazy rows");
        }
        let tag = if len == 10_000 { "1e4" } else { "1e5" };
        {
            let (f, rows, start, cs) = (f.clone(), rows.clone(), start.clone(), cs.clone());
            let mut acc = vec![f.wide_zero(); len];
            let mut out = vec![f.zero(); len];
            arm(
                arms,
                &backend,
                prime,
                &format!("axpy_lazy_{tag}"),
                (K_ROWS * len) as u64,
                move |it| {
                    for _ in 0..it {
                        f.row_load(&start, &mut acc);
                        for k in 0..K_ROWS {
                            f.row_sub(&mut acc, cs[k], &rows[k % 4]);
                        }
                        f.row_reduce(&acc, &mut out);
                        black_box(&out);
                    }
                },
            );
        }
        {
            let (f, rows, start, cs) = (f.clone(), rows.clone(), start.clone(), cs.clone());
            let mut dst = start.clone();
            let ncs: Vec<F::Elem> = cs.iter().map(|&c| f.neg(c)).collect();
            arm(
                arms,
                &backend,
                prime,
                &format!("axpy_eager_{tag}"),
                (K_ROWS * len) as u64,
                move |it| {
                    for _ in 0..it {
                        dst.copy_from_slice(&start);
                        for k in 0..K_ROWS {
                            f.row_axpy_eager(&mut dst, ncs[k], &rows[k % 4]);
                        }
                        black_box(&dst);
                    }
                },
            );
        }
    }
    {
        let len = 10_000;
        let u = rand_vec(len);
        let v = rand_vec(len);
        let naive = u
            .iter()
            .zip(&v)
            .fold(f.zero(), |s, (&x, &y)| f.add(s, f.mul(x, y)));
        assert_eq!(f.dot(&u, &v), naive);
        {
            let (f, u, v) = (f.clone(), u.clone(), v.clone());
            arm(
                arms,
                &backend,
                prime,
                "dot_lazy_1e4",
                len as u64,
                move |it| {
                    for _ in 0..it {
                        black_box(f.dot(black_box(&u), &v));
                    }
                },
            );
        }
        {
            let f = f.clone();
            arm(
                arms,
                &backend,
                prime,
                "dot_eager_1e4",
                len as u64,
                move |it| {
                    for _ in 0..it {
                        let s = black_box(&u)
                            .iter()
                            .zip(&v)
                            .fold(f.zero(), |s, (&x, &y)| f.add(s, f.mul(x, y)));
                        black_box(s);
                    }
                },
            );
        }
    }
}

/// The f4_fp_tower kernels and the matching Fp32 kernels at a 31/32-bit prime.
fn tower_arms(arms: &mut Vec<Arm>, p: u64, prime: &str, seed: u64) {
    let tfp = TowerFp::new(p);
    let kern = TowerKernel::new(tfp);
    let f = Fp32::new(p);
    let mut rng = StdRng::seed_from_u64(seed);
    let a: Vec<u32> = (0..N_ARR).map(|_| rng.gen_range(0..p) as u32).collect();
    let b: Vec<u32> = (0..N_ARR).map(|_| rng.gen_range(0..p) as u32).collect();
    {
        let (a, b) = (a.clone(), b.clone());
        let mut out = vec![0u32; N_ARR];
        arm(
            arms,
            "f4_fp_tower",
            prime,
            "mul_tp",
            N_ARR as u64,
            move |it| {
                for _ in 0..it {
                    for i in 0..N_ARR {
                        out[i] = tfp.mul(u64::from(black_box(a[i])), u64::from(b[i])) as u32;
                    }
                    black_box(&out);
                }
            },
        );
    }
    // Dense: the elimination kernel (Kernel::sub, then the reduce that
    // builds `bred`).
    for &len in &[10_000usize, 100_000] {
        let tag = if len == 10_000 { "1e4" } else { "1e5" };
        let rows: Vec<Vec<u32>> = (0..4)
            .map(|_| (0..len).map(|_| rng.gen_range(0..p) as u32).collect())
            .collect();
        let start: Vec<u32> = (0..len).map(|_| rng.gen_range(0..p) as u32).collect();
        let cs: Vec<u32> = (0..K_ROWS).map(|_| rng.gen_range(0..p) as u32).collect();
        // Both must agree.
        {
            let mut acc: Vec<u64> = start.iter().map(|&x| u64::from(x)).collect();
            for k in 0..K_ROWS {
                kern.sub(&mut acc, u64::from(cs[k]), &rows[k % 4]);
            }
            let want: Vec<u32> = acc.iter().map(|&x| tfp.reduce(x) as u32).collect();
            let mut wacc = vec![0u64; len];
            f.row_load(&start, &mut wacc);
            for k in 0..K_ROWS {
                f.row_sub(&mut wacc, cs[k], &rows[k % 4]);
            }
            let mut out = vec![0u32; len];
            f.row_reduce(&wacc, &mut out);
            assert_eq!(out, want, "Fp32 vs f4_fp_tower dense");
        }
        {
            let (rows, start, cs) = (rows.clone(), start.clone(), cs.clone());
            let mut acc = vec![0u64; len];
            let mut out = vec![0u32; len];
            arm(
                arms,
                "f4_fp_tower",
                prime,
                &format!("axpy_dense_{tag}"),
                (K_ROWS * len) as u64,
                move |it| {
                    for _ in 0..it {
                        for (x, &s) in acc.iter_mut().zip(&start) {
                            *x = u64::from(s);
                        }
                        for k in 0..K_ROWS {
                            kern.sub(&mut acc, u64::from(cs[k]), &rows[k % 4]);
                        }
                        for (o, &x) in out.iter_mut().zip(&acc) {
                            *o = tfp.reduce(x) as u32;
                        }
                        black_box(&out);
                    }
                },
            );
        }
        {
            let f = f.clone();
            let mut acc = vec![0u64; len];
            let mut out = vec![0u32; len];
            arm(
                arms,
                "Fp32",
                prime,
                &format!("axpy_dense_{tag}"),
                (K_ROWS * len) as u64,
                move |it| {
                    for _ in 0..it {
                        f.row_load(&start, &mut acc);
                        for k in 0..K_ROWS {
                            f.row_sub(&mut acc, cs[k], &rows[k % 4]);
                        }
                        f.row_reduce(&acc, &mut out);
                        black_box(&out);
                    }
                },
            );
        }
    }
    // Sparse scatter (Acc::sub_tail): 256 rows of 64 entries over 8192
    // columns; eager Barrett per entry against the lazy accumulator.
    let n_cols = 8192usize;
    let sparse: Vec<(Vec<u32>, Vec<u32>, u32)> = (0..256)
        .map(|_| {
            let mut cols: Vec<u32> = (0..64).map(|_| rng.gen_range(0..n_cols as u32)).collect();
            cols.sort_unstable();
            cols.dedup();
            let vals = cols.iter().map(|_| rng.gen_range(1..p) as u32).collect();
            (cols, vals, rng.gen_range(1..p) as u32)
        })
        .collect();
    let nnz: u64 = sparse.iter().map(|r| r.0.len() as u64).sum();
    {
        let mut vals = vec![0u64; n_cols];
        let mut bits = vec![0u64; n_cols / 64];
        for (cols, v, c) in &sparse {
            tower_sub_tail(tfp, &mut vals, &mut bits, u64::from(*c), cols, v);
        }
        let mut wv = vec![0u64; n_cols];
        for (cols, v, c) in &sparse {
            let nc = f.neg(*c);
            for (&col, &x) in cols.iter().zip(v) {
                f.wide_mul_add(&mut wv[col as usize], nc, x);
            }
        }
        for i in 0..n_cols {
            assert_eq!(f.wide_reduce(&wv[i]), tfp.reduce(vals[i]) as u32, "sparse");
        }
    }
    {
        let sparse = sparse.clone();
        let mut vals = vec![0u64; n_cols];
        let mut bits = vec![0u64; n_cols / 64];
        arm(
            arms,
            "f4_fp_tower",
            prime,
            "scatter_sparse",
            nnz,
            move |it| {
                for _ in 0..it {
                    vals.iter_mut().for_each(|x| *x = 0);
                    for (cols, v, c) in &sparse {
                        tower_sub_tail(tfp, &mut vals, &mut bits, u64::from(*c), cols, v);
                    }
                    black_box(&vals);
                }
            },
        );
    }
    {
        let sparse = sparse.clone();
        let f = f.clone();
        let mut vals = vec![0u64; n_cols];
        let mut bits = vec![0u64; n_cols / 64];
        arm(arms, "Fp32", prime, "scatter_sparse", nnz, move |it| {
            for _ in 0..it {
                vals.iter_mut().for_each(|x| *x = 0);
                for (cols, v, c) in &sparse {
                    f.row_scatter_axpy(&mut vals, f.neg(*c), cols, v);
                    for &col in cols {
                        bits[col as usize / 64] |= 1u64 << (col % 64);
                    }
                }
                black_box(&vals);
            }
        });
    }
    {
        let mut vals = vec![0u64; n_cols];
        let mut bits = vec![0u64; n_cols / 64];
        arm(
            arms,
            "Fp32",
            prime,
            "scatter_sparse_scalar",
            nnz,
            move |it| {
                for _ in 0..it {
                    vals.iter_mut().for_each(|x| *x = 0);
                    for (cols, v, c) in &sparse {
                        let nc = f.neg(*c);
                        for (&col, &x) in cols.iter().zip(v) {
                            let ci = col as usize;
                            f.wide_mul_add(&mut vals[ci], nc, x);
                            bits[ci / 64] |= 1u64 << (ci % 64);
                        }
                    }
                    black_box(&vals);
                }
            },
        );
    }
}

/// ct_bignum and BigUint at a wide prime, against FpWide on the same data.
fn reference_arms<const L: usize>(arms: &mut Vec<Arm>, p: &BigUint, prime: &str, seed: u64) {
    let ctx = MontgomeryContext::<L>::new(Uint::<L>::from_biguint(p)).unwrap();
    let mut rng = StdRng::seed_from_u64(seed);
    let xs: Vec<BigUint> = (0..N_ARR).map(|_| rng.gen_biguint_below(p)).collect();
    let ys: Vec<BigUint> = (0..N_ARR).map(|_| rng.gen_biguint_below(p)).collect();
    let am: Vec<Uint<L>> = xs
        .iter()
        .map(|x| ctx.to_montgomery(&Uint::from_biguint(x)))
        .collect();
    let bm: Vec<Uint<L>> = ys
        .iter()
        .map(|x| ctx.to_montgomery(&Uint::from_biguint(x)))
        .collect();
    for i in 0..64 {
        let r = ctx
            .from_montgomery(&ctx.mont_mul(&am[i], &bm[i]))
            .to_biguint();
        assert_eq!(r, (&xs[i] * &ys[i]) % p, "ct_bignum check");
    }
    {
        let (ctx, am, bm) = (ctx.clone(), am.clone(), bm.clone());
        let mut out = vec![Uint::<L>([0; L]); N_ARR];
        arm(
            arms,
            "ct_bignum",
            prime,
            "mul_tp",
            N_ARR as u64,
            move |it| {
                for _ in 0..it {
                    for i in 0..N_ARR {
                        out[i] = ctx.mont_mul(black_box(&am[i]), &bm[i]);
                    }
                    black_box(&out);
                }
            },
        );
    }
    {
        let (ctx, x0, y) = (ctx.clone(), am[1], bm[2]);
        arm(arms, "ct_bignum", prime, "mul_lat", 1024, move |it| {
            let mut x = black_box(x0);
            for _ in 0..it * 1024 {
                x = ctx.mont_mul(&x, &y);
            }
            black_box(x);
        });
    }
    {
        let (ctx, x0) = (ctx.clone(), am[3]);
        arm(arms, "ct_bignum", prime, "sqr_lat", 1024, move |it| {
            let mut x = black_box(x0);
            for _ in 0..it * 1024 {
                x = ctx.mont_sqr(&x);
            }
            black_box(x);
        });
    }
    {
        let (p, xs, ys) = (p.clone(), xs.clone(), ys.clone());
        arm(arms, "BigUint", prime, "mul_tp", 256, move |it| {
            for _ in 0..it {
                for i in 0..256 {
                    black_box((black_box(&xs[i]) * &ys[i]) % &p);
                }
            }
        });
    }
    {
        let p = p.clone();
        let e = &p - BigUint::from(2u32);
        let xs: Vec<BigUint> = xs[..8].to_vec();
        arm(arms, "BigUint", prime, "inv", 8, move |it| {
            for _ in 0..it {
                for x in &xs {
                    black_box(x.modpow(&e, &p));
                }
            }
        });
    }
}

fn prime_fast_arms(arms: &mut Vec<Arm>, p: u64, prime: &str, seed: u64) {
    let f = prime_fast::Fp64::new(p);
    let mut rng = StdRng::seed_from_u64(seed);
    let a: Vec<u64> = (0..N_ARR).map(|_| f.to_mont(rng.gen_range(1..p))).collect();
    let b: Vec<u64> = (0..N_ARR).map(|_| f.to_mont(rng.gen_range(1..p))).collect();
    {
        let (a, b) = (a.clone(), b.clone());
        let mut out = vec![0u64; N_ARR];
        arm(
            arms,
            "prime_fast::Fp64",
            prime,
            "mul_tp",
            N_ARR as u64,
            move |it| {
                for _ in 0..it {
                    for i in 0..N_ARR {
                        out[i] = f.mul(black_box(a[i]), b[i]);
                    }
                    black_box(&out);
                }
            },
        );
    }
    {
        let (x0, y) = (a[1], b[2]);
        arm(
            arms,
            "prime_fast::Fp64",
            prime,
            "mul_lat",
            1024,
            move |it| {
                let mut x = black_box(x0);
                for _ in 0..it * 1024 {
                    x = f.mul(x, y);
                }
                black_box(x);
            },
        );
    }
    {
        let a: Vec<u64> = a[..256].to_vec();
        arm(arms, "prime_fast::Fp64", prime, "inv", 256, move |it| {
            for _ in 0..it {
                for &x in &a {
                    black_box(f.inv(black_box(x)));
                }
            }
        });
    }
}

/// Comparison arms for the algorithm choices inside `FpWide`: the Fermat
/// and plain binary inverses the binary-GCD inverse replaced, Euler's
/// criterion against the binary Jacobi symbol, and plain Tonelli–Shanks
/// against the windowed one (when `s ≥ 16`).
fn choice_arms<const N: usize, const S: bool>(
    arms: &mut Vec<Arm>,
    f: FpWide<N, S>,
    prime: &str,
    seed: u64,
) {
    let backend = f.name().to_string();
    let p = f.modulus().clone();
    let mut rng = StdRng::seed_from_u64(seed);
    let xs: Vec<[u64; N]> = (0..32)
        .map(|_| f.from_biguint(&rng.gen_biguint_below(&p)))
        .filter(|&x| !f.is_zero(x))
        .collect();
    for &x in &xs {
        assert_eq!(f.inv_fermat(x), f.inv(x));
        assert_eq!(f.inv_plain_binary(x), f.inv(x));
        assert_eq!(f.legendre_euler(x), f.legendre(x));
    }
    let n = xs.len() as u64;
    {
        let (f, xs) = (f.clone(), xs.clone());
        arm(arms, &backend, prime, "inv_fermat", n, move |it| {
            for _ in 0..it {
                for &x in &xs {
                    black_box(f.inv_fermat(black_box(x)));
                }
            }
        });
    }
    {
        let (f, xs) = (f.clone(), xs.clone());
        arm(arms, &backend, prime, "inv_plain_binary", n, move |it| {
            for _ in 0..it {
                for &x in &xs {
                    black_box(f.inv_plain_binary(black_box(x)));
                }
            }
        });
    }
    {
        let (f, xs) = (f.clone(), xs.clone());
        arm(arms, &backend, prime, "legendre_euler", n, move |it| {
            for _ in 0..it {
                for &x in &xs {
                    black_box(f.legendre_euler(black_box(x)));
                }
            }
        });
    }
    if f.sqrt_pre().chunks.is_some() {
        let sq: Vec<[u64; N]> = xs.iter().map(|&x| f.sqr(x)).collect();
        for &x in &sq {
            let r = sqrt_tonelli_shanks_plain(&f, x).unwrap();
            assert_eq!(f.sqr(r), x);
        }
        arm(arms, &backend, prime, "sqrt_ts_plain", n, move |it| {
            for _ in 0..it {
                for &x in &sq {
                    black_box(sqrt_tonelli_shanks_plain(&f, black_box(x)));
                }
            }
        });
    }
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    let pos = q * (sorted.len() - 1) as f64;
    let (lo, hi) = (pos.floor() as usize, pos.ceil() as usize);
    sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo as f64)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .map(|i| args[i + 1].clone())
    };
    let reps: usize = get("--reps").map_or(7, |s| s.parse().unwrap());
    let target_ms: f64 = get("--target-ms").map_or(30.0, |s| s.parse().unwrap());
    let out_dir = get("--out").unwrap_or_else(|| ".".into());
    let filter = get("--filter");
    std::fs::create_dir_all(&out_dir).unwrap();

    let p192 = hex("fffffffffffffffffffffffffffffffeffffffffffffffff");
    let p224 = hex("ffffffffffffffffffffffffffffffff000000000000000000000001");
    let p256 = hex("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff");
    // The random 200-bit prime of the tests: first prime above 2^199 + 0x1234567.
    let mut p200 = (BigUint::one() << 199u32) + BigUint::from(0x1234567u64);
    while !is_probable_prime(&p200) {
        p200 += 1u32;
    }

    let mut arms: Vec<Arm> = Vec::new();
    let p31 = 2_013_265_921u64; // 15·2^27 + 1, the PKM test prime p₁
    let p32 = 4_294_967_291u64; // largest prime below 2^32
    tower_arms(&mut arms, p31, "p31", 1);
    tower_arms(&mut arms, p32, "p32", 2);
    field_arms(&mut arms, Fp32::new(p31), "p31", 3);
    field_arms(&mut arms, Fp32::new(p32), "p32", 4);
    field_arms(&mut arms, Fp64::new(p31), "p31", 5);
    field_arms(&mut arms, Fp64::new((1 << 61) - 1), "p61", 6);
    field_arms(&mut arms, Fp64::new((1 << 63) - 25), "p63", 7);
    prime_fast_arms(&mut arms, (1 << 61) - 1, "p61", 8);
    field_arms(&mut arms, Fp192::new(&p192), "P-192", 9);
    field_arms(&mut arms, FpP192::new(), "P-192", 10);
    field_arms(&mut arms, Fp256::new(&p200), "p200", 11);
    field_arms(&mut arms, Fp256::new(&p224), "P-224", 12);
    field_arms(&mut arms, Fp224::new(&p200), "p200", 19);
    field_arms(&mut arms, Fp224::new(&p224), "P-224", 20);
    field_arms(&mut arms, Fp256::new(&p256), "P-256", 13);
    reference_arms::<3>(&mut arms, &p192, "P-192", 14);
    reference_arms::<4>(&mut arms, &p224, "P-224", 15);
    reference_arms::<4>(&mut arms, &p256, "P-256", 16);
    choice_arms(&mut arms, Fp192::new(&p192), "P-192", 17);
    choice_arms(&mut arms, Fp256::new(&p224), "P-224", 18);
    if let Some(flt) = &filter {
        arms.retain(|a| a.name.contains(flt.as_str()));
    }
    eprintln!(
        "{} arms, {} reps, target {} ms",
        arms.len(),
        reps,
        target_ms
    );

    // Calibrate: grow iterations until a call takes ≥ target/4, then scale.
    for a in arms.iter_mut() {
        let mut it = 1u64;
        loop {
            let t = Instant::now();
            (a.run)(it);
            let ms = t.elapsed().as_secs_f64() * 1e3;
            if ms >= target_ms / 4.0 || it > 1 << 40 {
                a.iters = ((it as f64) * target_ms / ms.max(1e-6)).ceil().max(1.0) as u64;
                break;
            }
            it *= 4;
        }
    }

    let mut raw = std::fs::File::create(format!("{out_dir}/raw.csv")).unwrap();
    writeln!(
        raw,
        "rep,slot,arm,backend,prime,op,iters,ops_per_iter,ns_per_op"
    )
    .unwrap();
    let n = arms.len();
    let mut samples: Vec<Vec<f64>> = vec![Vec::new(); n];
    for rep in 0..reps {
        for slot in 0..n {
            // Rotate the order each repetition; odd repetitions run reversed.
            let base = (slot + rep * 7) % n;
            let i = if rep % 2 == 1 { n - 1 - base } else { base };
            let a = &mut arms[i];
            let t = Instant::now();
            (a.run)(a.iters);
            let ns = t.elapsed().as_secs_f64() * 1e9 / (a.iters * a.ops) as f64;
            samples[i].push(ns);
            writeln!(
                raw,
                "{rep},{slot},{},{},{},{},{},{},{:.4}",
                a.name, a.backend, a.prime, a.op, a.iters, a.ops, ns
            )
            .unwrap();
        }
        eprintln!("rep {} done", rep + 1);
    }

    let mut sum = std::fs::File::create(format!("{out_dir}/summary.csv")).unwrap();
    writeln!(
        sum,
        "arm,backend,prime,op,n,median_ns,q1_ns,q3_ns,min_ns,max_ns"
    )
    .unwrap();
    for (a, s) in arms.iter().zip(samples.iter_mut()) {
        s.sort_by(|x, y| x.partial_cmp(y).unwrap());
        let line = format!(
            "{},{},{},{},{},{:.4},{:.4},{:.4},{:.4},{:.4}",
            a.name,
            a.backend,
            a.prime,
            a.op,
            s.len(),
            quantile(s, 0.5),
            quantile(s, 0.25),
            quantile(s, 0.75),
            s[0],
            s[s.len() - 1]
        );
        writeln!(sum, "{line}").unwrap();
        println!("{line}");
    }
}
