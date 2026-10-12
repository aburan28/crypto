//! Rust bindings for the GPU root-finding 2-decomposition oracle
//! (`gpu/rootfind/`), plus the CPU reference it is validated against.
//!
//! The oracle answers, for a batch of targets `R_k` and a factor base
//! `F = {F_0 … F_{N−1}}` on `y² = x³ + ax + b`, every
//! `R_k = ε·F_i + τ·F_j` (`j ≥ i`) by the root-finding identity of
//! [`crate::cryptanalysis::prime_fast::collect_relations`]: the roots of
//! Semaev's `S₃(x_R, x_i, X)` are `x(R ± F_i)`, so each `(k, i)` pair is
//! two slope computations and two table probes.  The C/CUDA side
//! implements that per `(target, slice)` thread; this module
//!
//! - exposes it through [`decompose64`] (any odd prime `p < 2^62`, the
//!   [`prime_fast::Fp64`] bound) and [`decompose_p256`] (NIST P-256), both
//!   available only with the `gpu-rootfind` cargo feature, which pulls in
//!   the `gpu-rootfind-sys` path crate (`gpu/rootfind/sys/`) whose build
//!   script compiles the host-emulation backend (no CUDA needed) or links
//!   a CUDA-enabled library when `GPU_ROOTFIND_CUDA_LIB_DIR` is set;
//! - provides [`cpu_reference64`], the identical sweep written directly
//!   on [`prime_fast`] primitives (`Fp64::batch_inv`,
//!   `FastFactorBase::lookup`), so a harness can demand set equality
//!   between the two on every target.
//!
//! Claim hygiene: this is the decomposition *subroutine* of index
//! calculus on public synthetic factor bases.  Nothing here changes the
//! asymptotics of the ECDLP or says anything about deployed curves; every
//! relation the oracle returns is meant to be re-added in the group by
//! the caller ([`verify_relation64`]), and the validation example
//! (`examples/gpu_rootfind_validate.rs`) fails closed if any is not.
//!
//! [`prime_fast`]: crate::cryptanalysis::prime_fast
//! [`prime_fast::Fp64`]: crate::cryptanalysis::prime_fast::Fp64

use std::collections::BTreeSet;
use std::fmt;

use crate::cryptanalysis::prime_fast::{FastCurve, FastFactorBase, Pt};

/// One decomposition `R_target = eps·F_i + tau·F_j`.  `tau == 0` marks
/// the degenerate `R = eps·F_i` (then `j == i`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Relation {
    pub target: u32,
    pub i: u32,
    pub j: u32,
    pub eps: i8,
    pub tau: i8,
}

/// Which backend ran the oracle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Backend {
    /// Let the library pick: CUDA when linked and a device answers,
    /// host emulation otherwise.
    Auto,
    /// The kernel's per-thread function executed on the host.
    Emulate,
    /// A CUDA device.
    Cuda,
}

/// Lookup structure for the factor-base x-coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Lookup {
    /// Hash when it fits the device's shared-memory budget, sorted otherwise.
    Auto,
    /// Sorted key array + binary search (global memory, any size).
    Sorted,
    /// Open-addressing hash (shared memory on the device when small).
    Hash,
}

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub backend: Backend,
    pub lookup: Lookup,
    /// Factor-base slice per thread: 16, 32 or 64 (0 = library default).
    pub chunk: u32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            backend: Backend::Auto,
            lookup: Lookup::Auto,
            chunk: 0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Report {
    pub relations: Vec<Relation>,
    pub backend_used: Backend,
    pub lookup_used: Lookup,
    pub chunk_used: u32,
    /// Oracle wall time as measured by the library (excludes table build).
    pub seconds: f64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    /// Crate built without the `gpu-rootfind` feature.
    FeatureDisabled,
    /// Output buffer too small even after retry; carries `n_found`.
    Overflow(u32),
    BadArgument,
    NoCuda,
    Cuda,
    BadChunk,
    NoMemory,
    Unknown(i32),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::FeatureDisabled => write!(f, "crate built without the `gpu-rootfind` feature"),
            Error::Overflow(n) => write!(f, "output overflow: kernel found {n} relations"),
            Error::BadArgument => write!(f, "bad argument"),
            Error::NoCuda => write!(f, "CUDA backend unavailable"),
            Error::Cuda => write!(f, "CUDA runtime error"),
            Error::BadChunk => write!(f, "chunk must be 16, 32 or 64"),
            Error::NoMemory => write!(f, "out of memory"),
            Error::Unknown(c) => write!(f, "unknown status {c}"),
        }
    }
}

impl std::error::Error for Error {}

// ── FFI ───────────────────────────────────────────────────────────────

/// Raw C ABI, declared (and the native library compiled or linked) by the
/// `gpu-rootfind-sys` path crate in `gpu/rootfind/sys/`.
#[cfg(feature = "gpu-rootfind")]
use gpu_rootfind_sys as ffi;

#[cfg(feature = "gpu-rootfind")]
fn backend_code(b: Backend) -> i32 {
    match b {
        Backend::Auto => 0,
        Backend::Emulate => 1,
        Backend::Cuda => 2,
    }
}

#[cfg(feature = "gpu-rootfind")]
fn backend_from(c: i32) -> Backend {
    match c {
        2 => Backend::Cuda,
        1 => Backend::Emulate,
        _ => Backend::Auto,
    }
}

#[cfg(feature = "gpu-rootfind")]
fn lookup_code(l: Lookup) -> i32 {
    match l {
        Lookup::Auto => 0,
        Lookup::Sorted => 1,
        Lookup::Hash => 2,
    }
}

#[cfg(feature = "gpu-rootfind")]
fn lookup_from(c: i32) -> Lookup {
    match c {
        2 => Lookup::Hash,
        1 => Lookup::Sorted,
        _ => Lookup::Auto,
    }
}

#[cfg(feature = "gpu-rootfind")]
fn status_to_error(rc: i32, found: u32) -> Error {
    match rc {
        ffi::RF_ERR_OVERFLOW => Error::Overflow(found),
        ffi::RF_ERR_BAD_ARG => Error::BadArgument,
        ffi::RF_ERR_NO_CUDA => Error::NoCuda,
        ffi::RF_ERR_CUDA => Error::Cuda,
        ffi::RF_ERR_BAD_CHUNK => Error::BadChunk,
        ffi::RF_ERR_NOMEM => Error::NoMemory,
        other => Error::Unknown(other),
    }
}

/// True when the linked library has a CUDA backend *and* a device answers.
pub fn cuda_available() -> bool {
    #[cfg(feature = "gpu-rootfind")]
    {
        unsafe { ffi::rf_cuda_available() != 0 }
    }
    #[cfg(not(feature = "gpu-rootfind"))]
    {
        false
    }
}

/// Human-readable description of the linked backend.
pub fn backend_info() -> String {
    #[cfg(feature = "gpu-rootfind")]
    {
        let p = unsafe { ffi::rf_backend_info() };
        if p.is_null() {
            return "unknown".into();
        }
        unsafe { std::ffi::CStr::from_ptr(p) }.to_string_lossy().into_owned()
    }
    #[cfg(not(feature = "gpu-rootfind"))]
    {
        "gpu-rootfind feature disabled".into()
    }
}

/// Run the oracle over a 64-bit field.  Coordinates are canonical
/// (`< p`); `fb_x.len() == fb_y.len()`, `tx.len() == ty.len()`.  The
/// output buffer is sized from `n_targets` and grown once on overflow.
#[allow(unused_variables)]
pub fn decompose64(
    p: u64,
    fb_x: &[u64],
    fb_y: &[u64],
    tx: &[u64],
    ty: &[u64],
    opts: Options,
) -> Result<Report, Error> {
    #[cfg(not(feature = "gpu-rootfind"))]
    {
        Err(Error::FeatureDisabled)
    }
    #[cfg(feature = "gpu-rootfind")]
    {
        if fb_x.len() != fb_y.len() || tx.len() != ty.len() || fb_x.is_empty() || tx.is_empty() {
            return Err(Error::BadArgument);
        }
        let call = |cap: u32| -> (i32, ffi::RfResult, Vec<ffi::RfRelation>) {
            let mut out = vec![ffi::RfRelation::default(); cap as usize];
            let opt = ffi::RfOptions {
                backend: backend_code(opts.backend),
                lookup: lookup_code(opts.lookup),
                chunk: opts.chunk,
            };
            let mut res = ffi::RfResult::default();
            let rc = unsafe {
                ffi::rf64_decompose(
                    p,
                    fb_x.as_ptr(),
                    fb_y.as_ptr(),
                    fb_x.len() as u32,
                    tx.as_ptr(),
                    ty.as_ptr(),
                    tx.len() as u32,
                    out.as_mut_ptr(),
                    cap,
                    &opt,
                    &mut res,
                )
            };
            (rc, res, out)
        };
        finish(call, tx.len())
    }
}

/// Run the oracle on NIST P-256.  Each coordinate is 4 little-endian
/// `u64` limbs (`fb_x.len() == 4 * n_fb`, etc.), canonical (`< p`).
#[allow(unused_variables)]
pub fn decompose_p256(
    fb_x: &[u64],
    fb_y: &[u64],
    tx: &[u64],
    ty: &[u64],
    opts: Options,
) -> Result<Report, Error> {
    #[cfg(not(feature = "gpu-rootfind"))]
    {
        Err(Error::FeatureDisabled)
    }
    #[cfg(feature = "gpu-rootfind")]
    {
        if fb_x.len() != fb_y.len()
            || tx.len() != ty.len()
            || fb_x.len() % 4 != 0
            || tx.len() % 4 != 0
            || fb_x.is_empty()
            || tx.is_empty()
        {
            return Err(Error::BadArgument);
        }
        let n_fb = (fb_x.len() / 4) as u32;
        let n_t = (tx.len() / 4) as u32;
        let call = |cap: u32| -> (i32, ffi::RfResult, Vec<ffi::RfRelation>) {
            let mut out = vec![ffi::RfRelation::default(); cap as usize];
            let opt = ffi::RfOptions {
                backend: backend_code(opts.backend),
                lookup: lookup_code(opts.lookup),
                chunk: if opts.chunk == 0 { 16 } else { opts.chunk },
            };
            let mut res = ffi::RfResult::default();
            let rc = unsafe {
                ffi::rf256_decompose(
                    fb_x.as_ptr(),
                    fb_y.as_ptr(),
                    n_fb,
                    tx.as_ptr(),
                    ty.as_ptr(),
                    n_t,
                    out.as_mut_ptr(),
                    cap,
                    &opt,
                    &mut res,
                )
            };
            (rc, res, out)
        };
        finish(call, n_t as usize)
    }
}

#[cfg(feature = "gpu-rootfind")]
fn finish<F>(call: F, n_targets: usize) -> Result<Report, Error>
where
    F: Fn(u32) -> (i32, ffi::RfResult, Vec<ffi::RfRelation>),
{
    // Expected hits are sparse; start with 4 per target plus slack and
    // retry once with the exact count if the kernel reports overflow.
    let cap0 = (4 * n_targets + 1024) as u32;
    let (mut rc, mut res, mut out) = call(cap0);
    if rc == ffi::RF_ERR_OVERFLOW {
        let (rc2, res2, out2) = call(res.n_found + 64);
        rc = rc2;
        res = res2;
        out = out2;
    }
    if rc != ffi::RF_OK {
        return Err(status_to_error(rc, res.n_found));
    }
    let relations = out[..res.n_written as usize]
        .iter()
        .map(|r| Relation {
            target: r.target,
            i: r.i,
            j: r.j,
            eps: r.eps,
            tau: r.tau,
        })
        .collect();
    Ok(Report {
        relations,
        backend_used: backend_from(res.backend_used),
        lookup_used: lookup_from(res.lookup_used),
        chunk_used: res.chunk_used,
        seconds: res.seconds,
    })
}

// ── CPU reference (prime_fast primitives) ─────────────────────────────

/// The same sweep as the kernel, written on [`prime_fast`] primitives:
/// one `Fp64::batch_inv` per target over `x_R − x_i`, both slopes, and a
/// `FastFactorBase::lookup` per `(i, sign)`.  This is the inner loop of
/// `collect_relations` with the output turned into index-level
/// relations, so agreement with the GPU is agreement with the production
/// CPU oracle.
///
/// [`prime_fast`]: crate::cryptanalysis::prime_fast
pub fn cpu_reference64(curve: &FastCurve, fb: &FastFactorBase, targets: &[Pt]) -> BTreeSet<Relation> {
    let f = curve.f;
    let n = fb.len();
    let mut d = vec![0u64; n];
    let mut degenerate = vec![false; n];
    let mut scratch = Vec::with_capacity(n);
    let mut out = BTreeSet::new();
    for (k, &r) in targets.iter().enumerate() {
        for i in 0..n {
            let dd = f.sub(r.x, fb.pts[i].x);
            degenerate[i] = dd == 0;
            d[i] = if dd == 0 { f.one } else { dd };
        }
        f.batch_inv(&mut d, &mut scratch);
        for i in 0..n {
            let fi = fb.pts[i];
            if degenerate[i] {
                let eps = if r.y == fi.y { 1i8 } else { -1 };
                out.insert(Relation {
                    target: k as u32,
                    i: i as u32,
                    j: i as u32,
                    eps,
                    tau: 0,
                });
                continue;
            }
            let inv = d[i];
            let lam_minus = f.mul(f.add(r.y, fi.y), inv);
            let lam_plus = f.mul(f.sub(r.y, fi.y), inv);
            for (eps, lam) in [(1i8, lam_minus), (-1i8, lam_plus)] {
                let x3 = f.sub(f.sub(f.sqr(lam), r.x), fi.x);
                if let Some(j) = fb.lookup(f.from_mont(x3)) {
                    if j >= i {
                        let y3 = f.sub(f.mul(lam, f.sub(r.x, x3)), r.y);
                        let tau = if y3 == fb.pts[j].y { 1i8 } else { -1 };
                        out.insert(Relation {
                            target: k as u32,
                            i: i as u32,
                            j: j as u32,
                            eps,
                            tau,
                        });
                    }
                }
            }
        }
    }
    out
}

/// Fail-closed check: re-add `eps·F_i + tau·F_j` in the group and compare
/// with the target.
pub fn verify_relation64(curve: &FastCurve, fb: &FastFactorBase, target: Pt, rel: &Relation) -> bool {
    let (i, j) = (rel.i as usize, rel.j as usize);
    if i >= fb.len() || j >= fb.len() || j < i {
        return false;
    }
    let sign = |p: Pt, s: i8| if s < 0 { curve.neg(p) } else { p };
    let pi = sign(fb.pts[i], rel.eps);
    let sum = if rel.tau == 0 {
        if i != j {
            return false;
        }
        Some(pi)
    } else {
        curve.add(Some(pi), Some(sign(fb.pts[j], rel.tau)))
    };
    sum == Some(target)
}

/// Canonical coordinate arrays for a factor base, in the layout the FFI
/// takes.
pub fn factor_base_canonical(curve: &FastCurve, fb: &FastFactorBase) -> (Vec<u64>, Vec<u64>) {
    let mut xs = Vec::with_capacity(fb.len());
    let mut ys = Vec::with_capacity(fb.len());
    for &p in &fb.pts {
        let (x, y) = curve.canonical(p);
        xs.push(x);
        ys.push(y);
    }
    (xs, ys)
}

/// Canonical coordinate arrays for a target batch.
pub fn targets_canonical(curve: &FastCurve, targets: &[Pt]) -> (Vec<u64>, Vec<u64>) {
    let mut xs = Vec::with_capacity(targets.len());
    let mut ys = Vec::with_capacity(targets.len());
    for &p in targets {
        let (x, y) = curve.canonical(p);
        xs.push(x);
        ys.push(y);
    }
    (xs, ys)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::prime_fast::find_a3_curve;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn instance(bits: u32, fb_size: usize, n_random: usize) -> (FastCurve, FastFactorBase, Vec<Pt>, usize) {
        let (curve, _g) = find_a3_curve(bits, 3).expect("curve");
        let fb = FastFactorBase::build(&curve, fb_size);
        let mut rng = StdRng::seed_from_u64(7);
        let mut targets = Vec::new();
        for _ in 0..40 {
            let i = rng.gen_range(0..fb.len());
            let j = rng.gen_range(0..fb.len());
            let pi = if rng.gen_bool(0.5) { curve.neg(fb.pts[i]) } else { fb.pts[i] };
            let pj = if rng.gen_bool(0.5) { curve.neg(fb.pts[j]) } else { fb.pts[j] };
            if let Some(r) = curve.add(Some(pi), Some(pj)) {
                targets.push(r);
            }
        }
        let n_ctl = targets.len();
        while targets.len() < n_ctl + n_random {
            let x = rng.gen_range(1..curve.f.p);
            if let Some(p) = curve.lift_x(x) {
                targets.push(p);
            }
        }
        targets.push(fb.pts[3]);
        (curve, fb, targets, n_ctl)
    }

    #[test]
    fn cpu_reference_relations_re_add_and_cover_controls() {
        let (curve, fb, targets, n_ctl) = instance(40, 200, 100);
        let rels = cpu_reference64(&curve, &fb, &targets);
        for r in &rels {
            assert!(verify_relation64(&curve, &fb, targets[r.target as usize], r), "{r:?}");
        }
        for k in 0..n_ctl {
            assert!(rels.iter().any(|r| r.target as usize == k), "control {k} not decomposed");
        }
        let last = (targets.len() - 1) as u32;
        assert!(rels.contains(&Relation {
            target: last,
            i: 3,
            j: 3,
            eps: 1,
            tau: 0
        }));
    }

    #[cfg(feature = "gpu-rootfind")]
    #[test]
    fn emulated_oracle_matches_cpu_reference_bit_for_bit() {
        for bits in [40u32, 61] {
            let (curve, fb, targets, _) = instance(bits, 300, 200);
            let (fbx, fby) = factor_base_canonical(&curve, &fb);
            let (tx, ty) = targets_canonical(&curve, &targets);
            let reference = cpu_reference64(&curve, &fb, &targets);
            for lookup in [Lookup::Sorted, Lookup::Hash] {
                for chunk in [16u32, 32, 64] {
                    let rep = decompose64(
                        curve.f.p,
                        &fbx,
                        &fby,
                        &tx,
                        &ty,
                        Options {
                            backend: Backend::Emulate,
                            lookup,
                            chunk,
                        },
                    )
                    .expect("oracle");
                    let got: BTreeSet<Relation> = rep.relations.iter().copied().collect();
                    assert_eq!(got, reference, "bits={bits} lookup={lookup:?} chunk={chunk}");
                    for r in &rep.relations {
                        assert!(verify_relation64(&curve, &fb, targets[r.target as usize], r));
                    }
                }
            }
        }
    }

    #[cfg(feature = "gpu-rootfind")]
    #[test]
    fn emulated_oracle_rejects_bad_input() {
        let (curve, fb, targets, _) = instance(40, 50, 10);
        let (fbx, fby) = factor_base_canonical(&curve, &fb);
        let (mut tx, ty) = targets_canonical(&curve, &targets);
        tx[0] = curve.f.p; // not canonical
        assert_eq!(
            decompose64(curve.f.p, &fbx, &fby, &tx, &ty, Options::default()).unwrap_err(),
            Error::BadArgument
        );
        let (tx, ty) = targets_canonical(&curve, &targets);
        assert_eq!(
            decompose64(
                curve.f.p,
                &fbx,
                &fby,
                &tx,
                &ty,
                Options {
                    chunk: 24,
                    ..Options::default()
                }
            )
            .unwrap_err(),
            Error::BadChunk
        );
        assert_eq!(
            decompose64(curve.f.p + 1, &fbx, &fby, &tx, &ty, Options::default()).unwrap_err(),
            Error::BadArgument
        );
    }

    #[cfg(feature = "gpu-rootfind")]
    #[test]
    fn emulated_p256_oracle_matches_crate_field() {
        use crate::ecc::p256_field::P256FieldElement;
        use num_bigint::BigUint;
        // a = -3, b = P-256's b; factor base of small x lifted with the crate field.
        let p = BigUint::parse_bytes(
            b"ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let b = P256FieldElement::from_biguint(
            &BigUint::parse_bytes(b"5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b", 16).unwrap(),
        );
        let a = P256FieldElement::from_biguint(&(&p - BigUint::from(3u32)));
        let rhs = |x: &P256FieldElement| x.sqr().mul(x).add(&a.mul(x)).add(&b);
        let sqrt = |v: &P256FieldElement| -> Option<P256FieldElement> {
            // p = 3 mod 4: v^((p+1)/4)
            let e: BigUint = (&p + BigUint::from(1u32)) >> 2u32;
            let mut acc = P256FieldElement::ONE;
            for i in (0..256u64).rev() {
                acc = acc.sqr();
                if e.bit(i) {
                    acc = acc.mul(v);
                }
            }
            if bool::from(acc.sqr().ct_eq(v)) {
                Some(acc)
            } else {
                None
            }
        };
        let mut fb: Vec<(P256FieldElement, P256FieldElement)> = Vec::new();
        let mut x = 1u64;
        while fb.len() < 48 {
            let xe = P256FieldElement::from_biguint(&BigUint::from(x));
            if let Some(y) = sqrt(&rhs(&xe)) {
                fb.push((xe, y));
            }
            x += 1;
        }
        let limbs = |e: &P256FieldElement| e.to_canonical().0;
        let mut fbx = Vec::new();
        let mut fby = Vec::new();
        for (xe, ye) in &fb {
            fbx.extend_from_slice(&limbs(xe));
            fby.extend_from_slice(&limbs(ye));
        }
        let add = |p1: &(P256FieldElement, P256FieldElement), p2: &(P256FieldElement, P256FieldElement)| {
            let lam = p2.1.sub(&p1.1).mul(&p2.0.sub(&p1.0).inv());
            let x3 = lam.sqr().sub(&p1.0).sub(&p2.0);
            let y3 = lam.mul(&p1.0.sub(&x3)).sub(&p1.1);
            (x3, y3)
        };
        let mut rng = StdRng::seed_from_u64(11);
        let mut tx = Vec::new();
        let mut ty = Vec::new();
        let mut expect = Vec::new();
        for k in 0..16u32 {
            let i = rng.gen_range(0..fb.len() - 1);
            let j = rng.gen_range(i + 1..fb.len());
            let neg = rng.gen_bool(0.5);
            let pj = if neg { (fb[j].0, fb[j].1.neg()) } else { fb[j] };
            let r = add(&fb[i], &pj);
            tx.extend_from_slice(&limbs(&r.0));
            ty.extend_from_slice(&limbs(&r.1));
            expect.push(Relation {
                target: k,
                i: i as u32,
                j: j as u32,
                eps: 1,
                tau: if neg { -1 } else { 1 },
            });
        }
        for lookup in [Lookup::Sorted, Lookup::Hash] {
            let rep = decompose_p256(
                &fbx,
                &fby,
                &tx,
                &ty,
                Options {
                    backend: Backend::Emulate,
                    lookup,
                    chunk: 16,
                },
            )
            .expect("p256 oracle");
            for e in &expect {
                assert!(rep.relations.contains(e), "{lookup:?}: missing {e:?}");
            }
        }
    }
}
