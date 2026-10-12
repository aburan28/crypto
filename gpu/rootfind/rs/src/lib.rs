//! Rust bindings for the GPU root-finding 2-decomposition oracle
//! (`gpu/rootfind/`): the raw C ABI of `rootfind_api.h` and safe wrappers
//! [`decompose64`] (any odd prime `p < 2^62`, the `prime_fast::Fp64`
//! bound) and [`decompose_p256`] (NIST P-256).
//!
//! The native library is compiled by this crate's `build.rs` — the
//! host-emulation backend by default (no CUDA needed), or a prebuilt CUDA
//! library when `GPU_ROOTFIND_CUDA_LIB_DIR` is set.  The pure-Rust half
//! (the [`Relation`] type, the CPU reference `cpu_reference64` written on
//! `prime_fast` primitives, and the group re-add check) lives in the root
//! crate at `crypto_lib::cryptanalysis::gpu_rootfind` and is re-exported
//! here so a harness can demand set equality between the two.
//!
//! Claim hygiene: this is the decomposition *subroutine* of index
//! calculus on public synthetic factor bases.  Nothing here changes the
//! asymptotics of the ECDLP or says anything about deployed curves; every
//! relation the oracle returns is meant to be re-added in the group by the
//! caller (`verify_relation64`), and `gpu_rootfind_validate` fails closed
//! if any is not.

use std::fmt;

pub use crypto_lib::cryptanalysis::gpu_rootfind::{
    cpu_reference64, factor_base_canonical, synthetic_instance, targets_canonical,
    verify_relation64, Relation,
};

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
    /// Factor-base slice per thread: 16, 32 or 64.  `0` is passed through
    /// to the library, whose default is per field: 32 for the 64-bit field
    /// (`RF_DEFAULT_CHUNK_64`) and 16 for P-256 (`RF_DEFAULT_CHUNK_256`,
    /// register pressure of four-limb elements).  `Report::chunk_used`
    /// reports what ran.
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

// ── raw C ABI (rootfind_api.h) ────────────────────────────────────────

pub mod ffi {
    //! Layouts mirror `rootfind_api.h` exactly: `rf_relation` is 16
    //! bytes, `rf_options` / `rf_result` are plain structs of `int` /
    //! `uint32_t` / `double`.
    use std::os::raw::{c_char, c_int};

    #[repr(C)]
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
    pub struct RfRelation {
        pub target: u32,
        pub i: u32,
        pub j: u32,
        pub eps: i8,
        pub tau: i8,
        pub pad: [u8; 2],
    }

    #[repr(C)]
    #[derive(Clone, Copy, Debug)]
    pub struct RfOptions {
        pub backend: c_int,
        pub lookup: c_int,
        pub chunk: u32,
    }

    #[repr(C)]
    #[derive(Clone, Copy, Debug, Default)]
    pub struct RfResult {
        pub n_written: u32,
        pub n_found: u32,
        pub backend_used: c_int,
        pub lookup_used: c_int,
        pub chunk_used: u32,
        pub seconds: f64,
    }

    pub const RF_OK: c_int = 0;
    pub const RF_ERR_OVERFLOW: c_int = 1;
    pub const RF_ERR_BAD_ARG: c_int = -1;
    pub const RF_ERR_NO_CUDA: c_int = -2;
    pub const RF_ERR_CUDA: c_int = -3;
    pub const RF_ERR_BAD_CHUNK: c_int = -4;
    pub const RF_ERR_NOMEM: c_int = -5;

    extern "C" {
        pub fn rf64_decompose(
            p: u64,
            fb_x: *const u64,
            fb_y: *const u64,
            n_fb: u32,
            tx: *const u64,
            ty: *const u64,
            n_targets: u32,
            out: *mut RfRelation,
            cap: u32,
            opt: *const RfOptions,
            res: *mut RfResult,
        ) -> c_int;
        pub fn rf256_decompose(
            fb_x: *const u64,
            fb_y: *const u64,
            n_fb: u32,
            tx: *const u64,
            ty: *const u64,
            n_targets: u32,
            out: *mut RfRelation,
            cap: u32,
            opt: *const RfOptions,
            res: *mut RfResult,
        ) -> c_int;
        pub fn rf_cuda_available() -> c_int;
        pub fn rf_backend_info() -> *const c_char;
        pub fn rf_strerror(status: c_int) -> *const c_char;
    }
}

fn backend_code(b: Backend) -> i32 {
    match b {
        Backend::Auto => 0,
        Backend::Emulate => 1,
        Backend::Cuda => 2,
    }
}

fn backend_from(c: i32) -> Backend {
    match c {
        2 => Backend::Cuda,
        1 => Backend::Emulate,
        _ => Backend::Auto,
    }
}

fn lookup_code(l: Lookup) -> i32 {
    match l {
        Lookup::Auto => 0,
        Lookup::Sorted => 1,
        Lookup::Hash => 2,
    }
}

fn lookup_from(c: i32) -> Lookup {
    match c {
        2 => Lookup::Hash,
        1 => Lookup::Sorted,
        _ => Lookup::Auto,
    }
}

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
    // SAFETY: no arguments, returns an int.
    unsafe { ffi::rf_cuda_available() != 0 }
}

/// Human-readable description of the linked backend.
pub fn backend_info() -> String {
    // SAFETY: returns a pointer to static storage inside the library.
    let p = unsafe { ffi::rf_backend_info() };
    if p.is_null() {
        return "unknown".into();
    }
    unsafe { std::ffi::CStr::from_ptr(p) }
        .to_string_lossy()
        .into_owned()
}

/// Run the oracle over a 64-bit field.  Coordinates are canonical
/// (`< p`); `fb_x.len() == fb_y.len()`, `tx.len() == ty.len()`.  The
/// output buffer is sized from `n_targets` and grown once on overflow.
pub fn decompose64(
    p: u64,
    fb_x: &[u64],
    fb_y: &[u64],
    tx: &[u64],
    ty: &[u64],
    opts: Options,
) -> Result<Report, Error> {
    if fb_x.len() != fb_y.len() || tx.len() != ty.len() || fb_x.is_empty() || tx.is_empty() {
        return Err(Error::BadArgument);
    }
    let n_fb = u32::try_from(fb_x.len()).map_err(|_| Error::BadArgument)?;
    let n_t = u32::try_from(tx.len()).map_err(|_| Error::BadArgument)?;
    let call = |cap: u32| -> (i32, ffi::RfResult, Vec<ffi::RfRelation>) {
        let mut out = vec![ffi::RfRelation::default(); cap as usize];
        let opt = ffi::RfOptions {
            backend: backend_code(opts.backend),
            lookup: lookup_code(opts.lookup),
            chunk: opts.chunk,
        };
        let mut res = ffi::RfResult::default();
        // SAFETY: every pointer/length pair describes a live slice; `out`
        // has exactly `cap` elements; `opt`/`res` outlive the call.
        let rc = unsafe {
            ffi::rf64_decompose(
                p,
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
    finish(call, tx.len())
}

/// Run the oracle on NIST P-256.  Each coordinate is 4 little-endian
/// `u64` limbs (`fb_x.len() == 4 * n_fb`, etc.), canonical (`< p`).
pub fn decompose_p256(
    fb_x: &[u64],
    fb_y: &[u64],
    tx: &[u64],
    ty: &[u64],
    opts: Options,
) -> Result<Report, Error> {
    if fb_x.len() != fb_y.len()
        || tx.len() != ty.len()
        || fb_x.len() % 4 != 0
        || tx.len() % 4 != 0
        || fb_x.is_empty()
        || tx.is_empty()
    {
        return Err(Error::BadArgument);
    }
    let n_fb = u32::try_from(fb_x.len() / 4).map_err(|_| Error::BadArgument)?;
    let n_t = u32::try_from(tx.len() / 4).map_err(|_| Error::BadArgument)?;
    let call = |cap: u32| -> (i32, ffi::RfResult, Vec<ffi::RfRelation>) {
        let mut out = vec![ffi::RfRelation::default(); cap as usize];
        let opt = ffi::RfOptions {
            backend: backend_code(opts.backend),
            lookup: lookup_code(opts.lookup),
            // 0 = library default (RF_DEFAULT_CHUNK_256 = 16 for this field).
            chunk: opts.chunk,
        };
        let mut res = ffi::RfResult::default();
        // SAFETY: as in decompose64; limb counts were checked above.
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

fn finish<F>(call: F, n_targets: usize) -> Result<Report, Error>
where
    F: Fn(u32) -> (i32, ffi::RfResult, Vec<ffi::RfRelation>),
{
    // Expected hits are sparse; start with 4 per target plus slack and
    // retry once with the exact count if the kernel reports overflow.
    // Near 2^30 targets that guess no longer fits in the 32-bit capacity,
    // while the relation counter wraps only after 2^32 hits.  Probe with
    // an empty buffer and size the retry from `n_found`.  `u32::MAX` here
    // would allocate ~64 GiB of relations before the oracle runs.
    let cap0 = u32::try_from(4 * n_targets + 1024).unwrap_or(0);
    let (mut rc, mut res, mut out) = call(cap0);
    if rc == ffi::RF_ERR_OVERFLOW {
        let (rc2, res2, out2) = call(res.n_found.saturating_add(64));
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn ffi_layouts_match_the_c_header() {
        assert_eq!(std::mem::size_of::<ffi::RfRelation>(), 16);
        assert_eq!(std::mem::size_of::<ffi::RfOptions>(), 12);
        assert_eq!(std::mem::size_of::<ffi::RfResult>(), 32);
        let ok = unsafe { std::ffi::CStr::from_ptr(ffi::rf_strerror(ffi::RF_OK)) };
        assert_eq!(ok.to_str().unwrap(), "ok");
        assert!(!backend_info().is_empty());
    }

    #[test]
    fn emulated_oracle_matches_cpu_reference_bit_for_bit() {
        for bits in [40u32, 61] {
            let (curve, fb, targets, _) = synthetic_instance(bits, 300, 200, 7);
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
                    assert_eq!(
                        got, reference,
                        "bits={bits} lookup={lookup:?} chunk={chunk}"
                    );
                    for r in &rep.relations {
                        assert!(verify_relation64(
                            &curve,
                            &fb,
                            targets[r.target as usize],
                            r
                        ));
                    }
                }
            }
        }
    }

    #[test]
    fn emulated_oracle_rejects_bad_input() {
        let (curve, fb, targets, _) = synthetic_instance(40, 50, 10, 7);
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
        // A CUDA request on a machine without the CUDA backend fails closed.
        if !cuda_available() {
            assert_eq!(
                decompose64(
                    curve.f.p,
                    &fbx,
                    &fby,
                    &tx,
                    &ty,
                    Options {
                        backend: Backend::Cuda,
                        ..Options::default()
                    }
                )
                .unwrap_err(),
                Error::NoCuda
            );
        }
    }

    #[test]
    fn emulated_p256_oracle_matches_crate_field() {
        use crypto_lib::ecc::p256_field::P256FieldElement;
        use num_bigint::BigUint;
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        // a = -3, b = P-256's b; factor base of small x lifted with the crate field.
        let p = BigUint::parse_bytes(
            b"ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let b = P256FieldElement::from_biguint(
            &BigUint::parse_bytes(
                b"5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
                16,
            )
            .unwrap(),
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
        let add = |p1: &(P256FieldElement, P256FieldElement),
                   p2: &(P256FieldElement, P256FieldElement)| {
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
            for chunk in [0u32, 16, 32] {
                let rep = decompose_p256(
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
                .expect("p256 oracle");
                // chunk 0 resolves to the library's per-field default
                // (RF_DEFAULT_CHUNK_256 = 16), not to the 64-bit field's 32.
                let want = if chunk == 0 { 16 } else { chunk };
                assert_eq!(rep.chunk_used, want, "{lookup:?}: chunk_used");
                for e in &expect {
                    assert!(rep.relations.contains(e), "{lookup:?}: missing {e:?}");
                }
            }
        }
    }

    #[test]
    fn oversized_batch_sizes_the_buffer_from_n_found() {
        // 4 * n + 1024 == 2^32, so the initial guess does not fit in u32.
        // A probe must not request u32::MAX relations; the retry uses the
        // count the kernel reports.
        let n_targets = (u32::MAX as usize - 1024) / 4 + 1;
        assert!(u32::try_from(4 * n_targets + 1024).is_err());
        let caps = std::cell::RefCell::new(Vec::new());
        let report = finish(
            |cap| {
                caps.borrow_mut().push(cap);
                assert!(cap < 1_000_000, "refusing a huge relation buffer ({cap})");
                let mut res = ffi::RfResult::default();
                res.n_found = 7;
                res.n_written = if cap >= 7 { 7 } else { 0 };
                res.backend_used = 1;
                res.lookup_used = 1;
                res.chunk_used = 32;
                let rc = if (cap as u64) < u64::from(res.n_found) {
                    ffi::RF_ERR_OVERFLOW
                } else {
                    ffi::RF_OK
                };
                let mut out = vec![ffi::RfRelation::default(); cap as usize];
                if cap >= 7 {
                    for (k, slot) in out.iter_mut().take(7).enumerate() {
                        slot.target = k as u32;
                        slot.i = 1;
                        slot.j = 2;
                        slot.eps = 1;
                        slot.tau = 1;
                    }
                }
                (rc, res, out)
            },
            n_targets,
        )
        .expect("sized from n_found");
        assert_eq!(caps.into_inner(), vec![0, 7 + 64]);
        assert_eq!(report.relations.len(), 7);
        assert_eq!(report.relations[0].target, 0);
        assert_eq!(report.chunk_used, 32);
        assert_eq!(report.backend_used, Backend::Emulate);
    }

    #[test]
    fn chunk_zero_resolves_per_field_default() {
        let (curve, fb, targets, _) = synthetic_instance(40, 64, 8, 3);
        let (fx, fy) = factor_base_canonical(&curve, &fb);
        let (tx, ty) = targets_canonical(&curve, &targets);
        let rep = decompose64(
            curve.f.p,
            &fx,
            &fy,
            &tx,
            &ty,
            Options {
                backend: Backend::Emulate,
                ..Options::default()
            },
        )
        .expect("oracle");
        assert_eq!(
            rep.chunk_used, 32,
            "64-bit field default is RF_DEFAULT_CHUNK_64"
        );
    }
}
