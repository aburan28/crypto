//! Calibration measurements on the host this runs on: the effective clock,
//! cycles per field operation, read bandwidth per memory level, the rho
//! iteration of [`super::kernel`] across fields and batch sizes, and its
//! scaling over threads.
//!
//! Every measurement is a set of *arms* run for `reps` repetitions; each
//! repetition runs every arm of its group once, in an order rotated by the
//! repetition index, so drift on a shared host spreads over all arms.  An
//! arm's inner loop count is calibrated once to take about `target_ms`
//! and then frozen.  Summaries are medians over repetitions.

use std::collections::BTreeMap;
use std::hint::black_box;
use std::sync::Barrier;
use std::time::Instant;

use num_bigint::{BigUint, RandBigInt};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde::{Deserialize, Serialize};

use super::kernel::{self, Arith, Counted, LowWord, OpTally, Walks};
use crate::cryptanalysis::ecc2k130_guard::is_probable_prime;
use crate::cryptanalysis::fp_backend::PrimeField;

pub const MEASUREMENTS_SCHEMA: &str = "hwcost.measurements.v1";

/// Distinguished-point test mask of the timed kernel (`2⁻¹⁶`); a hit is
/// counted, the walk is not restarted.
pub const DP_MASK: u64 = (1 << 16) - 1;

/// The fields measured.  `w1 … w4` are one implementation (`FpWide<N>`,
/// Montgomery CIOS) at 1–4 limbs, with primes filling every limb, so cost
/// can be fitted against the limb count; `fp64` is the dedicated one-word
/// backend and `fp64bf` the same field and representation through
/// [`kernel::Fp64Bf`] (no data-dependent branch in subtraction or
/// reduction); `pf64` is the production prime-field rho engine's
/// arithmetic ([`kernel::PrimeFast`], `p < 2⁶²`).
#[derive(Clone, Copy, Debug)]
pub struct FieldSpec {
    pub name: &'static str,
    pub backend: &'static str,
    pub limbs: usize,
    pub prime_hex: &'static str,
    pub prime_label: &'static str,
}

pub const FIELDS: [FieldSpec; 8] = [
    FieldSpec {
        name: "fp64",
        backend: "Fp64",
        limbs: 1,
        prime_hex: "7fffffffffffffe7",
        prime_label: "2^63-25",
    },
    FieldSpec {
        name: "fp64bf",
        backend: "Fp64, branch-free sub/redc (hwcost::kernel::Fp64Bf)",
        limbs: 1,
        prime_hex: "7fffffffffffffe7",
        prime_label: "2^63-25",
    },
    FieldSpec {
        name: "fp64cm",
        backend: "Fp64, sub/redc by select_unpredictable (hwcost::kernel::Fp64Cm)",
        limbs: 1,
        prime_hex: "7fffffffffffffe7",
        prime_label: "2^63-25",
    },
    FieldSpec {
        name: "pf64",
        backend: "prime_fast::Fp64",
        limbs: 1,
        prime_hex: "3fffffffffffffc7",
        prime_label: "2^62-57",
    },
    FieldSpec {
        name: "w1",
        backend: "FpWide<1>",
        limbs: 1,
        prime_hex: "ffffffffffffffc5",
        prime_label: "2^64-59",
    },
    FieldSpec {
        name: "w2",
        backend: "FpWide<2>",
        limbs: 2,
        prime_hex: "ffffffffffffffffffffffffffffff61",
        prime_label: "2^128-159",
    },
    FieldSpec {
        name: "w3",
        backend: "FpWide<3>",
        limbs: 3,
        prime_hex: "fffffffffffffffffffffffffffffffeffffffffffffffff",
        prime_label: "P-192",
    },
    FieldSpec {
        name: "w4",
        backend: "FpWide<4>",
        limbs: 4,
        prime_hex: "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
        prime_label: "P-256",
    },
];

pub fn field_spec(name: &str) -> Result<FieldSpec, String> {
    FIELDS
        .iter()
        .copied()
        .find(|f| f.name == name)
        .ok_or_else(|| {
            format!("unknown field `{name}` (have fp64, fp64bf, fp64cm, pf64, w1, w2, w3, w4)")
        })
}

impl FieldSpec {
    pub fn prime(&self) -> BigUint {
        BigUint::parse_bytes(self.prime_hex.as_bytes(), 16).expect("prime hex")
    }
    pub fn bits(&self) -> u64 {
        self.prime().bits()
    }
    pub fn elem_bytes(&self) -> usize {
        8 * self.limbs
    }
}

/// Run `$body` with `$f` bound to the field named by `$spec` (its
/// `PrimeField`, for setup and counting) and, in the three-name form, `$a`
/// to the owned arithmetic that is timed ([`kernel::Owned`] of the field,
/// or [`kernel::Fp64Bf`] for `fp64bf`).
#[macro_export]
macro_rules! hwcost_with_field {
    ($spec:expr, $f:ident => $body:expr) => {
        $crate::hwcost_with_field!($spec, $f, _unused_arith => $body)
    };
    ($spec:expr, $f:ident, $a:ident => $body:expr) => {{
        use $crate::cryptanalysis::fp_backend::{Fp64, FpWide};
        use $crate::cryptanalysis::hwcost::kernel::{Fp64Bf, Owned};
        let spec: $crate::cryptanalysis::hwcost::measure::FieldSpec = $spec;
        let p = spec.prime();
        match spec.name {
            "fp64" => {
                let $f = Fp64::new(num_traits::ToPrimitive::to_u64(&p).unwrap());
                let $a = Owned($f.clone());
                $body
            }
            "fp64bf" => {
                let $f = Fp64::new(num_traits::ToPrimitive::to_u64(&p).unwrap());
                let $a = Fp64Bf::new(&$f);
                $body
            }
            "fp64cm" => {
                let $f = Fp64::new(num_traits::ToPrimitive::to_u64(&p).unwrap());
                let $a = $crate::cryptanalysis::hwcost::kernel::Fp64Cm::new(&$f);
                $body
            }
            "pf64" => {
                let p = num_traits::ToPrimitive::to_u64(&p).unwrap();
                let $f = Fp64::new(p);
                let $a = $crate::cryptanalysis::hwcost::kernel::PrimeFast(
                    $crate::cryptanalysis::prime_fast::Fp64::new(p),
                );
                $body
            }
            "w1" => {
                let $f = FpWide::<1>::new(&p);
                let $a = Owned($f.clone());
                $body
            }
            "w2" => {
                let $f = FpWide::<2>::new(&p);
                let $a = Owned($f.clone());
                $body
            }
            "w3" => {
                let $f = FpWide::<3>::new(&p);
                let $a = Owned($f.clone());
                $body
            }
            "w4" => {
                let $f = FpWide::<4>::new(&p);
                let $a = Owned($f.clone());
                $body
            }
            other => panic!("unknown field {other}"),
        }
    }};
}

// ── Records ─────────────────────────────────────────────────────────

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Sample {
    /// `clock`, `op`, `mem`, `kernel`, `scaling`.
    pub group: String,
    /// The arm, e.g. `w4.mul.throughput`, `mem.1t.65536`, `kernel.w4.1024`.
    pub arm: String,
    pub rep: usize,
    pub value: f64,
    /// `ghz`, `ns_per_op`, `gbps`, `ns_per_iter`.
    pub unit: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CacheInfo {
    pub level: u32,
    pub kind: String,
    pub size_bytes: u64,
    pub shared_cpu_list: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct KernelCounts {
    pub field: String,
    pub batch: usize,
    pub mul: f64,
    /// Multiplications on the batch inversion's serial chain.
    pub mul_dep: f64,
    pub sqr: f64,
    pub add: f64,
    pub inv: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Params {
    pub reps: usize,
    pub target_ms: f64,
    pub batches: Vec<usize>,
    pub mem_sizes: Vec<usize>,
    pub threads: Vec<usize>,
    pub scaling_batch: usize,
    pub fields: Vec<String>,
    pub seed: u64,
    /// When set, also time `prime_fast::rho_parallel` end to end on a
    /// prime-order `a = −3` curve of this many bits (arm `pfrho.<bits>`).
    #[serde(default)]
    pub prime_fast_rho_bits: Option<u32>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Measurements {
    pub schema: String,
    pub params: Params,
    /// `ecbench`'s host capsule (CPU, features, topology, build facts).
    pub host: serde_json::Value,
    pub caches: Vec<CacheInfo>,
    /// Target features the measured binary was compiled with.
    pub compiled_features: Vec<String>,
    pub counts: Vec<KernelCounts>,
    pub samples: Vec<Sample>,
}

impl Measurements {
    /// Median of an arm's samples.
    pub fn median(&self, arm: &str) -> Option<f64> {
        let mut v: Vec<f64> = self
            .samples
            .iter()
            .filter(|s| s.arm == arm)
            .map(|s| s.value)
            .collect();
        median(&mut v)
    }
    /// First and third quartile of an arm's samples.
    pub fn quartiles(&self, arm: &str) -> Option<(f64, f64)> {
        let mut v: Vec<f64> = self
            .samples
            .iter()
            .filter(|s| s.arm == arm)
            .map(|s| s.value)
            .collect();
        if v.is_empty() {
            return None;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let q = |t: f64| {
            let x = t * (v.len() - 1) as f64;
            let (i, fr) = (x.floor() as usize, x.fract());
            if i + 1 < v.len() {
                v[i] + fr * (v[i + 1] - v[i])
            } else {
                v[i]
            }
        };
        Some((q(0.25), q(0.75)))
    }
}

pub fn median(v: &mut [f64]) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    Some(if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    })
}

pub fn compiled_features() -> Vec<String> {
    let mut v = Vec::new();
    macro_rules! feat {
        ($($f:literal),*) => {$(
            if cfg!(target_feature = $f) {
                v.push($f.to_string());
            }
        )*};
    }
    feat!(
        "sse4.2",
        "avx",
        "avx2",
        "fma",
        "bmi1",
        "bmi2",
        "adx",
        "avx512f",
        "avx512ifma",
        "pclmulqdq",
        "vpclmulqdq",
        "neon"
    );
    v
}

pub fn caches() -> Vec<CacheInfo> {
    let mut out = Vec::new();
    for i in 0..8 {
        let dir = format!("/sys/devices/system/cpu/cpu0/cache/index{i}");
        let read = |f: &str| {
            std::fs::read_to_string(format!("{dir}/{f}"))
                .ok()
                .map(|s| s.trim().to_string())
        };
        let (Some(level), Some(kind), Some(size)) = (read("level"), read("type"), read("size"))
        else {
            continue;
        };
        let bytes = parse_size(&size).unwrap_or(0);
        out.push(CacheInfo {
            level: level.parse().unwrap_or(0),
            kind,
            size_bytes: bytes,
            shared_cpu_list: read("shared_cpu_list").unwrap_or_default(),
        });
    }
    out
}

fn parse_size(s: &str) -> Option<u64> {
    let s = s.trim();
    let (num, mult) = if let Some(n) = s.strip_suffix('K') {
        (n, 1024)
    } else if let Some(n) = s.strip_suffix('M') {
        (n, 1024 * 1024)
    } else if let Some(n) = s.strip_suffix('G') {
        (n, 1024 * 1024 * 1024)
    } else {
        (s, 1)
    };
    num.parse::<u64>().ok().map(|n| n * mult)
}

// ── Arms ────────────────────────────────────────────────────────────

/// Runs `n` units of work and returns (elapsed ns, units done).
type Run = Box<dyn FnMut(u64) -> (f64, f64)>;

struct Arm {
    group: &'static str,
    name: String,
    unit: &'static str,
    run: Run,
    /// Converts (elapsed ns, units) into the recorded value.
    value: fn(f64, f64) -> f64,
    n: u64,
}

fn ns_per_unit(ns: f64, units: f64) -> f64 {
    ns / units
}

fn calibrate(arm: &mut Arm, target_ns: f64) {
    let mut n = 1u64;
    loop {
        let (ns, _) = (arm.run)(n);
        if ns >= target_ns / 4.0 || n >= 1 << 40 {
            arm.n = ((n as f64) * target_ns / ns.max(1.0)).ceil().max(1.0) as u64;
            return;
        }
        n *= 2;
    }
}

fn run_group(
    arms: &mut [Arm],
    reps: usize,
    target_ms: f64,
    samples: &mut Vec<Sample>,
    log: &mut dyn FnMut(String),
) {
    let target_ns = target_ms * 1e6;
    for a in arms.iter_mut() {
        calibrate(a, target_ns);
    }
    for rep in 0..reps {
        let k = arms.len();
        for i in 0..k {
            let a = &mut arms[(i + rep) % k];
            let (ns, units) = (a.run)(a.n);
            samples.push(Sample {
                group: a.group.into(),
                arm: a.name.clone(),
                rep,
                value: (a.value)(ns, units),
                unit: a.unit.into(),
            });
        }
        log(format!("  rep {}/{reps} of {k} arms", rep + 1));
    }
}

// Clock: dependent chains of register-register `add` (one cycle each) and
// `imul` (three cycles each on every x86-64 core of the last decade, as a
// cross-check).  A chain of `add reg, imm` does not measure the clock:
// Golden Cove-class cores fold dependent immediate adds at rename and run
// the chain several times faster than one per cycle.

#[cfg(target_arch = "x86_64")]
macro_rules! chain8 {
    ($ins:literal) => {
        concat!(
            $ins, "\n", $ins, "\n", $ins, "\n", $ins, "\n", $ins, "\n", $ins, "\n", $ins, "\n",
            $ins
        )
    };
}

#[cfg(target_arch = "x86_64")]
fn add_chain(n: u64) -> u64 {
    let mut x: u64 = 0;
    let y: u64 = black_box(1);
    for _ in 0..n {
        // SAFETY: register-only arithmetic.
        unsafe {
            std::arch::asm!(
                chain8!("add {x}, {y}"), chain8!("add {x}, {y}"),
                chain8!("add {x}, {y}"), chain8!("add {x}, {y}"),
                x = inout(reg) x, y = in(reg) y,
                options(nomem, nostack)
            );
        }
    }
    x
}

#[cfg(target_arch = "x86_64")]
fn imul_chain(n: u64) -> u64 {
    let mut x: u64 = 3;
    let y: u64 = black_box(5);
    for _ in 0..n {
        // SAFETY: register-only arithmetic.
        unsafe {
            std::arch::asm!(
                chain8!("imul {x}, {y}"), chain8!("imul {x}, {y}"),
                chain8!("imul {x}, {y}"), chain8!("imul {x}, {y}"),
                x = inout(reg) x, y = in(reg) y,
                options(nomem, nostack)
            );
        }
    }
    x
}

#[cfg(target_arch = "x86_64")]
fn clock_arms() -> Vec<Arm> {
    vec![
        Arm {
            group: "clock",
            name: "clock.add_chain".into(),
            unit: "ghz",
            run: Box::new(|n| {
                let t = Instant::now();
                black_box(add_chain(black_box(n)));
                (t.elapsed().as_nanos() as f64, (32 * n) as f64)
            }),
            value: |ns, cycles| cycles / ns,
            n: 0,
        },
        Arm {
            group: "clock",
            name: "clock.imul_chain".into(),
            unit: "ghz",
            run: Box::new(|n| {
                let t = Instant::now();
                black_box(imul_chain(black_box(n)));
                (t.elapsed().as_nanos() as f64, (3 * 32 * n) as f64)
            }),
            value: |ns, cycles| cycles / ns,
            n: 0,
        },
    ]
}

#[cfg(not(target_arch = "x86_64"))]
fn clock_arms() -> Vec<Arm> {
    Vec::new()
}

fn random_elems<F: PrimeField>(f: &F, k: usize, rng: &mut StdRng) -> Vec<F::Elem> {
    (0..k)
        .map(|_| f.from_biguint(&rng.gen_biguint_below(f.modulus())))
        .collect()
}

/// Microbenchmarks of `arith`'s operations on random elements of `f`.
fn op_arms<F, A>(f: &F, arith: A, name: &str, seed: u64) -> Vec<Arm>
where
    F: PrimeField,
    F::Elem: 'static,
    A: Arith<E = F::Elem> + Clone + 'static,
{
    let mut rng = StdRng::seed_from_u64(seed);
    let a0 = random_elems(f, 8, &mut rng);
    let b0 = random_elems(f, 8, &mut rng);
    let mut arms = Vec::new();
    let mk = |op: &str, mode: &str, run: Run| Arm {
        group: "op",
        name: format!("{name}.{op}.{mode}"),
        unit: "ns_per_op",
        run,
        value: ns_per_unit,
        n: 0,
    };
    macro_rules! chains {
        ($k:expr, |$f:ident, $a:ident, $b:ident| $e:expr) => {{
            let ff = arith.clone();
            let (mut a, b) = (a0.clone(), b0.clone());
            Box::new(move |n: u64| {
                let $f = &ff;
                let t = Instant::now();
                for _ in 0..n {
                    for k in 0..$k {
                        let ($a, $b) = (a[k], b[k]);
                        a[k] = $e;
                    }
                }
                black_box(&a);
                (t.elapsed().as_nanos() as f64, (n * $k) as f64)
            }) as Run
        }};
    }
    arms.push(mk("mul", "throughput", chains!(8, |f, x, y| f.mul(x, y))));
    arms.push(mk("mul", "latency", chains!(1, |f, x, y| f.mul(x, y))));
    arms.push(mk("sqr", "throughput", chains!(8, |f, x, _y| f.sqr(x))));
    arms.push(mk("sqr", "latency", chains!(1, |f, x, _y| f.sqr(x))));
    arms.push(mk("add", "throughput", chains!(8, |f, x, y| f.sub(x, y))));
    arms.push(mk("add", "latency", chains!(1, |f, x, y| f.sub(x, y))));
    // The `+ 1` keeps the chain from alternating between two values.
    arms.push(mk(
        "inv",
        "throughput",
        chains!(4, |f, x, _y| f.inv(f.add(x, f.one())).unwrap_or(x)),
    ));
    arms.push(mk(
        "inv",
        "latency",
        chains!(1, |f, x, _y| f.inv(f.add(x, f.one())).unwrap_or(x)),
    ));
    arms
}

/// The production prime-field rho, `prime_fast::rho_parallel`, end to end:
/// one known-answer solve per call on the prime-order `a = −3` curve of
/// `bits` bits that [`prime_fast::find_a3_curve`] returns, one thread with
/// 1024 lock-step walkers, the negation map, a 2¹⁰-entry r-adding table and
/// 6 distinguished bits; wall time (table set-up included) per walk step.
///
/// [`prime_fast::find_a3_curve`]: crate::cryptanalysis::prime_fast::find_a3_curve
fn prime_fast_rho_arm(bits: u32, seed: u64) -> Result<Arm, String> {
    use crate::cryptanalysis::prime_fast::{find_a3_curve, rho_parallel, RhoConfig};
    use rand::Rng;
    let (curve, g) =
        find_a3_curve(bits, 3).ok_or_else(|| format!("no {bits}-bit a = -3 prime-order curve"))?;
    let cfg = RhoConfig {
        threads: 1,
        walkers: 1024,
        table_bits: 10,
        dp_bits: 6,
        negation: true,
        seed,
        max_steps: 0,
    };
    let mut rng = StdRng::seed_from_u64(seed ^ u64::from(bits));
    Ok(Arm {
        group: "pfrho",
        name: format!("pfrho.{bits}"),
        unit: "ns_per_step",
        run: Box::new(move |n: u64| {
            let (mut ns, mut steps) = (0.0, 0.0);
            for _ in 0..n {
                let x = rng.gen_range(1..curve.n);
                let q = curve.mul(Some(g), x).expect("x < n gives a finite point");
                let cfg = RhoConfig {
                    seed: rng.gen(),
                    ..cfg.clone()
                };
                let t = Instant::now();
                let r = rho_parallel(&curve, g, q, &cfg).expect("unlimited rho solves");
                ns += t.elapsed().as_nanos() as f64;
                assert_eq!(r.log, x, "pfrho: wrong logarithm");
                steps += r.total_steps as f64;
            }
            (ns, steps)
        }),
        value: ns_per_unit,
        n: 0,
    })
}

fn read_sum(buf: &[u64]) -> u64 {
    buf.iter().fold(0u64, |a, &x| a.wrapping_add(x))
}

fn mem_arm(bytes: usize, threads: usize) -> Arm {
    let words = (bytes / 8).max(1);
    let run: Run = if threads == 1 {
        let buf: Vec<u64> = (0..words as u64).collect();
        Box::new(move |n| {
            let t = Instant::now();
            let mut s = 0u64;
            for _ in 0..n {
                s = s.wrapping_add(read_sum(black_box(&buf)));
            }
            black_box(s);
            (
                t.elapsed().as_nanos() as f64,
                (n * (words as u64) * 8) as f64,
            )
        })
    } else {
        Box::new(move |n| {
            let barrier = Barrier::new(threads);
            let times: Vec<f64> = std::thread::scope(|sc| {
                let hs: Vec<_> = (0..threads)
                    .map(|_| {
                        let barrier = &barrier;
                        sc.spawn(move || {
                            let buf: Vec<u64> = (0..words as u64).collect();
                            black_box(read_sum(&buf));
                            barrier.wait();
                            let t = Instant::now();
                            let mut s = 0u64;
                            for _ in 0..n {
                                s = s.wrapping_add(read_sum(black_box(&buf)));
                            }
                            black_box(s);
                            t.elapsed().as_nanos() as f64
                        })
                    })
                    .collect();
                hs.into_iter().map(|h| h.join().unwrap()).collect()
            });
            let wall = times.iter().cloned().fold(0.0, f64::max);
            (wall, (n * (words as u64) * 8 * threads as u64) as f64)
        })
    };
    Arm {
        group: "mem",
        name: format!("mem.{threads}t.{bytes}"),
        unit: "gbps",
        run,
        value: |ns, bytes| bytes / ns,
        n: 0,
    }
}

/// `b` walks for a timed run: random points for `b ≤ 4096`; above that,
/// sums `Q_{i mod 4096} + S_{⌊i/4096⌋}` of two random pools, computed with
/// one batch inversion (distinct with overwhelming probability), so a
/// million walks cost a million additions rather than a million square
/// roots.
pub fn walks_for_timing<F: PrimeField>(
    f: &F,
    b: usize,
    seed: u64,
) -> (Walks<F::Elem>, kernel::Table<F::Elem>)
where
    F::Elem: LowWord,
{
    const POOL: usize = 4096;
    if b <= POOL {
        return kernel::setup(f, b, seed);
    }
    let (pool, t) = kernel::setup(f, POOL, seed);
    let mut rng = StdRng::seed_from_u64(seed ^ 0x9e37_79b9_7f4a_7c15);
    let shifts: Vec<_> = (0..b.div_ceil(POOL))
        .map(|_| kernel::random_point(f, &mut rng))
        .collect();
    let mut dx: Vec<F::Elem> = (0..b)
        .map(|i| f.sub(shifts[i / POOL].0, pool.xs[i % POOL]))
        .collect();
    assert!(
        dx.iter().all(|d| !f.is_zero(*d)),
        "x collision building walks"
    );
    let mut scratch = Vec::new();
    f.batch_inv(&mut dx, &mut scratch);
    let mut xs = Vec::with_capacity(b);
    let mut ys = Vec::with_capacity(b);
    for i in 0..b {
        let (x1, y1) = (pool.xs[i % POOL], pool.ys[i % POOL]);
        let (x2, y2) = shifts[i / POOL];
        let l = f.mul(f.sub(y2, y1), dx[i]);
        let x3 = f.sub(f.sub(f.sqr(l), x1), x2);
        ys.push(f.sub(f.mul(l, f.sub(x1, x3)), y1));
        xs.push(x3);
    }
    (Walks::new(f, xs, ys), t)
}

/// The rho kernel over `f`'s walks, timed through `arith` (which must
/// compute exactly what `f` computes; `run` checks it).
fn kernel_arm<F, A>(f: F, arith: A, name: &str, b: usize, seed: u64) -> (Arm, KernelCounts)
where
    F: PrimeField + 'static,
    F::Elem: LowWord + 'static,
    A: Arith<E = F::Elem> + 'static,
{
    let (mut w, t) = walks_for_timing(&f, b, seed);
    // One counted step on the very walks that are timed.
    let tally = {
        let c = Counted::new(&f);
        kernel::step(&c, &mut w, &t, DP_MASK);
        c.tally.get()
    };
    assert_eq!(w.degenerate, 0, "{name} B={b}: degenerate step");
    let per = |x: u64| x as f64 / b as f64;
    let counts = KernelCounts {
        field: name.into(),
        batch: b,
        mul: per(tally.mul),
        mul_dep: per(tally.mul_dep),
        sqr: per(tally.sqr),
        add: per(tally.add),
        inv: per(tally.inv),
    };
    let arm_name = format!("kernel.{name}.{b}");
    let run: Run = Box::new(move |n| {
        let t0 = Instant::now();
        for _ in 0..n {
            kernel::step(&arith, &mut w, &t, DP_MASK);
        }
        let ns = t0.elapsed().as_nanos() as f64;
        assert_eq!(w.degenerate, 0, "degenerate step while timing");
        black_box(&w.xs);
        (ns, (n as f64) * b as f64)
    });
    (
        Arm {
            group: "kernel",
            name: arm_name,
            unit: "ns_per_iter",
            run,
            value: ns_per_unit,
            n: 0,
        },
        counts,
    )
}

fn scaling_arm<F, A>(f: F, arith: A, name: &str, b: usize, threads: usize, seed: u64) -> Arm
where
    F: PrimeField + 'static,
    F::Elem: LowWord + 'static,
    A: Arith<E = F::Elem> + Sync + 'static,
{
    // Each thread owns its walks, built once, outside the timing.
    let mut sets: Vec<_> = (0..threads)
        .map(|i| walks_for_timing(&f, b, seed + i as u64))
        .collect();
    let run: Run = Box::new(move |n| {
        let barrier = Barrier::new(threads);
        let ar = &arith;
        let times: Vec<f64> = std::thread::scope(|sc| {
            let hs: Vec<_> = sets
                .iter_mut()
                .map(|(w, t)| {
                    let barrier = &barrier;
                    sc.spawn(move || {
                        barrier.wait();
                        let t0 = Instant::now();
                        for _ in 0..n {
                            kernel::step(ar, w, t, DP_MASK);
                        }
                        let ns = t0.elapsed().as_nanos() as f64;
                        assert_eq!(w.degenerate, 0);
                        ns
                    })
                })
                .collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        });
        let wall = times.iter().cloned().fold(0.0, f64::max);
        (wall, (n as f64) * (b * threads) as f64)
    });
    Arm {
        group: "scaling",
        name: format!("scaling.{name}.{threads}t"),
        unit: "ns_per_iter",
        run,
        value: ns_per_unit,
        n: 0,
    }
}

// ── The run ─────────────────────────────────────────────────────────

pub fn default_params(quick: bool) -> Params {
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let mut threads = vec![1];
    let mut t = 2;
    while t < cores {
        threads.push(t);
        t *= 2;
    }
    if cores > 1 {
        threads.push(cores);
    }
    if quick {
        Params {
            reps: 2,
            target_ms: 5.0,
            batches: vec![16, 256, 4096],
            mem_sizes: vec![16 << 10, 1 << 20, 64 << 20],
            threads,
            scaling_batch: 256,
            fields: vec!["fp64".into(), "w2".into()],
            seed: 20261016,
            prime_fast_rho_bits: None,
        }
    } else {
        Params {
            reps: 7,
            target_ms: 100.0,
            batches: vec![
                16,
                64,
                256,
                1024,
                4096,
                16384,
                65536,
                262144,
                1 << 20,
                1 << 22,
            ],
            mem_sizes: vec![
                16 << 10,
                32 << 10,
                128 << 10,
                512 << 10,
                1 << 20,
                4 << 20,
                16 << 20,
                64 << 20,
                256 << 20,
                1 << 30,
            ],
            threads,
            scaling_batch: 1024,
            fields: FIELDS.iter().map(|f| f.name.to_string()).collect(),
            seed: 20261016,
            prime_fast_rho_bits: None,
        }
    }
}

pub fn run(params: Params, log: &mut dyn FnMut(String)) -> Result<Measurements, String> {
    for name in &params.fields {
        let spec = field_spec(name)?;
        if !is_probable_prime(&spec.prime()) {
            return Err(format!("{name}: {} is not prime", spec.prime_label));
        }
        let r = crate::hwcost_with_field!(spec, f, a => {
            kernel::check(&f, 33, 4, params.seed)
                .and_then(|_| kernel::check_same_walk(&f, &a, 33, 4, params.seed))
        });
        r?;
        if spec.name == "fp64bf" {
            let p = num_traits::ToPrimitive::to_u64(&spec.prime()).unwrap();
            kernel::check_wrappers(&crate::cryptanalysis::fp_backend::Fp64::new(p), params.seed)?;
        }
        log(format!(
            "{name} ({}, {}): batched walk == unbatched reference, on curve; timed arithmetic == field",
            spec.backend, spec.prime_label
        ));
    }
    let host = crate::cryptanalysis::ecbench::host::capture()
        .map(|c| serde_json::to_value(c).unwrap_or(serde_json::Value::Null))
        .unwrap_or_else(|e| serde_json::json!({ "error": e }));
    let mut samples = Vec::new();
    let (reps, target) = (params.reps, params.target_ms);

    log("clock".into());
    let mut arms = clock_arms();
    run_group(&mut arms, reps, target, &mut samples, log);

    log("field operations".into());
    let mut arms = Vec::new();
    for name in &params.fields {
        let spec = field_spec(name)?;
        let mut v = crate::hwcost_with_field!(spec, f, a => op_arms(&f, a, name, params.seed));
        arms.append(&mut v);
    }
    run_group(&mut arms, reps, target, &mut samples, log);
    drop(arms);

    for &t in &params.threads {
        if t != 1 && t != *params.threads.last().unwrap() {
            continue;
        }
        log(format!("memory read bandwidth, {t} thread(s)"));
        let mut arms: Vec<Arm> = params.mem_sizes.iter().map(|&s| mem_arm(s, t)).collect();
        run_group(&mut arms, reps, target, &mut samples, log);
    }

    log("rho kernel".into());
    let mut counts = Vec::new();
    let mut arms = Vec::new();
    for name in &params.fields {
        let spec = field_spec(name)?;
        for &b in &params.batches {
            let (arm, c) = crate::hwcost_with_field!(spec, f, a => {
                kernel_arm(f, a, name, b, params.seed ^ b as u64)
            });
            arms.push(arm);
            counts.push(c);
        }
        log(format!(
            "  {name}: {} batch sizes set up",
            params.batches.len()
        ));
    }
    run_group(&mut arms, reps, target, &mut samples, log);
    drop(arms);

    log("thread scaling".into());
    let mut arms = Vec::new();
    for name in &params.fields {
        let spec = field_spec(name)?;
        for &t in &params.threads {
            arms.push(crate::hwcost_with_field!(spec, f, a => {
                scaling_arm(f, a, name, params.scaling_batch, t, params.seed)
            }));
        }
    }
    run_group(&mut arms, reps, target, &mut samples, log);
    drop(arms);

    if let Some(bits) = params.prime_fast_rho_bits {
        log(format!(
            "prime_fast::rho_parallel, {bits}-bit curve, one thread"
        ));
        let mut arms = vec![prime_fast_rho_arm(bits, params.seed)?];
        run_group(&mut arms, reps, target, &mut samples, log);
    }

    Ok(Measurements {
        schema: MEASUREMENTS_SCHEMA.into(),
        params,
        host,
        caches: caches(),
        compiled_features: compiled_features(),
        counts,
        samples,
    })
}

/// Exact per-iteration counts of the kernel at batch `b`, by running it
/// through the counting field (for workloads on fields not timed here).
pub fn counts_for(name: &str, b: usize, seed: u64) -> Result<KernelCounts, String> {
    let spec = field_spec(name)?;
    let (t, iters): (OpTally, f64) = crate::hwcost_with_field!(spec, f => {
        let (w, tbl) = walks_for_timing(&f, b, seed);
        let mut w = w;
        let c = Counted::new(&f);
        kernel::step(&c, &mut w, &tbl, DP_MASK);
        assert_eq!(w.degenerate, 0);
        (c.tally.get(), b as f64)
    });
    Ok(KernelCounts {
        field: name.into(),
        batch: b,
        mul: t.mul as f64 / iters,
        mul_dep: t.mul_dep as f64 / iters,
        sqr: t.sqr as f64 / iters,
        add: t.add as f64 / iters,
        inv: t.inv as f64 / iters,
    })
}

/// Medians of every arm, for tables.
pub fn medians(m: &Measurements) -> BTreeMap<String, f64> {
    let mut arms: Vec<&str> = m.samples.iter().map(|s| s.arm.as_str()).collect();
    arms.sort();
    arms.dedup();
    arms.into_iter()
        .filter_map(|a| m.median(a).map(|v| (a.to_string(), v)))
        .collect()
}
