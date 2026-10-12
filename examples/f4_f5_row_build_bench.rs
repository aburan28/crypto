//! Native allocation audit of complete Boolean matrix-F4/F5 calls.
//! No wall-clock measurements: ordinary hosts do not establish isolation.
//! Run with `cargo run --release --example f4_f5_row_build_bench -- [repeats]`.
//! `profile [repeats]` adds frozen holdouts, peak-live-byte accounting and
//! per-call Callgrind dumps; it makes no wall-time measurement.

use crypto_lib::cryptanalysis::koblitz_groebner::matrix_f4_f2_counted;
use crypto_lib::cryptanalysis::matrix_f5_f2::{matrix_f5_f2, F5Criterion};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use serde_json::json;
use std::alloc::{GlobalAlloc, Layout, System};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

struct CountingAllocator;

static ENABLED: AtomicBool = AtomicBool::new(false);
static ALLOCATIONS: AtomicU64 = AtomicU64::new(0);
static BYTES: AtomicU64 = AtomicU64::new(0);
static LIVE: AtomicU64 = AtomicU64::new(0);
static PEAK: AtomicU64 = AtomicU64::new(0);

fn live_increase(bytes: usize) {
    let live = LIVE.fetch_add(bytes as u64, Ordering::Relaxed) + bytes as u64;
    if ENABLED.load(Ordering::Relaxed) {
        PEAK.fetch_max(live, Ordering::Relaxed);
    }
}

fn record(bytes: usize) {
    if ENABLED.load(Ordering::Relaxed) {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        BYTES.fetch_add(bytes as u64, Ordering::Relaxed);
    }
}

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        record(layout.size());
        let ptr = unsafe { System.alloc(layout) };
        if !ptr.is_null() {
            live_increase(layout.size());
        }
        ptr
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        record(layout.size());
        let ptr = unsafe { System.alloc_zeroed(layout) };
        if !ptr.is_null() {
            live_increase(layout.size());
        }
        ptr
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        LIVE.fetch_sub(layout.size() as u64, Ordering::Relaxed);
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        record(size);
        let ptr = unsafe { System.realloc(ptr, layout, size) };
        if !ptr.is_null() {
            if size >= layout.size() {
                live_increase(size - layout.size());
            } else {
                LIVE.fetch_sub((layout.size() - size) as u64, Ordering::Relaxed);
            }
        }
        ptr
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

#[derive(Default)]
struct Audit {
    allocations: u64,
    bytes: u64,
    live_at_entry: u64,
    peak_live: u64,
}

// The standard Valgrind amd64 client-request ABI: rotate RDI by a full
// 128 bits, then its recognized RBX marker. Native execution is a no-op.
#[cfg(target_arch = "x86_64")]
fn callgrind_request(request: usize, argument: usize) {
    let args = [request, argument, 0, 0, 0, 0];
    let mut result = 0usize;
    // SAFETY: the array lives through the request. RAX carries its address,
    // RDX the default/result. RDI and condition flags are declared clobbered;
    // the rotations and RBX self-exchange preserve other registers.
    unsafe {
        std::arch::asm!(
            "rol rdi, 3", "rol rdi, 13", "rol rdi, 61", "rol rdi, 51",
            "xchg rbx, rbx",
            in("rax") args.as_ptr(), inout("rdx") result, lateout("rdi") _,
            options(nostack),
        );
    }
    std::hint::black_box(result);
}

#[cfg(not(target_arch = "x86_64"))]
fn callgrind_request(_: usize, _: usize) {
    panic!("the profile client-request ABI is currently supported on x86_64 only");
}

fn audited<T>(label: Option<&std::ffi::CStr>, run: impl FnOnce() -> T) -> (T, Audit) {
    ALLOCATIONS.store(0, Ordering::Relaxed);
    BYTES.store(0, Ordering::Relaxed);
    let live_at_entry = LIVE.load(Ordering::Relaxed);
    PEAK.store(live_at_entry, Ordering::Relaxed);
    if label.is_some() {
        callgrind_request(0x4354_0001, 0); // ZERO_STATS (discard prior stop request)
        callgrind_request(0x4354_0004, 0); // START_INSTRUMENTATION
    }
    ENABLED.store(true, Ordering::Relaxed);
    let result = run();
    ENABLED.store(false, Ordering::Relaxed);
    let audit = Audit {
        allocations: ALLOCATIONS.load(Ordering::Relaxed),
        bytes: BYTES.load(Ordering::Relaxed),
        live_at_entry,
        peak_live: PEAK.load(Ordering::Relaxed),
    };
    if let Some(label) = label {
        callgrind_request(0x4354_0003, label.as_ptr() as usize); // DUMP_STATS_AT
        callgrind_request(0x4354_0005, 0); // STOP_INSTRUMENTATION
    }
    (result, audit)
}

fn profile_label(
    profile: bool,
    case: &str,
    seed: u64,
    repeat: usize,
    arm: &str,
) -> Option<std::ffi::CString> {
    profile.then(|| std::ffi::CString::new(json!([case, seed, repeat, arm]).to_string()).unwrap())
}

fn emit(mut row: serde_json::Value, audit: Audit, profile: bool) {
    row["allocations"] = json!(audit.allocations);
    row["allocated_bytes"] = json!(audit.bytes);
    if profile {
        row["live_at_entry_bytes"] = json!(audit.live_at_entry);
        row["peak_live_bytes"] = json!(audit.peak_live);
    }
    println!("{row}");
}

fn criterion_digest(criterion: &F5Criterion, m: usize, n_vars: usize, degree: u32) -> String {
    let mut multipliers = vec![0u64];
    let mut level = vec![0u64];
    for _ in 0..degree {
        let mut next = Vec::new();
        for &mask in &level {
            for bit in (64 - mask.leading_zeros())..n_vars as u32 {
                next.push(mask | (1u64 << bit));
            }
        }
        multipliers.extend_from_slice(&next);
        level = next;
    }
    let mut hash = blake3::Hasher::new();
    for i in 0..m {
        for &multiplier in &multipliers {
            hash.update(&[u8::from(criterion.prunes(i, multiplier))]);
        }
    }
    hash.finalize().to_hex().to_string()
}

fn digest(polys: &[F2BoolPoly]) -> String {
    let mut hash = blake3::Hasher::new();
    for poly in polys {
        hash.update(&(poly.n_vars as u64).to_le_bytes());
        hash.update(&(poly.terms.len() as u64).to_le_bytes());
        for term in &poly.terms {
            hash.update(&term.mask.to_le_bytes());
        }
    }
    hash.finalize().to_hex().to_string()
}

fn system(vars: &[usize], m: usize, seed: u64) -> Vec<F2BoolPoly> {
    let mut state = seed;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let n_vars = vars.iter().max().unwrap() + 1;
    let planted = next();
    let mut monos = vec![0];
    for (i, &v) in vars.iter().enumerate() {
        monos.push(1u64 << v);
        for &w in &vars[i + 1..] {
            monos.push((1u64 << v) | (1u64 << w));
        }
    }
    (0..m)
        .map(|_| {
            let mut terms: Vec<F2BoolMono> = monos
                .iter()
                .filter(|_| next() & 1 != 0)
                .map(|&mask| F2BoolMono::from_mask(mask))
                .collect();
            if terms.iter().filter(|t| t.mask & planted == t.mask).count() % 2 == 1 {
                terms.push(F2BoolMono::from_mask(0));
            }
            let poly = F2BoolPoly::from_monos(terms, n_vars);
            assert!(!poly.is_zero());
            poly
        })
        .collect()
}

fn compare(baseline: &str, candidate: &str) {
    let read = |path: &str| {
        let mut rows = BTreeMap::new();
        for line in std::fs::read_to_string(path).unwrap().lines() {
            let row: serde_json::Value = serde_json::from_str(line).unwrap();
            let key = json!([row["case"], row["seed"], row["repeat"], row["arm"]]).to_string();
            assert!(rows.insert(key, row).is_none(), "duplicate audit row");
        }
        rows
    };
    let baseline = read(baseline);
    let candidate = read(candidate);
    assert_eq!(
        baseline.len(),
        72,
        "all eight systems and repetitions required"
    );
    assert_eq!(
        baseline.keys().collect::<Vec<_>>(),
        candidate.keys().collect::<Vec<_>>()
    );
    let mut totals: BTreeMap<String, [u64; 4]> = BTreeMap::new();
    let mut comparisons = Vec::new();
    for (key, a) in baseline {
        let b = &candidate[&key];
        let mut expected = a.clone();
        let mut actual = b.clone();
        for field in ["allocations", "allocated_bytes"] {
            expected.as_object_mut().unwrap().remove(field);
            actual.as_object_mut().unwrap().remove(field);
        }
        assert_eq!(expected, actual, "output or accounting mismatch: {key}");
        let counts = [
            a["allocations"].as_u64().unwrap(),
            b["allocations"].as_u64().unwrap(),
            a["allocated_bytes"].as_u64().unwrap(),
            b["allocated_bytes"].as_u64().unwrap(),
        ];
        let total = totals
            .entry(a["arm"].as_str().unwrap().to_owned())
            .or_default();
        for (sum, count) in total.iter_mut().zip(counts) {
            *sum += count;
        }
        comparisons.push(json!({"key": key, "counts": counts}));
    }
    let acceptance = totals["f4"][1] < totals["f4"][0]
        && totals["criterion"][1] * 5 <= totals["criterion"][0] * 4;
    let arms: Vec<_> = totals
        .into_iter()
        .map(|(arm, c)| {
            json!({
                "arm": arm, "calls": 24,
                "baseline_allocations": c[0], "candidate_allocations": c[1],
                "allocation_reduction_percent": 100.0 * (1.0 - c[1] as f64 / c[0] as f64),
                "baseline_allocated_bytes": c[2], "candidate_allocated_bytes": c[3],
                "allocated_byte_reduction_percent": 100.0 * (1.0 - c[3] as f64 / c[2] as f64),
            })
        })
        .collect();
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "scope": "Boolean matrix-F4/F5 allocation diagnostic",
            "systems": 8, "repeats": 3, "exact_outputs_and_counts": true,
            "acceptance_met": acceptance, "wall_time_speedup": null,
            "arms": arms, "comparisons": comparisons,
        }))
        .unwrap()
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).is_some_and(|s| s == "compare") {
        assert_eq!(args.len(), 4, "compare BASELINE.jsonl CANDIDATE.jsonl");
        compare(&args[2], &args[3]);
        return;
    }
    let profile = args.get(1).is_some_and(|s| s == "profile");
    let repeats = args
        .get(if profile { 2 } else { 1 })
        .map(|s| s.parse::<usize>().expect("positive repeat count"))
        .unwrap_or(3);
    assert!(repeats > 0);
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .unwrap();
    let mut workloads = vec![
        ("n8_d3", (0..8).collect::<Vec<_>>(), 8, 3),
        ("n8_d4", (0..8).collect::<Vec<_>>(), 8, 4),
        ("n10_d4", (0..10).collect::<Vec<_>>(), 10, 4),
        ("n10_d5", (0..10).collect::<Vec<_>>(), 10, 5),
    ];
    let mut cases = Vec::new();
    for (case, vars, m, degree) in workloads.drain(..) {
        for seed in [0x2545_f491_4f6c_dd1d, 0x9e37_79b9_7f4a_7c15] {
            cases.push((case, vars.clone(), m, degree, vars.len(), seed));
        }
    }
    if profile {
        for (case, vars, m, degree) in [
            ("holdout_n12_d3", (0..12).collect::<Vec<_>>(), 12, 3),
            ("holdout_n12_d4", (0..12).collect::<Vec<_>>(), 12, 4),
            ("holdout_shifted_d4", vec![52, 54, 57, 60, 62, 63], 6, 4),
        ] {
            cases.push((
                case,
                vars,
                m,
                degree,
                if case == "holdout_shifted_d4" { 6 } else { 12 },
                0xd1b5_4a32_d192_ed03,
            ));
        }
    }
    for (case, vars, m, degree, n_vars, seed) in cases {
        let polys = system(&vars, m, seed);
        let mask = (1u64 << n_vars) - 1;
        let input_digest = digest(&polys);
        // Warm lazy configuration and worker startup outside the audit.
        let expected_f4 = matrix_f4_f2_counted(&polys, n_vars, degree).unwrap();
        let expected_f5 = matrix_f5_f2(&polys, n_vars, degree).unwrap();
        assert_eq!(expected_f4.0, expected_f5.0);
        let expected_criterion = F5Criterion::new(&polys, n_vars, degree, mask);
        for repeat in 0..repeats {
            let label = profile_label(profile, case, seed, repeat, "f4");
            let (f4, audit) = audited(label.as_deref(), || {
                matrix_f4_f2_counted(&polys, n_vars, degree).unwrap()
            });
            assert_eq!(f4, expected_f4);
            emit(
                json!({
                    "case": case, "seed": seed, "repeat": repeat, "arm": "f4",
                    "input_digest": input_digest, "output_digest": digest(&f4.0),
                    "rank": f4.0.len(), "word_ops": f4.1,
                }),
                audit,
                profile,
            );
            let label = profile_label(profile, case, seed, repeat, "f5");
            let (f5, audit) = audited(label.as_deref(), || {
                matrix_f5_f2(&polys, n_vars, degree).unwrap()
            });
            assert_eq!(f5.0, expected_f5.0);
            assert_eq!(f5.1, expected_f5.1);
            emit(
                json!({
                    "case": case, "seed": seed, "repeat": repeat, "arm": "f5",
                    "input_digest": input_digest, "output_digest": digest(&f5.0),
                    "report": f5.1,
                }),
                audit,
                profile,
            );
            let label = profile_label(profile, case, seed, repeat, "criterion");
            let (criterion, audit) = audited(label.as_deref(), || {
                F5Criterion::new(&polys, n_vars, degree, mask)
            });
            assert_eq!(criterion.pruned_count(), expected_criterion.pruned_count());
            assert_eq!(criterion.word_ops(), expected_criterion.word_ops());
            let mut record = json!({
                "case": case, "seed": seed, "repeat": repeat, "arm": "criterion",
                "input_digest": input_digest,
                "pruned": criterion.pruned_count(), "word_ops": criterion.word_ops(),
                "lower_rows": criterion.lower_level_rows(),
                "pruned_by_part": criterion.pruned_by_part(),
            });
            if profile {
                record["criterion_digest"] = json!(criterion_digest(&criterion, m, n_vars, degree));
            }
            emit(record, audit, profile);
        }
    }
}
