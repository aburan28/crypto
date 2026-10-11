//! Native allocation audit of complete Boolean matrix-F4/F5 calls.
//! No wall-clock measurements: ordinary hosts do not establish isolation.
//! Run with `cargo run --release --example f4_f5_row_build_bench -- [repeats]`.

use crypto_lib::cryptanalysis::koblitz_groebner::matrix_f4_f2_counted;
use crypto_lib::cryptanalysis::matrix_f5_f2::{matrix_f5_f2, F5Criterion};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use serde_json::json;
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

struct CountingAllocator;

static ENABLED: AtomicBool = AtomicBool::new(false);
static ALLOCATIONS: AtomicU64 = AtomicU64::new(0);
static BYTES: AtomicU64 = AtomicU64::new(0);

fn record(bytes: usize) {
    if ENABLED.load(Ordering::Relaxed) {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        BYTES.fetch_add(bytes as u64, Ordering::Relaxed);
    }
}

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        record(layout.size());
        unsafe { System.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        record(layout.size());
        unsafe { System.alloc_zeroed(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        record(size);
        unsafe { System.realloc(ptr, layout, size) }
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

fn audited<T>(run: impl FnOnce() -> T) -> (T, u64, u64) {
    ALLOCATIONS.store(0, Ordering::Relaxed);
    BYTES.store(0, Ordering::Relaxed);
    ENABLED.store(true, Ordering::Relaxed);
    let result = run();
    ENABLED.store(false, Ordering::Relaxed);
    (
        result,
        ALLOCATIONS.load(Ordering::Relaxed),
        BYTES.load(Ordering::Relaxed),
    )
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

fn main() {
    let repeats = std::env::args()
        .nth(1)
        .map(|s| s.parse::<usize>().expect("positive repeat count"))
        .unwrap_or(3);
    assert!(repeats > 0);
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .unwrap();
    for (case, vars, m, degree) in [
        ("n8_d3", (0..8).collect::<Vec<_>>(), 8, 3),
        ("n8_d4", (0..8).collect::<Vec<_>>(), 8, 4),
        ("n10_d4", (0..10).collect::<Vec<_>>(), 10, 4),
        ("n10_d5", (0..10).collect::<Vec<_>>(), 10, 5),
    ] {
        for seed in [0x2545_f491_4f6c_dd1d, 0x9e37_79b9_7f4a_7c15] {
            let polys = system(&vars, m, seed);
            let n_vars = polys[0].n_vars;
            let mask = (1u64 << n_vars) - 1;
            let input_digest = digest(&polys);
            // Warm lazy configuration and worker startup outside the audit.
            let expected_f4 = matrix_f4_f2_counted(&polys, n_vars, degree).unwrap();
            let expected_f5 = matrix_f5_f2(&polys, n_vars, degree).unwrap();
            assert_eq!(expected_f4.0, expected_f5.0);
            let expected_criterion = F5Criterion::new(&polys, n_vars, degree, mask);
            for repeat in 0..repeats {
                let (f4, allocations, bytes) =
                    audited(|| matrix_f4_f2_counted(&polys, n_vars, degree).unwrap());
                assert_eq!(f4, expected_f4);
                println!("{}", json!({
                    "case": case, "seed": seed, "repeat": repeat, "arm": "f4",
                    "input_digest": input_digest, "output_digest": digest(&f4.0),
                    "allocations": allocations, "allocated_bytes": bytes,
                    "rank": f4.0.len(), "word_ops": f4.1,
                }));
                let (f5, allocations, bytes) =
                    audited(|| matrix_f5_f2(&polys, n_vars, degree).unwrap());
                assert_eq!(f5.0, expected_f5.0);
                assert_eq!(f5.1, expected_f5.1);
                println!("{}", json!({
                    "case": case, "seed": seed, "repeat": repeat, "arm": "f5",
                    "input_digest": input_digest, "output_digest": digest(&f5.0),
                    "allocations": allocations, "allocated_bytes": bytes,
                    "report": f5.1,
                }));
                let (criterion, allocations, bytes) =
                    audited(|| F5Criterion::new(&polys, n_vars, degree, mask));
                assert_eq!(criterion.pruned_count(), expected_criterion.pruned_count());
                assert_eq!(criterion.word_ops(), expected_criterion.word_ops());
                println!("{}", json!({
                    "case": case, "seed": seed, "repeat": repeat, "arm": "criterion",
                    "input_digest": input_digest,
                    "allocations": allocations, "allocated_bytes": bytes,
                    "pruned": criterion.pruned_count(), "word_ops": criterion.word_ops(),
                    "lower_rows": criterion.lower_level_rows(),
                    "pruned_by_part": criterion.pruned_by_part(),
                }));
            }
        }
    }
}
