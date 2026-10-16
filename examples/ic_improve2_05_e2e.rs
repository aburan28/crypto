//! IC improvement round 2, item 05: one end-to-end Koblitz index-calculus
//! run (relation collection, linear algebra, verification) through the
//! production entry point with the Gröbner decomposition oracle.  The arm is
//! selected by the environment (`KIC_LINEAR_HARVEST`, see
//! `koblitz_groebner::linear_harvest`), so each arm is its own process;
//! the arm label printed is the variable's value (`unset` when absent).
//!
//! Prints one JSON line: the arm, the verified logarithm, relations, trials,
//! the per-attempt disposition sequence (as a hash and counts), the solver
//! reductions, and CPU / relation-collection time.
//!
//! ```bash
//! KIC_LINEAR_HARVEST=0 cargo run --release --example ic_improve2_05_e2e -- --n 23 --ell 12 --m 2
//! cargo run --release --example ic_improve2_05_e2e -- --n 23 --ell 12 --m 2
//! ```
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, koblitz_index_calculus_dlp_with_factor_base_and_progress,
    DecompositionStrategy, KoblitzCurve, KoblitzIcEvent, KoblitzIcOptions,
    KoblitzRelationAttemptDisposition,
};
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

fn tcpu() -> f64 {
    let mut ts = libc::timespec {
        tv_sec: 0,
        tv_nsec: 0,
    };
    unsafe {
        libc::clock_gettime(libc::CLOCK_THREAD_CPUTIME_ID, &mut ts);
    }
    ts.tv_sec as f64 + ts.tv_nsec as f64 * 1e-9
}

/// Process CPU (user + system, every thread) in seconds.
fn pcpu() -> f64 {
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    unsafe {
        libc::getrusage(libc::RUSAGE_SELF, &mut ru);
    }
    let tv = |t: libc::timeval| t.tv_sec as f64 + t.tv_usec as f64 * 1e-6;
    tv(ru.ru_utime) + tv(ru.ru_stime)
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    let flag = |k: &str| {
        argv.iter()
            .position(|x| x == k)
            .map(|i| argv[i + 1].clone())
    };
    let n: u32 = flag("--n").map(|v| v.parse().unwrap()).unwrap_or(17);
    let ell: u32 = flag("--ell").map(|v| v.parse().unwrap()).unwrap_or(6);
    let m: usize = flag("--m").map(|v| v.parse().unwrap()).unwrap_or(3);
    let seed: u64 = flag("--seed")
        .map(|v| v.parse().unwrap())
        .unwrap_or(20261015);
    let max_trials: usize = flag("--max-trials")
        .map(|v| v.parse().unwrap())
        .unwrap_or(20_000);
    let arm = std::env::var("KIC_LINEAR_HARVEST").unwrap_or_else(|_| "unset".into());
    let (h0, _) = crypto_lib::cryptanalysis::koblitz_groebner::linear_harvest_counts_thread();
    let kc = KoblitzCurve::new(0, n)
        .or_else(|| KoblitzCurve::new(1, n))
        .expect("curve");
    let fb = build_standard_subspace_factor_base(&kc, ell).expect("factor base");
    let ord = kc.subgroup_order.to_u64_digits()[0];
    let d = StdRng::seed_from_u64(seed ^ n as u64).gen_range(1..ord);
    let q = kc.mul(kc.generator(), &BigUint::from(d));
    let opts = KoblitzIcOptions {
        m,
        strategy: DecompositionStrategy::Groebner,
        seed,
        max_trials,
        ..Default::default()
    };
    let mut dispositions: Vec<KoblitzRelationAttemptDisposition> = Vec::new();
    let c0 = tcpu();
    let p0 = pcpu();
    let report =
        koblitz_index_calculus_dlp_with_factor_base_and_progress(&kc, &q, &fb, &opts, &mut |e| {
            if let KoblitzIcEvent::RelationAttemptFinished { disposition, .. } = e {
                dispositions.push(disposition);
            }
        });
    let cpu = tcpu() - c0;
    let proc_cpu = pcpu() - p0;
    let (h1, _) = crypto_lib::cryptanalysis::koblitz_groebner::linear_harvest_counts_thread();
    let harvests = h1 - h0;
    let mut h = DefaultHasher::new();
    for x in &dispositions {
        format!("{x:?}").hash(&mut h);
    }
    let count = |want: KoblitzRelationAttemptDisposition| {
        dispositions.iter().filter(|&&x| x == want).count()
    };
    let (log_ok, relations, trials, reductions, rc_ns) = match &report {
        Some(r) => (
            r.log
                .as_ref()
                .map(|l| *l == BigUint::from(d) % &kc.subgroup_order),
            r.relations,
            r.trials,
            r.reductions,
            r.relation_collection_ns,
        ),
        None => (None, 0, 0, 0, 0),
    };
    println!(
        r#"{{"arm":"{arm}","n":{n},"ell":{ell},"m":{m},"seed":{seed},"d":{d},"log_ok":{},"relations":{relations},"trials":{trials},"attempts":{},"found":{},"refuted":{},"unknown":{},"dispositions_hash":"{:016x}","reductions":{reductions},"harvests":{harvests},"relation_collection_s":{:.4},"cpu_s":{cpu:.4},"proc_cpu_s":{proc_cpu:.4}}}"#,
        match log_ok {
            Some(b) => b.to_string(),
            None => "null".into(),
        },
        dispositions.len(),
        count(KoblitzRelationAttemptDisposition::RelationFound),
        count(KoblitzRelationAttemptDisposition::Refuted),
        count(KoblitzRelationAttemptDisposition::Unknown),
        h.finish(),
        rc_ns as f64 * 1e-9,
    );
}
