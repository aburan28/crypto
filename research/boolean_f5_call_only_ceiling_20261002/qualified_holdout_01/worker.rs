//! Native complete-call phase-cost screen for Boolean matrix F5.
//! This computes an optimistic conditional ceiling, not a candidate speedup.

use crypto_lib::cryptanalysis::koblitz_groebner::matrix_f4_f2;
use crypto_lib::cryptanalysis::matrix_f5_f2::{
    canonical_row_space_fingerprint, matrix_f5_f2_with_form_timed, F5OutputForm, F5Report,
    F5Timings,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, BufWriter, Write};
use std::process::Command;
use std::time::Instant;

const PROTOCOL_BYTES: &[u8] = include_bytes!("protocol.json");
const SOURCE_BYTES: &[u8] = include_bytes!("worker.rs");
const VERIFY_BYTES: &[u8] = include_bytes!("verify.rs");
const FROZEN_PROTOCOL_SHA256: &str =
    "8ec427456cc5dce1bcd3964adc7755030d976693ee6470b8999f06417706b098";

mod verify;

#[derive(Deserialize)]
struct Protocol {
    schema_version: u32,
    degree_bound: u32,
    output_form: String,
    direct_pack_fallback_allowed: bool,
    require_direct_unpack: bool,
    env: BTreeMap<String, String>,
    variables: Vec<u8>,
    batches: Vec<usize>,
    families: Vec<String>,
    discovery_seeds: Vec<u64>,
    holdout_seeds: Vec<u64>,
    repetitions: usize,
    bootstrap_resamples: usize,
    bootstrap_seed: u64,
    primary_variables: Vec<u8>,
    primary_batch: usize,
    primary_groups: usize,
    upper_bound_greater_than: f64,
    worker_seconds: u64,
    evidence_cap_bytes: usize,
    reservation: String,
}

#[derive(Clone)]
struct System {
    n: u8,
    quadratic: Vec<Vec<u64>>,
    affine: Vec<u64>,
}

#[derive(Clone, Serialize, Deserialize)]
struct Reference {
    report: F5Report,
    output_digest: String,
    returned_rows: usize,
    returned_terms: usize,
    direct_pack_used: bool,
    direct_unpack_used: bool,
}

#[derive(Serialize)]
struct CallSample {
    input_index: usize,
    call_ns: u128,
    destruction_ns: u128,
    total_ns: u128,
    criterion_ns: u64,
    build_ns: u64,
    reduce_ns: u64,
    unpack_ns: u64,
    validation_ns: u128,
    returned_rows: usize,
    returned_terms: usize,
    report: F5Report,
    output_digest: String,
    direct_pack_used: bool,
    direct_unpack_used: bool,
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn next(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn fixture_with_state(n: u8, seed: u64) -> (System, u64) {
    assert!([12, 16, 20, 24].contains(&n));
    let mut state = seed;
    let mut quadratic = Vec::with_capacity(n as usize);
    for _ in 0..n {
        let mut terms = BTreeSet::new();
        while terms.len() < usize::from(2 * n) {
            let a = next(&mut state) % u64::from(n);
            let b = next(&mut state) % u64::from(n);
            if a != b {
                terms.insert((1u64 << a) | (1u64 << b));
            }
        }
        quadratic.push(terms.into_iter().collect());
    }
    (
        System {
            n,
            quadratic,
            affine: vec![0; n as usize],
        },
        state,
    )
}

fn assignments(n: u8, seed: u64, batch: usize, family: &str) -> Vec<System> {
    assert!([2, 8, 32].contains(&batch));
    assert!(["independent_affine", "walk_affine"].contains(&family));
    let (base, mut state) = fixture_with_state(n, seed);
    let mut current = base.clone();
    let mut out = Vec::with_capacity(batch);
    out.push(base.clone());
    for _ in 1..batch {
        let mut input = base.clone();
        if family == "independent_affine" {
            for affine in &mut input.affine {
                for slot in 0..=n {
                    if next(&mut state) & 1 != 0 {
                        *affine |= 1u64 << slot;
                    }
                }
            }
        } else {
            for affine in &mut current.affine {
                let slot = next(&mut state) % u64::from(n + 1);
                *affine ^= 1u64 << slot;
            }
            input.affine.clone_from(&current.affine);
        }
        out.push(input);
    }
    out
}

fn polynomials(system: &System) -> Vec<F2BoolPoly> {
    system
        .quadratic
        .iter()
        .zip(&system.affine)
        .map(|(quadratic, &affine)| {
            let mut terms = quadratic
                .iter()
                .copied()
                .map(F2BoolMono::from_mask)
                .collect::<Vec<_>>();
            if affine & 1 != 0 {
                terms.push(F2BoolMono::from_mask(0));
            }
            for variable in 0..system.n {
                if affine & (1u64 << (variable + 1)) != 0 {
                    terms.push(F2BoolMono::from_mask(1u64 << variable));
                }
            }
            F2BoolPoly::from_monos(terms, system.n as usize)
        })
        .collect()
}

fn output_digest(rows: &[F2BoolPoly]) -> (String, usize) {
    let mut hasher = Sha256::new();
    hasher.update(b"boolean-f5-returned-rows-v1");
    hasher.update((rows.len() as u64).to_le_bytes());
    let mut terms = 0;
    for row in rows {
        hasher.update((row.n_vars as u64).to_le_bytes());
        hasher.update((row.terms.len() as u64).to_le_bytes());
        for monomial in &row.terms {
            hasher.update(monomial.mask.to_le_bytes());
            terms += 1;
        }
    }
    (format!("{:x}", hasher.finalize()), terms)
}

fn route() -> (BTreeMap<String, String>, Protocol) {
    assert_eq!(sha256(PROTOCOL_BYTES), FROZEN_PROTOCOL_SHA256);
    let protocol: Protocol = serde_json::from_slice(PROTOCOL_BYTES).expect("protocol");
    assert_eq!(protocol.schema_version, 1);
    assert_eq!(protocol.degree_bound, 4);
    assert_eq!(protocol.output_form, "Echelon");
    assert!(protocol.direct_pack_fallback_allowed);
    assert!(protocol.require_direct_unpack);
    assert_eq!(protocol.variables, [12, 16, 20, 24]);
    assert_eq!(protocol.batches, [2, 8, 32]);
    assert_eq!(protocol.families, ["independent_affine", "walk_affine"]);
    assert_eq!(protocol.repetitions, 7);
    assert_eq!(protocol.bootstrap_resamples, 4000);
    assert_eq!(protocol.primary_variables, [24]);
    assert_eq!(protocol.primary_batch, 32);
    assert_eq!(protocol.primary_groups, 4);
    assert_eq!(protocol.upper_bound_greater_than, 2.0);
    assert_eq!(protocol.reservation, "all_smt_siblings_of_logical_cpu_2");
    (protocol.env.clone(), protocol)
}

fn check_route_environment(route: &BTreeMap<String, String>) {
    for (key, value) in route {
        assert_eq!(
            std::env::var(key).as_deref(),
            Ok(value.as_str()),
            "route {key}"
        );
    }
    for (key, _) in std::env::vars() {
        if (key.starts_with("KIC_F5_") || key.starts_with("KIC_GF2_")) && !route.contains_key(&key)
        {
            panic!("unregistered F5/GF2 route option {key}");
        }
    }
}

fn check_host() {
    #[cfg(target_arch = "x86_64")]
    {
        assert_eq!(std::env::consts::OS, "linux", "Linux required");
        assert!(std::arch::is_x86_feature_detected!("avx2"), "AVX2 required");
    }
    #[cfg(not(target_arch = "x86_64"))]
    panic!("qualified phase requires Linux x86-64 AVX2");
}

fn call_f5(polys: &[F2BoolPoly], n: u8) -> (Vec<F2BoolPoly>, F5Report, F5Timings) {
    matrix_f5_f2_with_form_timed(polys, n as usize, 4, F5OutputForm::Echelon)
        .expect("F5 call must complete")
}

fn reference(polys: &[F2BoolPoly], n: u8) -> Reference {
    let (rows, report, timings) = call_f5(polys, n);
    assert!(timings.direct_unpack_used, "n={n} direct unpack required");
    let (output_digest, returned_terms) = output_digest(&rows);
    Reference {
        report,
        output_digest,
        returned_rows: rows.len(),
        returned_terms,
        direct_pack_used: timings.direct_pack_used,
        direct_unpack_used: timings.direct_unpack_used,
    }
}

fn f4_cross_check(polys: &[F2BoolPoly], n: u8) -> Option<u64> {
    if n > 16 {
        return None;
    }
    let f4 = matrix_f4_f2(polys, n as usize, 4).expect("small F4 reference");
    let (f5, _, timings) = call_f5(polys, n);
    assert!(timings.direct_unpack_used, "n={n} direct unpack required");
    let f4_fingerprint = canonical_row_space_fingerprint(&f4).expect("F4 row space");
    let f5_fingerprint = canonical_row_space_fingerprint(&f5).expect("F5 row space");
    assert_eq!(f4_fingerprint, f5_fingerprint, "small F4/F5 row spaces");
    Some(f4_fingerprint)
}

fn timed_call(polys: &[F2BoolPoly], n: u8, index: usize, expected: &Reference) -> CallSample {
    let start = Instant::now();
    let (rows, report, timings) = call_f5(polys, n);
    let call_ns = start.elapsed().as_nanos();
    let check = Instant::now();
    assert!(timings.direct_unpack_used);
    assert_eq!(timings.direct_pack_used, expected.direct_pack_used);
    assert_eq!(report, expected.report);
    let (digest, terms) = output_digest(&rows);
    assert_eq!(digest, expected.output_digest);
    assert_eq!(rows.len(), expected.returned_rows);
    assert_eq!(terms, expected.returned_terms);
    let validation_ns = check.elapsed().as_nanos();
    let returned_rows = rows.len();
    let destruction_start = Instant::now();
    drop(rows);
    let destruction_ns = destruction_start.elapsed().as_nanos();
    let total_ns = call_ns + destruction_ns;
    let phase_sum = u128::from(timings.criterion_ns)
        + u128::from(timings.build_ns)
        + u128::from(timings.reduce_ns)
        + u128::from(timings.unpack_ns);
    assert!(phase_sum <= call_ns, "exclusive call-phase accounting");
    assert!(u128::from(timings.build_ns) + u128::from(timings.reduce_ns) < total_ns);
    CallSample {
        input_index: index,
        call_ns,
        destruction_ns,
        total_ns,
        criterion_ns: timings.criterion_ns,
        build_ns: timings.build_ns,
        reduce_ns: timings.reduce_ns,
        unpack_ns: timings.unpack_ns,
        validation_ns,
        returned_rows,
        returned_terms: terms,
        report,
        output_digest: digest,
        direct_pack_used: timings.direct_pack_used,
        direct_unpack_used: timings.direct_unpack_used,
    }
}

fn sample_batch(
    polys: &[Vec<F2BoolPoly>],
    references: &[Reference],
    n: u8,
    repetition: usize,
    arm: usize,
) -> Value {
    let mut calls = Vec::with_capacity(polys.len());
    let mut call_total = 0u128;
    let mut destruction_total = 0u128;
    let mut total = 0u128;
    let mut criterion = 0u128;
    let mut build = 0u128;
    let mut reduce = 0u128;
    let mut unpack = 0u128;
    let mut validate = 0u128;
    for position in 0..polys.len() {
        let index = (position + repetition + arm) % polys.len();
        let call = timed_call(&polys[index], n, index, &references[index]);
        call_total += call.call_ns;
        destruction_total += call.destruction_ns;
        total += call.total_ns;
        criterion += u128::from(call.criterion_ns);
        build += u128::from(call.build_ns);
        reduce += u128::from(call.reduce_ns);
        unpack += u128::from(call.unpack_ns);
        validate += call.validation_ns;
        calls.push(call);
    }
    let residual = total - build - reduce;
    let ceiling = total as f64 / residual as f64;
    json!({
        "kind":"sample", "repetition":repetition,
        "arm":if arm==0 {"aa_a"} else {"aa_b"},
        "call_ns":call_total,"destruction_ns":destruction_total,
        "total_ns":total,"criterion_ns":criterion,
        "build_ns":build,"reduce_ns":reduce,"unpack_ns":unpack,
        "validation_ns":validate,"optimistic_ceiling":ceiling,
        "calls":calls
    })
}

fn cell(n: u8, seed: u64, batch: usize, family: &str, repetitions: usize, split: &str) {
    let systems = assignments(n, seed, batch, family);
    let polys = systems.iter().map(polynomials).collect::<Vec<_>>();
    let references = polys
        .iter()
        .map(|input| reference(input, n))
        .collect::<Vec<_>>();
    let f4_fingerprint = f4_cross_check(&polys[0], n);
    let cell_id = format!("n{n}-{split}-{seed}-{family}-b{batch}");
    println!(
        "{}",
        json!({
            "kind":"fixture","cell":cell_id,"n":n,"seed":seed,
            "family":family,"batch":batch,
            "quadratic":systems[0].quadratic,
        "affine":systems.iter().map(|s|&s.affine).collect::<Vec<_>>(),
        "references":references,
        "small_f4_fingerprint":f4_fingerprint
        })
    );
    for repetition in 0..repetitions {
        for arm in 0..2 {
            let mut sample = sample_batch(&polys, &references, n, repetition, arm);
            sample["cell"] = json!(cell_id);
            println!("{sample}");
        }
    }
}

fn campaign(phase: &str, protocol_path: &str) {
    let (route_env, protocol) = route();
    check_host();
    check_route_environment(&route_env);
    let supplied = std::fs::read(protocol_path).expect("read protocol");
    assert_eq!(supplied, PROTOCOL_BYTES, "changed protocol bytes");
    let (split, seeds) = match phase {
        "discovery" => ("discovery", protocol.discovery_seeds),
        "holdout" => ("holdout", protocol.holdout_seeds),
        _ => panic!("phase"),
    };
    let mut writer = BufWriter::new(io::stdout().lock());
    let header = json!({
        "kind":"campaign","phase":phase,
        "protocol_sha256":sha256(PROTOCOL_BYTES),
        "source_sha256":sha256(SOURCE_BYTES),
        "verifier_sha256":sha256(VERIFY_BYTES),
        "host_arch":std::env::consts::ARCH,
        "route":route_env
    });
    writeln!(writer, "{header}").unwrap();
    let mut bytes = header.to_string().len() + 1;
    let started = Instant::now();
    let executable = std::env::current_exe().expect("current executable");
    for n in protocol.variables {
        for &seed in &seeds {
            for family in &protocol.families {
                for &batch in &protocol.batches {
                    let output = Command::new(&executable)
                        .args([
                            "--cell",
                            split,
                            &n.to_string(),
                            &seed.to_string(),
                            family,
                            &batch.to_string(),
                            &protocol.repetitions.to_string(),
                        ])
                        .output()
                        .expect("launch cell");
                    assert!(
                        output.status.success(),
                        "cell failed: {}",
                        String::from_utf8_lossy(&output.stderr)
                    );
                    bytes += output.stdout.len();
                    assert!(bytes <= protocol.evidence_cap_bytes, "evidence cap");
                    writer.write_all(&output.stdout).expect("cell records");
                    assert!(
                        started.elapsed().as_secs() <= protocol.worker_seconds,
                        "worker cap; incomplete campaign censored"
                    );
                }
            }
        }
    }
    writer.flush().expect("flush campaign");
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    if args.len() == 3 && args[1] == "--wait-quiet" {
        verify::wait_quiet(&args[2]);
        return;
    }
    if args.len() == 3 && args[1] == "--seal" {
        verify::seal(&args[2]);
        return;
    }
    if args.len() == 3 && args[1] == "--verify-bundle" {
        verify::verify_bundle(&args[2]);
        return;
    }
    if args.len() == 4 && args[1] == "--failure" {
        verify::failure(&args[2], &args[3]);
        return;
    }
    if args.len() == 4 && args[1] == "--check-discovery" {
        verify::check_discovery(&args[2], &args[3]);
        return;
    }
    if args.len() == 5 && args[1] == "--verify" {
        verify::verify_run(&args[2], &args[3], &args[4]);
        return;
    }
    if args.len() == 4 && args[1] == "--campaign" {
        campaign(&args[2], &args[3]);
        return;
    }
    if args.len() == 8 && args[1] == "--cell" {
        let (route_env, _) = route();
        check_host();
        check_route_environment(&route_env);
        cell(
            args[3].parse().unwrap(),
            args[4].parse().unwrap(),
            args[6].parse().unwrap(),
            &args[5],
            args[7].parse().unwrap(),
            &args[2],
        );
        return;
    }
    panic!("usage: --campaign PHASE PROTOCOL | --cell SPLIT N SEED FAMILY BATCH REPS | --verify PHASE RAW RESULTS | --check-discovery BUNDLE BINDING | --wait-quiet OUTPUT | --failure BUNDLE REASON | --seal BUNDLE | --verify-bundle BUNDLE");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_route_and_full_core_reservation_parse() {
        let (env, protocol) = route();
        assert_eq!(env["KIC_F5_DIRECT_PACK"], "1");
        assert_eq!(env["KIC_F5_UNPACK_DIRECT"], "1");
        assert_eq!(protocol.reservation, "all_smt_siblings_of_logical_cpu_2");
    }

    #[test]
    fn fixture_and_affine_walk_keep_quadratic_core() {
        for n in [12, 16, 20, 24] {
            let systems = assignments(n, 17, 8, "walk_affine");
            for system in &systems {
                assert_eq!(system.quadratic, systems[0].quadratic);
            }
            for pair in systems.windows(2) {
                for (&a, &b) in pair[0].affine.iter().zip(&pair[1].affine) {
                    assert_eq!((a ^ b).count_ones(), 1);
                }
            }
        }
    }

    #[test]
    fn optimistic_ceiling_is_at_least_one() {
        let total = 1000u128;
        let build = 100u128;
        let reduce = 300u128;
        assert_eq!(total as f64 / (total - build - reduce) as f64, 5.0 / 3.0);
    }

    #[test]
    fn development_route_matches_small_f4_and_phase_accounting() {
        if std::env::var("KIC_F5_DIRECT_PACK").as_deref() != Ok("1")
            || std::env::var("KIC_F5_UNPACK_DIRECT").as_deref() != Ok("1")
        {
            return;
        }
        for n in [12, 16] {
            let inputs = assignments(n, 17, 2, "independent_affine");
            let polys = polynomials(&inputs[0]);
            assert!(f4_cross_check(&polys, n).is_some());
            let expected = reference(&polys, n);
            let sample = timed_call(&polys, n, 0, &expected);
            assert!(sample.call_ns > u128::from(sample.build_ns + sample.reduce_ns));
            assert_eq!(sample.total_ns, sample.call_ns + sample.destruction_ns);
            assert_eq!(sample.output_digest, expected.output_digest);
        }
    }
}
