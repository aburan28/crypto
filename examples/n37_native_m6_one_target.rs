#![recursion_limit = "256"]
//! Cold, one-public-target n37 descendant-native bounded residual search.
//! Protocol: research/notes/ecc2k130/n37_native_one_target_20261003/PROTOCOL.md
//!
//! `cargo run --release --example n37_native_m6_one_target -- INDEX RESULT.json`

use crypto_lib::binary_ecc::{F2mElement, F2mPoly, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_field_basis::BinaryFieldBasis;
use crypto_lib::cryptanalysis::binary_velu::{
    elt, isogeny_from_kernel, to_u64, velu_point_map, Curve,
};
use crypto_lib::cryptanalysis::ec_index_calculus::gaussian_eliminate_mod_n;
use crypto_lib::cryptanalysis::koblitz_fast::FastPoint;
use crypto_lib::cryptanalysis::koblitz_relation_solver::{IncrementalRelationSolver, RowStatus};
use crypto_lib::cryptanalysis::native_signed_mitm::{MitmCounts, NativeSignedMitm, SignedWitness};
use crypto_lib::hash::sha256::sha256;
use flate2::read::GzDecoder;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::cell::Cell;
use std::collections::HashSet;
use std::fs;
use std::io::Read;
use std::path::Path;
use std::time::Instant;

const NOTE: &str = "research/notes/ecc2k130/n37_native_one_target_20261003";
const BASE_PATH: &str = "research/notes/ecc2k130/n37_native_basis_bridge_20261002/NATIVE42.json";
const BASE_SHA: &str = "bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c";
const FROZEN_PATH: &str = "research/notes/ecc2k130/disjoint_cold_v2_20261001/FROZEN.json";
const FROZEN_SHA: &str = "da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d";
const POINTS_PATH: &str =
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b02.points.jsonl";
const POINTS_SHA: &str = "248d405e834b2ee7c27752e5742ca568d21156d63e981fa49f362cdaab2f9cd9";
const ARCHIVE_PATH: &str = "research/koblitz_isogeny_descent_37_results_20260925/raw.json.gz";
const ARCHIVE_SHA: &str = "eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90";
const ORDER: u64 = 137_439_487_532;
const R: u64 = 230_603_167;
const K: usize = 42;
const MAX_PROBES: usize = 256;
const SCALAR_SEED: u64 = 0x6e33_376d_365f_3031;
const SHIFT_SEED: u64 = 0x6e33_375f_7265_7331;
const SHIFT_CAP: usize = 16;

fn read_verified(path: &str, expected_sha: &str) -> Result<Vec<u8>, String> {
    let bytes = fs::read(path).map_err(|error| format!("read {path}: {error}"))?;
    let actual = hex::encode(sha256(&bytes));
    if actual != expected_sha {
        return Err(format!("SHA-256 mismatch for {path}: {actual}"));
    }
    Ok(bytes)
}

fn archive_field_element(expression: &str) -> Result<u64, String> {
    let mut value = 0u64;
    for monomial in expression.trim_matches(['(', ')']).split(" + ") {
        let bit = match monomial {
            "1" => 0,
            "a" => 1,
            power => power
                .strip_prefix("a^")
                .ok_or_else(|| format!("invalid archived field monomial {power}"))?
                .parse::<u32>()
                .map_err(|error| format!("archived field exponent: {error}"))?,
        };
        if bit >= 37 {
            return Err("archived field coefficient exceeds degree 37".into());
        }
        value ^= 1u64 << bit;
    }
    Ok(value)
}

fn archive_x_degree(suffix: &str) -> Result<usize, String> {
    if suffix.is_empty() {
        Ok(1)
    } else {
        suffix
            .strip_prefix('^')
            .ok_or("invalid archived x power".to_string())?
            .parse()
            .map_err(|error| format!("archived x exponent: {error}"))
    }
}

fn archive_kernel(expression: &str) -> Result<F2mPoly, String> {
    let mut terms = Vec::new();
    let mut start = 0;
    let mut depth = 0i32;
    for (index, symbol) in expression.char_indices() {
        match symbol {
            '(' => depth += 1,
            ')' => depth -= 1,
            '+' if depth == 0 => {
                terms.push(expression[start..index].trim());
                start = index + 1;
            }
            _ => {}
        }
        if depth < 0 {
            return Err("unbalanced archived kernel".into());
        }
    }
    if depth != 0 {
        return Err("unbalanced archived kernel".into());
    }
    terms.push(expression[start..].trim());
    let mut coeffs = vec![F2mElement::zero(37); 37];
    for term in terms {
        let (coefficient, degree) = if let Some((coefficient, suffix)) = term.rsplit_once("*x") {
            (
                archive_field_element(coefficient)?,
                archive_x_degree(suffix)?,
            )
        } else if let Some(suffix) = term.strip_prefix('x') {
            (1, archive_x_degree(suffix)?)
        } else {
            (archive_field_element(term)?, 0)
        };
        if degree > 36 {
            return Err("archived kernel exceeds degree 36".into());
        }
        coeffs[degree] = elt(to_u64(&coeffs[degree]) ^ coefficient, 37);
    }
    let kernel = F2mPoly::from_coeffs(coeffs, 37);
    if kernel.degree() != Some(36) || to_u64(&kernel.lead()) != 1 {
        return Err("archived kernel is not monic degree 36".into());
    }
    Ok(kernel)
}

fn point_from_value(value: &Value) -> Result<FastPoint, String> {
    let words: [u64; 2] = serde_json::from_value(value.clone())
        .map_err(|error| format!("point coordinates: {error}"))?;
    Ok(FastPoint::affine(words[0], words[1]))
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut value = *state;
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn row_from_witness(witness: &SignedWitness, rhs: u64) -> Vec<u64> {
    let mut row = Vec::with_capacity(K + 2);
    for &coefficient in &witness.coefficients {
        row.push((coefficient as i64).rem_euclid(R as i64) as u64);
    }
    row.push(0); // The incremental solver's extra target column is unused.
    row.push(rhs);
    row
}

fn target_log(witness: &SignedWitness, base_logs: &[u64]) -> u64 {
    witness
        .coefficients
        .iter()
        .zip(base_logs)
        .map(|(&coefficient, &log)| coefficient as i128 * log as i128)
        .sum::<i128>()
        .rem_euclid(R as i128) as u64
}

fn accumulate(total: &mut MitmCounts, added: MitmCounts) {
    total.query_adds += added.query_adds;
    total.query_batches += added.query_batches;
    total.lookups += added.lookups;
    total.witness_adds += added.witness_adds;
}

fn elapsed_ms(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1_000.0
}

fn run(output: &Path, target_index: usize) -> Result<(), String> {
    let total_start = Instant::now();
    let setup_start = Instant::now();
    let mut setup_scalar_muls = 0u64;
    let base: Value = serde_json::from_slice(&read_verified(BASE_PATH, BASE_SHA)?)
        .map_err(|error| format!("base JSON: {error}"))?;
    let frozen: Value = serde_json::from_slice(&read_verified(FROZEN_PATH, FROZEN_SHA)?)
        .map_err(|error| format!("frozen JSON: {error}"))?;
    let point_bytes = read_verified(POINTS_PATH, POINTS_SHA)?;
    let archive_bytes = read_verified(ARCHIVE_PATH, ARCHIVE_SHA)?;
    let spec = &frozen["specs"]["n37_L1024"];
    if base["subgroup_order"].as_u64() != Some(R)
        || base["group_order"].as_u64() != Some(ORDER)
        || spec["subgroup_order"].as_u64() != Some(R)
        || base["source_generator_image"].as_u64() != Some(10_156_182_909)
    {
        return Err("frozen group parameters changed".into());
    }
    let source_modulus = IrreduciblePoly {
        degree: 37,
        low_terms: serde_json::from_value(spec["field_modulus_low_terms"].clone())
            .map_err(|error| format!("source modulus: {error}"))?,
    };
    let archive_modulus = IrreduciblePoly {
        degree: 37,
        low_terms: vec![0, 1, 2, 3, 4, 5],
    };
    let bridge =
        BinaryFieldBasis::from_generator_image(&source_modulus, &archive_modulus, 10_156_182_909)?;
    let control = Curve::with_order(37, &source_modulus, 0, 1, ORDER).ok_or("control curve")?;
    let archived_source =
        Curve::with_order(37, &archive_modulus, 0, 1, ORDER).ok_or("archive curve")?;
    let mut decompressed = String::new();
    GzDecoder::new(&archive_bytes[..])
        .read_to_string(&mut decompressed)
        .map_err(|error| format!("archive gzip: {error}"))?;
    let archive: Value =
        serde_json::from_str(&decompressed).map_err(|error| format!("archive JSON: {error}"))?;
    let certificate = &archive["isogeny_certificate"];
    if certificate["degree"].as_u64() != Some(73)
        || certificate["field_modulus"] != "x^37 + x^5 + x^4 + x^3 + x^2 + x + 1"
    {
        return Err("archived degree-73 certificate changed".into());
    }
    let kernel = archive_kernel(
        certificate["kernel_polynomial"]
            .as_str()
            .ok_or("archived kernel polynomial")?,
    )?;
    let isogeny = isogeny_from_kernel(&archived_source, kernel, 73);
    let leaf_a6 = archive_field_element(
        certificate["target_normalized_coefficients"]["a6"]
            .as_str()
            .ok_or("archived leaf a6")?,
    )?;
    if isogeny.a6_codomain != leaf_a6 {
        return Err("reconstructed leaf coefficient mismatch".into());
    }
    let leaf = Curve::with_order(37, &archive_modulus, 0, leaf_a6, ORDER).ok_or("leaf curve")?;
    let transport_calls = Cell::new(0u64);
    let transport = |point: FastPoint| -> Result<FastPoint, String> {
        transport_calls.set(transport_calls.get() + 1);
        velu_point_map(
            bridge.map_point(point).ok_or("basis map point")?,
            &isogeny.kernel,
            37,
            &archive_modulus,
        )
        .ok_or("undefined degree-73 point map".into())
    };
    let source_generator = point_from_value(&spec["generator"])?;
    setup_scalar_muls += 1;
    if !control.fast.is_on_curve(source_generator)
        || !control.fast.mul_u64(source_generator, R).infinity
    {
        return Err("invalid control generator".into());
    }
    let leaf_generator = transport(source_generator)?;
    setup_scalar_muls += 1;
    if !leaf.fast.is_on_curve(leaf_generator) || !leaf.fast.mul_u64(leaf_generator, R).infinity {
        return Err("invalid leaf generator".into());
    }
    let manifest_rows = base["rows"].as_array().ok_or("base rows")?;
    if manifest_rows.len() != K {
        return Err("native base is not K=42".into());
    }
    let mut points = Vec::new();
    let mut scanned_x = 0u64;
    for x in 1..4096 {
        scanned_x = x;
        let Some(point) = leaf.points_with_x(x).into_iter().next() else {
            continue;
        };
        setup_scalar_muls += 1;
        let native = leaf.fast.mul_u64(point, ORDER / R);
        if native.infinity || points.iter().any(|other: &FastPoint| other.x == native.x) {
            continue;
        }
        let row = &manifest_rows[points.len()];
        if row["index"].as_u64() != Some(points.len() as u64)
            || point_from_value(&row["leaf"])? != native
        {
            return Err("native base selection differs from frozen manifest".into());
        }
        let archived_pullback = point_from_value(&row["archive_pullback"])?;
        let control_pullback = point_from_value(&row["control_pullback"])?;
        setup_scalar_muls += 1;
        if !archived_source.fast.is_on_curve(archived_pullback)
            || !control.fast.is_on_curve(control_pullback)
            || bridge.map_point(control_pullback) != Some(archived_pullback)
            || transport(control_pullback)? != native
            || !leaf.fast.mul_u64(native, R).infinity
        {
            return Err(format!("invalid native base row {}", points.len()));
        }
        points.push(native);
        if points.len() == K {
            break;
        }
    }
    if points.len() != K {
        return Err("native base selection exhausted scan".into());
    }
    // Input loading is not target-dependent arithmetic. The producer never
    // opens the separate b02 fixture containing known-answer scalars.
    let public_text =
        String::from_utf8(point_bytes).map_err(|error| format!("target text: {error}"))?;
    let point_line = public_text
        .lines()
        .nth(target_index)
        .ok_or("target file is shorter than requested index")?;
    let public_words: [u64; 2] =
        serde_json::from_str(point_line).map_err(|error| format!("target point: {error}"))?;
    let setup_ms = elapsed_ms(setup_start);

    let table_start = Instant::now();
    let oracle = NativeSignedMitm::build(&leaf.fast, &points)?;
    let build_counts = oracle.build_counts();
    let table_ms = elapsed_ms(table_start);

    let relations_start = Instant::now();
    let mut state = SCALAR_SEED;
    let mut seen_scalars = HashSet::new();
    let mut ranker =
        IncrementalRelationSolver::new(K, &BigUint::from(R)).ok_or("incremental ranker")?;
    let mut independent_rows: Vec<(Vec<u64>, u64)> = Vec::new();
    let mut relations = Vec::new();
    let mut relation_query_counts = MitmCounts::default();
    let mut relation_scalar_muls = 0u64;
    let mut external_witness_adds = 0u64;
    for _ in 0..MAX_PROBES {
        let scalar = loop {
            let candidate = 1 + splitmix64(&mut state) % (R - 1);
            if seen_scalars.insert(candidate) {
                break candidate;
            }
        };
        let target = leaf.fast.mul_u64(leaf_generator, scalar);
        relation_scalar_muls += 1;
        let (witness, counts) = oracle.query(target, 3)?;
        accumulate(&mut relation_query_counts, counts);
        let status = if let Some(witness) = &witness {
            external_witness_adds += witness.l1_norm() as u64;
            if !oracle.verify_coefficients(target, &witness.coefficients) {
                return Err("relation coefficient replay failed".into());
            }
            let row = row_from_witness(witness, scalar);
            let status = ranker.add_row(row.clone());
            if status == RowStatus::Inconsistent {
                return Err("verified relation contradicted modular system".into());
            }
            if status == RowStatus::Independent {
                independent_rows.push((row[..K].to_vec(), scalar));
            }
            format!("{status:?}")
        } else {
            "Miss".to_string()
        };
        relations.push(json!({
            "probe_index":relations.len(), "scalar":scalar,
            "target":[target.x,target.y], "witness":witness,
            "oracle_counts":counts, "row_status":status,
            "rank_after":ranker.rank()
        }));
        if ranker.rank() == K {
            break;
        }
    }
    let relation_ms = elapsed_ms(relations_start);

    let solve_start = Instant::now();
    let mut base_logs = None;
    let mut base_verify_muls = 0u64;
    if ranker.rank() == K {
        let mut matrix: Vec<Vec<BigUint>> = independent_rows
            .iter()
            .map(|(row, _)| row.iter().copied().map(BigUint::from).collect())
            .collect();
        let mut rhs: Vec<BigUint> = independent_rows
            .iter()
            .map(|(_, scalar)| BigUint::from(*scalar))
            .collect();
        let values = gaussian_eliminate_mod_n(&mut matrix, &mut rhs, &BigUint::from(R))
            .ok_or("full-rank dense solve failed")?;
        let logs: Vec<u64> = values
            .iter()
            .map(|value| value.to_u64_digits().first().copied().unwrap_or(0))
            .collect();
        for (&point, &log) in points.iter().zip(&logs) {
            base_verify_muls += 1;
            if leaf.fast.mul_u64(leaf_generator, log) != point {
                return Err("solved base logarithm failed point verification".into());
            }
        }
        base_logs = Some(logs);
    }
    let solve_ms = elapsed_ms(solve_start);

    let shifts_start = Instant::now();
    let mut shift_state = SHIFT_SEED;
    let mut seen_shifts = HashSet::new();
    let mut shifts = Vec::with_capacity(SHIFT_CAP);
    while shifts.len() < SHIFT_CAP {
        let candidate = 1 + splitmix64(&mut shift_state) % (R - 1);
        if seen_shifts.insert(candidate) {
            shifts.push(candidate);
        }
    }
    let shift_points: Vec<FastPoint> = shifts
        .iter()
        .map(|&shift| leaf.fast.mul_u64(leaf_generator, shift))
        .collect();
    let shifts_ms = elapsed_ms(shifts_start);

    if ranker.rank() != K {
        return Err(format!(
            "target-blind relation rank {}/{}",
            ranker.rank(),
            K
        ));
    }
    let logs = base_logs.as_ref().ok_or("missing solved base logarithms")?;
    let mut direct_counts = MitmCounts::default();
    let mut residual_counts = MitmCounts::default();
    let mut shift_adds = 0u64;
    let online_start = Instant::now();
    let source = FastPoint::affine(public_words[0], public_words[1]);
    setup_scalar_muls += 1;
    if !control.fast.is_on_curve(source) || !control.fast.mul_u64(source, R).infinity {
        return Err("target is off curve or outside the prime subgroup".into());
    }
    let image = transport(source)?;
    setup_scalar_muls += 1;
    if !leaf.fast.is_on_curve(image) || !leaf.fast.mul_u64(image, R).infinity {
        return Err("target image is off curve or outside the prime subgroup".into());
    }
    let after_target_query = Instant::now();

    let mut attempts = Vec::new();
    let mut hit: Option<(usize, u64, FastPoint, SignedWitness)> = None;
    for attempt_index in 0..=SHIFT_CAP {
        let (shift, query) = if attempt_index == 0 {
            (0, image)
        } else {
            let shift = shifts[attempt_index - 1];
            shift_adds += 1;
            (shift, leaf.fast.add(image, shift_points[attempt_index - 1]))
        };
        let (witness, counts) = oracle.query(query, 3)?;
        if attempt_index == 0 {
            accumulate(&mut direct_counts, counts);
        } else {
            accumulate(&mut residual_counts, counts);
        }
        if let Some(ref found) = witness {
            hit = Some((attempt_index, shift, query, found.clone()));
        }
        attempts.push(json!({
            "attempt_index":attempt_index, "shift":shift,
            "query_leaf":[query.x,query.y], "witness":witness,
            "oracle_counts":counts
        }));
        if hit.is_some() {
            break;
        }
    }
    let after_target_pdp = Instant::now();

    if let Some((attempt_index, _, query, witness)) = &hit {
        external_witness_adds += witness.l1_norm() as u64;
        if !oracle.verify_coefficients(*query, &witness.coefficients) {
            return Err(format!(
                "target {target_index} attempt {attempt_index} coefficient replay failed"
            ));
        }
    }
    let after_relation_check = Instant::now();

    let recovered = hit.as_ref().map(|(_, shift, _, witness)| {
        let shifted_log = target_log(witness, logs);
        (shifted_log as i128 - *shift as i128).rem_euclid(R as i128) as u64
    });
    let after_target_descent = Instant::now();

    let target_verify_muls = if let Some(log) = recovered {
        if leaf.fast.mul_u64(leaf_generator, log) != image
            || control.fast.mul_u64(source_generator, log) != source
        {
            return Err(format!(
                "target {target_index} logarithm failed full-point verification"
            ));
        }
        2u64
    } else {
        0u64
    };
    let online_stop = Instant::now();
    let phase_ns = [
        (after_target_query - online_start).as_nanos() as u64,
        (after_target_pdp - after_target_query).as_nanos() as u64,
        (after_relation_check - after_target_pdp).as_nanos() as u64,
        (after_target_descent - after_relation_check).as_nanos() as u64,
        (online_stop - after_target_descent).as_nanos() as u64,
    ];
    let online_ns = (online_stop - online_start).as_nanos() as u64;
    if phase_ns.iter().sum::<u64>() != online_ns {
        return Err("online phase times do not sum exactly".into());
    }
    let target_row = json!({
        "index":target_index, "source":public_words, "leaf":[image.x,image.y],
        "attempts":attempts, "recovered_log":recovered
    });
    let direct_hits = usize::from(hit.as_ref().is_some_and(|(i, _, _, _)| *i == 0));
    let shifted_hits = usize::from(hit.as_ref().is_some_and(|(i, _, _, _)| *i > 0));
    let verified_logs = usize::from(recovered.is_some());
    let attempts_total = target_row["attempts"].as_array().ok_or("attempts")?.len();
    let result = json!({
        "schema":"n37-native-m6-one-target-v1",
        "protocol":format!("{NOTE}/PROTOCOL.md"),
        "claim_boundary":"One frozen n37 one-target IC-versus-rho panel; n41, n53, m83 and n131 are separate gates. S and speedup require calibrated complete costs.",
        "mode":"one_target",
        "inputs":[
            {"path":BASE_PATH,"sha256":BASE_SHA},
            {"path":FROZEN_PATH,"sha256":FROZEN_SHA},
            {"path":POINTS_PATH,"sha256":POINTS_SHA},
            {"path":ARCHIVE_PATH,"sha256":ARCHIVE_SHA}
        ],
        "group_order":ORDER, "subgroup_order":R, "cofactor":ORDER/R,
        "columns":K, "base_selection_scanned_x_through":scanned_x,
        "degree_73_transport_calls":transport_calls.get(),
        "target_count":1, "target_index":target_index,
        "relation_probe_seed":format!("0x{SCALAR_SEED:016x}"),
        "relation_probe_cap":MAX_PROBES,
        "rank":ranker.rank(), "relation_probes":relations.len(),
        "dependent_rows":ranker.dependent_rows(),
        "relation_misses":relations.iter().filter(|row| row["witness"].is_null()).count(),
        "base_logs":base_logs,
        "shift_seed":format!("0x{SHIFT_SEED:016x}"),
        "shift_cap":SHIFT_CAP,
        "shifts":shifts.iter().zip(&shift_points).map(|(&scalar,&point)| json!({
            "scalar":scalar,"leaf_point":[point.x,point.y]
        })).collect::<Vec<_>>(),
        "direct_hits":direct_hits, "shifted_hits":shifted_hits,
        "unresolved_targets":1-verified_logs,
        "verified_target_logs":verified_logs,
        "total_target_attempts":attempts_total,
        "maximum_target_attempts":attempts_total,
        "complete_one_target":verified_logs==1,
        "build_counts":build_counts,
        "relation_oracle_counts":relation_query_counts,
        "direct_oracle_counts":direct_counts,
        "residual_oracle_counts":residual_counts,
        "residual_shift_adds":shift_adds,
        "external_witness_replay_adds":external_witness_adds,
        "scalar_multiplications":{
            "setup_subgroup_and_base_selection":setup_scalar_muls-2,
            "relation_probes":relation_scalar_muls,
            "base_log_verification":base_verify_muls,
            "shift_precomputation":SHIFT_CAP,
            "target_validation":2,
            "target_log_verification":target_verify_muls,
            "total":setup_scalar_muls+relation_scalar_muls+base_verify_muls+SHIFT_CAP as u64+target_verify_muls
        },
        "phase_wall_ms":{
            "setup_transport":setup_ms, "half_table":table_ms,
            "relations":relation_ms, "solve_verify_base":solve_ms,
            "shift_precomputation":shifts_ms,
            "internal_excluding_output":elapsed_ms(total_start)
        },
        "online_start_ns":(online_start-total_start).as_nanos() as u64,
        "online_stop_ns":(online_stop-total_start).as_nanos() as u64,
        "online_ns":online_ns,
        "online_phase_ns":{
            "target_query":phase_ns[0],
            "target_pdp":phase_ns[1],
            "target_relation_check":phase_ns[2],
            "target_descent":phase_ns[3],
            "recovery_check":phase_ns[4]
        },
        "matched_rho_total":Value::Null,
        "S":Value::Null,
        "end_to_end_speedup":Value::Null,
        "relations":relations,
        "target":target_row
    });
    let bytes = serde_json::to_vec_pretty(&result).map_err(|error| error.to_string())?;
    fs::write(output, &bytes).map_err(|error| format!("write {}: {error}", output.display()))?;
    println!(
        "rank={}/{} index={} direct={} shifted={} verified={} attempts={} online_ns={} output={} bytes={}",
        ranker.rank(),
        K,
        target_index,
        direct_hits,
        shifted_hits,
        verified_logs,
        attempts_total,
        online_ns,
        output.display(),
        bytes.len()
    );
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let (index, output) = match args.as_slice() {
        [_, index, output] => (
            index.parse::<usize>().expect("target index"),
            output.as_str(),
        ),
        _ => panic!("usage: n37_native_m6_one_target INDEX RESULT.json"),
    };
    assert!(index < 1024);
    if let Err(error) = run(Path::new(output), index) {
        eprintln!("n37_native_m6_one_target: {error}");
        std::process::exit(1);
    }
}
