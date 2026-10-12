//! Exact private-column gate: research/f6_n83_private_quartic_core_20261005.
use std::io::Write;
use std::time::Instant;

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_standard_subspace_factor_base, cofactor_project_factor_base, KoblitzCurve,
};
use crypto_lib::cryptanalysis::wide_groebner::TwoWordFieldStructure;
use crypto_lib::cryptanalysis::wide_sixsum::{
    Mono512, Poly512, RootReduction, System512, MAX_ROOT_COLS,
};
use num_bigint::BigUint;
use serde_json::json;

const CANDIDATES_PER_ROW: usize = 32;
const RSS_LIMIT_BYTES: u64 = 7 * 1024 * 1024 * 1024;

fn peak_rss_bytes() -> Option<u64> {
    #[cfg(any(target_os = "macos", target_os = "linux"))]
    {
        let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
        if unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) } != 0 {
            return None;
        }
        let usage = unsafe { usage.assume_init() };
        #[cfg(target_os = "macos")]
        return Some(usage.ru_maxrss as u64);
        #[cfg(target_os = "linux")]
        return Some((usage.ru_maxrss as u64) * 1024);
    }
    #[cfg(not(any(target_os = "macos", target_os = "linux")))]
    None
}

fn x_of(point: &BinaryPoint) -> &F2mElement {
    match point {
        BinaryPoint::Affine { x, .. } => x,
        BinaryPoint::Infinity => panic!("affine point required"),
    }
}

fn set_element(bits: &mut Mono512, offset: usize, width: usize, element: &F2mElement) {
    let words = element.raw_bits();
    for i in 0..width {
        if words[i / 64] >> (i % 64) & 1 == 1 {
            bits.0[(offset + i) / 64] |= 1u64 << ((offset + i) % 64);
        }
    }
}

fn point_json(point: &BinaryPoint) -> serde_json::Value {
    match point {
        BinaryPoint::Infinity => json!({"infinity":true}),
        BinaryPoint::Affine { x, y } => json!({
            "x":x.to_biguint().to_str_radix(16), "y":y.to_biguint().to_str_radix(16)
        }),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert!(
        args.len() == 3 || args.len() == 4 || args.len() == 5,
        "usage: probe planted K | ordinary OFFSET(0..3) K [core16|core16wide|cubic90|cubic90all]"
    );
    let run_core16 = args.len() == 5 && (args[4] == "core16" || args[4] == "core16wide");
    let run_cubic90 = args.len() == 5 && (args[4] == "cubic90" || args[4] == "cubic90all");
    let all_cubic = args.len() == 5 && args[4] == "cubic90all";
    let wide_cap = args.len() == 5 && (args[4] == "core16wide" || run_cubic90);
    assert!(args.len() != 5 || run_core16 || run_cubic90);
    let k_arg = if run_core16 || run_cubic90 {
        &args[3]
    } else {
        args.last().unwrap()
    };
    let k: usize = k_arg.parse().expect("integer k");
    assert!([16, 90].contains(&k));
    assert!(!run_core16 || k == 16);
    assert!(!run_cubic90 || k == 90);
    let source_variables: Vec<usize> = (0..18)
        .flat_map(|bit| (0..5).map(move |summand| bit + 18 * summand))
        .collect();
    assert_eq!(source_variables.len(), 90);
    let start = Instant::now();
    let kc = KoblitzCurve::known_n83_k0().expect("pinned n83 K0");
    assert_eq!(kc.label(), "icv1-f2m83-tm6151469093347-debefd74");
    let source = build_standard_subspace_factor_base(&kc, 18).expect("standard dimension-18 base");
    let projected = cofactor_project_factor_base(&kc, &source).expect("usable subgroup base");
    let usable_points = projected.points.len();
    let signed_columns = projected.signed_orbits.len();
    drop(projected);
    let table = TwoWordFieldStructure::new(83, &kc.curve.irreducible).expect("field table");
    let four = BigUint::from(4u32);
    let inv_four = (&kc.subgroup_order * BigUint::from(3u32) + BigUint::from(1u32)) / &four;
    let one = F2mElement::one(83);
    let zero = F2mElement::zero(83);
    let torsion = [
        BinaryPoint::Infinity,
        BinaryPoint::Affine {
            x: zero.clone(),
            y: one.clone(),
        },
        BinaryPoint::Affine {
            x: one.clone(),
            y: zero,
        },
        BinaryPoint::Affine {
            x: one.clone(),
            y: one,
        },
    ];
    for point in &torsion {
        assert!(kc.curve.is_on_curve(point));
        assert_eq!(kc.mul(point, &four), BinaryPoint::Infinity);
    }
    let public_target = BinaryPoint::Affine {
        x: F2mElement::from_hex("355fb5df7a905f16921eb", 83),
        y: F2mElement::from_hex("5900a390f42d290f1bbe", 83),
    };
    assert!(kc.curve.is_on_curve(&public_target));
    assert_eq!(
        kc.mul(&public_target, &kc.subgroup_order),
        BinaryPoint::Infinity
    );
    let setup_ns = start.elapsed().as_nanos();

    let (target, planted, offset) = match args[1].as_str() {
        "planted" => {
            let chosen = [0usize, 2, 4, 6, 8];
            let sum = chosen.iter().fold(BinaryPoint::Infinity, |acc, &i| {
                kc.add(&acc, &source.points[i])
            });
            assert_ne!(sum, BinaryPoint::Infinity);
            (sum, Some(chosen), None)
        }
        "ordinary" => {
            let i: usize = args
                .get(2)
                .expect("offset required")
                .parse()
                .expect("integer offset");
            assert!(i < 4);
            let subgroup_preimage = kc.mul(&public_target, &inv_four);
            assert_eq!(kc.mul(&subgroup_preimage, &four), public_target);
            (kc.add(&subgroup_preimage, &torsion[i]), None, Some(i))
        }
        _ => panic!("unknown mode"),
    };

    let build = Instant::now();
    let system = System512::build(
        &source.subspace_basis,
        x_of(&target),
        &kc.curve.b,
        5,
        &table,
    )
    .expect("339-variable system admitted");
    let build_ns = build.elapsed().as_nanos();
    assert_eq!(system.n_vars, 339);
    let original_equations = system.equations.len();
    assert_eq!(original_equations, 332);
    let mut planted_control = None;
    if let Some(chosen) = planted {
        let mut assignment = Mono512::default();
        let mut prefix = BinaryPoint::Infinity;
        for (i, &index) in chosen.iter().enumerate() {
            let point = &source.points[index];
            set_element(&mut assignment, i * 18, 18, x_of(point));
            prefix = kc.add(&prefix, point);
            if (1..=3).contains(&i) {
                set_element(&mut assignment, 90 + (i - 1) * 83, 83, x_of(&prefix));
            }
        }
        assert_eq!(prefix, target);
        assert!(system.all_vanish(&assignment));
        let mut products_checked = 0usize;
        for &variable in &source_variables[..k] {
            let multiplier = Poly512::var(variable);
            for equation in &system.equations {
                assert!(!equation.mul(&multiplier).eval(&assignment));
                products_checked += 1;
            }
        }
        assert_eq!(products_checked, original_equations * k);
        let projected = kc.mul(&target, &four);
        let subgroup_preimage = kc.mul(&projected, &inv_four);
        assert_eq!(kc.mul(&subgroup_preimage, &four), projected);
        let matching_offsets: Vec<_> = torsion
            .iter()
            .enumerate()
            .filter(|(_, t)| kc.add(&subgroup_preimage, t) == target)
            .map(|(i, _)| i)
            .collect();
        assert_eq!(matching_offsets.len(), 1);
        planted_control = Some(json!({
            "indices":chosen, "source_target":point_json(&target),
            "projected_target":point_json(&projected), "all_equations_zero":true,
            "group_sum_replayed":true, "torsion_offset":matching_offsets[0],
            "products_checked":products_checked
        }));
    }
    println!(
        "{}",
        json!({
            "phase":"constructed", "mode":args[1], "offset":offset,
            "curve_id":kc.label(), "source_points":source.points.len(),
            "usable_points":usable_points, "signed_columns":signed_columns,
            "target":point_json(&target), "n_vars":system.n_vars,
            "k":k, "source_variables":&source_variables[..k],
            "original_equations":original_equations,
            "original_term_occurrences":system.monomial_count(),
            "max_original_degree":system.max_degree(),
            "setup_ns":setup_ns, "build_ns":build_ns,
            "peak_rss_bytes":peak_rss_bytes(), "planted_control":planted_control,
            "claim_scope":"algebraic_feasibility_only"
        })
    );
    std::io::stdout()
        .flush()
        .expect("flush construction receipt");
    if offset.is_none() {
        return;
    }
    if peak_rss_bytes().is_some_and(|bytes| bytes > RSS_LIMIT_BYTES) {
        println!(
            "{}",
            json!({
                "phase":"private_certificate", "offset":offset, "k":k,
                "status":"rss_cap", "rss_limit_bytes":RSS_LIMIT_BYTES,
                "peak_rss_bytes":peak_rss_bytes(), "claim_scope":"algebraic_feasibility_only"
            })
        );
        return;
    }
    let begin = Instant::now();
    let certificate = system
        .private_degree_four(&source_variables[..k], CANDIDATES_PER_ROW)
        .expect("valid exact source-variable certificate");
    let certificate_ns = begin.elapsed().as_nanos();
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"f6-private-degree-four-v1");
    for (row, monomial) in &certificate.private_witnesses {
        hasher.update(&(*row as u64).to_le_bytes());
        for word in monomial.0 {
            hasher.update(&word.to_le_bytes());
        }
    }
    let witness_digest = hasher.finalize().to_hex().to_string();
    let unresolved_count = certificate.unresolved_rows.len();
    if run_core16 {
        assert_eq!(unresolved_count, 1_928, "frozen core size changed");
        assert_eq!(
            witness_digest, "ba12cc6bb138b135f514f6bb6198db21286101b4a55eeb6f42a5af12e4a49d74",
            "frozen private witness set changed"
        );
    }
    if run_cubic90 {
        assert_eq!(unresolved_count, 15_822, "frozen quartic core size changed");
        assert_eq!(
            witness_digest, "153342de74d4e40a2ea5f96ba35d24b7bdf2726700ed9b04714c9c89ca49dbdd",
            "frozen quartic witness set changed"
        );
    }
    println!(
        "{}",
        json!({
            "phase":"private_certificate", "offset":offset, "k":k,
            "status":"complete", "prolonged_rows":certificate.prolonged_rows,
            "rows_with_degree_four":certificate.rows_with_degree_four,
            "degree_four_occurrences":certificate.degree_four_occurrences,
            "candidate_columns":certificate.candidate_columns,
            "certified_rows":certificate.private_witnesses.len(),
            "unresolved_rows":&certificate.unresolved_rows,
            "witness_digest_blake3":witness_digest,
            "certificate_ns":certificate_ns, "peak_rss_bytes":peak_rss_bytes(),
            "claim_scope":"algebraic_feasibility_only"
        })
    );
    std::io::stdout()
        .flush()
        .expect("flush certificate receipt");
    let cubic_certificate = if run_cubic90 {
        let begin = Instant::now();
        let cubic = system
            .private_degree_three_after_quartic(
                &source_variables[..k],
                &certificate.unresolved_rows,
                if all_cubic {
                    usize::MAX
                } else {
                    CANDIDATES_PER_ROW
                },
            )
            .expect("ordered exact quartic-unresolved rows");
        let mut cubic_hasher = blake3::Hasher::new();
        cubic_hasher.update(b"f6-private-degree-three-after-four-v1");
        for (row, monomial) in &cubic.private_witnesses {
            cubic_hasher.update(&(*row as u64).to_le_bytes());
            for word in monomial.0 {
                cubic_hasher.update(&word.to_le_bytes());
            }
        }
        println!(
            "{}",
            json!({
                "phase":"private_cubic_certificate", "offset":offset, "k":k,
                "candidate_limit":if all_cubic {"all"} else {"32"},
                "status":"complete", "input_rows":cubic.input_rows,
                "rows_with_degree_three":cubic.rows_with_degree_three,
                "degree_three_occurrences":cubic.degree_three_occurrences,
                "candidate_columns":cubic.candidate_columns,
                "certified_rows":cubic.private_witnesses.len(),
                "unresolved_rows":&cubic.unresolved_rows,
                "witness_digest_blake3":cubic_hasher.finalize().to_hex().to_string(),
                "certificate_ns":begin.elapsed().as_nanos(),
                "peak_rss_bytes":peak_rss_bytes(),
                "claim_scope":"algebraic_feasibility_only"
            })
        );
        std::io::stdout().flush().expect("flush cubic certificate");
        Some(cubic)
    } else {
        None
    };
    let core_rows = cubic_certificate
        .as_ref()
        .map_or(&certificate.unresolved_rows, |cubic| &cubic.unresolved_rows);
    if (k != 90 && !run_core16) || core_rows.len() > 2_000 {
        return;
    }
    let core_build = Instant::now();
    let core = system
        .source_prolongation_core(&source_variables[..k], core_rows)
        .expect("ordered unresolved row IDs");
    let core_build_ns = core_build.elapsed().as_nanos();
    println!(
        "{}",
        json!({
            "phase":"core_constructed", "offset":offset, "k":k,
            "equations":core.equations.len(),
            "term_occurrences":core.monomial_count(),
            "core_build_ns":core_build_ns, "peak_rss_bytes":peak_rss_bytes(),
            "claim_scope":"algebraic_feasibility_only"
        })
    );
    std::io::stdout().flush().expect("flush core receipt");
    if peak_rss_bytes().is_some_and(|bytes| bytes > RSS_LIMIT_BYTES) {
        println!(
            "{}",
            json!({
                "phase":"core_reduction", "offset":offset, "k":k,
                "status":"rss_cap", "rss_limit_bytes":RSS_LIMIT_BYTES,
                "peak_rss_bytes":peak_rss_bytes(), "claim_scope":"algebraic_feasibility_only"
            })
        );
        return;
    }
    let reduce = Instant::now();
    let result = if wide_cap {
        core.root_reduce_with_column_cap(6_500_000)
    } else {
        core.root_reduce()
    };
    let (status, columns, rank, linear, linear_equations, source_only, contradiction, xor_ops) =
        match result {
            RootReduction::ColumnLimit { columns } => {
                ("column_limit", columns, None, None, None, None, None, None)
            }
            RootReduction::Reduced {
                columns,
                rank,
                linear,
                linear_equations,
                contradiction,
                xor_ops,
            } => {
                let encoded: Vec<_> = linear_equations
                    .into_iter()
                    .map(|equation| {
                        let source_bits = equation
                            .variables
                            .iter()
                            .filter(|&&v| v < core.summand_bits)
                            .count();
                        json!({
                            "variables":equation.variables,
                            "constant":equation.constant,
                            "source_only":source_bits == equation.variables.len()
                        })
                    })
                    .collect();
                let source_only = encoded
                    .iter()
                    .filter(|row| {
                        row["source_only"].as_bool() == Some(true)
                            && row["variables"]
                                .as_array()
                                .is_some_and(|variables| !variables.is_empty())
                    })
                    .count();
                (
                    "reduced",
                    columns,
                    Some(rank),
                    Some(linear),
                    Some(encoded),
                    Some(source_only),
                    Some(contradiction),
                    Some(xor_ops),
                )
            }
        };
    println!(
        "{}",
        json!({
            "phase":"core_reduction", "offset":offset, "k":k,
            "column_cap":if wide_cap {6_500_000} else {MAX_ROOT_COLS},
            "status":status, "columns":columns, "rank":rank,
            "linear":linear, "linear_equations":linear_equations,
            "source_only_nonconstant":source_only,
            "contradiction":contradiction, "xor_ops":xor_ops,
            "reduce_ns":reduce.elapsed().as_nanos(), "peak_rss_bytes":peak_rss_bytes(),
            "claim_scope":"algebraic_feasibility_only"
        })
    );
}
