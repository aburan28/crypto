fn logical_json(l: &Logical) -> String {
    format!("{{\"screen_second_checks\":{},\"screen_second_rejected\":{},\"screen_calls\":{},\"screen_points\":{},\"screen_batches\":{},\"screen_passes\":{},\"screen_full_checks\":{},\"screen_full_rejected\":{},\"screen_fallback_points\":{},\"screen_fallback_batches\":{},\"fiber_filter_rounds\":{},\"fiber_selected\":{},\"fiber_dimension\":{},\"fiber_selection_states\":{},\"fiber_prefixes\":{},\"fiber_batches\":{},\"fiber_queries\":{},\"fiber_zero_rejected\":{},\"fiber_linear_rejected\":{},\"fiber_extensions_rejected\":{},\"fiber_rank_sum\":{},\"projection_rows\":{},\"projection_xors\":{},\"projection_certified\":{},\"nodes\":{},\"kernel_calls\":{},\"decisions\":{},\"forced\":{},\"affine_eliminated\":{},\"derived_rows\":{},\"enumeration_points\":{},\"enumeration_batches\":{},\"enumeration_leaves\":{},\"specialized_terms\":{},\"source_rows\":{},\"source_columns\":{},\"max_depth\":{}}}",l.screen_second_checks,l.screen_second_rejected,l.screen_calls,l.screen_points,l.screen_batches,l.screen_passes,l.screen_full_checks,l.screen_full_rejected,l.screen_fallback_points,l.screen_fallback_batches,l.fiber_filter_rounds,l.fiber_selected,l.fiber_dimension,l.fiber_selection_states,l.fiber_prefixes,l.fiber_batches,l.fiber_queries,l.fiber_zero_rejected,l.fiber_linear_rejected,l.fiber_extensions_rejected,l.fiber_rank_sum,l.projection_rows,l.projection_xors,l.projection_certified,l.nodes,l.kernel_calls,l.decisions,l.forced,l.affine_eliminated,l.derived_rows,l.enumeration_points,l.enumeration_batches,l.enumeration_leaves,l.specialized_terms,l.source_rows,l.source_columns,l.max_depth)
}
const QUIET_ARMS: &[&str] = &[
    "packed_untraced",
    "gray_quiet",
    "wide64_quiet",
    "wide64_unrolled",
    "word16_unrolled",
    "word_dispatch",
    "leaf16_quiet",
    "leaf16_word_unrolled",
    "byte_scalar",
    "byte_simd",
    "byte_planes",
    "byte_unrolled",
    "byte_single_quiet",
    "leaf16_byte_scalar",
    "leaf16_byte_simd",
    "leaf16_byte_planes",
    "leaf16_byte_unrolled",
    "leaf16_single_quiet",
];
const DIAGNOSTIC_ARMS: [&str; 6] = [
    "construction16_rebound",
    "construction16_profile",
    "construction64_rebound",
    "construction64_profile",
    "construction_dispatch_rebound",
    "construction_dispatch_profile",
];
pub fn diagnostic_main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 7);
    let n: u8 = args[1].parse().unwrap();
    let seed: u64 = args[2].parse().unwrap();
    let family = &args[3];
    let repetitions: usize = args[4].parse().unwrap();
    let limit: u64 = args[5].parse().unwrap();
    let order_seed: u64 = args[6].parse().unwrap();
    assert!([12, 16, 20, 24].contains(&n) && repetitions == 8 && limit == 200000);
    let (system, witness) = fixture(n, seed, family);
    let reference = solve(&system, n, "search", limit);
    let original16 = solve(&system, n, "word16_unrolled", limit);
    let original64 = solve(&system, n, "wide64_unrolled", limit);
    let mut clocks = Vec::with_capacity(4096);
    for _ in 0..4096 {
        let a = Instant::now();
        let b = Instant::now();
        clocks.push(b.duration_since(a).as_nanos());
    }
    clocks.sort_unstable();
    let witness = witness.map_or("null".into(), |v| v.to_string());
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"family\":\"{family}\",\"polys\":{:?},\"planted_witness\":{witness},\"search_reference\":\"{}\",\"clock_samples\":4096,\"clock_median_ns\":{},\"clock_p99_ns\":{},\"clock_max_ns\":{}}}",system,status(&reference.outcome),clocks[2048],clocks[4055],clocks[4095]);
    let names: Vec<_> = FINAL_ARMS
        .iter()
        .chain(DIAGNOSTIC_ARMS.iter())
        .copied()
        .collect();
    let mut stable: Vec<Option<(Outcome, Logical, u64)>> =
        std::iter::repeat_with(|| None).take(names.len()).collect();
    for rep in 0..repetitions {
        for (order, index) in paired_order(names.len(), order_seed, rep)
            .into_iter()
            .enumerate()
        {
            let arm = names[index];
            let begin = Instant::now();
            let (got, phases) = if arm.starts_with("construction") {
                let width = if arm.contains("16_") || (arm.contains("dispatch") && n <= 20) {
                    16
                } else {
                    64
                };
                if arm.ends_with("profile") {
                    construction_solve::<true>(black_box(&system), n, limit, width)
                } else {
                    construction_solve::<false>(black_box(&system), n, limit, width)
                }
            } else {
                (solve(black_box(&system), n, arm, limit), None)
            };
            let solve_ns = begin.elapsed().as_nanos();
            let validation = Instant::now();
            match got.outcome {
                Outcome::Sat(point) => {
                    assert!(point < (1u64 << n));
                    assert!(satisfies(&system, point));
                }
                Outcome::Unsat => {
                    assert_eq!(reference.outcome, Outcome::Unsat);
                }
                Outcome::Unknown(_) => (),
            }
            if !matches!(got.outcome, Outcome::Unknown(_))
                && !matches!(reference.outcome, Outcome::Unknown(_))
            {
                assert_eq!(status(&got.outcome), status(&reference.outcome));
            }
            if arm.starts_with("construction") {
                let old = if arm.contains("16_") || (arm.contains("dispatch") && n <= 20) {
                    &original16
                } else {
                    &original64
                };
                assert_eq!(got.outcome, old.outcome);
                assert_eq!(got.logical, old.logical);
                assert_eq!(got.trace, old.trace);
            }
            let signature = (got.outcome.clone(), got.logical.clone(), got.trace);
            if let Some(prior) = &stable[index] {
                assert_eq!(prior, &signature);
            } else {
                stable[index] = Some(signature);
            }
            let verified = matches!(got.outcome, Outcome::Sat(_))
                || (got.outcome == Outcome::Unsat && reference.outcome == Outcome::Unsat);
            let validation_ns = validation.elapsed().as_nanos();
            let total_ns = begin.elapsed().as_nanos();
            let model = match got.outcome {
                Outcome::Sat(v) => v.to_string(),
                _ => "null".into(),
            };
            let reason = match got.outcome {
                Outcome::Unknown(r) => format!("\"{r}\""),
                _ => "null".into(),
            };
            let trace = if QUIET_ARMS.contains(&arm) || arm.starts_with("construction") {
                "null".into()
            } else {
                got.trace.to_string()
            };
            let phases=phases.map_or("null".into(),|p| {
                assert!(p.encoding_ns+p.construction_ns+p.scan_ns+p.release_ns<=solve_ns);
                format!("{{\"encoding_ns\":{},\"construction_ns\":{},\"scan_ns\":{},\"release_ns\":{}}}",p.encoding_ns,p.construction_ns,p.scan_ns,p.release_ns)
            });
            println!("{{\"type\":\"sample\",\"rep\":{rep},\"order\":{order},\"variant\":\"{arm}\",\"outcome\":\"{}\",\"model\":{model},\"reason\":{reason},\"verified\":{verified},\"solve_ns\":{solve_ns},\"validation_ns\":{validation_ns},\"total_ns\":{total_ns},\"trace\":{trace},\"logical\":{},\"phases\":{phases}}}",status(&got.outcome),logical_json(&got.logical));
        }
    }
}
