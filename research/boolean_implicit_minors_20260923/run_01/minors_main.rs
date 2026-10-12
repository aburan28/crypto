const MINOR_ARMS: [&str; 3] = ["numeric_minor4", "symbolic_minor4", "affine_filter"];
fn minor_output_json(output: &MinorOutput) -> String {
    format!(
        "{{\"accepted\":{:?},\"minors\":{:?}}}",
        output.accepted, output.minors
    )
}
pub fn minors_main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 7);
    let n: u8 = args[1].parse().unwrap();
    let seed: u64 = args[2].parse().unwrap();
    let family = &args[3];
    let repetitions: usize = args[4].parse().unwrap();
    let limit: u64 = args[5].parse().unwrap();
    let order_seed: u64 = args[6].parse().unwrap();
    assert!([12, 16].contains(&n) && repetitions == 8 && limit == 200000);
    let (system, witness) = fixture(n, seed, family);
    let oracle = minor_oracle(&system, n);
    let witness = witness.map_or("null".into(), |x| x.to_string());
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"family\":\"{family}\",\"polys\":{:?},\"planted_witness\":{witness},\"original_solutions\":{},\"rows\":{:?},\"oracle_minor\":{},\"oracle_affine\":{}}}",system,oracle.original_solutions,oracle.rows,minor_output_json(&oracle.minor),minor_output_json(&oracle.affine));
    let names: Vec<_> = FINAL_ARMS
        .iter()
        .chain(PROJECTED_ARMS.iter())
        .chain(MINOR_ARMS.iter())
        .copied()
        .collect();
    let mut stable: Vec<Option<String>> = vec![None; names.len()];
    for rep in 0..repetitions {
        for (order, index) in paired_order(names.len(), order_seed, rep)
            .into_iter()
            .enumerate()
        {
            let arm = names[index];
            let begin = Instant::now();
            let (solved, filtered, projected) = if MINOR_ARMS.contains(&arm) {
                (
                    None,
                    Some(minor_build(black_box(&system), n, arm, 50000000)),
                    None,
                )
            } else if PROJECTED_ARMS.contains(&arm) {
                let (got, work) = match arm {
                    "projected4" => solve_projected::<4>(black_box(&system), n, 1 << 24, 1 << 24),
                    "projected5" => solve_projected::<5>(black_box(&system), n, 1 << 24, 1 << 24),
                    "projected6" => solve_projected::<6>(black_box(&system), n, 1 << 24, 1 << 24),
                    _ => unreachable!(),
                };
                (Some(got), None, Some(work))
            } else {
                (Some(solve(black_box(&system), n, arm, limit)), None, None)
            };
            let compute_ns = begin.elapsed().as_nanos();
            let validation = Instant::now();
            let verified = if let Some(got) = &solved {
                match got.outcome {
                    Outcome::Sat(point) => {
                        assert!(
                            point < 1u64 << n
                                && satisfies(&system, point)
                                && oracle.original_solutions > 0
                        );
                        true
                    }
                    Outcome::Unsat => {
                        assert_eq!(oracle.original_solutions, 0);
                        true
                    }
                    Outcome::Unknown(_) => false,
                }
            } else {
                let got = filtered.as_ref().unwrap();
                assert!(got.setup_ns + got.construction_ns + got.evaluation_ns <= compute_ns);
                if let Some(output) = &got.output {
                    assert_eq!(
                        output,
                        if arm == "affine_filter" {
                            &oracle.affine
                        } else {
                            &oracle.minor
                        }
                    );
                    assert!(got.reason.is_none());
                    true
                } else {
                    assert!(got.reason.is_some());
                    false
                }
            };
            let validation_ns = validation.elapsed().as_nanos();
            let total_ns = begin.elapsed().as_nanos();
            // Diagnostic serialization and repeat-signature checks are outside arm clocks.
            let (kind, outcome, model, reason, trace, logical, work, filter, phases) =
                if let Some(got) = solved {
                    let model = match got.outcome {
                        Outcome::Sat(x) => x.to_string(),
                        _ => "null".into(),
                    };
                    let reason = match got.outcome {
                        Outcome::Unknown(r) => format!("\"{r}\""),
                        _ => "null".into(),
                    };
                    let trace = if QUIET_ARMS.contains(&arm) || PROJECTED_ARMS.contains(&arm) {
                        "null".into()
                    } else {
                        got.trace.to_string()
                    };
                    (
                        "solver",
                        status(&got.outcome),
                        model,
                        reason,
                        trace,
                        logical_json(&got.logical),
                        "null".into(),
                        "null".into(),
                        "null".into(),
                    )
                } else {
                    let got = filtered.unwrap();
                    let outcome = if got.output.is_some() {
                        "COMPLETE"
                    } else {
                        "CAPPED"
                    };
                    let reason = got.reason.map_or("null".into(), |r| format!("\"{r}\""));
                    let filter = got.output.as_ref().map_or("null".into(), minor_output_json);
                    let phases = format!(
                        "{{\"setup_ns\":{},\"construction_ns\":{},\"evaluation_ns\":{}}}",
                        got.setup_ns, got.construction_ns, got.evaluation_ns
                    );
                    (
                        "filter",
                        outcome,
                        "null".into(),
                        reason,
                        "null".into(),
                        "null".into(),
                        got.work.json(),
                        filter,
                        phases,
                    )
                };
            let projected = projected.map_or("null".into(), |w| w.json());
            let signature =
                format!("{outcome}/{model}/{reason}/{trace}/{logical}/{work}/{filter}/{projected}");
            if let Some(previous) = &stable[index] {
                assert_eq!(previous, &signature);
            } else {
                stable[index] = Some(signature);
            }
            println!("{{\"type\":\"sample\",\"rep\":{rep},\"order\":{order},\"variant\":\"{arm}\",\"kind\":\"{kind}\",\"outcome\":\"{outcome}\",\"model\":{model},\"reason\":{reason},\"verified\":{verified},\"compute_ns\":{compute_ns},\"validation_ns\":{validation_ns},\"total_ns\":{total_ns},\"trace\":{trace},\"logical\":{logical},\"projected_work\":{projected},\"work\":{work},\"filter\":{filter},\"phases\":{phases}}}");
        }
    }
}
