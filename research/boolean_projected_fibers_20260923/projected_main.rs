const PROJECTED_ARMS: [&str; 3] = ["projected4", "projected5", "projected6"];
pub fn projected_main() {
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
    let witness = witness.map_or("null".into(), |v| v.to_string());
    println!("{{\"type\":\"fixture\",\"n\":{n},\"seed\":{seed},\"family\":\"{family}\",\"polys\":{:?},\"planted_witness\":{witness},\"search_reference\":\"{}\"}}",system,status(&reference.outcome));
    let names: Vec<_> = FINAL_ARMS
        .iter()
        .chain(PROJECTED_ARMS.iter())
        .copied()
        .collect();
    let mut stable: Vec<Option<(Outcome, Logical, u64, Option<ProjectedWork>)>> =
        vec![None; names.len()];
    for rep in 0..repetitions {
        for (order, index) in paired_order(names.len(), order_seed, rep)
            .into_iter()
            .enumerate()
        {
            let arm = names[index];
            let begin = Instant::now();
            let (got, work) = if arm.starts_with("projected") {
                let (got, work) = match arm {
                    "projected4" => solve_projected::<4>(black_box(&system), n, 1 << 24, 1 << 24),
                    "projected5" => solve_projected::<5>(black_box(&system), n, 1 << 24, 1 << 24),
                    "projected6" => solve_projected::<6>(black_box(&system), n, 1 << 24, 1 << 24),
                    _ => unreachable!(),
                };
                (got, Some(work))
            } else {
                (solve(black_box(&system), n, arm, limit), None)
            };
            let solve_ns = begin.elapsed().as_nanos();
            let validation = Instant::now();
            match got.outcome {
                Outcome::Sat(point) => {
                    assert!(point < 1u64 << n && satisfies(&system, point));
                }
                Outcome::Unsat => assert_eq!(reference.outcome, Outcome::Unsat),
                Outcome::Unknown(_) => (),
            }
            if !matches!(got.outcome, Outcome::Unknown(_))
                && !matches!(reference.outcome, Outcome::Unknown(_))
            {
                assert_eq!(status(&got.outcome), status(&reference.outcome));
            }
            let signature = (
                got.outcome.clone(),
                got.logical.clone(),
                got.trace,
                work.clone(),
            );
            if let Some(previous) = &stable[index] {
                assert_eq!(&signature, previous);
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
            let trace = if QUIET_ARMS.contains(&arm) || work.is_some() {
                "null".into()
            } else {
                got.trace.to_string()
            };
            let work = work.map_or("null".into(), |w| w.json());
            println!("{{\"type\":\"sample\",\"rep\":{rep},\"order\":{order},\"variant\":\"{arm}\",\"outcome\":\"{}\",\"model\":{model},\"reason\":{reason},\"verified\":{verified},\"solve_ns\":{solve_ns},\"validation_ns\":{validation_ns},\"total_ns\":{total_ns},\"trace\":{trace},\"logical\":{},\"projected\":{work}}}",status(&got.outcome),logical_json(&got.logical));
        }
    }
}
