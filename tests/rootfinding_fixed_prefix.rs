//! Frozen equal-probe coverage control; see the study's FIXED_PREFIX_PROTOCOL.md.
use crypto_lib::cryptanalysis::ic_boundary::{
    collect_and_solve_relations, prime_factor_base_strided, prime_instance_for, CompletionRule,
    CountedGroup, FieldOps, GroupOps, IncrementalGauss, PairTable, PrimePoint, RestartPool,
    RowStatus, TargetSource,
};
use serde_json::{json, Value};

const PREFIX: usize = 8_192;
const TRIAL_CAP: u64 = 20_000;
const MILESTONES: [usize; 5] = [512, 1_024, 2_048, 4_096, PREFIX];

// coverage-holdout ECBS1h1b585c89ebd1, measured round 1, ordinary B48.
// records.jsonl SHA-256: 8305acc015dc77437df966126bb915bdeae244a1b7c559bc52d7f3413e02601a
const TARGETS: [(u64, u64, u64); 8] = [
    (0xcba6, 0x36367, 13085358730139183541),
    (0x149f, 0x80ddd, 4888437509440305614),
    (0x3e54, 0x5e39e, 6804247407036939956),
    (0x406b4, 0x2c935, 10931747248848763017),
    (0x2f9fb, 0x76b6, 17647644974317840469),
    (0x13e09, 0x80b60, 16023983659542281631),
    (0x40f2, 0x5e73e, 11298333846333399661),
    (0x48082, 0x296dd, 12176717737306572033),
];

fn plus_mod(left: u64, right: u64, modulus: u64) -> u64 {
    ((left as u128 + right as u128) % modulus as u128) as u64
}

#[test]
fn fixed_public_prefix_gives_equal_probe_budgets_for_three_bases() {
    let inst = prime_instance_for(20, 59297);
    assert_eq!(inst.r, 562427);
    let curve = &inst.curve;
    let generator = inst.generator_point();
    let stream_base = prime_factor_base_strided(&inst, 48, 1);
    let mut targets_json = Vec::<Value>::new();

    for (target_index, (x, y, seed)) in TARGETS.into_iter().enumerate() {
        let target = PrimePoint::affine(x, y);
        assert!(curve.is_on_curve(target));
        let mut collector_matrix = IncrementalGauss::new(stream_base.columns + 1, inst.r);
        let mut stream: Vec<(PrimePoint, u64, u64, u64, FieldOps)> = Vec::new();
        let collected = collect_and_solve_relations(
            curve,
            generator,
            target,
            inst.r,
            inst.cofactor,
            &stream_base,
            seed,
            TRIAL_CAP,
            TargetSource::WalkPrimeLanes {
                lanes: 8,
                batched: true,
            },
            RestartPool::Lazy,
            &mut collector_matrix,
            CompletionRule::TargetPinned,
            |ops, _, point, a, b| {
                stream.push((point, a, b, ops.adds + ops.doubles, ops.field_ops));
                Vec::new()
            },
        );
        assert_eq!(collected.trials, TRIAL_CAP);
        assert!(
            stream.len() >= PREFIX,
            "target {target_index} had only {} surviving probes",
            stream.len()
        );
        stream.truncate(PREFIX);
        let mut arms = Vec::<Value>::new();

        for stride in [1, 7, 23] {
            let fb = prime_factor_base_strided(&inst, 48, stride);
            let table = PairTable::build_negation_folded_prime_batch(curve, &fb).unwrap();
            let mut matrix = IncrementalGauss::new(fb.columns + 1, inst.r);
            let mut probe_ops = GroupOps::default();
            let mut lookups = 0u64;
            let mut hits = 0u64;
            let mut milestones = Vec::<Value>::new();
            for (index, (point, a, b, walk_gae, walk_field_ops)) in stream.iter().enumerate() {
                if let Some((first, second)) =
                    table.probe(curve, &fb, &mut probe_ops, &mut lookups, *point)
                {
                    hits += 1;
                    let mut row = vec![0u64; fb.columns + 1];
                    for summand in [first, second] {
                        let col = fb.col_of[summand];
                        row[col] = plus_mod(row[col], fb.coef_of[summand], inst.r);
                    }
                    row[fb.columns] = (inst.r - b % inst.r) % inst.r;
                    assert_ne!(matrix.add_row(row, a % inst.r), RowStatus::Inconsistent);
                }
                if MILESTONES.contains(&(index + 1)) {
                    let cold_gae_lower_bound = *walk_gae as f64
                        + fb.cost.group_ops.gae()
                        + table.build_ops.gae()
                        + probe_ops.gae();
                    milestones.push(json!({
                        "probes": index + 1,
                        "hits": hits,
                        "independent_rank": matrix.rank(),
                        "dependent_rows": matrix.dependent,
                        "row_ops": matrix.row_ops,
                        "cold_gae_lower_bound": cold_gae_lower_bound,
                        "group_field_invs_only": walk_field_ops.invs + table.build_ops.field_ops.invs + probe_ops.field_ops.invs,
                    }));
                }
            }
            assert_eq!(lookups as usize, PREFIX);
            let structural_bytes_floor = fb.points.len() * std::mem::size_of::<PrimePoint>()
                + table.entries as usize * (std::mem::size_of::<u64>() + 2 * std::mem::size_of::<u32>());
            arms.push(json!({
                "stride": stride,
                "unique_pair_keys": table.entries,
                "uniform_hit_probability": 2.0 * table.entries as f64 / (inst.r - 1) as f64,
                "structural_bytes_floor": structural_bytes_floor,
                "milestones": milestones,
            }));
        }
        targets_json.push(json!({
            "target_index": target_index,
            "target": [format!("0x{x:x}"), format!("0x{y:x}")],
            "algorithm_seed": seed,
            "attempted_points": TRIAL_CAP,
            "surviving_points": collected.relations.get("target_guard_probes") - collected.relations.get("repeated_targets_skipped"),
            "fixed_prefix": PREFIX,
            "arms": arms,
        }));
    }
    println!(
        "FIXED_PREFIX {}",
        json!({
            "schema": "rootfinding.fixed_prefix/v1",
            "curve": "icv1-fp20-t727-cd198a38",
            "subgroup_r": inst.r,
            "prefix_surviving_probes": PREFIX,
            "trial_cap": TRIAL_CAP,
            "lane_width": 8,
            "target_source": "coverage-holdout measured round 1",
            "targets": targets_json,
        })
    );
}
