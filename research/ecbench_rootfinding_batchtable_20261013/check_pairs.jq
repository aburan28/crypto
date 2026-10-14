# Mechanical equality check for the two implementations of the same
# two-summand search. Read the sealed records with jq -s -f check_pairs.jq.
[.[] | select(.arm == "pair-scalar" or .arm == "pair-batch")]
| group_by([.workload.workload_id, .round, .warmup])
| map(
    sort_by(.arm)
    | if length != 2 then
        {pair_count: length, mismatch: "missing arm"}
      else
        .[0] as $batch | .[1] as $scalar
        | {
            workload: $batch.workload.workload_id,
            round: $batch.round,
            warmup: $batch.warmup,
            same_seed: ($batch.algorithm_seed == $scalar.algorithm_seed),
            same_status: ($batch.outcome.status == $scalar.outcome.status),
            same_recovered: ($batch.outcome.recovered == $scalar.outcome.recovered),
            same_target_check: ($batch.outcome.matches_target == $scalar.outcome.matches_target),
            same_fb: ($batch.factor_base.fb_sha256 == $scalar.factor_base.fb_sha256),
            same_trials: ($batch.counters.targets_tried == $scalar.counters.targets_tried),
            same_relation_count: ($batch.counters.relations_found == $scalar.counters.relations_found),
            same_rank: ($batch.counters.matrix_rank == $scalar.counters.matrix_rank),
            same_total_gae: ($batch.cost.total_gae == $scalar.cost.total_gae),
            same_counted_phases: (
              [$batch.phases[] | {name, adds, doubles, scalar_mults, gae}]
              == [$scalar.phases[] | {name, adds, doubles, scalar_mults, gae}]
            )
          }
      end
  )
| {
    schema: "rootfinding.batchtable.pair-check/v1",
    pairs: length,
    passed: (map(select([.same_seed, .same_status, .same_recovered,
      .same_target_check, .same_fb, .same_trials, .same_relation_count,
      .same_rank, .same_total_gae, .same_counted_phases] | all)) | length),
    mismatches: [ .[] | select(( [.same_seed, .same_status,
      .same_recovered, .same_target_check, .same_fb, .same_trials,
      .same_relation_count, .same_rank, .same_total_gae,
      .same_counted_phases] | all) | not) ]
  }
