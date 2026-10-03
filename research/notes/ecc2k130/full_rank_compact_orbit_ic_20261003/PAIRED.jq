def phase($name): (.phases | map(select(.name == $name)) | first // {});

def arm_record:
  {
    arm,
    algorithm_seed,
    scalar: .outcome.recovered,
    verified: (.outcome.status == "verified" and .outcome.matches_target),
    gae: .cost.total_gae,
    s: .cost.s,
    lower_bound: .cost.lower_bound,
    factor_base_gae: phase("factor_base").gae,
    oracle_setup_gae: phase("oracle_setup").gae,
    relations_gae: phase("relations").gae,
    linalg_gae: phase("linear_algebra").gae,
    verify_gae: phase("verify").gae,
    trials: .counters.targets_tried,
    rows: .counters.matrix_rows,
    rank: .counters.matrix_rank,
    dependent: .counters.matrix_dependent,
    first_pin_rank: phase("linear_algebra").native.first_target_pinned_rank,
    first_pin_trial: phase("linear_algebra").native.first_target_pinned_trial,
    first_pin_relation: phase("linear_algebra").native.first_target_pinned_relation,
    full_rank_reached: phase("linear_algebra").native.full_rank_reached,
    factor_base_id: .factor_base.fb_id,
    binary_sha256
  };

[.[] | select(.warmup == false) | {
  workload: .workload.workload_id,
  target_index: .workload.target_index,
  round,
  arm
} + arm_record]
| group_by([.workload, .round])
| map(
    . as $g
    | ($g | map(select(.arm == "rho-strong")) | first) as $rho
    | ($g | map(select(.arm == "ic-early-pin")) | first) as $early
    | ($g | map(select(.arm == "ic-full-rank")) | first) as $full
    | {
        workload: $g[0].workload,
        target_index: $g[0].target_index,
        round: $g[0].round,
        rho: $rho,
        early: $early,
        full: $full,
        full_over_rho: ($full.gae / $rho.gae),
        full_over_early: ($full.gae / $early.gae),
        full_online_over_rho: (
          ($full.relations_gae + $full.linalg_gae + $full.verify_gae) / $rho.gae
        ),
        checkpoint_agrees: (
          $full.algorithm_seed == $early.algorithm_seed
          and $full.scalar == $early.scalar
          and $full.factor_base_id == $early.factor_base_id
          and $full.factor_base_gae == $early.factor_base_gae
          and $full.oracle_setup_gae == $early.oracle_setup_gae
          and $full.first_pin_rank == $early.rank
          and $full.first_pin_trial == $early.trials
          and $full.first_pin_relation == $early.rows
        )
      }
  )
