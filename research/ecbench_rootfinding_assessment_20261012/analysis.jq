[ .[] | select(.warmup == false)]
| group_by(.arm)
| map(. as $rows | {
    session_id: $rows[0].session_id,
    arm: $rows[0].arm,
    curve: $rows[0].workload.curve.slug,
    n: ($rows | length),
    verified: ([$rows[].outcome.status] | map(select(. == "verified")) | length),
    mean_s: ([$rows[].cost.s] | add / length),
    mean_gae: ([$rows[].cost.total_gae] | add / length),
    mean_max_rss_kib: ([$rows[].time.max_rss_kib] | add / length),
    max_max_rss_kib: ([$rows[].time.max_rss_kib] | max),
    lower_bound: ([$rows[].cost.lower_bound] | any),
    unpriced: ([$rows[].cost.unpriced[]] | unique),
    fb_id: $rows[0].factor_base.fb_id,
    signed_points: $rows[0].factor_base.signed_points,
    columns: $rows[0].factor_base.columns,
    mean_trials: ([$rows[].counters.targets_tried] | map(select(. != null)) | if length == 0 then null else add / length end),
    phases: ([$rows[].phases[].name] | unique | map(. as $p | {
      name: $p,
      mean_gae: ([$rows[].phases[] | select(.name == $p) | .gae] | add / length)
    }))
  })
