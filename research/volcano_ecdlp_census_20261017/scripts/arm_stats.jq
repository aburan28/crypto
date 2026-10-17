# input: records.jsonl slurped; output: per-arm stats over measured (non-warmup) records.
def mean: add/length;
def sd: (mean) as $m | (map(.-$m|.*.)|add/((length-1)|if .<1 then 1 else . end))|sqrt;
[ .[] | select(.warmup==false) ] | group_by(.arm) | map(
  . as $rs |
  ($rs | group_by(.workload.workload_id) | map({s:(map(.cost.s)|mean), w:(map(.online.wall_ns // 0)|mean)})) as $wl |
  { arm: $rs[0].arm, solves: ($rs|length), workloads: ($wl|length),
    verified: ([$rs[]|select(.outcome.status=="verified")]|length),
    verify_rate: (([$rs[]|select(.outcome.status=="verified")]|length)/($rs|length)),
    mean_S: ($rs|map(.cost.s)|mean),
    cv_S_over_workloads: (($wl|map(.s)) as $x | ($x|mean) as $m | if $m==0 then null else (($x|sd)/$m) end),
    cv_S_over_records: (($rs|map(.cost.s)) as $x | (($x|sd)/($x|mean))),
    mean_wall_s: (($rs|map(.online.wall_ns // 0)|mean)/1e9),
    levels: ($rs|map(.isolation.level // "?")|group_by(.)|map({(.[0]):length})|add) })
