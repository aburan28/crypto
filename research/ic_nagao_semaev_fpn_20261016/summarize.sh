#!/usr/bin/env bash
# Per (cell, arm) summary of a run's rows.jsonl: medians and ranges of the
# recorded counters, completion counts, failures by reason, and the
# correctness gate (equal relation sets, planted found). jq only; no new
# computation beyond order statistics and the registered ratio.
#   summarize.sh <run-dir>  → <run-dir>/summary.json
set -eu
RUN=${1:?run dir}
jq -s '
  def med: sort | if length == 0 then null else .[(length - 1) / 2 | floor] end;
  def rng: if length == 0 then null else [min, max] end;
  (map(select(.row == "target")) | map(. as $t | .arms[] | . + {n: $t.n, m: $t.m, kind: $t.kind, target: $t.target})) as $ok
  | (map(select(.row == "failure"))) as $bad
  | ([$ok[], $bad[]] | map({n, m}) | unique) as $cells
  | [ $cells[] as $c
      | ["sem", "rr"][] as $arm
      | ($ok | map(select(.n == $c.n and .m == $c.m and .arm == $arm))) as $r
      | ($r | map(select(.status != "undetermined"))) as $done
      | ($bad | map(select(.n == $c.n and .m == $c.m and .arm == $arm))) as $f
      | {
          n: $c.n, m: $c.m, arm: $arm,
          rows: ($r | length), completed: ($done | length),
          undetermined: ($r | map(select(.status == "undetermined")) | length),
          failures: ($f | group_by(.reason) | map({key: .[0].reason, value: length}) | from_entries),
          n_vars: ($r | first | .n_vars), n_eq: ($r | first | .n_eq), input_degree: ($r | first | .input_degree),
          solving_degree: ($done | map(.solving_degree) | rng),
          max_cols: ($done | map(.max_cols) | rng),
          delta: ($done | map(.delta) | rng),
          f4_muls_median: ($done | map(.f4_muls) | med), f4_muls_range: ($done | map(.f4_muls) | rng),
          lin_muls_median: ($done | map(.lin_muls) | med), lin_muls_range: ($done | map(.lin_muls) | rng),
          total_median: ($done | map(.f4_muls + .lin_muls) | med),
          total_range: ($done | map(.f4_muls + .lin_muls) | rng),
          f4_ms_median: ($done | map(.f4_ms) | med), lin_ms_median: ($done | map(.lin_ms) | med),
          planted_found: ($r | map(select(.planted_found != null)) | [map(select(.planted_found)) | length, length]),
          relations_total: ($done | map(.relations | length) | add)
        } ]
  | . as $s
  | {
      arms: $s,
      ratio_rr_over_sem: [ $cells[] as $c
        | ($s | map(select(.n == $c.n and .m == $c.m))) as $p
        | ($p | map(select(.arm == "sem")) | first) as $a
        | ($p | map(select(.arm == "rr")) | first) as $b
        | {n: $c.n, m: $c.m,
           complete: ($a.rows > 0 and $b.rows > 0 and $a.completed == $a.rows and $b.completed == $b.rows and ($a.failures | length) == 0 and ($b.failures | length) == 0),
           rho_total: (if $a.total_median and $b.total_median and $a.total_median > 0 then $b.total_median / $a.total_median else null end),
           rho_f4: (if $a.f4_muls_median and $b.f4_muls_median and $a.f4_muls_median > 0 then $b.f4_muls_median / $a.f4_muls_median else null end)} ],
      gate_relation_sets_equal: (
        [ ($ok | group_by([.n, .m, .kind, .target])[]
            | select(length == 2 and all(.status != "undetermined"))
            | (.[0].relations == .[1].relations)) ] | {checked: length, equal: map(select(.)) | length})
    }
' "$RUN/rows.jsonl" > "$RUN/summary.json"
echo "wrote $RUN/summary.json"
