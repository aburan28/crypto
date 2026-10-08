#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo"
out=research/f6_ic_geometric_closure_20261003/small_cold
base=research/f6_ic_geometric_closure_20261003/three_way
binary=${IC_WORKER_BINARY:-/private/tmp/f6-ic-target/release/examples/ic_tournament_worker}
digest() { shasum -a 256 "$1" | awk '{print $1}'; }
test ! -e "$out/freeze.txt"
test "$(digest "$binary")" = "$(sed -n 's/^binary_sha256=//p' "$base/freeze.txt")"
test "$(jq -r .fixture.target_scalar_constructed "$out/fixture.json")" = false
test "$(jq -r .fixture.target_seeds[0] "$out/fixture.json")" = 20261004041
test "$(jq -r .status "$out/inventory.json")" = inventory
test "$(jq -r .effective_columns "$out/usable-inventory.json")" = "$(jq -r .columns "$out/inventory.json")"
test "$(jq -r '.usable_points|length' "$out/usable-inventory.json")" = 14
test "$(jq -r '.geometric_points|length' "$out/usable-inventory.json")" = 21
geometric_from_worker=$(jq -c '[.factor_base[]|map(tonumber)]|sort_by(.[0],.[1])' "$out/inventory.json")
geometric_from_helper=$(jq -c '[.geometric_points[]|map(tonumber)]|sort_by(.[0],.[1])' "$out/usable-inventory.json")
test "$geometric_from_worker" = "$geometric_from_helper"
geometric_sha=$(jq -jcS .geometric_points "$out/usable-inventory.json" | shasum -a 256 | awk '{print $1}')
usable_sha=$(jq -jcS .usable_points "$out/usable-inventory.json" | shasum -a 256 | awk '{print $1}')

jq -cnS --slurpfile fixture "$out/fixture.json" \
  '{field:{p:2,n:9,basis:"polynomial",modulus:515,
     element_encoding:"unsigned-integer-polynomial-bits"},
    curve:{model:"y^2+x*y=x^3+a2*x^2+a6",tag:"kb1",
      coefficients:{a1:1,a2:1,a3:0,a4:0,a6:1},
      order:($fixture[0].fixture.group_order|tonumber),
      trace:(513-($fixture[0].fixture.group_order|tonumber)),
      r:($fixture[0].fixture.subgroup_order|tonumber),
      cofactor:($fixture[0].fixture.cofactor|tonumber),
      generator:($fixture[0].fixture.generator|map(tonumber)),
      target_group:"prime-order-subgroup-generated-by-G"}}' \
  > "$out/curve-hash-input.json"
curve_sha=$(jq -jcS . "$out/curve-hash-input.json" | shasum -a 256 | awk '{print $1}')
curve_id="EC1N9Ckb1h$(printf '%.12s' "$curve_sha")"
jq -cS --arg id "$curve_id" '.curve.curve_id=$id' "$out/curve-hash-input.json" > "$out/curve-record.json"

for arm in inherited_f4 f5 f6_ic; do
  case "$arm" in
    inherited_f4) stage=f4; split=highest-free ;;
    f5) stage=f5; split=lowest-free ;;
    f6_ic) stage=f6; split=highest-free ;;
  esac
  jq -cS --slurpfile curve "$out/curve-record.json" \
    --slurpfile inv "$out/usable-inventory.json" \
    --arg geometric_sha "$geometric_sha" --arg usable_sha "$usable_sha" \
    --arg split "$split" --arg arm "$arm" --arg binary_sha "$(digest "$binary")" \
    '.field=$curve[0].field | .curve=$curve[0].curve |
      .factor_base={
        construction:{recipe:{kind:"standard_subspace",dimension:4},
          ordering:"constructor-abscissa-order-then-half-trace-root-and-negative",
          subgroup_map:"cofactor-multiplication-remove-identity-deduplicate"},
        nominal_bound:4,
        inventory:{geometric_point_count:($inv[0].geometric_points|length),
          geometric_set_sha256:$geometric_sha,
          geometric_order_sha256:null,
          identity_images:$inv[0].identity_images,
          duplicate_nonidentity_images:$inv[0].duplicate_nonidentity_images,
          usable_point_count:($inv[0].usable_points|length),
          usable_set_sha256:$usable_sha,
          effective_columns:$inv[0].effective_columns,
          subgroup_map:"cofactor-multiplication",quotient:"sign-and-Frobenius"}} |
      .point_decomposition.solver=$arm |
      .point_decomposition.limits.split_rule=$split |
      .relation_collection={collector:"sample",query_distribution:"target-independent",
        query_rule:"trial-keyed-sample",filtering:"group-verified-witness-only",
        verification:"group-readd-every-relation",
        duplicates:"scalar-and-sorted-base-indices",
        dependencies:"retain-unique-even-if-dependent",
        stop_rule:"first-verified-full-column-log-table-or-4096-ordinary-trials",
        source_sha256:.relation_collection.source_sha256} |
      .relation_linear_algebra.modulus=37 |
      .target_descent.stop_rule.max_trials=4096 |
      .implementation.components |= map(select(.role!="prepared-log-state")) |
      .implementation.flags={binary_sha256:$binary_sha,rayon_threads:1,
        exclusive_phases:true,preparation_mode:"native-cold-relation-collection",
        reusable_template:"process-local-source-defaults",f6_ic:($arm=="f6_ic")}' \
    "$base/$arm-record.json" > "$out/$arm-record.json"
  sha=$(jq -jcS . "$out/$arm-record.json" | shasum -a 256 | awk '{print $1}')
  id="IC1N9Ckb1fb14PDP3${stage}RCsampleLAgaussTDpdpISO0h$(printf '%.12s' "$sha")"
  jq -cnS --slurpfile record "$out/$arm-record.json" --arg sha "$sha" --arg id "$id" \
    '{record:$record[0],record_sha256:$sha,candidate_id:$id}' > "$out/$arm-candidate.json"
done

jq -cnS --slurpfile fixture "$out/fixture.json" --arg curve_id "$curve_id" \
  '{algorithm_seed:20261004042,cache_policy:"cold",
    curve_id:$curve_id,input_law:"blake3-ic-workflow-public-target-v1; fixture-construction-excluded",
    point_was_previously_supplied:false,
    question:"full-cold-n9-f4-f5-f6-ic-one-public-target-diagnostic",
    resource_envelope:{cpu_workers:1,host_class:"physical-macos-arm64",
      memory_limit_bytes:null,target_count:1,total_wall_limit_seconds:180},
    target_count:1,target_seeds:[20261004041],
    targets:[$fixture[0].fixture.targets[0]|map(tonumber)]}' \
  > "$out/workload-record.json"
workload_sha=$(jq -jcS . "$out/workload-record.json" | shasum -a 256 | awk '{print $1}')
jq -cnS --slurpfile record "$out/workload-record.json" --arg sha "$workload_sha" \
  '{record:$record[0],record_sha256:$sha,workload_id:$sha[0:12]}' > "$out/workload.json"

printf '%s\n' "source_commit=$(sed -n 's/^source_commit=//p' "$base/freeze.txt")" \
  "binary_sha256=$(digest "$binary")" \
  "inventory_helper_sha256=$(digest examples/ic_f6_small_inventory.rs)" \
  "fixture_sha256=$(digest "$out/fixture.json")" \
  "inventory_sha256=$(digest "$out/inventory.json")" \
  "usable_inventory_sha256=$(digest "$out/usable-inventory.json")" \
  "input_inherited_f4_sha256=$(digest "$out/input-inherited_f4.json")" \
  "input_f6_ic_sha256=$(digest "$out/input-f6_ic.json")" \
  "input_f5_sha256=$(digest "$out/input-f5.json")" \
  "curve_id=$curve_id" "workload_id=$(jq -r .workload_id "$out/workload.json")" \
  "inherited_f4=$(jq -r .candidate_id "$out/inherited_f4-candidate.json")" \
  "f6_ic=$(jq -r .candidate_id "$out/f6_ic-candidate.json")" \
  "f5=$(jq -r .candidate_id "$out/f5-candidate.json")" \
  > "$out/freeze.txt"
