#!/usr/bin/env bash
# Thin orchestration for PROTOCOL.md: everything measured is native
# (`ecbench exec`, Valgrind Callgrind); this script only sequences runs.
#
#   run.sh build              baseline (main fa80835a) and candidate binaries
#   run.sh inputs             export every frozen child input
#   run.sh gate               identical-output gate on every input
#   run.sh ir [JOBS]          Callgrind total Ir per (binary, input)
#   run.sh wall [ROUNDS]      pinned ABAB timing with an A/A pair
#
# Outputs go to $OUT (default results/run-1); nothing is overwritten.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
OUT=${OUT:-$HERE/results/run-1}
BASE_REV=fa80835af5f7c18bd400ccb2002f5e6bd8f6c058
BIN=$OUT/bin
LADDER=$ROOT/research/f6_ic_ecbench_ladder_20261005

# set, spec, sequences
SETS=(
  "L13 $LADDER/SPEC-n13.json 24 25 27 28 30 31 33 34 36 37 39 40 42 43 45 46"
  "L19 $LADDER/SPEC-n19.json 24 25 27 28 30 31 33 34 36 37 39 40 42 43 45 46"
  "L23 $LADDER/SPEC-n23.json 24 25 27 28"
  "H13 $HERE/SPEC-holdout-n13.json 0 1 2 3 4 5 6 7"
  "H19 $HERE/SPEC-holdout-n19.json 0 1 2 3"
)
IR_INPUTS="L13-* H13-* L19-24 L19-25 L19-27 L19-28 H19-*"
WALL_INPUTS="L19-25 L19-24 L23-25"

strip() { jq -S -f "$HERE/strip.jq"; }

case "${1:-}" in
build)
  mkdir -p "$BIN"
  [ -e "$BIN/ecbench.base" ] && { echo "refusing to overwrite $BIN"; exit 1; }
  wt=$(mktemp -d)
  git -C "$ROOT" worktree add --detach "$wt" "$BASE_REV" >/dev/null
  CARGO_TARGET_DIR="$ROOT/target/ab-base" cargo build --release --bin ecbench --manifest-path "$wt/Cargo.toml"
  cp "$ROOT/target/ab-base/release/ecbench" "$BIN/ecbench.base"
  cp "$BIN/ecbench.base" "$BIN/ecbench.base-copy"
  git -C "$ROOT" worktree remove --force "$wt"
  cargo build --release --bin ecbench --manifest-path "$ROOT/Cargo.toml"
  cp "$ROOT/target/release/ecbench" "$BIN/ecbench.cand"
  {
    echo "candidate_revision $(git -C "$ROOT" rev-parse HEAD) (+ working tree: $(git -C "$ROOT" status --porcelain src | wc -l) changed src files)"
    echo "baseline_revision $BASE_REV"
    rustc --version
    uname -srm
    grep -m1 'model name' /proc/cpuinfo
    echo "logical_cpus $(nproc)"
    grep MemTotal /proc/meminfo
    echo "features $(grep -m1 '^flags' /proc/cpuinfo | grep -o -w 'popcnt\|avx2\|avx512f\|bmi2\|pclmulqdq' | sort -u | tr '\n' ' ')"
    valgrind --version
  } > "$OUT/host.txt"
  (cd "$BIN" && sha256sum ecbench.base ecbench.base-copy ecbench.cand) > "$OUT/binaries.sha256"
  ;;
inputs)
  mkdir -p "$OUT/inputs"
  for s in "${SETS[@]}"; do
    read -r name spec seqs <<<"$s"
    for q in $seqs; do
      "$BIN/ecbench.base" profile-input --spec "$spec" --seq "$q" > "$OUT/inputs/$name-$q.json"
    done
  done
  (cd "$OUT/inputs" && sha256sum *.json) > "$OUT/inputs.sha256"
  ;;
gate)
  f=$OUT/gate.tsv
  [ -e "$f" ] && { echo "refusing to overwrite $f"; exit 1; }
  printf 'input\tverdict\tstripped_sha256\trecovered\tsolver_ops\ttotal_gae\n' > "$f"
  for in in "$OUT"/inputs/*.json; do
    id=$(basename "$in" .json)
    a=$("$BIN/ecbench.base" exec < "$in" | strip)
    b=$("$BIN/ecbench.cand" exec < "$in" | strip)
    v=DIFFER; [ -n "$a" ] && [ "$a" == "$b" ] && v=IDENTICAL
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$v" "$(sha256sum <<<"$a" | cut -c1-64)" \
      "$(jq -r '.report.recovered' <<<"$a")" \
      "$(jq -r '[.report.phases[]|select(.name=="relations")|.native.solver_ops][0]' <<<"$a")" \
      "$(jq -r '.report.total_gae' <<<"$a")" | tee -a "$f"
    [ "$v" == IDENTICAL ] || diff <(echo "$a") <(echo "$b") > "$OUT/gate-$id.diff" || true
  done
  ;;
ir)
  jobs=${2:-4}
  f=$OUT/ir.tsv
  [ -e "$f" ] && { echo "refusing to overwrite $f"; exit 1; }
  mkdir -p "$OUT/ir"
  list=()
  for pat in $IR_INPUTS; do for in in "$OUT"/inputs/$pat.json; do list+=("$(basename "$in" .json)"); done; done
  printf '%s\n' "${list[@]}" | while read -r id; do for b in base cand; do echo "$id $b"; done; done |
    xargs -P "$jobs" -L 1 bash -c '
      id=$0; b=$1
      valgrind --tool=callgrind --cache-sim=no --branch-sim=no \
        --callgrind-out-file="'"$OUT"'/ir/$id.$b.cg" "'"$BIN"'/ecbench.$b" exec \
        < "'"$OUT"'/inputs/$id.json" > /dev/null 2> "'"$OUT"'/ir/$id.$b.err"
      ir=$(grep -o "Collected : [0-9]*" "'"$OUT"'/ir/$id.$b.err" | awk "{print \$3}")
      rm -f "'"$OUT"'/ir/$id.$b.cg"*
      printf "%s\t%s\t%s\n" "$id" "$b" "$ir"' > "$OUT/ir.raw.tsv"
  { printf 'input\tbase_ir\tcand_ir\tratio\n'
    sort "$OUT/ir.raw.tsv" | awk -F'\t' '{v[$1,$2]=$3; ids[$1]=1}
      END {for (i in ids) printf "%s\t%s\t%s\t%.4f\n", i, v[i,"base"], v[i,"cand"], v[i,"cand"]/v[i,"base"]}' | sort
  } > "$f"
  cat "$f"
  ;;
wall)
  rounds=${2:-5}
  f=$OUT/wall.tsv
  [ -e "$f" ] && { echo "refusing to overwrite $f"; exit 1; }
  printf 'input\tround\tbinary\tcpu_ns\tsolve_wall_ns\trun_delay_ns\ttimeslices\n' > "$f"
  for id in $WALL_INPUTS; do
    for r in $(seq 1 "$rounds"); do
      for b in base cand base-copy cand; do
        o=$(taskset -c 3 "$BIN/ecbench.$b" exec < "$OUT/inputs/$id.json")
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$r" "$b" \
          "$(jq -r '.. | objects | select(has("cpu_ns")) | .cpu_ns' <<<"$o" | head -1)" \
          "$(jq -r '.. | objects | select(has("solve_wall_ns")) | .solve_wall_ns' <<<"$o" | head -1)" \
          "$(jq -r '.. | objects | select(has("run_delay_ns")) | .run_delay_ns' <<<"$o" | head -1)" \
          "$(jq -r '.. | objects | select(has("timeslices")) | .timeslices' <<<"$o" | head -1)" | tee -a "$f"
      done
    done
  done
  ;;
*)
  sed -n '2,12p' "$0"; exit 1 ;;
esac
