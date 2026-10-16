#!/usr/bin/env bash
# One job of the N83 precomputed-rho control: IC, plain fold rho, and the two
# precomputed-table rho arms on the same runner and the same 16 public points.
# Thin orchestration only; every measurement is a native producer's record and
# the analysis is the Rust analyzer.
#
# usage: run_job.sh <k3300|k600> <ic-bin> <rho-bin> <bl-bin> <analyze-bin> <out-dir>
set -euo pipefail
label=$1 ic_bin=$2 rho_bin=$3 bl_bin=$4 analyze_bin=$5 out=$6
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
[ ! -e "$out" ] || { echo "output exists: $out" >&2; exit 2; }
mkdir -p "$out"/{ic,rho_plain,bl_a,bl_b,host}
cp "$here/targets.jsonl" "$out/targets.jsonl"
test "$(sha256sum "$out/targets.jsonl" | cut -d' ' -f1)" = "$(cut -d' ' -f1 "$here/targets.sha256")"
rho_config="$repo/research/n83_dense_root_20261013/rho_config.json"
test "$(sha256sum "$rho_config" | cut -d' ' -f1)" = cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b
stamp() { date -u +%Y-%m-%dT%H:%M:%S.%NZ; }

# Host receipt.
{ lscpu; echo; nproc; echo; uname -a; } > "$out/host/lscpu.txt"
cat /proc/meminfo > "$out/host/meminfo.txt"
cat /sys/fs/cgroup/cpu.max > "$out/host/cgroup_cpu_max.txt" 2>/dev/null || true
sha256sum "$ic_bin" "$rho_bin" "$bl_bin" "$analyze_bin" > "$out/host/binaries.sha256"
echo "job=$label started=$(stamp)" > "$out/host/timeline.txt"

# 1. Index calculus: one process, all 16 targets; setup once, per-target online.
echo "ic start $(stamp)" >> "$out/host/timeline.txt"
if [ "$label" = k3300 ]; then
  # Exact environment of research/n83_tier1_refine_reserve_20261015/run.py (root11).
  env -i PATH="$PATH" KIC_PACKED_ROOT_TABLE=1 KIC_PACKED_ROOT_TABLE_DENSE=0 \
    KIC_PACKED_ROOT_TABLE_LOAD75=0 KIC_STATE28=0 KIC_STATE26=0 KIC_ROOT_ONLY22=0 \
    KIC_ROOT_ONLY21=0 KIC_ROOT_ONLY11=1 KIC_INDEXED_ROOT5=1 KIC_INDEXED_ROOT4=0 \
    KIC_ROOT_BLOOM=0 KIC_SYMMETRIC_S3=1 KIC_INDEXED5_LOAD75=0 KIC_INDEXED5_FIXED30=1 \
    KIC_INDEXED4_FIXED31=0 KIC_SPARSE_CANONICAL=0 KIC_RANK_THREADS=4 KIC_TARGET_THREADS=1 \
    KIC_DUMP_BASE="$out/ic/base_dump.jsonl" \
    /usr/bin/time -v -o "$out/ic/time.txt" taskset -c 0-3 \
    "$ic_bin" construct:83:1:3300 "$out/targets.jsonl" 29 "$out/ic/records.jsonl" \
    > "$out/ic/stdout.log" 2> "$out/ic/stderr.log" || echo "ic exit $?" >> "$out/host/timeline.txt"
  sha256sum "$out/ic/base_dump.jsonl" > "$out/ic/base_dump.sha256"
  rm -f "$out/ic/base_dump.jsonl"
else
  # Exact command of hosted run gha-38597929540-1 (rank seed 7, one pinned CPU).
  env RAYON_NUM_THREADS=1 /usr/bin/time -v -o "$out/ic/time.txt" taskset -c 3 \
    "$ic_bin" construct:83:1:600 "$out/targets.jsonl" 7 "$out/ic/records.jsonl" \
    > "$out/ic/stdout.log" 2> "$out/ic/stderr.log" || echo "ic exit $?" >> "$out/host/timeline.txt"
fi
echo "ic end $(stamp)" >> "$out/host/timeline.txt"

# 2. Plain fold rho (empty table per process), the Tier 1 reference worker.
echo "rho_plain start $(stamp)" >> "$out/host/timeline.txt"
i=0
while IFS= read -r line; do
  [ -n "$line" ] || continue
  w=$(printf '%02d' "$i")
  printf '%s\n' "$line" | tr -d ' ' | sed -E 's/\[([0-9]+),([0-9]+)\]/["\1","\2"]/' > "$out/rho_plain/w$w.point.jsonl"
  taskset -c 3 "$rho_bin" "$out/rho_plain/w$w.point.jsonl" $((2026101600 + i)) "$rho_config" \
    "$out/rho_plain/w$w.json" > /dev/null 2> "$out/rho_plain/w$w.stderr.log" || true
  i=$((i + 1))
done < "$out/targets.jsonl"
echo "rho_plain end $(stamp)" >> "$out/host/timeline.txt"

# 3. Precomputed rho, time-matched: one thread, table walks stop 650 s after setup start.
echo "bl_a start $(stamp)" >> "$out/host/timeline.txt"
cp "$here/configs/bl_a_time_matched.json" "$out/bl_a/config.json"
/usr/bin/time -v -o "$out/bl_a/time.txt" taskset -c 3 \
  "$bl_bin" "$out/bl_a/config.json" "$out/targets.jsonl" "$out/bl_a/records.jsonl" \
  > "$out/bl_a/stdout.log" 2> "$out/bl_a/stderr.log" || true
echo "bl_a end $(stamp)" >> "$out/host/timeline.txt"

# 4. Precomputed rho, memory-matched: table bytes = this job's IC max RSS.
ic_rss_kib=$(awk -F': ' '/Maximum resident set size/ {print $2}' "$out/ic/time.txt")
ic_rss=$((ic_rss_kib * 1024))
sed "s/IC_MAX_RSS_BYTES/$ic_rss/" "$here/configs/bl_b_memory_matched.template.json" > "$out/bl_b/config.json"
echo "bl_b start $(stamp) table_bytes=$ic_rss" >> "$out/host/timeline.txt"
/usr/bin/time -v -o "$out/bl_b/time.txt" taskset -c 0-3 \
  "$bl_bin" "$out/bl_b/config.json" "$out/targets.jsonl" "$out/bl_b/records.jsonl" \
  > "$out/bl_b/stdout.log" 2> "$out/bl_b/stderr.log" || true
echo "bl_b end $(stamp)" >> "$out/host/timeline.txt"

# 5. Independent replay and analysis.
cp "$here/reference_scalars.txt" "$out/reference_scalars.txt"
set +e
"$analyze_bin" "$out" "$label" "$out/analysis.json" > "$out/analysis.log" 2>&1
echo "$?" > "$out/analysis_exit.txt"
set -e
echo "finished $(stamp)" >> "$out/host/timeline.txt"
cat "$out/analysis.log"
