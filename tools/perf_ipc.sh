#!/usr/bin/env bash
# Hardware-counter measurement of one command with Linux `perf stat`:
# user-space cycles, instructions, IPC, branches and branch misses, and task
# clock, over N repeats, written as JSON.  It never estimates a counter the
# host cannot count: on a host without a hardware PMU (most cloud VMs,
# including Claude Code's cloud containers) it writes a record with
# "status": "unsupported" and the reason, and exits 3.
#
#   tools/perf_ipc.sh [--repeat N] [--label L] [--out FILE] [--isolate]
#                     [--events LIST] [--steps K] -- COMMAND [ARGS...]
#
#   --repeat N    runs of COMMAND (default 5); one record per run, then a
#                 summary with medians, minima and maxima
#   --label L     free-text label stored in every record
#   --out FILE    append the JSON lines there (default: stdout)
#   --isolate     run each repeat under isolated_bench (AGENTS.md section 10):
#                 $ISOLATED_BENCH, else target/release/isolated_bench;
#                 its isolation record goes to FILE.isolation.jsonl
#   --events L    extra comma-separated perf events, counted beside the core set
#   --selftest    software events only (task-clock, faults, switches): checks
#                 the parsing and JSON on a host without counters; never IPC
#   --steps K     the command's unit of work (walk steps, relations, ...);
#                 adds cycles and instructions per unit to the summary
#
# Environment: PERF (perf binary), MAX_OTHER_CPU (isolated_bench, default 0.1).
# Exit status: 0 measured, 2 usage error, 3 no hardware counters,
# 4 the command failed, 5 counts were multiplexed (scaled) in some run.
set -euo pipefail

repeat=5 label="" out="" isolate=0 extra="" steps="" selftest=0
while [ $# -gt 0 ]; do
    case $1 in
        --repeat) repeat=$2; shift 2 ;;
        --label) label=$2; shift 2 ;;
        --out) out=$2; shift 2 ;;
        --isolate) isolate=1; shift ;;
        --events) extra=$2; shift 2 ;;
        --steps) steps=$2; shift 2 ;;
        --selftest) selftest=1; shift ;;
        --) shift; break ;;
        *) echo "perf_ipc: unknown option $1" >&2; exit 2 ;;
    esac
done
[ $# -gt 0 ] || { echo "perf_ipc: no command (put it after --)" >&2; exit 2; }
command -v jq > /dev/null || { echo "perf_ipc: needs jq" >&2; exit 2; }

emit() { if [ -n "$out" ]; then printf '%s\n' "$1" >> "$out"; else printf '%s\n' "$1"; fi; }

# A perf binary that runs on this kernel.  Ubuntu's /usr/bin/perf is a wrapper
# that refuses when no linux-tools package matches the running kernel; the
# versioned binaries under /usr/lib/linux-tools-*/ work on other kernels.
find_perf() {
    local p
    for p in "${PERF:-}" "$(command -v perf || true)" $(ls -d /usr/lib/linux-tools-*/perf /usr/lib/linux-tools/*/perf 2>/dev/null | sort -rV); do
        [ -n "$p" ] && [ -x "$p" ] && "$p" --version > /dev/null 2>&1 && { echo "$p"; return 0; }
    done
    return 1
}

host_json() {
    local cpu pmus
    cpu=$(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ *//' || true)
    [ -n "$cpu" ] || cpu=$(grep -m1 -i 'CPU part' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ *//' || echo unknown)
    pmus=$(ls /sys/bus/event_source/devices 2>/dev/null | tr '\n' ' ')
    jq -nc --arg cpu "$cpu" --arg kernel "$(uname -r)" --arg arch "$(uname -m)" \
        --arg pmus "$pmus" --arg paranoid "$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unknown)" \
        --argjson cpus "$(nproc)" --arg perf "${PERF_VERSION:-}" \
        '{cpu_model: $cpu, kernel: $kernel, arch: $arch, logical_cpus: $cpus,
          event_sources: ($pmus | split(" ") | map(select(. != ""))),
          perf_event_paranoid: $paranoid, perf_version: $perf}'
}

unsupported() { # reason
    emit "$(jq -nc --arg label "$label" --arg reason "$1" --argjson host "$(host_json)" \
        --args '{schema: "perf-ipc/v1", kind: "unsupported", status: "unsupported", label: $label,
                 reason: $reason, command: $ARGS.positional, host: $host}' -- "${cmd[@]}")"
    echo "perf_ipc: hardware counters unavailable: $1" >&2
    exit 3
}

cmd=("$@")
PERF_BIN=$(find_perf) || unsupported "no working perf binary (Ubuntu: apt-get install linux-tools-common linux-tools-generic)"
PERF_VERSION=$("$PERF_BIN" --version 2>/dev/null)

if [ "$selftest" = 1 ]; then
    # Exercise the parsing and summary paths on any host: software events
    # only, labelled so, never a measurement of cycles or IPC.
    label="selftest${label:+ $label}"
    extra="page-faults,context-switches"
fi
# Probe: can this process count user-space cycles and instructions at all?
probe=$("$PERF_BIN" stat -x, -e cycles:u,instructions:u -- true 2>&1 || true)
if [ "$selftest" = 0 ] && grep -q '<not supported>' <<< "$probe"; then
    unsupported "perf reports cycles:u/instructions:u <not supported>: no hardware PMU is exposed (event sources: $(ls /sys/bus/event_source/devices | tr '\n' ' '))"
fi
if grep -qiE 'permission|access to performance monitoring' <<< "$probe"; then
    unsupported "perf_event_paranoid=$(cat /proc/sys/kernel/perf_event_paranoid) forbids counting; set it to 2 or lower"
fi

events="cycles:u,instructions:u,branches:u,branch-misses:u,task-clock"
[ "$selftest" = 1 ] && events="task-clock"
[ -n "$extra" ] && events="$events,$extra"
bench=""
if [ "$isolate" = 1 ]; then
    bench=${ISOLATED_BENCH:-target/release/isolated_bench}
    [ -x "$bench" ] || { echo "perf_ipc: --isolate needs \$ISOLATED_BENCH or $bench" >&2; exit 2; }
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
runs=()
status=0
for i in $(seq 1 "$repeat"); do
    csv="$tmp/run$i.csv"
    set +e
    if [ -n "$bench" ]; then
        "$bench" run --wait --cpus 3 --max-other-cpu "${MAX_OTHER_CPU:-0.1}" \
            --label "${label:-perf_ipc}-r$i" --out "${out:-perf_ipc}.isolation.jsonl" -- \
            "$PERF_BIN" stat -x, --no-big-num -e "$events" -o "$csv" -- "${cmd[@]}" > "$tmp/stdout$i" 2> "$tmp/stderr$i"
    else
        "$PERF_BIN" stat -x, --no-big-num -e "$events" -o "$csv" -- "${cmd[@]}" > "$tmp/stdout$i" 2> "$tmp/stderr$i"
    fi
    rc=$?
    set -e
    [ "$rc" = 0 ] || status=4
    # perf -x, rows: value,unit,event,run-time,percent-running,...
    rec=$(grep -v '^#' "$csv" | grep -v '^$' | jq -Rsc --arg label "$label" --argjson run "$i" --argjson rc "$rc" '
        split("\n") | map(select(length > 0) | split(",")) |
        map({event: .[2], raw: .[0], unit: .[1], running_pct: (.[4] | tonumber? // null)}) as $rows |
        ($rows | map({key: (.event | sub(":u$"; "")), value: (.raw | tonumber? // null)}) | from_entries) as $v |
        {schema: "perf-ipc/v1", kind: "run", label: $label, run: $run, exit_status: $rc,
         counters: $v,
         not_counted: ($rows | map(select(.raw | startswith("<"))) | map(.event)),
         multiplexed: ($rows | map(select(.running_pct != null and .running_pct < 100)) | map(.event)),
         ipc: (if ($v.cycles // 0) > 0 and $v.instructions != null then $v.instructions / $v.cycles else null end),
         branch_miss_rate: (if ($v.branches // 0) > 0 and $v["branch-misses"] != null then $v["branch-misses"] / $v.branches else null end)}')
    emit "$rec"
    runs+=("$rec")
    [ "$(jq '.multiplexed | length' <<< "$rec")" = 0 ] || { [ "$status" = 0 ] && status=5; }
done

summary=$(printf '%s\n' "${runs[@]}" | jq -sc --arg label "$label" --argjson host "$(host_json)" \
    --arg steps "$steps" --arg bench "${bench:+isolated_bench}" --args '
    def med: sort | if length == 0 then null elif length % 2 == 1 then .[length/2|floor] else (.[length/2-1] + .[length/2]) / 2 end;
    def st(f): (map(f) | map(select(. != null))) as $x | {median: ($x | med), min: ($x | min), max: ($x | max), n: ($x | length)};
    {schema: "perf-ipc/v1", kind: "summary", label: $label, command: $ARGS.positional,
     repeats: length, failed_runs: (map(select(.exit_status != 0)) | length),
     multiplexed_runs: (map(select(.multiplexed | length > 0)) | length),
     isolation: (if $bench == "" then "none (not under isolated_bench)" else $bench end),
     ipc: st(.ipc), cycles: st(.counters.cycles), instructions: st(.counters.instructions),
     branch_miss_rate: st(.branch_miss_rate), task_clock_ms: st(.counters["task-clock"]),
     per_step: (if $steps == "" then null else ($steps | tonumber) as $k |
        {steps: $k, cycles: st(if .counters.cycles == null then null else .counters.cycles / $k end),
                  instructions: st(if .counters.instructions == null then null else .counters.instructions / $k end)} end),
     host: $host}' -- "${cmd[@]}")
emit "$summary"
jq -r '"perf_ipc: IPC median \(.ipc.median) [\(.ipc.min)-\(.ipc.max)] over \(.ipc.n) runs; cycles median \(.cycles.median); instructions median \(.instructions.median)"' <<< "$summary" >&2
exit "$status"
