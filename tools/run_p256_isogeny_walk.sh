#!/usr/bin/env bash
# Run the frozen P-256 2^12-edge feasibility gate end to end.

set -euo pipefail

readonly STEPS=4096
readonly EXPECTED_JSONL_RECORDS=4100

usage() {
    cat <<'EOF'
Usage: tools/run_p256_isogeny_walk.sh [options]

Build the native Rust walker, cross-check its independent division-polynomial/
Velu fixtures, generate exactly 4,096 deterministic P-256 isogeny edges, and
replay-verify every edge. No AWS, S3, Cairn, or network access is used.

Options:
  --out-dir PATH    New directory for the certificate and run records.
                    Default: ./p256-isogeny-walk-run-<UTC timestamp>
  --cpus LIST       Linux CPU list reserved for each timed command, such as
                    "4" or "2-3". Default: choose one complete physical core.
  --no-isolation    Run without CPU reservation. Correctness is still checked,
                    but wall time and throughput are not rigorous evidence.
  --skip-build      Use an existing release binary in CARGO_TARGET_DIR.
  -h, --help        Show this help.

The gate takes about 16 minutes on the host used for the committed reference
run. It stops immediately on any fixture, edge, twist, digest, or replay error.
EOF
}

fail() {
    printf 'p256-isogeny-gate: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || fail "required command not found: $1"
}

expand_cpu_list() {
    local list=$1
    local part
    local first
    local last
    local cpu
    local -a parts

    IFS=',' read -r -a parts <<<"$list"
    for part in "${parts[@]}"; do
        if [[ $part == *-* ]]; then
            first=${part%-*}
            last=${part#*-}
            for ((cpu = first; cpu <= last; cpu++)); do
                printf '%s\n' "$cpu"
            done
        else
            printf '%s\n' "$part"
        fi
    done
}

# Select a full SMT sibling set while leaving at least one logical CPU for
# other work. tools/isolated_bench.py performs the authoritative validation.
choose_isolated_cpus() {
    local allowed_list
    local cpu
    local sibling_list
    local sibling
    local all_allowed
    local -a allowed
    local -a siblings
    local -A allowed_set=()

    [[ -r /proc/self/status ]] || return 1
    allowed_list=$(awk '/^Cpus_allowed_list:/ { print $2 }' /proc/self/status)
    [[ -n $allowed_list ]] || return 1
    mapfile -t allowed < <(expand_cpu_list "$allowed_list")
    ((${#allowed[@]} >= 2)) || return 1
    for cpu in "${allowed[@]}"; do
        allowed_set[$cpu]=1
    done

    for cpu in "${allowed[@]}"; do
        sibling_list=$cpu
        if [[ -r /sys/devices/system/cpu/cpu"$cpu"/topology/thread_siblings_list ]]; then
            sibling_list=$(<"/sys/devices/system/cpu/cpu${cpu}/topology/thread_siblings_list")
        fi
        mapfile -t siblings < <(expand_cpu_list "$sibling_list")
        ((${#siblings[@]} < ${#allowed[@]})) || continue
        all_allowed=1
        for sibling in "${siblings[@]}"; do
            if [[ -z ${allowed_set[$sibling]+present} ]]; then
                all_allowed=0
                break
            fi
        done
        if ((all_allowed)); then
            printf '%s\n' "$sibling_list"
            return 0
        fi
    done
    return 1
}

out_dir=
cpus=auto
isolate=1
skip_build=0

while (($#)); do
    case "$1" in
        --out-dir)
            (($# >= 2)) || fail '--out-dir requires a path'
            out_dir=$2
            shift 2
            ;;
        --cpus)
            (($# >= 2)) || fail '--cpus requires a CPU list'
            cpus=$2
            shift 2
            ;;
        --no-isolation)
            isolate=0
            shift
            ;;
        --skip-build)
            skip_build=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            fail "unknown argument: $1"
            ;;
    esac
done

script_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
readonly script_dir
repo_root=$(CDPATH= cd -- "$script_dir/.." && pwd)
readonly repo_root
readonly isolation_tool="$repo_root/tools/isolated_bench.py"

require_command cargo
require_command rustc
require_command git
require_command gzip
require_command sha256sum
require_command wc
require_command cmp
require_command grep
require_command awk
require_command tee

if ((isolate)); then
    require_command python3
    [[ -f $isolation_tool ]] || fail "missing isolation tool: $isolation_tool"
    if [[ $cpus == auto ]]; then
        cpus=$(choose_isolated_cpus) || fail \
            'cannot reserve a complete physical core; pass --cpus LIST or --no-isolation'
    fi
elif [[ $cpus != auto ]]; then
    fail '--cpus and --no-isolation cannot be used together'
fi

if [[ -z $out_dir ]]; then
    out_dir=$PWD/p256-isogeny-walk-run-$(date -u +%Y%m%dT%H%M%SZ)
fi
[[ ! -e $out_dir ]] || fail "output path already exists: $out_dir"
mkdir -p -- "$out_dir"
out_dir=$(CDPATH= cd -- "$out_dir" && pwd)
readonly out_dir

target_dir=${CARGO_TARGET_DIR:-$repo_root/target}
if [[ $target_dir != /* ]]; then
    target_dir=$repo_root/$target_dir
fi
export CARGO_TARGET_DIR=$target_dir
export RAYON_NUM_THREADS=1
readonly binary="$CARGO_TARGET_DIR/release/p256_isogeny_walk"
readonly certificate="$out_dir/certificate.jsonl.gz"

printf 'P-256 isogeny gate: output directory %s\n' "$out_dir"
if ((isolate)); then
    printf 'P-256 isogeny gate: reserving CPU(s) %s\n' "$cpus"
else
    printf '%s\n' 'P-256 isogeny gate: WARNING: CPU isolation disabled; timing is diagnostic only' >&2
fi

if ((!skip_build)); then
    printf '%s\n' 'P-256 isogeny gate: building the release binary'
    if ((isolate)); then
        python3 "$isolation_tool" busy --wait -- \
            cargo build --locked --release --bin p256_isogeny_walk \
            --manifest-path "$repo_root/Cargo.toml"
    else
        cargo build --locked --release --bin p256_isogeny_walk \
            --manifest-path "$repo_root/Cargo.toml"
    fi
fi
[[ -x $binary ]] || fail "release binary not found: $binary"

{
    printf 'schema: p256-isogeny-gate-environment/v1\n'
    printf 'commit: %s\n' "$(git -C "$repo_root" rev-parse HEAD)"
    printf 'rustc: %s\n' "$(rustc --version)"
    printf 'cargo: %s\n' "$(cargo --version)"
    printf 'kernel: %s\n' "$(uname -srmo)"
    printf 'steps: %s\n' "$STEPS"
    printf 'rayon_threads: %s\n' "$RAYON_NUM_THREADS"
    if ((isolate)); then
        printf 'isolation: tools/isolated_bench.py\n'
        printf 'reserved_cpus: %s\n' "$cpus"
    else
        printf 'isolation: none\n'
    fi
} >"$out_dir/environment.txt"

printf '%s\n' 'P-256 isogeny gate: checking independent division-polynomial/Velu fixtures'
"$binary" fixture | tee "$out_dir/fixture.json"

printf 'P-256 isogeny gate: generating and internally verifying %s edges\n' "$STEPS"
if ((isolate)); then
    python3 "$isolation_tool" run --wait --cpus "$cpus" \
        --out "$out_dir/isolated-generate.jsonl" --label p256-isogeny-generate -- \
        "$binary" generate --steps "$STEPS" --output "$certificate" \
        >"$out_dir/generate-summary.json"
else
    "$binary" generate --steps "$STEPS" --output "$certificate" \
        >"$out_dir/generate-summary.json"
fi

printf '%s\n' 'P-256 isogeny gate: replay-verifying the completed certificate'
if ((isolate)); then
    python3 "$isolation_tool" run --wait --cpus "$cpus" \
        --out "$out_dir/isolated-verify.jsonl" --label p256-isogeny-verify -- \
        "$binary" verify --input "$certificate" \
        >"$out_dir/verify-summary.json"
else
    "$binary" verify --input "$certificate" >"$out_dir/verify-summary.json"
fi

gzip -t "$certificate"
record_count=$(gzip -dc "$certificate" | wc -l)
record_count=${record_count//[[:space:]]/}
[[ $record_count == "$EXPECTED_JSONL_RECORDS" ]] || fail \
    "certificate has $record_count JSONL records; expected $EXPECTED_JSONL_RECORDS"
edge_count=$(gzip -dc "$certificate" | grep -c '^{"record":"edge",')
[[ $edge_count == "$STEPS" ]] || fail \
    "certificate has $edge_count edge records; expected $STEPS"
cmp -s "$out_dir/generate-summary.json" "$out_dir/verify-summary.json" || fail \
    'generation and standalone replay summaries differ'
grep -q '"emitted_steps": 4096' "$out_dir/verify-summary.json" || fail \
    'replay summary does not report 4,096 emitted edges'
grep -q '"verified_steps": 4096' "$out_dir/verify-summary.json" || fail \
    'replay summary does not report 4,096 verified edges'

printf '%s\n' PASS >"$out_dir/STATUS"

(
    cd "$out_dir"
    sha256sum STATUS certificate.jsonl.gz environment.txt fixture.json \
        generate-summary.json verify-summary.json >SHA256SUMS
    if ((isolate)); then
        sha256sum isolated-generate.jsonl isolated-verify.jsonl >>SHA256SUMS
    fi
)

printf '\n%s\n' 'P-256 isogeny gate: PASS'
cat "$out_dir/verify-summary.json"
printf 'certificate: %s\n' "$certificate"
printf 'certificate_sha256: %s\n' "$(sha256sum "$certificate" | awk '{ print $1 }')"
printf 'artifacts: %s\n' "$out_dir"
