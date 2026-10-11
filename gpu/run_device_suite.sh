#!/usr/bin/env bash
# Run every GPU arithmetic family that has a Modal entry point on one rented
# GPU, and keep the logs as evidence.
#
#   MODAL_TOKEN_ID=... MODAL_TOKEN_SECRET=... gpu/run_device_suite.sh [GPU] [FAMILIES]
#
#   GPU        a Modal GPU type: RTX-PRO-6000 (default), H100, H200, B200,
#              L40S, A100, L4, T4
#   FAMILIES   comma-separated subset of: ecc-secp256k1, ecc-curve25519, nist
#              (default: all three)
#
# What runs (all of it is the existing per-directory Modal code; this script
# only sequences it and files the output):
#
#   gpu/ecc,  secp256k1   selftest, bench, tune, bsgs --wbits 44
#   gpu/ecc,  curve25519  selftest, bench, bsgs --wbits 44
#   gpu/nist              validate, tune --repeats 3
#
# Output goes to research/gpu_device_suite_<UTC date>/<family>/<step>.log plus
# a provenance.json (git revision, GPU type, start/end times, exit codes).
# Nothing here reads the logs back or claims a number: a run is evidence only
# once a person has read it against the family's README and recorded it.
#
# Exit status is non-zero if any step failed.  A failed selftest stops that
# family's later steps, since a throughput number from a kernel that
# disagrees with the host reference is not a number.
set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
GPU=${1:-${GPU:-RTX-PRO-6000}}
FAMILIES=${2:-ecc-secp256k1,ecc-curve25519,nist}
STAMP=$(date -u +%Y%m%d)
OUT=${OUT:-$ROOT/research/gpu_device_suite_$STAMP}

if [[ -z "${MODAL_TOKEN_ID:-}" || -z "${MODAL_TOKEN_SECRET:-}" ]]; then
  echo "MODAL_TOKEN_ID / MODAL_TOKEN_SECRET are not set; nothing runs." >&2
  echo "Tokens: https://modal.com/settings/tokens" >&2
  exit 2
fi
if ! python3 -c 'import modal' 2>/dev/null; then
  python3 -m pip install -q 'modal>=1.1,<2' || { echo "pip install modal failed" >&2; exit 2; }
fi

mkdir -p "$OUT"
rev=$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)
dirty=$(git -C "$ROOT" status --porcelain 2>/dev/null | wc -l | tr -d ' ')
started=$(date -u +%Y-%m-%dT%H:%M:%SZ)
declare -a results=()
rc_all=0

run_step() {   # family dir step cmd...
  local family=$1 dir=$2 step=$3; shift 3
  mkdir -p "$OUT/$family"
  local log="$OUT/$family/$step.log"
  echo "== [$family] $step: $* (log: ${log#$ROOT/})"
  local t0=$(date -u +%s)
  ( cd "$ROOT/$dir" && "$@" ) 2>&1 | tee "$log"
  local rc=${PIPESTATUS[0]}
  local t1=$(date -u +%s)
  results+=("{\"family\":\"$family\",\"step\":\"$step\",\"command\":\"$*\",\"exit\":$rc,\"seconds\":$((t1-t0))}")
  [[ $rc -ne 0 ]] && rc_all=1
  return $rc
}

IFS=, read -ra fams <<< "$FAMILIES"
for fam in "${fams[@]}"; do
  case $fam in
    ecc-secp256k1)
      export ECC_GPU=$GPU ECC_CURVE=secp256k1
      run_step "$fam" gpu/ecc selftest modal run modal_app.py::selftest || continue
      run_step "$fam" gpu/ecc bench    modal run modal_app.py::bench
      run_step "$fam" gpu/ecc tune     modal run modal_app.py::tune
      run_step "$fam" gpu/ecc bsgs     modal run modal_app.py::bsgs --wbits 44
      ;;
    ecc-curve25519)
      export ECC_GPU=$GPU ECC_CURVE=curve25519
      run_step "$fam" gpu/ecc selftest modal run modal_app.py::selftest || continue
      run_step "$fam" gpu/ecc bench    modal run modal_app.py::bench
      run_step "$fam" gpu/ecc bsgs     modal run modal_app.py::bsgs --wbits 44
      ;;
    nist)
      export NIST_GPU=$GPU
      run_step "$fam" gpu/nist validate modal run modal_app.py::validate || continue
      run_step "$fam" gpu/nist tune     modal run modal_app.py::tune --repeats 3
      ;;
    *) echo "unknown family $fam" >&2; rc_all=1 ;;
  esac
done

ended=$(date -u +%Y-%m-%dT%H:%M:%SZ)
{
  echo "{"
  echo "  \"git_revision\": \"$rev\", \"dirty_files\": $dirty,"
  echo "  \"gpu\": \"$GPU\", \"families\": \"$FAMILIES\","
  echo "  \"started\": \"$started\", \"ended\": \"$ended\","
  echo "  \"steps\": [$(IFS=,; echo "${results[*]}")]"
  echo "}"
} > "$OUT/provenance.json"
echo "== done; provenance: ${OUT#$ROOT/}/provenance.json (exit $rc_all)"
exit $rc_all
