#!/bin/bash
# Build an independent public replay without replacing any running library.
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/.." && pwd)
if [[ $# -ne 1 ]]; then echo 'Usage: build_field_public.sh NEW_OUTPUT_DIRECTORY' >&2; exit 2; fi
out=$(realpath -m -- "$1")
if [[ -e "$out" ]]; then echo 'Refuse output overwrite' >&2; exit 2; fi
mkdir -p "$out"
modules=(modarith util group_common group_zp group_ec ec_tau rho bsgs probe)
inputs=()
for module in "${modules[@]}"; do inputs+=("$root/field_backend/src/$module.c"); done
cc -O3 -march=native -ftls-model=initial-exec -std=c11 -D_POSIX_C_SOURCE=200809L -fPIC -shared \
  -I"$root/field_backend/include" -I"$root/field_backend/src" "${inputs[@]}" \
  -pthread -lm -o "$out/libcryptanalysis.so"
g++ -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -I"$root/field_backend/include" \
  "$root/source/field_replay_public.cpp" -L"$out" -lcryptanalysis \
  -Wl,-rpath,'$ORIGIN' -pthread -o "$out/field_replay_public"
sha256sum "$out/libcryptanalysis.so" "$out/field_replay_public" > "$out/BUILD_SHA256SUMS"
echo "Built independent answer-free replay in $out; original binaries and frozen bindings are unchanged."
