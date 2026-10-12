#!/bin/sh
set -u
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
cd "$repo" || exit 1
out=research/ic_solver_online_20261003/confirmation
binary=${IC_WORKER_BINARY:-/private/tmp/ic-solver-e2e-target/release/examples/ic_tournament_worker}
expected=$(sed -n 's/^binary_sha256=//p' "$out/freeze.txt")
actual=$(shasum -a 256 "$binary" | awk '{print $1}')
if [ "$actual" != "$expected" ]; then
  echo "worker binary differs from frozen digest" >&2
  exit 1
fi
mkdir -p "$out/runs"
uname -srm > "$out/runs/host.txt"
rustc --version > "$out/runs/rustc.txt"
printf '%s\n' '1 inherited_f4' '1 f5' '1 rho' \
  '2 f5' '2 inherited_f4' '2 rho' \
  '3 inherited_f4' '3 f5' '3 rho' |
while read -r repeat arm; do
  stem="$out/runs/R${repeat}-${arm}"
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.start-utc"
  RAYON_NUM_THREADS=1 IC_ARTIFACT_CACHE=off IC_F2_BACKEND=cpu \
    gtimeout 600 "$binary" < "$out/input-$arm.json" \
    > "$stem.stdout.json" 2> "$stem.stderr.txt"
  rc=$?
  printf '%s\n' "$rc" > "$stem.exit-status"
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.end-utc"
done
