#!/bin/sh
set -eu
repo=/Volumes/SSD990-2/crypto/isogeny-degree-challenge-20261015
base=$repo/research/isogeny_world_frontier_20261016/replay-optimization-v1/custody-approved-20261016
out=$base/root-custody
payload=$base/root-gate-3ce9694.tar.gz
store_exe=/Volumes/SSD990-2/crypto/.claude/worktrees/agent-affd575914239114f/target/release/crypto
logical=evidence/isogeny-fixed-modulus-20261016/015fa7bb504e/merged-root-gate-3ce9694.tar.gz
export CRYPTO_STORE_INDEX=$out/store-index.json
export CRYPTO_STORE_REMOTE=s3://crypto-autoresearcher/store
export CRYPTO_STORE_CACHE=/Volumes/SSD990-2/crypto/isogeny-fixed-modulus-builds-20261016/root-gate-custody-readback-cache
export CRYPTO_STORE_AWS=/opt/homebrew/bin/aws
unset CRYPTO_STORE_ENDPOINT
mkdir "$out"
trap 'upload_rc=$?; printf "%s\n" "$upload_rc" > "$out/transfer.exit"' EXIT
date -u +%Y-%m-%dT%H:%M:%SZ > "$out/started-utc.txt"
shasum -a 256 "$store_exe" > "$out/store-client.sha256"
"$store_exe" store push "$payload" --name "$logical" --description 'Complete mandatory release library gate after unrelated main merge changes' --provenance 'TASK-20261016-de577e; approved build/control evidence custody; root source 3ce9694072773badcedc815ddbb953339888e433; frozen scientific A/B remains source015fa7bb504e' --no-cache --json > "$out/push.json" 2> "$out/push.stderr"
"$store_exe" store pull "$logical" --json > "$out/pull.stdout" 2> "$out/pull.stderr"
"$store_exe" store verify --remote --json > "$out/verify.stdout" 2> "$out/verify.stderr"
readback=$(cat "$out/pull.stdout")
cmp "$payload" "$readback"
shasum -a 256 "$readback" > "$out/readback.sha256"
date -u +%Y-%m-%dT%H:%M:%SZ > "$out/ended-utc.txt"
printf 'Merged-source root gate custody verified.\n'
