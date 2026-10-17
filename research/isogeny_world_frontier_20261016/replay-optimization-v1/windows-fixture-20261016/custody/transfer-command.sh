#!/bin/sh
set -eu
repo=/Volumes/SSD990-2/crypto/isogeny-degree-challenge-20261015
out=$repo/research/isogeny_world_frontier_20261016/replay-optimization-v1/windows-fixture-20261016/custody
store_exe=/Volumes/SSD990-2/crypto/.claude/worktrees/agent-affd575914239114f/target/release/crypto
payload=$repo/research/isogeny_world_frontier_20261016/replay-optimization-v1/windows-fixture-20261016/native-controls.tar.gz
logical=evidence/isogeny-fixed-modulus-20261016/015fa7bb504e/windows-fixture-native-controls-20261016.tar.gz
export CRYPTO_STORE_INDEX=$out/store-index.json
export CRYPTO_STORE_REMOTE=s3://crypto-autoresearcher/store
export CRYPTO_STORE_CACHE=/Volumes/SSD990-2/crypto/isogeny-fixed-modulus-builds-20261016/windows-custody-readback-cache
export CRYPTO_STORE_AWS=/opt/homebrew/bin/aws
unset CRYPTO_STORE_ENDPOINT
mkdir "$out"
trap 'upload_rc=$?; printf "%s\n" "$upload_rc" > "$out/transfer.exit"' EXIT
date -u +%Y-%m-%dT%H:%M:%SZ > "$out/started-utc.txt"
shasum -a 256 "$store_exe" > "$out/store-client.sha256"
"$store_exe" store push "$payload" --name "$logical" --description 'Frozen native valid/cache-negative retake after exact Windows fixture checkout repair' --provenance 'TASK-20261016-de577e; approved control evidence custody; native source 015fa7bb504ef63953e0e04e464bb266f4f207bf; full member manifest committed in PR1991' --no-cache --json > "$out/push.json" 2> "$out/push.stderr"
"$store_exe" store pull "$logical" --json > "$out/pull.stdout" 2> "$out/pull.stderr"
"$store_exe" store verify --remote --json > "$out/verify.stdout" 2> "$out/verify.stderr"
readback=$(cat "$out/pull.stdout")
cmp "$payload" "$readback"
shasum -a 256 "$readback" > "$out/readback.sha256"
date -u +%Y-%m-%dT%H:%M:%SZ > "$out/ended-utc.txt"
printf 'Supplemental native control custody verified.\n'
