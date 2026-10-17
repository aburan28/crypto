#!/bin/sh
# Authorized public raw-control custody. No mathematical executable is invoked.
set -eu
study=/Volumes/SSD990-2/crypto/isogeny-degree-challenge-20261015/research/isogeny_world_frontier_20261016/replay-optimization-v1/interruption-retake-015-v1
out=$study/custody
payload=$study/native-controls.tar.gz
store_exe=/Volumes/SSD990-2/crypto/.claude/worktrees/agent-affd575914239114f/target/release/crypto
logical=evidence/isogeny-fixed-modulus-20261016/015fa7bb504e/interruption-retake-v1
[ ! -e "$out" ]
mkdir "$out"
trap 'transfer_rc=$?; printf "%s\n" "$transfer_rc" > "$out/transfer.exit-code.txt"' EXIT
date -u +%Y-%m-%dT%H:%M:%SZ > "$out/started-utc.txt"
actual=$(/usr/bin/shasum -a 256 "$store_exe")
printf '%s\n' "$actual" > "$out/store-client-before.sha256"
actual_hash=$(printf '%s\n' "$actual" | /usr/bin/awk '{print $1}')
[ "$actual_hash" = aab28a781f807f4f3ef364336b04d2ff18a2bb8483e1c6d1494baf5c1c3a561a ]
export CRYPTO_STORE_INDEX=$out/store-index.json
export CRYPTO_STORE_REMOTE=s3://crypto-autoresearcher/store
export CRYPTO_STORE_CACHE=/Volumes/SSD990-2/crypto/isogeny-fixed-modulus-builds-20261016/interruption-retake-015-v1/custody-readback-cache
export CRYPTO_STORE_AWS=/opt/homebrew/bin/aws
unset CRYPTO_STORE_ENDPOINT
"$store_exe" store push "$payload" --name "$logical" --description 'Prospective native pure-logging interruption control raw evidence' --provenance 'TASK-20261016-7b4402; new reference/candidate logging invocations on source-bound current unit executable hashes; no historical raw recovery or mathematical work' --no-cache --json > "$out/push.json" 2> "$out/push.stderr.log"
"$store_exe" store pull "$logical" --json > "$out/pull.stdout.log" 2> "$out/pull.stderr.log"
"$store_exe" store verify --remote --json > "$out/verify.stdout.log" 2> "$out/verify.stderr.log"
readback=$(cat "$out/pull.stdout.log")
[ -f "$readback" ]
cmp "$payload" "$readback"
/usr/bin/shasum -a 256 "$readback" > "$out/readback.sha256"
/usr/bin/shasum -a 256 "$store_exe" > "$out/store-client-after.sha256"
cmp "$out/store-client-before.sha256" "$out/store-client-after.sha256"
date -u +%Y-%m-%dT%H:%M:%SZ > "$out/ended-utc.txt"
printf 'Authorized pure-logging archive push, fresh pull, remote verify and byte comparison completed.\n'
