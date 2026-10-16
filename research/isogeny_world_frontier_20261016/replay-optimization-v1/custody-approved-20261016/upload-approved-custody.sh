#!/bin/sh
set -eu
repo=/Volumes/SSD990-2/crypto/isogeny-degree-challenge-20261015
custody_dir=$repo/research/isogeny_world_frontier_20261016/replay-optimization-v1/custody-approved-20261016
store_exe=/Volumes/SSD990-2/crypto/.claude/worktrees/agent-affd575914239114f/target/release/crypto
export CRYPTO_STORE_INDEX=$custody_dir/store-index.json
export CRYPTO_STORE_REMOTE=s3://crypto-autoresearcher/store
export CRYPTO_STORE_CACHE=/Volumes/SSD990-2/crypto/isogeny-fixed-modulus-builds-20261016/custody-approved-20261016/readback-cache
export CRYPTO_STORE_AWS=/opt/homebrew/bin/aws
unset CRYPTO_STORE_ENDPOINT
mkdir "$custody_dir/push" "$custody_dir/pull"
trap 'upload_rc=$?; printf "%s\n" "$upload_rc" > "$custody_dir/transfer.exit"' EXIT
cd "$repo"
date -u +%Y-%m-%dT%H:%M:%SZ > "$custody_dir/transfer-started-utc.txt"
shasum -a 256 "$store_exe" > "$custody_dir/store-client.sha256"
jq -r '.inputs[] | [.id,.path,.store_name,.sha256] | @tsv' "$custody_dir/inputs.json" > "$custody_dir/transfer-inputs.tsv"
while IFS="$(printf '\t')" read -r ident input logical expected; do
 actual=$(shasum -a 256 "$input" | cut -d ' ' -f 1)
 test "$actual" = "$expected"
 "$store_exe" store push "$input" --name "$logical" --description 'Frozen matched isogeny build/control custody; engineering evidence only' --provenance 'TASK-20261016-de577e; native source 015fa7bb504ef63953e0e04e464bb266f4f207bf; explicit approved transfer 2026-10-16; input binding inputs.json' --no-cache --json > "$custody_dir/push/$ident.json" 2> "$custody_dir/push/$ident.stderr"
 printf 'Stored %s\n' "$ident"
done < "$custody_dir/transfer-inputs.tsv"
"$store_exe" store push "$custody_dir/inputs.json" --name evidence/isogeny-fixed-modulus-20261016/015fa7bb504e/inputs.json --description 'Frozen input identities and authorization scope for matched native build custody' --provenance 'TASK-20261016-de577e; approved custody inputs 2026-10-16' --no-cache --json > "$custody_dir/push/input-manifest.json" 2> "$custody_dir/push/input-manifest.stderr"
while IFS="$(printf '\t')" read -r ident input logical expected; do
 "$store_exe" store pull "$logical" --json > "$custody_dir/pull/$ident.json" 2> "$custody_dir/pull/$ident.stderr"
 printf 'Downloaded %s\n' "$ident"
done < "$custody_dir/transfer-inputs.tsv"
"$store_exe" store pull evidence/isogeny-fixed-modulus-20261016/015fa7bb504e/inputs.json --json > "$custody_dir/pull/input-manifest.json" 2> "$custody_dir/pull/input-manifest.stderr"
"$store_exe" store verify --remote --json > "$custody_dir/remote-verification.json" 2> "$custody_dir/remote-verification.stderr"
date -u +%Y-%m-%dT%H:%M:%SZ > "$custody_dir/transfer-ended-utc.txt"
printf 'Approved custody transfer and native remote verification completed.\n'
