#!/bin/zsh
set -euo pipefail

readonly P61_ADJ_RUNNER="/Users/adamburan/Library/Application Support/CryptoP61Campaign/run_adjacency_campaign-bd5ee733.py"
readonly P61_ADJ_ROOT="/Volumes/data/ec-hardness-atlas-20261012/p61-selected-adjacency-v1"
readonly P61_VERTEX_ROOT="/Volumes/data/ec-hardness-atlas-20261012/p61-class-streaming-v1"
readonly P61_ADJ_LOG="/Users/adamburan/Library/Logs/CryptoP61Campaign/adjacency-v2-from-x000000016704.log"

if [[ ! -d /Volumes/data ]]; then
  print -u2 "required UGREEN mount is absent: /Volumes/data"
  exit 1
fi
if [[ ! -f "${P61_ADJ_RUNNER}" ]]; then
  print -u2 "committed adjacency supervisor is absent: ${P61_ADJ_RUNNER}"
  exit 1
fi

exec /usr/bin/caffeinate -disu \
  /Library/Frameworks/Python.framework/Versions/3.13/Resources/Python.app/Contents/MacOS/Python \
  -u "${P61_ADJ_RUNNER}" \
  --producer "${P61_ADJ_ROOT}/bin/p61_isogeny_adjacency-darwin-arm64-f9b48fcf9" \
  --verifier "${P61_ADJ_ROOT}/bin/p61_isogeny_adjacency-darwin-arm64-92a62be64" \
  --vertex-manifest-dir "${P61_VERTEX_ROOT}/manifests" \
  --root "${P61_ADJ_ROOT}" \
  --mount-root /Volumes/data \
  --producer-source-commit f9b48fcf916bfa2b0b6fbec06e8b7c4137ceef96 \
  --verifier-source-commit 92a62be6471358edb3f9becaeb56c94c24a90b7e \
  --class-number-receipt-sha256 44e202b8eb9cdad3ef09fb50e416cc8bdcf0942a6f140e24db587a3c4d44493f \
  --window-vertices 4096 \
  --threads 4 \
  --prefix adjacency-v2 \
  --minimum-x 16704 \
  --maximum-x 19894128 \
  --min-free-bytes 1500000000000 \
  --max-logical-bytes 6000000000000 \
  --follow \
  --poll-seconds 60 \
  >> "${P61_ADJ_LOG}" 2>&1
