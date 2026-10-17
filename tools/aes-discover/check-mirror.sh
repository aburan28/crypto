#!/usr/bin/env bash
# Read-only software parity check. Does not execute an experiment.
set -euo pipefail
if [ "$#" -ne 1 ]; then
  echo 'usage: bash check-mirror.sh /path/to/cryptanalysis' >&2
  exit 2
fi
local_crate=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
reference=$(CDPATH= cd -- "$1/tools/aes-discover" && pwd)
sources() {
  (cd -- "$1" && { printf '%s\n' Cargo.toml Cargo.lock; find src tests examples -type f; } | LC_ALL=C sort)
}
diff -u <(sources "$local_crate") <(sources "$reference")
count=0
while IFS= read -r path; do
  cmp -- "$local_crate/$path" "$reference/$path"
  count=$((count + 1))
done < <(sources "$local_crate")
printf 'Byte-identical companion sources, tests, example and Cargo files: %s\n' "$count"
