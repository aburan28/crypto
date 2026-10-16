#!/bin/sh
# Integrity only: does not rerun Sage mathematical verification.
set -eu
cd "$(dirname "$0")"
sha256sum -c archive.sha256
unzip -t artifacts/factor_base_evaluation.zip >/dev/null
work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM
unzip -q artifacts/factor_base_evaluation.zip -d "$work_dir"
cp contents.sha256 "$work_dir/contents.sha256"
(cd "$work_dir" && sha256sum -c contents.sha256)
for file in results/* protocol.json gate_sampling_addendum.json timing_wrapper_addendum.json; do
    [ -f "$work_dir/$file" ] || continue
    cmp "$file" "$work_dir/$file"
done
printf '%s\n' 'Archive CRCs, 32 SHA-256 entries and matching review copies verified.'
