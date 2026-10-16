#!/usr/bin/env bash
# ReleaseMe creates the release. Upload archives as bytes, then verify downloads.
set -euo pipefail
if [[ $# != 4 ]]; then
  echo "usage: $0 REPOSITORY TAG ASSET_DIRECTORY FRESH_VERIFICATION_DIRECTORY" >&2
  exit 64
fi
repository=$1
release_tag=$2
asset_dir=$3
verification_dir=$4
[[ -s "$asset_dir/SHA256SUMS" ]] || { echo "missing SHA256SUMS" >&2; exit 65; }
[[ ! -e "$verification_dir" ]] || { echo "verification directory already exists" >&2; exit 73; }
check_sums() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum -c SHA256SUMS
  else shasum -a 256 -c SHA256SUMS; fi
}
(cd "$asset_dir" && check_sums)
assets=("$asset_dir"/*)
# No --clobber: an existing asset is an error, never silently replaced.
gh release upload "$release_tag" "${assets[@]}" --repo "$repository"
mkdir "$verification_dir"
gh release download "$release_tag" --repo "$repository" --dir "$verification_dir"
# Bind downloaded bytes to the locally computed manifest, including the manifest itself.
cmp "$asset_dir/SHA256SUMS" "$verification_dir/SHA256SUMS"
(cd "$verification_dir" && check_sums)
# The release stays a draft until all downloaded bytes match the staged manifest.
gh release edit "$release_tag" --repo "$repository" --draft=false --latest
echo "Published release assets match the staged SHA256SUMS."
