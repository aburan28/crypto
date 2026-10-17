#!/bin/sh
# Metadata-only reproduction recipe; use a verified published executable.
set -eu
release_exe=$1
out=$2
mkdir "$out"
shasum -a 256 "$release_exe" > "$out/binary.sha256"
"$release_exe" --version > "$out/version.json" 2> "$out/version.stderr"
printf '%s\n' "$?" > "$out/version.exit"
"$release_exe" curves > "$out/curves.json" 2> "$out/curves.stderr"
printf '%s\n' "$?" > "$out/curves.exit"
"$release_exe" challenge policy > "$out/policy.json" 2> "$out/policy.stderr"
printf '%s\n' "$?" > "$out/policy.exit"
"$release_exe" --help > "$out/help.txt" 2> "$out/help.stderr"
printf '%s\n' "$?" > "$out/help.exit"
