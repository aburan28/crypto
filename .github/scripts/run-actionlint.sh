#!/usr/bin/env bash
# GitHub workflow semantics only; Python/Rust suites own language validation.
set -euo pipefail

ci_lint_temp=$(mktemp -d)
trap 'rm -rf "$ci_lint_temp"' EXIT
archive="$ci_lint_temp/actionlint_1.7.12_linux_amd64.tar.gz"
curl --fail --silent --show-error --location --retry 2 \
  --connect-timeout 10 --max-time 60 \
  https://github.com/rhysd/actionlint/releases/download/v1.7.12/actionlint_1.7.12_linux_amd64.tar.gz \
  --output "$archive"
# Official release checksum, pinned with the tool version.
printf '%s  %s\n' \
  8aca8db96f1b94770f1b0d72b6dddcb1ebb8123cb3712530b08cc387b349a3d8 \
  "$archive" | sha256sum --check --status
tar -xzf "$archive" -C "$ci_lint_temp" actionlint
"$ci_lint_temp/actionlint" -color -shellcheck='' -pyflakes='' "$@"
