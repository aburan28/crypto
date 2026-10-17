#!/usr/bin/env bash
# Exercise byte custody and failed publication with a local gh transport fixture.
set -euo pipefail
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
fixture=$(mktemp -d)
trap 'rm -rf -- "$fixture"' EXIT
mkdir "$fixture/bin" "$fixture/payload" "$fixture/staged"
printf '\000\001\177\200\201\376\377\n' > "$fixture/payload/bytes"
tar -czf "$fixture/staged/isogeny-control.tar.gz" -C "$fixture/payload" bytes
(cd "$fixture/payload" && zip -q "$fixture/staged/isogeny-control.zip" bytes)
sums() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$@"
  else shasum -a 256 "$@"; fi
}
(cd "$fixture/staged" && sums isogeny-*) > "$fixture/SHA256SUMS"
mv "$fixture/SHA256SUMS" "$fixture/staged/"
cat > "$fixture/bin/gh" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "$MOCK_ROOT/calls"
[[ "$1" = release ]] || exit 90
operation=$2
shift 3
case "$operation" in
  upload)
    [[ ${MOCK_UPLOAD_FAIL:-0} = 0 ]] || exit 31
    mkdir "$MOCK_ROOT/remote"
    for argument in "$@"; do
      [[ "$argument" != --clobber ]] || exit 91
      [[ "$argument" != --repo ]] || break
      cp "$argument" "$MOCK_ROOT/remote/"
    done
    ;;
  download)
    [[ ${MOCK_DOWNLOAD_FAIL:-0} = 0 ]] || exit 32
    destination=''
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --repo) shift 2 ;;
        --dir) destination=$2; shift 2 ;;
        *) exit 92 ;;
      esac
    done
    cp "$MOCK_ROOT/remote/"* "$destination/"
    if [[ ${MOCK_CORRUPT:-0} = 1 ]]; then
      printf '\377' >> "$destination/isogeny-control.tar.gz"
    fi
    if [[ ${MOCK_REWRITE_SUMS:-0} = 1 ]]; then
      (cd "$destination"; if command -v sha256sum >/dev/null 2>&1; then
        sha256sum isogeny-*; else shasum -a 256 isogeny-*; fi) > "$destination/SHA256SUMS"
    fi
    ;;
  edit)
    [[ ${MOCK_EDIT_FAIL:-0} = 0 ]] || exit 33
    [[ "$*" = '--repo fixture/repo --draft=false --latest' ]] || exit 94
    touch "$MOCK_ROOT/published"
    ;;
  *) exit 93 ;;
esac
MOCK
chmod +x "$fixture/bin/gh"
export PATH="$fixture/bin:$PATH"
export MOCK_ROOT
for case_name in success corrupted_download replaced_manifest upload_failure download_failure publication_failure local_mismatch existing_output existing_assets; do
  MOCK_ROOT="$fixture/$case_name"
  mkdir "$MOCK_ROOT"
  export MOCK_CORRUPT=0 MOCK_REWRITE_SUMS=0 MOCK_UPLOAD_FAIL=0 MOCK_DOWNLOAD_FAIL=0 MOCK_EDIT_FAIL=0
  staged="$fixture/staged"
  case "$case_name" in
    corrupted_download) MOCK_CORRUPT=1 ;;
    replaced_manifest) MOCK_CORRUPT=1; MOCK_REWRITE_SUMS=1 ;;
    upload_failure) MOCK_UPLOAD_FAIL=1 ;;
    download_failure) MOCK_DOWNLOAD_FAIL=1 ;;
    publication_failure) MOCK_EDIT_FAIL=1 ;;
    local_mismatch)
      cp -R "$fixture/staged" "$MOCK_ROOT/staged"
      staged="$MOCK_ROOT/staged"
      printf '\377' >> "$staged/isogeny-control.zip"
      ;;
    existing_output) mkdir "$MOCK_ROOT/verification" ;;
    existing_assets) mkdir "$MOCK_ROOT/remote" ;;
  esac
  if bash "$script_dir/publish_release_assets.sh" fixture/repo v0.0.0 "$staged" "$MOCK_ROOT/verification" > "$MOCK_ROOT/stdout" 2> "$MOCK_ROOT/stderr"; then
    [[ "$case_name" = success ]] || { echo "unexpected success: $case_name" >&2; exit 1; }
    for archive in isogeny-control.tar.gz isogeny-control.zip SHA256SUMS; do
      cmp "$staged/$archive" "$MOCK_ROOT/verification/$archive"
    done
  else
    [[ "$case_name" != success ]] || { cat "$MOCK_ROOT/stderr" >&2; exit 1; }
  fi
  if [[ "$case_name" = local_mismatch || "$case_name" = existing_output ]]; then
    [[ ! -e "$MOCK_ROOT/calls" ]] || { echo "published before local admission" >&2; exit 1; }
  fi
  if [[ "$case_name" = success ]]; then
    [[ -e "$MOCK_ROOT/published" ]] || { echo "release was not published" >&2; exit 1; }
  else
    [[ ! -e "$MOCK_ROOT/published" ]] || { echo "published despite failure: $case_name" >&2; exit 1; }
  fi
  printf 'PASS %s\n' "$case_name"
done
