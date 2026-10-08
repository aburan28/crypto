#!/usr/bin/env bash
# Build the disclosed post-buffer CMS variant only after the sole natural
# preparation has completed and passed its original source-bound audit.
# This creates source/build evidence; it never admits SAT target execution.
set -euo pipefail

if (( $# != 4 )); then
  echo 'usage: build-marked-cms.sh ACCEPTED_ASSETS_TAR PREPARATION_TERMINAL PREPARATION_AUDIT NEW_OUTPUT_DIR' >&2
  exit 2
fi
bundle=$1
preparation_terminal=$2
preparation_audit=$3
output=$4
for path in "$bundle" "$preparation_terminal" "$preparation_audit" "$output"; do
  if [[ $path != /* ]]; then
    echo 'all paths must be absolute' >&2
    exit 2
  fi
done
if [[ $(uname -s) != Darwin || $(uname -m) != arm64 ]]; then
  echo 'this source-bound build recipe is for the declared physical macOS ARM64 host' >&2
  exit 2
fi
seal=3b768aa18127c2aaafb261770cfafe8eed1c7631e9aaf2a42db211b678c3e470
jq -e --arg seal "$seal" '
  .schema_version == 1 and
  .registration_sha256 == $seal and
  .source_gate_passed == true and
  .worker_drain_passed == true and
  .role_drain_passed == true and
  .source_bound_execution_admitted == false
' "$preparation_terminal" >/dev/null
jq -e --arg seal "$seal" '
  .schema_version == 1 and
  .registration_sha256 == $seal and
  .status == "PASS_NATIVE_SOURCE_BOUND_ORDINARY_PREPARATION_AUDIT" and
  .source_bound_execution_admitted == true and
  .mathematics.declared_family == "cryptominisat" and
  .mathematics.planned_queries == 512 and
  .mathematics.audited_queries == 512 and
  .mathematics.panel_complete == true and
  .mathematics.mathematical_preparation_complete == true and
  .mathematics.rank == 29 and
  .mathematics.folded_columns == 29 and
  .mathematics.usable_points == 62
' "$preparation_audit" >/dev/null

sha() { shasum -a 256 "$1" | awk '{print $1}'; }
require_sha() {
  if [[ $(sha "$1") != "$2" ]]; then
    echo "SHA-256 mismatch: $1" >&2
    exit 1
  fi
}
require_sha "$bundle" a1fd5bd49c80076f3b64fd5cb51d891b278afde765b4d3e39853692ac318bd96

script_dir=$(cd "$(dirname "$0")" && pwd -P)
patch_file="$script_dir/cms-stdin-ready-postbuffer.patch"
require_sha "$patch_file" b335c52f6673fd6a17d32d0c5506699058b3ee04f0089200305f2c89fe9fde82
source_root=$(cd "$script_dir/../../../.." && pwd -P)
output_parent=$(cd "$(dirname "$output")" && pwd -P)
preparation_parent=$(cd "$(dirname "$preparation_terminal")" && pwd -P)
bundle_parent=$(cd "$(dirname "$bundle")" && pwd -P)
if [[ $bundle == */immutable/native-archive/assets.tar.gz ]]; then
  capsule_root=$(cd "$bundle_parent/../.." && pwd -P)
else
  capsule_root=$bundle_parent
fi
for protected in "$source_root" "$preparation_parent" "$capsule_root"; do
  if [[ $output_parent == "$protected" || $output_parent == "$protected/"* ]]; then
    echo 'build output must be outside source, registration and original execution evidence' >&2
    exit 2
  fi
done
mkdir "$output"
phase=extract
finish() {
  code=$?
  trap - EXIT
  if (( code == 0 )); then verdict=BUILT_UNVALIDATED; else verdict=BUILD_FAILED; fi
  jq -n --arg status "$verdict" --arg phase "$phase" --arg seal "$seal" \
    --argjson exit_code "$code" \
    '{schema_version:1,status:$status,phase:$phase,exit_code:$exit_code,
      preparation_registration_sha256:$seal,source_bound_execution_admitted:false,
      transport_parity_passed:false,target_execution_admitted:false,
      online_speedup:null}' > "$output/terminal.json"
  exit "$code"
}
trap finish EXIT

tar -xOzf "$bundle" cms/source.tar > "$output/source.tar"
tar -xOzf "$bundle" cms/cadical.tar > "$output/cadical.tar"
tar -xOzf "$bundle" cms/cadiback.tar > "$output/cadiback.tar"
require_sha "$output/source.tar" 467b1c3d00a7d6e893332b4d8b42c6326301974d22885aa745f1f926da050323
require_sha "$output/cadical.tar" 8264713f3dc1c4455162d2912238712bd8030fceabec0f4b430d106b5a58058d
require_sha "$output/cadiback.tar" e0aa8f5d67c04527135dde5fe5f943e672d6af11a6f7f924e9f0c07ad3bffba0
cp "$patch_file" "$output/postbuffer.patch"
mkdir "$output/source" "$output/cadical" "$output/cadiback"
tar -xf "$output/source.tar" -C "$output/source"
tar -xf "$output/cadical.tar" -C "$output/cadical"
tar -xf "$output/cadiback.tar" -C "$output/cadiback"
(cd "$output/source" && patch -p1 -i "$output/postbuffer.patch") > "$output/patch.stdout" 2> "$output/patch.stderr"
require_sha "$output/source/src/main.cpp" b6c65e963442cc63df10b0d4888737693aeebccb94b6d2379b3fb2bef37160b2
require_sha "$output/source/src/dimacsparser.h" a3919c1761e77243ed526c3be38b620f111da071a7ee098595e88ab3980a2355
require_sha "$output/source/src/streambuffer.h" e58dce81a60271641ab6f2c8923e3a1c02d945a75463c9d666d791ec0a077c0a

phase=toolchain
export LANG=C LC_ALL=C TZ=UTC OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1
export PATH=/usr/bin:/bin:/opt/homebrew/bin:/usr/sbin:/sbin
cmake=/opt/homebrew/bin/cmake
[[ -x $cmake && -x /usr/bin/cc && -x /usr/bin/c++ ]]
"$cmake" --version > "$output/cmake.version" 2>&1
/usr/bin/cc --version > "$output/cc.version" 2>&1
/usr/bin/c++ --version > "$output/cxx.version" 2>&1
jq -n --arg cmake "$(sha "$cmake")" --arg cc "$(sha /usr/bin/cc)" \
  --arg cxx "$(sha /usr/bin/c++)" --arg patch "$(sha "$output/postbuffer.patch")" \
  --arg source "$(sha "$output/source.tar")" --arg cadical "$(sha "$output/cadical.tar")" \
  --arg cadiback "$(sha "$output/cadiback.tar")" --arg script "$(sha "$0")" \
  --arg bundle "$(sha "$bundle")" --arg audit "$(sha "$preparation_audit")" \
  --arg terminal "$(sha "$preparation_terminal")" \
  '{schema_version:1,source_archive_sha256:$source,cadical_archive_sha256:$cadical,
    cadiback_archive_sha256:$cadiback,postbuffer_patch_sha256:$patch,
    build_script_sha256:$script,accepted_asset_bundle_sha256:$bundle,
    preparation_audit_sha256:$audit,preparation_terminal_sha256:$terminal,
    toolchain_sha256:{cmake:$cmake,cc:$cc,cxx:$cxx},
    build_profile:"Release",generator:"Unix Makefiles",parallel_jobs:2,
    source_bound_execution_admitted:false}' > "$output/build-inputs.json"

phase=configure
"$cmake" -S "$output/source" -B "$output/build" -G 'Unix Makefiles' \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DSTATIC_BINARY=ON \
  -DENABLE_TESTING=OFF -DNOMPI=ON -DNOBREAKID=ON \
  -DFETCHCONTENT_FULLY_DISCONNECTED=ON \
  -DCMAKE_C_COMPILER=/usr/bin/cc -DCMAKE_CXX_COMPILER=/usr/bin/c++ \
  "-DFETCHCONTENT_SOURCE_DIR_CADICAL=$output/cadical" \
  "-DFETCHCONTENT_SOURCE_DIR_CADIBACK=$output/cadiback" \
  > "$output/configure.stdout" 2> "$output/configure.stderr"
phase=compile
"$cmake" --build "$output/build" --target cryptominisat5-bin --parallel 2 \
  > "$output/build.stdout" 2> "$output/build.stderr"
phase=pin-binary
mkdir "$output/bin"
cp "$output/build/cryptominisat5" "$output/bin/prepared-cms"
chmod 755 "$output/bin/prepared-cms"
jq -n --arg binary "$(sha "$output/bin/prepared-cms")" \
  --arg configure_stdout "$(sha "$output/configure.stdout")" \
  --arg configure_stderr "$(sha "$output/configure.stderr")" \
  --arg build_stdout "$(sha "$output/build.stdout")" \
  --arg build_stderr "$(sha "$output/build.stderr")" \
  '{schema_version:1,status:"BUILT_UNVALIDATED",prepared_cms_sha256:$binary,
    logs_sha256:{configure_stdout:$configure_stdout,configure_stderr:$configure_stderr,
      build_stdout:$build_stdout,build_stderr:$build_stderr},
    transport_parity_passed:false,source_bound_execution_admitted:false}' \
  > "$output/build-receipt.json"
