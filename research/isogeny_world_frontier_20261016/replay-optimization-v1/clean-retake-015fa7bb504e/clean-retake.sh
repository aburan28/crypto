#!/bin/sh
set -eu
repo=/Volumes/SSD990-2/crypto/isogeny-degree-challenge-20261015
build_root=/Volumes/SSD990-2/crypto/isogeny-fixed-modulus-builds-20261016/clean-retake-v1
source_commit=$1
cd "$repo"
[ "$(git rev-parse HEAD)" = "$source_commit" ]
[ -z "$(git status --porcelain=v1)" ]
mkdir "$build_root"
printf '%s\n' "$source_commit" > "$build_root/source-commit.txt"
git rev-parse 'HEAD^{tree}' > "$build_root/source-tree.txt"
rustc --version --verbose > "$build_root/rustc.txt"
rustc --print cfg > "$build_root/compiler-cfg.txt"
rustc --print target-features > "$build_root/compiler-target-features.txt"
sw_vers > "$build_root/os.txt"
sysctl hw.model hw.ncpu hw.memsize machdep.cpu.brand_string hw.optional.arm.FEAT_AES hw.optional.arm.FEAT_PMULL > "$build_root/host.txt"
printf 'RUSTFLAGS=%s\nCARGO_ENCODED_RUSTFLAGS=%s\nCARGO_BUILD_TARGET=%s\nRAYON_NUM_THREADS=1\n' "${RUSTFLAGS-<unset>}" "${CARGO_ENCODED_RUSTFLAGS-<unset>}" "${CARGO_BUILD_TARGET-<unset>}" > "$build_root/flags.txt"
run_arm() {
  arm=$1
  shift
  arm_root="$build_root/$arm"
  mkdir "$arm_root" "$arm_root/install" "$arm_root/controls"
  git status --porcelain=v1 > "$arm_root/git-status-before.txt"
  [ ! -s "$arm_root/git-status-before.txt" ]
  [ "$(git rev-parse HEAD)" = "$source_commit" ]
  printf '%s\n' "$@" > "$arm_root/feature-arguments.txt"
  printf '%s\n' 'cargo build --release --locked --manifest-path tools/isogeny-cli/Cargo.toml --bins' > "$arm_root/build-command.txt"
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$arm_root/build-started.txt"
  set +e
  /usr/bin/time -l -o "$arm_root/build-os-time.txt" env RAYON_NUM_THREADS=1 CARGO_TARGET_DIR="$arm_root/target" cargo build --release --locked --manifest-path tools/isogeny-cli/Cargo.toml --bins "$@" > "$arm_root/build.stdout" 2> "$arm_root/build.stderr"
  status=$?
  set -e
  printf '%s\n' "$status" > "$arm_root/build.exit"
  [ "$status" = 0 ]
  date -u '+%Y-%m-%dT%H:%M:%SZ' > "$arm_root/build-ended.txt"
  for binary in isogeny fixed-modulus-bench; do
    cp "$arm_root/target/release/$binary" "$arm_root/install/$binary"
    chmod 555 "$arm_root/install/$binary"
    cmp "$arm_root/target/release/$binary" "$arm_root/install/$binary"
    "$arm_root/install/$binary" version > "$arm_root/$binary-version.json"
    shasum -a 256 "$arm_root/target/release/$binary" "$arm_root/install/$binary" > "$arm_root/$binary-sha256.txt"
  done
  # Every test output is outside the source worktree, so both builds see clean Git.
  set +e
  /usr/bin/time -l -o "$arm_root/tests-os-time.txt" env RAYON_NUM_THREADS=1 CARGO_TARGET_DIR="$arm_root/target" ISOGENY_CHALLENGE_BINARY="$arm_root/install/isogeny" ISOGENY_CONTROL_BINARY="$arm_root/install/isogeny" ISOGENY_CHALLENGE_EVIDENCE="$arm_root/controls/challenge" ISOGENY_TRACE_EVIDENCE="$arm_root/controls/trace" ISOGENY_FIXED_MODULUS_EVIDENCE="$arm_root/controls/cache" ISOGENY_COMPANION_EVIDENCE="$arm_root/controls/companion" ISOGENY_NEWTON_IDENTITY_OUTPUT="$arm_root/controls/newton-identity.json" cargo test --release --locked --manifest-path tools/isogeny-cli/Cargo.toml "$@" -- --test-threads=1 --nocapture > "$arm_root/tests.stdout" 2> "$arm_root/tests.stderr"
  status=$?
  set -e
  printf '%s\n' "$status" > "$arm_root/tests.exit"
  [ "$status" = 0 ]
  for binary in isogeny fixed-modulus-bench; do
    cmp "$arm_root/target/release/$binary" "$arm_root/install/$binary"
    set +e
    env RAYON_NUM_THREADS=1 CARGO_TARGET_DIR="$arm_root/target" ISOGENY_MONITOR_EVIDENCE="$arm_root/controls/monitor-$binary" cargo test --release --locked --manifest-path tools/isogeny-cli/Cargo.toml --bin "$binary" "$@" supervisor::tests::missing_and_zero_rss_fail_closed_before_and_after_acknowledgement -- --exact --nocapture > "$arm_root/monitor-$binary.stdout" 2> "$arm_root/monitor-$binary.stderr"
    monitor_status=$?
    env RAYON_NUM_THREADS=1 CARGO_TARGET_DIR="$arm_root/target" ISOGENY_TIMEOUT_EVIDENCE="$arm_root/controls/timeout-$binary" cargo test --release --locked --manifest-path tools/isogeny-cli/Cargo.toml --bin "$binary" "$@" companion_shape_retains_deadline_after_positive_acknowledgement -- --nocapture > "$arm_root/timeout-$binary.stdout" 2> "$arm_root/timeout-$binary.stderr"
    timeout_status=$?
    set -e
    printf '%s\n' "$monitor_status" > "$arm_root/monitor-$binary.exit"
    printf '%s\n' "$timeout_status" > "$arm_root/timeout-$binary.exit"
    [ "$monitor_status" = 0 ] && [ "$timeout_status" = 0 ]
  done
  git status --porcelain=v1 > "$arm_root/git-status-after.txt"
  [ ! -s "$arm_root/git-status-after.txt" ]
}
run_arm reference --no-default-features
run_arm candidate --no-default-features --features cached-fixed-modulus-reduction
printf 'completed\n' > "$build_root/completed.txt"
