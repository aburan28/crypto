#!/usr/bin/env bash
# Thin execution sequencing only. Native Rust performs every evidence check.
# Run from the repository root; the first argument is the active census PID.
set -u
study=research/ecbench_conductor_gap_ir_20261014
archive=research/ecbench_conductor_gap_20261012
if [ "$#" -ne 1 ]; then
  printf 'usage: bash %s CENSUS_PID\n' "$0" >&2
  exit 2
fi
tail --pid="$1" -n 0 -f "$study/census-v2.log" > /dev/null
target/release/examples/conductor_gap_ir analyze "$study" "$archive" \
  > "$study/native-validation.log" 2>&1
native_status=$?
printf 'native_analysis_exit=%s\n' "$native_status" > "$study/finish-local.status"
# This workspace uses cone-mode sparse checkout. The CI helper's leading-slash
# patterns require non-cone mode, so use its dependency graph with native dirs.
node --input-type=module > "$study/materialize-final-v2.log" 2>&1 <<'JS'
import {existsSync} from 'node:fs';
import path from 'node:path';
import {execFileSync} from 'node:child_process';
import {sourceFiles, embeddedInputs} from './.github/scripts/ci-scope.mjs';
const cwd=process.cwd();
const sources=sourceFiles(cwd);
const missing=[...new Set([...sources,...embeddedInputs(cwd,sources)])]
  .filter(file=>!existsSync(path.join(cwd,file)));
const dirs=[...new Set(missing.map(file=>path.posix.dirname(file)))];
if(dirs.includes('.')) throw new Error('unexpected missing root compile input');
if(dirs.length) execFileSync('git',['sparse-checkout','add','--skip-checks',...dirs],{stdio:'inherit'});
console.log(`materialized ${missing.length} missing source/input files in ${dirs.length} directories`);
JS
materialize_status=$?
printf 'materialize_exit=%s\n' "$materialize_status" >> "$study/finish-local.status"
if [ "$materialize_status" -ne 0 ]; then
  exit "$materialize_status"
fi
# Keep release opt-level=3 and all default features/tests. More codegen units
# bound monolithic compilation work; this executable supplies no measurements.
CARGO_HOME=/workspace/toolchains/cargo \
RUSTUP_HOME=/workspace/toolchains/rustup \
CARGO_TARGET_DIR=/workspace/toolchains/volcano-ir/test-target \
CARGO_PROFILE_RELEASE_CODEGEN_UNITS=256 \
  python3 tools/isolated_bench.py busy -- \
  /workspace/toolchains/cargo/bin/cargo test --locked --release --lib -j3 \
  -- --test-threads=3 > "$study/library-tests-final.log" 2>&1
library_status=$?
printf 'library_tests_exit=%s\n' "$library_status" >> "$study/finish-local.status"
printf 'native_analysis_exit=%s library_tests_exit=%s\n' "$native_status" "$library_status"
if [ "$native_status" -ne 0 ]; then
  exit "$native_status"
fi
exit "$library_status"
