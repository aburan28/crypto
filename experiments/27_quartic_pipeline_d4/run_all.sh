#!/usr/bin/env bash
# E19 (RESEARCH_GLV_INVARIANT_FACTOR_BASES.md §9.5): one curve per
# invocation, in registered order; a curve whose JSON exists is skipped, so
# the script resumes after a restart of the host.
set -u
cd "$(dirname "$0")/../.."
out=experiments/27_quartic_pipeline_d4
for p in 53 101 149 197 269; do
  for s in 1 2; do
    f=$out/p${p}s${s}.json
    [ -s "$f" ] && continue
    ./target/release/examples/quartic_pipeline_d4 --sizes "$p" --first-seed "$s" --seeds 1 \
      --rho-runs 16 --json "$f.part" >>"$out/run.log" 2>&1 && mv "$f.part" "$f"
  done
done
echo "finished $(date -u +%FT%TZ)" >>"$out/run.log"
