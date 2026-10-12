# Separate F4/F5 history plots

This native, read-only tool prepares graph data from existing `perf.history/v1` snapshots and the frozen Boolean allocation audit. It does not run a benchmark. Python is used only to lay out recorded values and compile PDFs.

```sh
cargo build --locked --release --manifest-path tools/groebner-history/Cargo.toml \
  --target-dir /tmp/groebner-history-target
/tmp/groebner-history-target/release/groebner-history prepare \
  --history research/f4_f5_history_20261012/records/history \
  --groups docs/perf/groebner-groups.json \
  --allocation-index docs/perf/f4-f5-allocation-history.json \
  --root . --out /tmp/f4-f5-plot-data.json
python3 -m venv /tmp/groebner-render
/tmp/groebner-render/bin/pip install -r scripts/perf/groebner-render-requirements.txt
/tmp/groebner-render/bin/python scripts/perf/render_groebner_history.py \
  --data /tmp/f4-f5-plot-data.json --out /tmp/f4-f5-graphs
```

Outputs: separate matrix-step and solver F4/F5 series, instruction and diagnostic-runtime graphs, complete-call and criterion allocation groups, eight SVGs, four group PDFs, a combined PDF, editable Typst sources, a Markdown report and machine-readable provenance.

Every historical source retains its SHA-256, revision, host, compiler, profiler, harness hash and recorded sample count. The producer's UTC-commit filenames establish revision ordering. Missing data remains null. Connections break at missing points, changed fingerprints or changed measurement conditions. Different algorithm fingerprints are expected; equality is checked only along each algorithm's own history.

Allocation import revalidates the original JSONL rows against the frozen native comparison: matched keys, exact output/accounting fields, shared inputs, per-call counts and aggregate totals. Append future accepted audits to `docs/perf/f4-f5-allocation-history.json`. The workload and repetition contract must match; a different contract requires a separate group. A/A and replication runs are controls, not extra production revisions.

The historical snapshots were Python-orchestrated. Their runtime has no stored spread or isolation eligibility and remains diagnostic. The historical profiler scope misses some worker-thread work. No wall-speedup, complete solver operation count, end-to-end ECDLP gain, peak-memory result or new benchmark is claimed.

Validation:

```sh
cargo test --locked --manifest-path tools/groebner-history/Cargo.toml
cargo clippy --locked --manifest-path tools/groebner-history/Cargo.toml \
  --all-targets -- -D warnings
```

The performance-history workflow prepares and renders these outputs beside the existing aggregate history, commits every generated form to `perf-history`, and uploads them together. The frozen initial rendering is in `research/f4_f5_history_20261012/`.
