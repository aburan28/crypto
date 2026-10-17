# End-to-end performance release gate

This is a release requirement, not a measured speedup claim. The initial
hard-panel performance corpus is GF(2); finite-field correctness tests have
broader coverage. Passing this panel would establish a result for these
workloads and this host, rather than a speed advantage on every finite field.

The gate requires all of the following:

* At least 20 distinct, preregistered, coupled nonlinear inputs, whose fastest
  completed F4/F5/F6 baseline takes at least 1 second end to end.
* Genuine, source-identified F4, F5, and F6 implementations executing the same
  all-roots contract. A first-relation search or a basis-only timing is not a
  comparable all-roots baseline. An F4 signature option is not silently F5.
  Multiple executable configurations can share a truthful algorithm identity:
  `F4_msolve` and `F4_native` each declare `algorithm: "F4"`; each worker
  still reports `algorithm: "F4"`. The gate retains per-configuration results
  and uses the faster median for the F4 family. Every registered configuration
  must complete successfully. A timeout cannot be discarded to manufacture a win.
* At least three cold runs per backend and instance (default five), with order
  alternating across rounds. The process interval includes input parsing,
  field construction and verification, domain construction or formulation,
  solving, extraction of every root, exact verification, and serialized output.
* Identical exact root sets, counts, finite-field representation, and variable
  order. Every returned tuple is independently checked against the original
  equations and domain constraint. Every backend must claim complete output;
  truncated output, an exact count without all tuples, or a timeout cannot pass.
* At least 2x geometric-mean speedup against the fastest baseline on each
  eligible instance. No instance may exceed 1.10 times that fastest baseline's
  median time. Stronger thresholds are configurable; release floors cannot be
  weakened to admit a toy benchmark.
* An A/A check of each executable before comparisons: the same command runs
  twice per round on the first input. Maximum pairwise relative spread must
  stay within 10%. A difference obscured by contention or noise is incomplete
  evidence.
* An uncontended completed `isolated_bench.py` record bound to the exact harness
  invocation and its pinned CPU. All runs use one thread. Isolation records and
  raw failed runs remain available for inspection.

`PASS` approves release under the declared benchmark policy. `FAIL` means
complete evidence violates a requirement. `INCOMPLETE` means evidence or a
matched baseline is missing; any observed errors or regressions are also listed.
An unavailable or failed baseline never counts as an infinite speedup.

## Root-report protocol

Each backend command receives the JSON input path through `{instance}`.
`{python}` expands to the harness's Python executable. It must write exactly
one JSON object to stdout, and diagnostics to stderr:

```json
{
  "schema_version": 1,
  "instance_id": "value of ROOTFINDER_INSTANCE_ID",
  "input_sha256": "SHA-256 of the original input bytes",
  "algorithm": "candidate or F4 or F5 or F6",
  "field": {"kind": "prime", "characteristic": 2},
  "variables": ["x", "y"],
  "roots": [[0, 0], [1, 1]],
  "root_count": 2,
  "complete": true,
  "verified": true,
  "pipeline": ["input", "field", "domain", "solve", "extract", "verify", "serialize"]
}
```

`field` and `variables` copy the original input objects exactly. Root
coordinates are canonical nonnegative integers or decimal/hexadecimal strings;
no floating-point conversion is permitted. Extension-field coordinates use
the input's polynomial basis. Duplicate roots fail validation. The protocol's
verification flag is checked alongside independent tuple membership checks;
completeness also depends on comparison with independently implemented
all-root algorithms.

## Running and preserving evidence

Use an otherwise idle host. Reserve a physical core including its SMT sibling,
and pin the solver subprocesses to one of those reserved logical CPUs. For
example, with CPUs 2 and 98 known to be siblings on that host:

```sh
python /workspace/generalized_roots/vendor/crypto/tools/isolated_bench.py reserve \
  --cpus 2,98 --out /workspace/generalized_roots/performance/isolation.jsonl \
  --period 1 --label rootfinder-all-roots \
  -- python /workspace/generalized_roots/performance_gate.py run \
  /workspace/generalized_roots/performance/release_manifest.json \
  --output /workspace/generalized_roots/performance/raw_measurements.json --cpu 2

python /workspace/generalized_roots/performance_gate.py evaluate \
  /workspace/generalized_roots/performance/release_manifest.json \
  /workspace/generalized_roots/performance/raw_measurements.json \
  --isolation /workspace/generalized_roots/performance/isolation.jsonl \
  --output /workspace/generalized_roots/performance/gate_report.json
```

Choose CPUs from the host topology; the example numbers are not portable. Do
not build, test, or run other benchmarks during the reservation. The isolation
tool refuses high pressure or busy hosts. A diagnostic subset may be measured
with a separate manifest; it cannot discharge the 20-instance release gate and
must remain `INCOMPLETE`.

The harness preserves exact commands, JSON stdout, stderr, wall and CPU time,
resource limits, host data, and failures. The per-child `wait4` RSS measurement
is the largest individual process high-water mark, including waited children;
it is not the sum of concurrent descendant RSS. The configured address-space
limit applies per process and is inherited by descendants. Timeout cleanup
terminates the entire subprocess group. Incremental measurement records are
written after each comparison round so interrupted runs remain inspectable.

## Manifest layout

Use `schema_version: 1`, `preregistration` containing `registered_at_utc` and
`selection_rule`, `policy`, `workdir`, `backends`, and `instances`. Every backend
needs a command argument array and a source/algorithm `identity`. An unavailable
backend has `command: null` and an explicit `reason`; keep it in the manifest.
Every instance needs `id`, a path relative to the manifest, original-byte
`sha256`, and `hard: true` if preregistered for the hard panel. Input content
hashes are checked before every run. Duplicate canonical JSON inputs cannot
inflate panel size. Resource budgets and identities must be frozen before
comparison; calibration failures are retained rather than removed after seeing
which algorithm wins.

The unit tests use explicitly synthetic timings to test gate logic. Their
passing results are not performance measurements.

The release and earlier bounded manifests remain as preregistration records.
Creating a manifest does not execute it. The combined bounded panel is
`full_bounded_combined_manifest.json`; its raw measurements and isolation
records, if present, establish which experiment actually ran. Recorded source,
binary, and harness hashes are checked before launching that experiment.
