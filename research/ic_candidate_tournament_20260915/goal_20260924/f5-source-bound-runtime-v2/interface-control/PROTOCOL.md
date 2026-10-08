# Actual native Job/Config deserialization control

Status: preregistered; execution pending. This is a source-bound interface
correctness control after the closed F5 v1 parse failure, not an IC run.
Commit this protocol, exporter, append, runner and eight input cases before
the one native control. Preserve the first build/native failure; no retry,
resume or budget extension. A new IC registration is forbidden until this
control passes and its retained proof and the v2 adapter are accepted.

The hypothesis is that the corrected v2 supplied-point job and its omitted-seed
equivalent deserialize through the actual historical private `Job` and `Config`
types with exactly the declared/default configuration. Six adversarial inputs
must fail: historical `[null]`, negative seed, algorithm-seed overflow, unknown
job field, unknown config field and wrong coordinate shape. These controls
generate no target, run no feasibility oracle and call no `execute_job`,
arithmetic, relation collector or PDP solver. Parsing does not prove a complete
IC implementation or measure natural yield.

Copy the exact source at commit `765c3c5f19032bd852163805f257c56babef2040`,
with full root/dependency manifest
`c64e4b3102bface63a2305efbff4bd85810cc112cb43546da2992ba48e9e85b7`.
Append only `kernel_append.rs` to the copied `examples/ic_tournament_worker.rs`
to expose a parser-only function from that module. Add `export.rs` as a new
example importing the copied module. Preserve the original source, exact append
and complete compiled source/dependency archives. Do not modify the historical
worktree or the admitted IC worker binary.

Use the existing controlled toolchain receipt path: offline, locked release,
no default features, no ambient compiler/solver overrides, one build job and
one Rayon thread. Bound compilation at 1800 seconds and the parser control at
30 seconds, each with owned process-group cleanup. Retain actual tool-version
probe outcomes, including unsupported version flags. Raw process intervals
are uncalibrated diagnostics, never IC or PDP costs. This is a physical macOS
ARM64 control; portable Python replay proves no additional hardware execution.

Success requires both accepted inputs to produce every original/default field
exactly and all six invalid inputs to retain a nonempty rejection error, with
unchanged source, input, executable and command gates. Independent archive
transport must reproduce the audit without running a native solver. Any
changed job, field, compiler overlay, output or binary must reject replay.
Complete DLP costs, online time, rho ratios, S/floor and promotion remain null.

Run once from the accepted repository:

```sh
PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/goal_20260924/f5-source-bound-runtime-v2/interface-control/run_control.py \
  --historical-source /Volumes/SSD990/crypto/worktrees/ic-generic-source-pinned-20260929 \
  --out /absolute/new/job-schema-control
```

The v2 protocol separately registers its new seed and one complete IC attempt
only after these prerequisites. All v1 and historical confirmation registrations
remain closed.
