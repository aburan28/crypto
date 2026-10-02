# Production F4/F5 specialization correctness control

Status: preregistered; execution pending. This is a new, single bounded stage
control following the closed Boolean/full-readback root control in PR #1048.
Neither that run nor the closed 256-query F5 development run is resumed.

Hypothesis: the historical solver's canonical substitutions preserve the
disclosed roots, and its production active-multiplier, linear-tail, flat-matrix
F4/F5 readback returns exactly the decisive constant/forced-variable
consequences of independently enumerated unfiltered Boolean F4 matrices at
the registered prefixes. Those kernels differ from the public full-readback
routine checked previously.

The frozen inputs are the same two already exposed n17a1 queries, trials 164
and 173, with their disclosed points and all six independently derived group
chain assignments. These are deliberately guided correctness controls, not
ordinary-query search, an estimate of natural relation yield, or actual
recursive solver traversal. Original and interleaved layouts give 24 guided
paths. Specialize variables in numerical order through all 35 bits, preserving
all 840 substitution steps. Reduce at depths 0, 6, 12, 18, 24, 30 and 35.
Deduplicate identical systems within each query before reducing them. There
are at most 336 matrix calls, including both F4 and F5. Use degree three and
the historical default production flags with no ambient KIC overrides.

The reference independently constructs the exact S3 Boolean equations and
specializes square-free GF(2) polynomials by set parity. It enumerates every
unfiltered F4 multiple and uses Python integer row reduction to test membership
of the constant and every single-variable assignment polynomial. Compare the
entire forced-variable set, not only rank, row count or witness preservation.
The independently checked model means no registered node may contain a
contradiction. This reference runs after the native producer and supplies no
matrix or consequence to that producer.

Build from the immutable source commit
`765c3c5f19032bd852163805f257c56babef2040`, with root/dependency manifest
`c64e4b3102bface63a2305efbff4bd85810cc112cb43546da2992ba48e9e85b7`.
Copy those sources into the new output directory. Append only the frozen
`kernel_append.rs` diagnostic access wrappers to the copied kernel and add
`export.rs` as a new example. Record the original kernel, exact append bytes,
base manifest, compiled manifest, complete source/dependency archives,
toolchain and interpreter receipts before compilation. No historical source
file, registered candidate ID or old worker is modified. This diagnostic
overlay has its own source identity and no complete-candidate ID.

Use one Rayon thread and one Cargo build job; offline locked release build,
no default features or inherited compiler/solver overrides. Stop after the
first build or native failure, 1,800-second build watchdog, 600-second native
watchdog, or the single registered native execution. Keep owned process-group
cleanup and preserve partial output. No retry, resume, budget extension or
post-outcome selection of paths is permitted. Version probes may fail and
their actual stdout/stderr/exit status are retained; that is not a failed
compilation. This is a physical macOS ARM64 correctness diagnostic. The
unbounded memory envelope and raw times establish no competitive ratio.

Success requires exact substitutions for all paths, exact decisive sets for
both engines at every unique node, unchanged source/input/executable gates,
and independent byte-verified transport replay. Any changed coefficient,
substitution, node input, consequence, binary, source archive or invocation
must reject replay. Preserve negative outcomes, including a control failure,
with null complete-DLP costs. Internal GF(2) PDP matrices remain distinct from
the final subgroup-prime relation LA.

Commit the protocol, inputs, code and hashes before execution and publish
the complete raw inventory, native binary, source archives, process records,
result and independent replay in this PR. The run command is:

```sh
PYTHONDONTWRITEBYTECODE=1 python3.12 research/ic_candidate_tournament_20260915/goal_20260924/f5-production-path-control-20260930/run_control.py \
  --historical-source /Volumes/SSD990/crypto/worktrees/ic-generic-source-pinned-20260929 \
  --out /absolute/new/production-control
```

A pass closes this stage's production-reduction/specialization check only.
Actual solver traversal and a complete source-bound one-target F4/F5 IC solve
remain unestablished. Fresh paired qualification with the incumbent and matched
one-target rho remains pending. All three old confirmation rounds stay closed.
