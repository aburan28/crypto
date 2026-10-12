# F5 coverage and budget accounting review

Classification: **accounting**. This is a source review and a pending design,
with `candidate_id: null` and every proposed measured cost null. No new solver
job, target or qualification registration executes here. The accepted
[terminal panel](results-20260929/README.md) remains unchanged and every
comparison speedup stays unknown.

## What the registered budget actually limited

The measured worker is bound to source commit
`765c3c5f19032bd852163805f257c56babef2040`. Its `SolveStats` retains separate
`reductions`, `splits`, `propagations`, `exhausted`, `unsupported`,
`max_degree_built` and `oversize` fields. The solver checks
`stats.reductions >= opts.node_budget`, both on entering an unfinished branch
and after a propagation/reduction iteration. Cached reductions also increment
the call counter; it is not simply a count of uncached matrix builds.

In contrast, the symmetrised plugin's documentation, parameter description and
human description called `node_budget` a split limit, and the CLI inventory
repeated that description. This PR corrects those descriptions to **algebraic
reduction calls**. It keeps the parameter name and solver behavior. Earlier
candidate manifests, binaries, receipts and measured limits are not rewritten.

The discrepancy exists in the reviewed historical source and the current
implementation. Verify the enforcement in
[koblitz_groebner.rs at the measured commit](https://github.com/aburan28/crypto/blob/765c3c5f19032bd852163805f257c56babef2040/src/cryptanalysis/koblitz_groebner.rs#L4313),
the [separate counters](https://github.com/aburan28/crypto/blob/765c3c5f19032bd852163805f257c56babef2040/src/cryptanalysis/koblitz_groebner.rs#L3665)
and [cached reduction accounting](https://github.com/aburan28/crypto/blob/765c3c5f19032bd852163805f257c56babef2040/src/cryptanalysis/koblitz_groebner.rs#L3696).
The geometric frontend reports `effort = stats.splits`, which is a different
counter from the enforced budget. Do not relabel that effort as reductions or
interpret `node_budget=4096` as 4096 visited branches.

Historical source bytes checked during this review:

| Source role | SHA-256 |
| --- | --- |
| `src/cryptanalysis/koblitz_groebner.rs` | `056fa31dc32f203e9efaa8dafc7672a632e13efab21340467ef98049ac81d89d` |
| `src/cryptanalysis/matrix_f5_f2.rs` | `4615c6c6d04342209002b571cb8ff4ff8ee0574d5547bd1e750a995304ef4cb1` |
| `src/cryptanalysis/koblitz_symmetrised.rs` | `383959eb66610db54893735f8f2274522d8611d60cbc87ad121f38c2ecd31974` |
| `src/cryptanalysis/ic_framework/plugins.rs` | `1dd4e87f25edd97c8211f68d9c4520503eedbcf8dd9885d3fffa9f4d0d51ee1d` |
| `src/bin/ic/bench.rs` | `12ec55b67aa19660e797199ec7737e2d90c4f6ba9019ca4e4c8e0c7c1e6cf3db` |

The current description correction has different source hashes. These listed
hashes identify the measured source under review, not a newly built worker.

## What the matrix evidence does and does not establish

Independent scalar-field elimination of the terminal F5 relations gives rank
28/29, with the nullspace unit vector at column 27. Two retained SAT witnesses
cover that direction on queries that F5 left incomplete. This demonstrates
lost PDP coverage before final relation LA. It does not establish that the
Boolean Macaulay row filter, row-reduction arithmetic or final LA is incorrect.

The source uses a fixed-degree matrix engine with an F5 row criterion. It is
not a full incremental Gröbner-basis algorithm. Its default resolved variable
rule for MatrixF5 is lowest-free, and its recursive branch values are false
then true. A bounded traversal can therefore have structured coverage; that
is a hypothesis to test, not a causal result derived from two missing rows.
Degree settings and actual built degree can also differ. Record observed
`max_degree_built` and oversize events rather than inferring them from a label.

## Pending development gates before another fresh panel

1. Extend the next versioned producer receipt to retain reductions, splits,
   propagations, exhausted/unsupported status, actual built degree, oversize
   events, effective variable and value order, and row-construction/reduction
   word counts per ordinary and target query. Keep the existing effort count
   with its definition. Cache hits and misses need their own counts.
2. Check the algebraic encoding independently against disclosed full-point
   decompositions before attributing lost witnesses to traversal. Include
   repeated points, sign choices, equal abscissae and identity intermediate
   states. A valid wide-S4 witness does not alone prove a chained encoding has
   represented every case. Then cross-check the F5-filtered matrix consequences against unfiltered F4 on
   identical disclosed Boolean systems, including forced assignments and
   contradictions. Use independent equation evaluation and full-point witness
   replay. Charge criterion construction and every matrix operation; internal
   GF(2) elimination remains PDP work, separate from final LA modulo 65587.
3. Register a bounded **disclosed-input** development study before executing it.
   Compare a small predeclared set of reduction-budget and variable/value-order
   policies on ordinary controls, including failed SAT/F5 cases. Declare the
   same source, system/base construction, input order, resources and stop rules.
   Old exposed queries are development controls only. Do not redispatch their
   old registrations or count this study as fresh qualification.
4. Retain every failed, timed-out, unsupported and zero-yield control. Report
   useful rank per ordinary query and cost per useful row as diagnostics, with
   uncertainty appropriate to the fixed input set. Select a complete pipeline
   policy using development evidence; do not add cheap phases from separate
   candidates or choose a policy on confirmation targets.
5. Require full-rank factor logs, verified target descent and independent scalar
   replay in bounded complete-DLP development controls before generating a new
   qualification target. Freeze the chosen policy and all actual source,
   candidate/workload identities and resource limits for the next panel. Exclude
   every previously exposed point, preserve the three sealed rounds, and use
   the source-complete SAT integration plus strong same-point rho on calibrated
   Linux. Missing completions or costs block a speedup claim.

No numeric new budget, branch-order winner, fresh target, candidate ID or
benchmark outcome is selected by this review. Those are subsequent registered
decisions. The persistent goal remains active.

## Validation scope

The correction changes descriptions and documentation. Source review verifies
the two budget guards and cached/uncached counter handling; the final diff
must preserve the solver functions and numeric defaults. Rust formatting and
`git diff --check` are the local checks. No new timing, kernel-performance or
complete-solve test is justified by a text correction. Hosted applicable CI
still gates the PR merge.
