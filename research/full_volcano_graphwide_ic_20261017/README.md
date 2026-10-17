# Full graph-wide index-calculus execution

The full-graph IC requirement is preserved in the committed [protocol](PROTOCOL.md).
It includes every vertex and conductor, direct EC relations and finite-field
reductions, factor-base grids, matched generic references, complete setup and
online accounting, raw yields, failures and min/max/mean coverage statistics.
The original A--E requirements remain open; the three additional complete
graphs provide an immediately available corpus.
The [complete graph explorer](graph-viewer/README.md) embeds every vertex and
recorded selected-degree map, with conductor filters and exact kernel inspection.

Native input certification has completed over **every one of the 24,823
vertices and 51 strata**. This is workload preparation, not a solver benchmark.

**Execution update:** [native control and pilot measurements](execution/README.md)
now exist: 60 attempts total, 38 verified, 18 exhausted, 4 timed out. All
successful and exhausted deterministic reports reproduce on replay. A separate
checkpointed bounded survey has a stopped, audited checkpoint of **99 vertices
across all 51 strata**, with 1,584 attempts (1,065 verified, 422 exhausted,
97 timed out). [Download the complete measured checkpoint and graphs](https://drive.google.com/file/d/1hLVWVad2SrgHaBn6U5iSyqnKbx5cAkeY/view).
The separate optimized continuation has frozen the complete 24,823-vertex queue
after 56 exact completed deterministic profile comparisons, with four censored
pairs excluded. Consult its launch receipt and ledger for execution status;
planning alone is not execution. The original full matrix and A--E remain open. The preparation
archive below is an unchanged historical NOT_RUN snapshot.

[Download the certified inputs and preparation report](https://drive.google.com/file/d/1SPpdcNjm7q6ChcMkaty0NgED8lyHofHq/view?usp=drivesdk).
The filename explicitly states NOT_RUN. Local ZIP CRC passed; Drive metadata
readback verified the archive name, 3,210,305-byte size and parent folder.
Local SHA-256: `8acb83dca7a6bca8736429e515e6e31a141d3dbadba8d40438ac707be8bf8f37`.
No remote whole-file hash was checked.

| Source model | Vertices | Strata | Largest prime subgroup r | Exact embedding degree |
| --- | ---: | ---: | ---: | ---: |
| icv1-fp57-t536870914-6cabb0c9 | 8192 | 14 | 61751 | 30875 |
| icv1-fp57-tm536870941-e6d0c766 | 6561 | 9 | 19768889568671 | 3953777913734 |
| icv1-fp59-t1073741954-4dd0b1ef | 10070 | 28 | 4253 | 2126 |

The native inventory checks original database SHA-256 values, complete graph
receipts, every conductor population, curve membership, nonzero subgroup
generators and [r]G=O. Complete native factorizations certify r and k=ord_r(p).
It records exact explicit ecbench curve inputs and a canonical digest per
vertex. In that frozen preparation snapshot, solver states are NOT_RUN;
timings, relation counts and total solver costs are null. Current measurements
are separate. Representation registration is explicitly pending.
No smaller subgroup replaced the actual 45-bit largest prime factor of the
degree-3 class. A largest-prime-subgroup experiment does not itself measure
full composite-order PH or all prime-power components.

Finite-field reductions must target F_(p^k). The existing supersingular k=2
MOV routine does not implement these ordinary-field destinations and uses
BSGS/PH rather than finite-field IC. That is a concrete backend gap, not a
proof that no reduction exists. The direct-EC path now has a checked subgroup
base, audited controls and a resumable bounded survey; full qualification is open.
See [the shared open problem](../../problems/full-volcano-class-storage-and-execution.md).

Preparation source is [inventory.rs](inventory.rs); command:

```sh
volcano-ic-inventory /path/to/crypto /path/to/graphs /fresh/output
```

It is linked to the repository library with the retained native command in
`BUILD.sh`; no Python arithmetic or research orchestration is used. Run from
the pinned code/library revision; no measured solver is invoked. Output has
one JSONL row per vertex plus a compact native inventory, hashes, byte counts
and per-conductor populations. Source commit `e1b0e302e` preceded execution.

The first link exhausted the 8.8 GiB temporary filesystem and aborted, and
the missing executable launch exited 127. Empty build output, disk evidence
and the failed launch were retained. Re-linking with symbol stripping and a
workspace temporary directory succeeded. All prior datasets were preserved.

Next acceptance gate: registered exact representations and an audited native
control/pilot with subgroup-correct relations, complete cold costs and matched
references; then cost/provision and execute the complete frozen matrix.
The original full matrix has not launched. The dev bounded survey is paused;
optimized continuation and current coverage are recorded separately. The workload
archive is labeled inputs/NOT_RUN and is not a set of benchmark results.
