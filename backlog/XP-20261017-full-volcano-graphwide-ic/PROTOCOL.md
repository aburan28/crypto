# Full graph-wide index-calculus campaign

User scope: execute direct elliptic-curve index calculus with factor bases
and finite-field reductions where supported, plus the previously requested
generic-method comparison. Preserve the original fixed A--E requirements.
The three additional completed graphs are an immediately enumerable corpus,
not a substitute for A/C/D/E. “Euclidean” is provisionally interpreted as
“elliptic-curve”; a question about a specific Euclidean-algorithm method is
pending and would change the applicable algorithm, not graph coverage.

Full graph coverage requires every stored vertex, every conductor stratum,
and an exact expected-workload ledger. Stratified sampling, estimates from
other vertices, and historical B samples do not count as full execution.
Run one-target workloads separately; retain failure, timeout, unsupported,
and not-run states. A row without successful verified scalar recovery is not
a completed solve. Distinguish complete coverage from successful recovery.

## Workload preparation (not solver timing)

Read the byte-frozen SQLite graphs from the prior delivery. Verify their
hashes, full-class flags, outgoing-map replay receipts and conductor counts.
Use each class's actual largest prime factor r and cofactor h=N/r, certified
from the native factorization. Deterministically find a nonzero G=[h]P per
vertex and check curve membership and [r]G=O. Do not replace the 45-bit r of
the degree-3 class with a cheaper subgroup. Bind the exact model, vertex,
conductor, r and G into a native canonical workload-input digest.

Compute k=ord_r(p) using a complete factorization of r-1; preserve its prime
divisibility/minimality evidence. The finite-field target is F_(p^k), not F_p
unless k=1. Compare with the actual supported pairing/field solver backend;
an unsupported large k is an implementation/resource status, not a proof
that no pairing reduction exists. No BSGS/rho fallback may be called IC.

## Execution contract (not yet launched)

Use native ecbench. Before launch, register exact model and subgroup/generator
representations, freeze candidate identities and all method parameters, and
pass a control that prevents factor-base points outside the selected subgroup
from being treated as having unqualified logs to G. Explicitly distinguish a
projected subgroup base from a quotient/fiber base and prove/check the
corresponding relation equation. Do not silently project/change coordinates
and reuse an earlier factor-base identity.

For every vertex use 8 independently generated public single-target workloads,
5 interleaved measurement rounds and 1 warmup. Primary per-target comparisons
share exact curve, r, G and Q across all arms. A secondary transported-target
experiment must verify edge maps preserve the subgroup; reusable setup,
transport and recovery remain separately charged. Keep it secondary.

Freeze a factor-base size grid 16,32,64,128; direct relation oracles include
explicit two-summand subtraction and meet-in-the-middle, plus algebraic
relations when implemented. References include rho.negation and BSGS; full
composite-order PH is a separate workload from a largest-prime-subgroup DLP.
Test applicable MOV/Frey--Rueck routes and actual finite-field IC; otherwise
record the exact missing backend/degree/cost and leave that requirement open.
Initial deterministic per-arm relation budget is 1,000,000 trials. Any staged
bounded diagnostic gets its own run id and does not replace this campaign.
Do not launch the full matrix until per-run storage/time estimates come from
qualified controls; those estimates do not count as measured full results.

The baseline is the existing unmodified native IC recipe on these exact
workloads; candidate variations change only declared factor-base/oracle
parameters. Matched rho is the reference. No hypothesis of faster-than-rho
is assumed. Success for measurement is exact coverage with verified answers
or explicit failures, complete phase accounting and replay; improvement
requires lower cold total cost against the verified IC baseline on the same
instances. Stop only per frozen operation/resource limits or correctness
failure; preserve all failed attempts and resume states.

Record factor-base build, decomposition setup, relation collection, matrix
filtering/rank/linear algebra, target descent/recovery, scalar verification,
map transport, failed attempts, memory and disk. Emit trials, valid/duplicate/
dependent relations and independent rank: relation yield is valid relations
divided by attempted decompositions, with numerator/denominator retained.
Separate cold setup from target-dependent online cost; no free reusable logs.
Record all unpriced counters and compare the same 0,0.1,1 GAE pricing band.
Wall time needs L2, A/A and matched controls; unadmitted time remains unknown
as a qualified result. Warm/amortized and cold totals are separate views.

Aggregate min/max/mean by graph, conductor, method, base, target and phase,
with denominators and success/failure/timeout counts. Do not drop censored
or failed rows into a success-only average or infer timings for unrun rows.
Verify every recovered scalar with [x]G=Q and replay all deterministic records.
Preserve exact candidate/workload/run identities, source/binary hashes, host,
isolation and a durable per-vertex completion ledger. Upload the completed
sessions, receipts, full coverage tables and source-linked visual/PDF reports
to the existing Drive study folder; keep original uploads unchanged.

Current status: workload inventory preparation only. No graph-wide solver
timing, relation yield or runtime distribution has yet been measured. The
existing graph archive remains graph data, not an IC benchmark dataset.
Conductor reservation attempted before edits; executable unavailable (127).
