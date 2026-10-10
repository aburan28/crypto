# Two-hour-total N83 pilot extension

The owner extended the local pilot to two hours on 2026-10-10. This protocol
interprets that as 7,200 seconds total, including the frozen prior charge of
3,593.092994294 seconds. The extension has at most 3,606.907005706 seconds;
worker rebuilds and retained runs are charged. No cloud workers are allocated.

The full requested objective remains a factor base minimizing complete,
single-target cold index-calculus runtime across the declared combinations.
This bounded extension first tests a prerequisite, not that objective itself.

| Requirement | Authoritative evidence | Status before launch |
| --- | --- | --- |
| All declared combinations | Frozen v1/v2 ordinal grid | Addressable; mostly unexecuted |
| Stored, exact factor bases | Panel replay and S3 download receipts | 54 v1 objects verified; v2 pending |
| Executable factored-S4 model for a retained base | Worker counts and original input receipts | Unmeasured |
| Best complete cold runtime | Matched, complete, verified single-target runs | Unmeasured |

## First gate

Hypothesis: the current factored S4 builder can construct the exact primary
N83 K=64 finite-domain system inside a 4 GiB, zero-swap Docker cgroup and a
600-second process-wall cap. The reference is the current committed Rust
builder, whose source-derived algebraic variable/clause bounds are documented
in `SOLVER_GATES.md`. This is a construction-only test; it performs no solver
search or lifting and supplies no natural relation, rank or total runtime.

Frozen input: the generically replayed and S3-verified `pilot-01` panel, its
public fixture corpus, and the primary K=64 hash-policy base selected by the
existing `load_panel` contract. Preserve the panel and every historical run.
Freeze HEAD, source hashes, Cargo.lock, compiler identity, binary SHA-256,
image ID, exact commands and container state before interpreting the result.

`run-capacity.sh` is thin shell orchestration of the native Rust worker. It
does not use a Python execution path. The container uses one pinned Docker
CPU, no network, read-only inputs/root, a 4 GiB memory/swap ceiling, a PID cap,
and an immutable locally present image. The Docker VM and shared host are not
a dedicated L2 benchmark; elapsed times are budget/informational L0 evidence.
No timing ratio, throughput comparison or speedup is admitted from this gate.

PASS requires worker exit zero, the exact construction-only schema, K=64,
matching memory limits, positive model counts, and false search/lifting flags.
Timeout or a resource exit is censored UNKNOWN, never UNSAT. Preserve stderr,
stdout, worker JSON, container state, launch/wait/cleanup status and outer JSON
for every launched outcome. A failed source-matched build is a producer
failure and remains charged, not an observation about decomposition existence.

## Conditional next gate

After the capacity result, and only with sufficient remaining allowance,
construct one declared signed-Frobenius v2 K=1,182 public-x-hash object on
the primary curve, seed 2026100801, under the same memory/no-network guard.
Use a fresh directory and separate bounded generic replay. Publication may
use only the existing uploader after the matching replay passes, into the
versioned S3 destination already declared in `SIZE_FRONTIER_V2.md`.
Construction/replay alone does not rank that base. If insufficient time
remains, retain the missing phase as pending and stop at the total cap.

Stop immediately on input/source mutation, cleanup uncertainty or the total
allowance. Preserve negative and inconclusive results. The runtime winner
and all incomplete end-to-end measurements remain null. Canonical runtime
graphs/leaderboards remain unchanged unless a graphable admitted runtime
result exists; a source-linked visual receipt will accompany this gate.
