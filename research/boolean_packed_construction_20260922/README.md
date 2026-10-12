# Packed Boolean construction: fixed candidate and strong control

This standalone bounded experiment continues the parameterized support-envelope
study. Its previous turn produced verified progress: a correct coefficient-aware
contract, a frozen negative timing result, and evidence that application cost,
not only compilation, exceeded direct construction. Those artifacts remain
unchanged. The present objective is a large measured construction improvement;
full solving and cryptanalytic costs are outside this worker and remain unknown.

## Mechanisms and correctness contract

The packed envelope compiles product parity groups to fixed ambient column slots.
For each possible current generator degree, it retains a view of eligible plans
in the same numeric multiplier order as direct construction. Concatenating degree
buckets would reorder the matrix and is deliberately avoided. Runtime evaluates
coefficient parity directly into packed rows, omits zero rows, and ORs surviving
rows into an occupancy bitmap. That bitmap defines the actual ascending column
support. A final remapping packs only occupied columns; when every ambient
column occurs, the packed rows can be returned without remapping.

All previous coefficient, cancellation, degree-drop, zero/constant-generator,
context-guard and direct-fallback rules still apply. Compilation and output caps
are separate. A wide ambient basis does not cause a narrow actual matrix to be
rejected. Row ordering, column ordering and every output bit must agree with the
independent set-parity oracle, including across 64-bit boundaries.

The new **packed direct** control reuses only `(n,D,active)` coordinates and
multiplier lists. It XORs each current term product directly into a row and uses
the same occupancy compactor. It has no generator-envelope restriction, stored
coefficient routing or cached polynomial products. This control distinguishes a
better output representation from an advantage due to symbolic schedule reuse.
The five previous arms remain in the same executable.

## Predeclared gates

`protocol.json` is frozen with the source before timing. The fixed grid includes
n=6/8/10/12, four families, batch sizes 1/4/16/64, two discovery seeds, **new**
holdout seeds 20260925/99083, and fourteen balanced repetitions of seven arms.
No candidate tuning is allowed after these holdouts are observed.

- Correctness: every output equals the independent oracle, and exact expected
  hits/fallbacks are verified. Both packed arms pass the entire 8,192-case small
  coefficient/degree/mask grid, plus word-boundary and current-cap tests.
- Incremental packed-envelope gate: all 48 coefficient/degree-cycle comparisons
  at batches 16 and 64 must have a 95% paired-bootstrap lower bound above 1.05
  against direct, layout and exact-matrix controls.
- **Dramatic construction gate:** each new candidate is evaluated separately.
  At batch 64, for both changing-coefficient families and all four n values, its
  95% paired-bootstrap lower bound must exceed **2.0** against the **pointwise
  minimum of all five legacy arms**, in every one of eight comparisons.
- Attribution: the packed envelope must separately beat packed direct by the
  1.05 lower-bound criterion in all 16 size/family/batch comparisons before a
  gain can be credited specifically to symbolic reuse.

Repeat and support-escape workloads are retained guardrails, not members of the
2x gate. The fixed cold workload includes compilation, all applications and
fallbacks, output allocation, full output validation and destruction. Envelopes,
fixtures and independent reference generation are common supplied-input costs
outside arm timing; fresh-process receipts include them. Retained capacity bytes
exclude allocator metadata; whole-worker RSS is not candidate-specific memory.
Bootstrap intervals describe this finite paired experiment, not a population.

The best correct measured constructor is the reference boundary. Output
materialization remains a necessary work term. There is no calibrated operation
floor or full-pipeline measurement here. A gate pass is an engineering result
for matrix construction, not a result about Gröbner solving or index calculus.
The WDSat/full IC suite is inapplicable to this standalone worker.

## Run and reproduce

```sh
python3 research/boolean_packed_construction_20260922/run.py --out research/boolean_packed_construction_20260922/run_01
python3 research/boolean_packed_construction_20260922/run.py --protocol research/boolean_packed_construction_20260922/confirmation_protocol.json --out research/boolean_packed_construction_20260922/run_02
python3 -m unittest discover -s research/boolean_packed_construction_20260922 -p 'test_*.py'
```

Every output directory is new. Source/protocol copies, executable hashes,
fixtures, raw measurements, process receipts and a file manifest are retained.
Evidence replay tests use temporary copies; frozen runs are never rewritten.

The confirmation protocol uses discovery seeds 211/1009 and holdout seeds
20260926/104729, all disjoint from the first run. The worker source and all
acceptance gates are unchanged. The runner gained only the optional protocol
path argument between runs; each run retains its actual runner source. A
confirmation hash binds its worker to the primary frozen source. This is a
second seed set on the same host, not external reproduction.
