# P-256 projective pair-table locality, round 35: protocol

Date registered: 2026-10-06

## Requested scope and evidence plan

| continuing requirement | round-35 evidence |
|:--|:--|
| pursue parity with rho | attempt to close round 34's measured cache-to-DRAM gap without removing any lookup, addition or branch from the round-33 selector |
| retain all 17 variable columns | import cutoff 219, eight signed-pair entries and the remaining seventeenth point unchanged |
| propagate structure before materialization | compare software-prefetched random access with a global radix-sorted request stream over whole batches |
| exact group replay | require every locality variant to reproduce the same P-256 group sum for every segment |
| preserve promotion gates | keep a DRAM working-set result separate from a fully materialized 3.318-TB table and from end-to-end relation coverage |

## Question

Round 34 found that directly addable 96-byte projective entries fit the
`0.334410508423`-addition segment budget only in a 96-KiB cache row.  At a
384-MiB working set, eight random lookups add `3.612595` group-addition
equivalents and move the stage to `1.014028` rho.  Can latency hiding or global
request ordering reduce that exact eight-entry access cost below the frozen
budget while preserving every segment result?

This is a locality experiment on the registered round-33 pair-table selector.
It does not change or prove the occurrence-based P-256 coverage bound.

## Frozen dependency and work unit

Import the canonical round-34 artifact by SHA-256
`df68204433a4abf408ac2c39f86d801da7db296056f030360fae8f719bca6a4b`.
Require:

- projective entry size 96 bytes;
- complete pair entries `34,562,148,612`;
- complete projected table `3,317,966,266,752` bytes;
- round-33 corrected ratio `0.998569034150286`;
- mean segment capacity `224.361249307511`;
- admissible overhead `0.334410508423` additions per segment.

One segment begins at the seventeenth point and adds exactly eight signed-pair
points.  The timing denominator is the median time for those same eight
complete P-256 additions with operands supplied from a fixed L1-resident
array.  Index generation is outside every timed variant, as in round 34.

## Registered locality variants

Evaluate the following exact variants:

1. `random`: the round-34 random table lookup and immediate addition control.
2. `prefetch-d`: before processing segment `s`, issue read prefetches for both
   cache lines of all eight entries of segment `s+d`, then process segment
   `s`.  Sweep `d in {1,2,4,8,16,32,64}`.  Charge every prefetch instruction;
   the first and last `d` segments are included, not discarded.
3. `radix-global`: form all eight `(table_index, segment_id)` requests for a
   batch, perform an exact stable two-pass 11-bit LSD radix sort over the
   22-bit measured table index, read projective entries in sorted order and
   scatter-add them into per-segment accumulators.  Charge request filling,
   counter clearing, both radix passes, accumulator initialization, sorted
   reads and scatter additions.  Buffers may be allocated once outside timed
   repetitions, but no filled or sorted state may be reused.

For radix ordering sweep batch sizes `2^8`, `2^10`, `2^12`, `2^14` and
`2^16` segments.  For random and prefetch variants use `2^16` segments.  Use
table working sets of `2^15`, `2^20` and `2^22` entries (3 MiB, 96 MiB and
384 MiB).  The `2^22` row is the registered decision row.  Smaller rows expose
cache transitions only.

If a platform lacks the x86 read-prefetch intrinsic, mark prefetch variants
unsupported rather than replacing them with an uncharged abstraction.

## Deterministic construction and timing

Generate 4,096 hash-selected nonzero scalar multiples of the P-256 generator
and repeat them to fill each physical working set.  Generate request indices
as independent SHA-256 domain-separated words masked to the table size.  Write
every table page before measurement.

Use one release binary and a single OS thread.  Run one untimed warmup and nine
timed repetitions.  Rotate variant order deterministically per repetition.
Record wall time, process CPU time, segment count, checksum, request/index
digests, bytes addressed, compiler, target and CPU model.  Use sorted ranks 1,
4 and 7 of nine repetitions as the registered interval.  For a variant with
layout time `T` and direct-eight-addition time `A`, report

```text
extra addition equivalents = 8*(T/A - 1), clamped at zero.
projected / rho = 0.998569034150286
                + 0.964336477130181*extra/(224.361249307511+1).
```

Select the smallest median projected ratio on the `2^22` decision row;
smaller prefetch distance, then smaller radix batch, breaks exact ties.  The
timing gate uses the selected interval's upper endpoint, not its median.

## Correctness controls

- On 4,096 deterministic segments, compare random, every prefetch distance
  and every radix batch against an independently accumulated reference point.
- Require exactly eight requests and eight additions per segment.
- Require radix output to be a permutation of the input request multiset and
  every segment to receive exactly eight entries.
- Record zero false positives, false negatives, missing requests, duplicate
  requests, out-of-range accesses and group mismatches.
- Hash the point pool, index stream, sorted request stream and canonical group
  outputs.  Repeated release runs must reproduce all semantic hashes; raw host
  timing rows remain volatile evidence.

## Gates and stop condition

A locality variant advances the projective stage only if its `2^22` row has:

- zero correctness and request-accounting failures;
- median and registered upper projected ratios no greater than one;
- all prefetch/sort/scatter work charged; and
- no claim that a 384-MiB resident table proves a 3.318-TB deployment.

Even a passing locality row is not end-to-end promotion.  Complete-time parity
still requires a materialized or defensible full-table hardware model; the
original gates still require proved P-256 relation probability, structured
residual degree at most 5, relation collection below `2^120`, cost per usable
row below `2^103`, storage below `2^50`, and a non-generic algorithm.  Attempt
no full-depth relation unless all gates pass.  If every `2^22` row fails,
reject this projective locality implementation and report the smallest
measured overhead without generalizing to all possible memory systems.
