# Speed round, 2026-10-11 (Conductor T-293)

Wall-time speedups of the pair-table family, and the multiplier-closure
prediction tested at a landed rung. Operation counts are unchanged by
everything in the first half: the speedups are implementation, not
algorithm, and are invisible to the ledger's unit by construction.

**Host and method.** Apple M4 Pro, load average 60–130 from other sessions
during every run, so every wall-time figure below is the **minimum** (or
paired median) over interleaved repeats in one process, never a single
timing, and none is an isolated-host receipt (`ICMS` level L0). Where two
binaries are compared they run alternately, three rounds, on the same
sweep and seed.

## 1. Batched Frobenius pair table and three-summand probe

`FrobeniusPairTable::build` formed one scalar affine addition (one field
inversion) per pair, and `decompose_mitm_frobenius` at `m = 3` did the same
for every base point of every trial, although `FastCurve::add_many` (one
shared inversion, Montgomery's trick, four interleaved chains) already
served other paths. Now:

- **build:** each representative's pair sums are one `add_many`;
- **three-summand probe** (`FrobeniusPairTable::probe_three`): blocks of
  `target − P_i` by `add_many`, abscissae canonicalised side by side, only
  table hits resolved. The block length starts at 16 and doubles to 2048: a
  fixed 2048 block made the probe **0.83×** (slower) on a base where 199 of
  200 targets decompose within the first few points, because it computed the
  whole base each time; the geometric schedule wastes at most about twice
  the points a scan needs and still reaches the full batch rate on long
  scans.

**Equivalence is a test, not a claim.**
`batched_frobenius_table_and_probe_equal_the_scalar_reference` holds the new
code to the scalar references (`build_scalar`,
`decompose_mitm_frobenius_scalar`, both kept): every table entry, all build
counters and charged operations, and for random and planted targets at
`m = 2, 3` the decomposition, charged operations, lookups and
canonicalisations, at block schedules (1,1) (3,3) (7,7) (64,64) (1,64)
(4,4096). A first version of the test found nothing to decompose at
`m = 3` and the test refused to pass (`hits > 0` guard); targets are now
sums of base points, which also places the first success at varied scan
offsets. Charging follows the scalar loop: additions and probes up to the
first success, not for the rest of the block.

### Micro-benchmark (`examples/ic_batch_bench.rs`, minimum of repeats, interleaved)

| rung, base | table build | 3-summand probe | identical |
|---|---:|---:|---|
| n = 23 `K_0`, subspace 645 pts | 5.7× | 4.1× (399 / 400 targets decompose early) | 4,412 entries, 5,475 charged adds |
| n = 41 `K_0`, 60 plain columns (4,920 pts) | 8.4× | 13.3× (min), 9.4× (median) | 147,420 entries, 150,060 adds |
| n = 41 `K_0`, 6 columns norm 16 (4,920 pts) | 8.4× | 12.7× | 141,780 entries |
| n = 41 `K_0`, **K = 600** (49,200 pts) | **6.9×** (9.6 s → 1.4 s) | — | **14,742,146 entries**, 14,784,600 adds |

### End to end through `ic bench` (n = 41, `m = 3`, old and new binary alternating, 3 rounds, 4 repeats)

Counted operations, trials, relations, lookups, canonicalisations, rank and
verification are **identical across the two binaries on all 24 paired
rows**. Wall time:

| base | pair-table build | relation phase | whole pipeline |
|---|---:|---:|---:|
| 6 columns, norm 16 | 8.9× | 12.2× | **10.0×** (0.29 s → 0.03 s) |
| 60 plain columns | 9.1× | 11.4× | **11.2×** (3.26 s → 0.29 s) |

### Where the batched build's time goes now (K = 600, 29.5 M pairs timed, three runs)

| part | ns per pair | share |
|---|---:|---:|
| `add_many` | 22–24 | 19 % |
| `canon_with_shift` (scalar) | 18.5–20 | 16 % |
| hash-map insert | 74–86 | **64–66 %** |

The insert figure is from the harness's `std::collections::HashMap` (SipHash)
and so overstates the project's hasher; the cache-missing random writes into
a 15-million-entry table are the cost either way. The sampling profile of the
benchmark agrees (the batched build function's own time is six times
`add_many` plus `batch_inv`; `Gf2::inv` at 67 % of samples is the scalar
reference arm).

## 2. The multiplier closure tested at a landed rung (n = 53, `K_0`, 52-bit r)

`ic bench`, `mitm-frobenius:m=3`, batched table, 2 repeats, all verified.
600 plain columns is the ledger's `K`; 60 columns with norm ≤ 16
multipliers gives the **same 63,600 points**, hence the same table and
yield:

| base | points | columns | trials per relation | relations needed | probe additions | table additions | relation-phase wall | whole pipeline |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 600 plain | 63,600 | 600 | 1.1–1.2 | 533–547 | 9.8–11.4 M | 19.1 M | 880 ms | 3.9–4.0 s |
| 60, norm 16 | 63,600 | 60 | 1.1–1.2 | 56–58 | 1.1–1.2 M | 19.1 M | **93 ms** | 2.7–3.6 s |

The prediction of the earlier rounds holds at the landed rung: equal yield
and table, relations needed and probe work **÷ 9.5 / ÷ 9.3**, relation-phase
wall **÷ 9.5**. The whole run improves much less because the identical
**2.9 s table build is now 80–100 % of it** (Amdahl): the pair-table build is
the bottleneck of the batched pipeline at this size. Two honest caveats:
trials per relation is about 1 here (the base is above the 3-summand family
optimum at n = 53, so the sizing sweep should be redone at this rung before
quoting a best total), and this is the framework's table family, not the
production compact-orbit rank fixture (`koblitz_rank_fixture`), which
needs the multiplier closure ported separately.

## 3. What the profile says to do next

1. **Shard and parallelise the table build** (65 % of the build is cache-
   missing inserts): partition pairs by a mixed hash of the key into
   cache-sized shards, generate pairs in parallel chunks, merge shards in
   chunk order so first-wins is preserved, and probe through
   `shard(key).get`. Expected ≥ 2× on one core (sequential appends and
   L2-resident shard maps) and near-linear in cores beyond that.
2. **Batched `canon_with_shift`** (16 %): the lane-parallel least-rotation
   exists for `canon_in_place`; a with-shift variant removes the scalar loop.
3. **aarch64 field kernel.** The folding reduction and the inlined
   carry-less multiply are x86-only (`fold_tail` is gated on `x86_64`;
   `add_many` has `pclmulqdq` twins only); this host reduces products
   through the byte table and calls PMULL out of line. Applies to every
   binary-field loop on Apple silicon, uncertain gain (`batch_inv` already
   runs four interleaved chains).
4. **Port the multiplier closure to the production rank fixture**, where the
   rank precompute at `K = 600` is hours (n = 73: 18,400 s single thread).
   The ÷ 9.5 above is the estimate of what it buys there.
