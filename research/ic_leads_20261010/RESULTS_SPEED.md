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

## 2b. Sharded, parallel table build, and the landed n = 53 rung end to end

The profile of §1 put 65 % of the build in cache-missing inserts, so the
table is now `2^bits` maps (shards of about 2^15 entries, at most 2^10)
selected by a hash of the key that is independent of the in-shard hash.
Pair sums are generated in parallel chunks of consecutive representatives,
each appending to one buffer per shard, and every shard is filled from the
chunks in chunk order, so the first pair to reach a key wins as in the
sequential build, **for any thread and shard count** (the equivalence test
forces 1, 8 and 64 shards on pools of 1 and 5 threads and requires every
entry, the charged operations and the `m = 3` probes to equal the scalar
reference). Threads are used only above 2^21 pairs: the first version
spun up chunks for any table and was slower than scalar in the median at
n = 23 (0.69×) on this loaded host, so a small build now runs the same
algorithm on one thread.

K = 600 at n = 41 (14.7 M entries, 14,784,600 charged additions, identical
in both):

| table build | min | median |
|---|---:|---:|
| scalar | 8.9 s | 9.2 s |
| batched, unsharded (§1) | 1.40 s | 1.45 s |
| batched, sharded, **one thread** | **0.78 s (11.4×)** | 0.79 s (11.9×) |
| batched, sharded, all cores | 0.29 s (30×) | 0.50 s (18×) |

Sharding alone, with no extra thread, halves the unsharded batched build.

**The landed n = 53 rung** (`K_0`, 52-bit r, 600 plain columns = the ledger's
`K`, 63,600 points, 19.1 M charged table additions, `m = 3`), old scalar
`ic` binary and new binary alternating, 2 rounds × 2 repeats, all verified,
**counted operations, trials, relations, lookups, canonicalisations and rank
identical on every paired row**:

| | old scalar | new, one thread (matched) | new, all cores |
|---|---:|---:|---:|
| pair-table build | 16.6–22.4 s | 1.08–1.22 s (15–19×) | 0.41–0.76 s (34×) |
| relation phase | 8.3–10.1 s | 0.61–0.76 s (13×) | 0.62–0.76 s (13×) |
| **whole pipeline** | **25.5–32.6 s** | **1.8–2.1 s (≈ 14×)** | 1.1–1.6 s (≈ 21×) |

and with the multiplier closure on top (60 columns, norm ≤ 16, the same
63,600 points and the same table): whole pipeline 1.3–1.4 s on one thread
(≈ 20× the scalar baseline for the ledger's shape) and 0.6–0.9 s on all
cores (≈ 30×); its relation phase is 78–94 ms against 608–758 ms for 600
plain columns, the ÷ 9 of §2. The single-thread column is the
matched-resource figure; the all-cores column used the host's threads for
the table build only.

## 3. What the profile says to do next

1. **Batched `canon_with_shift`** (16 % of the unsharded build, now a larger
   share of what is left): the lane-parallel least-rotation exists for
   `canon_in_place`; a with-shift variant removes the scalar loop.
2. **The relation phase at `m = 3`** is now the larger cost at n = 53
   (0.6–0.76 s against the 1.1 s table build on one thread): per probe it is
   an `add_many` element, a lane-parallel canonicalisation and a sharded
   lookup. The lazy variant of `add_many` (no ordinate) saves one of about
   five field multiplications per element.
3. **aarch64 field kernel.** The folding reduction and the inlined
   carry-less multiply are x86-only (`fold_tail` is gated on `x86_64`;
   `add_many` has `pclmulqdq` twins only); this host reduces products
   through the byte table and calls PMULL out of line. Applies to every
   binary-field loop on Apple silicon, uncertain gain (`batch_inv` already
   runs four interleaved chains).
4. **Port the multiplier closure to the production rank fixture**
   (`koblitz_rank_fixture`), where the rank precompute at `K = 600` is
   hours (n = 73: 18,400 s single thread). The ÷ 9 above is the estimate of
   what it buys there; the files are free in Conductor since T-261 was
   cancelled, but PR 1605 still touches them.
