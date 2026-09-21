# Round 0015 pre-registration: the scan, and what the native clock can see

Written before any measured stage of this round was run.

## Where the IC arm's own cost sits after round 0014

A line-level profile of the round-0014 winner puts the IC arm's own work on
the largest cell in two phases: the block scan that decomposes each probe
(1.10M–1.76M instructions on `n23a0`, depending on how many probes fail)
and the pair-table build (0.82M). Inside one scanned rest, the cost divides
as: the batched inversion of the block's denominators 25%, naming the rest's
orbit 29% (the coordinate change 5%, the rotation minimum 24%), the field's
reduction 15%, and the rest table probe and bookkeeping.

The pair table's shape is not the lever. With `K = 8` orbits and a build
cost of about 680 instructions per stored sum against about 380 per scanned
rest — both measured, where the code's model assumes them equal — the row
count that minimises the expected total is still the 2 the model picks
(801,840 / 796,010 / 962,880 expected instructions for 1, 2 and 3 rows on
`n23a0`). The model's assumption is wrong and its answer is right.

Two changes were tried in development and **dropped**:

- *Reduction by folding a sparse modulus* instead of by byte tables. Every
  irreducible this pipeline uses is a trinomial or a pentanomial, so `z^n`
  folds into at most five low terms. It is **22% to 47% worse**: the table
  path needs only `⌈(n−1)/8⌉` lookups, which is two or three at these
  degrees, while folding a pentanomial costs four shift-mask-xor triples
  twice over. The tables were already the better answer.
- *Indexing those tables through a slice of exactly 256 words* so no bound
  is checked per lookup. A wash: 0.874 against 0.884 on `n23a0` but 1.009
  against 0.994 on `n13a0`, so it loses where the campaign's small cells
  weigh most.

## What changes, and for whom

**Nothing shared changes.** The baseline is the round-0014 winner source
exactly as it was sealed, so the incumbent's ratios should reproduce round
0014's up to fresh-fixture variation, and whatever the challengers show is
their own.

**Challengers**, both IC-only, both byte-identical to the incumbent on all
304 frozen round-0014 fixtures:

- `scan` ([round15-scan.patch](round15-scan.patch)): the orbit name of an
  abscissa found from the longest circular run of zeros — the runs by
  doubling, their maximal length by a binary search over those masks, and
  only the one or two rotations that put such a run at the top compared —
  instead of trying all `n − 1`, with the rotation loop kept as the test
  reference; and the block scan keeping only each rest's **abscissa**,
  which is all the orbit lookup needs, with the ordinate computed for the
  one rest in a few hundred whose orbit is stored. Development: 0.874–0.994
  of the incumbent's instructions by cell, strongest on `n23a0`.
- `scan_io` ([round15-fastio.patch](round15-fastio.patch) on `scan`): adds
  the report serialisation change. Development: 0.852–0.942, **the largest
  IC-only step this campaign has measured**.

## What this round expects, and why it may not be enough

Both challengers are expected to pass the instruction gate. On the native
clock the pre-registration expects them to sit at the edge of what the
measurement can resolve, and the reasoning is stated here so the outcome
can be checked against it rather than explained afterwards:

- The IC arm's own phases are about two thirds of its instructions, so a
  10% cut in them is about 6.5% of the job's instructions.
- Round 0014's instrumentation found that the job's own execution is about
  60% of the measured process wall; the rest is the spawn, the image
  page-in and the reap the protocol charges to every trial. So 6.5% of
  instructions is about 4% of the measured wall.
- Round 0014's challenger carried a paired native interval of
  [0.9724, 1.0163] about a point of 0.9947 — a half-width of 2.2%.

A 4% effect against a 2.2% half-width passes only if it lands in full. The
same serialisation change has now failed this gate in rounds 0011, 0012,
0013 and 0014 at 3.5–8% of instructions, which is the same arithmetic seen
four times. If `scan_io` fails too, the finding this round records is not
about the change but about the gate: **an IC-only change worth less than
about a tenth of the job cannot be resolved by a clock that spends a third
of its measurement on process creation**, and future IC-only work should
be judged on instructions with the native gate read as a no-regression
check, or else pursued in steps large enough to clear the floor.

## Objective and gates

`--objective rho`, as rounds 0007–0014. Incumbent: the round-0014 winner
configuration (`batch_trials: 1`) on the unchanged round-0014 source. A
challenger is promoted only if it passes the no-regression gate against
the incumbent (instruction ratio at most 0.98 with the paired upper limit
below one, native upper limit below one, no cell more than 10% worse) and
the strict rho gate on confirmation and replay. The decision records the
winner's `winner_over_rho`, `beats_rho_strict` and `rho_parity`.

## Parent, seed, budget

Fresh seed 2026091615; target count 1; pilot profile; 2,400 paired-job
budget; one pinned CPU; 8 GiB cap; 60-second watchdog; Valgrind 3.22.0
`Ir`; native progress recorded; the round-0010 evaluator (spawn without a
fork, caps before job delivery).

## Boundary, floor, class, honesty

Unit and boundary unchanged. Base support, `m = 3`, no direct relations,
every verification obligation, the worker's phase dumps and report fields
are unchanged; the final check is still `[d]G = Q` in the general
arithmetic on every recovered scalar. The decomposition searched, the
relations accepted and the certificates written are identical — the
challengers change how a rest is named and how much of it is computed
before it is discarded, not which rests decompose. Class: engineering.
Absolute times and instruction counts are comparable with rounds
0010–0014 (same evaluator, same host, same compiler) and not with earlier
rounds. No arithmetic-complexity, family-wide or cryptographic-size claim
follows. Fresh fixtures from the new seed; every failure retained; panels
stay separate.
