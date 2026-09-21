# Earning on cairn: what a distinguished point has to carry

**Route B is implemented.** The walk now carries the eight per-branch step
counts, they reach disk in a v2 corpus, and a claim built from one is accepted
by cairn's own checker.

**`main` has since grown Route A.** `src/witness.cpp` (`make witness`,
`test_witness.py`) emits the same claim artifact by *replaying* each corpus
record on the CPU. That is the other half of the table in §4 and it was the
right thing to build first: it needs no device change and it made the format
testable before any of this existed. This branch is the Route B half, so the
two are now reconciled rather than duplicated — `witness.cpp` reads both corpus
formats and uses the carried counts when they are there, replaying only a v1
corpus. The redundant Python emitter this branch had added is deleted.

What is measured and what is not:

| | |
|---|---|
| bitsliced backend (CPU and the default CUDA path) | **implemented and tested here**: counters verified against `Solver::rewalk` on 300 distinguished points of a `GF(2^83)` run and on runs on `GF(2^23)` and `GF(2^41)`, planted logarithms still recovered, and the witness algebra checked in `--test`. The comparison is not vacuous: a one-bit change to the branch bucket in `bumpCounts` is caught on the first record (§7.1). On `GF(2^131)` a reference re-walk is 2^25.27 scalar steps, so there the check is cairn's, below |
| packed backend (the RTX preset, and the campaign) | **implemented, not compiled**: there is no CUDA toolchain on the machine this was written on, so `make gpu` has not been run. The edits mirror the bitsliced ones site for site and are listed in §5 |
| end to end against cairn | **done**: 128 real weight-34 ECC2K-130 orbits, emitted by `build/witness` from the carried counters with **0 steps replayed**, accepted by `examples/certicom-ecdlp/checkers/ecc2k130_orbit_batch.py` in the merged objective and by `orbit_dp.py verify` against the pinned job |
| what the witness costs the bitsliced walk | **measured: +6.0%** median on the CPU backend (spread +5.4% to +7.1%), against +33.8% predicted (§5) |
| what the witness costs the packed walk | **not measured** — it needs the GPU, and §5 says what to run |

Under `AGENTS.md` §3 this change is still none of the four classes. It does not
lower `S` and it does not move the ratio to any floor: it *raises* the cost of
the walk to buy a property that is not in the cost model at all — that a
distinguished point becomes checkable by a stranger in `2^7.8` group operations
instead of the `2^25.27` it cost to make. It therefore does not belong on
`docs/index-calculus-scoreboard.html`, which prices attacks against rho; saying
so is §7's requirement met, not skipped.

## 1. What cairn now pays for

[`aburan28/cairn#159`](https://github.com/aburan28/cairn/pull/159) added a
piecework objective that pays per **orbit** of this walk
(`docs/design/orbit-piecework.md` there). One claim carries up to 64 elements:

```json
{"dps": [{"x": "<orbit>", "seed": "<64-bit walk seed>", "j": [n3, …, n10]}]}
```

- `x` is the orbit under `⟨−1⟩ × ⟨σ⟩`, named by the least cyclic rotation of the
  abscissa's coordinates in a **normal basis the objective pins**;
- `seed` is this client's walk seed, unchanged;
- `j` is the **witness**: how many steps took each of the eight branches. The
  step is `R ↦ R + σ^j(R) = [1 + s^j]R` and the endomorphism ring is
  commutative, so a trail of any length is `[μ]R₀` with
  `μ = ∏_j (1 + s^j)^{n_j}` whatever order the steps came in. A checker rebuilds
  `μ`, forms `[μ·α₀]P + [μ]Q`, and requires it to land on the claimed orbit —
  one double scalar multiplication.

Without `j` an orbit name is a low-weight bit string anyone can type, so the
objective would pay for nothing and a forged name would take the payment of
whoever reaches that orbit honestly. The witness is not decoration; it is the
whole reason the objective can exist.

## 2. Three things that already work

The scope is smaller than it looks, because this client already computes almost
all of it.

- **The counters exist.** `Solver<Cfg>::rewalk` in `include/solver.h` already
  returns `counts[8]`, and `Solver<Cfg>::multiplier` already builds
  `μ = ∏_j (1 + s^j)^{n_j}`. They were written for collision resolution and are
  exactly the witness. `--verify N` already calls `rewalk` on sampled records.
- **The orbit name is recoverable from the corpus we already write.** A
  `DpFileRecord` is `(seed, canon[3])` and carries no coefficients — but `canon`
  determines the orbit, and the orbit determines cairn's name. Checked in §3.
  **The 32-byte corpus format does not have to change for `x`.**
- **`s` and `ℓ` already agree.** `generated/eccF131.h` carries
  `S_DEC = 196511074115861092422032515080945363956` and cairn's job document
  carries the same value, derived independently there as the root of
  `s² + s + 2 ≡ 0 (mod ℓ)` that acts as σ on `⟨P⟩`.

So the only genuinely new thing is getting `j` out of a walk that does not
currently count.

## 3. The basis reconciliation (measured)

This was the risk worth retiring first, because getting it wrong is silent:
both sides produce well-formed 131-bit strings and simply disagree about which
orbit a point is on.

cairn pins a normal element γ and indexes coordinate `i` by `γ^(2^i)`, so σ is a
cyclic rotation. This client uses a *permuted* type-II ONB where σ is
`SQ_PERM`. Same basis **set**, different index order. Measured:

| check | result |
|---|---|
| cairn's γ expressed in this client's ONB coordinates | popcount **1** — it *is* one of our basis elements |
| its client index | **39** |
| `perm[i] = SQ_PERM^i(39)` is a bijection on `0…130` | yes |
| `perm` agrees with γ's own conjugate matrix on `P.x, P.y, Q.x, Q.y` | all four, values and weights |
| the same four against cairn's own implementation, from the other side | agree |
| client `canonical()` and cairn's least rotation pick the same representative | **no** — different member of the same orbit |
| cairn's name recovered from the client's `canon[3]` alone | **yes** |

Reproduce with `python3 cairn_basis.py`, which derives the index from
`Z_TO_ONB` and `SQ_PERM` in `generated/eccF131.h`, checks the map is a
bijection, and cross-checks it against a conversion that shares none of its
machinery — the challenge points in *polynomial* basis from `codegen/gen.py`,
taken through γ's own conjugate matrix. That independence is the point:
a cross-check built out of the permutation would pass no matter what the
permutation was. The first draft inverted `Z_TO_ONB` transposed, and this
check is what caught it.
`--canon` names the orbit of a corpus record: on `P.x` it returns
`3115f4e0f72a1d2c6ff5eca554d4f0bf`, which is what cairn's own tools return for
the same point from the other side.

Two consequences, and the second is the trap:

1. `x` costs no device state and no corpus change: it is a pure function of
   `canon[3]`.
2. **Do not ship index 39 as a constant.** Derive it at build time from
   `Z_TO_ONB` and the job's γ, and assert the popcount is 1. If the generator's
   basis ever changes, a hardcoded 39 keeps producing valid-looking orbit names
   for the wrong orbits, and every claim silently stops matching anyone else's.

## 4. Two routes, one table

Unit: cost per payable orbit, and the fraction of a campaign that becomes
payable. Both rows produce an identical artifact; they differ only in who pays
for the counters.

| | device change | cost per orbit | tranche of 2^21 orbits | payable fraction of a 128-GPU fleet |
|---|---|---|---|---|
| **A — host re-walk** (`src/witness.cpp`, on `main`) | none | one full trail again: 3.21 CPU core-s, or 2.87 ms of a GPU | +77.9 CPU core-days, or +1.67 GPU-h on top of the 1.67 GPU-h that produced it | linear in CPU budget: one core buys **0.089%** of one GPU; the fleet would need **143,418 cores** for all of it |
| **B — counters on device** (this branch) | packed + bitsliced walks, checkpoint, DP record | ~0 extra walk steps; the cost is throughput: **+6.0%** measured on the bitsliced CPU walk, unmeasured on the packed one (§5) | ~0 | **100%** |

Route A is not a fallback to be embarrassed about, and it is not hypothetical:
it is `src/witness.cpp` on `main`. It needs no GPU change, no checkpoint bump
and no corpus change, it makes the format testable **today**, and at the
tranche the objective actually funds (`2^21` orbits, about one part in 24,000
of the search) it is 1.67 extra GPU-hours. It does not scale to a campaign, and
the table says so in the column that matters.

Route B is what makes a campaign payable. It is built, and on the one backend
that can be measured here it costs 6.0% of the walk against Route A's 100%.
What still needs hardware is the backend the campaign actually runs.

The two are one tool, not two. `witness.cpp` frames the corpus by its magic and
takes the counts from a v2 record or replays a v1 one, so a campaign that
switches the counters on does not switch emitters, and every corpus written
before they existed still pays out. The measured difference on the same 128
ECC2K-130 records: **5,039,383 steps carried and 0 replayed**, 30.8 s wall for
128 claims; Route A on the same corpus would have had to walk all 5.04M steps
back. That ratio is the reason §4's last column reads the way it does.

### Route A, made cheaper

Two things were costing more than they had to, and both are fixed.

**The witness check was two scalar multiplications where one says the same
thing.** It rebuilt the endpoint as `[mu·alpha_0]P + [mu]Q`. But `R_0` *is*
`[alpha_0]P + Q`, and the walk already builds it, so `[mu]R_0` is the identical
statement for 196 point operations instead of 392. `WalkResult` now carries
`startPt` and both the replay check and `fromCounts` use it. The carried path
gets this too: emitting 128 ECC2K-130 claims went from **30.8 s to 18.2 s**
single-threaded, which is the whole of that 1.69x.

**The replay was serial across records that do not depend on each other.** It
now runs under OpenMP, with the emission left serial so the artifact bytes are
whatever a one-record-at-a-time run would have produced. That is checked rather
than asserted: the same corpus at `OMP_NUM_THREADS` 1, 2 and 4 gives
byte-identical output on both paths.

| | before | after, 1 thread | after, 4 threads |
|---|---:|---:|---:|
| replay, 8 records / 20,394 steps | 15.18 s | 12.78 s | **4.61 s** |
| carried, 128 records | 30.8 s | 18.18 s | **4.90 s** |

So 3.3x on the replay and 6.3x on the carried path, on four cores.

**What this does not change is the shape of the cost**, and that is the honest
limit. A serial replay costs the SUM of the trail lengths; parallelism divides
that sum by the core count, but a single trail still cannot be split. On the
128-record ECC2K-130 corpus the sum is 5,039,383 steps and the longest single
trail is 77,146 -- a factor of 65 between what a replay must pay and what a
*batched* replay bounded by the longest trail would pay. Closing that gap means
replaying through the bitsliced walk, 512 lanes at a time, which is the walk
this client already has and the counters Route B already added: a replay kernel
is Route B's kernel pointed at a corpus's seeds instead of at fresh ones. It
needs seed injection into the reseed path, which is the hottest code in the
campaign, so it is written up here rather than rushed.

## 5. Route B, priced

The two backends pay for the witness by different mechanisms, and the reason is
one line of each. `packedkernels.cuh` computes `const int j = 3 + ((hw >> 1) & 7);`
— an `int`, per worker — so incrementing one of eight counters is an ordinary
indexed read-modify-write of four bytes. `kernel.h` carries the same selector as
three *bitsliced words* `hb[1..3]`, so every lane in the word wants a different
counter and all eight must be touched with a ripple-carry every step.

This note first said the bitsliced one would therefore be expensive. Measured,
it is not: **+6.0%** median on the CPU backend (three paired rounds, +5.96%,
+7.13%, +5.43%; `WITNESS=0` median 40.026 M it/s against `WITNESS=1` 37.762 M
it/s at m = 131, weight 34), against +33.8% predicted. The estimate is left in
the table beside the measurement rather than quietly corrected, because the gap
is the useful part.

An earlier revision of this note measured +4.2% against the pre-merge tree. The
number moved because the tree did, not because the counting changed: `main`'s
walk is faster now, so the same per-step counting work is a larger fraction of
it. Both numbers are from the same paired method; the one above is the one that
describes the code in this branch.

| | added per step | against | predicted | measured |
|---|---:|---:|---:|---:|
| packed: one 4-byte RMW | ~2–6 instructions | 2,189.75 instr/update (`THROUGHPUT-30B.md`) | +0.09% to +0.27% | *needs the GPU* |
| bitsliced, 32-bit counters | `8 × (1 + 2×32)` = 520 lane-ops | ~1,540 lane-ops/update (`THROUGHPUT-CEILING.md`) | +33.8% | **+6.0%** |
| bitsliced, 16-bit + flush | `8 × (1 + 2×16)` = 264 lane-ops | same | +17.1% | not built |

`main` has since made the same correction from the other side
(*Correct the counter cost: it was priced against the wrong backend*): the
campaign runs packed, so the ~105-slot / ~4% / ~1,800 GPU-hour figure that
once argued against carrying counters at all was pricing the *bitsliced*
layout against the *packed* preset's budget. Its conclusion — a scalar counter
on the packed path is ~+0.09%, about 36 GPU-hours across the campaign — agrees
with the first row above, which is the reassuring part: two independent
routes, one number.

One difference worth naming, since the two descriptions are not identical.
That correction costs the scalar counter at two instructions by packing eight
8-bit fields into one `u64` and flushing before a field can overflow. What is
implemented here is a single indexed `+= 1u` into eight 32-bit fields per walk:
8 bytes of per-walk state against 32, in exchange for needing no flush. Same
cost class, different trade, and **neither is measured** — that still needs a
GPU.

**Two distinguishing weights appear in this note and they are different
things.** The cairn job pins weight **34**, where a trail averages `2^25.27`
steps; the campaign's own `dpWeight` is **32**, where it is `2^28.41`. Every
`2^25.27` here is the cairn-objective figure and is correct for it. Campaign
accounting uses 32.

The measurement, so the number can be argued with: `make cpu WITNESS=0` and
`WITNESS=1`, four threads on a four-core host, `--curve 131 --dp-weight 34
--launches 900 --verify 0`, three rounds with the two binaries interleaved
within each round so drift hits both. Control 40.822 / 40.026 / 39.811 M it/s,
witness 38.526 / 37.361 / 37.762, giving +5.96% / +7.13% / +5.43% per round and
**+6.00%** on the medians. Three paired samples on a shared four-core box is a
thin measurement and the spread says so; what it is good enough to establish is
that the cost is single-digit percent and not the third of the walk that was
predicted.

`make` does not relink on a flag change alone, so both builds were preceded by
`touch src/main.cu` and the two binaries were compared with `cmp` before the
first round. Without that the control is a copy of the treatment and the result
is a perfect 0.00%, which is exactly what success looks like.

**The bitsliced prediction was wrong, and by a lot.** It assumed a branch-free
ripple over all `ECC_COUNT_BITS` of all eight counters. The implementation's
carry loop stops as soon as no lane in the word is still carrying, and a lane
carries past bit `b` only when its low `b` bits are all set — so on a 64-lane
word it runs about seven bits, not 32. The other half of the gap is the
denominator: a walk step is not 1,540 lane-ops of counting-comparable work, it
is a batched field inversion amortised over `BATCH` slots, and the counters are
a much smaller share of that than the lane-op ratio suggested.

Two things the CPU number does *not* transfer to the GPU. The early exit is a
data-dependent loop bound, so on a warp it costs the maximum across 32 threads
rather than each thread's own; and the packed backend is a different mechanism
entirely, one scalar read-modify-write against a kernel that is issue-bound at
95% of the part's best measured rate. Neither is predicted here from the other.

The packed kernel is at **74.1 lane-instructions per SM-clock, 95% of the best
mixed rate ever measured on the part**, so it is issue-bound and added
instructions convert to lost throughput nearly one for one. +0.27% is the
pessimistic end of a rounding error. That is the claim to *measure*, not to
believe: it assumes issue, not state traffic, stays binding — and
`THROUGHPUT-30B.md` records that the binding constraint has already swapped
once in this kernel's history.

State is the real question, because that is what swapped it last time:

| counter layout | hot bytes/slot | hot footprint at 385,024 workers × batch 16 | cold |
|---|---:|---:|---:|
| (a) 8 × `u32` | 32 | 188.0 MiB | — |
| (b) 7 × `u32`, eighth from `iters` | 28 | 164.5 MiB | — |
| (c) 8 × `u16` hot, flushed to 8 × `u32` once per launch | 16 | **94.0 MiB** | 188.0 MiB, touched once per 1,024 steps |
| *current coordinate state* | | *399.5 MiB* | |
| *current metadata (`dead`/`seed`/`startIter`)* | | *117.5 MiB* | |

Against a 128 MB L2 that `BATCH-TUNING.md` already shows the kernel sitting
right against — two clients on one GPU "land back in the regime batch 32
measured at 8.5 B/s" — (a) is a 36% increase in hot state and (c) is 18%.
(b) is free relative to (a) because `iters` is already derivable from `now` and
`startIter`, both of which the record already holds.

### What changed

- `include/walk.h` — `DpRecord` gains `counts[ECC_JCOUNT]`; `ECC_WITNESS`,
  `ECC_COUNT_BITS` and the scalar counter layout the packed path indexes with.
- `include/kernel.h` — `WalkParams` gains the counter array; `bumpCounts`
  (bitsliced, with a carry loop that stops when no lane is still carrying),
  `clearCounts` on start and on restart, and the per-lane read-out in
  `handleDistinguished`.
- `include/packedkernels.cuh` — the same three points, scalar: zero in `init`
  (which also serves restart), one `+= 1` after `j` is computed, copy at the
  report.
- `include/packedengine.cuh`, `src/main.cu` — allocation, `bytesPerThread`,
  and the checkpoint. Counters are walk state, so a checkpoint written with
  them cannot be read by a build without them and the version has to say so:
  `ECC_CKPT_BUMP` adds **16**, not 1. The versions already in use are 1
  (bitsliced), 2 (packed) and 3 (`main`'s packed table walk), so a bump of one
  would have put the witnessed packed format on 3 and given two different
  layouts the same number — the one thing a format version exists to prevent.
  `WITNESS=0` leaves every version at what it was.
- `include/walk.h`, `Makefile` — the table walk picks its addend from a table
  instead of computing `σ^j(R)`, so a trail there is not `[∏_j (1+s^j)^{n_j}]R₀`
  and eight counters say nothing about it. `WITNESS=1` with `WALK_TABLE=1` is a
  compile error rather than a build that writes meaningless — or, worse,
  uninitialised — counts, since that path's report sites never fill them. The
  *default* follows the walk (`WALK_TABLE=1` implies `WITNESS ?= 0`), because
  several recipes select the table walk and say nothing about the witness; they
  want a binary, not a diagnostic. Only an explicit request for both fails.
- `include/durable.h` — `CorpusOutput::openFile` takes a header size and
  records the size it saw **under the lock**, which is what decides whether a
  fresh corpus gets its v2 header. The append handle's own `ftell` cannot
  answer that: a stream opened `"ab"` has an implementation-defined position
  until the first write, and on a library that reports 0 it would put a header
  in the middle of an existing corpus.
- `include/solver.h` — `Solver::fromCounts` rebuilds a walk result from a
  carried witness with no replay. What it checks is narrow and deliberate: the
  counts sum to the claimed trail length and the endpoint is distinguished.
  Comparing the endpoint against `μ` would be a tautology there, because the
  endpoint *is* built from `μ`; the check that binds is the caller's orbit
  comparison, which is the same one the payer makes.
- `src/main.cu` — corpus v2 behind a magic, `--verify` compares the device's
  counters against `Solver::rewalk`, and a `testWitness` in `--test`.
- `src/witness.cpp` (Route A, from `main`) — reads either corpus format, takes
  the carried counts when they are there and replays only a v1 record. The
  `μ`-reproduces-the-endpoint check stays where it is an independent check (the
  replay) and is skipped where it would be a tautology.
- `test_witness.py` — frames the corpus from its magic rather than assuming
  32-byte records, emits across several batches, and verifies each batch
  separately, since `orbit_dp.py verify` reads one artifact per file and the
  emitter writes one per line. With ≤ `max_batch` records those look the same,
  which is why it was worth writing out. It also derives a v1 corpus from the
  v2 one and requires the two to produce **identical** `j` vectors — the direct
  test that a counter carried on the device is the number a replay recovers.
- `aws/merge.py` — reads either format. The witness is deliberately *not*
  carried into the bucket files: merging matches two seeds against one orbit
  key and nothing else, and widening every bucket record to carry something
  the merge never reads would cost the pass its margin.
- deleted: `cairn_artifacts.py`. It turned a v2 corpus into claims in Python,
  which is what `src/witness.cpp` now does for both formats. Two emitters that
  must agree byte for byte is a liability, not redundancy.
- kept: `cairn_basis.py`, which is not an emitter. It derives the basis index
  independently and cross-checks it against the challenge constants by a route
  that shares no machinery with the permutation (§3) — the check that caught a
  transposed matrix inverse once already.

One thing deliberately not done: `Solver::Entry` still holds only `(seed,
iters)`, so resolving a collision still re-walks both trails rather than
reading their witnesses. The witness would save those two walks, but `Entry`
is the per-orbit store that `--load-max` exists to bound, and growing it by 40
bytes to save work that happens once per campaign is the wrong trade.

### What to measure, and against what

`TUNING.md`'s paired comparison, not a before/after across builds: same CUDA,
same geometry, same arithmetic options, `WITNESS=1` against `WITNESS=0`, three
benchmark and three DP34 collection samples each, reporting the median and the
spread. Report instructions/update from the receipt alongside B/s, because if
throughput falls by more than the instruction count did, state traffic became
binding again and layout (c) is the answer rather than a smaller counter.

`make` does not rebuild on a flag change alone, so `touch src/main.cu` between
the two builds and `cmp` them before trusting a round, or the control is a copy
of the treatment. That mistake was made once already while measuring the CPU
number above, and it reports a perfect 0.00% difference, which is exactly what
it looks like when it works.

## 6. Operational constraints

- **A checkpoint cannot survive this.** `aws/campaign.json` already says
  geometry "must never change while checkpoints exist"; counters are walk state,
  so every in-flight checkpoint is invalidated. This lands at a campaign
  boundary or it strands work.
- **The corpus grew and the format moved.** v2 records are 72 bytes against
  v1's 32, so the full `2^35.5`-orbit corpus goes from about 1.6 TB to 3.6 TB.
  The format is told apart by an `ECC2KDP2` magic rather than by size, and a
  build refuses to append one to a non-empty file of the other. **Every
  uploaded delta carries the header**, because `aws/worker.py` ships byte
  ranges as standalone objects and the merge frames each on its own — without
  that, only a slot's first delta would announce the format and every later one
  would be read as v1, which mis-frames every record in it and is invisible
  until the merge reports orbits nobody walked.
- **Bucket files stay v1.** `aws/merge.py` reads either corpus format but
  writes the old 32-byte shape into its buckets: the merge matches two seeds
  against one orbit key and never reads a witness, so widening every bucket
  record to carry one would cost the pass its margin. Whatever wants the
  witness reads the corpus.
- **`dpWeight` must stay 34 campaign-wide** (campaign.json already warns), and
  the cairn job pins 34 too, so the two agree by construction — but a future
  cutoff change forks the cairn job id and is a new objective, not a tweak.

## 7. Validation

What was run, in the order each became worth running:

1. **The device against the oracle.** `--verify N` now compares the walk's
   counters with `Solver::rewalk`'s, not just the endpoint — a wrong count
   still produces a well-formed record, a well-formed claim, and a `mu` that
   lands on somebody else's orbit, so nothing else would catch it. 300
   distinguished points of a `GF(2^83)` run, plus every distinguished point of
   runs on `GF(2^23)` and `GF(2^41)` at several distinguishing weights,
   including trails of 448 steps where the carry actually propagates, and the
   counts sum to `iters` on all of them. Not on `GF(2^131)`: one reference
   re-walk there is 2^25.27 scalar steps, which is why step 4 exists and does
   the same job from the other side.

   **And the oracle is not vacuous.** Inverting one bit of the branch selector
   in `bumpCounts` — `(k & 1) ? jb[1] : ~jb[1]` to its negation, which permutes
   counts between buckets `k` and `k^1` and leaves their sum alone — is caught
   on the first verified record (`witness[2] = 4, the reference walk took 3
   steps on that branch`, exit 3). The first mutation tried was not caught, and
   that was informative rather than alarming: dropping the `| dp` term from the
   live mask changes nothing, because `handleDistinguished` reads the counters
   *before* the bump and `clearCounts` zeroes them on revive. That term is
   defensive, not load-bearing, and now the note says so instead of the reader
   having to work it out.
2. **The algebra, in `--test`.** `testWitness` walks a bounded few steps and
   requires `[mu]R_0` to be the point it reached, on every curve the suite
   covers including `GF(2^131)`. The identity holds at every step, which is
   what makes it a unit test rather than a 2^25.27-step one; the endpoint half
   runs on one trail per curve rather than four, because the reference's
   scalar multiplication is `O(m^3)` and four of them at m = 131 turned a
   21-second suite into a six-minute one.
3. **The planted logarithms still come out.** `make break-small` on both
   backends, which is `AGENTS.md` §6's "a run whose answer was not checked
   against the planted secret".
4. **Interop on the real instance.** 128 genuine weight-34 ECC2K-130 orbits
   collected by this client, written to a v2 corpus, turned into two 64-element
   claims by `build/witness` **straight from the carried counters — 5,039,383
   steps carried, 0 replayed** — and accepted by
   `examples/certicom-ecdlp/checkers/ecc2k130_orbit_batch.py`, the checker
   whose hash is pinned inside the merged objective, and independently by
   `orbit_dp.py verify` against `jobs/ecc2k130.json`. Verification is 1.3 s per
   64-point batch against the ~2^25.27 steps each of those points cost, which
   is the asymmetry the whole design is for. Bumping one counter is refused
   (`the witness does not reach this orbit`), as is tampering with the seed,
   the orbit name, or the order of the counters.
5. **The corpus formats.** A fresh v2 corpus is `16 + 72N` bytes with exactly
   one header; appending to a real ECC2K-130 v2 corpus grows it 128 → 134
   records, still aligned, still one header. A v2 file round-trips through
   reload (and a reload collision still solves); a v1 file written by a
   `WITNESS=0` build still reads. Both append directions are refused with exit
   2 and the corpus untouched: a `WITNESS=1` build onto a v1 corpus, and a
   `WITNESS=0` build onto a v2 one.
6. **The suite.** `make test`'s targets pass except `test-clmad` and
   `test-clmad-square`, which fail identically on a pristine `origin/main`
   worktree (`codegen/prove_native_square.py`: *native branch is outside the
   proved source contract*, plus a guard test). They are not this branch's and
   are not touched by it. `test-table-walk-host` *was* broken by the guard
   above and is fixed here, not exempted.

What has *not* been run: `make gpu`. There is no CUDA toolchain here, so the
packed edits are uncompiled — see the status table at the top.

## 8. What would make this not worth doing

State it before measuring, per `AGENTS.md` §4. Abandon Route B, and fall back
to Route A for tranches only, if the paired comparison shows the packed walk
losing **more than 2%** of its collecting rate under layout (c) with 16-bit
counters. Below 2% the cost is smaller than the spread between CUDA releases
this kernel has already absorbed; above it, the witness is buying an external
property at a price that shows up in the campaign's wall-clock, and that trade
should be made deliberately rather than by default.

That target stands for the packed walk and is still unmeasured.

The exemption written here for the bitsliced path does not: it was granted on
an estimate of +17%, and the path measures **+6.0%** on the CPU
backend. So the bitsliced walk carries the witness by default too, and
`WITNESS=0` is there for whoever wants the old rate rather than for whoever
wants correctness. What that does not settle is the same path under CUDA,
where the carry loop's early exit costs the warp's maximum rather than each
thread's own; if that turns out to matter, the 16-bit-plus-flush layout is the
next thing to build, and it is already costed above.
