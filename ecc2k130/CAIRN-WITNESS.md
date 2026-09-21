# Earning on cairn: what a distinguished point has to carry

**Route B is implemented.** The walk now carries the eight per-branch step
counts, they reach disk in a v2 corpus, and a claim built from one is accepted
by cairn's own checker. What is measured and what is not:

| | |
|---|---|
| bitsliced backend (CPU and the default CUDA path) | **implemented and tested here**: counters verified against `Solver::rewalk` on every distinguished point of runs on `GF(2^23)` and `GF(2^41)`, planted logarithms still recovered, and the witness algebra checked in `--test`. On `GF(2^131)` a reference re-walk is 2^25.27 scalar steps, so there the check is cairn's, below |
| packed backend (the RTX preset, and the campaign) | **implemented, not compiled**: there is no CUDA toolchain on the machine this was written on, so `make gpu` has not been run. The edits mirror the bitsliced ones site for site and are listed in §5 |
| end to end against cairn | **done**: 53 real weight-34 ECC2K-130 orbits, emitted by `cairn_artifacts.py` and accepted by `examples/certicom-ecdlp/checkers/ecc2k130_orbit_batch.py` in the merged objective |
| what the witness costs the bitsliced walk | **measured: +4.2%** on the CPU backend, against +33.8% predicted (§5) |
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
| **A — host re-walk** | none | one full trail again: 3.21 CPU core-s, or 2.87 ms of a GPU | +77.9 CPU core-days, or +1.67 GPU-h on top of the 1.67 GPU-h that produced it | linear in CPU budget: one core buys **0.089%** of one GPU; the fleet would need **143,418 cores** for all of it |
| **B — counters on device** | packed + bitsliced walks, checkpoint, DP record | ~0 extra walk steps; the cost is throughput: **+4.2%** measured on the bitsliced CPU walk, unmeasured on the packed one (§5) | ~0 | **100%** |

Route A is not a fallback to be embarrassed about: it needs no GPU change, no
checkpoint bump and no corpus change, it makes the format testable **today**,
and at the tranche the objective actually funds (`2^21` orbits, about one part
in 24,000 of the search) it is 1.67 extra GPU-hours. It does not scale to a
campaign, and the table says so in the column that matters.

Route B is what makes a campaign payable. It is built, and on the one backend
that can be measured here it costs 4.2% of the walk against Route A's 100%.
What still needs hardware is the backend the campaign actually runs.

## 5. Route B, priced

The two backends pay for the witness by different mechanisms, and the reason is
one line of each. `packedkernels.cuh` computes `const int j = 3 + ((hw >> 1) & 7);`
— an `int`, per worker — so incrementing one of eight counters is an ordinary
indexed read-modify-write of four bytes. `kernel.h` carries the same selector as
three *bitsliced words* `hb[1..3]`, so every lane in the word wants a different
counter and all eight must be touched with a ripple-carry every step.

This note first said the bitsliced one would therefore be expensive. Measured,
it is not: **+4.2%** on the CPU backend, against +33.8% predicted. The estimate
is left in the table beside the measurement rather than quietly corrected,
because the gap is the useful part.

| | added per step | against | predicted | measured |
|---|---:|---:|---:|---:|
| packed: one 4-byte RMW | ~2–6 instructions | 2,189.75 instr/update (`THROUGHPUT-30B.md`) | +0.09% to +0.27% | *needs the GPU* |
| bitsliced, 32-bit counters | `8 × (1 + 2×32)` = 520 lane-ops | ~1,540 lane-ops/update (`THROUGHPUT-CEILING.md`) | +33.8% | **+4.2%** |
| bitsliced, 16-bit + flush | `8 × (1 + 2×16)` = 264 lane-ops | same | +17.1% | not built |

The measurement, so the number can be argued with: `make cpu WITNESS=0` and
`WITNESS=1`, four threads on a four-core host, `--curve 131 --steps 256
--launches 25 --verify 0`, runs interleaved. Control 8.066 and 8.098 M it/s,
witness 7.779 and 7.711, so 4.2% on the means and 3.9–4.4% pairwise. Two
samples each on a shared box is a thin measurement and the spread says so;
what it is good enough to establish is that the cost is single-digit percent
and not the third of the walk that was predicted.

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
  and the checkpoint. Counters are walk state, so `checkpointVersion` moves to
  2 for the bitsliced engine and 3 for the packed one, and `WITNESS=0` leaves
  both at what they were.
- `src/main.cu` — corpus v2 behind a magic, `--verify` compares the device's
  counters against `Solver::rewalk`, and a `testWitness` in `--test`.
- `aws/merge.py` — reads either format. The witness is deliberately *not*
  carried into the bucket files: merging matches two seeds against one orbit
  key and nothing else, and widening every bucket record to carry something
  the merge never reads would cost the pass its margin.
- new: `cairn_artifacts.py` turns a v2 corpus into claims, and `cairn_basis.py`
  is the basis map it needs.

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
the two builds or the control is a copy of the treatment. That mistake was
made once already while measuring the CPU number below, and it reports a
perfect 0.00% difference, which is exactly what it looks like when it works.

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
   lands on somebody else's orbit, so nothing else would catch it. Every
   distinguished point of runs on `GF(2^23)` and `GF(2^41)` at several
   distinguishing weights, including trails of 448 steps where the carry
   actually propagates, and the counts sum to `iters` on all of them. Not on
   `GF(2^131)`: one reference re-walk there is 2^25.27 scalar steps, which is
   why step 4 exists and does the same job from the other side.
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
4. **Interop on the real instance.** 53 genuine weight-34 ECC2K-130 orbits
   collected by this client, written to a v2 corpus, turned into a claim by
   `cairn_artifacts.py`, and accepted by
   `examples/certicom-ecdlp/checkers/ecc2k130_orbit_batch.py` — the checker
   whose hash is pinned inside the merged objective — and independently by
   `orbit_dp.py verify`. Tampering with one counter, the seed, the orbit name,
   or the order of the counters is refused by the checker in every case.
5. **The corpus formats.** A v2 file round-trips through reload; a v1 file
   written by a `WITNESS=0` build still reads; and a v2 build refuses to
   append to a non-empty v1 corpus rather than producing a file neither reader
   can frame.

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
an estimate of +17%, and the path measures **+4.2%** on the CPU
backend. So the bitsliced walk carries the witness by default too, and
`WITNESS=0` is there for whoever wants the old rate rather than for whoever
wants correctness. What that does not settle is the same path under CUDA,
where the carry loop's early exit costs the warp's maximum rather than each
thread's own; if that turns out to matter, the 16-bit-plus-flush layout is the
next thing to build, and it is already costed above.
