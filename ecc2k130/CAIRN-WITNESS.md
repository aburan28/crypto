# Earning on cairn: what a distinguished point has to carry

**This is a scope, not a result.** Nothing here is measured on a GPU, and the
one number that decides the design — what eight counters cost the packed walk —
is named below as a measurement to run, not estimated and then quoted. What *is*
measured here is the basis reconciliation in §3, because it was the one thing
that could have made the whole idea impossible.

Under `AGENTS.md` §3 this change is none of the four classes. It does not lower
`S` and it does not move the ratio to any floor: it *raises* the cost of the
walk to buy a property that is not in the cost model at all — that a
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
| **B — counters on device** | packed + bitsliced walks, checkpoint, DP record | ~0 extra walk steps; cost is state and issue slots (§5) | ~0 | **100%** |

Route A is not a fallback to be embarrassed about: it needs no GPU change, no
checkpoint bump and no corpus change, it makes the format testable **today**,
and at the tranche the objective actually funds (`2^21` orbits, about one part
in 24,000 of the search) it is 1.67 extra GPU-hours. It does not scale to a
campaign, and the table says so in the column that matters.

Route B is what makes a campaign payable, and it is the one that needs hardware.

## 5. Route B, priced

**The packed backend is cheap and the bitsliced backend is not**, and the reason
is one line of each. `packedkernels.cuh` computes `const int j = 3 + ((hw >> 1) & 7);`
— an `int`, per worker — so incrementing one of eight counters is an ordinary
indexed read-modify-write of four bytes. `kernel.h` carries the same selector as
three *bitsliced words* `hb[1..3]`, so every lane in the word wants a different
counter and all eight must be touched with a ripple-carry every step.

| | added per step | against | share |
|---|---:|---:|---|
| packed: one 4-byte RMW | ~2–6 instructions | 2,189.75 instr/update (`THROUGHPUT-30B.md`) | **+0.09% to +0.27%** |
| bitsliced, 32-bit counters | `8 × (1 + 2×32)` = 520 lane-ops | ~1,540 lane-ops/update (`THROUGHPUT-CEILING.md`) | **+33.8%** |
| bitsliced, 16-bit + flush | `8 × (1 + 2×16)` = 264 lane-ops | same | **+17.1%** |

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

### What to change

- `include/walk.h` — `DpRecord` gains `counts[8]`; `WalkParams` gains the
  counter array.
- `include/kernel.h` — bitsliced increment in `run`, counters zeroed in `init`
  and in `reseed`, counters copied in `handleDistinguished`.
- `include/packedkernels.cuh`, `include/packedengine.cuh` — the same for the
  packed path, plus the allocation and the tile layout if (c) is chosen.
- `src/main.cu` — `DpFileRecord` gains the counters (32 → 64 bytes with (a));
  `CkptHeader` payload grows, so `checkpointVersion()` **must** go to 2 for the
  bitsliced engine and 3 for the packed one, and `ckptPayloadIsWhole` follows.
- `include/solver.h` — `rewalk` stays as the oracle the device is checked
  against; `solve` can then take counters directly and skip re-walking on a
  collision, which is a small speedup it gets for free.
- `aws/merge.py` — `RECORD` dtype and `RECORD_BYTES`; bucket files hold the old
  32-byte shape, so this is a corpus-format break (§6).
- new: the artifact emitter and the §3 permutation, with the popcount assertion.

### What to measure, and against what

`TUNING.md`'s paired comparison, not a before/after across builds: same CUDA,
same geometry, same arithmetic options, counters on and off, three benchmark
and three DP34 collection samples each, reporting the median and the spread.
The control is the same binary with the counter writes compiled out. Report
instructions/update from the receipt alongside B/s, because if throughput falls
by more than the instruction count did, state traffic became binding again and
layout (c) is the answer rather than a smaller counter.

## 6. Operational constraints

- **A checkpoint cannot survive this.** `aws/campaign.json` already says
  geometry "must never change while checkpoints exist"; counters are walk state,
  so every in-flight checkpoint is invalidated. This lands at a campaign
  boundary or it strands work.
- **The corpus format breaks.** DP files are headerless fixed-size records, so
  v1 and v2 files are indistinguishable by content. Either write v2 under a new
  extension or give the file a header; do not overload the size.
- **`dpWeight` must stay 34 campaign-wide** (campaign.json already warns), and
  the cairn job pins 34 too, so the two agree by construction — but a future
  cutoff change forks the cairn job id and is a new objective, not a tweak.

## 7. Validation

In this order, because each step is only worth running if the previous passed:

1. **Host, no GPU.** Emit an artifact from `rewalk` on the `GF(2^23)` and
   `GF(2^41)` instances and run cairn's own checker
   (`examples/certicom-ecdlp/checkers/ecc2k130_orbit_batch.py` shape, its toy
   twin for the small fields) over it. These instances have a planted `k`, so a
   collision must recover it — `AGENTS.md` §6's "a run whose answer was not
   checked against the planted secret" applies.
2. **Device against the oracle.** Extend `--verify N` to compare the device's
   counters with `rewalk`'s, not just the endpoint. A mismatch is the only
   failure mode this change can introduce that nothing else would catch.
3. **Cross-implementation.** Feed the emitted artifact to cairn's checker *and*
   to `orbit_dp.py verify`; both must accept, and `orbit_dp.py audit` must
   re-walk it clean.
4. **Interop on the real instance.** One artifact from an ECC2K-130 trail,
   accepted by the checker. This costs one trail and is the end-to-end proof.

## 8. What would make this not worth doing

State it before measuring, per `AGENTS.md` §4. Abandon Route B, and fall back
to Route A for tranches only, if the paired comparison shows the packed walk
losing **more than 2%** of its collecting rate under layout (c) with 16-bit
counters. Below 2% the cost is smaller than the spread between CUDA releases
this kernel has already absorbed; above it, the witness is buying an external
property at a price that shows up in the campaign's wall-clock, and that trade
should be made deliberately rather than by default.

The bitsliced path is exempt from that target and is expected to fail it: at
+17% it should stay witness-free and keep using Route A, which is what it is
for — it is the correctness oracle and the small-field path, not the production
walk.
