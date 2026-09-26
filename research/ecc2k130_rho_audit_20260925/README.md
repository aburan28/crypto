# ECC2K-130 rho audit: walk length, duplicate points, seed resets, fruitless cycles

2026-09-25. An audit of the live ECC2K-130 campaign's walk, its seeds and
its distinguished-point store, and the Rust guards that come with it
(`ecc2k-guard`, `src/cryptanalysis/ecc2k130_guard.rs`). Nothing here is a
speed or cost-per-solve claim; under `AGENTS.md` §3 every item is
**accounting** (what the campaign already does, measured and written down)
or tooling. Campaign figures are from the frozen public feed
`live_status_20260925T050509Z.json` unless another source is named.

## Answers

### 1. Average steps per walk

**2^28.41 ≈ 3.57 × 10^8 steps per completed walk** at the live cutoff
`HW(x) ≤ 32` (`aws/campaign.json` `dpWeight` 32). That is the per-client
measurement the repository already relies on (`CAMPAIGN_ITER_PER_DP_LOG2`
in `aws/dp_ingest.py`, `benchmarks/dp-interval/`). The live feed agrees
within its noise: 3.65 B it/s over a 3,242 s window against about 39,500
points an hour gives 2^28.3.

Two numbers are easy to mistake for it:

- **2^27.93**, total iterations over stored points (2^56.38 / 367,123,410).
  It is low because the corpus still holds the points from the `HW(x) ≤ 34`
  era (2^25.27 steps a point). The repository retired this fleet-wide ratio
  once already.
- **2^28.34**, the mean length of a walk that *ended*. Walks cut by the
  `maxIters` guard count here, and that guard is finding F1 below.

The Rust counted-rho reference (`ic rho`) now reports it too: a
`steps / walk` column and `steps`, `walks`, `mean_steps_per_walk`, plus a
per-run `walks` and `steps_per_walk`. On a 16/20/24-bit ladder the tuned
walks sit on their design length (15.8, 34.2, 62.0 steps against 16, 32,
64), and the frozen walk runs 1.1–1.5× long. That excess is its
documented `θ²/r` share of walks that close a cycle and burn the cap.

### 2. Are distinguished points duplicated by the ingester?

**They are deduplicated, and the duplicates are counted.** The DP ingester
is `aws/dp_ingest.py`. There are three layers, and each one keeps a point
once per orbit key:

| layer | same orbit, same seed | same orbit, different seed |
|:--|:--|:--|
| client (`Solver::insert`, `include/solver.h`) | counted, not stored | collision, solved in process |
| ingester (`distinguished_points` primary key `(campaign_id, point_key)`, `ON CONFLICT DO NOTHING`; `dp_ingest_progress` makes re-ingesting an object a no-op) | counted in `dp_ingest_progress.duplicates`, published as `ingest.duplicate_records` | written to `rho_collisions`, page state `COLLISION_RECORDED` |
| `aws/merge.py` over the S3 corpus (the canonical detector) | dropped per bucket | reported and solved |

Live: 367,123,410 distinct points stored, **9,953,494 re-reports dropped
over the campaign's life, 0 in the last 24 hours**, and 0 collisions. The
~10 M are nearly all from 2026-09-19/20, when Modal runs 1–4 walked the
seeds of AWS slots 0–3 (see F3).

One blind spot, low priority: two records *in one upload object* with the
same orbit and different seeds are skipped by the ingester's collision
check (`findCollisions`, `worker == mine`). The client's own table sees
that pair first, since both came from one process, and `merge.py` sees it
regardless.

### 3. How often is the seed reset?

**Every time a walk ends, not on a clock.** A walk's seed is `run id (16)
‖ walk index (32) ‖ restart counter (16)`, and the kernel starts the
next walk on `seed + 1` (`packedkernels.cuh` `init`, `kernel.h`
`reseed`) in two cases:

- **The walk reaches a distinguished point.** This happens about every
  2^28.41 steps.
- **The `maxIters` guard fires.** The guard is checked every 4,096 steps
  and fires once a walk has gone 2^30 steps without a point. That happens
  to about 4.9% of walks.

In numbers:

- **Per walk.** On the audited RTX PRO 6000 preset (6,160,384 walks at
  about 14.6 B it/s, so about 2,370 steps/s per walk) a walk reseeds about
  every **40 hours**.
- **Fleet-wide, now.** 11.2 resets/s at a distinguished point plus about
  0.6/s at the guard, so about **11.7/s**. At the 2026-09-23 peak of
  842,000 points an hour it was about 245/s.
- **Scheduled restarts don't reseed.** The 48-hour process restart
  (`restartHours`) resumes from the checkpoint with the same run id and
  the same per-walk seeds and counters.
- **Run ids are permanent.** Each slot keeps its run id; retired slots
  keep theirs, so their seeds are never walked again.
- **The counter can't wrap.** It has 65,535 restarts of headroom, about
  300 years at 40 hours each, and the client exits 9 rather than wrap
  into another walk's seeds.

Nothing reseeds from fresh entropy, and that is deliberate. Every point is
replayed from its recorded 64-bit seed, and that replay is how a collision
gets solved (see Seeds).

### 4. How are fruitless cycles avoided?

**By construction, and now by a certificate that CI runs.**

- **Construction.** The walk `R ← σʲ(R) + R`, with
  `j = 3 + ((HW(x_R)/2) mod 8)`, is equivariant under `⟨σ, −1⟩`. It
  descends to the orbit set with no canonical representative in the loop,
  so there is no sign for a cycle to feed back through.
- **What is left to rule out.** Only an arithmetic coincidence:
  `∏(λ^{j_t} + 1) ≡ ±λⁱ (mod ℓ)` for some multiset of steps.
- **The certificate.** `ecc2k-guard cycles` rules that out exhaustively.
  At m = 131 **no multiset of up to 32 steps closes an orbit**
  (76,904,684 multisets, 2.96 × 10⁻²⁹ expected by chance). The Python
  check it replaces stopped at 8 steps.
- **Backstop.** The `maxIters` guard releases any walk that does not reach
  a point.
- **The alternative walk.** The table walk (`WALK_TABLE=1`) does leak
  fruitless cycles and is correctly not deployed (`WALK-CONSTANT.md` §5).

## What is costing money now

**F1. The guard throws away about 15.6% of the fleet's steps.** At
`dpWeight` 32 a trail averages 2^28.41 steps, so `maxIters` = 2^30 is only
2^1.59 trail lengths. The guard cuts 4.9% of honest trails before their
point, and those are the longest ones, so they hold 15.6% of all steps
(`WALK-CONSTANT.md` §6, a model on measured rates). No point is ever
recorded for those steps, so they cannot take part in a collision.
`WALK-CONSTANT.md` §9 recommends `maxIters ≥ 2^32` (a 0.01% loss).
`aws/campaign.json` still says 2^30. Raising it starts a new corpus,
because `maxIters` is part of the campaign contract, and
`src/witness.cpp --max-iters` has to move with it. That makes it the
campaign owner's decision, so this audit does not make the change.

**F2. The next fleet rebuild would feed the store points it cannot read.**
This is latent today. The Makefile defaults `WITNESS=1` for the σ walk
(since `0cb45d1e`, 2026-09-21), and `aws/build.sh` does not override it.
A rebuilt client therefore writes the 72-byte v2 corpus behind a 16-byte
`ECC2KDP2` header. `worker.py` and `merge.py` frame v2; two other paths do
not:

- **The strict upload marker.** `protocol.envelope` counts
  `size // 32` records and raises `partial DP record` whenever
  `16 + 72N` is not a multiple of 32, which is three deltas in four. The
  deltas that pass carry the wrong record count (N = 2 is recorded as 5).
- **The ingester.** `dp_ingest.py` decodes every object as 32-byte
  records. It would store header bytes as a "point" (the magic read as a
  seed), about 2.25 rows per real point, and none of those rows could
  deduplicate or collide.

The live slot is still on v1: its points arrive at the weight-32 interval
that 32-byte decoding predicts, and no slot is flagged off-weight. The fix
is either `WITNESS=0` in `build.sh`'s knobs, or v2 framing in
`protocol.py` and `dp_ingest.py`. The second one is best done as part of
the ingester's move to Rust.

**F3. Run-id separation between AWS and Modal is enforced on one side
only.** A run id names a seed space. Modal refuses run ids outside
8000–9999 (`modal_app.py checkCampaignRunId`). All three AWS allocators
(`DynamoSlots`, `S3Slots`, controlplane `PostgresSlots`) hand out
`max(slot) + 1` up to 65,534, with run id = slot + 1, and none of them
reserves the Modal range. The fleet has 195 slots, so the overlap is
7,804 slots away. This is the same kind of failure as the 9.95 M
duplicates, where two workers walked one seed space. The fix is a
reserved range in the AWS allocators, or a campaign-wide seed-space
registry.

## Formal validation of the cycle property

`ecc2k-guard cycles` decides a finite question exhaustively. It is not a
proof-assistant proof. For each supported curve it does four things:

1. **Derives ℓ.** It computes `#E` from the Frobenius trace recurrence
   (`V₀ = 2, V₁ = −1, V_{k+1} = −V_k − 2V_{k−1}`), requires `#E = 4ℓ`
   with ℓ prime (Miller–Rabin, 20 bases), and derives λ from `√−7`,
   taking the root of order dividing m.
2. **Cross-checks ℓ and λ at m = 131.** Both are compared with
   `ELL_DEC`/`S_DEC` in `generated/eccF131.h`. `codegen/gen.py` chose that
   λ by checking `[λ]P = σ(P)` on the challenge point itself, so the λ in
   the certificate is the one that actually acts on the campaign's points.
   The Python check never made this comparison.
3. **Checks both degeneracies.** `λʲ = 1` would make a step a doubling,
   and `λʲ = −1` would send it to infinity.
4. **Searches every multiset.** It walks every multiset of the eight step
   exponents up to the length bound and reports any product that lands on
   one of the 2m orbit scalars `±λⁱ`. The multipliers commute, so a walk
   that returns to its own orbit after k steps must produce exactly such a
   product. **A clean result therefore rules out every cycle of length ≤
   the bound, whichever path a walk takes.**

A hit counts as real when chance does not explain it (expected ≤ 0.01).
The certificate fails on any real hit, any degeneracy, or any mismatch
with the generated constants.

| m | multisets (≤ 8 steps) | expected by chance | result |
|--:|--:|--:|:--|
| 23 | 12,869 | 0.28 | clean |
| 41 | 12,869 | 1.9 × 10⁻⁶ | clean |
| 83 | 12,869 | 8.8 × 10⁻¹⁹ | clean |
| 97 | 12,869 | 6.3 × 10⁻²³ | clean |
| 131 | 12,869 | 5.0 × 10⁻³³ | clean, ℓ and λ equal the generated constants |
| **131, ≤ 32 steps** | **76,904,684** | **3.0 × 10⁻²⁹** | **clean** |

The L = 8 rows reproduce `codegen/cycles.py` exactly: the same ℓ, λ,
counts and expectations, pinned in the unit tests. Two controls show the
search can see cycles and classify them:

- **m = 19, a positive control.** Two of its eight multipliers share an
  orbit, and the search finds the same two length-5 cycles the Python
  found, classed as coincidence (3.7 expected).
- **m = 23 at 24 steps.** The search finds 229 closing multisets against
  231 expected. That is the rate random products give, which is what a
  correct search should report on a toy curve.

**What this does not cover.** It does not cover a trail that crosses
itself at random, which is not a fruitless cycle but the ordinary
self-intersection of a random mapping. Modelled as `T² · m / ℓ` with
`T = 2^28.41`, that is about 2^−65 per trail and about 10⁻¹⁰ over the
2^32.5 trails of the expected campaign. This is an estimate, not a
measurement, and the guard releases such a walk. The step exponents are
declared in the Rust source, mirroring `include/walk.h`; they are not
parsed from the kernel. The new CI job (`no-fruitless-cycles` in
`.github/workflows/ecc2k130-certification.yml`) runs on every change
under `ecc2k130/`, but a change to the schedule must also change
`STEP_EXPONENTS`.

## Seeds

**What the live client uses.**

- **The seed.** A walk seed is `run id (16) ‖ walk index (32) ‖ restart
  counter (16)`.
- **The start point.** A splitmix64 finaliser expands the seed into 128
  bits `c`, and the start point is `Q + Σ cᵢσⁱ(P)` (`include/walk.h`
  `eccSeedFor`, `eccPrf`).
- **Where it is stored.** Every point record carries its 64-bit seed, and
  every checkpoint carries each walk's seed and start iteration. So the
  live corpus already stores exactly what it used, and
  `ecc2k-guard seed-decode` now prints the run, walk, restart count and
  start bits for any recorded seed.

**Is it good enough?** For rho, the starts must be distinct orbits, must
not be correlated with the walk, and must be reproducible. The finaliser
is a bijection on the seed, so distinct seeds give distinct start
strings, and two starts sharing an orbit has a chance of about 2^−121
per pair. It is therefore adequate. It is still a mixer rather than a
keyed PRF, and nothing binds the seed space to a committed value anyone
can check.

**What `ecc2k-guard seed-generate` adds.**

- **A master seed.** It is 256 bits:
  `BLAKE3-derive_key("crypto ecc2k130 2026-09-25 rho campaign master seed
  v1", anchor ‖ entropy)`.
- **The anchor.** The anchor is the published Certicom instance (curve,
  field polynomial, ℓ, P, Q in polynomial basis), which nobody running the
  campaign chose.
- **The entropy.** The entropy is 32 bytes from the OS CSPRNG
  (`getrandom(2)`), or a recorded public beacon such as a block hash
  announced in advance.
- **The record.** The provenance file records every input, the anchor's
  hash and a commitment over the whole record, and it is written once:
  the tool refuses to overwrite it. `ecc2k-guard seed-verify` re-derives
  every field and fails on any change; a single flipped entropy bit fails
  three independent checks.
- **Test vectors.** The record carries vectors for a keyed-BLAKE3 start
  derivation (`walk_start_bits`). It keeps the same seed layout, so every
  point still replays from its 64-bit seed.

**Not wired into the client.** Changing how starts are derived changes
the corpus identity, like `maxIters`. The right time is the next corpus
boundary, and F1 already argues for one. At that point the master seed is
generated, committed in the pull request that opens the corpus, and the
kernel's `eccPrf` is replaced by `walk_start_bits` against the recorded
vectors.

**Why this construction.** It follows the same principle as the
verifiably random constructions in the challenge literature: every input
is either the published instance or recorded randomness, passed through
a standard KDF, so anyone can re-derive it. I did not find Bailey et
al.'s own start-point procedure in this repository, so none of it is
copied here. Grinding a rho seed buys nothing, because no seed makes a
collision come sooner. What a bad seed can do is overlap another seed
space (F3) or produce degenerate starts, and the construction prevents
both.

## Python → Rust

This PR converts `codegen/cycles.py`, which is deleted, and makes
`make check-cycles` call `ecc2k-guard`. It adds the seed tooling and the
`ic rho` walk statistics in Rust. About 35,000 lines of Python remain.
Rust rather than Go, because the repository is already a Rust crate with
a pinned toolchain and CI clippy. The order below keeps a parity
check in front of each live component:

| next | lines | why this order |
|:--|--:|:--|
| `aws/merge.py` | 398 | pure computation over files, and the canonical collision detector; byte-for-byte parity against the Python on a rehearsal corpus |
| `aws/dp_ingest.py` | 1,903 | fixes F2 as it moves; needs a Postgres parity harness first |
| `aws/status.py` | 152 | small; gains per-slot steps per walk |
| `aws/protocol.py` | 118 | with the ingester (v2 framing) |
| `aws/worker.py`, `aws/controlplane/` | 3,990 | last: they drive the live fleet, and a regression costs GPU hours |
| `scripts/rho_status/` | 2,896 | the public dashboard's publisher |
| Modal and RunPod launchers | 3,303 | their SDKs are Python-first; the thinnest shim that can call a Rust binary |
| `codegen/` | 14,261 | offline generators of the CUDA and C++ field code, not on the walk's path |

Rewriting the ingester or the supervisor mid-campaign without a parity
harness risks the failure recorded on 2026-09-17: five hours of points
invisible to the store because of one unrecognised key shape.

## Evidence

| file | sha256 |
|:--|:--|
| `live_status_20260925T050509Z.json` (public feed, generated 2026-09-25T05:05:09Z) | `4824dff5f5ce10a7323004d2165b2ba1e1a239d8fec23f36873bdefb07621a9f` |
| `no_fruitless_cycle_L8.json` | `35b4f8c4aedf19dc7d7283cbd93c31c016044bf573d099cad4c866efd35f0de9` |
| `no_fruitless_cycle_m131_L32.json` | `d526c62e803e3dac34d5e4cb7364befd4e6464ee128046f7a582f63467ed1c76` |

The `.txt` files beside each certificate are the tool's printed output.
To reproduce:

```
cargo test --lib ecc2k130_guard
cargo run --release --bin ecc2k-guard -- cycles --max-length 8 --json no_fruitless_cycle_L8.json
cargo run --release --bin ecc2k-guard -- cycles --curves 131 --max-length 32 --json no_fruitless_cycle_m131_L32.json
cargo run --release --bin ecc2k-guard -- seed-generate --out campaign-seed.json
cargo run --release --bin ecc2k-guard -- seed-verify campaign-seed.json
make -C ecc2k130 check-cycles CYCLE_LENGTH=24
```
