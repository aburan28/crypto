# P-256 S17 batched collision selector, round 20: result

Date run: 2026-10-05

**EXACT SAMPLED-LIST REPLAY PASSES; FINITE-DOMAIN ACCOUNTING REJECTS GENERIC
BATCHING 16.203 BITS ABOVE RHO.**  The measured width ladder is unchanged, but
the preregistered symmetric projection is not attainable: it asks for
`2^136.538` distinct signed-eight images while the entire left domain contains
only `2^128.734`.  Repeating left representations cannot create independent
rows.  The initially reported `2^144.625` direct cost and `2^141.625` generic
boundary are therefore withdrawn, not promoted.

The corrected asymmetric schedule materialises the complete left domain and
streams target-adjusted signed-nine images.  After charging the exact
finite-domain relation population and duplicate occupancy, direct batching
costs `2^151.885` P-256 field-multiplication equivalents (FME), 19.419 bits
above Pollard rho.  Even an optimistic boundary charging one group addition
per emitted image costs `2^148.668` FME, 16.203 bits above rho.  This remains a
13.271-bit projected improvement over round 19, but it is not parity.

## One boundary table, one unit

One FME is one P-256 field multiplication.  A P-256 group addition is charged
as 17 FME.  The batch rows include 138,031 independent relation rows and the
frozen sparse-linear-algebra model.  Sorting and byte movement remain
unpriced, so both batch figures are lower bounds.

| variant | total FME, log2 | log2(total / rho) | per usable row, log2 FME | materialised bytes, log2 | correctness | class |
|:--|--:|--:|--:|--:|:--:|:--|
| Pollard rho | 132.466 | 0 | n/a | small | generic reference | reference |
| round-19 exhaustive variable 8+9 | 165.156 | 32.690 | 148.081 | 148.926 structured frontier | exact prefixes | superseded projection |
| withdrawn symmetric batching | ~~144.625~~ | ~~12.159~~ | ~~127.551~~ | ~~143.860~~ | exceeds finite left domain | withdrawn accounting |
| direct finite-domain signed 8+9 | **151.885** | **19.419** | **134.811** | **135.056** left table | exact sampled lists | projected engineering |
| one-addition-per-image finite-domain boundary | **148.668** | **16.203** | **131.594** | 135.056 left table | generic lower boundary | boundary |

The corrected direct candidate remains 31.885 bits above the `2^120`
collection gate, 31.811 bits above the `2^103` per-row gate, and 85.056 bits
above the `2^50` materialisation gate.  Sparse linear algebra contributes only
634,220,698,496 optimistic FME and does not affect the displayed exponent.

## Measured P-256 image ladder

Each cell uses actual points from `FB1h2f8621cda105`, 256 deterministic
hash-selected unplanted targets, and one separately labelled planted control.
`N` left signed-eight and `N` target-adjusted signed-nine images are generated,
where `N = 2^(w/2+3)`.  The expected number of random projected-key matches is
64 in every cell.  Every projected match is checked against the complete
257-bit compressed point key and every full-key match receives fresh group
replay.

| projected bits | N per side | projected matches | random matches / 64 | full-key matches | unplanted relations | FME | logical list bytes | unranked wall |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 20 | 8,192 | 72 | 71 / 64 | 1 planted | 0 | 2,311,203 | 1,310,720 | 361 ms |
| 22 | 16,384 | 65 | 64 / 64 | 1 planted | 0 | 4,622,113 | 2,621,440 | 342 ms |
| 24 | 32,768 | 57 | 56 / 64 | 1 planted | 0 | 9,243,933 | 5,242,880 | 638 ms |
| 26 | 65,536 | 57 | 56 / 64 | 1 planted | 0 | 18,487,573 | 10,485,760 | 1,325 ms |
| 28 | 131,072 | 58 | 57 / 64 | 1 planted | 0 | 36,974,853 | 20,971,520 | 2,637 ms |

The five cells contain 304 random projected matches versus 320 expected, a
ratio of 0.95.  All 304 fail the complete-key comparison and are retained as
truncated screening candidates, not false relations.  The five forced planted
matches are the only full-key matches and all replay exactly.

The counted-work slope is 0.499979 FME bits per projected key bit; logical
memory is exactly 0.5.  Both agree with the frozen birthday-ladder prediction.
Peak RSS was 208,617,472 bytes, dominated by the complete rebuilt factor base,
and disk traffic was zero.  Target generation was separately charged at
32,891 projective additions and 65,064 doublings.

## Exactness boundary

- Every cell's independently full-key-sorted reference returns exactly the
  selector's relation set.
- The 20-bit cell additionally compares all 67,108,864 left/right pairs and
  returns the same set.
- False positives / false negatives: `0 / 0` in every cell.
- Duplicate relations: zero.
- Every planted relation uses 17 distinct columns and passes fresh P-256 group
  replay.
- No unplanted exact relation was observed.  That is expected at these sample
  counts and is not a full-domain nonexistence claim.

Completeness applies only to each frozen pair of sampled lists.  Projected-key
collisions are never treated as exact group collisions, and omitted samples
are never described as searched branches.

## Full-size projection

For `B = 131458`, the exact probability that independently drawn distinct
eight- and nine-column samples are disjoint is
`0.9994524131979801`.  The complete signed representation domains are

```text
L = C(B,8) * 2^8 = 2^128.734424
R = C(B,9) * 2^9 = 2^143.568654.
```

The original balanced birthday formula requested `2^136.537712` images on
each side.  Its left request exceeds `L` by 7.803 bits and is impossible
without repeats, so its `2^144.625` / `2^141.625` projections are withdrawn.
This was found by checking the projection against the exact domains after the
first ladder run; `ACCOUNTING_CORRECTION.md` freezes the replacement.

The exact signed-17 mean is 3.324241 relations per target, corresponding to
the inherited 0.964000182 success probability.  The corrected projection
allows 143,185.658 targets for 138,031 usable rows, giving 475,983.691
distinct relations in expectation.  Uniform occupancy then requires
163,013.787 collision events, or 1.180994 events per usable row.  Materialise
all `L` left images and stream

```text
N_right = collision_events*n/(L*p_disjoint) = 2^144.581001 images.
```

This is only `2^127.454` right images per allowed target, below the complete
signed-nine domain.  Seven additions per left image, nine per streamed right
image, measured batch normalisation, target setup, relation replay, and sparse
linear algebra total `2^151.885227` FME.  Charging only one addition per image
and giving key extraction and replay away produces the `2^148.668488` generic
lower boundary.  Sorting, comparison, and byte movement remain unpriced.

The complete 80-byte left table alone needs `2^135.056352` bytes.  Streaming
the right side avoids materialising it but cannot pass the storage gate.  A
distinguished-point or external-memory schedule may trade storage for replay
and I/O, but cannot lower the generic image count; no unimplemented schedule
is credited.

## Dickson and degree evidence

The native harness rebuilds all 131,458 factor-base columns and recomputes all
18 Dickson steps for every leaf.  Rejections are zero and survival remains
1.0.  As in round 19, Dickson membership is a separable leaf property; it does
not alter the generic collision probability.

The frozen degree receipt remains:

- local split join degree: 2;
- structured residual maxima at depths 1, 2, 3: 3, 3, 4;
- unsplit S17 degree of regularity: unknown.

Degree at most 5 passes, but local degree does not close the 16.203-bit generic
floor gap.

## Decision and next target

Do not run the full-depth unplanted P-256 collision search.  Exactness and the
degree gate pass; cost, per-row, storage, and rho-parity gates fail.

The next credible promotion target is now numeric: demonstrate and prove an
exact cross-column invariant or representation reuse worth more than
`2^16.203` over the one-addition-per-image boundary, and more than `2^19.419`
over direct random signed sums, without discarding probabilistic branches.
Bucket rearrangements, more local degree reductions, and memory-only
distinguished-point schedules cannot meet that target.

## Reproduction and custody

Source commit used for the run:
`4aa8ae7268d8b5db248ba3ee8e289e00613da786`.

```bash
cargo test --bin p256_s17_batch_collision
cargo clippy --bin p256_s17_batch_collision -- -D warnings \
  -A clippy::mismatched_bit_width_type
cargo build --release --bin p256_s17_batch_collision
target/release/p256_s17_batch_collision \
  --round19 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_multilevel_selector_round19_20261005/selector-result.json \
  --round6 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_batch_collision_round20_20261005/collision-result.json
```

The run used Rust 1.99.0 (`b940084d7`), Linux x86-64, five exposed cores of
an AMD EPYC 9V74 VM, and made no wall-time performance claim.
`collision-result.json` is 31,729 bytes with SHA-256
`72f269c401d5192fda1e2a9f8c7e4b7bd537d92fe29a18b851882b610fa5c459`.
