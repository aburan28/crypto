# P-256 S17 batched collision selector, round 20: result

Date run: 2026-10-05

**EXACT SAMPLED-LIST REPLAY PASSES; GENERIC BATCHING RECOVERS 20.531 BITS,
BUT PARITY STILL FAILS.**  Sharing one two-sided birthday process across
138,031 relation rows lowers round 19's optimistic projection from
`2^165.156` to `2^144.625` P-256 field-multiplication equivalents (FME).
That is a real projected algorithmic improvement over one exhaustive 8+9
frontier per target.  It remains 12.159 bits above Pollard rho.

More importantly, an optimistic boundary that emits every new signed-sum
image for only one group addition still costs `2^141.625` FME, 9.159 bits
above rho.  Generic batching alone therefore cannot reach parity.  A successor
needs a non-generic Dickson/curve invariant worth more than `2^9.159`, or a
mechanism that obtains multiple independent rows from one generic collision.

## One boundary table, one unit

One FME is one P-256 field multiplication.  A P-256 group addition is charged
as 17 FME.  The batch rows include 138,031 independent relation rows and the
frozen sparse-linear-algebra model.  Sorting and byte movement remain
unpriced, so both batch figures are lower bounds.

| variant | total FME, log2 | log2(total / rho) | per usable row, log2 FME | materialised bytes, log2 | correctness | class |
|:--|--:|--:|--:|--:|:--:|:--|
| Pollard rho | 132.466 | 0 | n/a | small | generic reference | reference |
| round-19 exhaustive variable 8+9 | 165.156 | 32.690 | 148.081 | 148.926 structured frontier | exact prefixes | superseded projection |
| direct batched random signed 8+9 | **144.625** | **12.159** | **127.551** | **143.860** | exact sampled lists | projected algorithmic advance |
| one-addition-per-image boundary | **141.625** | **9.159** | **124.551** | schedule-dependent | generic lower boundary | boundary |

The direct candidate remains 24.625 bits above the `2^120` collection gate,
24.551 bits above the `2^103` per-row gate, and 93.860 bits above the `2^50`
direct-materialisation gate.  Sparse linear algebra contributes only
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
| 20 | 8,192 | 72 | 71 / 64 | 1 planted | 0 | 2,311,203 | 1,310,720 | 354 ms |
| 22 | 16,384 | 65 | 64 / 64 | 1 planted | 0 | 4,622,113 | 2,621,440 | 326 ms |
| 24 | 32,768 | 57 | 56 / 64 | 1 planted | 0 | 9,243,933 | 5,242,880 | 687 ms |
| 26 | 65,536 | 57 | 56 / 64 | 1 planted | 0 | 18,487,573 | 10,485,760 | 1,362 ms |
| 28 | 131,072 | 58 | 57 / 64 | 1 planted | 0 | 36,974,853 | 20,971,520 | 2,642 ms |

The five cells contain 304 random projected matches versus 320 expected, a
ratio of 0.95.  All 304 fail the complete-key comparison and are retained as
truncated screening candidates, not false relations.  The five forced planted
matches are the only full-key matches and all replay exactly.

The counted-work slope is 0.499979 FME bits per projected key bit; logical
memory is exactly 0.5.  Both agree with the frozen birthday-ladder prediction.
Peak RSS was 208,498,688 bytes, dominated by the complete rebuilt factor base,
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
`0.9994524131979801`.  To obtain `R = 138031` useful cross-side collisions in
the P-256 group requires

```text
N_left = N_right = sqrt(R*n/p_disjoint) = 2^136.537712 images.
```

The directly measured sampling rule uses seven additions per left image and
nine per right image, giving `16*N = 2^140.537712` group additions and
`2^144.625174` FME.  The optimistic generic boundary replaces all sixteen
additions by only two, one per emitted image, and still costs
`2^141.625174` FME.

Directly retaining both lists as the measured 80-byte records needs
`2^143.859640` bytes.  A distinguished-point or external-memory schedule may
trade storage for replay and I/O, but cannot lower the generic image count.
No such unimplemented schedule is credited with passing the storage gate.

## Dickson and degree evidence

The native harness rebuilds all 131,458 factor-base columns and recomputes all
18 Dickson steps for every leaf.  Rejections are zero and survival remains
1.0.  As in round 19, Dickson membership is a separable leaf property; it does
not alter the generic collision probability.

The frozen degree receipt remains:

- local split join degree: 2;
- structured residual maxima at depths 1, 2, 3: 3, 3, 4;
- unsplit S17 degree of regularity: unknown.

Degree at most 5 passes, but local degree does not close the 9.159-bit generic
floor gap.

## Decision and next target

Do not run the full-depth unplanted P-256 collision search.  Exactness and the
degree gate pass; cost, per-row, storage, and rho-parity gates fail.

The next credible promotion target is now numeric: demonstrate and prove an
exact cross-column invariant or representation reuse worth at least
`2^9.159` over the one-addition-per-image boundary, and at least `2^12.159`
over direct random signed sums, without discarding probabilistic branches.
Bucket rearrangements, more local degree reductions, and memory-only
distinguished-point schedules cannot meet that target.

## Reproduction and custody

Source commit used for the run:
`1aa5d4ad4cfe89cd406c4ce53ba989fe285e363a`.

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
`collision-result.json` is 31,060 bytes with SHA-256
`503f80a84b8af9e3690d78470cae058641f2626e3c753348aceed623c2cf6bc9`.
