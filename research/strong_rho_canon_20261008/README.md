# Strong-rho canonicalization and batched-inversion scratch: result

[`PROTOCOL.md`](PROTOCOL.md) (commit `4c243e69`) was declared before the
evidence ran. It measures two changes, `23ed8d3a` and `a8f4bb60`, against their
parent `7ce52fa2`.

**Class: engineering.** Every canonical form, walk, relation, logarithm and
verification is identical across the arms. The strong signed-Frobenius rho is
1.40× faster per capped walk at m = 83 and 1.46–1.61× faster at n = 41–61. No
native time effect on the index-calculus pipeline is established, in either
direction.

## The changes

- **`23ed8d3a`: canonicalization from the longest zero runs.**
  `StrongRhoG::canonicalize` used to try all `n − 1` rotations of the
  normal-basis `x` coordinate, plus a negation at each. That is two tuple
  compares per rotation, on `u128` words at m = 83.
  - The least rotation starts with a longest run of zero coordinates. So the
    new code finds the longest run's length by AND-ing rotations of the
    complement, then compares only the rotations that start at those runs.
  - The sign is then read off one `y` rotation.
  - A tie happens only for a periodic `x`, which includes 0 and all-ones. A tie
    still takes the old exhaustive search, so the sign rule is unchanged.
  - A new unit test checks the search against the exhaustive scan for
    n = 2–127.
- **`a8f4bb60`: no scratch zero-fill per call.**
  - `Gf2::batch_inv_clmul` and `batch_inv_inline` used to clear their scratch
    and zero-fill it to the batch length on every call. They now truncate and
    resize, so only a batch longer than the last one zero-fills, and only its
    new tail.
  - The AVX-512 `add_many_lazy` grows its accumulator only when it is too
    short.
  - The scratch's post-state is unchanged, which
    `gf2_matches_b072fcf5_code_on_every_operation` pins.

## Results

This host is one x86-64 cloud container with AVX-512 and 4 logical cores.
Binary hashes and the toolchain are in [`runs/host.json`](runs/host.json).

**Isolation.** 176 declared processes ran under `isolated_bench --cpus 3
--max-other-cpu 0.3` with `RAYON_NUM_THREADS=1`. None was contended and none
exited nonzero. `summary.json` holds every figure below; `summarise.sh`
derives it from `runs/`.

### Identity gate: passed

- **M1.** One canonical digest across all 16 runs of both arms, and every
  walk reached its cap.
- **M2.** Every output line is identical once the timing fields are set aside:
  charges, restarts, fruitless cycles and the recovered scalar.
- **M3.** `counts`, `recovered` and `all_verified` are identical in every run
  of both arms. They are also identical in the post-hoc runs and the Callgrind
  runs below.

**Note on the M2 rule.** The first pass of `summarise.sh` dropped only
`peak_rss_bytes` and reported M2 as not identical. The fields that differed
were `setup_ms`, `target_generation_ms`, `total_ms`, `validation_ms` and
`walk_ms`. They differ between two runs of the same binary in the A/A pass as
well. The protocol's rule is "identical output lines, timing aside", so
`summarise.sh` now drops every `*_ms` key. No other field differed in any run.

### M1: m = 83 reference

Curve `icv1-f2m83-tm6151469093347-debefd74`. Medians of 5 interleaved A/B
rounds, with the range in brackets; the A/A spread is max/min − 1 over 6
baseline runs.

| measurement | baseline | candidate | ratio | A/A spread |
|:--|--:|--:|--:|--:|
| canonicalize + hash, ns per state | 486 [454–507] | 220 [211–281] | 2.21× | 3.1% |
| walk capped at 10^6 steps, ms | 981 [931–1018] | 701 [686–733] | 1.40× | 4.9% |

- **Against the declared criterion.** Met: the capped walk is 1.40× faster,
  well outside the 4.9% A/A spread, and the two ranges do not overlap.
- The ns-per-state figure includes hashing each result with blake3, which
  both arms pay equally. A shared cost lowers a ratio above 1, so the
  canonicalization's own ratio is above 2.21×.
- **Observation.** The walk saves about 280 ns per step and canonicalization
  saves about 266 ns per state. That fits the walk's saving being
  canonicalization, but this round does not separate it further.

### M2: narrow strong rho

`koblitz_rho_fixture n a signed_frobenius 4 strong`: 4 fixtures per process,
each solved to the end. The columns are:

- wall: process wall time from the isolation record, in ms; median [range] of
  5 A/B rounds;
- spread: the A/A spread of the process wall over 6 baseline runs;
- walk ratio: the ratio of the summed `walk_ms` of the four fixtures, with its
  A/A spread in brackets.

| curve | baseline wall | candidate wall | ratio | A/A spread | walk ratio |
|:--|--:|--:|--:|--:|--:|
| `icv1-f2m41-tm2308219-7f48b14a` | 68 [62–71] | 47 [45–59] | 1.46× | 23.2% | 1.50× [22.5%] |
| `icv1-f2m47-t22705043-f4e44623` | 41 [40–58] | 27 [25–32] | 1.55× | 47.4% | 1.66× [45.6%] |
| `icv1-f2m53-tm56619371-dac20a85` | 645 [626–662] | 412 [405–421] | 1.57× | 2.9% | 1.58× [2.8%] |
| `icv1-f2m57-tm747311035-c1f545af` | 65 [63–77] | 41 [37–58] | 1.58× | 3.3% | 1.73× [3.4%] |
| `icv1-f2m59-tm943548413-98844ecc` | 555 [501–619] | 345 [323–365] | 1.61× | 8.5% | 1.63× [8.5%] |
| `icv1-f2m61-t158598901-ab42b6c5` | 810 [785–823] | 506 [474–511] | 1.60× | 4.6% | 1.62× [4.6%] |

- **Against the declared criterion.** Every gain is outside its A/A spread.
  - The n = 41 and n = 47 processes last 40–70 ms, so start-up noise sets
    their wide spreads.
  - At n = 53, 59 and 61 the walk dominates and the spread is 3–9%.
- Both changes act here. This round does not measure how the gain splits
  between them.

### M3: index-calculus pipeline

`ic price --repeats 1 --repeats-fast 1` on the frozen §23 parameter files,
single thread. Medians of 5 A/B rounds, in ms; ratio = baseline / candidate.

| curve | path | collect baseline | collect candidate | collect ratio | A/A spread | total ratio | A/A spread |
|:--|:--|--:|--:|--:|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | AVX-512 | 445 [438–462] | 437 [429–450] | 1.017 | 6.8% | 1.044 | 6.6% |
| `icv1-f2m53-tm56619371-dac20a85` | portable | 594 [564–610] | 587 [576–618] | 1.013 | 5.9% | 1.028 | 5.5% |
| `icv1-f2m61-t158598901-ab42b6c5` | AVX-512 | 2458 [2394–2552] | 2465 [2375–2566] | 0.997 | 6.0% | 0.996 | 3.8% |
| `icv1-f2m61-t158598901-ab42b6c5` | portable | 3235 [3199–3415] | 3346 [3203–3438] | **0.967** | 2.6% | 0.978 | 2.5% |

**Against the declared criterion:**

- Three cells are inside their A/A spreads. No gain is established in them,
  and no loss.
- **Portable n61 collection is a regression by the declared rule.** The
  candidate's median is 3.3% slower, against a 2.6% A/A spread. Its total is
  inside its spread (2.2% against 2.5%).

### M3 follow-ups (post hoc, not declared)

Both follow-ups ran after the declared runs, because of the portable n61
cell. [`posthoc.sh`](posthoc.sh) has their exact commands.

**Instruction counts.** Callgrind, portable n61, one run per arm, one price
round:

| | baseline | candidate | change |
|:--|--:|--:|--:|
| all instructions | 46,661,818,722 | 46,030,572,942 | −631 M (−1.35%) |
| `memset` | 1,060,989,932 | 432,911,463 | −628 M |
| `memset` called from `batch_inv_clmul` | 687,598,351 | 59,519,699 | −628 M |
| `batch_inv_clmul` itself | 5,062,174,955 | 5,060,238,107 | −2 M |

Nearly the whole instruction saving is the zero-fill `a8f4bb60` removes.

**Three-arm screen.** Portable n61, 5 interleaved rounds of three arms, all 15
processes uncontended, outputs identical. The `item1` arm is `23ed8d3a` alone.
It runs the same index-calculus source as the baseline: that commit changes
only `koblitz_strong_rho.rs`, which `ic price` does not call without
`--rho-seed`.

| arm | collect ms, median [range] | total ms, median |
|:--|--:|--:|
| baseline (`7ce52fa2`) | 3265 [3218–3409] | 3807 |
| item1 (`23ed8d3a`) | 3358 [3086–3411] | 3931 |
| candidate (`a8f4bb60`) | 3193 [3118–3379] | 3752 |

**What was observed:**

- The declared loss did not repeat. In this screen the candidate had the
  lowest median, and baseline / candidate = 1.023.
- Two builds of the same index-calculus source, baseline and item1, differ in
  median by 2.8%. That is about the size of the declared loss.
- All three ranges overlap.

**Hypothesis, not tested:** the portable n61 cell measures build-to-build code
layout, which an A/A pass of one binary cannot capture. No layout-randomized
build was run.

**Status of `a8f4bb60`:**

- The declared verdict on portable n61 stays as recorded: a regression.
- Taken with the screen, no native time effect of `a8f4bb60` is established,
  in either direction.
- What is established is 1.35% fewer instructions.
- The change was kept, because it removes work and changes no output.

## What this changes downstream

- **Counted-unit records do not change.**
  `docs/bounds/records/koblitz-rho-strong.json`, its pair-claw sibling and
  the ecbench frontier count `ecbench.gae` units, not time. Speed does not
  enter them.
- **Wall-time and instruction comparisons against the strong rho get harder
  for index calculus.**
  - Index-calculus/rho ratios measured from now on against
    `koblitz_rho_fixture … strong` price rho with the new canonicalization,
    which M2 shows is 1.5–1.6× faster per process at n = 53–61.
  - Committed sessions keep their numbers. A dashboard row that cites a
    strong-rho wall or instruction count measured before `23ed8d3a` prices rho
    with the slower canonicalization.
  - This round re-measures none of those rows.
- **The admissibility bar rises.** BOUNDARY_TARGETS admits a different rho only
  at a per-step cost no worse than the strong one, and the strong one is now
  cheaper.
  - The batch references `koblitz_rho_batch_ks_strong` and
    `koblitz_rho_batch_ks_v3` carry their own exhaustive rotation scans, so
    `23ed8d3a` does not reach them. Whether each still meets the bar is not
    measured here.
- **Table A's rho reference gets only the second change.** `ParallelRho`
  canonicalizes through `FrobeniusCanon::canon_with_shift`, so `23ed8d3a` does
  not reach it. It does invert through `Gf2::batch_inv`, via
  `FastCurve::add_many`, so `a8f4bb60` acts on it. Its time was not measured
  here; M2 is a different walk.

## Limits

- **m = 83 is the reference only.** There is no pair-table index calculus at
  m = 83, and M1 measures a walk capped at 10^6 steps, not a solve.
- **One host class.** No Arm64, AMD or GPU host was measured, and the GPU
  ECC2K-130 rho is separate code. The portable M3 figures are this host with
  AVX-512 switched off, not a measurement of a host without it.
- **Leftover zero-fill.** The candidate still spends 433 M instructions in
  `memset`:
  - 362 M are `FrobeniusCanon::canon_in_place` zeroing a 16-word array on
    each of 10 M calls;
  - 60 M are `batch_inv_clmul` growing its scratch over 1,959 calls.

  Neither is changed here.
- **No split of M2.** The gain at n ≤ 61 is attributed to both changes
  together.

## Reproduce

Build the three arms (the baseline is `7ce52fa2` plus this round's
measurement example), then run the declared measurements, the follow-ups and
the summary:

```sh
for arm in baseline:7ce52fa2 item1:23ed8d3a candidate:a8f4bb60; do
  name=${arm%%:*}; rev=${arm#*:}; w=/tmp/arm-$name
  git worktree add --detach "$w" "$rev"
  git show a8f4bb60:examples/strong_rho_step_rate.rs > "$w/examples/strong_rho_step_rate.rs"
  (cd "$w" && cargo build --release --bin ic --example koblitz_rho_fixture --example strong_rho_step_rate)
  mkdir -p /tmp/arms
  cp "$w/target/release/ic" /tmp/arms/ic-$name
  cp "$w/target/release/examples/koblitz_rho_fixture" /tmp/arms/koblitz_rho_fixture-$name
  cp "$w/target/release/examples/strong_rho_step_rate" /tmp/arms/strong_rho_step_rate-$name
done
cargo build --release --bin isolated_bench
BIN=/tmp/arms BENCH=target/release/isolated_bench MAX_OTHER_CPU=0.3 \
  bash research/strong_rho_canon_20261008/run.sh          # declared M1-M3; skips files already in runs/
BIN=/tmp/arms CG=/tmp bash research/strong_rho_canon_20261008/posthoc.sh callgrind
BIN=/tmp/arms BENCH=target/release/isolated_bench \
  bash research/strong_rho_canon_20261008/posthoc.sh screen
bash research/strong_rho_canon_20261008/summarise.sh
```
