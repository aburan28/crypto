# Laned Frobenius-orbit keys: result

Declared in [`PROTOCOL.md`](PROTOCOL.md) (commit `45ebdcb6`) before the evidence
ran. The change is `a98ba95a`, measured against its parent `7424539b`.

**Class: engineering.** One key path is faster by a constant. Every key, every
relation, every logarithm and every verification is identical.

## Where the cost was

In `ic price` on the frozen §23 n61 parameter file, collection is 90% of the
cold pipeline. The pair-table build is most of the rest, and linear algebra is
0.06%.

Under Callgrind, which sees the portable path because Valgrind hides AVX-512,
`PairSumTable::keys_of` was 66.8% of all instructions. That is about 586
instructions per key over 74.4 M keys, from the exhaustive loop that tries all
`n − 1` rotations of eight words.

A floor run bounds the key's share natively. It used a deliberately wrong key
(`x + 1`) and was capped at the same 19 collection units and 66,576 trials,
3 runs. It showed the AVX-512 key costing about 0.85 s of 2.55 s collection.
On the portable path the key was about 60%.

The scalar search, which the portable bulk path did not use, has two problems
on blocks of keys:

- **Run-length exit.** It grows the zero-run length bit by bit and exits when
  growth stops. Run lengths spread over 3–9 bits, so nearly every key
  mispredicts that exit.
- **Candidate loop.** 27% of words have two or more longest runs, so the
  candidate loop mispredicts too. Each mispredict stops the out-of-order core
  overlapping the next independent key.

Distribution of 2·10⁶ random words (`n = 61`):

| longest-run candidates | share |
|:--|--:|
| one | 72.6% |
| two | 17.9% |
| three or more | 9.5% |

## The change

`FrobeniusCanon::least_rotations::<8>` returns the same least rotation in a
fixed number of steps per word:

1. `Z_{2^k}` by doubling.
2. The longest run by binary lifting, with selects instead of branches.
3. The two highest run tops compared with a select.
4. A loop only for a third run or more.

The eight lanes advance in lockstep, so their chains overlap.

- `canon_in_place` uses it on the portable path, so the scan keys and the
  folded build keys use it.
- `point_keys` now goes through `canon_in_place` everywhere.
- The scalar `least_rotation` is unchanged, and so are `canon_with_shift` (the
  rho walks) and the AVX-512 kernel.

## Results

All figures are for one target per size, under
`isolated_bench run --cpus 3 --max-other-cpu 0.3`, with `RAYON_NUM_THREADS=1`.
Each cell is the median of 5 A/B rounds. The raw reports are in `runs/ab/`, and
the medians, minima, maxima and quotients are `summary.json`, from
`summarise.sh`.

The reports were first written under their §23 run-directory names and renamed
to the ICV1 slugs before publication (AGENTS.md §11); no content changed. The
`isolation.jsonl` labels keep the names the runs were launched with.

- **Columns.** Times are ms. The ratio columns are baseline / candidate.
  "units" is the report's `total_units` (batched additions, measured in the
  same process).
- **S.** "S" is the batch pricer's `s_per_target` for this one target. It is
  not Table A's `S`, which is priced with `--single-target` over 64 targets on
  an older commit.

| size | path | collect, base | collect, cand | collect × | build, base | build, cand | total, base | total, cand | total × | units × | S | identical |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| icv1-f2m41-tm2308219-7f48b14a | portable | 94.8 | 85.6 | 1.108 | 21.1 | 19.3 | 132.9 | 121.7 | 1.092 | 1.093 | 16.37 → 14.98 | yes |
| icv1-f2m47-t22705043-f4e44623 | portable | 21.7 | 18.5 | 1.174 | 12.0 | 11.1 | 45.5 | 41.6 | 1.094 | 1.121 | 12.83 → 11.45 | yes |
| icv1-f2m53-tm56619371-dac20a85 | portable | 781.4 | 629.0 | 1.242 | 309.2 | 280.7 | 1160.4 | 985.1 | 1.178 | 1.229 | 22.99 → 18.70 | yes |
| icv1-f2m57-tm747311035-c1f545af | portable | 52.8 | 40.6 | 1.298 | 19.1 | 17.6 | 91.5 | 79.7 | 1.149 | 1.194 | 15.73 → 13.18 | yes |
| icv1-f2m59-tm943548413-98844ecc | portable | 868.5 | 626.4 | 1.387 | 324.8 | 287.8 | 1275.9 | 999.5 | 1.277 | 1.259 | 22.94 → 18.23 | yes |
| icv1-f2m61-t158598901-ab42b6c5 | portable | 4733.1 | 3507.3 | 1.349 | 534.6 | 485.3 | 5377.9 | 4106.7 | 1.310 | 1.311 | 38.18 → 29.12 | yes |
| icv1-f2m41-tm2308219-7f48b14a | avx512 | 57.8 | 57.7 | 1.002 | 15.8 | 16.5 | 89.4 | 90.6 | 0.986 | 0.993 | 10.99 → 11.06 | yes |
| icv1-f2m47-t22705043-f4e44623 | avx512 | 12.8 | 12.9 | 0.993 | 9.3 | 9.7 | 33.5 | 35.5 | 0.945 | 0.993 | 9.52 → 9.58 | yes |
| icv1-f2m53-tm56619371-dac20a85 | avx512 | 461.4 | 473.9 | 0.974 | 230.0 | 227.9 | 749.9 | 761.8 | 0.984 | 0.990 | 14.89 → 15.04 | yes |
| icv1-f2m57-tm747311035-c1f545af | avx512 | 29.3 | 29.5 | 0.992 | 14.7 | 14.7 | 62.3 | 62.8 | 0.993 | 1.014 | 10.76 → 10.62 | yes |
| icv1-f2m59-tm943548413-98844ecc | avx512 | 472.1 | 469.7 | 1.005 | 241.0 | 240.6 | 788.2 | 788.2 | 1.000 | 1.043 | 14.36 → 13.77 | yes |
| icv1-f2m61-t158598901-ab42b6c5 | avx512 | 2624.8 | 2609.4 | 1.006 | 406.3 | 414.0 | 3132.6 | 3121.5 | 1.004 | 1.020 | 22.29 → 21.85 | yes |

### A/A spread

The A/A spread is max/min − 1 over the baseline's 6 A/A runs:

| size | portable collect | portable total | avx512 collect | avx512 total |
|:--|--:|--:|--:|--:|
| icv1-f2m41-tm2308219-7f48b14a | 56.9% | 43.2% | 10.9% | 8.9% |
| icv1-f2m47-t22705043-f4e44623 | 10.3% | 8.8% | 14.2% | 8.2% |
| icv1-f2m53-tm56619371-dac20a85 | 14.1% | 10.0% | 8.4% | 6.5% |
| icv1-f2m57-tm747311035-c1f545af | 15.3% | 8.8% | 7.8% | 6.1% |
| icv1-f2m59-tm943548413-98844ecc | 6.4% | 4.8% | 4.5% | 5.7% |
| icv1-f2m61-t158598901-ab42b6c5 | 11.6% | 9.9% | 5.9% | 4.7% |

### Against the declared criteria

- **Portable.** At five of the six sizes, collection and total fall by more
  than the A/A spread. The n47 total, 1.094 against 8.8%, clears it narrowly.
  - On `icv1-f2m41-tm2308219-7f48b14a` **the gain is not established**: 1.108
    and 1.092 sit inside a 56.9% and 43.2% spread.
  - That spread comes from two slow baseline A/A runs on a 100 ms process,
    148 and 119 ms against 94–102.
  - The ratio grows with `n`, as the key's share of the scan does.
- **avx512 control.** Every size stays within its A/A spread, from 0.945 to
  1.006. The kernel this path runs is unchanged.
- **Identity.** Counts, recovered logarithms and verification are identical in
  all 12 cells and in every run.
- **Isolation.** 192 runs. One was contended (`aa/k1n47-avx512-baseline-r1`,
  other processes 0.11 CPU s), in the A/A set of the control, and is kept and
  flagged. None exited non-zero.
- **Refused attempt.** The first launch was refused before any run, because the
  session harness itself uses about 0.1 CPU. `runs/refused-attempt-1.log` keeps
  it. The allowance was then raised to 0.3 and recorded in `runs/host.json`.

### Instruction counts

These are Callgrind counts for the whole process, portable path, one run each.
They are deterministic. The reports are in `runs/callgrind/`.

| size | baseline Ir | candidate Ir | ratio | outputs |
|:--|--:|--:|--:|:--|
| icv1-f2m53-tm56619371-dac20a85 | 13,354,911,993 | 11,082,479,113 | 1.205 | identical |
| icv1-f2m61-t158598901-ab42b6c5 | 65,219,299,765 | 46,663,161,498 | 1.398 | identical |

`canon_in_place` is still 55% of the candidate's instructions at n61, about
340 instructions a key. The generated code holds the eight lanes in memory: it
uses 706 `mov`s and zeroes the level array with `memset` per block. Three
codegen rewrites were tried; none measurably beat it (below).

### Rho reference

The run was `ic price --single-target --rho-seed 7` on icv1-f2m61-t158598901-ab42b6c5, 3 interleaved
rounds per arm, isolated, with none contended (`runs/rho/`).

| arm | walk operations | rho online (ms) | IC online (ms) |
|:--|--:|:--|:--|
| baseline | 1,068,075 | 113.7, 114.8, 113.3 | 9.78, 9.97, 9.85 |
| candidate | 1,068,075 | 121.5, 111.2, 113.2 | 10.10, 9.65, 9.81 |

The rho walk calls the scalar search, which is unchanged.

## Tried and rejected

These were exploratory screens: n61, 5 interleaved rounds each, pinned with
`taskset` but not under `isolated_bench`. Each cell is the median collect in
ms. Every variant gave identical outputs.

| variant | path | what it was | baseline | variant | why rejected |
|:--|:--|:--|--:|--:|:--|
| A | portable | per key, existing scalar `canon` | 4617 | 4108 | superseded by the laned search |
| B | portable | A plus binary lifting in the scalar search | 4617 | 4186 | slower than A |
| C | portable | lifting plus two-candidate selects, per key | 4713 | 3958 | superseded by the laned search |
| F | AVX-512 | lifting kernel, scalar fix-up for 3+ candidates | 2528 | 2702 | slower; same-binary exhaustive control 2576 |
| G | AVX-512 | F plus three refinement steps | 2564 | 2479 | 3%, within noise; same-binary control 2528 |
| G | portable | laned plus refinement | 3397 (laned) | 3503 | slower than laned alone |
| K | portable | one word at a time, fully unrolled, mask selects | 3514 (laned) | 3843 | the lanes' overlap is worth more than registers |
| J | portable | lanes unrolled for `n ≥ 33`, mask selects | 3441 / 3453 (laned) | 3325 / 3433 | 3% in one screen and 0% in the next |

Two more variants were rejected for the scalar search itself, which the rho
walks call once per step on a serial chain:

| variant | rho online at n61, seed 7 (ms) |
|:--|--:|
| baseline | 113 |
| lifting plus two-candidate selects | 124 |
| lifting plus refinement | 128 |

Its latency is what counts there, so it stays as it was.

## Limits

- **AVX-512 hosts.** Measured, and no gain: the default path there is the
  AVX-512 kernel, untouched. This host and §23's host are that class.
- **Non-AVX-512 hosts.** Apple silicon, Graviton and AMD/consumer x86 run the
  portable path but were not measured. The portable figures here are this x86
  host with the AVX-512 path switched off, not a measurement of those hosts.
- **Field width.** Only fields with `n ≤ 63` are affected, the only ones
  `FrobeniusCanon` exists for. The m = 83 gate (§8a) and ECC2K-130 (m = 131)
  use wide fields and are untouched, so nothing here transfers to them.
- **Dashboards.** No IC/rho ratio is measured or changed. The scoreboard,
  leaderboard and progress chart cite §23's AVX-512-class measurements, and
  this round adds no reference arm on the portable path.
- **WDSat suite.** The frozen WDSat regression suite (§8) measures a SAT stage
  this change does not touch, so it was not run. The matched suite is the six
  frozen §23 parameter files.

## Reproduce

```sh
git show 7424539b:src/cryptanalysis/koblitz_fast.rs > src/cryptanalysis/koblitz_fast.rs
cargo build --release --bin ic --bin isolated_bench && cp target/release/ic /tmp/ic-baseline
git checkout src/cryptanalysis/koblitz_fast.rs
cargo build --release --bin ic && cp target/release/ic /tmp/ic-candidate
BASE=/tmp/ic-baseline CAND=/tmp/ic-candidate BENCH=target/release/isolated_bench \
  MAX_OTHER_CPU=0.3 bash research/canon_laned_keys_20261008/run.sh   # refuses to overwrite runs/
bash research/canon_laned_keys_20261008/summarise.sh
```
