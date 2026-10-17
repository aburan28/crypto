# A15: each scanned key's filter slot hashed once a block

**What this is.** A timed exploration on v8, R12's accepted baseline. It ran on 2026-10-17 from 03:50 to 04:08 UTC on
the Emerald Rapids host R12 ran on. It is a record, not a round's
measurement; none of its runs is pooled with any round's.

**Why.** The stage probes at v8 ([record](../scan-stages-v8-20261017/README.md))
put the filter stage at 5.0–6.6 ns a scanned summand, the largest of the four.
In v8 the filter loop hashes each key twice:
- once to prefetch its filter word 32 keys ahead;
- once more to test it.

The hash is two 64-bit multiplies.

**The outcome: negative, no round.** Hashing once a block into two slot
arrays makes cold time 2–4% *slower* at all three target sizes, and
collection 2–6.5% slower. In v8 the second hash overlaps the filter's
memory wait, so removing it saves nothing measurable. The extra pass and its
16 KiB of slot arrays cost more than the hash they replace.

## The candidate

`7ec4e6cb`, one commit on v8 (`23138d20`), frozen as
[`filter-slots.patch`](filter-slots.patch):
- Before the filter loop, one pass fills `slot_words` and `slot_bits` with
  each key's `filter_slot(pair_filter_hash(key), …)`.
- The prefetch ahead and the test read those slots.
- The values tested are the same; only where they are computed changes.

The tests were not run on its tree, because no round was declared.

## What ran

- **The rows:** `M1`'s two rows at the three target sizes:
  - `icv1-f2m53-tm56619371-dac20a85` (2^44.3);
  - `icv1-f2m59-tm943548413-98844ecc` (2^44.5);
  - `icv1-f2m61-t158598901-ab42b6c5` (2^47.2).
- **The arms:** `v8` (`ic-r12-v7-23138d20`) and `a15`
  (`ic-a15-v8-7ec4e6cb`), hashed in [`binaries.sha256`](binaries.sha256).
- **The schedule:** six rounds with the order alternating, 72 processes, all completed:
  - each process ran `ic price --single-target` with the suite's rho seed and `RAYON_NUM_THREADS=1`;
  - each was isolated on CPU 2
    ([`explore.sh`](explore.sh), [`explore.log`](explore.log)).
- **Refusals:** the busy gate refused a run 39 times before it started. Each
  refused run waited 15 s and was asked for again.
- **Contended runs:** one run was marked contended (`a15`, 2^44.3) and is left out.
- **Logarithms:** every pair recovered the same logarithm in both arms, all
  36 ([`pairs.tsv`](pairs.tsv)).

## Results

Cold time, `v8` over `a15`; above 1 means the candidate is faster. The
intervals are 95% t intervals of the geometric mean over the clean pairs
([`summarise.sh`](summarise.sh), [`summary.tsv`](summary.tsv)):

| curve | `log₂ r` | clean pairs | cold time | 95% interval | all pairs |
|:--|--:|--:|--:|:--|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 11 | **0.958** | [0.911, 1.007] | 0.957 |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 12 | **0.964** | [0.935, 0.993] | 0.964 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 12 | **0.983** | [0.967, 1.000] | 0.983 |

The phases are the geometric mean over all 12 pairs of each pair's ratio of median repetitions,
`v8` over `a15` ([`phases.sh`](phases.sh), [`phases.tsv`](phases.tsv)):

| curve | collection | build | selection |
|:--|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 0.935 | 0.985 | 1.045 |
| `icv1-f2m59-tm943548413-98844ecc` | 0.944 | 0.999 | 1.018 |
| `icv1-f2m61-t158598901-ab42b6c5` | 0.980 | 0.997 | 1.050 |

The median cold time over each arm's 12 processes:
- `v8`: 412.1, 422.7 and 1647.7 ms;
- `a15`: 423.7, 441.4 and 1670.0 ms.

## What it shows

- **The filter stage is not bound by the hash.** Its 5–6.6 ns is mostly the
  filter word's memory access. That is about 11–14 cycles at 2.1 GHz with 32
  keys of lookahead.
- **Splitting the loop costs more than it saves.** The split adds a pass
  over the block and 16 KiB of stores and loads, and it moves the hashes out from
  under the memory wait they used to fill.

## Graph and report

- **Graph:** [`report/cold-ratio.svg`](report/cold-ratio.svg).
- **Report:** [`report/report.pdf`](report/report.pdf), with its source
  [`report/report.typ`](report/report.typ).
- **Copies:** in
  [`benchmarks/20261017-ic-a15-filter-slots/`](../../../../benchmarks/20261017-ic-a15-filter-slots/).

## Scope

One host and one core, the three target sizes, two rows a size.
