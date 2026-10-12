# R02b's go/no-go on v3: the wide-tail kernel against main's reduction

**What this is.** The exploration that R02b's amendment 3 declared
([protocol](../../rounds/R02b-wide-tail-retest/PROTOCOL.md)), and its
decision. It ran on 2026-10-06 from 05:40 to 05:52 UTC, after R07's last
step, on the newest baseline, v3 (main's head `995ea207`). The base ran
against the base plus R02's
[`candidate.patch`](../../rounds/R02-wide-tail-kernel/candidate.patch).
Its rule decides whether R02b runs or is withdrawn. It is not a round's
measurement, and none of its runs is pooled with any round's.

**The decision: R02b is withdrawn, and the kernel is retired.** On v3
the candidate is slower than the base at both target sizes.

## The arms

| arm | binary | SHA-256 | tree |
|:--|:--|:--|:--|
| base | `ic-main-995ea207` | `9c332039…` | v3, `995ea207` |
| candidate | `ic-r02b-on-995ea207-8b146d1a` | `b9cc60aa…` | `8b146d1a`: `995ea207` plus R02's `candidate.patch` |

- **No port is disclosed, because none was needed.** The patch applied
  unchanged: `8b146d1a`'s diff against `995ea207` adds and removes the
  patch's 140 lines, line for line.
- **The tests are R02b's:** `cargo test --release --lib -- koblitz_fast`
  on `8b146d1a`'s tree, 30 passed ([`tests.log`](tests.log), re-run for
  this record after the timed runs).
- **The full hashes** are in [`binaries.sha256`](binaries.sha256).

## What ran

- **The rows:** `M1`'s two rows, `M1-T01` and `M1-T02`, at R02b's two
  target sizes, `icv1-f2m59-tm943548413-98844ecc` and
  `icv1-f2m61-t158598901-ab42b6c5`.
- **A control size:** the same rows at `icv1-f2m53-tm56619371-dac20a85`.
  The narrow 8-lane kernel runs there in both arms, and nothing the patch
  changes runs.
- **Three rounds,** with the order alternating (base first, candidate
  first, base first).
- **Every process** was `ic price --single-target` with the suite's rho
  seed, `RAYON_NUM_THREADS=1`, isolated on CPU 2 by `isolated_bench run`:
  36 processes ([`explore.sh`](explore.sh)).
- **Two of the 36 were marked contended:** the base's `M1-T02` in
  round 3 at `2^47.2` and the candidate's `M1-T01` in round 1 at
  `2^44.5`. The script did not retry them.
  - The decision is read on the clean pairs, since contended runs are
    not pooled with clean ones.
  - The figure on all six pairs is shown beside it, and the verdict is
    the same.

## Results

Cold time (the median repetition's set-up plus online interval), base
over candidate, so a figure above 1 means the candidate is faster. The
intervals are 95% t intervals of the geometric mean.

| curve | `log₂ r` | clean pairs | cold time, base over candidate | 95% interval | all six pairs |
|:--|--:|--:|--:|:--|--:|
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 5 | **0.962** | [0.865, 1.069] | 0.960 |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 5 | **0.918** | [0.858, 0.982] | 0.922 |
| `icv1-f2m53-tm56619371-dac20a85` (control) | 44.3 | 6 | 1.002 | [0.917, 1.095] | 1.002 |

- **Collection is where the kernel acts, and it is slower.** Over the six
  processes of each arm, the median collection phase was 468 ms (base)
  against 525 ms (candidate) at `2^44.5`, and 2,604 against 2,821 ms at
  `2^47.2`. At the control size it was 465 against 459 ms.
- **Every pair recovered the same logarithm** in both arms, on all 18
  pairs ([`pairs.tsv`](pairs.tsv)).
- **The raw runs** are in [`runs.tar.xz`](runs.tar.xz): each process's
  report, isolation record and logs, checked by
  [`SHA256SUMS`](SHA256SUMS). The archive keeps the suite's own directory
  names.

## The decision

The rule (amendment 3, item 3): if the cold-time ratio's geometric mean
is below 1.05 at both target sizes, R02b is withdrawn without running.
- **It is below 1.05 at both,** at 0.962 and 0.918. Indeed it is below
  1: the candidate's cold time is about 4% and 8% longer.
- **So R02b is withdrawn,** and the kernel is retired, as a rejection
  would have retired it.
- **The withdrawal is neither an acceptance nor a failed round** on the
  scan (plan §11).

## Why

- **Main's #1242 removed the premium the kernel was built for.** It
  reduces every Gf2 product with two carry-less folds by the modulus's
  sparse tail.
- **The scan-stage diagnostic on v3** puts the scalar subtraction at the
  two wide-tail sizes at 11.3 and 12.1 ns a scanned summand
  ([record](../scan-stages-v3-20261006/README.md)). On v0′, R04
  measured 32.5 and 33.3 ns.
- **The kernel now costs more than the path it replaces.** R02's kernel
  carries `H·t`'s bits above `z^63` into a second fold, lane by lane.
  On v3 that is slower than main's scalar path: collection is 12% and 8%
  slower.

**What is left in the subtraction:** it is still about a third of the
scan on v3. A kernel built on main's two-fold reduction, vectorised, is
the lever there, and it is explored separately.
