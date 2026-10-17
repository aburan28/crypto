# Sixty distinct targets: full-image setup reuse and whole-process price

**Runs:** V3 began 2026-10-15 19:49:52 UTC; V4 began 19:54:56 UTC.
**Host:** Darwin arm64, macOS 26.6, shared L0. The
[V1 protocol](../../PROTOCOL_V1.md) and its
[pre-run target-order correction](../../TARGET_ORDER_CORRECTION_V1B.md)
preceded the native measurement. [V3](../../RUN_SPEC_V3.json) and
[V4](../../RUN_SPEC_V4.json) each pinned the solver source, binary,
public/checker input and script hashes before their run. V4 reused
the identical solver binary and 60 point pairs, adding a separately
hash-pinned [Rust process meter](../../../../../examples/different_j_process_meter.rs).

## Input and correctness

The explicit degree-three map connects `j=7` to `j=11` curves over
`F_(67²)` with a source subgroup of order 61. The target file contains
all **60 distinct nonidentity** source subgroup points and their mapped
destination points, sorted by a hash of public point coordinates.
The solver child received only the [public file](../../PUBLIC_BATCH_V2.json).
The separate [checker file](../../CHECKER_BATCH_V2.json) holds known
scalar correspondences for replay. A deliberately answer-bearing input
was [rejected with exit 101](../../CUSTODY_REJECT_V2.stderr) before
timing. The freeze process and shell hash gate accessed the checker
file; the measured solver process did not.

Each of the seven arms solved all 60 targets in each of three rotated
rounds. The [V3 replay](../../REPLAY_V3.jsonl) and
[V4 replay](../../REPLAY_V4.jsonl) each independently checked
**1,260/1,260** recovered scalars, all 60 distinct source/destination
point pairs, eight full-image preimage witnesses, and **10,960**
relation witnesses with separate field and point-add formulas. The
four relation arms reached rank nine on every target; no arm reached
its 1,220-trial per-target cap or 512 MiB ceiling. The native
[V4 receipt checker](../../SUMMARY_V4.json) matched all 21 process
meter records to their raw stdout, stderr, exit codes and replay.

The earlier [V2 script](../../RUN_V2_FAILURE.md) stopped after one
source arm because zsh reserves `status` as a read-only parameter.
Its raw stdout and wrapper stderr remain versioned. It supplied no
multi-arm comparison. The V3 wrapper changed that shell variable and
used new output paths; V4 then added the whole-process meter.

## One-unit cost table

All table entries are **microseconds of CPU**, medians of three rounds.
V3 columns are cumulative current-thread CPU through each prefix;
they charge JSON input, curve construction, validation of all 60
public points, one arm setup, every query, linear recovery and final
scalar check. The input validation is charged once even for the
shorter prefixes because the process reads the full frozen batch.
The final column is V4 child **user + system CPU** for the complete
60-target process, including startup, final JSON serialization and
stdout. V3 and V4 were separate shared-host runs; values across them
are descriptive, not paired differences.

| Arm | V3 N=1 CPU (µs) | V3 N=4 CPU (µs) | V3 N=16 CPU (µs) | V3 N=60 CPU (µs) | V4 whole process N=60 CPU (µs) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Source projected | 661.666 | 829.875 | 1,704.292 | 5,872.041 | 15,199 |
| Isogeny-image projected | 727.875 | 1,108.708 | 2,711.000 | 6,900.125 | 12,735 |
| Destination direct | 613.417 | 907.042 | 1,732.584 | 5,389.917 | 14,021 |
| Full rational image | 4,538.750 | 4,752.500 | 6,528.750 | 11,049.667 | 15,351 |
| Direct source BSGS | 407.833 | 424.500 | 468.833 | 643.208 | 6,148 |
| Direct negation rho | 624.125 | 729.459 | 999.167 | 2,431.542 | 6,908 |
| Rho A/A stream | 832.792 | 961.625 | 1,377.625 | 2,392.959 | 8,028 |

The **paired** V3 full-image/BSGS ratio medians at `N=1,4,16,60`
were **11.039, 11.104, 13.926, 18.103**; paired full-image/rho
medians were **7.213, 6.462, 6.534, 5.284**. The paired V4
whole-process full-image/BSGS ratios were **2.196, 2.483, 3.473**
(median **2.483**), and full-image/rho ratios were **2.222, 2.166,
2.906** (median **2.222**). These paired medians are computed from
same-round totals, so they need not equal ratios of table medians.

The V3 full-image `N=60` median common validation, image/table setup,
per-target work and timing/accounting residual were **438.875,
4,095.917, 5,873.456 and 58.669 µs**, respectively. Phase medians
come from different rounds and should not be summed to recreate the
median total. Full-image setup reuse reduced its displayed CPU per
target from 4,538.750 µs at `N=1` to 184.161 µs at `N=60`;
direct BSGS was 10.720 µs per target at `N=60`. The 60-target
full-image relation trials were 1,145, 987 and 976 across rounds,
versus 1,022, 972 and 915 for source projected. Source and
isogeny-image arms had identical trial and relation counts in every
round, consistent with their identical frozen support sets.

In the V4 run, median child-minus-native-timer CPU was 6,330.083 µs
for full-image, 5,120.250 µs for BSGS and 4,788.666 µs for rho.
That difference includes process startup, output and timer-scope
differences; it is not an isolated startup estimate. V4 median full
process wall times were 26.344, 10.822 and 10.621 ms for those arms.
Peak child RSS for full-image was at most 9,764,864 bytes across its
three V4 rounds. The [native V3 summary](../../SUMMARY_V3.json) and
[metered V4 summary](../../SUMMARY_V4.json) retain every round,
prefix, phase, primitive count, wall time and memory receipt.

![V3 cold prefix CPU from the native timer](../../GRAPH_V3.svg)

![V4 whole-child CPU and native timer in the same process](../../GRAPH_V4.svg)

## Decision and limits

This bounded setup-reuse experiment changed the full-image cost
structure as intended: one 1,525-point rational-image construction
served 60 genuinely distinct target solves. On this frozen order-61
group, both complete local CPU views remained above direct source
BSGS and rho, and direct BSGS stayed lower at every tested prefix.
The full-image base's different support membership did not produce a
lower complete solve cost here. This observation does not decide a
larger-order curve, another factor-base selection or a direct
χ-coordinate resultant.

The V3 prefix CPU excludes process startup and final output. V4
charges those for `N=60` only, so no whole-process claim is made for
`N=1,4,16`. These are shared-host L0 CPU and wall measurements with
three rounds, not isolated performance intervals or a portable
ecbench `S` in counted group-operation equivalents. Field/group
primitive ledgers and raw results remain available for a later
counted-unit conversion and second-host replay. No exponent or
deployed-curve claim follows from this fixed instance.

## Reproduction and hashes

- [V3 run script](../../run_v3.zsh), [raw 21-row run](../../RUN_V3.jsonl)
  SHA-256 `ac85539489c417981adcc41ab5f1cca00539f55269c8ec83ebe4d7231a683e10`,
  [replay](../../REPLAY_V3.jsonl) SHA-256
  `47d6ee15964c390968fa9636068e4e0b7a9cb4b519a548622c17053dbff85d1f`,
  and [native summary](../../SUMMARY_V3.json) SHA-256
  `f41d682b5c77441295d79d8be33e0f101e12daf5483e4a9c080792acc2bcf89c`.
- [V4 script](../../run_v4.zsh), [raw solver output](../../RUN_V4.jsonl)
  SHA-256 `7d5a055d348cb6fa8edb5fc85b648a4850236adef10292506533a3ef902ff421`,
  [process meters](../../METERS_V4.jsonl) SHA-256
  `72890c5526949f069eb6afbb0857ff20d63b0ec9fde680e13ae7ceb0adbd7a0e`,
  [replay](../../REPLAY_V4.jsonl) SHA-256
  `1b81efec41a26c6d6791fb3f32e9fd6f12ee2b106f7169da7a409ed8efdaa09a`,
  and [V4 summary](../../SUMMARY_V4.json) SHA-256
  `3599e7fb89fdfdd039a703ba3b260be19ce5fb0bac4ca0b92e815650294b9ffc`.
- The [V3 raw directory](../../raw-v3/) holds 21 stdout/stderr/exit
  triples. The [V4 raw directory](../../raw-v4/) holds each solver
  stdout/stderr, meter JSON/stdout/stderr and exit code. The
  [V2 failure](../../RUN_V2_FAILURE.md) and
  [custody rejection](../../CUSTODY_REJECT_V2.json) remain separate.
