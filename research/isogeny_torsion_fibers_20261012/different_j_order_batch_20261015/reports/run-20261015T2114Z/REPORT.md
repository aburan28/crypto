# Order-229 different-\(j\) distinct-target run

**Run:** V1, 2026-10-15 21:14:48–21:14:52 UTC. **Host:** Darwin
25.6.0 arm64 T6041, shared L0. **Execution:** 84/84 children exited
zero, with 120-second wall ceilings and 512 MiB peak-RSS acceptance
ceilings. The separate Rust checker passed on all 84 cells and all
1,785 scalar rows, including 13,141 relation witnesses. These are
same-host correctness checks; the whole-child CPU values below are
descriptive shared-host measurements.

## Frozen question and instance

The [V1 protocol](../../PROTOCOL_V1.md) and [V2 prefix
correction](../../PROTOCOL_V2_PREFIX_CORRECTION.md) preceded this run.
The [run spec](../../RUN_SPEC_V1.json), SHA-256
`53a64bf1f1b4a55a07e6472ab9783158e43a2997dfc16e6df4cd03eec9eb5570`,
pins native source and binary hashes, the exact commands, public and
checker inputs, limits and seven arms. The runner and input were committed
at `d1253790ef3e3c98de170612d526fea21d161be9` before measurement.

The field is `F_211[s]/(s²−2)`. The degree-three map sends
`E:y²=x³−21x+37` with `j=77` to `E':y²=x³−81x−243` with
`j=160`. Its rational kernel is `{O,(3,1),(3,−1)}`. The source and
image generators both have order `r=229`. The earlier
[full fiber audit](../../../different_j_order_screen_20261015/FIBER_REPORT.md)
checked the 44,655 source points, 14,885 image points, every three-point
fiber and every one of the 52,441 subgroup-pair map identities.

The [solver-public input](../../PUBLIC_V1.json), SHA-256
`f6c3126b4d95e3c457c9e552471bc2cdff322de0e0a9de9fbba701d7c05897ba`,
contains 64 distinct source and image targets in a frozen coordinate-hash
order, with no scalar answers. The separate [checker-only input](../../CHECKER_V1.json),
SHA-256 `a3faedc54e9fa6d033e9fde9d972c4a1f233507f97edd1f6322a637963aac236`,
contains their scalar correspondences. Each arm received the same public
prefixes `N=1,4,16,64` in three rotated rounds. The direct controls were
reusable source BSGS, source negation rho and a repeated rho A/A arm.

The measured command was `zsh
research/isogeny_torsion_fibers_20261012/different_j_order_batch_20261015/run_v1.zsh`.
The script checked all source, binary and input hashes before a
three-round, 100,000-addition conversion calibration and the 84 fresh
metered child processes. Every child includes input/model construction,
setup, target queries, failed attempts, linear recovery, scalar
verification, startup and output. The [raw solver records](../../RUN_V1.jsonl)
have SHA-256 `d0e4616236aeb8978f5ac1549dde373be294bdf52503c954a6f964904f68faf8`;
the [whole-child meters](../../METERS_V1.jsonl) have SHA-256
`8cbb7e5d5737792231fa703d1bdf35535230132bff24062746abced60daaf776`.
The versioned script retained the 84 per-child stdout, stderr, exit and
meter files under [`raw-v1/`](../../raw-v1/), along with the wrapper,
calibration and checker exits.

## Observations

The single table uses **whole-child user+system CPU, microseconds**.
Each entry is the median of three rounds on the same public prefix.
The native [summary](../../SUMMARY_V1.json), SHA-256
`ad4b969aec7e6cd5aebf7f7d2ca0dfd851d22a96acfaef492838fc194d4e8af5`,
also records every round value, minima, maxima, field counters, peak RSS,
relation counts and paired ratios. The [graph](../../GRAPH_V1.svg) shows
the full-image arm and two direct controls with min–max whiskers; its
[150 dpi PNG](../../GRAPH_V1.png) was rendered from that SVG.

| Arm | N=1 | N=4 | N=16 | N=64 |
| --- | ---: | ---: | ---: | ---: |
| Source projected relations | 4,771 | 4,759 | 9,535 | 20,774 |
| Isogeny-image relations | 5,999 | 5,778 | 8,575 | 20,042 |
| Destination-direct relations | 5,039 | 6,229 | 11,332 | 22,054 |
| Full rational-image relations | 56,482 | 64,045 | 53,560 | 72,774 |
| Reusable source BSGS | 4,529 | 5,330 | 6,328 | 6,383 |
| Source negation rho | 5,030 | 6,333 | 4,797 | 7,653 |
| Rho A/A repeat | 5,248 | 4,980 | 6,088 | 7,357 |

At `N=64`, the full-image child CPU values were 70,847, 72,774 and
74,986 µs; BSGS values were 6,360, 6,383 and 6,794 µs; rho values
were 7,368, 7,653 and 8,686 µs. Full-image peak child RSS was at most
14,794,752 bytes, below the common 536,870,912-byte acceptance ceiling.
The full-image setup alone used 3,667,277 base-field multiplications and
89,636 inversions in each `N=64` child. Median complete `N=64` counters
were 5,138,577 multiplications and 143,637 inversions for full-image,
43,313 and 1,641 for BSGS, and 160,515 and 6,345 for rho.

The destination-direct two-sum table covered 90 point keys; the
full-image table covered 98. They shared 40 keys, with 50 direct-only
and 58 full-only keys. At `N=64` over three rounds, destination-direct
recorded 2,628 accepted relations in 6,776 trials; full-image recorded
2,329 in 5,389 trials. The source and transported-image arms each
recorded 2,433 relations in 5,763 trials. The separate
[checker receipt](../../REPLAY_V1.json), SHA-256
`c4572cc1dda31ddbfaec37190433eac60e8862769f77333de0d6ea6940d92037`,
reconstructed all target scalars, base witnesses, relation rows and
rank-nine solves independently.

The median calibration was 197.50791 CPU ns per repeated source-group
addition. The summary's CPU-derived local `S` is a conversion of these
local measurements, **not** a portable counted-operation bound. No L2
isolation, independent-host timing, polynomial/asymptotic result or
growing-order conclusion is claimed by this run report. See the
[results assessment](../../RESULTS.md) for the cost interpretation and
remaining gates.
