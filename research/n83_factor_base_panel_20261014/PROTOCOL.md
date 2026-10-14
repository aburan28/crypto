# N83 secondary factor-base selection panel

This panel compares seven complete factor-base/index configurations on the
same five archived public points. It follows the single-target capacity
measurements for each included layout. The selection variable is factor-base
size and state layout; every arm uses the same N83 curve, rank seed 29,
symmetric S3 traversal, exact five-byte indexed root table, four rank
threads, one target thread, and 14-GiB observed peak RSS envelope.

## Frozen input and selection rule

The target file is
`research/n83_evaluation_20261012/optimization/fixtures/holdout_5_targets.jsonl`,
SHA-256 `b96672b4d75df4b90a337c1e471d54dc147edbf236ef8d0054adb716d2314871`.
The five scalar sidecars were independently replayed in the archived
`verification/final_receipt.json`; this panel's verifier independently
multiplies the generator by each scalar again. The producer receives only
the five public points. The frozen source SHA-256 is
`07ce5094f0e561a75bd65ca143f2301830447742268936805c0b920f7a93f6ae`.
The point file is a **secondary** selection workload, run after the primary
one-point gates. It cannot replace a fresh one-target IC/rho timing panel.

| Order | Arm | Usable points | State layout | Table rule | Child cap |
| ---: | --- | ---: | --- | --- | ---: |
| 1 | K3400 root11 | 564,400 | one 83-bit root, second reconstructed | fixed `2^30` slots | 1500 s |
| 2 | K3200 root11 | 531,200 | one 83-bit root, second reconstructed | fixed `2^30` slots | 1300 s |
| 3 | K3240 root21 | 537,840 | two retained roots | fixed `2^30` slots | 1300 s |
| 4 | K3114 root22 | 516,924 | two retained roots | `ceil(4S/3)` request | 900 s |
| 5 | K2940 root22 | 488,040 | two retained roots | `ceil(4S/3)` request | 800 s |
| 6 | K2543 root22 | 422,138 | two retained roots | `ceil(4S/3)` request | 700 s |
| 7 | K2543 root11 | 422,138 | one 83-bit root, second reconstructed | fixed `2^30` slots | 700 s |

Each arm builds its base, final relation rank, factor logs and root index
once in a fresh process, then queries the five points in file order.
Preparation is reported separately. The target-dependent interval and its
five exclusive phases are reported **per point**; no five-target average is
used as a primary speedup. The verifier independently replays every complete
base, subgroup, four-point relation, recovered scalar, phase sum, exact
state count and table allocation. The K2543 controls must have identical
base hash, exact root-key count and target probe counts.

Among arms with all five verified targets and a final `wait4` peak at or
below 14 GiB, the descriptive selection rule chooses the lowest arithmetic
mean of the five online intervals. Report all per-target times, median,
maximum, target-block bootstrap intervals, setup and resource peaks, and
retain timeouts, guard stops, errors
and failed checks. The selected arm is validated on a **new** five-point
Tier 1 alternating same-point IC/rho panel after this comparison; panel
selection alone is an exploratory measurement.
