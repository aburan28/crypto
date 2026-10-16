# N83 100-point factor-base selection panel

This secondary panel compares twelve complete N83 factor-base and root-index
configurations on the same 100 archived public points. The targets were
generated and hash-frozen before the four-byte capacity result was read.
The source is the exact four-byte producer with SHA-256
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.
The 100-target file has SHA-256
`fa0044f7cad85185812ac4c79f4feabb300ddcdc1732b6a0890778d5b8e72af6`;
its independent answer sidecar has SHA-256
`fcb167978cda4a72732c322442a3b77d3f0a997b74676a2a4d5767eb6c5262bf`.
The fixture excludes every earlier capacity, five-point selection, and fresh
Tier 1 public point. The producer receives only the public targets.

## Frozen arms and resource envelope

Every arm uses the same exact curve, rank seed 29, symmetric S3 traversal,
four rank threads, one target thread, complete factor-log setup, and
14-GiB observed `wait4` peak guard. Each arm starts in a fresh process and
processes all 100 targets after its own reusable base/index/log setup.
Order and caps are fixed below. Root4 uses four-byte exact slots and
`2^31` slots; root5/root11/root21 use five-byte exact slots and `2^30`
slots; root22 uses the five-byte load-75 request rule.

| Order | Arm | Usable points | Child cap |
| ---: | --- | ---: | ---: |
| 1 | K3240 root21 | 537,840 | 1300 s |
| 2 | K3400 root4 | 564,400 | 1500 s |
| 3 | K3300 root11 | 547,800 | 1200 s |
| 4 | K3500 root4 | 581,000 | 1500 s |
| 5 | K3220 root21 | 534,520 | 1300 s |
| 6 | K3550 root4 | 589,300 | 1500 s |
| 7 | K3114 root22 | 516,924 | 1100 s |
| 8 | K3596 root4 | 596,936 | 1500 s |
| 9 | K3200 root11 | 531,200 | 1100 s |
| 10 | K3400 root5 | 564,400 | 1300 s |
| 11 | K3190 root21 | 529,540 | 1300 s |
| 12 | K2940 root22 | 488,040 | 1000 s |

The K3400 root4/root5 controls must have identical complete base hash,
exact root-key count, scalar, and target probe counts on every point.
The independent verifier replays each complete base, subgroup membership,
four-point relation, recovered scalar, five exclusive online phases,
state count, and table allocation. A failed arm retains its raw row and
cannot be selected.

## Selection and validation

Among arms with all 100 target replays and observed peak at or below
14 GiB, select the lowest arithmetic mean target-online interval. Report
per-target times and probes, mean, median, maximum, target-block bootstrap
intervals, setup, memory peaks, and the paired mean difference from the
runner-up. This panel estimates target-search cost across the frozen input
law; each target's one-point online interval remains the unit of comparison.

Five reserve public points were also frozen before this panel, with target
SHA-256 `644c56a7bae2c92085704bb6f58ab70b38b278710de90152ece7631058c41dc9`
and answer sidecar SHA-256
`b30206982765ba649dcdf6f1eb48ded2cbe0fd3c8ec6a226a1fdfb261e5873b4`.
If this panel changes the selected candidate, those separate points supply
a new alternating same-point Tier 1 IC/rho validation workload.
