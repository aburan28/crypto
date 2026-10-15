# N83 factor-base refinement on 1,000 public points

This secondary panel refines the 100-point selection around K2940 and tests
the largest four-byte root index. Its 1,000 distinct N83 subgroup points and
independent answer sidecar were frozen before any refinement arm ran. The
public-point file has SHA-256
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`;
the answer sidecar has SHA-256
`f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71`.
The producer receives only the public points. They are disjoint from the
earlier capacity, 100-point selection, and fresh Tier 1 workloads.

## Frozen arms and resource envelope

Every arm uses the same curve, rank seed 29, symmetric S3 traversal, four
rank threads, one target thread, and complete factor-log setup. A fresh
process per arm receives an observed 14-GiB `wait4` peak guard and the child
cap below. Root4 uses exact four-byte slots and `2^31` fixed slots; root11
and root21 use exact five-byte slots and `2^30` fixed slots; root22 uses the
five-byte load-75 request rule. The executable source is the merged root4
producer with SHA-256
`5e9fd826d22816b429a82e658c517d9387916aabe7660f708e4c29901cce647e`.

| Order | Arm | Usable points | Child cap |
| ---: | --- | ---: | ---: |
| 1 | K2940 root22 | 488,040 | 1200 s |
| 2 | K2700 root22 | 448,200 | 1200 s |
| 3 | K3220 root21 | 534,520 | 1200 s |
| 4 | K3060 root22 | 507,960 | 1200 s |
| 5 | K2800 root22 | 464,800 | 1200 s |
| 6 | K3300 root11 | 547,800 | 1200 s |
| 7 | K2600 root22 | 431,600 | 1200 s |
| 8 | K3000 root22 | 498,000 | 1200 s |
| 9 | K3596 root4 | 596,936 | 1200 s |
| 10 | K2900 root22 | 481,400 | 1200 s |
| 11 | K3100 root22 | 514,600 | 1200 s |
| 12 | K2850 root22 | 473,100 | 1200 s |
| 13 | K2980 root22 | 494,680 | 1200 s |
| 14 | K3020 root22 | 501,320 | 1200 s |

The independent verifier checks the complete base and its subgroup points,
four-point relations, recovered scalars, five exclusive online phases, state
counts, and table allocations. K2940 root22 must reproduce the previously
sealed base hash
`12c164dd4bbc72df19d8b0c638f1d808c534bd464397c9033822dace74cca92b`.
Failed, timed-out, over-cap, or unverified arms remain in the result and are
ineligible for selection.

## Selection and fresh validation

Among arms with all 1,000 point replays and observed peak at or below 14 GiB,
select the lowest arithmetic mean target-online interval. Report every
point's time and probe count, mean, median, maximum, target-block bootstrap
intervals, setup, memory peak, and paired mean difference from the runner-up.
This panel measures selection cost under one frozen input law. A candidate
change then uses a separate five-point Tier 1 same-point IC/rho workload.
Those five points were frozen before this panel, with public-point SHA-256
`2ccaac2baf367598e0d4e724b500ba66c5b4a711f0be8207565ce3acd6b1d0b4`
and answer-sidecar SHA-256
`5b170d9bc4d34fdd21c032990ceee3a62af335c9af7fcaae2596335b8699dd8c`.
