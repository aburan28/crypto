# N83 four-byte root slots with 21-byte states

This implementation pairs exact four-byte root slots with the direct
21-byte two-root state layout. A new `KIC_INDEXED4_FIXED30=1` mode requests
`2^29` entries, producing `2^30` four-byte slots (4 GiB). Its admission
checks require N83, symmetric S3 traversal, `KIC_ROOT_ONLY21=1`, a 29-bit
state ID, and fewer than `2^30` possible root insertions. Existing
four-byte/11-byte fixed31 and five-byte modes remain available. The small
exact-root parity test checks both present and absent root keys against the
five-byte table with the same first witness.

## Frozen comparison

Six arms use the same curve, rank seed 29, four rank threads, one target
thread, and one fresh process per arm under a 14-GiB observed RSS guard.
Each child has a 1,200-second cap. The IC source SHA-256 is
`68f9b22972ee09511e15ae74a998ba74b82e74eed36e3014bd0460c918b6cdcb`.
The common 1,000 public points are the previously frozen refinement fixture
from commit `b1a4131558b3121ed562ca830abce972aa1dc580`, SHA-256
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`;
their independent answer sidecar SHA-256 is
`f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71`.
This is a secondary paired layout comparison on known common inputs. The
producer receives only public points.

| Order | Arm | Usable points | Index slots | State bytes |
| ---: | --- | ---: | ---: | ---: |
| 1 | K2940 root4f30 | 488,040 | 2^30 × 4 | 21 |
| 2 | K2940 root22 | 488,040 | load-75 × 5 | 22 |
| 3 | K3200 root4f30 | 531,200 | 2^30 × 4 | 21 |
| 4 | K3200 root21 | 531,200 | 2^30 × 5 | 21 |
| 5 | K3300 root4f30 | 547,800 | 2^30 × 4 | 21 |
| 6 | K3300 root11 | 547,800 | 2^30 × 5 | 11 |

The independent verifier replays every full base, subgroup point, four-point
relation, recovered scalar and exclusive target-online phase sum. Within
each same-K pair, base hash, exact root-key count and all 1,000 target probe
counts must agree. Failed, over-cap, and parity-failed arms remain in the
receipt. Among arms with all target replays and peak RSS at or below 14 GiB,
select the lowest arithmetic mean of 1,000 target-online intervals. Report
per-target times and probes, setup, memory, target-block bootstrap intervals,
and paired difference to the runner-up. This comparison does not contain rho.

Five new public points were also frozen before any run of this mode. Their
target SHA-256 is
`fef03b92aef0352e2e21ff7362878b8bed8f4826fe8d90503d70174508d631f1`;
their answer-sidecar SHA-256 is
`a1979e2bbcd3d82d47e29f38d66056228e3733faf125b0d28f8e4c32831b5665`.
They are disjoint from the common 1,000, prior selection and validation
points, and the other refinement reserve set. A changed selected candidate
uses those five points for a separate alternating same-point Tier 1 IC/rho
validation.
