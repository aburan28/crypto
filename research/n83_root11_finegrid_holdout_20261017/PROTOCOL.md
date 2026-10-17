# Prospective N83 root11 fine-grid confirmation on 24 public points

The twenty-four public points in `targets.jsonl` were sealed before the local
fine-grid batch produced any outcome. `answers.jsonl` is a verification-only
sidecar. The curve is `EC1N83Ckb1a1haacb69ee9f37` over GF(2^83), with
polynomial integer `9671406556917033397649557`, equation
`y² + xy = x³ + x² + 1`, generator and subgroup order fixed by the canonical
curve record. The public-target SHA-256 is
`aa4453240d8c79d66a9cc3525d2bc9237140f04b4b5bb39b5dc5259356693cb3`;
the answer-sidecar SHA-256 is
`cf4ffaeef7efb2a7abf3334323f8e1ccb67c8ad3d9173cb86b0925ed52f73a1f`.
All 24 points are distinct and disjoint from both 128-point local batches,
the 20-point intermediate holdout, the 12-point final online holdout and the
10-point five-arm panel. The executed source SHA-256 must be
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.

## Arm rule fixed before the batch outcome

The fixed arms are K2543/root11, K2816/root11 and K3072/root11, using their
exact current-source IC1 candidate manifests. The fourth arm is the K value
with the lowest geometric mean of 128 target-online intervals in the separate
local fine-grid batch among K2752, K2816, K2880, K2944 and K3072. Break an
exact tie by the smaller K. Require complete independent scalar, base,
four-point-relation and source replay for every point and arm. If that gate is
incomplete, leave the fourth-arm selection null and do not dispatch this
holdout. Deduplicate by exact candidate ID, leaving three or four arms. The
local batch is a secondary shared-setup arm-selection diagnostic. Its
[batch plan](local_batch_plan.json), [verifier plan](local_verify_plan.json)
and [selection plan](local_selection_plan.json) were frozen before the batch.

The three fixed current-source candidate IDs are:

| K | Actual usable points | Folded columns | Exact candidate ID |
| ---: | ---: | ---: | --- |
| 2543 | 422,138 | 2,543 | `IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0h7a345a837c55` |
| 2816 | 467,456 | 2,816 | `IC1N83Ckb1a1fb467456PDP4rootRCguidedLAgaussTDdirectISO0hf693b1b33c8b` |
| 3072 | 509,952 | 3,072 | `IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a` |

The exact manifest for a selected new K must be materialized from its
independently replayed single-point base dump before dispatch. The unresolved
fourth arm has `candidate_id: null` at this fixture stage.

## Four Tier 1 host strata

Assign point indices 0–5, 6–11, 12–17 and 18–23 to four distinct cloud
instance or Pod allocations. Run every chosen arm on its six assigned public
points in fresh IC processes with empty base, root-index and log caches. Pair
each arm-point run with a separately seeded strong rho solve of the same
public point. Rotate arm order by point index and alternate `IC,rho` and
`rho,IC` order by point and arm. Use a frozen exposed-point rho calibration
before and after each pair. Record source and binary hashes, commands,
topology, affinity, cgroup limits, throttling, steal time, memory peak,
process and online intervals, all attempts, and complete raw outputs.
Admit a stratum only when all pair, calibration, scalar, base, relation,
resource, source, execution-order and provider-identity gates pass. Tier 1
does not require physical host isolation; report each cloud allocation's
observed envelope. An IC process has a 900-second cap and each rho or
calibration process a 120-second cap; use a 300-minute producer cap and
330-minute job cap per host. Preserve every timeout and noisy or incomplete
pair as a censored row.

## Frozen comparison

For each point and candidate, report the independently replayed one-target
online IC and same-point rho intervals, their ratio, exclusive target phases,
complete cold process, rank, base count, folded columns, probes and memory.
The online IC clock starts at target-dependent work after reusable preparation
and ends after scalar recovery and the producer's group check. The rho clock
uses the same public point and stops after its independent group check.
Report the external Python scalar replay separately. Base/index/log
preparation and launch belong to the supplementary complete cold interval.

For each unordered pair of distinct candidates, compare 24 paired log ratios
of online IC intervals and separately 24 paired cold-process differences.
Use the exact two-sided sign-flip test over `2^24` assignments and the
Bonferroni threshold `0.05 / C(m,2)` for m distinct arms. Select an online
winner only if its observed mean log interval is lower than every other arm
and every corresponding pair passes that threshold; apply the same rule
separately to cold process costs. Otherwise retain a null selection and all
intervals and p values. A verified Tier 1 one-target speedup against rho
requires independently solved matching public points and IC online time
below rho on all six points of that host.
