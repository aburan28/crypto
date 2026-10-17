# Prospective N83 root11 intermediate-size confirmation

The twenty public points in `targets.jsonl` and the verification-only
answer sidecar were sealed before the 128-point batch screen finished. Their
hashes, generation law, exact curve and subgroup are in `manifest.json`.
These points have no overlap with the ten-point five-arm panel, the twelve-point
final online holdout, or the 128-point batch screen.

## Question and arm rule

The primary question is which measured root11 factor base gives the shortest
verified one-target online interval for one previously unseen N83 public point
under the frozen resource envelope. The fixed arms are K2304, K2543 and K3072.
The fourth arm is the K value with the lowest geometric mean of the 128
independently verified target-online intervals in the secondary shared-setup
batch among K2304, K2543, K2688, K2816 and K3072. Break an exact tie by the
smaller K. If any of the five batch arms lacks complete independent replay of
all 128 answers, use K2543 as the fallback. Deduplicate by exact candidate ID,
leaving three or four arms. Batch timing is only an arm-selection input; it is
not a one-target cloud speedup result. The separate twelve-point holdout
outcome plays no role in this arm rule.

Every executed arm must use the exact current-source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and the root11 dynamic indexed-five layout. The base points were independently
replayed from the single-target dumps and materialized as exact `IC1` records
before the full batch outcome:

| K | Actual usable points | Exact candidate ID |
| ---: | ---: | --- |
| 2304 | 382,464 | `IC1N83Ckb1a1fb382464PDP4rootRCguidedLAgaussTDdirectISO0hf7c27ecf3b60` |
| 2543 | 422,138 | `IC1N83Ckb1a1fb422138PDP4rootRCguidedLAgaussTDdirectISO0h7a345a837c55` |
| 2688 | 446,208 | `IC1N83Ckb1a1fb446208PDP4rootRCguidedLAgaussTDdirectISO0hc92369996197` |
| 2816 | 467,456 | `IC1N83Ckb1a1fb467456PDP4rootRCguidedLAgaussTDdirectISO0hf693b1b33c8b` |
| 3072 | 509,952 | `IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a` |

The complete manifests retain the enumerated base digest and exact source,
curve and stage configuration; their bytes must be pinned before dispatch.

## Four Tier 1 strata

Run points 0–4, 5–9, 10–14 and 15–19 on four distinct cloud instance or Pod
allocations. Each host runs every chosen arm on its five assigned points, with
a fresh process and empty base, index and log caches for every arm-point run.
For each arm-point pair, solve and independently verify the same public point
with a separately seeded strong rho process. Rotate arm order by global point
index; alternate `IC,rho` and `rho,IC` by point and arm. Run the frozen
exposed-point rho calibration before and after each pair. Record source and
binary hashes, exact command and environment, topology, affinity, cgroup
limits, throttling, steal time, memory peak, process and online intervals,
all failed attempts, and the full raw outputs. A host stratum is admitted only
when every pair, calibration, base and relation, scalar, resource envelope,
and event-order check passes. Physical host isolation is not required for
Tier 1; preserve its actual cloud-instance envelope.

Set each IC process cap at 900 seconds and each rho/calibration cap at 120
seconds, with a 225-minute producer cap and 260-minute job cap per host. A
global-cap stop remains censored and prevents aggregate admission. The same
source and factor-base hashes must appear across all four strata for each
candidate. Distinct concurrent provider runner identities must match the
archived host receipts before pooling.

## Frozen analysis

The primary metric is `rho_online_ms / IC_online_ms` for each exact point and
arm, with online clocks starting at target-dependent work and stopping after
producer scalar recovery and its group check. Independent Python scalar replay
is outside both online intervals and its time is reported separately.
Base/index/log construction and process launch are outside the online
interval; report complete cold process time separately. Report every point
and every arm, per-host means and median
paired rho/IC ratios, phase costs, rank, probes, memory and correctness.

For each unordered pair of distinct chosen arms, compare the twenty paired
log ratios of target-online intervals and the twenty paired differences of
complete cold process costs. Use the exact two-sided sign-flip test over all
`2^20` assignments. Select an online winner only if its observed mean log
interval is lower than every other arm and every corresponding p value is at
most `0.05 / C(m,2)` for m distinct arms. Apply the same rule separately to
the cold metric. Otherwise retain a null selection and the full interval
and p-value table. A verified one-target Tier 1 speedup against rho requires
an independently solved same point and IC online time below rho on all five
points of that host. Any failed or noisy pair remains a row and blocks the
relevant aggregate rather than being silently dropped.
