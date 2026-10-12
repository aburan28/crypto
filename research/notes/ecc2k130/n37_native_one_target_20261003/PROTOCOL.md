# Frozen n37 one-target descendant-native IC versus signed-Frobenius rho

Status: **preregistered, before opening or solving any b02 target**. This is
the first primary one-target comparison following the b00 direct-support and
b01 bounded-residual batch diagnostics (merged PRs #1254 and #1261). Neither
of those batches supplies a one-target attack-speed estimate.

## Question, inputs, and stop rules

For each index `i=0..31` of the orbit-disjoint b02 **public point** file,
solve Q independently from a fresh process and a fresh target-dependent
state. The IC hypothesis is that the fixed 42-point descendant-native factor
base, target-blind full-rank relation log, exact signed 3+3 oracle, and at
most 16 frozen target-independent shifts recover every Q. A single missing
or incorrect answer, incomplete oracle decision, relation rank below 42,
or failed independent replay falsifies that correctness hypothesis. The
separate performance hypothesis is lower **cold total operations**, in one
calibrated unit and with correct answers on every point, than the matched
signed-Frobenius rho arm. Leave `S` and speedup unset until that accounting,
the complete one-target phases, and a matched resource envelope exist.

The frozen file is
`research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b02.points.jsonl`,
SHA-256 `248d405e834b2ee7c27752e5742ca568d21156d63e981fa49f362cdaab2f9cd9`.
It contains 1024 public points from corpus
`compact-disjoint-cold-v2-n37-L1024-b02-20261001`, seed `2026100110102`.
The producer may read **only** that point file. Independent replay may read
the separate b02 scalar fixture, SHA-256
`0d23bb74f85fe3654c55a502c868c9491ab2d2d28417e3d35739c233e128b11c`,
only after both arms have written their recovered answers. The common
`FROZEN.json` SHA-256 is
`da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d`.
Use the degree-73 archive SHA-256
`eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90`
and `NATIVE42.json` SHA-256
`bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c`.
The source curve is Koblitz `y²+xy=x³+1` over the frozen GF(2^37)
representation; subgroup order `r=230603167`.

## Fixed algorithms and accounting

IC rebuilds and verifies the degree-73 bridge, all 42 leaf base points, the
signed half table, and the target-blind relation matrix in **each** process.
The rank probe uses SplitMix64 state `0x6e33376d365f3031`, cap 256, and
stops at rank 42. Solve and full-point-check all base logs before online
timing. The ordered list of 16 distinct nonzero residual shifts uses state
`0x6e33375f72657331`; precompute their points before online timing. These
settings are fixed from the prior b01 protocol, with no b02 tuning.

After input loading, the IC online interval begins at Q validation and
degree-73 transport and stops after source and leaf full-point scalar
verification. Its consecutive, exclusive phases are `target_query`
(validation, subgroup checks, transport), `target_pdp` (direct and, only
after complete misses, ordered shifted oracle queries),
`target_relation_check` (external full-point witness replay),
`target_descent` (combine base logs, subtract shift), and
`recovery_check` (both full-point scalar checks). Failed attempts and their
operation counts remain in the raw row. The five phase durations must sum
exactly to the charged online duration. Reusable setup, input/output,
whole-process time, and all group-operation counters are reported separately.
No known-answer scalar may enter either producer.

The comparator is `koblitz_rho_batch_ks_strong_online` at rung 3, 32 lanes,
8 distinguished-point bits, command
`37 0 signed_frobenius 1 <seed>`, with
`KIC_RHO_TARGET_POINT=<public x>,<public y>` and
`KIC_RHO_BATCH_CORPUS=compact-disjoint-cold-v2-n37-L1024-b02-20261001`.
For index `i` and repetition `j=0..4`, fix the rho seed to
`2026100110102 + 1000*i + j`. Run only one point per process: no shared
distinguished-point table or cross-target collision. Rho online starts at
its first target-dependent walk and ends after its own `[d]G=Q` check.
Report exclusive `walk`, `collision`, and `recovery_check` durations and
all operation counts. Keep cold setup and process costs separate.

Record one row per `(candidate_id, workload_id, run_id)`; the workload is
the indexed b02 Q, and each repetition is a fresh execution of both arms.
Run an A/A control (at least five paired repetitions) before timing claims,
then interleave IC/rho order across five A/B repetitions per point. Fix the
host, binary, thread count (`RAYON_NUM_THREADS=1`), and CPU/resource envelope;
exclude contended runs from runtime inference while retaining them as raw
evidence. Report all 32 per-point outcomes and failures. A wall-time claim
additionally needs a paired 95% interval excluding no improvement and a
difference beyond the A/A spread. The algorithmic comparison needs a
calibrated common operation unit and all cold costs; otherwise wall times
are descriptive and `S`, speedup, and a crossover verdict stay null.

An independent source-curve replay must reconstruct the frozen base and
shift schedule, check every queried point and oracle support/miss decision
using its own reachable-set closure, solve the relation logs independently,
check every recovered scalar on the source curve, and then compare the b02
fixture. Preserve replay certificate hashes and all failed runs. The n37
panel is a feasibility gate only; n41, n53, m=83, and ECC2K-130 require
their own matched evidence. The generic counting floor and the accepted
rho reference remain boundaries; a stage improvement alone is not an
algorithmic advance.
