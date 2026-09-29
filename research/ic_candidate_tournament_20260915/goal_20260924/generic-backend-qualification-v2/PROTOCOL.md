# Generic F4/F5/SAT qualification, second registration

This protocol is frozen before any measurement under seed **2026092902**. The
first registration, seed 2026092901, ran for six hours and lost its raw bundle
when direct upload traversed the result tree. Its solver outcome is unknown.
The 25 public targets that registration could have generated are reconstructed
in [lost-campaign-exposures.json](lost-campaign-exposures.json), SHA-256
`a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41`.
The reproducer checks the sealed round archives and five retained fixtures
against the pinned worker. Treat all 25 points as exposed, whether or not a
timed-out process actually reached them. The first registration is never
redispatched or counted as a failed F4/SAT solve.

## Question and success rule

Can at least one source-bound generic F4/F5 arm and one generic SAT arm complete
independently verified, one-target IC solves on **every** smoke and development
point in the five-cell panel, with natural-query evidence and exclusive costs?
For each admitted arm, require all 5 smoke and all 15 development jobs verified,
with no censored query audit, checked usable base and folded columns, observed
PDP dispatch, stored relation matrix and rank, final relation LA, target
descent, scalar replay, complete profile/native phase accounting and matching
source/build identity. The [family gate](../../generic_backend_gate_v2.py) reads
the complete tournament summary and independently audited natural-query rows.
One scalar produced by a different path, a planted decomposition, or a
successful stage-only microbenchmark does not qualify an arm.

The accepted incumbent and the two rho roles remain paired on the same public
point and resource envelope. A candidate that fails, times out, OOMs, or has
unverified work retains a row and null competitive totals. The whole schedule
must finish for a family verdict; a runner timeout is an operationally censored
campaign, not a negative solver result. Do not widen limits, substitute a new
point, or retry seed 2026092902 after dispatch. A later experiment needs its
own seed and exclusions.

This is a development qualification, with no held-out confirmation and no
familywise promotion or ECC2K-130 extrapolation. Preserve each complete
variant, per-cell specialist, Pareto tradeoff and exploration candidate rather
than choosing only the local winner. The accepted three improvement rounds and
their confirmation/replay sets remain sealed.

## Frozen panel and source

[panel.json](panel.json) has SHA-256
`d283a869b0412228d1c66260fdfd8f387d7243bd15456c7febf3c46ee5da27a8`.
Its algorithm arms, exact configs, references, five curve cells, and timeout
are unchanged from the first registration. Seed 2026092902 generates fresh
public-hash targets only after the original history, three sealed improvement
archives, seven supplemental fixture corpora, and the lost first-run exposures
are checked and excluded by exact curve ID and point. The n29a1 holdout is not
generated or inspected.

| Stage | Distinct public points | Paired arms | Processes per point | Trial slots |
| --- | ---: | ---: | ---: | ---: |
| A/A | 5 | incumbent, identical control | 1 | 10 |
| Smoke | 5 | 10 IC plus 2 rho | 1 | 60 |
| Development | 15 | 10 IC plus 2 rho | 1 | 180 |
| Total | 25 | — | — | 250 |

The former three process repetitions on each point are replaced by one. They
were neither independent targets nor additional natural-yield samples. All 25
distinct points and all five curve cells remain. No repeated-process precision
claim follows from this schedule; if a qualified arm merits confirmation, run
fresh point and process controls under a separate registration. Cap the runner
at 300 trial slots. Each child has 300 seconds, 8 GiB, one pinned CPU and one
Rayon thread; Linux amd64, Rust 1.94.1 and Valgrind 3.22.0 are required.

The generic worker is built from commit
`765c3c5f19032bd852163805f257c56babef2040` with the same pinned `src`,
example, `Cargo.toml`, calibration input and CI `Cargo.lock` objects as the
first registration. `KIC_F5_AVX512_UNPACK=0`. The new evaluator, schedule and
packaging source are retained in the result bundle. Prepared pairinv/both
reference sources come from the accepted reference archive. `generic_build`
records the actual compiler, dependency source trees, flags, executable hash
and build policy; a source mismatch aborts before fresh fixture generation.

## Cost and uncertainty

Primary metric: verified one-target native online wall time, from first
target-dependent work after reusable preparation through independent scalar
replay. The paired rho solve uses the **same supplied public point** and its
own target-dependent walk-to-replay interval. Report rho/IC online wall ratio
only for verified complete pairs. Supplementary metrics: full cold process
time and complete user-space Valgrind guest instructions, including setup,
failed attempts, factor-base construction, query generation, PDP, checking,
matrix build, final relation LA, descent, replay and report overhead. Report
`S = total Ir / sqrt(r)`, ratios to matched rho and the declared collector
floor with their boundaries fixed. Missing phase costs stay unknown.

Audit every ordinary sampled query and failed attempt. Keep proved UNSAT,
unresolved, timeout and absent report distinct; no missing report is zero
yield. Report actual subgroup-usable factor-base points before sign/Frobenius
folding, effective columns, accepted rows, novel rank, target descent,
verified witnesses and final rank. Bootstrap natural rates by **distinct
point** within fixed cells; report the conservative point-level interval too.
One process per point gives no within-point timing variance estimate, so
paired cost intervals are descriptive. Preserve all five cells and zero-yield
cells; do not select a favorable subset.

## Run and evidence retention

The measured command has a 300-minute **step** cap inside a 360-minute job,
reserving at least 60 minutes for an `if: always()` packaging and upload path.
The packer turns complete or partial raw output into one `.tar.zst` file with
SHA-256, byte count and capture status. A separate PR check uploads a synthetic
partial result through the same path before campaign dispatch. The GitHub run
attempt must be one. If the measured step times out, preserve its `state.json`,
receipts and logs, mark the panel incomplete and register a new seed for any
later measurement. Never infer a solve from upload success alone.

Only after `tournament.py verify`, the independent natural-query audit and the
family gate succeed may an evidence PR report qualification and paired costs.
Publish exact artifact link, archive hash, count of scheduled/verified/censored
jobs, source/build hashes, target census and the complete comparison table.
The repo scoreboard remains pending until that evidence is reviewed. This
registration and its workflow are a separate PR before any dispatch.
