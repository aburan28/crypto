# Frozen n37 native online wall screen

The [preceding target-only instruction panel](../ecbench_n37_online_ir_20261004/RESULT.md)
prioritized the K16 folded base for an isolated one-target online wall
comparison, while retaining K8 as the cold-setup control. Callgrind Ir is not
native wall time. This panel measures the same two complete IC pipelines and
the strong signed-Frobenius rho reference on **new public points** using the
native `ecbench` timer. It is a host-screen and reproduction step toward the
L2 comparison, not a claim that a hosted VM supplies L2 isolation.

`SPEC.json` freezes 16 one-target n37 Koblitz `a=0` workloads, public-point
law and seed `202610046101`, four arms, five measured rounds plus one warmup,
interleaved order seed `202610046102`, identical K16 A/A control, and a
120-second per-child cap. The methods retain the exact K8/K16 factor-base,
rank and target-descent policies of the preceding panel. Target index zero
is designated the **primary one-target row** before any run; indices 1–15
are secondary same-target replications for uncertainty and host-noise
assessment. The committed `PLAN.json` must show zero target-point overlap
with both earlier n37 K8/K16 panels before any execution. Do not replace a
failed point or change a cap after seeing outcomes.

The first native panel uses Ubuntu 24.04 on a GitHub-hosted runner, a single
release binary, `--cpus auto --allow-busy --settle 0.25`, and the harness's
fixed per-run core reservation where available. Preserve the host preflight,
CPU topology, affinity and isolation level, noise counters, executable and
source hashes, exact spec/plan/records, every child result, failures and
timeouts. Audit all measured records locally and replay every measured child
on a different host class. Keep the raw hosted session as an immutable
artifact, then commit compact raw records and independent replay receipts.

The **primary metric** is the target-zero verified online interval:
`rho_online_ms / IC_online_ms` for K8 and K16, paired on the same point,
algorithm seed, round and resource envelope. The IC timer starts after
target-independent factor-base, pair-table and rank-log preparation and
contains target query, PDP (including failed attempts), relation check,
descent, and scalar recovery check. Rho starts at its target-dependent walk
and stops after checked scalar recovery. Fixture point generation and process
startup are outside both timers. Check that the five exclusive IC phase
durations sum to its online duration, and retain per-run scalar replay and
the complete cold solve duration separately. The other 15 one-target rows
are reported individually; an across-target estimate is explicitly a
secondary panel diagnostic, never a substitute for the primary row.

Before execution, fix the secondary panel statistic to a ratio of sums of
paired online nanoseconds, with a 20,000-resample target-block percentile
interval at seed `202610046103`. Report K8/K16 and each IC/rho ratio, their
paired failures and timeout counts, the K16 A/A ratio, and the distribution
of actual earned isolation levels. A verified scalar and full independent
replay are required for any descriptive ratio; otherwise the affected cell
and overall ratio remain unknown. The phase table must retain base size,
folded columns, build and pair-table costs, rank trials and useful rows,
failed target attempts, online phase costs, cold solve, and memory peak.
The counted-GAE cold totals remain lower bounds while native costs are
unpriced and must not be divided into a wall-speedup claim.

The L2 gate from `AGENTS.md` remains strict: a speedup can be promoted only
if the *same frozen comparison* runs on an auditable physically isolated host
with all pairs verified, all required isolation/noise fields passing, and an
independent replay. A GitHub-hosted VM's L0/L1 results are exploratory even
if its native ratio is large. For resource prioritization only, a secondary
K8/K16 online ratio interval wholly above 1.10 and K16 A/A maximum relative
deviation below 5% keeps K16 first for an L2 run; an interval wholly below
0.90 selects K8; otherwise carry both. This does not certify IC/rho speed.
Any n41/n53 or ECC2K-130 transfer still requires separate matched full-rank
workloads and charged cold costs.
