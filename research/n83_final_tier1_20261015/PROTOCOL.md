# Five same-point N83 Tier 1 IC/rho pairs after finalist selection

The five public subgroup points in `targets.jsonl` were fixed before the
final factor-base tournament ran. Their SHA-256 is
`34c530851eb135b7e348ac976caaae37c5883c52b52962ad496028c058b1dcda`.
The verifier-only scalar sidecar SHA-256 is
`8705b909479f4ee1eb8a19fe327c0622cd890e0d8443e08bdbdf04a4c6b99548`.
`generate.py` excludes every prior N83 panel and all 4,000 tournament points,
checks curve membership and subgroup scalar construction, and refuses a
changed fixture. The curve is `EC1N83Ckb1a1haacb69ee9f37` over polynomial
basis `GF(2^83)`, equation `y² + xy = x³ + x² + 1`, subgroup order
`8569786107849059`.

Execute this panel only if the final tournament's independently replayed
analysis selects one candidate under `PROTOCOL.md` there. Pin the selected
candidate ID, its source manifest and base hash, the tournament raw artifact
and analysis digests, and the complete `admission.json` digest before starting
a Tier 1 allocation. The IC producer source SHA-256 is
`f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc`;
the paired strong-rho source SHA-256 is
`fb14072f4fa5d08e18d8b79130b6d58a220adae0f76369cc9be949691c83eaa0`,
with rho configuration SHA-256
`cf76ebb5372acfc1355b88374f875980be66cdffba1a82cb6e4c2f28bcfa723b`.
Both use four visible CPUs within the same affinity and memory envelope.
The versioned `selection.json` uses schema `n83.final-tier1-selection/v1` and
retains that candidate's pipeline, ID, K, variant, Bloom profile, base hash,
root-key count, manifest digest, the tournament ZIP and analysis digests,
its artifact ID, and the admission digest. The runner checks these against the admitted
candidate before executing.

Use one cloud VM or Pod allocation for all five points and record its
provider, allocation ID, region, start/end, CPU/core/SMT/NUMA topology,
affinity, cgroup CPU and memory limits, OS/kernel, compiler and binary hashes.
Require at least four visible CPUs, 14.5 GiB observed RSS capacity plus
256 MiB host headroom, and 2 GiB free disk. Pin the four execution CPUs
before any timings. Each IC child has a 1,500-second cap and every rho child
a 120-second cap. Retain exact start/stop times, stderr, operations, peak
memory, all target-dependent attempts, five exclusive IC online phases,
setup costs separately, cgroup event deltas, system steal and interrupts,
frequency if exposed, and process context switches.

For point indices 0 through 4, run `IC,rho`, `rho,IC`, `IC,rho`, `rho,IC`,
`IC,rho` in fresh processes. Both arms receive the same public point. Finish
and independently replay each arm's scalar and target relation before the
pair passes. The target-online interval starts with the first target-dependent
operation after reusable IC base/index/log setup, or the first rho walk step,
and ends when the scalar is recovered and verified. Record external replay
duration separately when outside the producer interval. Fresh process setup
is retained but excluded from this online clock for both arms.

Before and after each pair, run three fixed one-million-step SHA-256
calibration samples. Require identical checksums, each triple's maximum
sample divided by its minimum at most 1.15, and the after/before median
ratio between 0.85 and 1.15. Also require zero new cgroup CPU-throttled
periods and OOM kills, observed peak RSS at most 14.5 GiB per child,
nonnegative system counter deltas, and steal ticks no more than 1% of total
CPU ticks over the pair. Any failed correctness, resource or noise condition
keeps that pair's raw receipt and makes the five-pair result incomplete.
There is no replacement point after timing. A whole-panel retry uses a new
allocation and retains the first panel.

For five passing pairs, report each point's verified `rho_online_ms /
IC_online_ms`, the median ratio, arithmetic means, every phase cost, setup
and resource peaks, and a 95% paired bootstrap interval for the median ratio
with 10,000 resamples at seed `2026101505`. A Tier 1 speedup
requires IC to be faster on all five points. Label the result with the exact
cloud allocation and observed envelope. Source-panel many-point rankings
remain separate from this primary single-target comparison.
