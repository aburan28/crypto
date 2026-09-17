# Two-jump rho and fixed Frobenius paths on G7

Preregistered engineering experiment on the local AWS g7.2xlarge / RTX PRO
4500, GPU UUID GPU-827e82b5-4739-d339-d95e-b7214337553e, at 165 W. The objective
is 12 billion complete ECC2K-130 scalar iterations/s on one GPU. The immutable
selected reference is build/ecc2k130-local-packed, SHA-256
8c4ed76a152cc135280a68fd3559c3bf71bef726ab73eae0e04014c2e27c6d02. Its generic
work boundary is sqrt(n/262), generic work ratio 1, and full-DLP S is null.

## Hypothesis and variants

The selected rho map chooses eight Frobenius jumps, j=3..10, from three bits of
half the normal-basis X weight. The scalar multipliers for jumps 3 and 4 alone
generate all of F_ell*, so choosing j=3 or 4 from one weight bit does not trap
the walk in a proper scalar subgroup. A preliminary 200-target planted-solve
study found a 170.330 mean collision time for jumps 3/4 versus 162.695 for the
selected eight jumps, a 1.04693 collision-work ratio. Freeze this as the initial
algorithmic cost estimate and reproduce it with a larger deterministic cohort.

Build two candidates from an isolated selected-source copy. two-jump-indexed
changes only the partition and still uses the selected indexed/shared Frobenius
network. two-jump-immediate emits separate constant-mask sigma^3 and sigma^4
paths, eliminating indexed mask reads and the shared sigma table. The selected
binary remains the timing reference. No curve arithmetic, DP predicate, seed,
state format, checkpoint format, canonicalization or endpoint record changes.

Before timing, require exhaustive basis-vector equivalence for both fixed sigma
paths, dense random pair equivalence, complete GPU state and CPU endpoint replay
for each candidate, bidirectional resume, DP admission with zero drops, and
memcheck, initcheck and synccheck of the changed walk. Require a deterministic
planted-solve cohort with every target recovered; report its collision-work
ratio separately from throughput.

Time three interleaved repetitions of selected and both distinct verified rows
at 524,288 workers, B16, 1,024 steps and four launches, or 34,359,738,368
complete updates per sample. Preserve all samples and hardware identity. A row
advances only if its paired log-ratio Student-t 95% interval is wholly above 1.
Also report throughput divided by the frozen collision-work ratio; never label
raw iteration gain as solution-time gain. Promotion requires five fresh matched
pairs in benchmark and DP34 workloads, both confidence intervals above 1, a
fresh larger collision cohort, complete DLP recovery, and the full validation
suite. Otherwise retain the selected runtime and the negative evidence. No
hardware, power, driver, service or cloud-resource changes.

## Result

The round closed without a qualifying candidate. All nine timing samples
completed, no external GPU process was observed, and no sample was excluded.
The selected median was 6.279988 B/s. The indexed two-jump map measured
6.273318 B/s, paired ratio 0.992430 with 95% CI [0.964207, 1.021479]. Its
collision-adjusted rate is 6.192320 B/s. Fixed immediate sigma paths measured
5.418604 B/s, ratio 0.844542 [0.778245, 0.916487], and adjusted rate 5.348642
B/s. Neither raw rate qualifies, and both adjusted rates trail the selected
reference. Retain the selected eight-jump executable.

The matched collision study recovered all 6,000 planted GF(2^23) DLP targets
with no invalid points. Two jumps required 172.637 mean iterations versus
170.408 for eight jumps, a 1.013080 collision-work ratio. The exact multiplier
orders jointly generate all of F_ell*. The fixed helpers passed 16,646 GPU pair
checks against independent basis routing. Both implementations produced
identical complete states and 38,885 DP records across IDs 0 and 139, with 64
CPU endpoint replays and bidirectional checkpoint continuation. Each candidate
passed memcheck, initcheck and synccheck. The immediate walk used 80 registers,
16 stack bytes, 16 spill-store bytes and 24 spill-load bytes, versus 8/8/8 for
the indexed walk. Generic work remains sqrt(n/262); measured candidate ratio is
1.013080 and full-DLP S remains null.
