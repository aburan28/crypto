# Full-point pair-sum query for the compact-orbit producer

Status: preregistered algorithm and development gate. No point-sum timing or
target outcome had been measured when this protocol was committed. The control
source is main `0e7f5b482da270227e114b9020d958b4ea3b5578`, whose
`examples/koblitz_orbit_dlp_s3_batch.rs` SHA-256 is
`702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38`.
The point-sum candidate, its tests and any timing result belong in this PR;
the initial protocol commit must precede implementation and measurement.

## Hypothesis and exact boundary

The current W64 query solves S3 for each indexed pair sum and public target,
although both are rational curve points. For an indexed full point
`R=(u,y_R)` and public `Q=(v,y_Q)` with `u != v`, let
`d=(u+v)^(-1)`, `lambda=(y_R+y_Q)d`, and `h=ud`. The two possible partner
x-coordinates are

`z_+ = lambda^2 + lambda + u + v + a`,
`z_- = z_+ + h^2 + h`.

These are `x(R+Q)` and `x(Q-R)` by the binary-curve group law. They must
equal the unordered S3 roots for every regular candidate; zero/equal-x and
any unmaterialized indexed point use the existing S3 fallback. Build full
pair points for the *existing* swap/Frobenius states; do not enlarge the
factor base, change state order, root table, rank policy or rho. In a W64
query, charge one Montgomery batch inversion, its prefix/reverse
multiplications, and both per-candidate products. Charge both signed
pair-point constructions, their memory and all exceptional fallbacks in the
same cold process. A lower field-multiplication count is a mechanism, not a
speed claim.

The reference is the unchanged W64 S3 query with blocked root prefilter on
the identical base and public Q. The matched signed-Frobenius normal-basis
batched rho v3 is the method boundary. The candidate is successful as an
engineering change only if (1) it retains full rank, every relation and
scalar independently replays, and no support hit is lost; (2) the paired
95% interval for **complete-process CPU candidate/control** lies below one;
and (3) cold peak RSS stays below 5 GiB. A method crossover requires the
same complete CPU candidate/rho interval below one, a disjoint-Q repeat,
and the still-missing calibrated operation-unit and n83/n131 transfer gates.
Do not infer that from an index or query-stage timer.

## Development and stop gate

First prove the algebra on every ordered pair of rational GF(2^5) points,
including zero/equal-x and both signs. On the retained, already-public
point-only development corpora, use n37/L1024/K42, n41/L1/K85, and
n53/L1/K220 from `compact_ir_cold_gap_20260930/CONFIG.json`. Their point
SHA-256 values are respectively
`3e249b1ddb41f7c3f6e73cbe893d12db87ce58bebcf92013c8b39d2ee0ec05d6`,
`fff82518e01fa70ef0facfc8273ed112eb2ebfc18ba5c9a71d472d3af1aeb739`,
and `d71cd294f88cc87f85f0da9f6f3bd6585072f9128b2542047127307fb0c043c5`.
Use rank seed 7 and identical W64 state order. Independently replay all
point sums, rank rows, target group equations and `[d]G=Q`; compare base
hash, index/root counts, full-rank attempts, hit/miss counts and scalars.
Differences in first witness or probe count must be retained and explained.
These disclosed Q are for correctness and stage diagnostics only.

Stop and retain a negative result if any regular root set differs, a point
sum fails replay, a formerly solved Q fails, a scalar differs, a cell
exceeds 900 seconds or 5 GiB, or extra point-index cost dominates the
query savings in all three development cells. Archive failures rather than
selecting a favorable cell. If correctness passes and at least one cell
has a plausible complete-cost reduction, freeze candidate source and
generate a fresh disjoint point-only Q stream *before* any confirmatory
timing. The follow-on run must use five rotated cold blocks per cell with
control A/candidate/rho/control B, one pinned core, a reserved-core
contention monitor, wait4 user+system CPU, full rank, linear solve and
every target recovery in each child. A/A median must lie in [0.9,1.1],
its paired 95% log-t interval must contain one, and monitor contention
must be zero. Otherwise that cell's timing is ineligible. Publish source,
input and binary hashes, host manifests, all raw arms and failures,
independent second-host replay, paired intervals, and the canonical
scoreboard update in an outcome PR. No new Q or performance conclusion is
claimed by this protocol alone.
