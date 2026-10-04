# Frozen n41/n53 one-target cold and online control

This experiment tests the current merged source-curve compact-orbit pipeline
on the public Koblitz `a=0` curves over `GF(2^41)` and `GF(2^53)`. Its primary
question is whether one fresh point can be recovered with full factor-base,
index, rank, final linear algebra, target, and scalar-verification costs
charged, compared with strong signed-Frobenius Pollard rho on that **same**
point. The online interval is reported separately after reusable IC setup.
Neither a solver-stage gain nor an online-only ratio is a cold DLP win.

The two frozen cells are `(n,K,rank_seed,rho_batch_seed,public_hash_seed)` =
`(41,85,410041,410041,41261004)` and
`(53,220,530053,530053,53261004)`. `K=85` and `K=220` are the published
one-target cold controls in
`research/notes/ecc2k130/compact_ir_cold_gap_20260930/RESULT.md`; this is a
new source snapshot and point, not a replay of that panel. The public point
comes from `koblitz_rho_fixture <n> 0 signed_frobenius 1 strong <seed>
hash:<public_hash_seed>`, using the fixture's hash-to-curve and cofactor
projection. Fix the seed now and use the first resulting nonidentity subgroup
point. The scalar is not an input to IC. Search the committed corpus for the
seed and, once Q exists, check that it is absent from the named prior panels;
an accidental duplicate is a retained failed eligibility check, not grounds
to choose another seed after measurement.

For each cell take exactly **six** fresh-process pairs on the one frozen Q.
Pair 1 runs rho then IC to create the point file; pair 2 runs IC then rho,
alternating through pair 6. Both programs use release binaries built from
the same frozen commit, one thread, default 32-lane/4-distinguished-bit strong
rho, and no cross-target cache. IC uses
`koblitz_orbit_dlp_fast_online construct:<n>:0:<K> <Q-file> <rank_seed>
<target-output>`. Record actual nonidentity base points and signed-Frobenius
columns, every rank attempt including failure, final rank, point witness,
all five target phases, rho walk steps and collision state, wall/CPU/RSS,
source/input/binary hashes, toolchain and raw stderr. Each arm has a
900-second wall cap and a 16 GiB observed-RSS acceptance gate. Preserve all
timeouts, nonzero exits and failed replays; do not exclude an unfavorable pair.

The IC online clock begins at `target_query_begin` after base/index/rank/log
preparation and ends at `recovery_check_end`. Its five exclusive phases must
sum to `online_ms`. The rho online interval is `walk_ms + validation_ms` from
the first target-dependent walk through scalar verification; fixture point
generation and process launch are excluded. For the supplementary in-process
cold interval, charge `IC setup_complete_ns / 10^6 + IC online_ms` and
`rho setup_ms + walk_ms + validation_ms`. Split IC cold setup into exclusive
base construction, normal-basis/S3 preparation, root-index build, rank query,
rank PDP, rank relation check, rank matrix work, and final relation LA;
include remaining target-independent input/control overhead explicitly so
the exclusive phases sum to the charged cold interval. Rho setup includes
its jump-table construction. Archive outer process wall separately; it is
not a substitute for the in-process phase sum.

Use the independent Rust general-curve replay to rebuild every rank row and
the final modular solve, check every base-point log, sum the four target
points, and verify both scalars on the same Q. It must check the frozen public
hash seed and resource gate and reject an altered scalar or relation row.
The replay is outside both timed intervals and its cost is reported
separately. An exact base/rank transcript may be reused across repetitions
only if its content hash is identical and every target/scalar is still
independently replayed.

This protocol is committed **before** Q or any new timing outcome is
generated. Implementation and binary/source hashes will be committed next;
then the fixed panel runs without target/K substitution. The acceptance gate
for a measured pair is full rank, both verified same-point scalars, complete
exclusive phase accounting, and both arms inside the resource cap. Report
all six paired online and cold ratios, median, spread, and failures for each
cell. These macOS wall times are exploratory until an auditable host-isolation
receipt exists. A cold loss at either fixed K constrains that policy only; it
cannot establish a global index-calculus no-go or an ECC2K-130 transfer.
No multi-target result is run or used before this one-target gate closes.
