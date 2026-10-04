# Frozen n37 one-target online control

This is a bounded, native Rust control for the compact-orbit index-calculus
pipeline on the public Koblitz `a=0` curve over `GF(2^37)`. It does not model
ECC2K-130 directly. The question is whether the existing source-curve
factor-log method can recover one previously unseen point after reusable
precomputation, with complete target-phase accounting, and how that interval
compares with the repository's strong signed-Frobenius rho on that **same**
point. The derived generic rho reference is on the order of `sqrt(r/(2n))`
group steps for signed Frobenius; no tuning decision will be made from one
noisy wall-time sample.

The frozen inputs are `n=37`, `a=0`, seven signed Frobenius orbit columns
(expected 518 distinct points before folding), deterministic public x-scan
factor-base construction, rank seed `3737001`, rho batch seed `370041`, and
public hash-to-curve seed `370413`. First run native
`koblitz_rho_fixture 37 0 signed_frobenius 1 strong 370041 hash:370413`
to materialize the public point `Q` and its independently recovered scalar.
The IC input is the resulting `[x,y]` point only. The IC producer is
`koblitz_orbit_dlp_fast_online construct:37:0:7 <point-file> 3737001
<target-output>`; no known scalar enters it. Both are built from the same
commit and release profile. A 15-minute wall cap and 16 GiB memory ceiling
apply to each arm; record a timeout or failure as a result, never as a win.

The IC online interval starts with the first target-dependent query-hash
operation after base/index/rank/log setup and ends after scalar replay. Its
five exclusive costs are query, PDP, relation check, descent, and recovery
check. Rho starts at its first target-dependent walk and includes recovered
scalar verification; record its setup separately. Preserve full-rank attempts,
failures, base size, folded columns, peak RSS, public point, scalar, exact
source and input hashes, binary hash, build toolchain, raw JSON, and every
phase. Independently replay the IC relation, base logs, and both recovered
scalars using the general binary-curve group law. The pairing passes only if
both programs solve the identical `Q`, both scalar replays pass, IC rank is
full, and the five phases sum to its online interval. A time ratio is a local
diagnostic unless host-level CPU isolation is auditable. The result cannot
by itself support an ECC2K-130 extrapolation or an end-to-end cold-cost claim.

Before measuring, commit this protocol. Any implementation issue discovered
during the run must be fixed in a new commit and the affected measurement
repeated with the changed source hash.
