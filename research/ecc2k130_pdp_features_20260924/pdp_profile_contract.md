# Experiment contract: product-span and affine-tail PDP screen

Status: HYPOTHESIS / TOY-EVIDENCE; diagnostic, not an index-calculus speed claim.
Frozen before running on 2026-09-24.

Candidate: factor-base selection and target routing should use the quadratic
coefficient rank and exact affine consequences of the specialized Boolean PDP,
while retaining verified natural yield. Small systems that are cheaply refuted
are not necessarily useful relation generators.

For two summands on K0, t=x(R), V=<u_i>, W=<u_i*u_j>, the cross-block coefficient
is L_t(u_i*u_j), where L_t(z)=z^2+t*z. For t nonzero,
rank(L_t|W)=dim(W)-1[t in W]. This is a directly checkable linear-algebra identity,
not a degree-of-regularity assertion. A basis change within V preserves this
rank. The full affine tail is computed by eliminating quadratic columns from
the Boolean coefficient rows, retaining constants.

Parameters: K0/F2^n for n=7,9,13,17; dimensions 3,4,5,6 respectively. Each size
uses the polynomial-coordinate subspace, two independently random subspaces,
and two random nonzero field multiples of the polynomial subspace. Seeds 17
and 937 generate bases; 32 natural nonzero prime-subgroup targets are shared
by all bases per size. Four planted nonzero prime-subgroup targets per base
are correctness controls and excluded from natural-yield estimates. Infinity
and x=0 exceptional targets are excluded explicitly. No isogeny is measured.
The largest prime divisor r and generator order are checked; no secret scalar
is an input to the profile or solver.

Reference: exhaustive evaluation of S3 on all ordered x-pairs in V^2, followed
by rational point lifting and group verification. Candidate/check: for each
x1, solve the exact linear system in x2 over F2, recovering all roots. Both
must return identical full root sets. A separate exhaustive group point-pair
census checks the lifted result. Fixed full-root enumeration intentionally
measures an audit workload, not time to the first relation.

Controls: planted targets must lift; random subspaces break low product span;
scaled polynomial subspaces preserve product-span dimension; an invertible
change of subspace basis must preserve the algebraic ranks and root x-pairs.
Every original polynomial solution must satisfy the extracted affine tail.
An inconsistent affine tail must have no algebraic roots.

Metrics: exact factor-base size, product-span dimension, squaring closure
defect, quadratic rank, affine rank/inconsistency, algebraic/rational/group
root counts, target membership in W, exhaustive assignments, number of fiber
systems and row XORs, CPU/wall diagnostics. Native units remain separate;
no fabricated conversion to group additions or rho ratio. Full-DLP costs,
FFD, solving degree and rho speedup are unmeasured/null.

Correctness success: zero disagreements on all cases, rank identity exact,
all planted cases group-valid, all basis-change invariants preserved. Research
success is only detecting structural variation that could support a later
blinded solver panel; no runtime-based policy is selected here. Null: the
features vary but natural useful yield/work does not improve, or easy instances
are primarily UNSAT. Failures remain artifacts and falsify the relevant claim.
Budget: 5 minutes, local CPU; no external service. Output into a new directory;
no edits to existing solver baselines, frozen evidence, or boundary ledgers.
Reproduction: python3 pdp_profile_probe.py --source-repo /Volumes/SSD990/autolab/ecc2k130-normalx-e2e-20260923 --out profile_run_01

## Amendment 01 (before rerun)
The first run stopped at positive-control construction: degree-7 polynomial
V has only three curve points and its full pair census contains no nonzero
prime-subgroup two-sum. Preserve this base as a negative domain, not discard it.
Control generation now attempts at most 512 draws and, if fewer than four
controls exist in that sample, enumerates the full point-pair support. Record
zero available controls if support is empty; otherwise draw controls from
that exact support. This changes no natural target or factor base.

Runtime requirement: Python >=3.10 (`int.bit_count`). Run 02 was a Python3.9 environment failure; run03 uses /opt/homebrew/bin/python3.11 without an algorithm change.
