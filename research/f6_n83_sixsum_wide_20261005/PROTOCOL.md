# n83 six-summand wide-system feasibility gate

Registered before code changes and execution on 2026-10-05. This is a
bounded algebraic feasibility experiment for the K0 curve
`icv1-f2m83-tm6151469093347-debefd74` over
`F2[z]/(z^83+z^45+z^2+z+1)`, subgroup order
`2417851639230796216685689`, cofactor 4. The base is the standard
polynomial subspace of dimension 16. Its frozen inventory has 64,907
geometric points, 64,904 usable projected subgroup points and 32,452
signed columns. A direct source-point six-summand S3 chain has
`6*16+4*83=428` Boolean variables and `5*83=415` coordinate equations.
The current 128-bit monomial mask rejects it. The experiment tests whether
a 512-bit mask can construct the *coupled* direct system and whether a
bounded native reduction gives useful constraints on an ordinary public
target. It makes no coverage or speed claim.

Baseline: commit `fa6171e1a4bbc4ca448974feb585de49d763d327` (PR #1409),
including the two-word n83 field table, indexed valid coordinates, and
cofactor-four source/subgroup bridge. No change to the existing 128-bit
solver is permitted in this gate. Freeze the public T001 target already
used by the n83 K0 gate and a deterministic planted six-point control
from geometric source-point indices `[0,2,4,6,8,10]`. Record exact
point encodings and source hashes with the result. A planted case checks
equations and group arithmetic; it is never evidence of ordinary-query
yield.

Acceptance and falsification conditions:

1. Construct 428-variable direct equations with a native fixed eight-word
   monomial representation. Check Boolean cancellation, multiplication,
   assignment and evaluation at word boundaries 63/64, 127/128 and
   383/384. Confirm each planted equation vanishes under the six source
   coordinates and four exact intermediate sum abscissae. Check the six
   source points sum to the planted source target and cofactor projection
   matches the subgroup target.
2. Construct the same system for T001's source preimage obtained by
   inverting cofactor 4 in the odd-order subgroup. Try each rational
   4-torsion offset deterministically. Record equation count, distinct
   monomials, degree, construction time, peak RSS, and every failure.
3. Run at most one own-degree linearisation/reduction on the ordinary
   target under a hard 120-second and 8-GiB process envelope, and a
   small fixed node budget only if a sound branch search is implemented.
   Record rank, linear consequences, contradiction, or budget stop.
   A spent budget or resource failure is **inconclusive**, never a
   negative decomposition proof. Any candidate witness must be lifted
   and verified by full curve arithmetic and cofactor projection.

Stop this approach as currently implemented if direct construction exceeds
the envelope, if the planted control fails, or if the ordinary root
exhausts the envelope without a useful constraint. In that case preserve
the failed row and analyze a different coupling strategy before more
enumeration. A success here requires a verified ordinary T001 relation,
not merely a feasible matrix. A complete IC candidate, its one-target
online cost, and matched one-target rho remain unknown until relation
collection, rank, final linear algebra, target descent, scalar recovery,
and independent verification exist.

The counting capacity `C(64904+5,6)/r = 42.950379` is an upper bound
before duplicate sums; it is not a measured success probability or a
speedup. The generic reference is rho on the same public subgroup target,
with expected work on the order of `sqrt(r)` group operations, but no
matched rho run is in this gate. The stage boundary is the 120-second,
8-GiB feasibility envelope, not an end-to-end complexity boundary. All
new wall times are exploratory until an auditable host-isolation receipt
is available. Use Rust/Cargo only, one worker, preserve raw outputs and
resource failures, and do not assign an `IC1` candidate ID to this
incomplete pipeline.
