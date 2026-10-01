# Frozen design: charged target-blind selection of degree-7 m3 bases

Status: **protocol only; no candidate, target, timing or outcome has been run**.
This protocol is committed before source implementation and before generation
of its new challenge Q or relation-target streams. A separate source-lock PR
must pin the exact producer/verifier and imported dependency hashes, merge,
then independently generate and verify the fixtures before any scored run.
The [machine-readable config](CONFIG.json) fixes seeds, labels, caps and score.
The earlier [four-policy result](../m3_four_policy_20260930/RESULT.md) is
negative for a descendant-native cold-cost advantage. The subsequent
[probe-order bound](../m3_probe_order_bound_20260930/RESULT.md) makes exact
third-factor reordering an unattractive next 1,024-target lever on those
bases. This experiment changes **base composition** instead.

## Hypothesis and reference

At equal useful size, choosing between two independently sampled,
quota-matched bases using exact target-independent three-sum support and
first-witness base-row rank may reduce the complete cold cost to recover one
previously unseen Q. Compare it with the first candidate base from the same
generator, not with an unmatched or smaller base. The boundary is a maximum
`C(10,3)=120` unordered triple positions for eight base points, hence at most
`120/420` coverage over nonzero targets before duplicate sums. The accepted
reference is the prior two-seed degree-7 four-policy panel, whose source
base was cheapest in all four pairs; the new same-Q control is the actual
cost reference for this experiment. This toy result cannot establish an
ECC2K-130 or Pollard-rho crossover.

Use the same `GF(2^21)` polynomial, source `E₀`, order-421 subgroup, cofactor,
degree-7 line and normalized forward/dual maps as the prior panel. Recompute
and verify all geometry; use no saved base coordinates as inputs. The
challenge scalar is `1 + SHA256(secret label) mod 420`; it is audit-only and
must not enter either solver. Abort and preserve a collision receipt if the
new Q equals the prior panel's public Q. It is impossible to demand that all
relation targets be disjoint from the old toy group of only 421 points, so
independence here means new labels and coefficient streams, not impossible
pointwise disjointness.

## Candidate generation and score

For each of four frozen seeds and roles `source` and `leaf`, hash the base
label with seed, role, candidate index 0 or 1, and trial to a 21-bit x.
Follow the prior producer's lift, cofactor projection, duplicate rejection,
source-orbit classification and dual pullback rules. Accept exactly two
useful points from each of the first four sorted nonzero source signed-
Frobenius orbits. Cap **each** candidate scan at 4,096 x trials; keep a
failed/incomplete candidate in the archive. Candidate 0 is the deterministic
quota-matched control. The selector constructs both candidates and pays for
both, even if it chooses candidate 0.

Score each candidate only on the entire 420-point nonzero **source** subgroup,
before reading the new Q or any held-out coefficient. For a leaf candidate,
score its normalized dual pullback; the degree-7 subgroup isomorphism
preserves its support and row structure. Construct the 36 unordered pair
sums. For each third-factor index `k` in base order and pair `(i,j)` in
lexicographic order, compute the 288 pair-plus-third sums. Store the first
witness for each supported target. Build its eight base-point multiplicity
row modulo 421, then compute row rank and distinct row count independently
of the challenge coefficient. Charge all **324 group additions per scored
candidate**, their field arithmetic, row operations and memory. A candidate
is eligible only if these first-witness base rows have rank eight. Choose
the eligible candidate with larger distinct support, then more distinct
rows, then lower candidate index. If neither is eligible, record a selector
failure; do not substitute a third candidate or tune the rule.

The control policies use candidate 0 directly and do not pay for candidate
1 or scoring. The selected policies pay for both base scans and both scores.
Map each selected source base forward to `transported` and each selected
native leaf base back to `pullback`; do the same for controls. All four
policies in both arms have eight useful points and identical orbit quotas.
Keep source/transported and native/pullback hit, first-witness, rank and
scalar trajectories paired on the same Q. No target-dependent data may enter
base selection.

## New held-out streams and charged recovery

Split the ten complete source signed-Frobenius orbit IDs into A (first five)
and B (last five), exactly as in the prior panel. For each holdout, hash the
new target label with holdout, trial and `u` or `v`; form
`T=[u]G+[v]Q`, with `v` in 1..420, accept only nonidentity targets in that
holdout, and stop at 512 accepted targets or 8,192 draws. Charge rejected
draws. Use identical accepted coefficients for source and leaf and for
control and selector. Preserve all 512 outcomes, including post-rank audit
tails, and verify forward covariance independently.

For each seed, holdout, policy and arm, build the same complete 36-entry
pair table and use the same first-lexicographic pair / first-third residual
solver as the prior panel. Record every miss and witness, tested thirds,
rank trajectory, and the first full rank-nine solve. Verify every base log
and `[d]G=Q` using an independent point law. Stop the **charged** cold cost
at the verified scalar, but preserve the later audit tail separately.
Charge field setup, degree-7 kernel/forward/dual work when that policy needs
it, all base scans and score computation, generator/Q and target draws,
pair table, failed residuals, row reduction, solve and verification. The
control and selector each receive all costs they would incur if run alone;
shared audit work is excluded from both. Report field multiplications,
squarings, inversion calls, group additions, modular row operations, CPU,
wall, peak RSS and phase costs; field multiplications are the primary **toy
operation diagnostic**, never an ECC2K-130 speedup or calibrated `S`.

An independent verifier must regenerate both candidate bases, all exact
support/first-row scores and choices, the new Q/streams, four-policy map
covariance, every first witness and rank step, and each recovered scalar.
It must use independently written support enumeration and row elimination
and a separate bit-polynomial point law on the generator, Q, every selected
base point and edge targets. Reject changed raw bytes and mismatched hashes.
Preserve every source/input/binary hash, command, host, raw case, failure,
verification receipt and charged ledger in the outcome PR.

## Decision, caps and transfer

There are eight seed/holdout comparisons for each of the four policies:
source original, transported, descendant native and pullback. Treat the
**eight source** and **eight native-leaf** pairs as the primary selector
comparisons, and show their mapped transported/pullback pairs as covariance
checks. A selector advantage
requires verified rank nine and Q in **every** control and selected cell,
strictly lower selected/control cold field multiplications and no later
first full rank or lower 512-case hit count in all eight source and all
eight native-leaf pairs. Otherwise the result is mixed or negative for this
frozen selector; keep every pair. A missing arm, cap exceedance or OOM is
`CENSORED`; a witness/map/scalar mismatch is `FAIL` and blocks any advantage
claim. Never rerun or choose a favorable seed after opening the outcome.

The producer has 240 seconds and 1 GiB, the independent verifier 360
seconds and 1 GiB. A complete negative at this toy size rules out only this
two-candidate exact-support selector under this budget. A positive is merely
permission to preregister a Koblitz-family m31 test, then the m83 confidence
gate and matched-rho full-cost comparisons before n131 transfer. The
canonical scoreboard receives the outcome and its narrow claim in the same
PR; method-level `S`, crossover and n131 yield remain null.
