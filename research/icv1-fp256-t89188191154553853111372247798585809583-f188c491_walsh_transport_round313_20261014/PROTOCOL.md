# P-256 Dickson--Walsh global transport, round 313: preregistration

Date preregistered: 2026-10-14

## Question and boundary

Test whether the complete depth-18 Dickson trace fibre underlying
`FB1h2f8621cda105` carries a *global* elliptic-log identity that Round 312's
local additive atoms could not see.  The registered curve remains
`icv1-fp256-t89188191154553853111372247798585809583-f188c491` and the
registered factor base remains unchanged: 131,458 negation-folded columns,
262,916 signed points, and relation arity 17.

One end-to-end unit remains

```text
S = charged P-256 group-addition equivalents / sqrt(n),
```

with Pollard rho fixed at `S=1.3`.  Before this experiment, the exact
free-relation-oracle lower projection is `S=512.752864`, or `394.425280`
times rho, because all 131,458 factor-base logarithm classes remain
independent.  The executable S17 projection is `2^148.080980` P-256
field-multiplication equivalents (FME) per usable row and `2^165.155613`
FME for 138,031 rows; projected materialized attack storage is
`2^148.926255` bytes.  Those are projections, not measurements.

Parity means a verified, cold, one-target total `S <= 1.3` with every phase
charged.  A static-rank reduction, a local degree reduction, a stage below
rho, or an asymptotic `O(sqrt(n))` label is not parity.

## Frozen inputs

Reject the run on any dependency, curve, factor-base, enumeration, point,
or schema mismatch.

- Round-2 factor-base result SHA-256:
  `dcaa1f66f5f89757a3670a4a8e158ca6d5c360b7ca119d384d2c4a0b2794c2c8`.
- Round-19 selector result SHA-256:
  `3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114`.
- Round-312 additive-rank result SHA-256:
  `0cdec91f0097b498d1642c80c358aadac279c8a0bae951ed433381e1bf53db8e`.
- Factor-base implementation SHA-256 at preregistration:
  `a8740201f5f5c7dc56489827eb6f6a98d946338a512eb943741979f43eaf58c5`.
- Full FB1 digest:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`.
- Point-stream digest:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`.
- Depth 18, fibre length `N=2^18`, root exponent
  `0x2b6fdc73dc04e7667129`, and terminal
  `0x5b17195299a3158b93389ad04c776fff2a8bb23ca8659b5b0c2b75f9009b65e5`.

The implementation may add a read-only indexed-enumeration API to the frozen
builder, but `build("dickson-torus:depth=18")` must remain byte-identical and
all existing builder tests must pass.

## Complete global family

Let the builder's canonical torus order be `i=0,...,N-1`.  Put the low-y
factor-base representative in slot `i` when that trace abscissa lifts to the
curve, and the group identity in every non-lifting slot.  For every binary
character `s in F_2^18`, compute the exact group point

```text
W_s = sum_(i=0)^(N-1) (-1)^<s,i> P_i.
```

Compute all `N` values by the complete group-valued Walsh--Hadamard transform.
No character, sign branch, or non-lifting slot may be sampled or discarded.
The implementation must verify that the active index set affinely spans
`F_2^18`; this certifies that distinct characters give distinct coefficient
rows on the actual columns (rather than merely distinct formal fibre rows).

Canonicalize outputs only modulo elliptic negation, retaining the orientation
needed to reconstruct the coefficient row.  Exhaustively report and replay:

1. every identity output `W_s=O`;
2. every collision `W_s=+/- W_t`, `s<t`; and
3. every match `W_s=+/-[m]P_j` for every actual column and
   `m in {1,2,3,4,8,16}`.

Normalize each nonzero coefficient row by its integer gcd and global sign,
deduplicate exact rows, and compute the exact rank modulo both `2^61-1` and
the P-256 subgroup order.  Promotion credit is the smaller rank.  Every
reported relation must be replayed in the P-256 group from its implicit Walsh
row; a hash/key match alone is never a relation.

The experiment is complete only for this registered Walsh family.  Absence
of a collision is not a proof that no other nonlocal identity exists.

## Controls and accounting

1. Exhaust complete Walsh transforms in small prime cyclic groups, compare
   the optimized collision/multiple detector with a direct coefficient-row
   reference, and require zero false positives and false negatives.
2. Run a scalar-labelled P-256 control with planted Walsh identities and
   small-multiple matches.  Replay every control relation in the actual group.
3. Rebuild and verify the complete registered factor base independently.
4. Emit transform additions, field operations used for normalization,
   key comparisons, matches, exact rows, rank operations, logical RAM, peak
   RSS, disk traffic, wall time, and deterministic ordered digests.

The Walsh transform is cold preprocessing and is charged to a one-target
attack.  A zero-rank result leaves every Round-312 boundary unchanged.  A
nonzero row count without nonzero exact rank receives no relation-count
credit.

## ecbench gate

`ecbench` currently executes prime curves only for `p < 2^62`; it cannot run
the registered P-256 group, and an analytic or projected record must not be
inserted into an `ecbench` session.  This round therefore has a fail-closed
promotion rule:

- If the complete Walsh family has zero certified rank, there is no new DLP
  method to benchmark.  Record `ecbench candidate = ineligible`, preserve the
  exact P-256 result, and do not disguise a stage screen as an ecbench solve.
- If it has nonzero rank, first implement a word-size Dickson--Walsh analogue
  as a new, replayable `ic.pipeline` candidate, freeze at least four registered
  prime-curve sizes, and run matched `rho.negation`, unmodified IC, candidate,
  and A/A arms through `ecbench`; every recovered scalar must replay.
- No scaled ecbench result transfers to P-256 automatically.  It supplies an
  exponent/control measurement only.  A P-256 parity claim still requires an
  executable wide path and a verified one-target solve, neither of which is
  assumed here.

The ecbench promotion session, if eligible, must use an `ecbench.spec/v1`, at
least eight targets per curve, deterministic seeds, cold setup, complete
phase accounting, `verify --replay-all`, and per-curve `S` ratios.  Fit total
time and memory exponents over at least four sizes.  Missing or unpriced work
makes parity unset, not true.

## Promotion and stop conditions

All of the following are required before an actual unplanted full-depth P-256
relation attempt:

- zero false positives and false negatives on every complete control;
- exact P-256 replay of every reported Walsh relation;
- certified remaining log rank and structured residual degree at most five;
- projected usable relation below `2^103` FME;
- projected 138,031-row collection below `2^120` FME, adjusted by certified
  rank and duplicate allowance;
- projected materialized storage below `2^50` bytes;
- a verified ecbench-eligible complete method whose cold one-target `S` is at
  most the matched rho reference on growing sizes; and
- no sampled, omitted, probabilistic, or hash-discarded branch counted as an
  exhaustive search.

If the Walsh rank is zero or the updated free-oracle boundary remains above
rho, stop this candidate and publish the negative result.  Do not tune an
output-dependent subset of characters.  The next hypothesis must be a
different complete nonlocal identity family or a mechanism that certifies
more than one independent primitive quotient equation per charged event.

## Required artifacts

Commit this protocol before implementation.  Commit the native Rust runner
before execution.  Preserve the exact command, revision, host capsule,
isolation receipts, deterministic result JSON, hashes, failed runs, tests,
requirement-to-evidence table, boundary table, ecbench eligibility decision,
dashboard update, and pull request.  No Python cryptographic or measurement
code is permitted.
