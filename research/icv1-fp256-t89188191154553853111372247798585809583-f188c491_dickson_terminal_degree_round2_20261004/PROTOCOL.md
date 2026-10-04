# P-256 Dickson terminal-fibre degree protocol, round 2

Date frozen: 2026-10-04

Round 1 rejected the affine-bitbox replacement: its paired `S3` solving
degree was 5 or 6 while the accepted terminal-zero Dickson chain stayed at 4.
This round keeps that chain and asks whether a different nondegenerate terminal
trace/coset lowers solving degree without reducing factor-base cardinality.
The registered curve remains
`icv1-fp256-t89188191154553853111372247798585809583-f188c491` (P-256).

## Hypothesis

For the membership chain

```text
z0 = x,
z_(j+1) = z_j^2 - 2,
z_(t-1)^2 - 2 = c,
```

some full Dickson fibre with `c != +/-2` and at least the terminal-zero lift
count has paired common-positive `S3` F4 solving degree at most 3 at both
depths 4 and 5.  The matched terminal-zero reference has measured degree 4.

This is a strict success gate.  A lower matrix width or operation count at
degree 4 is engineering, not a degree result.  A smaller factor base is
inadmissible even if its degree is lower.

## Frozen inputs and sweep

- dependency: the two commits of round 1, ending at local commit
  `65330510b` and remote PR `#1329`;
- toy prime: `p = 1151`, with `p + 1 = 2^7 * 9`;
- toy curve: `y^2 = x^3 - 3x + (b_P256 mod 1151)`;
- depths: `t = 4, 5`;
- terminals: every `c in F_p` for which `D_(2^t)(x) = c` has exactly
  `2^t` distinct roots in `F_p`, excluding `c = +/-2`;
- cardinality gate: retain only terminals whose liftable-root count is at
  least the count for `c = 0` at the same depth;
- relation: `S3(x1, x2, xR) = 0`;
- factor membership: the triangular Dickson chain plus the curve-lift
  equation for each summand;
- monomial order: graded reverse lexicographic;
- F4 degree bound: 12;
- per-system stop: 30 seconds;
- targets: the first three affine points in increasing `(x,y)` order whose
  x-coordinate is a two-summand decomposition for both the candidate fibre
  and terminal-zero reference, using exhaustive signed-point addition.

The candidate and terminal-zero system are run on exactly the same target.
Repeated summands and doubling are allowed in both the exhaustive reference
and `S3`.  A terminal with no common positive target is retained as
inconclusive and cannot win.

The winner at each depth is the terminal minimizing, in order:

1. median completed `solving_degree_max`;
2. median columns to the last productive step;
3. median field multiplications;
4. negative lift count (more columns wins);
5. the integer terminal.

The terminal-zero reference is never eligible to be reported as a new
candidate.

## Correctness and boundary table

Exhaustive signed-point addition establishes every expected positive verdict.
A completed F4 run must remain consistent; a certified inconsistency is a
falsification.  Timeouts and degree truncations are preserved and excluded
from the winning median.

The primary unit and boundary are paired F4 solving degree.  The result table
has one row per `(depth, terminal)`:

| depth | terminal | full roots | liftable roots | common targets | correct | complete | solving degree min/median/max | columns median | field ops median | degree / c=0 | classification |
|---:|---:|---:|---:|---:|:--:|:--:|:--|---:|---:|---:|:--|

No end-to-end `S`, rho ratio, or ECDLP exponent is reported.  This is a solver
stage diagnostic over a small-prime analogue.

## Frozen P-256 coset search

Independently of the toy-degree outcome, construct eight depth-18 P-256
Dickson fibres.  For `k = 0,...,7`, define the root exponent as the low 78
bits of

```text
SHA256("icv1-fp256-t89188191154553853111372247798585809583-f188c491/dickson-terminal/v1/" || decimal(k))
```

interpreted big-endian, then set the low bit to one.  The odd exponent makes
the terminal nondegenerate and fixes a full `2^18` fibre.  Select the fibre
with the most liftable abscissae, ties by smaller `k`.

The accepted terminal-zero `FB1h2ea06bef7f7a` remains the reference.  Emit the
selected candidate with the existing `dickson-torus` FB1 family, registered
ICV1 slug, root exponent, terminal trace, and
`prime-affine-x-plus-one-shift-sign/be33/v1` point keys.  Rebuild it, verify
every point row, and record all eight counts, FB1 preimage/digests, and the
exact `m=17` signed domain and Poisson diagnostic.

The P-256 cardinality winner inherits no toy-prime degree measurement.  It may
be called a larger factor base if it exceeds 131,239 columns, but it may be
called a lower-degree candidate only if the toy sweep meets the strict degree
gate and only with an explicit small-prime qualifier.

## Stop conditions and evidence

Stop after all eligible toy terminals and all eight P-256 exponents.  Preserve
every regression, timeout, missing target cell, identity mismatch, and failed
hypothesis.  Commit compact raw JSON, a human-readable decision table,
artifact hashes, native commands, and verification receipts.  The full P-256
point dump remains deterministic derived data and need not be committed.

All implementation and experiment execution is native Rust.  No P-256
summation-polynomial Gröbner basis is attempted or implied.
