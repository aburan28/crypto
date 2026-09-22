# Tiered multi-target collision equations

Date: 2026-09-21

Status: exact ideal-model derivation, scoped prior-art audit, and finite
public-synthetic legal-chain controls.  The square-root batch scaling is known
prior art and asymptotically optimal in the generic model.  The tested storage
realization has a dependent-equation tail and no cryptographic-scale storage
feasibility.

## [T1] Prior-art result

Kuhn and Struik explicitly extend Pollard-style collision search to multiple
discrete logarithms and report roughly `sqrt(2NK)` group operations for `K`
instances in a prime-order group of size `N`.  They also identify the
`Omega(sqrt(KN))` generic lower bound as the natural open problem.
[Kuhn--Struik](https://doi.org/10.1007/3-540-45537-X_17)

Yun subsequently proves a tight generic lower bound for the multiple discrete
logarithm problem, showing that known `O(sqrt(KN))` generic algorithms are
asymptotically optimal.
[Yun, *Generic Hardness of the Multiple Discrete Logarithm Problem*](https://eprint.iacr.org/2014/637.pdf)

Fouque, Joux, and Mavromati analyze multi-user collisions and obtain
soft-`O(sqrt(NK))` attacks, including a matrix/rank interpretation of collision
relations and the effect of row sparsity.
[Fouque--Joux--Mavromati](https://eprint.iacr.org/2013/761.pdf)

Thus neither the batch square-root law nor collision equations for several
logs are novel.  This branch audits one exact archive-backed implementation and
its resource costs.

## [T2] Exact ideal distribution

Let `G=<P>` have prime order `N`, and let public targets be

    Q_j = x_j P,  j=1,...,K,

where the logs `x_j` are independently uniform in `F_N`.  A coefficient record
is

    c_i = (a_i,b_i1,...,b_iK)

and represents

    X_i = a_i P + sum_j b_ij Q_j
        = (a_i + sum_j b_ij x_j)P.

Under the ideal reference, coefficient vectors are independently uniform in
`F_N^(K+1)`.  Group labels are therefore uniform in `G`.  A collision
`X_i=X_l` gives

    sum_j (b_ij-b_lj)x_j = -(a_i-a_l).

Conditioned on a collision, the coefficient difference is uniform in the
`K`-dimensional kernel orthogonal to `(1,x_1,...,x_K)`.  Its `b` projection is
uniform in `F_N^K`, with `a` uniquely determined.  Independent collision pairs
therefore supply uniform linear equations.  The probability that `K` such rows
have full rank is

    product_(r=0)^(K-1) (1-N^(r-K)),

which is close to one for cryptographic prime `N`.

Among `M` uniform labels, the expected birthday-pair count is approximately
`M^2/(2N)`.  Reaching `K` equations therefore requires

    M approximately sqrt(2KN).

This establishes the claimed scaling only for the independent-uniform model.
It does not show that arbitrary legal addition-chain records behave that way.
The iid coefficient reference in the checker is labeled as samples, not group
additions, because materializing an arbitrary vector could cost many group
operations.

## [T3] Legal target-independent archive chain

The finite construction uses the public additive group `F_65537`.  Its initial
archive contains

    O, P, Q_1, ..., Q_K.

Those are public inputs and require zero group additions, but all `K+2` exact
records, coefficient vectors, and collision-directory entries are charged.

The current record adds one prior archive record selected by a deterministic
public seed/time function.  Each step therefore performs exactly one legal
group addition and one `(K+1)`-component coefficient-vector addition.  Every
output is appended and probed in the exact collision directory.  The selector
does not depend on any secret log or group encoding.

This construction is memory-full and time-dependent.  It is not a coalescing
rho function on the projected group state.  Its purpose is to test whether a
cheap legal recurrence supplies the iid-like collision ranks assumed by the
ideal argument.

For every exact group collision, the checker:

1. fetches the prior exact coefficient record;
2. recomputes the homogeneous relation against every public synthetic log;
3. incrementally row-reduces the `K` coefficients;
4. labels zero rows as dependent;
5. stops only at rank `K`; and
6. back-substitutes and verifies the complete public log vector.

## [T4] Finite scaling results

The experiment uses 32 fixed public trials for each `K` in
`{1,2,4,8,16}`.  Every full-rank run is capped at `24 sqrt(KN)` additions.
The independent comparator runs the same legal chain separately for every
target.  The iid reference draws uncharged independent coefficient samples.

| K | Legal multi mean | Legal median | iid sample mean | Independent total mean | Independent / multi mean |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 335.94 | 296 | 339.56 | 347.38 | 1.03 |
| 2 | 469.25 | 465 | 480.22 | 679.28 | 1.45 |
| 4 | 683.09 | 702 | 712.88 | 1,309.97 | 1.92 |
| 8 | 1,041.59 | 1,004 | 1,035.63 | 2,711.22 | 2.60 |
| 16 | 2,413.06 | 1,468 | 1,455.13 | 5,316.34 | 2.20 |

For `K<=8`, mean legal work is 1.30--1.44 times `sqrt(KN)` and closely tracks
the iid sample reference.  The speedup over independent work approaches the
expected `sqrt(K)` trend.

At `K=16`, the median remains close to the iid mean, but the legal-chain mean
is distorted by a long dependence tail.  Completion ranges from 1,093 to
18,067 additions.  Dependent collision rows have mean 104.125 and range
0--2,197.  The worst run reaches rank 15 by addition 1,618, then processes
2,197 dependent rows before obtaining the last pivot at addition 18,067.

This tail contradicts treating every legal collision as an independent uniform
row.  The asymptotic theorem remains valid for the established algorithms; the
negative result concerns this particular current-plus-history recurrence.

## [T5] Charged rank-deficiency fallback

The checker also stops shared collection at `4 sqrt(KN)`.  If rank is deficient,
it solves the free-column target logs with the charged single-target control and
uses the shared equations to recover pivot variables.  Shared work is retained;
the fallback is added rather than substituted for it.

| K | Rank-deficient shared trials / 32 | Shared plus fallback mean | Independent mean | Ratio |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0 | 335.94 | 347.38 | 0.967 |
| 2 | 0 | 469.25 | 679.28 | 0.691 |
| 4 | 0 | 683.09 | 1,309.97 | 0.521 |
| 8 | 0 | 1,041.59 | 2,711.22 | 0.384 |
| 16 | 3 | 1,730.34 | 5,316.34 | 0.325 |

For the 18,067-addition outlier, the cutoff has rank 15.  Solving its one free
target and retaining all shared work costs 4,309 additions, eliminating most of
the dependence tail while preserving exact correctness.

The cutoff is a finite engineering choice, not an optimized theorem.  Stronger
prior-art multi-walk constructions use tailored starting points and collision
graphs; this simple recurrence should not replace them.

## [T6] Rank and replay costs

Incremental row reduction charges finite-field additions, multiplications, and
inversions.  A typical independent `K=16` completion uses approximately 1,752
field multiplications for rank handling; the 18,067-addition tail uses 324,646
because 2,197 dependent rows are reduced before rejection.

Every primary record stores its complete coefficient vector, so rank testing
does not require replay.  Parent and addend indices are retained as well.  A
full audit replay costs one group addition and `K+1` coefficient-field additions
per generated record.  The worst finite run therefore charges 18,067 replay
group additions and 307,139 coefficient additions if full reconstruction is
requested.

A lineage-only compressed format could reduce coefficient bytes but would move
work into dependent random reads and replay.  It is not evaluated and receives
no claimed saving.

## [T7] 256-bit record and capacity envelope

The hypothetical exact record contains:

- a 33-byte canonical group encoding;
- `K+1` 32-byte coefficients;
- 8-byte parent and addend indices;
- 3 bytes of version/flags and a 4-byte transport CRC; and
- padding to a 64-byte slot boundary.

The external collision filter/directory adds 34 final bytes per record.  Peak
disk conservatively permits a complete second copy during recovery or merge.
The existing active-walker allocation is enlarged to hold `K+1` coefficients
per walker, and rank workspace stores `K` rows of `K+1` coefficients.

| K | Coefficient bytes | Exact slot | Peak bytes/record | Maximum records under 1 TB RAM / 100 TB peak disk |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 64 | 128 | 324 | 308,441,358,024 |
| 2 | 96 | 192 | 452 | 221,095,132,743 |
| 4 | 160 | 256 | 580 | 172,301,724,137 |
| 8 | 288 | 384 | 836 | 119,539,473,684 |
| 16 | 544 | 640 | 1,348 | 74,135,756,676 |

Peak disk binds in every row.  At `K=16`, final disk at the maximum is
49,968,499,999,624 bytes and conservative peak disk is
99,999,999,999,248 bytes.

Fewer group additions do not imply fewer stored bytes.  At `K=16`, a shared
record carries 544 coefficient bytes, while an independent single-target record
carries 64.  Using finite mean record counts, the shared coefficient payload is
several times the combined independent coefficient payload even though shared
group additions are lower.

For a 256-bit group, the asymptotic `sqrt(KN)` record requirement is still on
the order of `2^128 sqrt(K)`, vastly beyond the tens of billions of records in
this envelope.  Capacity accounting therefore does not establish a practical
256-bit attack.

## [T8] Verification and decision

The checker verifies:

- every record coefficient update and complete lineage replay;
- every collision relation against all public synthetic logs;
- every incremental rank increase and dependent row;
- every full-rank solve by exact back-substitution;
- every independent and residual-fallback completion;
- exact group-addition, directory-probe, verification-read, parent-read,
  field-operation, coefficient-byte, and peak-storage counters; and
- hash bindings for public logs, coefficient vectors, and lineage.

Artifacts:

- `check_multi_target_linear_collisions.py`
- `multi_target_linear_collision_checks.json`

The requested `O(sqrt(KN))` versus `K O(sqrt(N))` distinction is valid under
the independent-uniform collision distribution and is already established
generic-DLP prior art with a matching lower bound.  The finite legal archive
chain supports the scaling for typical small-`K` runs but exhibits severe
dependent-rank tails at `K=16`; a charged cutoff and residual fallback repairs
the finite tail.  Coefficient storage and cryptographic-scale record counts
remain prohibitive.  No novelty, private-target recovery, concrete ECDLP
speedup, or hardware-feasibility claim is made.
