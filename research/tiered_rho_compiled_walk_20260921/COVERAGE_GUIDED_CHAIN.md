# Coverage-guided legal addition chain

Date: 2026-09-21

Status: explicit target-independent finite mechanism with exact controls and a
parametric storage-cost audit.  The finite mathematical signal is positive;
the complete planning cost is quadratic and prevents a scalability or speedup
claim.

## [G1] Prior-art and novelty boundary

Bernstein and Lange already search legal addition chains by greedily maximizing
the number of distinct finite slopes after each multiplication for small prime
orders.  Their work also explains unavoidable anti-collisions and warns that a
good slope set without a cheap addition chain is insufficient.
[Bernstein--Lange, Sections 2 and 5](https://eprint.iacr.org/2012/294.pdf)

Chateauneuf, Ling, and Stinson established the slope-packing and slope-covering
connection to generic discrete-logarithm algorithms and performed computational
searches for small fields.
[Chateauneuf--Ling--Stinson](https://doi.org/10.1002/jcd.10033)

The approximate membership component is the classical Bloom-filter tradeoff:
space is reduced by permitting false positives while retaining no false
negatives under the data-structure assumptions.
[Bloom, 1970](https://doi.org/10.1145/362686.362692)

Consequently, coverage-guided chain search, slope coverage, legal addition
chains, and probabilistic membership are established ideas.  The bounded
combination below is a concrete research mechanism, not an established novel
algorithm.  Scoped searches did not locate this exact storage-tiered planner,
but no exhaustive novelty review was performed.

## [G2] Legal target-independent planner

The archive begins with the public coefficient points

    O=(0,0), A=(1,0), B=(0,1).

At step `t`, the planner deterministically proposes at most 32 unordered parent
pairs from the archive.  Four proposal lanes mix:

- the most recent record with a global archived record;
- a record from the most recent 64 with a global record;
- two global records; and
- two records from the most recent 64.

When the archive is very small, every unordered pair is proposed.  Candidate
coefficients are computed as `Z=P_i+P_j` without a group oracle call.  This is
legal because a generic algorithm knows the coefficient transcript and may
decide which archived operands to submit to the group oracle.  Only the selected
pair consumes one charged group addition.

For each candidate `Z`, the planner forms every useful finite slope from `Z` to
the existing archive.  It deduplicates these slopes locally and queries a Bloom
sketch containing every previously covered finite slope.  The score is the
number of Bloom negatives.  Deterministic parent indices break ties.  After the
selected group addition, all new exact slopes are inserted into the sketch and
the exact coefficient record is appended.

The primary experiment is non-inverting: subtraction candidates are disabled,
so no free negation is hidden.  All decisions use coefficients, the public seed,
step number, and archive length.  They do not hash or inspect a target group
element, making the coefficient chain target-independent.

## [G3] Bloom semantics and exact verification

The finite sketch has 1,048,576 bits and seven deterministic hash positions.
It never deletes.  Because every new point-to-history slope is inserted, a
correct sketch has no false negatives for the exact covered-slope set.  False
positives can only underestimate a candidate's gain.

The checker simultaneously maintains an exact slope set as an observer.  That
set validates false-negative absence, records exact marginal gains, and supplies
the scientific metric; it is not available to the Bloom planner.  A second
planner uses the exact set for scoring but receives exactly the same proposed
parent pairs.  In all eight tested seeds, Bloom and exact scoring selected the
same chain.  Bloom queries produced 176 false positives across the eight runs,
but none changed a winning proposal.

This agreement is finite evidence about the declared sketch load.  It does not
prove agreement at another budget, field, proposal width, or filter occupancy.

## [G4] Finite public-synthetic results

The field is `F_65537`; every method receives exactly 427 charged additions and
retains all 430 comparison events.  Useful slopes exclude the vertical
projective direction.  For a target-independent coefficient chain with `d`
finite slopes, exact generic collision success before a final guess is `d/p`.

Eight public seeds are used for each guided or random method:

| Method | Final finite slopes, range | Mean slopes | Mean exact success |
| --- | ---: | ---: | ---: |
| Bloom-guided | 51,557–51,700 | **51,634** | **0.787860** |
| Exact-score planner | 51,557–51,700 | **51,634** | **0.787860** |
| Best charged three-line endpoint | 51,210 | 51,210 | 0.781391 |
| Random accumulator | 48,992–49,152 | 49,051 | 0.748447 |

Every guided seed beats the strongest charged three-line endpoint and every
matched random accumulator in final slope count.

The complete prefix is also better.  For coverage `d_t` after `t` additions,
the checker computes the exact uniform-target censored mean

    E[min(T,427)]
      = (sum_t t(d_t-d_(t-1)) + 427(p-d_427))/p.

| Method | Censored mean additions |
| --- | ---: |
| Bloom-guided range | **284.290–284.567** |
| Bloom-guided mean | **284.415** |
| Random-accumulator mean | 292.429 |
| Best charged three-line prefix | 311.990 |

The finite result therefore survives both endpoint and prefix comparisons in
this model.  It is a result about generic slope success per charged group
addition while planning work remains separately counted.

## [G5] Finite planning and locality counters

One representative 427-addition run performs:

| Counter | Value |
| --- | ---: |
| Candidate evaluations | 13,584 |
| Parent-record queries | 27,168 |
| Bloom membership queries | about 2.942 million |
| Bloom hash probes | about 12.45 million |
| Sequential archive records scanned | 92,232 |
| Logical sequential scan bytes at 128 bytes/record | 11,805,696 |
| Parent page-cache misses | 0 |

The zero parent miss count is not a scalability result: 430 records fit inside
the modeled eight pages of 64 records.  Parent requests are known before
scoring and are grouped by page.  The complete coefficient archive is then
scanned once per step, allowing all 32 candidates to be evaluated in one
sequential pass rather than 32 independent scans.

The one-pass batching does not remove the asymptotic cost.  With `N` retained
records and 128-byte slots, the logical archive-scan traffic is

    64 N (N+5) bytes

for `N-3` planned additions.  Sketch queries and candidate arithmetic are also
quadratic in `N` for fixed proposal width.

## [G6] 1 TB RAM / 100 TB disk envelope

Starting from the external-directory Mode B envelope, the planner reserves an
additional 128,000,000,000 bytes of RAM for the all-slope Bloom sketch.  Fixed
RAM becomes 497,310,378,496 bytes.  With the same two bytes of per-record RAM
filter allocation and the same conservative disk peak equation, the byte
inequalities permit:

| Bound | Records |
| --- | ---: |
| Maximum by RAM | 251,344,810,752 |
| Maximum by peak disk | 308,441,358,024 |
| Capacity-envelope maximum | 251,344,810,752 |

This record capacity is misleading for the planner.  Under independent uniform
hash assumptions, a seven-hash 128-GB Bloom sketch reaches its nominal 1%
false-positive design load after approximately 106,745,005,079 distinct slope
insertions.  If pair slopes are nearly unique, that load is reached after only
about 462,050 records.

By 462,050 records, the cumulative sequential-scan formula already gives
13,663,520,816,000 bytes.  At the byte-envelope record maximum it gives
4,043,149,689,165,786,806,992,896 bytes of cumulative logical scans.  These are
traffic totals over time, not resident storage, but they show why satisfying
the RAM and peak-disk inequalities does not establish a usable planner.

Both guided and comparator methods still receive the same exact output archive,
collision filter/directory, external merge headroom, and exact certificate
verification.  The guided method's additional costs are candidate-parent reads,
full archive scans, slope hashing, and sketch updates.

## [G7] Verification and decision

The checker verifies:

- every selected output is the sum of two already archived parents;
- every method uses 427 charged additions and 430 retained records;
- the Bloom sketch has zero false negatives in all finite runs;
- exact slope sets equal independently recomputed final profiles;
- Bloom and exact scoring receive identical deterministic candidate schedules;
- the charged three-line controls are exhaustive over every feasible integer
  parameter under the budget; and
- all event and selected-parent transcripts are hash-bound in the JSON.

Artifacts:

- `check_coverage_guided_chain.py`
- `coverage_guided_chain_checks.json`

The mechanism is a genuine finite candidate: it improves exact generic endpoint
success and censored prefix mean over the charged three-line and random controls.
It does not yet support an end-to-end speedup claim.  The planner performs
millions of cheap operations and quadratic archive traffic to save a few hundred
group additions' worth of coverage.  A scalable successor would need a sketch
whose updates and gain estimates avoid enumerating all new point-to-history
slopes, or a proof that a small sampled anchor set preserves the ordering of
candidate gains.  Neither is supplied here.  Novelty, hardware benefit, and
cryptographic-scale feasibility remain unestablished.
