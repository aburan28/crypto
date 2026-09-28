# Completion audit: theoretical tiered rho design

Date: 2026-09-21

Audited objective: design a theoretical hybrid Pollard-rho-style approach using
1 TB RAM and 100 TB disk tiering that induces a speedup.

## Requirement 1: concrete new mechanism

**Evidence:** `COMPILED_RESERVOIR_GUIDED_WALK.md` specifies both compiler and
online executor.  The online recurrence is

    X_(t+1) = X_t + X_(t+1-d_t),   1 <= d_t <= 64.

There is one active register, one hot archived addend, one charged group
addition, one exact append, and one collision probe per step.  Six-bit offsets
are compiled from target-independent coefficient geometry using a bottom-`k`
reservoir and Bloom sampled-slope sketch.

**Judgment:** complete as a theoretical mechanism.  It is an archive-dependent,
time-inhomogeneous rho-style walk rather than a classical fixed `F:G->G` walk.

## Requirement 2: exact correctness

**Evidence:** every record retains exact coefficients and canonical group
encoding.  An equality gives

    x = -(a_i-a_j)(b_i-b_j)^(-1) mod n

when the denominator is nonzero, followed by exact `xP=Q` verification.  The
certification checker exhaustively replays all 65,537 public target logs.

**Judgment:** complete for the declared finite prime-order model.

## Requirement 3: storage tiering and peak budgets

**Evidence:** RAM holds the current record, 64 hot records, collision filter,
directory/cache, and compiler sketch.  NVMe holds 128-byte exact records,
external directory runs, temporary rewrite copies, and checkpoints.  Exact
integer capacity checks include maximum-plus-one failures.

Compiler ceilings under decimal limits:

- fixed reservoir, 800 GB RAM sketch: order scale floor 67 bits;
- fixed reservoir, phase-reused 100 TB external sketch: floor 75 bits;
- scaled reservoir, 800 GB RAM sketch: floor 41 bits; and
- scaled reservoir, external sketch: floor 48 bits.

**Judgment:** complete for the declared formats and order ceilings.  Capacity is
not a throughput or 256-bit feasibility result.

## Requirement 4: induced speedup

**Evidence:** the pinned 321-byte schedule with SHA-256
`3a07aeaf38483b9facabe21ef8415be95617534af03a6328d05d846a2efbccd3`
is replayed for every target in `F_65537`.  Exact capped means are

    guided = 19027631/65537,
    control = 19114086/65537.

The exact improvement is

    86455/65537 = 1.3189... group additions,

a 0.4523% online reduction after compilation.  Archive appends and directory
probes are one per executed addition, so their expected online counts decrease
with the stopping count as well.

**Judgment:** complete as an exact finite theoretical online speedup.  It is not
an asymptotic or measured wall-clock speedup.  Compilation must be charged once
or amortized over targets.

## Requirement 5: probabilistic data structure

**Evidence:** exact and one-megabit/seven-hash Bloom compilation are compared
under identical public seeds.  The retained walk audit records four false-
positive queries, zero false negatives, zero changed choices, and 64 matching
exact/Bloom transcripts.  Bloom results only choose legal steps; exact collision
verification never relies on approximate membership.

**Judgment:** complete for the declared finite filter load.

## Requirement 6: strong controls and falsifiers

**Evidence:** controls include matched unguided hot-addend walks, strong
pseudorandom-table HIAA controls, charged three-line constructions, legal random
accumulators, exact/Bloom pairs, fixed/scaled reservoirs, duty cycling, static
and epochal probes, and three-order schedule compilation.  Negative branches
and rank-deficient outcomes remain in the evidence package.

**Judgment:** complete for a theoretical candidate.  Scaling evidence shows the
fixed-reservoir gain shrinking, while scaled compilation becomes linear in the
group order; these are retained limitations rather than hidden failures.

## Requirement 7: novelty boundary

**Evidence:** primary-source review covers Bernstein--Lange greedy slope chains,
Teske/additive walks, Chateauneuf--Ling--Stinson slope packing, Bloom membership,
bottom-`k` sampling, delayed duplicate detection, HIAA-like memory dependence,
and multi-user collision methods.  Scoped searches did not locate the exact
combination of a bottom-`k` sampled-slope compiler, Bloom guidance, bounded hot
addends, six-bit target-independent replay, and tiered exact archive.

**Judgment:** complete as a distinct proposed combination after a scoped review.
Independent publication, patent, or exhaustive-literature novelty is not
claimed.

## Requirement 8: claim discipline

The following are proved or exactly checked:

- finite coefficient and collision correctness;
- the pinned all-target stopping law;
- Bloom/exact transcript equality at the tested load;
- exact byte ceilings; and
- matched finite control improvements.

The following remain unproved and are not part of the completion claim:

- a nonvanishing asymptotic constant improvement;
- positive wall-clock amortization on real hardware;
- external-sketch random-probe service;
- 256-bit feasibility;
- fixed-state rho coalescence; and
- exhaustive independent originality.

## Completion judgment

The requested theoretical design is complete: it is concrete, storage-tiered,
correct in its finite model, and has an exact finite online speedup after
target-independent compilation.  The broader scientific promotion gates remain
open and are explicitly excluded from the achieved claim.

Authoritative artifacts:

- `COMPILED_RESERVOIR_GUIDED_WALK.md`
- `check_certified_compiled_walk.py`
- `certified_compiled_walk.json`
- `certified_compiled_walk_schedule.bin`
- `check_walk_reservoir_guided.py`
- `walk_reservoir_guided_checks.json`
- `check_walk_reservoir_bloom.py`
- `walk_reservoir_bloom_checks.json`
- `check_compiled_walk_scaling.py`
- `compiled_walk_scaling_checks.json`
- `check_compiled_walk_capacity.py`
- `compiled_walk_capacity.json`
