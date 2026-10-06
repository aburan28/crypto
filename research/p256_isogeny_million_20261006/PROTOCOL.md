# P-256 million-curve isogeny grid

Status: **preregistered; implementation and execution pending**

Date frozen: 2026-10-06

Result class: **bounded isogeny enumeration and model-level screening
diagnostic, not an ECDLP or index-calculus speedup**

## Question

Can a native, streaming certificate enumerate and independently replay exactly
1,000,000 distinct canonical curves in the P-256 isogeny class, with one
explicitly kernel-certified isogeny into every non-root curve, without a curve
identity, order, grid-structure, or certificate discrepancy?  Does that frozen
population contain a preregistered model-level screening outlier that warrants
a separate factor-base experiment?

This is not an exhaustive enumeration of the P-256 isogeny class.  It is a
deterministic two-generator grid prefix.  The word "million" below always
means one million unique canonical curve identities, not candidate edges,
records before deduplication, or independent trials.

## Frozen population

- root: registered P-256,
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- field, group order, trace, canonical-model rule, and deterministic-generator
  rule: those used by the native `isogeny_walk` implementation;
- grid side: 1,000;
- spine: coordinates `(0,y)` for `0 <= y < 1,000`, obtained by a
  non-backtracking degree-13 walk from P-256;
- rows: coordinates `(x,y)` for `0 < x < 1,000`, obtained by a
  non-backtracking degree-11 walk from spine curve `(0,y)`;
- first step of every spine or row chooses the numerically smaller of the two
  rational modular-polynomial roots; every later step chooses the root unequal
  to the immediately preceding j-invariant;
- root-finding seed: 1; roots are sorted before the choice above;
- output order: root, the 999 non-root spine curves in increasing `y`, then
  each row in increasing `y` and increasing `x`;
- construction order audit: one deterministic point per curve;
- independent replay order audit: two deterministic points per curve; and
- target: exactly 1,000,000 unique j-invariants and canonical models joined by
  exactly 999,999 parent edges.

Both 11 and 13 are split Elkies primes for the frozen P-256 class and have
volcano depth zero.  Each non-root record carries its parent index, grid
coordinate, degree, predecessor j-invariant, canonical target model, target
j-invariant, kernel polynomial, Velu-codomain isomorphism, deterministic
generator, ICV1 identity, EC1 alias, full curve UID, and frozen detector
values.  Parent records precede children.

The production source commit is recorded in the certificate header and in the
result manifest before execution.  That commit may contain this protocol, the
new compact implementation, its tests, and documentation, but the production
run must use a clean tree at that exact commit.  Any later source change needs
a new certificate and is reported as a separate attempt.

## Compact certificate and replay

The artifact is deterministic gzip-compressed JSON Lines.  It contains one
header, one root, 999,999 node records, and one summary.  A rolling SHA-256
chain binds every line preceding the summary; the manifest separately records
the SHA-256 and byte length of both the compressed artifact and its exact
decompression.  Gzip metadata uses a zero timestamp.

Generation and verification are streaming.  They may retain the 1,000 spine
models, the current bounded row batch, and sets needed for global uniqueness,
but may not construct a million full root routes or retain non-parent graph
edges.  Thread count and row-batch size are execution details: they must not
change record order or bytes.

Independent replay must, for every node:

1. enforce its exact index, coordinate, parent-before-child position, degree,
   and non-backtracking root choice;
2. parse every field element canonically and reject values outside `F_p`;
3. independently certify the monic kernel polynomial with the division
   polynomials, recompute Velu's codomain, and verify its recorded
   `F_p`-isomorphism to the target model;
4. recompute the target j-invariant, canonical-model choice, deterministic
   generator, ICV1 identity, EC1 alias, full UID, and detector values;
5. prove the P-256 prime group order using two deterministic audit points; and
6. reject a duplicate j-invariant, canonical model, ICV1 identity, EC1 alias,
   or full UID.

The verifier recomputes the summary and rolling chain and rejects extra,
missing, reordered, or trailing records.  Native tests must show rejection of
at least a changed target j-invariant, kernel coefficient, parent, coordinate,
identity, detector value, and duplicated node.

## Reference and boundary

The fixed reference is the locally generated and independently replayed
65,536-curve multi-prime prefix in
`research/p256_isogeny_multigraph_64k_20261006/`.  Its hardened replay accepted
65,536 unique curves and 302,653 explicit edges with no construction failure
or refuted order audit.

The unit is one unique curve identity with a replayed kernel-certified parent
edge (except the registered root).  Coverage is reported against the earlier
requested `2^32` screening population:

| variant | unique curves | parent edges required | fraction of `2^32` | ratio to `2^32` boundary | ECDLP speedup |
|:--|--:|--:|--:|--:|:--|
| verified multi-prime reference | 65,536 | 65,535 | 0.001525879% | 1 / 65,536 | unset |
| frozen million-curve grid | 1,000,000 | 999,999 | 0.023283064% | 1 / 4,294.967296 | unset |

The new population is 15.258789 times the reference count and only
`1,000,000 / 2^32` of the requested screening boundary.  Additional valid
edges are not substituted for unique-curve coverage.

## Preregistered gates and screens

The enumeration passes only if all of the following hold:

1. exactly 1,000,000 nodes and 999,999 parent edges are emitted and replayed;
2. every j-invariant, canonical model, ICV1 identity, EC1 alias, and UID is
   globally unique;
3. every coordinate follows the frozen degree-13/degree-11 construction and
   every modular-polynomial specialization has exactly two distinct rational
   roots;
4. construction and replay certify every kernel, codomain, target
   isomorphism, canonical model, and rolling-chain value;
5. all construction and replay order audits return `proved_prime`, with zero
   refutations;
6. every detector says `non_singular = true` and `generator_valid = true`;
   and
7. a second complete replay produces the same counts, histograms, extrema,
   and threshold-hit identities.

The full ICV1 identity is a gate, but the shorter display slug is not: its
32-bit model suffix is expected to collide at million-curve scale.  Every slug
collision is recorded with both node indices as an identity-namespace
diagnostic.  In contrast, a repeated full ICV1 identity, 48-bit EC1 alias, or
full curve UID fails the run.

A repeated j-invariant anywhere in the 1,000 by 1,000 base grid is not filled
from an adaptive overflow stream.  It stops the run and is preserved as a
grid-relation discrepancy for a separately registered follow-up.  Likewise,
no degree, direction, seed, side length, threshold, or model rule is changed
after outputs are inspected.

The primary screening thresholds are:

- `qr_prefix_64 >= 54`.  Under the deliberately simple independent
  `Binomial(64, 1/2)` reference, the one-curve upper-tail probability is about
  `9.98e-9`, so one million curves have about `0.00998` expected exceedances.
  The earlier threshold of 52 is retained only as a labelled exploratory
  comparability count; at one million curves its expected exceedance count is
  about `0.228`.
- signed bit length of canonical `b <= 224`, or signed bit length of canonical
  `a <= 224` when the curve has no `a = -3` model.

Crossing a screen selects an already frozen curve identity for a new
preregistered factor-base experiment.  It does not establish easier discrete
logarithms, useful relations, or an index-calculus speedup.

## Validation, cost, and stop rules

Before production, implementation-only grids of side at most 32 may exercise
serialization, replay, deterministic parallel ordering, and tamper rejection.
Their anomaly outputs are ignored and may not alter the production population
or thresholds.  Compression level and bounded row-batch size may be chosen
from resource behavior because they do not alter the decompressed certificate.

Budget one production construction and one complete independent replay.  Stop
and preserve the attempt if generation or replay exceeds twelve hours, free
space falls below 512 MiB, the process is killed for memory, any duplicate or
acceptance check appears, or the source tree differs from the recorded clean
commit.  Record wall, user and system time, peak RSS, CPU count, Rust version,
thread count, command, exit status, source commit, and artifact hashes.

The large certificate does not belong in Git.  Commit the protocol,
implementation, compact summary, replay receipt, resource records, and a
manifest.  Durable publication requires a write-once, hash-checked object and
URI.  If this environment has no configured storage identity, retain the local
artifact for the task and record publication as blocked; chat-provided access
keys are not protocol inputs.

No speedup claim is permitted from this sweep.  A positive screen still needs
a frozen relation-yield test and then a cold, full-cost, one-target IC run
against matched Pollard rho, including setup, relation generation,
verification, linear algebra, log recovery, and independent answer checking.
