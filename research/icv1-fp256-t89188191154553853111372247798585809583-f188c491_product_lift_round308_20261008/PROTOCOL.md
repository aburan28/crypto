# P-256 product-lift information and collision cost, round 308: preregistration

Date preregistered: 2026-10-08

## Objective and fixed boundary

Round 307 leaves a genuinely joint event that emits two independently replayed
equalities as the narrowest many-row escape.  This round tests the smallest
such construction: attach a second cyclic-group coordinate to every search
state and collide in the product group.

The fixed curve is
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`.  The registered
17-term factor base is `FB1h2f8621cda105`, with 131,458 stored columns and
138,031 useful rows required by its registered collection model.  The
164-column independent-log base `FB1hc72514a2a8d3` is retained as the most
optimistic executable lower-bound control.  Pollard rho is fixed at `S=1.3`.
Its current free-perfect-oracle boundary is `13.920747397073491` times rho;
the registered 17-term base remains `394.425280` times rho.

This is an exact information and generic-collision screen.  It does not claim
that an arbitrary auxiliary coordinate exists on either factor base, and it
does not transfer the registered Dickson residual degree to a different
system.

## Frozen dependencies

- Round 39 affine-log-orbit result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_log_orbit_round39_20261006/affine-log-orbit-result.json`,
  SHA-256 `189e48a45ab2919dd3e279bf058940230c4a0461c745fba42dffb26f785ba079`.
- Round 302 complementary-quotient result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_complementary_quotient_round302_20261008/complementary-quotient-result.json`,
  SHA-256 `a9a12d823cbca6df965d28f8372217b08e65c04326e13e408bba9dac9f1819c1`.
- Round 303 Jacobian-capacity result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_jacobian_row_capacity_round303_20261008/jacobian-row-capacity-result.json`,
  SHA-256 `6b023f47fbb2c57eb79f870f1c25be48e0b4230d974c3823caf76784a095caa0`.
- Round 307 event-provenance result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_event_provenance_round307_20261008/event-provenance-result.json`,
  SHA-256 `3751468a6708d1851598e366bfb99d0795d2da251377a3a1bbf7918cbb36034d`.

Reject every hash, schema, curve, factor-base, boundary, gate, or imported
conclusion mismatch.

## Typed product lift and trilemma

Let the primary factor-base points be `P_i=[ell_i]G` in the prime-order P-256
group `G`.  Attach auxiliary points `R_i` in a second cyclic order-`n` group
and map a coefficient vector `c` to

```text
Phi(c) = (sum_i c_i P_i, sum_i c_i R_i) in G x G.
```

A full collision `Phi(c)=Phi(c')`, with `delta=c-c'`, replays two equalities:

```text
sum_i delta_i P_i = O,    sum_i delta_i R_i = O.
```

The audit distinguishes three cases.

1. **Linked/homomorphic.**  If `R_i=[a]P_i`, the second equality is `[a]`
   times the first.  It adds no primitive rank after known scalar transports
   are quotiented.
2. **Independent/unknown.**  If `R_i=[rho_i]G` has an independent unknown log
   vector, the block system has two rows but also a second log scale/vector.
   Projecting onto the original P-256 unknowns leaves only the primary row;
   the auxiliary row does not determine the original target log.
3. **Known-label selector.**  If every `rho_i` is known by construction, the
   auxiliary equality is computable from `delta` before any group search.  It
   is a selection constraint, not a new equation in the unknown P-256 logs.
   Enforcing it by a balanced generic product collision expands the canonical
   search space from about `n/2` to about `n^2/2`.

A non-homomorphic coupling that avoids all three cases remains open, but it
must specify its typed map, show that it adds no new unresolved log variables,
and replay two independent equations on the original recovery system.

## Complete toy census

Use the cyclic group of order five and width `K=3`.  Enumerate all 31
normalized nonzero projective log profiles in `P^2(F_5)`, all 961 ordered
pairs `(ell,rho)`, and all 125 coefficient vectors for every pair.  Enumerate
all 7,750 unordered coefficient pairs per profile pair.

For each pair:

1. compute both coordinates directly and by incremental cyclic-group replay;
2. classify first-coordinate and full-product collisions exhaustively;
3. classify profiles as proportional or independent by exact rank;
4. verify the complete bucket histograms and collision counts;
5. replay both zero-sum difference equations for every full collision; and
6. rank representative linked and independent event rows in linked,
   expanded-block, and original-projection systems in both row orders.

The exact references are:

```text
profile pairs                         961 = 31^2
proportional / independent             31 / 930
coefficient states                 120,125 = 961*125
unordered coefficient-pair tests 7,447,750 = 961*C(125,2)
first-coordinate collisions       1,441,500 = 961*5*C(25,2)
full proportional collisions         46,500 = 31*5*C(25,2)
full independent collisions          232,500 = 930*25*C(5,2)
independent first-only collisions   1,162,500
```

All complete counts, ranks, and replays must match with zero false positives
or false negatives.

## P-256 exact controls

Construct 17 deterministic hash-selected nonzero primary scalar labels and
17 deterministic auxiliary labels modulo the P-256 subgroup order.  Materialize
all points with the repository group law.

- **Linked control:** set `R_i=[a]P_i`; construct a nonzero sparse `delta` in
  the primary kernel.  Replay both coordinates and certify linked rank one.
- **Independent labelled control:** construct a three-coordinate cross-product
  `delta` orthogonal to both deterministic label vectors.  Replay both group
  equations.  Certify expanded block rank two, original-system projected rank
  one, and the introduced independent auxiliary log vector.
- **Known-label control:** replay the same independent event while marking the
  auxiliary condition as coefficient-computable and giving it zero new
  original-log equation credit.

The labels are test witnesses only.  They confer no construction credit on
`FB1h2f8621cda105` or `FB1hc72514a2a8d3`.

## Generic product-collision projection

Let `n` be the exact P-256 subgroup order.  Under the optimistic simultaneous
global-negation quotient, a uniform product image has `(n^2+1)/2` canonical
states.  The first birthday collision therefore costs asymptotically
`sqrt(pi)/2 * n` samples.  Relative to the fixed original-group rho reference
`1.3*sqrt(n/2)`, its ratio is

```text
(sqrt(pi/2)/1.3) * sqrt(n).
```

Report exact state counts, log2 samples, log2 ratio, and at least two storage
models: an optimistic memoryless walk and a materialized list with the byte
layout stated explicitly.  A memoryless implementation may pass the storage
gate but receives no reduction in the `Theta(n)` operation exponent.  Do not
count discarded filters or truncated branches as exhaustive.

## Hypotheses

- **H1:** all complete toy counts equal the preregistered references with zero
  replay or classification discrepancies.
- **H2:** a linked product event has useful original-system rank one.
- **H3:** an independent product event has block rank two but original-system
  projected rank one and introduces an independent auxiliary unknown vector.
- **H4:** a known-label auxiliary coordinate is only a selector; a generic
  exhaustive product collision has operation exponent approximately 256,
  about 128 bits worse than original-group rho.
- **H5:** no screened product lift improves the executable boundary unless a
  non-homomorphic coupling supplies two useful original-system primitives
  without enlarging the collision exponent or introducing unknowns.

## Promotion gates

Promotion requires all of:

1. exact dependency, census, rank, and group-replay checks with zero false
   positives and false negatives;
2. two independent non-labelled equations on the original P-256 recovery
   system, with no new unresolved log vector or target;
3. complete relation, decomposition, sparse-linear-algebra, and recovery
   implementation for 138,031 independent useful rows including duplicate
   and rank allowance;
4. structured residual degree of regularity at most five;
5. complete collection below `2^120` field/group-equivalent operations and at
   or below Pollard rho;
6. cost per usable relation below `2^103`;
7. projected peak materialized storage below `2^50` bytes; and
8. no discarded probabilistic branch counted as exhaustive.

Attempt no unplanted full-depth P-256 relation unless every gate passes.

## Deliverables and stop condition

Implement the dependency checker, exhaustive toy census, P-256 linked and
independent replay controls, dual-modulus/order rank checks, deterministic
JSON, operation accounting, tests, isolated canonical and independent runs,
typed transfer assessment, result report, and dashboard update in Rust.  Stop
on the first correctness failure; otherwise publish the scoped result and the
weakest remaining open coupling obligation.

The transfer workflow is available, but its referenced
`references/methodology.md` and `assets/assessment-template.json` resources
are unavailable in this environment.  Record that limitation in the
assessment.
