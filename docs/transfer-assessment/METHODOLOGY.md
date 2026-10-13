# Diagnostic methodology for ECDLP transfers

## Purpose and verdicts

Use this methodology to assess whether moving an elliptic-curve discrete-log
problem through a cover, restriction of scalars, isogeny, Jacobian, abelian
variety, quotient surface, or lattice representation creates a useful
reduction. A new equation or representation is not itself a reduction.

Every obligation has exactly one state: `unknown`, `supported`, `refuted`, or
`not_applicable`. Missing evidence is `unknown`, never failure. Verdicts are
cumulative:

1. **valid transfer**: map correctness, subgroup preservation, and scalar
   recovery are supported;
2. **algorithmically relevant**: a named destination solver applies to the
   transported subgroup;
3. **modeled advantage**: complete modeled cost beats the strongest applicable
   source baseline under stated assumptions;
4. **measured advantage**: paired end-to-end measurements support that model
   in the tested scope;
5. **novel reduction**: literature and equivalence checks support novelty
   inside a stated boundary.

A refuted hard gate rejects the dependent verdict. Do not collapse the gates
into a score. An abelian surface has dimension two; a genus-`g > 2` Jacobian is
a `g`-dimensional abelian variety. A Kummer surface or ordinary Euclidean
lattice is not itself an abelian group.

## 1. Freeze the exact problem

Record characteristic, exact field and basis, curve equation and model,
dimension or genus, field of definition, group order, prime subgroup order
`r`, cofactor, generator `P`, and target `Q=[x]P`. Bind the field, curve,
representation, candidate, and workload identities using the repository's
ICV1/EC1 rules. State classical or quantum scope, cold-start or amortized
accounting, hardware, memory, success probability, and target count.

Type every node and edge. For example,

```text
E(F_q^n) --rational-point correspondence--> Res(E)(F_q)
         --isogeny/correspondence?---------> J(C)(F_q)
         --named destination solver--------> x mod r.
```

The first edge changes representation but supplies no speedup. A dashed edge
remains an open obligation until its exact kind and certificate are known.

| Transformation | Required certificate |
|---|---|
| Coordinate isomorphism | Forward/inverse formulas, exceptional points, unchanged group |
| Isogeny | Homomorphism, kernel, rationality, degree, subgroup restriction |
| Field extension | Embedding and identification of the original subgroup |
| Restriction of scalars | Dimension, rational-point correspondence, arithmetic cost |
| Curve cover | Curve morphism and induced pullback/pushforward on divisor classes |
| Jacobian correspondence | Induced group maps, kernel, polarization, working-field descent |
| Quotient/Kummer | Lost sign or phase, compatible operations, scalar recovery |
| Lattice representation | Torsion location, membership, coordinates, reconstruction |

## 2. Hard gate: executable subgroup transfer

Require an explicitly computable homomorphism

\[
\iota:\langle P\rangle\longrightarrow A(\mathbb F_q)
\]

such that `iota(P)` and `iota(Q)` are computable without `x`, scalar
multiplication commutes with the map, `ker(iota) intersect <P>={0}`, and a
destination answer recovers `x mod r` and verifies `[x]P=Q` on the source.
For prime `r`, a proved homomorphism with nonzero `iota(P)` is injective.
Composite orders require every prime-power component and recombination.

A certificate `R o iota=[m]` with `gcd(m,r)=1` proves subgroup injectivity.
For a finite separable degree-`d` cover `pi:C->E`, the induced maps satisfy

\[
\pi_*\pi^*=[d].
\]

Thus `gcd(d,r)=1` is sufficient for injectivity of pullback on the order-`r`
subgroup. If `r` divides `d`, compute the actual kernel intersection. A curve
map `C->E` is not a pointwise group embedding `E->C`; the group is `J(C)`.
Random identities test code but do not replace the homomorphism proof.

## 3. Cheap algebraic screens

Run necessary filters before expensive construction or benchmarks.

| Screen | Required data | Refutation |
|---|---|---|
| Dimension | Every abelian/Jacobian dimension | Claimed full factor is larger than destination |
| Group order | `r` and destination `r`-primary part | `r` does not divide destination order |
| Frobenius/Weil | Exact polynomials over the working field | No compatible common isogeny factor |
| Field of definition | Galois action on formulas/kernel | Claimed base-field map exists only after unaccounted extension |
| Riemann--Hurwitz | Degree, different/ramification, genus | Data are inconsistent |
| Polarization/Jacobian locus | Polarization and curve realization | Unpolarized class is asserted to be a usable Jacobian |

For `Res_(F_q^n/F_q)(E)`, dimension is `n`. A compatible Weil polynomial only
supports an isogeny-class filter; it does not construct a map. For a separable
cover of an elliptic curve, `2g(C)-2=deg R`; use the different in wild or
inseparable cases.

## 4. Destination-solver gate

Name the exact algorithm and verify its hypotheses on the image subgroup:
object class, field regime, genus/dimension, model, group assumptions,
factor-base identity, smoothness model, relation validity, matrix rank modulo
`r`, linear algebra, individual-log descent, recovery, and source replay.
Dimension, genus, a smaller base field, or faster known-scalar arithmetic does
not imply a faster unknown-scalar solver.

If a cited solver costs `q^(alpha+o(1))`, define `kappa=log_q(r)` and use
`alpha < kappa/2` only as a first classical exponent screen. Recompute for the
actual subgroup and strongest baseline; equality leaves constants, memory,
setup, and success probability unresolved.

Record the named solver and its image-subgroup interface under
`solver_assessment`. A family label such as “index calculus” is insufficient:
bind the factor base or search space, relation/decomposition rule, linear
algebra, individual-log procedure, implementation revision, recovery
interface, and working field. Use the three-arm and nine-phase protocol in
[`MEASUREMENT.md`](MEASUREMENT.md) to distinguish a destination-only stage
improvement from a net transfer advantage.

## 5. Representation-specific rejection tests

- **Restriction of scalars:** require a second solver-useful edge. If
  `J(C) ~ E x B` and the image remains in `E`, appending `B` does not weaken
  that factor.
- **Abelian/Kummer surfaces:** name the Albanese, `Pic^0`, Jacobian, Prym, or
  abelian variety carrying the group. `A/{+1,-1}` loses sign and supplies only
  pseudomultiplication until collision and sign recovery are given.
- **Lattices:** `Z^m` is torsion-free, so every homomorphism `Z/rZ -> Z^m` is
  zero. A nontrivial proposal must use torsion in a finite quotient, torus,
  endomorphism-kernel lattice, or relation lattice. If explicit coordinates
  already reveal the log, audit circularity. GLV/GLS lattices accelerate
  known-scalar multiplication; that alone does not recover an unknown scalar.

## 6. Complete cost and controls

Charge cold-start cost

\[
C_{new}=C_{search}+C_{construction}+C_{map}+C_{preprocessing}
+C_{relations}+C_{linear\ algebra}+C_{individual\ log}
+C_{recovery}+C_{verification}.
\]

Include failed attempts, duplicates, timeouts, conversion, rank updates, data
movement, memory, hardware, parallelism, and success probability. Separately
report `C_amortized(N)=C_reusable/N+C_per-target` only with a reuse certificate
and exact `N`. Keep modeled and measured columns separate. Compare with every
applicable direct baseline, including Pohlig--Hellman, optimized parallel rho
with source automorphisms, pairing transfer where applicable, existing Weil
descent, and direct decomposition methods.

At minimum use: strongest direct source baseline; the same generic solver on
the transported subgroup; a coordinate-isomorphic control; a matched random
destination; transported versus destination-reselected data; and null or
known-false structural controls. Toy arithmetic validates formulas, not
scaling or a cryptographic advantage.

The primary comparison has three arms: the strongest direct source solver,
the same generic solver on the destination image, and the specialized
destination solver with the full route charged. The admissible net ratio is
`source_reference_total / destination_candidate_total`. Store cold
one-target evidence in `advantage_measurement`; keep stage-only, online, and
amortized ratios explicitly secondary. The validator checks phase totals,
unpriced work, ratio arithmetic, replay evidence, confidence bounds, and
scaling-size requirements before an advantage verdict can be supported.

## 7. Workflow and skill routing

Use `transfer` to type and audit any proposed change of algebraic object. Use
the repository's `ecdlp-surfaces` skill for its supported surface, binary,
cover, Jacobian, gluing, polarization, or Frobenius tools. Use
`curve-algebra` when CM, Frobenius, automorphisms, GLV/GLS, or endomorphism
actions determine the candidate or baseline. Once an executable subgroup map
and named destination solver exist, use `ecbench-measure` for the matched
cross-method experiment and `ecbench-bounds` to challenge an established
frontier. Apply `report-evidence` to every implementation or run and
`research-visuals` to every substantive search, run, or result PR.

The detailed trigger and handoff table is in
[`MEASUREMENT.md`](MEASUREMENT.md#which-repository-skill-to-use). Tool
selection does not relax a mathematical or evidence gate.

## 8. Novelty boundary

Fingerprint `(fields, Weil polynomials, dimensions, map type, degree, kernel,
polarization, target model, solver)`. Check twists, field automorphisms,
birational models, dual/conjugate isogenies, factor permutations, reverse
paths, and alternative equations for the same polarized variety. State only
`not found within <sources, families, parameters, date>` unless independent
review supports more.

## Required output

Return the typed graph and exact identities; obligation table with state,
evidence, and scope; preservation and recovery certificate; named solver
assessment; three-arm, nine-phase cost table; controls and exploration
boundary; narrowest verdict; limitations; and one checkable next obligation.
The native validator in
[`tools/transfer-assessment`](../../tools/transfer-assessment/README.md) checks
the machine-readable gates. It validates consistency, not the truth of a
citation and not a performance, novelty, or security claim.

A bounded unsuccessful search means `not found within this family and
boundary`, never universal nonexistence. This diagnostic profile makes no
official research-state transition.
