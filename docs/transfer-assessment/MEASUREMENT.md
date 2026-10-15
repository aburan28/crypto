# Measuring destination-solver and net transfer advantage

## The question this protocol answers

A transfer can be mathematically valid while making the discrete logarithm no
easier. This protocol separates two questions:

1. **Destination question:** does a named solver exploit the destination
   representation more effectively than a generic solver on the same image
   subgroup?
2. **Net question:** after every transport and recovery cost is charged, does
   the complete route beat the strongest direct solver on the source
   subgroup?

Only the second question establishes a transfer advantage. A faster relation
collector, smaller field, larger dimension, lower asymptotic exponent, or
faster destination group law is a stage diagnostic until the complete route
is priced and the recovered scalar verifies on the source.

This document complements [the transfer methodology](METHODOLOGY.md). First
establish the typed map and subgroup certificate there; then apply this
protocol. The assessment schema stores the result under `solver_assessment`
and `advantage_measurement`.

## Three-arm comparison

Every candidate is compared with three arms on the same order-`r` subgroup and
target workload.

| Arm | Purpose | Required scope |
|---|---|---|
| `source_reference` | Strongest applicable direct source attack | Exact source curve, subgroup, target, automorphisms, hardware and resource envelope |
| `destination_generic_control` | Same class of generic attack after transport | Exact image subgroup; exposes representation and arithmetic overhead without crediting a specialized solver |
| `destination_candidate` | Named destination solver plus the complete transfer and recovery pipeline | Construction, maps, solver setup, solve, recovery and source verification |

The destination control prevents a coordinate change or appended factor from
being mistaken for an algorithmic improvement. It is still only a control:
the primary ratio is

\[
R_{net}=\frac{C_{source\ reference}}{C_{destination\ candidate}}.
\]

`R_net > 1` means lower cost in the declared unit. Also report the equivalent
work-factor difference

\[
\Delta_{bits}=\log_2 C_{source\ reference}
              -\log_2 C_{destination\ candidate}.
\]

Positive `Delta_bits` favors the candidate. A destination-only ratio may be
reported as a stage diagnostic, but it must not populate `R_net` unless the
candidate total includes the complete source-to-destination-to-source route.

For repository ECDLP experiments, use the normalized operation count
`S = total_operations / sqrt(r)` and show its ratio to the generic floor and
the matched reference. Foreign operation types need measured conversion
weights or must remain separate axes. A missing price is not zero.

## Name the actual destination solver

Record a solver identity, implementation and revision. “Index calculus,”
“Weil descent,” “lattice reduction,” or “higher genus” is a family, not a
solver specification. The specification must say what the implementation
consumes and produces and which hypotheses are used.

For an index-calculus candidate, record at least:

- exact destination group and image subgroup;
- field, model, genus or dimension, and infinity convention;
- factor-base identity and membership rule;
- relation-generation algorithm and its success model;
- relation verification, duplicate filtering and rank policy modulo `r`;
- sparse or dense linear algebra and all conversions;
- target-dependent individual-log or descent procedure;
- how the solver output is projected or normalized to `x mod r`;
- memory, parallelism, success probability and failed attempts; and
- the implementation revision and reproducible invocation.

The relation database must determine logarithms in the transported
order-`r` component. Relations that live only in a complementary factor do
not solve the source DLP. If `J(C)` is isogenous to `E x B`, explicitly show
how relations, linear algebra and individual logarithms reach the `E` image;
the existence of the factorization alone is not an advantage. A Kummer solver
must also account for sign/lift recovery.

## Fast triage before implementation

Use the following screens to decide whether a candidate deserves a measured
trial. Passing a screen is not a promotion.

1. **Applicability:** the published or implemented solver applies to the exact
   destination class, field regime, genus/dimension and subgroup.
2. **Image access:** the factor base and relation map can reach the transported
   order-`r` component; the image is not merely carried as an opaque elliptic
   factor.
3. **Exponent screen:** express both algorithms in one size variable. If a
   destination estimate is `q^(alpha+o(1))` and `r=q^kappa`, compare `alpha`
   with `kappa/2`, then include constants, memory and setup. This is only a
   model screen.
4. **Recovery:** individual logs, projection/norm and any modular
   normalization recover `x mod r` without already knowing `x`.
5. **No cheaper equivalent route:** compare twists, isomorphic models,
   isogenous representatives and direct source methods under the same
   resource envelope.

Typical interpretations are:

| Destination | What would make it solver-relevant | Common false positive |
|---|---|---|
| Coordinate image or isomorphic curve | A solver uses structure absent from the source description and all conversion is priced | Calling cheaper encoding arithmetic a DLP reduction |
| `E x B` or `E x E'` | Relations and recovery act on the image component with a lower complete cost | Appending an unrelated factor and quoting the ambient dimension |
| Genus-2 or genus-3 Jacobian | Executable subgroup map plus a compatible Jacobian solver, relation matrix modulo `r`, individual log and norm/projection | Treating a cover equation or Frobenius factor as a solver |
| Weil restriction | A solver exploits the descended equations after field and dimension growth are charged | Counting the smaller base field but not the added variables |
| Kummer surface | A compatible pseudogroup solver with lift/sign recovery | Assuming an ordinary group law on the quotient |
| Lattice representation | A typed finite quotient or relation lattice whose solution reconstructs the scalar non-circularly | Mapping finite torsion into ordinary `Z^m` |

## Complete accounting

Use these nine exclusive phase names in the machine record:

1. `route_search` — instance-specific discovery or reconstruction of the
   usable route; zero only when a cited reusable route is already supplied;
2. `construction` — curve/variety, map, kernel and solver-data construction;
3. `forward_map` — transport of generator and target, including conversions;
4. `solver_setup` — factor base, tables, precomputation and reusable relations;
5. `relation_or_walk` — relation collection or the destination search;
6. `linear_algebra` — rank work and solution modulo `r`;
7. `individual_log_or_descent` — target-dependent descent or decomposition;
8. `source_recovery` — norm/projection, ambiguity resolution and scalar
   normalization; and
9. `verification` — destination checks and independent `[x]P=Q` replay.

The candidate total is the sum of all nine phases. A phase may be zero only
when its role is genuinely absent and the evidence explains why. Preserve
failed attempts, duplicates, timeouts and out-of-memory outcomes. Track peak
memory and unpriced counters beside operations; do not silently exchange
unbounded memory for a lower operation count.

Report two accounting views:

- **cold one-target:** every reusable and target-dependent phase charged to
  one previously unseen target; this is the primary comparison;
- **amortized `N` targets:** `C_reusable/N + C_per_target`, only after stating
  exactly which bytes and computations are reusable and fixing `N`.

An amortized result cannot promote the cold verdict. It is a separate
workload. Public route-reconstruction cost, designer setup with a retained
witness, per-target transport and destination solving are also separate
threat-model quantities.

## Modeled advantage gate

The machine validator permits `modeled_advantage=supported` only when:

- the transfer and destination-solver verdicts are already supported;
- the primary scope is `cold_end_to_end` and the basis is `modeled` or
  `measured`;
- all three arms use the same nonempty unit and exact workload;
- every candidate phase is present, finite and included in the total;
- no work is listed as unpriced;
- the recorded total equals the phase sum;
- `R_net` equals `source_total/candidate_total` within numerical tolerance and
  is greater than one; and
- evidence is bound to the assessment's declared sources.

This is an internally consistent model, not a measurement or security claim.
State uncertainty ranges and sensitivity to smoothness, rank, memory and
success assumptions in the accompanying report.

## Measured advantage gate

Run source reference, destination generic control and complete candidate on
identical targets and seeds. A measured promotion additionally needs:

- native `ecbench` records for cross-method work;
- correct answers for every admitted pair and independent source replay;
- an audit receipt reproducing the counted work;
- a ratio-of-totals 95% interval whose lower endpoint is greater than one;
- at least five matched verified pairs for a runtime statement;
- A/A timing controls and the required isolation for any wall-clock claim;
- at least four sizes, with separate per-size rows, before claiming scaling;
- cold totals as the headline, with online or amortized views secondary; and
- the same nonempty resource envelope for all compared arms.

Wall time is a practicality axis. Counted operations are primary. A measured
stage improvement stays a stage result even with a narrow confidence interval.
Use `ecbench-bounds` only after a replayable end-to-end session exists; it
classifies an admitted change as exponent, constant, primitive, trade, match
or regression relative to the established frontier.

## Experimental ladder

Advance one rung at a time and preserve failures.

| Rung | Experiment | What it can establish |
|---|---|---|
| A0 | Arithmetic and equation fixtures | Formula implementation works on bounded examples |
| A1 | Map, kernel, image order and scalar-recovery replay | Valid transfer on the tested subgroup |
| A2 | Destination solver on tiny known-answer images, including false relations | Solver compatibility and complete pipeline correctness |
| A3 | Three-arm cold comparison on at least four toy sizes | Measured toy scaling and phase crossover, not cryptographic advantage |
| A4 | Faithful family gate such as the repository's pinned `m=83` Koblitz case where applicable | Evidence on that exact family and size; not an automatic extrapolation |
| A5 | Independent replay and frontier challenge | Reproducible scoped comparison against the incumbent |

The first practical higher-genus experiment should be small enough for both
generic and specialized solvers to finish, but it must use the same cover/map
family intended for later sizes. Changing genus, field family, subgroup
embedding or infinity model between rungs creates a new candidate.

## Which repository skill to use

| Situation | Skill | Stop condition / handoff |
|---|---|---|
| A proposed cover, isogeny, descent, product, Jacobian, quotient or lattice route | `transfer` | Produce typed maps, subgroup/recovery certificates, solver obligations and a cost boundary before any advantage claim |
| Quartic surfaces, pinned binary products, Frobenius filters, hyperelliptic covers, Mumford arithmetic, gluing or polarization candidates | `ecdlp-surfaces` | Use the matching native constructor/checker; hand a valid subgroup map to the solver stage |
| CM, Frobenius, endomorphisms, GLV/GLS or automorphism effects on cost | `curve-algebra` | Establish exact map/action and distinguish scalar-arithmetic gains from unknown-scalar recovery |
| A new source-versus-destination ECDLP measurement | `ecbench-measure` | Produce a sealed session, exact replay, three-arm table, controls and scoped result |
| Asking whether a measured method beats the established record | `ecbench-bounds` | Fit/check a bound and run a frozen challenge; do not infer a frontier move from an ad hoc comparison |
| Any experiment, benchmark or completion report | `report-evidence` | Preserve requested scope, raw outcomes, failures, limits and requirement status |
| Any research run, substantive search or result PR | `research-visuals` | Add the source-linked report, editable visual and inspected PDF; performance results also go under `benchmarks/` |

Use `ecdlp-surfaces` and `transfer` while finding and certifying routes. Use
`ecbench-measure` only after the candidate has an executable map and named
solver adapter. Use `ecbench-bounds` after, not instead of, the matched
measurement. The evidence and visual workflows apply throughout.

## Candidate record and decision

Each assessment must leave one of these narrow outcomes:

- `valid transfer, destination solver absent`;
- `solver applicable, net cost unmodeled`;
- `destination stage improvement only`;
- `modeled cold net advantage`;
- `measured cold net advantage in the tested scope`;
- `refuted by a named mathematical or measurement control`; or
- `unknown with the next discriminating obligation`.

Never replace an unknown with zero cost, extrapolate a toy exponent to a
cryptographic size without a separate model, or describe failure within a
bounded family as universal nonexistence.
