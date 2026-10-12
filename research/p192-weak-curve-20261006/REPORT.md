---
title: "P-192 weak-object study: singular validation and bounded CM relations"
author: "crypto research harness"
date: "2026-10-06"
lang: en-US
toc: true
toc-depth: 3
geometry: margin=0.72in
fontsize: 10pt
mainfont: "DejaVu Sans"
monofont: "DejaVu Sans Mono"
colorlinks: true
linkcolor: blue
urlcolor: blue
---

# Status and verdict

**Evidence cutoff: 2026-10-06. Both native producer runs and the independent
review are complete; the review verdict is `BREAKS`, with separate subclaim
dispositions.** The validator independently reproduced the singular lane's
fixed public synthetic scalar recovery, so the tested unchecked interface is
`IMPLEMENTATION_WEAK`, not weak P-192. In the valid CM lane the order facts,
enumeration counts, theorem bound, and all 103 emitted scalar relations were
reproduced. However, the frozen protocol did not specify the canonical-set
digest serialization, so exact set identity was not independently bindable.
Explicit maps, HOT/COLD calibration, and charged payoff also remain
unavailable. The complete-set negative and full CM experiment are therefore
`INDETERMINATE`; no P-192 structural weakness or speedup was found.

The work uses one fixed, public, synthetic scalar. It does not target deployed
keys, live services, or private data. Its purpose is to distinguish an
implementation-validation failure from a property of a valid elliptic curve.

Two hypotheses are deliberately kept separate:

1. **Singular implementation-validation lane.** A nodal cubic shares the
   P-192 field and coefficient `a=-3`, but it is not an elliptic curve and is
   not isogenous to P-192. The experiment asks whether an unchecked low-level
   multiplication surface can turn x-only confirmation tags into recovery of
   the fixed synthetic scalar. The strongest admissible verdict is
   `IMPLEMENTATION_WEAK`.
2. **Valid CM/isogeny-class lane.** An exact, bounded class-relation census is
   run in the endomorphism order attached to the valid P-192 isogeny class.
   A form relation is not a weakness: global orientation, explicit maps,
   subgroup transport, point replay, and fully charged payoff must all pass.
   The preregistered negative verdict is only
   `NO_WEAKNESS_FOUND_WITHIN_SCOPE` for the frozen box.

![Two lanes with separate identity and evidence gates. Blue is the valid elliptic/CM lane; purple/red is the singular implementation-validation lane.](object-attack-map.svg)

# Evidence labels

Statements in this report use these labels:

- **FROZEN:** fixed by approved protocol before execution.
- **DERIVED:** algebraically derived from stated inputs; not a run measurement.
- **INDEPENDENTLY VERIFIED — PRIOR:** replayed evidence from an earlier
  committed/archived run.
- **PRODUCER-MEASURED:** hash-bound native producer output that was not itself
  independently rederived.
- **VALIDATED:** independently rederived against the content snapshot. This
  label applies to the singular scientific joints and named CM subclaims only;
  it does not override the review's overall `BREAKS` verdict.
- **OPEN:** required evidence or capability is explicitly unavailable.

# Exact object identities

## Valid elliptic curve

The valid object is P-192, registered as
`icv1-fp192-t31607402316713927207482677199-52e4af59`.

- **Field prime `p` (FROZEN):**

  `6277101735386680763835789423207666416083908700390324961279`.
- **Short-Weierstrass coefficients (FROZEN):** `a=p-3`; exact `b`:

  `b=2455155546008943817740293915197451784769108058161191238065`.
- **Prime group order `n` (FROZEN):**

  `6277101735386680763835789423176059013767194773182842284081`.
- **Trace `t=p+1-n` (DERIVED):** `31607402316713927207482677199`.
- **Compact representation:** `EC1P192Cp192h5531c4a08bdb`.
- **Full representation UID:** prefix `urn:ec-record:1:sha256:` followed by

  `5531c4a08bdb64b6e86a6e30e9a08aa57edef7af15ac5f6d4d2a83a53bf2f646`.

The exact registry source is
[`docs/curves/registry.json`](../../docs/curves/registry.json) at base revision:

`0dc6b7f4a03253e7c7641f99700f7f423f8d07e0`.

Its source SHA-256 is recorded in the artifact index.

## Singular companion: a typed non-curve object

**DERIVED.** Define over the same field

$$
S: y^2=x^3-3x-2=(x+1)^2(x-2).
$$

For a short-Weierstrass model, the discriminant is
`-16(4a^3+27b^2)`. At `a=-3,b=-2`, the inner term is
`4(-27)+27(4)=0`; therefore `S` is singular, with node `(-1,0)`. It is
not assigned ICV1, EC1, or a curve UID. The exact non-curve preimage is stored
in [`singular-model.canonical.json`](singular-model.canonical.json) with type
`singular-model/v1`; its SHA-256 is
`305129485b8a281f5cd73a7a0ec4d82a7f84f432847c1b9d5eb1734983db72b3`.

The probe `R=(66,536)` satisfies
`536^2=66^3-3*66-2=287296`, is away from the node, and is not on the valid
P-192 model because the valid `b` differs. Sharing `p` and `a` does not define
an isogeny, an isomorphism, or a curve-family member.

# Experiment A — native singular-companion recovery

Protocol: `EXP-SCURVE-29040c`, approved by `DEC-20261006-af8214`.
Exact inspected protocol SHA-256:
`ae25388412acb9fbc8a0387f8d4782fe215a493b6bccb6f1f0fa49f8f1b3fde6`.
Protocol v2 amendment SHA-256 is
`c5be1d42be520cbadec68c76928f579559f8e1ebfb7fa9aa93c58fa0853512a0`;
the x-only encoding correction is `CORR-20261006-500337`. Producer run
`RUN-SCURVE-af5caf` used code commit `6d9d1ed7674cead991549d2060de3e91be1d9a3c`
and is bound by the content snapshot described in
[`ARTIFACT_INDEX.md`](ARTIFACT_INDEX.md).

## Frozen algebra and attack path

On the smooth locus set `t=y/(x+1)`. The inverse parameterization and group
law used by the independent arithmetic control are

$$
P(t)=(t^2+2,\;t(t^2+3)), \qquad
t_1\star t_2=\frac{t_1t_2-3}{t_1+t_2}.
$$

The identity is `t=∞`, inverse is `t -> -t`, and `R=P(8)=(66,536)`.
Because `p mod 12 = 11` and `-3` is a nonsquare, this is the nonsplit nodal
torus with order

$$
N=p+1=
2^{64}\cdot3\cdot5\cdot17\cdot257\cdot641\cdot65537\cdot274177
\cdot6700417\cdot67280421310721.
$$

The frozen point-order certificate must prove `[N]R=O` and
`[N/q]R != O` for every distinct prime factor `q`. It may not be inferred from
the factorization alone.

The largest factor and its cofactor are

```text
qmax = 67280421310721             (log2 ≈ 45.93525)
M    = 93297598515282145091939986791922449178951680
```

The x-only oracle is sign-invariant. Nested composite-order probes therefore
retain one coherent pair `{r,-r mod m}` rather than independently choosing a
sign for every factor. The frozen sequence contains 64 lifts by two and eight
odd factors, so the expected primary-oracle count is exactly 72. Candidate
comparisons are counted separately and must not be relabelled as oracle calls.

After recovery modulo `M`, the legitimate public key `Q=[d]G` on P-192 is
required to solve the remaining interval. One shared negation-BSGS table is
frozen at `table_bits=24`, `m=5,800,019`, stride `11,600,037`, no more than
`5,800,018` finite baby entries, and no more than `5,800,018` giant positions
per orientation. Every fingerprint hit requires exact replay.

## Fixed public synthetic fixture

- **Seed:** `P192-SINGULAR-V1`; hex
  `503139322d53494e47554c41522d5631`.
- **Scalar mapping:** `1 + OS2IP_BE(SHA256(domain || seed)) mod (n-1)`.
- **Expected `d` in hex:**

  `8b70f85429c81223e88be1a60590627a3c8e558a24dd5180`.
- **Expected `d` in decimal:**

  `3419090462492416794041408622244871729251768983841006244224`.
- **Unordered residues modulo `M`:**

  - `78222222987683236539393919745863586265518464`.
  - `15075375527598908552546067046058862913433216`.
- **True residual `k`:** `36647143301682`. Check only after recovery; never
  use it to prune the search.

## Success, falsification, and controls

An `IMPLEMENTATION_WEAK` result requires all of the following, not merely
successful torus arithmetic:

- the exact P-192 and singular-point/order certificates pass;
- the actual unchecked `scalar_mul_secret` surface is used for every primary
  query and exposes only the frozen x-derived HMAC tag;
- one coherent plus/minus class is recovered after exactly 72 primary calls;
- the legitimate P-192 residual solve returns the fixed `d`;
- independent final replay proves `[d]G=Q`;
- the safe high-level `ecdh_raw` path rejects every singular probe before
  multiplication; and
- exhaustive `p=59`, wrong-parameter, mutated-transcript, full-point positive,
  and independent-torus controls all give their preregistered verdicts.

If instrument controls pass but the x-only transcript cannot recover the
coherent class, or neither residual orientation replays to `Q`, the hypothesis
is falsified. A timeout, OOM, or missing receipt is `INDETERMINATE`, not evidence
that the valid curve is strong.

## Producer result and independent validation

**VALIDATED AT THE STATED INTERFACE; OVERALL REVIEW `BREAKS` ON ARCHIVE
METADATA.** The native run recovered the exact fixed scalar listed in the
fixture above; a source-blind checker independently reproduced all 72 tags,
1,245,111 candidate comparisons, the coherent global sign, the BSGS equation,
the residual `k`, the recovered `d`, and `[d]G=Q`. Order, torus, transcript,
mutation, `p=59`, safe-path rejection, and full-point-model controls passed.
Both singular scientific joints are `HOLDS`. No speedup is asserted.

- **Valid P-192 reference:** approximately `sqrt(n)` group work; a matched
  native reference was not run.
- **Nested residues:** 72 primary x-only queries and 1,245,111 candidate
  comparisons, kept as different counted units.
- **Shifted-negation-v2 BSGS:** 5,800,018 baby entries. The exact hit has
  orientation `-1`, giant index 3,159,226, baby index 989,698, and sign `-1`.
- **Whole pipeline:** success; 20.809205507 seconds wall time, 133,616 KiB
  maximum RSS, and every frozen control passed.

Correctness is certificate-backed at the producer layer. A ratio to a matched
reference remains unavailable and is not meaningful until every phase is
expressed in the same counted unit.

Wall time is secondary engineering evidence. A performance claim additionally
requires the repository's isolated A/A and interleaved paired-measurement
protocol; this single recovery run is a deterministic feasibility result, not
a timing confidence interval or prevalence estimate. The supported object is
the exposed implementation interface. The singular cubic remains neither an
elliptic curve nor a P-192 isogenous representative.

The original snapshot receipt misspelled the singular producer-origin commit.
All 62 content hashes and both manifest chains passed, but strict provenance
therefore broke. The immutable receipt remains unchanged; additive correction
`CORR-20261006-cfc18a` overlays the actual commit
`8b81f8efa613766a76eaf46465892f22b37fad5d`. This administrative defect does
not erase the independently reproduced arithmetic, and the correction does not
turn the overall review verdict into a pass.

# Experiment B — exact weighted CM class-relation census

Protocol: `EXP-SCURVE-647ade`, approved by `DEC-20261006-af8214`.
Exact inspected protocol SHA-256:
`0fe122a0c6be3e2118b4d9c276a1085285c7e6a2cee38d15bbeed8099fa17145`.
Protocol v2 amendment SHA-256 is
`029fcf0ed147af8eb590e770bed9c2a023eadfd3bca8567d0297caa6130d3a33`.
Producer run `RUN-SCURVE-81c5e4` used corrected code commit
`78ba2da296c505621feb0b05c11507f3b5ae6669` and is bound by the content
snapshot described in [`ARTIFACT_INDEX.md`](ARTIFACT_INDEX.md).

## Frozen order and generator set

**PRODUCER-MEASURED; certificate emitted.** The Frobenius discriminant is

$$
D=t^2-4p=-24109379060336110122544161233113975664949272517896865359515
=-5\cdot11\cdot31\cdot C,
$$

where
`C=14140398275856956083603613626459809774163796198179979683`.
The run emitted and replayed a recursive Pocklington certificate for `C`, then
certified that `D` is squarefree and fundamental, the Frobenius conductor is
one, and `End(E)=Z[pi]=O_D`. The native p/n checks remain fixed-base
probable-prime checks with pinned SEC 2/FIPS provenance, not newly constructed
recursive proofs.

The frozen ramified roots are `5:[2]`, `11:[3]`, and `31:[7]`. The frozen
split roots are:

| `ell` | roots of `X^2-tX+p mod ell` |
|--:|:--|
| 13 | 2, 5 |
| 23 | 21, 22 |
| 37 | 12, 35 |
| 43 | 8, 26 |
| 73 | 60, 67 |
| 89 | 6, 83 |
| 101 | 17, 70 |
| 103 | 5, 36 |
| 107 | 56, 68 |
| 113 | 26, 42 |

Each root and orientation must be recomputed. A positive split generator uses
the smaller root; global conjugates are canonicalized by requiring the first
nonzero split exponent to be positive. Ramified exponents are unsigned.

## Exact boundary and theorem control

The census includes every canonical exponent vector with
`L=product ell^|e_ell| <= 2^48`; empirical weights may rank but may not prune.
The preregistered cardinalities are 9,948,061 oriented vectors, 4,974,348
canonical vectors including zero, and 4,974,347 canonical nonzero vectors,
under a hard cap of 5,000,000 canonical states. A differing count invalidates
the run.

For a non-scalar `alpha=u+v*pi`, the protocol must independently certify

$$
N(\alpha)\ge \left\lceil\frac{|D|}{4}\right\rceil
=6027344765084027530636040308278493916237318129474216339879
\approx2^{191.941}.
$$

If all premises pass, the exact `2^48` lane can contain only scalar ramified
relations. This is a proof-controlled scope statement, not a claim about all
relations or all attack mechanisms.

![Derived CM degree/norm values. The plot contains no empirical run or timing data.](cm-degree-frontier.svg)

The controls are `g5^2`, `g11^2`, and `g31^2` on P-192; the positive
non-scalar toy relation `g2^3=1` at `D=-23`; and the symbolic P-192 relation
`g5*g11*g31*gC=1` with `alpha=sqrt(D)=2*pi-t`. The last relation is genuine
algebraically but is outside the buildable `ell<=113` set: a `C`-degree edge
is unavailable, and its action on `E(F_p)[n]` does not by itself improve a
DLP. It must remain a rejected symbolic control.

## Result and payoff gates

**PRODUCER-MEASURED.** Complete exact algebra census; no non-scalar relation in
the frozen box. Both enumeration orders returned 9,948,061 oriented vectors,
635 conjugation-fixed vectors, 4,974,348 canonical vectors including zero, and
4,974,347 canonical nonzero vectors, with identical set digests. The exact
`L=2^48` layer is empty with the canonical empty-input SHA-256. There are 104
scalar-principal states including zero, hence 103 emitted nonzero relations;
all are scalar ramified, all 103 HNF/form replays passed, and zero are
non-scalar.

**INDEPENDENT REVIEW.** The order/discriminant/conductor facts, generator
roots, theorem lower bound, all enumeration cardinalities, all 103 emitted
relation certificates, and the mutation controls were independently
reproduced. The 103 emitted relations are therefore validated as scalar.
The review could not reproduce the producer's canonical-set digest from the
frozen protocol, because the protocol names a digest but never defines its
domain separator, vector serialization, XOR/sum chunking, or byte order. The
producer source uses a domain-separated construction, but learning that after
blind sealing cannot retroactively repair the protocol. Thus matching counts
do not independently certify exact set identity.

- **Order/CM proof — pass:** pinned tuple, Hasse uniqueness, recursive `C`
  proof, fundamental `D`, and conductor one.
- **Exact membership — inconclusive:** forward/reverse producer counts and
  digests agree, and counts independently reproduce, but the digest encoding
  was not frozen and exact set identity is not independently bound.
- **Algebraic replay — validated:** 103 scalar-ramified records accepted and zero
  rejected.
- **Global map realization — OPEN:** no non-scalar candidate exists, and
  globally oriented explicit maps and terminal isomorphisms are not implemented.
- **HOT/COLD accounting — OPEN:** zero of 368 required calibration records;
  missing tuples are null, not zero.
- **Whole payoff — not passed:** no charged attack or canonicalization
  inequality exists.

The exact search took 28.239448942 seconds and the full verifier replay took
11.462448644 seconds; aggregate phase wall time was 39.733819654 seconds with
639,024 KiB peak RSS. These are single-run engineering observations, not a
speed comparison. The degree-linear proxy was not substituted for missing
HOT/COLD operation costs.

Two preserved attempts found defects in the research harness rather than in
P-192: an incorrect Bézout coefficient plus duplicated factor in generic
quadratic-form composition, and a one-ULP JSON float round-trip mismatch that
the exact verifier correctly rejected. Commits `bcd59e1` and `78ba2da` repaired
them and added full-census, one-ULP, and NaN regressions before attempt 3.

The producer's bounded negative remains a useful observation, but independent
review does not license `NO_WEAKNESS_FOUND_WITHIN_SCOPE` for the complete set.
Because exact set identity, map/calibration/payoff evidence, and map-specific
controls are unavailable, the CM experiment is `INDETERMINATE`, not evidence
that P-192 is globally strong or weak. No non-scalar relation was found in the
producer census, and none appears among the independently replayed emissions.

# What will count as “weak”

The word *weak* is reserved for an auditable statement with an exact object,
threat model, boundary, and verified advantage.

- **`INVALID_INSTANCE`:** a purported elliptic-curve candidate is singular,
  malformed, or otherwise fails identity preconditions. The singular companion
  has this object-level classification even though it remains useful for an
  implementation-validation test.
- **`IMPLEMENTATION_WEAK`:** the actual unchecked surface, fixed x-only oracle,
  end-to-end planted-key recovery, independent replay, and every control pass.
  The singular cubic is still not a weak elliptic or isogenous curve.
- **`CLASS_WEAK`:** a certified attack applies across the declared valid
  isogeny class, with subgroup transport, replay, and material advantage.
- **`WEAK_REPRESENTATIVE_EXISTS`:** at least one valid class representative has
  a certified weakness, without claiming every representative shares it.
- **`SOURCE_TRANSFER_WEAK`:** a certified map transfers a concrete weakness
  back to the named source curve with all conversion costs charged.
- **`NO_WEAKNESS_FOUND_WITHIN_SCOPE`:** the frozen CM box is exhausted with
  exact counts, independently bindable set identity, proofs, and replay, and no
  candidate passes the buildability and payoff gates. This gate was not met in
  the present CM review; it is a criterion, not the result.
- **`INDETERMINATE`:** mandatory map, cost, control, or replay evidence is
  unavailable. Resource exhaustion or unsupported tooling is not strength
  evidence.

The safe public ECDH path's rejection of every singular probe is recorded as a
passing control, not as a separate final verdict.

P-192's roughly 96-bit generic security level and legacy status are baseline
properties, not discoveries of this campaign.

# Confidence model

## Deterministic claims

Primality, factorization use, point order, class-group products, HNF identities,
enumeration membership, map equations, and final scalar replay are deterministic
claims. Statistical confidence is the wrong tool for them. Confidence comes
from self-contained certificates, exact counts/digests, a second enumeration
order, mutation-negative controls, and an independent verifier that does not
trust the producer's intermediate state.

The CM producer artifacts support a finite-census observation and independently
reproduced counts, but the review does not support the stronger complete-set
statement because the digest encoding was not frozen. No confidence interval
can repair an underspecified deterministic identity check; a new additive
protocol must define the bytes and rerun or re-digest them independently.

## Timing and implementation variability

Any wall-clock comparison must first measure A/A spread, then interleave
baseline and candidate for at least five paired rounds under the repository's
isolated runner. Report counted native operations as primary, plus median,
minimum, contention flags, and a 95% paired confidence interval. A runtime
improvement is supported only if the interval excludes no improvement and all
outputs/digests match. The recorded 20.809-second recovery and 39.734-second CM
pipeline are single executions, so no timing comparison or confidence interval
is asserted.

## Generalizing beyond one fixture

One fixed synthetic scalar can demonstrate a mechanism and exact recovery; it
cannot estimate a prevalence rate across keys, implementations, or curves. A
later population claim would require a preregistered sampling frame, independent
holdouts, multiple keys per implementation, explicit failure retention, exact
binomial intervals, and multiplicity control across curve/candidate searches.
No such population inference is made here.

# Prior evidence and non-duplication

The prior committed P-192 isogeny walk
[`research/isogeny_walk_p192_20261004/README.md`](../isogeny_walk_p192_20261004/README.md)
is **INDEPENDENTLY VERIFIED — PRIOR**, not a result of these experiments. Its
20,000-curve run is `p192-71c04205135fcf8f`, archived under
`s3://crypto-autoresearcher/isogeny-walk/runs/p192-71c04205135fcf8f`.
It reported 20,000 valid curve records, 128,700 edges, zero failures, and replay
pass, but explicitly established no DLP weakness or per-curve IC difference.
The new exact CM census asks a different, relation-focused question and must not
rewrite that frozen evidence.

The repository's typed-link rules
([`docs/curves/ic/curve-links/README.md`](../../docs/curves/ic/curve-links/README.md);
source hash in the artifact index)
require exact endpoints, kernels/maps, dual status, subgroup transport, and
replay before an isogeny route can carry a DLP claim. That rule is the basis of
the CM payoff gate and the red “NOT an isogeny” separation in the diagram.

# Open obligations

- Independently accept or reject additive custody correction
  `CORR-20261006-cfc18a`; do not edit the immutable snapshot receipt or erase
  the historical review break.
- Freeze the CM canonical-vector byte encoding, domain separator, digest
  aggregation, and byte order in a prospective amendment, then rerun or
  independently re-digest the complete set.
- Keep explicit global maps, terminal isomorphisms, HOT/COLD calibration, and
  fully charged payoff open. Do not fabricate them from degree or wall time.
- Preserve separate evidence records and the post-result `revise` decision;
  never merge the singular implementation conclusion into the valid
  CM/isogeny conclusion.
- Regenerate and visually inspect this report PDF after the final archive
  hashes are incorporated.
- Do not extend the conclusion to sect113r1, other standardized curves, other
  APIs, or a population of implementations without separately frozen inputs
  and evidence.

# Sources and provenance

## Repository and protocol sources

1. P-192 registry entry: [`docs/curves/registry.json`](../../docs/curves/registry.json),
   exact EC1/UID above; source hash and revision recorded in the artifact index.
2. Repository curve parameters: [`src/ecc/curve_zoo.rs`](../../src/ecc/curve_zoo.rs),
   source hash recorded in the artifact index.
3. Prior P-192 walk: [`research/isogeny_walk_p192_20261004/README.md`](../isogeny_walk_p192_20261004/README.md),
   with report, run, and store hashes in the artifact index.
4. Frozen experiment records `EXP-SCURVE-29040c` and `EXP-SCURVE-647ade`,
   their v2 amendments, archive correction `CORR-20261006-9ef467`, and the
   additive origin overlay `CORR-20261006-cfc18a`, approved or recorded under
   the campaign ledger; exact hashes are recorded above and in the artifact
   index.
5. Independent review `TASK-20261006-7d34e7`, evidence records
   `EV-SCURVE-44f051` and `EV-SCURVE-8e5c70`, and post-result decision
   `DEC-20261006-581065`; their exact hashes and archive commit are recorded in
   the artifact index.
6. Evidence and accounting rules: [`AGENTS.md`](../../AGENTS.md), SHA-256
   recorded in the artifact index.
7. The canonical `audit-curve` skill and `KN-TECH-6a2ef9` weak-curve audit
   synthesis in `crypto-autoresearcher`; these define the verdict taxonomy and
   keep invalid-input, representative, class-wide, transfer, and implementation
   claims separate.

## External parameter and terminology sources

- NIST, [FIPS 186-4](https://csrc.nist.gov/pubs/fips/186-4/final), archived
  standard containing the legacy NIST prime-curve parameters.
- SECG, [SEC 2 version 2.0](https://www.secg.org/sec2-v2.pdf), secp192r1
  domain parameters.
- Andrew V. Sutherland, [“Isogeny Volcanoes”](https://msp.org/obs/2013/1-1/obs-v1-n1-p25-s.pdf),
  for isogeny-volcano terminology. The report does not treat terminology as a
  certificate for any new edge.

The complete artifact set and remaining technical obligations are in
[`ARTIFACT_INDEX.md`](ARTIFACT_INDEX.md).
