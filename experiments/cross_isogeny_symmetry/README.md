# Cross-isogeny symmetry: experimental design and execution

Companion to CryptoAutoResearcher PR #2059. This PR provides a **runnable Sage toy pilot** and a staged, falsifiable experiment design, not just research goals.

## Run pilot
```sh
sage -python experiments/cross_isogeny_symmetry/pilot.py --prime 7 --extension 2 --ell 3 --seed 11 --output /tmp/fiber-pilot.json
sage -python experiments/cross_isogeny_symmetry/pilot.py --prime 11 --extension 2 --ell 5 --seed 17
```
Requires a SageMath installation. No Sage execution or CI result is claimed by this PR.

## Pilot protocol (implemented)
- Generate ordinary short-Weierstrass E/F_(p^n) with coefficients in F_p and a degree-ell isogeny phi.
- Exhaustively enumerate small curve points and check phi commutes with both p-Frobenius (valid here because coefficients/maps are over F_p) and q-Frobenius.
- Enumerate Frobenius orbits, compare matched-size factor bases defined by trace(x(P))=0 versus trace(x(phi(P)))=0.
- For identical sampled target points, exhaustively count two-term decompositions P+R=Q in each factor base; report density, hit counts and timings.
- JSON results carry parameters, seed, curve and codomain, measured counts, timing, and an explicit warning that this is *not* a full index-calculus speedup.

## Next experiments: precise independent variables and controls
| ID | Treatment | Matched control | Endpoint |
|---|---|---|---|
| X01 | Phase-aware joint Frobenius/isogeny orbit keys | Frobenius-only keys | construction + relation cost |
| X02 | SAT symmetry-breaking constraints with phase payload | unbroken SAT | verified independent relations/s |
| X03 | Cross-isogeny pullback factor base | same-cardinality ordinary factor base | relation probability and total cost |
| X04 | Invariant-ring and equivariant Macaulay blocks | direct F4/F5 | degree of regularity and elimination |
| X05 | Syzygy pullback with saturation | direct syzygies | reusable identities and net solver cost |
| X06 | Curve-model switching per attack stage | best single model | end-to-end cost incl transfers |
| X07 | Cached isogeny-path elimination templates | cold uncached solver | amortization break-even |
| X08 | Trace/norm with phase recovery | explicit coordinate equations | reconstruction-adjusted cost |
| X09 | Torsion-coset slicing | unsliced relation variety | independent relation yield |
| X10 | Weil restriction and fiber-product decomposition | standard Semaev formulation | total verified relation cost |
| X11 | Solver portfolio trained on other classes | fixed tuned solver | held-out class wall-clock |
| X12 | Conductor-stratum discontinuity search | isomorphic and same-level controls | replicated d_reg/yield change |

## Experiment protocol for X01–X12
1. Generate certified ordinary isogeny classes and `ell`-isogeny paths, including conductor valuations and trace/order certificates.
2. Start with toy p^n and binary m=small for exact exhaustive verification; advance to m=31,51,53,83 when mathematically suitable.
3. Generate fixed target sets and pair each treatment with same field, trace, subgroup, factor-base size, solver version, monomial order, seeds, time budget and hardware.
4. Record preprocessing, map construction, fiber construction, factor-base membership, misses, timeouts, duplicates, relation verification, independence/rank, linear algebra, peak RAM and final solve.
5. Report both *cost per new independent verified relation* and *full verified ECDLP solve cost*. Do not equate two-term toy relation counts with independent relations.
6. Use multiple runs, bootstrap confidence intervals, and hold out entire isogeny classes to prevent leakage.
7. Reject claimed speedups that disappear after map transfer, failed targets, preprocessing, or density effects.

## Mathematical validation
- Degree of separable isogeny equals geometric kernel size; rational fibers may be smaller.
- For prime-order subgroup r coprime to ell, isogeny restriction is injective, so no generic r-fold compression.
- Frobenius action must preserve the curve/model and map; distinguish p-Frobenius and q-Frobenius.
- If rational substitutions introduce denominators, saturate and check exceptional loci before transporting polynomial ideals.
- Equivariant block decompositions require a group-stable ideal and compatible linear action.

## Acceptance gates
Pilot: deterministic JSON, exact commutation checks, nonzero comparable factor bases, reproducible small-field instances. Next: Sage tests for known curves, verified kernel cosets and affine actions. Main: >=2x improvement in total cost per independent relation on at least two larger valid settings, replicated and with uncertainty intervals; otherwise publish negative result.

## Reproducibility and limitations
This branch includes a pilot script and this design. It does **not** include F4/F5/SAT integration, populated results, or claims of running Sage. The second command is illustrative and may fail to find an isogeny for a particular seed; the pilot raises an explicit error in that case.
