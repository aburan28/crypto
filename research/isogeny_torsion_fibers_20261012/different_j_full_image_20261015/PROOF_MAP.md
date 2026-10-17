# Full-image finite support proof map

The [native support census](SUPPORT_V2.json) records scalar coordinates
for the eight selected base points in each arm. The readable argument is
finite: in a cyclic group of order 61, take each base scalar and its
negative, add every unordered pair modulo 61, discard zero, then compare
the resulting support sets. The [Lean module](lean/FullImageSupport.lean)
checks this exact finite calculation, using all ordered signed pairs;
ordered and unordered pair sets agree because addition is commutative.

| Report statement | Lean declaration | Hypotheses and status |
| --- | --- | --- |
| Full-image support consists of every nonzero residue except 1, 5, 56 and 60. | `FullImageSupport.fullSupport` | Proved for the eight literal scalar coordinates in `SUPPORT_V2.json`, modulo 61, with signed pair sums. |
| Projected-image support misses 2, 28, 33 and 59; destination-direct misses 24, 30, 31 and 37. | `projectedSupport`, `directSupport` | Proved for their literal scalar coordinates under the same finite definition. |
| Source and projected-image supports coincide. | `sourceEqualsProjected` | Proved because their frozen scalar-coordinate lists are identical. |
| Full-image and projected-image have 52 common residues and the specified four residues unique to each. | `fullProjectedDifference` | Proved over the two frozen scalar-coordinate lists. |
| Full-image and destination-direct have 52 common residues and the specified four residues unique to each. | `fullDirectDifference` | Proved over the two frozen scalar-coordinate lists. |
| The listed scalar coordinates correspond to the selected `F_(67²)` curve points, their signs and the explicit degree-three image map. | Open | The separate native arithmetic replay checks the bounded points and 178 relation witnesses. A faithful Lean `F_(67²)` curve/group model, point-to-scalar bridge, Vélu map and equality proof are not yet present. The finite Lean theorems do not assert these curve-level facts. |

The project pins Lean `v4.34.1` (toolchain commit
`5045d0056413266e57c625dcd7c365b10e377c52`) in `lean/lean-toolchain`
and has no external
dependencies. `lake build FullImageSupport` completed with the exact
output in [LEAN_BUILD_V1.stdout](LEAN_BUILD_V1.stdout) and
[LEAN_BUILD_V1.stderr](LEAN_BUILD_V1.stderr). `#print axioms` reports
no axioms for all six named theorems. `lean4checker --fresh` is unavailable
in this environment. The proof-map gap does not affect the separate
native scalar replay; it limits the formal claim to finite residue
combinatorics, leaving the geometric bridge open.
