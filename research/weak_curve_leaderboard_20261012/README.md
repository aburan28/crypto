# WeakCurves evidence snapshot 03

This snapshot joins **814 exact elliptic models** from the source-locked
weak-curve database, **105 named curves** from the repository's curve registry
and the pinned pairing roster,
and **36 published odd-characteristic covering cases**. The named roster
includes Oakley Groups 3 and 4, with exact GHS invariant candidates and
separate WID1 finding IDs. It adds three exhaustively checked prime-field
examples and a twenty-four-criterion
registry for anomalous, smooth-order, MOV, Frey–Rück, supersingular,
binary-GHS, odd-characteristic GHS, cubic norm-one, quadratic-torus and
proper-subfield mechanisms, plus direct extension-field decomposition,
cover-and-decomposition, relation-support capacity and useful-rank checks, and cubic Jacobian-isogeny routes. The first [WeakCurves site](../../website/weakcurves/)
uses this exact export.

The new [pairing roster](pairing_curves.json) contains the exact Ethereum
BN_SNARK1 variant and BLS12-381, with separate G1, G2, and GT parameters.
The [arithmetic replay](validate_pairing_roster.py) checks both published
generators and the minimal root-of-unity degree. A [SIKEp434 protocol
case](protocol_cases.json) records the published auxiliary-torsion key-recovery
setting outside the ordinary-curve isogeny graph. DDH and static-DH criteria
state their challenge and oracle requirements; a curve-only record does not
assert those protocol inputs.

The [isogeny graph](ISOGENY_GRAPH.md) adds 170 verified-walk curve identities
to the 922-model atlas, grouping 1,092 nodes into 676 finite-field isogeny
classes and displaying 232 certified directed maps. Its criterion checks
power the site's vulnerability filters; check stages distinguish exact
conditions, geometry, explicit maps, candidates and structural screens.

The [WAM1 attack map](attack_map.json) adds 38 source-linked research checks
to the 24 existing criteria, giving the website 62 searchable checks across
curve structure, selected subgroups, auxiliary inputs, pairings, index calculus,
hyperelliptic divisor and Jacobian decomposition, endomorphisms,
implementations and isogeny protocols. Each WAM1 row states an
input prerequisite and a concrete evidence gate. These 38 rows are an
unassessed research queue: they create no curve-level WID1 finding. The
[validator](validate_attack_map.py) checks ID uniqueness, source fields,
source/site byte parity and certified graph edge consistency.

The new `ic_relation_support_rank` criterion separates geometric eligibility
from actual index-calculus relation yield. For an unordered positive
\(m\)-sum policy on \(B\) distinct usable classes, at most
\(\binom{B+m-1}{m}\) targets can be represented; a uniform accepted query set
of size \(A\) has per-query hit probability at most this support divided by
\(A\). The [P-192 Type I certificate](../p192_extension_cover_20261012/BRANCH_CONTINUATION.md)
supplies an exact instance: all 64 accepted guided queries had degree-three
prime-field divisor pairs, while the fixed 16-class three-sum base produced
zero verified relations and rank zero. The registry keeps the support bound,
accepted-query law, verified relation yield, and useful rank as separate gates.
The generic support and conditioning bounds are checked in
[Lean](lean/RelationSupport.lean); the [proof map](PROOF_MAP.md) states their
exact hypotheses and the remaining curve-specific bridge.

The [genus-128 GHS cost audit](ghs_genus128_cost_audit.json) resolves the
nominal factor-base and relation-query stages behind the two exact-invariant
`GF(4)` formula-advantage routes. Its independent integer calculation selects
degree bound 15 and 48,921,052 nominal half-split base places for both
`c2pnb176w1` and `c2pnb272w1`; actual cover places, relations, and rank are
separate witness gates in the binary GHS criterion.
The [run report and cost graph](reports/run-20261015T0854Z/REPORT.md)
show the full conditional sweep and its missing measurement stages.
The [source-lock refresh](reports/run-20261015T0929Z/REPORT.md) records the
same quantitative sweep against the newer main-branch registries.
The [cover-source alignment](reports/run-20261015T0950Z/REPORT.md) records the
current-main covers inventory used by the decomposition property index.

The [decomposition property index](decomposition_properties.json) makes thirteen
hyperelliptic checks filterable for all 922 atlas curves. It binds 917 exact
same-field cover certificates by ICV1 identity and full model hash, and
replays the explicit prime-field formula for the remaining five models under
their source field assumptions. The cubic and quadratic genus-three screens,
binary GHS invariants, Joux–Vitse intermediate-field and quadratic-sieve
screens, projected genera, divisor relation measurements,
Jacobian factorization, selected-subgroup transport and complete-cost gates
have separate states. Open or excluded branch states remain visible; the
same-field cover is a geometry property, not a weakness finding. The
[property validator](validate_decomposition_properties.py) audits the source
joins, status counts and site copy.
The paper's exact 149-bit worked example is separately preserved as a
[source-attributed case](published_cover_cases.json), with no atlas identity
or local replay assigned to it.

The default leaderboard contains 29 finding rows: two prime-field models
meeting a published DLP criterion, 24 exact small-field covering-model
signals, one embedding-degree candidate, and two Oakley GHS invariant
candidates. Its full-corpus view includes the 814 prior models, 105 named
curves, and all three prime-field examples. The 537 prior
models over 192–252-bit fields remain searchable with their actual cubic and
quadratic screens: 247 pass the cubic necessary condition, 137 fail its
depth screen, and 153 are outside the genus-three domain. The 36 published
cases are parameter families rather
than exact curve models.

Run from this directory:

```sh
cargo test --release --locked
cargo run --release --locked -- build snapshot.json
cargo run --release --locked -- check snapshot.json
python3 validate_pairing_roster.py
node --check ../../website/weakcurves/dist/app.js
python3 build_isogeny_graph.py --check
node --check ../../website/weakcurves/dist/graph.js
node --check ../../website/weakcurves/dist/attack-map.js
python3 validate_attack_map.py
cargo run --release --locked --bin ghs_genus128_cost_audit -- --check
python3 build_decomposition_properties.py --check
python3 validate_decomposition_properties.py
```

The Rust exporter reopens both frozen SQLite inputs and the JSON
registries, hashes their bytes, checks the exact import counts, enumerates
the three prime-field models, and emits stable JSON with no time- or
host-dependent fields. The graph builder binds the snapshot to verified
P-192/P-224 walk records and their receipt. Copy `snapshot.json`,
`criteria.json`, `attack_map.json`, `isogeny_graph.json`, and
`decomposition_properties.json`, `ghs_genus128_cost_audit.json`, and
`IDENTIFIERS.md` into
`../../website/weakcurves/dist/` when they change.
The site is static and reads the JSON files in the browser.

The [report](REPORT.md) gives the ranking contract, Oakley source checks,
exact example certificates, mechanism recognition gates, and next research
steps. The [WID1 specification](IDENTIFIERS.md) defines the weakness marker.
The [coverage audit](COVERAGE_AUDIT.md) maps prime-degree extensions,
hyperelliptic covers, and other Jacobian routes to their exact open witnesses.
The [hyperelliptic decomposition guide](HYPERELLIPTIC_DECOMPOSITION.md)
separates divisor relations from isogeny factorization and the elliptic/Prym
split induced by a cover. The coverage audit also records missing twist/input,
auxiliary-input, trace-zero, pairing-cost,
endomorphism, prime-field decomposition, and implementation-context checks.
The [proof map](PROOF_MAP.md)
specifies the Lean boundary. [VALIDATION.md](VALIDATION.md) records the build,
PARI replay, and site checks. This release does not infer a single security
score from field size or cover existence.
