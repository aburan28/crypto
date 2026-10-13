# WeakCurves evidence snapshot 03

This snapshot joins **814 exact elliptic models** from the source-locked
weak-curve database, **105 named curves** from the repository's curve registry
and the pinned pairing roster,
and **36 published odd-characteristic covering cases**. The named roster
includes Oakley Groups 3 and 4, with exact GHS invariant candidates and
separate WID1 finding IDs. It adds three exhaustively checked prime-field
examples and a twenty-three-criterion
registry for anomalous, smooth-order, MOV, Frey–Rück, supersingular,
binary-GHS, odd-characteristic GHS, cubic norm-one, quadratic-torus and
proper-subfield mechanisms, plus direct extension-field decomposition,
cover-and-decomposition, and cubic Jacobian-isogeny routes. The first [WeakCurves site](../../website/weakcurves/)
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

The [WAM1 attack map](attack_map.json) adds 36 source-linked research checks
to the 23 existing criteria, giving the website 59 searchable checks across
curve structure, selected subgroups, auxiliary inputs, pairings, index calculus,
hyperelliptic divisor and Jacobian decomposition, endomorphisms,
implementations and isogeny protocols. Each WAM1 row states an
input prerequisite and a concrete evidence gate. These 36 rows are an
unassessed research queue: they create no curve-level WID1 finding. The
[validator](validate_attack_map.py) checks ID uniqueness, source fields,
source/site byte parity and certified graph edge consistency.

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
```

The Rust exporter reopens both frozen SQLite inputs and the JSON
registries, hashes their bytes, checks the exact import counts, enumerates
the three prime-field models, and emits stable JSON with no time- or
host-dependent fields. The graph builder binds the snapshot to verified
P-192/P-224 walk records and their receipt. Copy `snapshot.json`,
`criteria.json`, `attack_map.json`, `isogeny_graph.json`, and
`IDENTIFIERS.md` into `../../website/weakcurves/dist/` when they change.
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
