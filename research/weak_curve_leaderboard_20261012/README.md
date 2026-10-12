# WeakCurves evidence snapshot 02

This snapshot joins **814 exact elliptic models** from the source-locked
weak-curve database, **103 named curves** from the repository's curve registry,
and **36 published odd-characteristic covering cases**. The named roster
includes Oakley Groups 3 and 4, with exact GHS invariant candidates and
separate WID1 finding IDs. It adds three exhaustively checked prime-field
examples and a sixteen-criterion
registry for anomalous, smooth-order, MOV, Frey–Rück, supersingular,
binary-GHS, odd-characteristic GHS, cubic norm-one, quadratic-torus and
proper-subfield mechanisms, plus direct extension-field decomposition,
cover-and-decomposition, and cubic Jacobian-isogeny routes. The first [WeakCurves site](../../website/weakcurves/)
uses this exact export.

The [isogeny graph](ISOGENY_GRAPH.md) adds 170 verified-walk curve identities
to the 920-model atlas, grouping 1,090 nodes into 674 finite-field isogeny
classes and displaying 232 certified directed maps. Its criterion checks
power the site's vulnerability filters; check stages distinguish exact
conditions, geometry, explicit maps, candidates and structural screens.

The default leaderboard contains 29 finding rows: two prime-field models
meeting a published DLP criterion, 24 exact small-field covering-model
signals, one embedding-degree candidate, and two Oakley GHS invariant
candidates. Its full-corpus view includes the 814 prior models, 103 named
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
node --check ../../website/weakcurves/dist/app.js
python3 build_isogeny_graph.py --check
node --check ../../website/weakcurves/dist/graph.js
```

The Rust exporter reopens both frozen SQLite inputs and both existing JSON
registries, hashes their bytes, checks the exact import counts, enumerates
the three prime-field models, and emits stable JSON with no time- or
host-dependent fields. The graph builder binds the snapshot to verified
P-192/P-224 walk records and their receipt. Copy `snapshot.json`,
`criteria.json`, `isogeny_graph.json`, and
`IDENTIFIERS.md` into `../../website/weakcurves/dist/` when they change.
The site is static and reads the JSON files in the browser.

The [report](REPORT.md) gives the ranking contract, Oakley source checks,
exact example certificates, mechanism recognition gates, and next research
steps. The [WID1 specification](IDENTIFIERS.md) defines the weakness marker.
The [coverage audit](COVERAGE_AUDIT.md) maps prime-degree extensions,
hyperelliptic covers, and other Jacobian routes to their exact open witnesses.
The [proof map](PROOF_MAP.md)
specifies the Lean boundary. [VALIDATION.md](VALIDATION.md) records the build,
PARI replay, and site checks. This release does not infer a single security
score from field size or cover existence.
