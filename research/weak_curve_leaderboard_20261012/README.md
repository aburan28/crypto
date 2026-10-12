# WeakCurves evidence snapshot 01

This snapshot joins **814 exact elliptic models** from the source-locked
weak-curve database to **36 published odd-characteristic covering cases**.
It adds three exhaustively checked prime-field examples and a ten-criterion
registry for anomalous, smooth-order, MOV, Frey–Rück, supersingular,
binary-GHS, odd-characteristic GHS, cubic norm-one, quadratic-torus and
proper-subfield mechanisms. The first [WeakCurves site](../../website/weakcurves/)
uses this exact export.

The default leaderboard contains 27 exact signals: two prime-field models
meeting a published DLP criterion, 24 exact small-field covering-model
signals, and one embedding-degree candidate. Its full-corpus view includes
the 814 prior models and all three new prime-field examples. The 537 prior
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
```

The exporter reopens both frozen SQLite inputs, hashes their bytes, checks
the exact import counts, enumerates the three prime-field models, and emits
stable JSON with no time- or host-dependent fields. Copy `snapshot.json` and
`criteria.json` into `../../website/weakcurves/dist/` when they change. The
site is static and reads these files in the browser.

The [report](REPORT.md) gives the ranking contract, exact example certificates,
mechanism recognition gates, and next research steps. The [proof map](PROOF_MAP.md)
specifies the Lean boundary. [VALIDATION.md](VALIDATION.md) records the build,
PARI replay, and site checks. This release does not infer a single security
score from field size or cover existence.
