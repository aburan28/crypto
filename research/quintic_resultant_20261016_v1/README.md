# Native inversion-free quintic predicates — 8x research target

This is a new exact **component**, not a complete generalized-rootfinding or
index-calculus pipeline and not a novelty claim. The immutable original 2x
release gate is preserved; the updated requested end-to-end target is **8x**.

For characteristic-two monic H=X^5+e1 X^4+e2 X^3+e3 X^2+e4 X+e5, compute
U=e3+e1 e2, V=e5+e1 e4, D=V^2+e2 U V+e4 U^2. H is squarefree iff D !=0.
This specializes a classical resultant and avoids field inversions, square
roots and generic GCDs. See PROOF.md for exceptional cases and the proof.
Other characteristics use a division-free polynomial remainder sequence.
Repeated-root support always retains the exact nilpotence fallback.

Native exhaustive controls independently test all4,392 monic quintics over
GF4/F3/F5 against degree1/2 square divisors and literal allowed-root products.
Supplemental GF9/GF25 controls test1,024 further monic inputs and336 repeated
irreducible-quadratic ghost cases. Native Rabin tests certify9 binary moduli.
The finite panel includes20 fields: five NIST primes, four small primes,
nine binary extensions through571 including83, and GF9/GF25. Wide coefficients
use full-width seeded field elements. No curve pipeline integration is claimed.

Run from this directory: `cargo run --release --locked -- controls`, or
`cargo run --release --locked -- native-quintic-resultant-calibration-20261016-v1 2026101601`.
Commands for the confirmation seed2026101602 were fixed before measurements.
This reproduces native predicate diagnostics; it is not a qualified timing run.

The original user-mandated lab scheduler only accepts Python module commands.
The retained `coordinator_receipts/launch.py.snapshot` contains only os.execv
into the frozen Rust executable; there is no Python mathematical, oracle,
timing or serialization execution path in this study. The separate original
lab's Python coordinator was used as explicitly instructed by the user.
All actual algorithm/harness/controls are native. Cargo sources/lock and binary
hashes are frozen; prior preparation mistakes remain quarantined.

Current boundary: supplied quintic classification, degree5,12 inputs/field,
three alternating internal rounds, at most90seconds/4GiB/one CPU per worker.
These are not three cold full-solver runs, not a matched F4/F5/F6 comparison,
not hard20-system evidence, and not a claim of uniform8x performance.

Novelty review remains open: fraction-free PRS and resultant identities are
classical, treated in von zur Gathen and Gerhard, Modern Computer Algebra,
3rd ed. (author-maintained book site https://www-math.upb.de/mca/).
No priority or asymptotic claim is made. No Lean proof is claimed; the exact
readable proof and independent controls remain the current proof evidence.

Next: price the complete support recurrence and root extraction with these
predicates before adoption; if generic PRS regresses reject unconditional use.
A separate hard-family branch must reduce the remaining search, not merely
classify its supplied candidates faster. Budget pauses do not close the goal.
