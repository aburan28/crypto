# Native endomorphism search

Turns the decision tree into exact bounded searches and verified explicit maps.
The execution path, correctness tests and replay verifier are Rust. The small
prime-field arithmetic is an independent exhaustive-check oracle; metadata
hashing reuses the repository's SHA-256 implementation. Build the standalone
Cargo package and invoke its binary directly:

```sh
cargo build --manifest-path tools/endomorphism-search/Cargo.toml --locked --release
tools/endomorphism-search/target/release/endomorphism-search scan --discriminant=-619 --degree-bound 1000
tools/endomorphism-search/target/release/endomorphism-search sweep --max-abs-discriminant 1000
tools/endomorphism-search/target/release/endomorphism-search probe --p 167 --a 25 --b 36 --count-ops
tools/endomorphism-search/target/release/endomorphism-search demo --count-ops --out /tmp/native-demo.json
cargo test --manifest-path tools/endomorphism-search/Cargo.toml --locked --release
bash tools/endomorphism-search/check.sh
```

`scan` enumerates both signs and conjugates of non-scalar elements of a declared
quadratic order through a degree bound, with exact integer square-root bounds.
It computes the conductor, geometric units and primitive reduced class forms.
Facts about the order remain conditional for any curve until bound to that
curve's full ring. For discriminant −619 the minimum non-scalar degree is 155,
there are two units and five ideal classes. Class representatives of leading
coefficient 5 or 7 are not self-endomorphisms of those degrees.

`sweep` searches all valid discriminants in a declared range. Work, output,
kernel, order and finite time limits produce explicit incomplete statuses.
An incomplete scan cannot establish absence. To keep native integer arithmetic
auditable, discriminant magnitude and degree bound are capped at 10^12; the
range sweeper caps its range at 100000. Limits reject NaN and infinity.

`probe` enumerates points on a nonsingular short Weierstrass model over a prime
field with `3 < p <= 16381`. It discovers scaling automorphisms and rational
cyclic Vélú kernels of degrees 2, 3, 5, 7, 11 and 13, then finds an isomorphism
back to the starting model. Every rational image is checked, 64 deterministic
homomorphism pairs are sampled, and the subgroup action is checked exhaustively.
The algebraic kernel, quotient and isomorphism conditions are separate from the
finite checks; samples alone are not a geometric proof.

For ordinary curves it enumerates conductors compatible with Frobenius. Special
j-invariants, a fundamental Frobenius discriminant, or the norm of a constructed
non-scalar prime-degree map can certify the full ring. A conductor-one case
also constructs `omega = pi - [c]`; its restriction to rational points is
`[1-c]`. This distinguishes geometric degree from the scalar action observed
on one subgroup. Supersingular curves keep their quaternion-ring status.

`--count-ops` verifies signed lattice decomposition and complete joint scalar
multiplication against binary and width-3/4 wNAF arithmetic for every scalar in
the selected subgroup. It records separate, uncalibrated counters including
map evaluation, per-call tables and group arithmetic, plus lattice setup once
per batch. Integer counters are named algorithm checkpoints, not instruction
counts. This is a scalar stage diagnostic: wall time, end-to-end ECDLP cost and
speedup remain null. No native performance claim is inferred from the archive's
Python pilot timings. The code is variable-time laboratory arithmetic.

`decision_tree` records conclusions supported by the actual evidence. Distinct
subgroup symmetry actions include negation and retain overlap exactly. Cheap
orbit handling, measured rho gain and index-calculus relations remain unresolved.
Non-rational kernels, longer isogeny paths, extension fields and binary fields
are unsearched. Verification and exhaustive scalar checks are separately
bounded by the field cap, outside the kernel-search time limit.

Native JSON uses `endomorphism-search/native-v1`; replay receipts use
`endomorphism-search-replay/v1`. A probe retains exact field, model, subgroup
and generator metadata, full EC1 curve UID and computed ICV1 identity. These
are representation identities, not isomorphism-class identities. Generated
models are marked `not_registered_by_this_tool`; register a model before citing
its slug or publishing it to the curve catalog. The runner does not change the
registry or infer cross-repository equivalence.

The [protocol and frozen evidence](../../research/endomorphism_search_20261005/PROTOCOL.md)
pin the original source ZIP, inputs, correctness acceptance and scope. The
native replay checks all 15 file hashes against a separately pinned manifest,
reproduces all 500 order summaries and all four archived map workloads, then
verifies every accelerated scalar. It does not execute any archived Python or
certify historical timings. Altered source, evidence or manifest bytes fail
closed without emitting a successful result.

Mathematical references: [ordinary endomorphism orders](https://math.mit.edu/classes/18.783/2023/LectureNotes13.pdf),
[norm and trace](https://math.mit.edu/classes/18.783/2013/LectureNotes14.pdf),
[Bisson–Sutherland ring certification](https://arxiv.org/abs/0902.4670),
[Sage Vélú reference fixture](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_curve_isogeny.html),
and [GLV decomposition](https://arxiv.org/abs/1106.5149).
