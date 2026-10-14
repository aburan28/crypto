---
title: "ECDLP audit execution checkpoint"
subtitle: "Native control plane, cover-transfer boundary, and draft workstream contracts"
date: "2026-10-14"
subject: "authorized public/synthetic audit; no live-target execution; no completion or solver claim"
geometry: margin=0.8in
fontsize: 10pt
---

# Result

This checkpoint moves the eight-workstream program from planning into a first
locally reproduced implementation slice. It adds a native, fail-closed validator for
the non-authoritative planning record; audits the present cover and Jacobian
surface; records a canonical-JSON identity conflict; and supplies draft
contracts and synthetic fixtures for evidence, route, corpus, accounting, and
governance work.

No scientific acceptance gate is promoted by this checkpoint. The defensible
mathematical result remains narrower: explicit cover construction has
source-level evidence for all 2,091 catalog rows, but catalog-wide Jacobian
transfer, exact subgroup preservation, an independent divisor-class oracle,
and an ECDLP advantage are not established.

![Current evidence boundary. Green is bounded source evidence, yellow is
partial or draft evidence, red is an open mathematical or release obligation,
and blue is a non-authoritative control-plane implementation.](../../../docs/ecdlp-audit/execution-boundary.png)

# Authorization and claim boundary

Work in this checkpoint is limited to repository inspection, public source
evidence, and synthetic or known-answer fixtures. It does not use secret key
material, probe a live target, recover a production key, or authorize another
party to do so. “Audit of illicit ECDLP” is treated precisely as an authorized
audit of publicly documented unauthorized recovery routes.

A passing schema or planning check establishes only its stated contract. A
failed bounded search is not proof of nonexistence; an explicit cover is not a
subgroup certificate; a relation collector is not automatically a complete
solver; and an elliptic norm projection is not a Jacobian index-calculus
solver.

# Local execution in this checkpoint

The repository now exposes:

```text
ecdlp-audit verify-plan \
  --input docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json \
  --repo-root .
```

The native Rust validator executes every keyword used by its embedded,
self-contained Draft 2020-12 schema through a fail-closed interpreter with
format checks; unknown keywords, patterns, and remote references are rejected.
It then applies strict duplicate-key-safe deserialization and semantic checks. It
binds the full planning record to its canonical identity, enforces the exact
WS1–WS8 and A–E skeletons, recomputes gate and obligation metrics, rejects
invalid promotion, unbound partial artifacts, provenance substitution, and
missing or optional obligations, detects active blocking cycles, confines the
input and referenced paths to the repository, verifies size-bounded pinned Git
blobs and current-file hashes, and requires every positive mock completion
counter and flag to remain zero or false. Its tests cover those boundaries plus
schema-invalid-but-typed records, unknown fields, malformed JSON, report/input
mismatch, and raw-UTF-8 identity behavior.

On Rust 1.90.0, a targeted local component run produced:

```text
cargo test --locked --offline --bin ecdlp-audit
45 passed; 0 failed

cargo run --locked --offline --bin ecdlp-audit -- verify-plan \
  --input docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json --repo-root .
ECDLP audit plan verification: PASS
workstreams=8 required_gates=20 open_required=20
authoritative=false promotion_allowed=false subjects=0 assessments=0

cargo test --locked --offline --lib \
  binary_ecc::cover_transfer::tests
7 passed; 0 failed; 4,126 filtered out; 104.68s
```

The repository-wide release-library gate was also exercised. An initial run in
the restricted sandbox produced:

```text
cargo test --locked --offline --release --lib -j1
4009 passed; 10 failed; 113 ignored
```

This was not a clean pass. Nine failures were loopback-socket tests blocked with
`EPERM`; all 37 `cryptanalysis::pollard_collab` tests passed when rerun with
network permission. The remaining failure was the sealed ECBench legacy-record
decimal-semantics test. Adding `jsonschema` had feature-unified `serde_json`
with `float_roundtrip`, altering that legacy decimal parse. The repair removed that
dependency and uses the audit binary's fail-closed embedded-schema interpreter,
avoiding the crate-wide parser change. A full post-repair replay in a
network-enabled sandbox, still with Cargo offline, produced:

```text
cargo test --locked --offline --release --lib -j1
test result: ok. 4019 passed; 0 failed; 113 ignored; 0 measured; 0 filtered out; finished in 646.70s
```

This is a repository-wide release-library pass for the stated tree and
environment. It is not a packaged-archive replay, cross-platform CI result, or
E5 reproduction receipt.

A fresh release-mode `ecdlp-audit` executable was then built from the final
source with Rust 1.90.0 and the local GNU bfd linker fallback. Its SHA-256 was
`f861e7fc73b9be9ca7d01b5ce750607d77e6bad22b0a84ad99cdedb18b0065e3`;
`--help`, `verify-plan`, and `verify-fixed-f8-order` all exited zero. This is a
local executable replay, not a published archive or catalog-wide package
receipt, so the WS1 packaging gates remain open.

The first build attempt exhausted the runner's temporary filesystem while
compiling the pre-existing library. The retry used a fresh build directory on
the workspace filesystem and passed. That capacity event is an environment
receipt, not a product or test failure. Clippy completed for the binary; it
reported existing library warnings plus the inherited unknown-lint warning and
no diagnostic in `src/bin/ecdlp_audit.rs` or its fixed-F8 verifier module.

These runs began from base revision
`1e447956d3616092a223264f1c6c31eff8ac9979` in a modified worktree. No immutable
raw command log or environment image is attached, so these local results are
not an E5 reproduction receipt and cannot promote a scientific claim.

# WS1 — release and packaged replay

- **Current source ledger — evidence present, but not package evidence.** The
  generator's byte-exact check passed in a clean checkout. Semantic totals
  reached 2,091 rows: 160 binary, 1,931 odd, 664 EC1 bindings, and 2,089 unique
  cover/Jacobian records. Three historical `git show` subprocesses were
  blocked by the audit sandbox; manual file hashes matched the pinned records.
- **Published standalone replay — open.** The retained release receipt binds
  release `v0.1.0-main.1371` / tag `bbd330a…` and a Linux glibc artifact, but
  bundles and replays 592 rows rather than the current 2,091-row denominator.
- **Artifact identity — open.** No Rust `BUILD.json`, embedded build identity,
  SBOM, dependency manifest, or full catalog/export replay receipt was found.
- **Package/documentation parity — open.** The published/base configuration
  omitted the main musl target while its documentation still listed it; the
  inspected archive contained six executables where documentation said five,
  and bundled catalogs where documentation said otherwise. This worktree now
  configures and documents seven native Rust executables, including
  `ecdlp-audit`, and stages non-authoritative audit examples, but no
  cross-platform CI run or new archive receipt has verified that correction.
- **Final CI and review — open.** The head checks observed during the audit
  reported failures in certificate, transfer-assessment, and Lean jobs. This
  checkpoint does not pin a raw CI receipt or check URL, no approving review
  was attached, and no polling was performed.

The next admissible WS1 run is an exact-artifact replay: bind the final source,
toolchain, lockfile, dependencies, SBOM, configuration, catalogs, archive and
binary hashes; extract the archive in a clean directory; execute its packaged
`hyperelliptic-cover catalog-export` against the bundled registry; retain all
2,091 rows including failures; compare canonical mathematical records; and run
malformed, singular, reducible, and explicitly unsupported controls.

# WS2 — exact supported mathematics

- **Binary covers.** Over polynomial-basis `GF(2^m)`,
  `E: y^2+xy=x^3+a*x^2+b` maps to genus-3
  `H: v^2+u^2*v=u^7+a*u^4+d*u`, `d^2=b`, by
  `(u,v) -> (u^3,u*v+d)`. The degree-3 map is separable and the source has one
  rational infinity. This covers 160 catalog rows. The exact F8 fixture
  described below now has an all-class norm; a catalog-general arbitrary-class
  norm is still missing.
- **Odd covers.** In characteristic greater than 3, short-Weierstrass curves
  map from genus-2 even sextics `v^2=F(u^2+c)` by
  `(u,v) -> (u^2+c,v)`. There are 1,931 rows: 577 direct-prime, 30 normalized,
  and 1,324 extension rows. Those sources have two infinities; the current odd
  Jacobian backend rejects even degree, and no induced-map implementation or
  certificate exists for these 1,931 odd catalog covers.
- **Binary Jacobian.** The implementation checks characteristic-two,
  one-infinity Mumford validity, canonicalization, equality, negation, Cantor
  addition/reduction, and scalar multiplication. A fixed F8 fixture exhausts
  486 classes and exercises the arithmetic's pair cases, including shared
  support and the ramified origin. A separate native Rust module counts
  the source over the first three extensions and certifies the cardinality
  `#Jac(H)(GF(8))=486`; it does not replay divisor arithmetic. All-triples
  associativity and an independent arithmetic oracle remain open.
- **Odd Jacobian.** Genuine odd-prime monic odd-degree `2g+1` one-infinity
  arithmetic exists. A genus-2 F3 fixture exhausts 24 classes, pairs, and
  triples. It accepts zero of the 1,931 catalog even sextics and has no generic
  `GF(p^k)` backend.
- **Transfer.** Binary rational-point pushforward, pullback of `[P]-[O]`,
  infinity/two-torsion fibers, and typed point-class pull-push certificates
  exist. For the exact polynomial-basis F8 model `a=1`, `b=3`, a cached
  486-class presentation now maps every valid reduced class. Native tests
  enumerate degree counts `1, 5, 56, 424`, check all 236,196 ordered
  homomorphism pairs, six fibers of size 81, all rational source points, all
  target-point pull-push compositions, explicit shared/nonsplit/exceptional
  classes, and malformed inputs. The implementation-separated order certificate confirms
  the ambient cardinality but does not validate the presentation, norm, or
  homomorphism, so this remains same-implementation transfer evidence and does
  not generalize to other fields or models. The existing conditional Lean
  group lemma does not formalize the fixture's field/curve model, `N1`–`N3`
  and zeta bridge, divisor presentation, or norm homomorphism; those proof
  obligations remain open and no Lean receipt is claimed for them. Odd transfer
  and a catalog-general binary norm remain open. The retained cryptanalysis
  parity receipt binds a pre-norm `cover_transfer.rs` digest, so a current
  mirror or integration of the fixed-F8 path is not established either.
  The target group has order 6, so `gcd(3,6)=3`; this fixture exercises norm
  mechanics, but the degree-coprime theorem does not make it an injective
  full-target-group ECDLP transfer.
- **Subgroups.** Results are conditional when an exact subgroup order `r` is
  supplied and `gcd(3,r)=1`; narrow order-29 overlay evidence exists.
  Factor/primality chains, exact target and image orders, kernel intersection,
  and recovery for all 664 bindings are not established.
- **Independent reference.** The pinned [fixed-F8 order certificate](../../../docs/curves/fixed-f8-genus3-order-certificate.json)
  and its [implementation-separated replay](../../../src/bin/ecdlp_audit/fixed_f8_order.rs)
  use direct finite-field enumeration, Newton identities, and the genus-three
  functional equation without importing production curve arithmetic. No independent
  divisor-class/transfer oracle receipt is attached. A separately implemented
  brute-force arithmetic/transfer oracle or a suitable pinned CAS/container is
  still required.
- **Solver boundary.** `hyperelliptic_index_calculus.rs` is a real bounded
  odd-prime, one-infinity relation/linear-algebra/recovery library for tiny
  enumerable domains. It does not consume the catalog's binary/even-sextic
  covers and lacks catalog-wide accounting and independent output evidence.

The correct sequencing is: freeze an independent small-field oracle protocol
and compare the implemented 486-class F8 norm to it; implement a balanced
two-infinity representation over one small odd field and compare it to the
oracle; bind proof-grade subgroup certificates; then generalize fields and
catalog coverage. A failed search at any rung records that bounded failure
only.

# WS3–WS8 control-plane progress

All artifacts listed here are deliberately non-authoritative. Where a draft
matches a named deliverable, its existence can change that deliverable's
progress from “not attempted” to “partial”; it does not satisfy or promote the
associated acceptance gates.

- **WS3 ledger:** native plan validator and canonicalization audit. Still open:
  authoritative unit/assessment schemas, frozen denominator, additive ledger
  generator, cross-language identity vectors, semantic verifier, and reviewed
  migration.
- **WS4 evidence:** E0–E7 profile schema, five route profiles, and transition
  fixtures. Still open: governed registry, native transition engine,
  cross-record freshness/history checks, and independent review.
- **WS5 routes:** five-lane route-registry and mock coverage contracts. Still
  open: frozen cited route population, authoritative denominator, full
  route-by-asset classification, and review.
- **WS6 corpus:** search protocol, strict source/claim schema, pointer-only
  seed, and fail-closed fixtures. Still open: frozen reviewed protocol,
  executed searches, acquired primary evidence, citation passes, comparator,
  and gap reproduction.
- **WS7 accounting:** cross-route cost schema, four valid synthetic route
  receipts, and four structural negative controls. Still open: native
  recomputation/claim checker, identity reconciliation, migration adapters,
  and independently replayed real receipts.
- **WS8 governance:** policy, role/ownership matrix, strict schema, and six
  workflow drills. Still open: approved policy, authenticated identities,
  protected decisions, append-only storage, cross-platform release/package
  verification, and external decisions.

Pinned AJV CLI 5 / Draft 2020-12 replay accepts the plan and every designated
valid WS4–WS8 draft record. It rejects all designated WS5, WS6, and WS7
structural negative fixtures. WS4 and WS8 rejection/demotion drills are
schema-valid event records by design; their expected decisions still require
the absent native semantic engines. Schema conformance does not establish the
truth, completeness, or authority of any fixture.

The WS3 audit found a concrete identity hazard. The documented EC1 convention
and historical reference emit sorted compact JSON with raw UTF-8, while the
historical ecbench codec escapes non-ASCII. Existing ASCII records conceal the
difference. The proposed new audit codec follows raw UTF-8 and pins a Unicode
digest vector, while existing ecbench identities remain immutable and require
typed migration receipts. WS3-G1 remains open pending a reviewed codec,
duplicate-key-safe ingest, frozen cross-language vectors, and domain-separated
identity schemas.

# Open acceptance boundary

At this checkpoint all 20 required planning gates remain open and all attack
route skeletons remain `unknown`. The scoped WS1 package-replay and WS2
bounded-transfer claims cannot advance until exact artifact replay and an
independently checked bounded Jacobian-transfer fixture exist; those two
conditions are not sufficient for any broader promotion, whose remaining
gates and WS8 approval also stay mandatory. Full Jacobian index calculus
remains a separate deliverable whose acceptance must include actual Jacobian
relations, verified linear algebra and target decomposition, independent
output checks, explicit parameter support, and total accounting for
preprocessing, failed searches, duplicates, memory, and runtime.

# Artifact index

- Human plan: `RESEARCH_COMPREHENSIVE_ECDLP_AUDIT_PLAN.md`
- Machine planning twin: `docs/ecdlp-audit/ECDLP_AUDIT_PLAN.mock.json`
- Native validator: `src/bin/ecdlp_audit.rs`
- Identity audit: `docs/ecdlp-audit/ledger/IDENTITY_CANONICALIZATION_AUDIT.md`
- Draft contracts: `docs/ecdlp-audit/{evidence,routes,corpus,accounting,governance}/`
- Execution-boundary diagram: `docs/ecdlp-audit/execution-boundary.{dot,svg,png}`

The PDF twin of this report is a rendered convenience copy. The Markdown,
machine records, and source revisions are reviewable. The local command output
is summarized above but is not yet retained as an immutable receipt.
