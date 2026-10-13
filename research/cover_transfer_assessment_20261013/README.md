---
title: "Exact curve-cover transfer assessment"
subtitle: "Current-catalog crosswalk for two explicit hyperelliptic covers"
date: "2026-10-13"
geometry: margin=0.8in
fontsize: 10pt
---

# Outcome

Two exact same-field hyperelliptic covers are constructively and
reproducibly identified. Only one has any implemented Jacobian transfer, and
that transfer is deliberately bounded.

- The binary degree-three map over \(\mathbf F_{2^7}\) has verified curve
  geometry, checked one-infinity genus-three Jacobian arithmetic, pullback of
  rational elliptic point classes, pushforward of rational source-point
  divisors and the implementation's own typed fibers, and a narrow
  order-29 subgroup certificate.
- Its arbitrary-class \(\operatorname{Jac}(H)\to E\) norm is **not**
  implemented. The public negative-control path reports it as unsupported.
- The prime degree-two map over \(\mathbf F_{827}\) has verified genus-two
  cover geometry. Its source has two rational infinity points, while the
  current prime Jacobian accepts only odd-degree, one-infinity models.
  Prime-cover Jacobian arithmetic and induced maps therefore remain absent.
- No measured computational advantage, complete Jacobian relation solver,
  linear algebra, or target decomposition is claimed.

The canonical machine assessment is
[`docs/curves/cover-transfer-assessments.json`](../../docs/curves/cover-transfer-assessments.json),
validated by
[`cover-transfer-assessment.schema.json`](../../docs/curves/cover-transfer-assessment.schema.json).
Its SHA-256 is:

```text
edd6101f4a9339d14ac02c259f7c323a718da6dee7918097f1053f016b4db9a8
```

The semantic verifier independently pins this complete digest, the catalog
base commit, and the current implementation and contract sources used below.

All claim states use exactly `unknown`, `supported`, `refuted`, and
`not_applicable`. `Refuted` always names a precise implementation or evidence
claim; it is never used to infer mathematical nonexistence from a failed or
missing construction.

# Why this is a continuation, not a rewrite

The [October 6 infrastructure audit](../hyperelliptic_cover_infrastructure_20261006/README.md)
is an immutable, internally consistent snapshot. It pins 121 cover rows, two
rich evidence records, exact source revisions, and the primary-source theorem
audit. Rewriting its counts would destroy provenance. This continuation
instead joins those same two exact records to the current HC1/CV1 graph and
pins the current catalog inputs.

The pinned binary transfer, binary Jacobian, and prime curve kernels still
match that snapshot. The catalog checker and standalone CLI have since
changed, so the sidecar marks wholesale freshness of the historical replay
`unknown` and uses current graph joins as separate evidence.

The primary-source conclusion remains:

- Xarles, *Hyperelliptic curves covering an elliptic curve twice*,
  [arXiv:1303.4220v1](https://arxiv.org/abs/1303.4220v1), Corollary 3, works
  over an arbitrary field \(K\) of characteristic different from 2 and 3.
- It gives a hyperelliptic \(H/K\) of genus at most 5 and two independent
  \(K\)-defined maps \(H\to E\). For \(j\ne0,1728\), the displayed
  genus-five construction has two degree-three maps, hence separable in the
  stated characteristic. The exceptional \(j\)-cases cite genus-two work.
- This is not a uniform executable construction for every input model, and
  it is not the trivial genus-one interpretation. The standard
  hyperelliptic usage here has genus at least two.

The repository's fixed prime and binary constructors are separate,
machine-replayed families. They establish neither all cases of Xarles's
two-map theorem nor a solver.

# Current catalog snapshot

Assessment base commit:

```text
6600c7c6c0cfbb6a9e4e596bf8c5cfa292fb269e
```

| Artifact | SHA-256 | Audited content |
|---|---|---|
| `registry.json` | `9e420a7f...19e872` | 2,091 ICV1 cards |
| `covers.json` | `baac4c80...928324` | 2,091 verified geometry rows |
| `cover-links.yaml` | `a662b445...8e145` | 2,091 maps, 2,089 unique covers, 664 EC1 representations |
| `jacobian-transfers.json` | `aec63380...e6ca40` | stale 321-curve formula screen; not imported as certificates |
| prior `evidence.json` | `b5860bea...9dff6` | two detailed exact records |

Full artifact digests:

```text
registry.json            9e420a7f2e69a1d72a2590899e6bffb8b3d45eff0f197c5c898ee66c1219e872
covers.json              baac4c80f89e3cd9c8ef14421b600cce98c9a642b54ed1f182e8dcc8a2928324
cover-links.yaml         a662b4458eab6d6c10fc7a51b74c0d2c7f6de3a9f817ea981ee898ada8f8e145
jacobian-transfers.json  aec6338011158c8c1aaddb829d4e90fb96eebc43cd4027d1067979bba3e6ca40
prior evidence.json      b5860beae57838245146c5183aa1a96a5c2bf78f2f1b7c669fb401cc8439dff6
```

Every `covers.json` row is `verified_over_declared_field`, every row still
has `subgroup_transfer: not_tested`, and every `dlp_advantage` is null. The
content-addressed graph likewise has 664 `not_tested` subgroup-transport
entries. This assessment crosswalks two maps and two EC1 representations;
the other 2,089 map rows and 662 EC1 representations remain unassessed for
transfer purposes.

The separate `jacobian-transfers.json` file pins registry digest
`f8297bb8...`, contains only 321 curves and 1,156 formula routes, and omits
1,770 current slugs. Its conditional labels and modeled formulas are useful
discovery leads, not cover-map, kernel, image-order, or recovery
certificates.

Browser impact was reviewed. This change adds no curve, map, solver result,
or performance datum and deliberately leaves the generated lab browser and
its geometry-only cards unchanged. The curves index links the new sidecar and
report; a future browser transfer badge must consume these four-state claims
rather than infer transfer from `verified_over_declared_field`.

# Typed correspondence graph

![](typed_correspondence.png)

*Figure 1: Typed correspondence graph.*

Editable source:
[`typed_correspondence.dot`](typed_correspondence.dot). Vector rendering:
[`typed_correspondence.svg`](typed_correspondence.svg).

Green edges are supported only on their labels' exact domains. The binary
pullback route is orange because its point-class constructor is supported
while an independent additivity oracle remains unknown. The red binary norm
edge refutes the claim that the current API handles arbitrary Jacobian
classes. The red prime arithmetic edge refutes only compatibility with the
existing one-infinity implementation; it says nothing against mathematical
existence of a balanced Jacobian representation.

# Exact records

| Item | Binary record | Prime record |
|---|---|---|
| Field | \(\mathbf F_2[z]/(z^7+z+1)\) | \(\mathbf F_{827}\) |
| ICV1 | `icv1-f2m7-t13-616700dd` | `icv1-fp10-t5-192cb216` |
| EC1 | `EC1N7Ce0hb6d297a2ca08` | `EC1P10Cfphcf15d87ab724` |
| Target \(E\) | \(y^2+xy=x^3+1\) | \(y^2=x^3+x+15\) |
| Source \(H\) | \(v^2+u^2v=u^7+u\) | \(v^2=u^6+u^2+15\) |
| Genus / infinity | 3 / one rational point | 2 / two rational points |
| Map | \((u,v)\mapsto(u^3,uv+1)\) | \((u,v)\mapsto(u^2,v)\) |
| Degree / field | 3 / same field | 2 / same field |
| Cover ID | `HC1G3h2530bd70db89` | `HC1G2h998582908336` |
| Map ID | `CV1D3h0dc60016b03d` | `CV1D2h5768ce11909c` |

The exact unabridged UIDs, model hashes, source models, target objects,
subgroup generators, obligations, and evidence references are retained in
the machine sidecar.

# Transfer obligations

| ID | Typed claim | Status | Exact boundary or remaining work |
|---|---|---|---|
| T1 | Binary source, target, and degree-three morphism exist over the declared field | `supported` | Current HC1/CV1 objects and native polynomial replay |
| T2 | Binary source smoothness, genus, separability, defining identity, infinity, and exceptional point | `supported` | Fixed-family verifier; see [`cover_transfer.rs`](../../src/binary_ecc/cover_transfer.rs#L195) |
| T3 | Checked binary Jacobian validity, reduction, equality, negation, addition, and scalar multiplication | `supported` | One-infinity model only; implementation in [`hyperelliptic.rs`](../../src/binary_ecc/hyperelliptic.rs#L351) |
| T4 | Binary pullback \(\phi^*\) | `supported` | Rational elliptic point classes only; infinity and exceptional fibers explicit |
| T5 | Independent pullback additivity oracle on the exact \(\mathbf F_{2^7}\) record | `unknown` | Internal consistency is not an independent CAS/divisor oracle |
| T6 | Binary pushforward \(\phi_*\) | `supported` | Single rational source-point divisors and private typed pullback witnesses only |
| T7 | Arbitrary reduced-class binary norm | `refuted` | Public API returns a typed unsupported result; see [`cover_transfer.rs`](../../src/binary_ecc/cover_transfer.rs#L707) |
| T8 | \(\phi_*\phi^*=[3]\) | `supported` | Typed rational point-class domain, exhaustively replayed across ordinary F8 models |
| T9 | Exact order-29 kernel, image order, and degree-inverse recovery | `supported` | Only `EC1N7Ce0hb6d297a2ca08`; not a catalog-wide result |
| T10 | Prime source, target, and degree-two morphism exist over \(\mathbf F_{827}\) | `supported` | Current HC1/CV1 objects and native geometry replay |
| T11 | Current prime one-infinity arithmetic accepts the even sextic source | `refuted` | It requires degree \(2g+1\); this source has degree \(2g+2\) and two infinities; see [`curve.rs`](../../src/prime_hyperelliptic/curve.rs#L7) |
| T12 | Prime pullback, norm, and \(\phi_*\phi^*=[2]\) | `unknown` | Balanced divisor representation and induced maps are not implemented |
| T13 | Order-823 subgroup preservation for the prime record | `unknown` | Exact EC1 identity exists; kernel, image, and recovery evidence does not |
| T14 | The remaining current cover-map rows preserve their catalog subgroups | `unknown` | Absence of assessment is not refutation |
| T15 | Any cover here yields a computational advantage | `unknown` | No full pipeline or matched end-to-end measurement |
| T16 | A failed bounded constructor proves nonexistence | `refuted` | Structured unsupported/search failure never proves no cover exists |

For a degree-\(d\) finite separable morphism, the expected identity is
\(\phi_*\phi^*=[d]\) on \(\operatorname{Jac}(E)\). Even a verified
identity does not by itself prove preservation of a selected subgroup. The
kernel intersection, exact image order, and recovery route must be bound to
the full EC1 object. The order-29 fixture does this only for its narrow typed
domain.

# Cost and benefit assessment

No calibrated benchmark was performed. Correctness-test elapsed time is not
an algorithmic cost, and a symbolic formula is not a measurement.

| Phase | Value | Status | Reason |
|---|---:|---|---|
| Cover construction | null | `unknown` | executable but uncalibrated |
| Failed searches | null | `unknown` | no catalog-wide search ledger |
| Duplicate candidates | null | `unknown` | no deduplication ledger |
| Geometry verification | null | `unknown` | replayable but uncalibrated |
| Jacobian arithmetic | null | `unknown` | correctness receipts only |
| Pullback / pushforward | null | `unknown` | bounded implementation, no measured conversion cost |
| Subgroup recovery | null | `unknown` | order-29 correctness witness only |
| Reusable preprocessing | null | `unknown` | no amortization boundary |
| Peak memory | null | `unknown` | not measured |
| Total runtime | null | `unknown` | no end-to-end transfer or solver run |
| Matched generic baseline | null | `unknown` | not run |
| Jacobian relations / linear algebra / target decomposition | null | `not_applicable` | separate solver was not implemented or run |
| End-to-end speedup | null | `unknown` | no complete numerator or denominator |

Benefit hypothesis: a cover helps only if a selected EC1 subgroup embeds in
an operationally usable source-Jacobian subgroup and the complete destination
algorithm plus both conversions is cheaper than the matched generic baseline.
This assessment neither confirms nor measures that hypothesis.

Falsification condition: any missing general norm, failed kernel/image
certificate, failed independent relation replay, or total measured cost at or
above the matched baseline rejects the claimed end-to-end benefit even when
the curve morphism itself is correct.

# Controls

| Control | Status | Evidence boundary |
|---|---|---|
| ICV1, EC1, HC1, and CV1 identity replay | `supported` | exact current hashes and UIDs |
| Cover geometry replay | `supported` | both fixed families |
| Binary arithmetic internal exhaustive control | `supported` | all 486 reduced F8 fixture classes and shared-support cases |
| Independent CAS/divisor-class operation oracle | `unknown` | not run |
| Transported-data composition | `supported` | binary typed rational point classes only |
| Coordinate-isomorphism commutation | `unknown` | not run |
| Exact subgroup kernel/image/recovery | `supported` | binary order-29 record only |
| Prime order-823 transfer | `unknown` | no induced maps |
| Malformed/tampered/unsupported negative controls | `supported` | structured invalid or unsupported results |
| Cross-repository source parity | `supported` | PR #377 merge content, pinned manifest, and native replay; not fresh exact-tip CI |
| Paired performance control | `not_applicable` | no benefit experiment occurred |

The mirrored implementation in `aburan28/cryptanalysis` is merged through
[PR #377](https://github.com/aburan28/cryptanalysis/pull/377). The
[replay receipt](release_replay.txt) pins merge `eadd838a...`, the parity
manifest digest, the no-diff path set, and a successful native parity replay
against this checkout. Its prior
`fb-refs` content defect was repaired by
[PR #352](https://github.com/aburan28/cryptanalysis/pull/352); a descendant
validator run recorded zero problems. Current exact-tip workflow freshness
is still absent because the repository's Actions allowlist blocks setup and
the workflow is manually disabled. This is reported as a CI-delivery gap,
not as a mathematical failure.

# Standalone release replay

Crypto [PR #1492](https://github.com/aburan28/crypto/pull/1492) merged as
`7636a3c8e12c210c4b7ab611c791d1979a177d9e`. Its immediate release run 999
failed on Windows and did not create tag `v0.1.0-main.999`. A later descendant
[release run 1011](https://github.com/aburan28/crypto/actions/runs/37617671402)
completed all eight jobs and published Linux GNU, Linux musl, macOS, and
Windows Rust bundles containing `hyperelliptic-cover`.

The later descendant
[`v0.1.0-main.1371`](https://github.com/aburan28/crypto/releases/tag/v0.1.0-main.1371)
was independently replayed for this assessment. The exact commands, inputs,
stdout hashes, and catalog comparison are retained in the
[`release_replay.txt`](release_replay.txt) receipt:

- tag commit `bbd330a9d8e9bb59346bb9f5f45416851d2af3db` contains the cover merge;
- the exact asset is
  `crypto-v0.1.0-main.1371-x86_64-unknown-linux-gnu.tar.gz`;
- its SHA-256 is
  `77e69e60a68ada60c7b5fe13503e326f50fce4594b96340eca402041e60dcb2d`,
  matching the published `SHA256SUMS`;
- the archive contains `hyperelliptic-cover`, `registry.json`, `covers.json`,
  and `cover-links.yaml`;
- direct `construct`, `verify`, and binary-Jacobian `identity` commands
  returned status `ok`; and
- `catalog-export` completed without `cargo run`. Its mathematical rows
  matched the bundled catalog; the exported file intentionally changed the
  producer provenance fields.

The latest current-main release attempt has unrelated native/isogeny job
failures, while all three Rust bundle jobs (Linux, macOS, Windows) passed.
That newer failure does not erase the published descendant bundles above.

# Exploration boundary

This work did not search production targets, recover keys, or run a discrete
logarithm campaign. It assessed supplied constructions and current repository
evidence only.

It does not:

- search all cover families, genera, degrees, or characteristics;
- assess the other 2,089 current map rows beyond their existing geometry
  status;
- turn missing evidence or failed bounded search into nonexistence;
- import stale formula screens as transfer certificates;
- claim a general binary norm or any prime induced map;
- claim subgroup preservation without exact EC1 evidence; or
- implement the full Jacobian index-calculus solver.

# Open obligations

1. Implement and independently verify an arbitrary-class binary norm.
2. Replay pullback additivity and full composition through a genuinely
   independent divisor-class implementation.
3. Add balanced/two-infinity prime Jacobian arithmetic and prime-family
   induced maps.
4. Bind kernel, image order, and recovery certificates to every EC1 subgroup
   for which transfer is claimed.
5. Add coordinate-isomorphism controls.
6. Refresh the formula-screen inventory against the current registry without
   presenting modeled routes as certificates or measurements.
7. Assess the remaining 2,089 maps and 662 EC1 representations explicitly.
8. Keep the full solver separate. Its acceptance gate must require actual
   Jacobian relations, independently verified linear algebra and target
   decomposition, checked outputs, explicit parameter limits, and complete
   preprocessing, failure, duplicate, memory, and total-runtime accounting.

# Reproduction

Structural and identity checks:

```text
jq -e . docs/curves/cover-transfer-assessment.schema.json
jq -e . docs/curves/cover-transfer-assessments.json
jsonschema \
  -i docs/curves/cover-transfer-assessments.json \
  docs/curves/cover-transfer-assessment.schema.json
sha256sum docs/curves/registry.json docs/curves/covers.json \
  docs/curves/cover-links.yaml \
  docs/curves/jacobian-transfers.json
node research/cover_transfer_assessment_20261013/verify_assessment.mjs
```

Visual and report derivation:

```text
dot -Tsvg typed_correspondence.dot -o typed_correspondence.svg
dot -Tpng typed_correspondence.dot -o typed_correspondence.png
pandoc README.md \
  --from=markdown+tex_math_single_backslash \
  --pdf-engine=pdflatex \
  -o cover_transfer_assessment_20261013.pdf
pdfinfo cover_transfer_assessment_20261013.pdf
```

Relevant sources are the [binary transfer](../../src/binary_ecc/cover_transfer.rs),
[binary Jacobian](../../src/binary_ecc/hyperelliptic.rs),
[prime Jacobian](../../src/prime_hyperelliptic/curve.rs), and
[catalog checker](../../src/bin/curve_cover_check/checker.rs), plus the
[CLI documentation](../../docs/curves/HYPERELLIPTIC_INFRASTRUCTURE.md) and
[prior detailed evidence](../hyperelliptic_cover_infrastructure_20261006/evidence.json).
