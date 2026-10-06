# Artifact index — P-192 weak-object study

Status: **report scaffold; experiment artifacts pending**

Date: 2026-10-06

This index separates immutable inputs, derived communication artifacts, and
future run evidence. A path marked `PENDING` is a required slot, not evidence
that a run occurred. No missing hash is to be inferred from a filename.

## Frozen protocol sources

| record | role | source repository path | source snapshot SHA-256 | archive status |
|:--|:--|:--|:--|:--|
| `EXP-SCURVE-29040c` | singular-companion recovery protocol | `experiments/EXP-SCURVE-29040c/specification.yaml` in `crypto-autoresearcher` | `5de49bcf72cc676290982314046c6579b10473372aed822c6955813660e9030e` | coordinator archive pending |
| `EXP-SCURVE-647ade` | exact CM relation protocol | `experiments/EXP-SCURVE-647ade/specification.yaml` in `crypto-autoresearcher` | `0fe122a0c6be3e2118b4d9c276a1085285c7e6a2cee38d15bbeed8099fa17145` | coordinator archive pending |
| `DEC-20261006-af8214` | approval decision | `ledger/decisions/DEC-20261006-af8214.yaml` in `crypto-autoresearcher` | PENDING archive manifest | coordinator archive pending |

The hashes bind the exact protocol bytes inspected while this scaffold was
written. Until those records are committed and archived, the hashes are local
snapshot provenance, not a durable source URL.

## Existing immutable or committed provenance

| artifact | immutable identifier or SHA-256 | evidentiary role |
|:--|:--|:--|
| P-192 representation | `EC1P192Cp192h5531c4a08bdb`; `urn:ec-record:1:sha256:5531c4a08bdb64b6e86a6e30e9a08aa57edef7af15ac5f6d4d2a83a53bf2f646` | exact valid curve/subgroup/generator identity |
| `docs/curves/registry.json` | `25083c9cf2ac96f8b117fa32492a3dc86e0ea42500ddcaf086a3b7b5d7948cb3` at base revision `0dc6b7f4a03253e7c7641f99700f7f423f8d07e0` | registered ICV1/EC1 crosswalk |
| `src/ecc/curve_zoo.rs` | `cc1377410a05a7e7c5806317043e8cfef50473303ed25fd73dff40b8610d0ee9` at the same base revision | repository P-192 parameters |
| `AGENTS.md` | `c8b8eb367d1153421fa440752e41982bf3e59605b467dbd376ef4a60c7202c23` at the same base revision | evidence, accounting, identity, and isolation rules |
| `docs/curves/ic/curve-links/README.md` | `3c00e1a03aca118b118ea2534669f52936cfb0bb2624bbecfe72120e26a2ac9c` at the same base revision | typed map/transport evidence requirements |
| prior P-192 walk | run `p192-71c04205135fcf8f`; store marker `s3://crypto-autoresearcher/isogeny-walk/runs/p192-71c04205135fcf8f/complete.json` | prior valid-isogeny-class coverage; not a weakness result |
| prior walk `README.md` | `6e807de8299dc9631e4cd6e48fdcc22829329c76fca76f36ff30e9e57c136f62` | prior conclusion and limitations |
| prior walk `STORE.json` | `58f14a65a28d6a3f3ba6f462a3633ed614f0b29bd9e1267f101ec57167541fe5` | archived object keys and digests |
| prior walk `walk.json` | `a706473796b58e03cf2a0bef46d828b2c0de686302dd73ebbadf4babed0323bc` | prior run/class metadata |
| singular model preimage | `singular-model/v1 sha256:305129485b8a281f5cd73a7a0ec4d82a7f84f432847c1b9d5eb1734983db72b3` | typed canonical preimage; deliberately not EC1/ICV1 |

## Derived communication artifacts

| path | format | status | source of truth |
|:--|:--|:--|:--|
| `REPORT.md` | Markdown | present; results remain PENDING | cited protocol and evidence records |
| `REPORT.pdf` | PDF | generated after source validation | `REPORT.md` plus current SVGs |
| `object-attack-map.dot` | Graphviz source | present | report claim boundaries |
| `object-attack-map.svg` | vector rendering | generated | `object-attack-map.dot` |
| `object-attack-map.pdf` | PDF-build vector rendering | generated | `object-attack-map.dot` |
| `cm-degree-frontier.csv` | exact/derived plot data | present; no empirical data | `EXP-SCURVE-647ade` |
| `render-cm-degree-frontier.mjs` | deterministic SVG renderer | present | CSV only |
| `cm-degree-frontier.svg` | vector rendering | generated | CSV and renderer |
| `cm-degree-frontier.pdf` | PDF-build vector rendering | generated | current SVG |
| `pdf-images.lua` | Pandoc image selector | present | chooses PDF siblings only for XeLaTeX |

## Required recovery-run evidence (`EXP-SCURVE-29040c`)

Fresh directory: `runs/singular-recovery/<RUN-ID>/` — **PENDING**.

| required artifact | status | SHA-256 | verification note |
|:--|:--|:--|:--|
| `run.yaml` | PENDING | PENDING | fresh RUN ID, exact commit/argv/environment/timestamps |
| `recovery.json` | PENDING | PENDING | transcript, residue stages, BSGS certificate, final replay |
| `order-and-point-certificates.json` | PENDING | PENDING | exact order and smooth-point checks |
| `controls.json` | PENDING | PENDING | safe-path rejection and mutation controls |
| `isolation.jsonl` | PENDING | PENDING | isolation/contended status |
| `phase-metrics.json` | PENDING | PENDING | counted native operations; time secondary |
| `stdout.log` / `stderr.log` | PENDING | PENDING | retained verbatim |
| `SHA256SUMS` | PENDING | PENDING | covers every run artifact, excluding itself |
| independent replay receipt | PENDING | PENDING | separate validator and final `[d]G=Q` replay |

## Required CM-run evidence (`EXP-SCURVE-647ade`)

Fresh directory: `runs/cm-relation/<RUN-ID>/` — **PENDING**.

| required artifact | status | SHA-256 | verification note |
|:--|:--|:--|:--|
| `run.yaml` | PENDING | PENDING | fresh RUN ID, exact commit/argv/environment/timestamps |
| `p192-pocklington.json` | PENDING | PENDING | recursive certificate for `C` |
| `order.json` | PENDING | PENDING | tuple, order, Hasse, discriminant/conductor claims |
| `generators.json` | PENDING | PENDING | roots and orientation fixtures through `ell=113` |
| `exact-boundary-manifest.json` | PENDING | PENDING | state counts and boundary digest |
| `weights.json` | PENDING | PENDING | if unsupported, preserve explicit unsupported receipt |
| `relations.json` | PENDING | PENDING | exact census and classifications |
| `verification.json` | PENDING | PENDING | independent forms/HNF/map replay status |
| `operation-counts.json` / `resource-metrics.json` | PENDING | PENDING | primary counted units; timing only if isolated |
| `stdout.log` / `stderr.log` | PENDING | PENDING | retained verbatim |
| `SHA256SUMS` | PENDING | PENDING | covers every run artifact, excluding itself |
| independent replay receipt | PENDING | PENDING | second enumeration order and certificate replay |

## Publication closeout checklist

- [ ] Replace every observed-value `PENDING` in `REPORT.md` from hashed run files only.
- [ ] Preserve unsupported map/calibration output as an explicit gap; do not
      synthesize HOT/COLD costs.
- [ ] Regenerate both SVGs and `REPORT.pdf`; inspect all pages.
- [ ] Record report/figure hashes without editing immutable run receipts.
- [ ] Link the crypto PR and the archived coordinator decision/run records.
- [ ] Keep the singular lane classified as an implementation-validation test,
      not as an isogeny or a weak elliptic curve.
