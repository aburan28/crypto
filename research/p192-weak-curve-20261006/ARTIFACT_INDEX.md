# Artifact index — P-192 weak-object study

Status: **independent review complete — overall `BREAKS`; singular scientific
joints `HOLDS`; CM complete-set identity `INCONCLUSIVE`**

Date: 2026-10-06

This index separates immutable inputs, producer evidence, archive metadata,
independent review, and derived communication artifacts. A `HOLDS` subclaim
does not override the independent review's overall `BREAKS` verdict.

## Frozen protocol sources

| record | role | source repository path | source snapshot SHA-256 | archive status |
|:--|:--|:--|:--|:--|
| `EXP-SCURVE-29040c` | singular-companion recovery protocol | `experiments/EXP-SCURVE-29040c/specification.yaml` | `ae25388412acb9fbc8a0387f8d4782fe215a493b6bccb6f1f0fa49f8f1b3fde6` | committed at `3cb96ec70` |
| singular v2 amendment | shifted-center BSGS coverage | `experiments/EXP-SCURVE-29040c/amendments/v1_to_v2.yaml` | `c5be1d42be520cbadec68c76928f579559f8e1ebfb7fa9aa93c58fa0853512a0` | committed at `3cb96ec70` |
| `CORR-20261006-500337` | exact x-only tag encoding correction | `ledger/corrections/CORR-20261006-500337.yaml` | `bad4ab248d8244cf56c3a021fb7e69c441deb13ec930be47b40031a56c053c94` | committed at `3cb96ec70` |
| `EXP-SCURVE-647ade` | exact CM relation protocol | `experiments/EXP-SCURVE-647ade/specification.yaml` | `0fe122a0c6be3e2118b4d9c276a1085285c7e6a2cee38d15bbeed8099fa17145` | committed at `3cb96ec70` |
| CM v2 amendment | exact boundary/strata | `experiments/EXP-SCURVE-647ade/amendments/v1_to_v2.yaml` | `029fcf0ed147af8eb590e770bed9c2a023eadfd3bca8567d0297caa6130d3a33` | committed at `3cb96ec70` |
| `DEC-20261006-af8214` | approval decision | `ledger/decisions/DEC-20261006-af8214.yaml` | `b374c1b25d1bfd27c94d2f616333f4fe116b68784dd32e9539133ac3798d897b` | committed at `3cb96ec70` |
| `CORR-20261006-9ef467` | content-at-commit and run-registration correction | `ledger/corrections/CORR-20261006-9ef467.yaml` | `d453b8b4d33eee6715ea86efb85aadc95544e2f81ef721eedc26c143c541968d` | committed at `2744a4363` |
| `CORR-20261006-51b7d2` | first singular-origin transcription correction; superseded for overbroad scope | `ledger/corrections/CORR-20261006-51b7d2.yaml` | `0fe2a2360afbf03fe145a5566542d36067f9cb3c8da9fa60b23147c0b914f8f1` | committed at `9422ebece`; superseded by `CORR-20261006-cfc18a` |
| `CORR-20261006-cfc18a` | authoritative additive singular-origin overlay | `ledger/corrections/CORR-20261006-cfc18a.yaml` | `675a56868dc0de1aeee69e445a289e208d47882ba06072e66c7ca1163c5ef051` | committed at `4a6830200` |

Snapshot receipt `TASK-20261006-11f0a2` has SHA-256
`37d779d50bfbdedb8288e31410babc1908aefc1e2977b42e7d4369462b3b12c8`
at commit `8fbdb72ba`. It binds 62 governed paths at content commit
`2744a436378489cb9dcb9808712edf6362e51e29`, while retaining the separate
producer-origin commits `8b81f8efa` and `9b2a47bca`.

The immutable receipt spells out a non-resolving full SHA for the first of
those abbreviated commits. `CORR-20261006-cfc18a` overlays the verified full
value `8b81f8efa613766a76eaf46465892f22b37fad5d`; the receipt itself remains
byte-for-byte identical to commit `8fbdb72ba`.

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
| `REPORT.md` | Markdown | validator disposition incorporated | cited protocol, run, and review artifacts |
| `REPORT.pdf` | PDF | regenerated and visually inspected, 10 pages | `REPORT.md` plus current SVGs |
| `object-attack-map.dot` | Graphviz source | present | report claim boundaries |
| `object-attack-map.svg` | vector rendering | generated | `object-attack-map.dot` |
| `object-attack-map.pdf` | PDF-build vector rendering | generated | `object-attack-map.dot` |
| `cm-degree-frontier.csv` | exact/derived plot data | present; no empirical data | `EXP-SCURVE-647ade` |
| `render-cm-degree-frontier.mjs` | deterministic SVG renderer | present | CSV only |
| `cm-degree-frontier.svg` | vector rendering | generated | CSV and renderer |
| `cm-degree-frontier.pdf` | PDF-build vector rendering | generated | current SVG |
| `pdf-images.lua` | Pandoc image selector | present | chooses PDF siblings only for XeLaTeX |

Final communication-artifact hashes after the independent disposition:

- `REPORT.md`: `babc8777503cee828c2251448b4a99c00290b4e14d6cef5d720513c12441bcac`
- `REPORT.pdf`: `060a3aad15f36c282a7dc2aef329a2edc09ccd3554cc8644e9bdf6a7cfba2e3a`
- `object-attack-map.dot`: `3a8407003a29082cee4a74c09c58a0eb129697989635fd1d040dc950a4f5a95f`
- `object-attack-map.svg`: `361d3e053a5a944370ba7b7c4e03ca7c37a54bd2749e34a57f096a263b21a0b6`
- `object-attack-map.pdf`: `74566a8484ee3c0b3e36e3146a8fb7af25d77193041e1fdbd0606dcbc8f39d83`
- `cm-degree-frontier.svg`: `9058cfb51bda4b87f5a11ea168258b268ac8c098990a7d4f0e43fea18be2a4bc`
- `cm-degree-frontier.pdf`: `9537a68d06450aafd55d1e37e04f317668fa5397e569299c55090480fa10a578`

## Recovery-run evidence (`RUN-SCURVE-af5caf`)

Source directory:
`experiments/EXP-SCURVE-29040c/runs/RUN-SCURVE-af5caf/` in
`crypto-autoresearcher`, producer commit `8b81f8efa`.

| artifact | status | SHA-256 | verification note |
|:--|:--|:--|:--|
| `run.yaml` | producer-complete | `1568b333dae66e1801e9f431d32ddaefcbb9e81bd292635ea542855de84d766f` | code commit, argv, environment, timing |
| `recovery.json` | producer-complete | `3f45d23151926a4d38086b1cce916891c13a0866af9a278d173d8a4470f2ab5c` | transcript, residues, BSGS, final replay |
| `bsgs-cover-certificate.json` | producer-complete | `808d79a8af586db6e5bf6cb0091c9dcd62cb38b520611007ea86c9fcfca7fd64` | shifted-negation-v2 interval coverage |
| `order-and-point-certificates.json` | producer-complete | `af55f158b7bf96b4d1f97ee0125123b1bc894d9a41dd58215ec9845a9017c5ef` | exact order and point witnesses |
| `controls.json` | producer-complete | `074a3f3d9b4dd95ea07e3174eee8d3f20085f3bb493a0ad5452fb26fae03feec` | safe-path rejection and mutation controls |
| `phase-metrics.json` / `isolation.jsonl` | producer-complete | `4b06868072e11d65723d76eb56629b2fd4b62d5cbf02cb9bb0c3564af05c42a1` / `fcba9c7dcebffd8e2cfa0211da68b1aa5200b5131b2af3b06480e1b15c38cc1c` | wall/CPU/RSS and uncontended flag |
| logs and `task-report.yaml` | producer-complete | manifest-bound | preflight refusal, run, and verifier logs retained |
| `manifest.sha256` | verified | `0c1485432aaf35f0be94f8f9c3446fc4cb83bfd90c2c44cbbcf6e34acd7873e4` | covers all 14 producer files other than itself |
| `registration.sha256` | verified | `d3271d644c75c2430e3b4a82927242ac9e5cb4cfb031dc95f3ad346d1db66b27` | additive ledger envelope; producer bytes unchanged |
| independent validation | **SCIENTIFIC JOINTS `HOLDS`; OVERALL `BREAKS`** | recorded in validator report | all 72 tags, linkage, BSGS equation, `k`, `d`, replay, and controls rederived; snapshot-origin typo preserved |

## CM-run evidence (`RUN-SCURVE-81c5e4`)

Source directory:
`experiments/EXP-SCURVE-647ade/runs/RUN-SCURVE-81c5e4/` in
`crypto-autoresearcher`, producer commit `9b2a47bca`.

| artifact | status | SHA-256 | verification note |
|:--|:--|:--|:--|
| `run.yaml` | producer-complete, full experiment incomplete | `3ccf11348eebe42ee7021d6ebe6431b115013f2e01c7c593217d0b3db85f99be` | all three attempts and dependency limitation |
| `p192-pocklington.json` | producer-complete | `7db2cb0d691739e7a862d3db99eda303633e82281f9170e6db43125271f102bc` | recursive certificate for `C` |
| `order.json` | producer-complete | `5ff98b29189055141dac92de8560b290c54dc1e10dd5dbe798d3f50fdf6c2a6f` | tuple/order/Hasse/discriminant/conductor |
| `generators.json` | producer-complete | `a213ca6aa2df6760b63b16d38ceab226f8ea9b1ef3e0c5cc75bfb0a5b988b72e` | 13 roots/orientations and algebra controls |
| `exact-boundary-manifest.json` | producer-complete | `019e8c90d7fa5cd6bd19706f4b1ea7a41f00a33e9aa5d10344c1f61d495e0689` | forward/reverse counts and digests |
| `weights.json` | explicit OPEN | `57ce6d32f879cddb00169f00679966d98e44f74496fa2127dd334857abe176f1` | `unsupported_open`, 0/368 records |
| `relations.json` | producer-complete | `e374872f15e43ad1bf890a1dfeb8aeb31d61a42449e873f4b11da73aae6f7227` | 103 scalar, zero non-scalar |
| `verification.json` | producer-complete | `997c1897c026e23e6e3b0fea661a2c80ddd7a33747f284a0ad7297884f324d84` | 103 accepted, zero rejected; payoff false |
| `operation-counts.json` | explicit OPEN | `e8aace30d5b2cb2a5711824db4dc2ecbed6138db552bd645dcf0f712c6e1c7ca` | missing map costs are null, never zero |
| `resource-metrics.json` / `isolation.jsonl` | producer-complete | `662715ef593ff6fe344c746e1a97bcb8eb90bbdaa06844867601ecded5c01f83` / `1e316aad8957223ed52e91decee975af7ed3ecd5155a18d82e1ba50a9c24db16` | five uncontended phase records |
| `report.md` / `report.pdf` | generated and inspected | `e33b0e627a21f7633681883cc9b872886858b683bc5a32d8564436e55ef84bd7` / `709cfa26dd49ec60b765c81cc01e0cd8b9e85ef8f6ff3a642f651c743014153b` | source-linked run report and 4-page PDF |
| failure-attempt and launcher logs | preserved | manifest-bound | composition and float-parser defects retained |
| `dependencies.txt` | producer-complete | `7d8f2e9d33548b6232682c884507e56f19eb0cac30f0c6f2c0625e037ba056a1` | resolved dependency/features snapshot; no committed lockfile |
| `manifest.sha256` | verified | `318534b0de05bb3b09aeff3ad76fdb660fd1a41172db294be57473e1d6a6c981` | covers all 35 producer files other than itself |
| `registration.sha256` | verified | `0b34db5e210132fb07bdeb3e860ae052ae1eabac56eb5b134ca9121405bbbd7e` | additive ledger envelope; producer bytes unchanged |
| independent validation | **BOUNDARY `INCONCLUSIVE`; EMITTED-RELATION REPLAY `HOLDS`; OVERALL `BREAKS`** | recorded in validator report | counts and 103 scalar relations reproduce; canonical-set digest bytes were not frozen; maps/cost remain open |

## Independent review and ledger disposition

The review and decision were archived together at evidence commit
`b445b0a3f` (`TASK-20261006-3d7ad9`), whose diff contains exactly the two
review artifacts, two evidence records, decision, and self-neutral receipt.

| artifact | SHA-256 | disposition |
|:--|:--|:--|
| `reviews/TASK-20261006-7d34e7/validation-report.yaml` | `257af9850daaf6bde68b63a7f73caa9ab4154e2bf012e13f8022163f3da2c256` | overall `BREAKS`; A joints `HOLDS`; B boundary `INCONCLUSIVE`; B emitted-relation replay `HOLDS` |
| `reviews/TASK-20261006-7d34e7/runtime-session-receipt.json` | `bdf73bf09a1bdf7aba6ffa503d2e4698672a851bdbf7cf61d726da0559f9ae3d` | fresh independent validator session |
| `ledger/evidence/EV-SCURVE-44f051.yaml` | `75d530e2436247b45536bd1a829bfa138b1029c4b6f2f6f8a55a7592501f261f` | supports scoped `IMPLEMENTATION_WEAK` result |
| `ledger/evidence/EV-SCURVE-8e5c70.yaml` | `bcc80fa963c21c21a4d503376b3c9685264b13b656fa0477143acfa56bd8cea4` | neutral/inconclusive whole CM result |
| `ledger/decisions/DEC-20261006-581065.yaml` | `25d7ba969ac52b52a8c2ae889bd9ed961735105ecc68061f89ab4d3b3c23e09c` | `revise`; conclusions kept separate |

## Publication closeout checklist

- [x] Fill producer-observed values from hash-bound run artifacts only.
- [x] Preserve unsupported map/calibration output as an explicit gap; do not
      synthesize HOT/COLD costs.
- [x] Incorporate the independent validator verdict and receipt.
- [x] Regenerate both SVGs and `REPORT.pdf`; inspect all 10 pages.
- [x] Record final report/figure hashes without editing immutable run receipts.
- [ ] Link the crypto and evidence PRs plus the post-result ledger decision.
- [x] Keep the singular lane classified as an implementation-validation test,
      not as an isogeny or a weak elliptic curve.
