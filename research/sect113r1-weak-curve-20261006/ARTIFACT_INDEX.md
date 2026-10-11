# sect113r1 pre-admission diagnostic artifact index

Status: **PRE-ADMISSION DIAGNOSTIC; NOT AN ADMITTED SCIENTIFIC RUN; NO
INDEPENDENT IMPLEMENTATION VALIDATION**

This index binds the exact files used for the deterministic `sect113r1`
diagnostic. The producer and replay verifier are two modes of the same native
Rust implementation. Passing replay therefore establishes same-implementation
determinism, not independent reproduction.

The diagnostic ran after the source paths below were committed, but before the
branch and protocol were pushed. It intentionally records
`admitted_scientific_run=false`. It consumed no admitted experiment run ID or
randomized-search budget, and it produced no timing or end-to-end speedup
measurement.

## Revisions and environment

| item | exact value |
|---|---|
| merged source-tree revision used to build | `a3bd6e7e1072a644415d922d5da4f84a38b6eaf5` |
| certificate-hardening commit | `c0dadee6e034d96efd764ecbe9aaec5f1d32f1c9` |
| initial diagnostic implementation commit | `5a871ae5546ea650c7a8f057317f1152c11e37dd` |
| source-path dirty check at execution | clean relative to `HEAD`; the overall worktree contained the report and generated-view outputs listed below |
| Rust compiler | `rustc 1.99.0 (b940084d7 2026-09-28)`, LLVM `23.1.1`, host `x86_64-unknown-linux-gnu` |
| Cargo | `cargo 1.99.0 (5f94df478 2026-08-27)` |
| OS | `Linux 6.18.44 x86_64` |
| CPU | `AMD EPYC 9V74 80-Core Processor`; 5 logical CPUs exposed; no timing claim |

Committed native-source hashes:

| path | SHA-256 |
|---|---|
| `src/cryptanalysis/sect113r1_audit.rs` | `0e3f9fcd41473d12feb0508da5bed0badfd97253bfa6a9f9a69693aa1d929105` |
| `src/cryptanalysis/binary_velu.rs` | `1e93ff1ca36c787792240b7c017297a96a34f056a6f52d3a8a3602804946b4c6` |
| `src/binary_ecc/curve.rs` | `175656465a90db03d5dca365c4764c9b2a3db326ec7b9637cf956095e9b8ddfb` |
| `src/bin/sect113r1_audit.rs` | `38c1ffa5971f9917ca7c04719f417a764b5c2126411b1d3107cda06dd14b947c` |
| `Cargo.lock` | `4ddc5ef76119493b41352c64a551fbd3fea3b5789644949ee775975b38150bca` |

The release producer was built at
`/workspace/crypto-target-sect-v2/release/sect113r1_audit`; its 2,996,320-byte
binary SHA-256 was
`2e36f3d3ffffac7327b8419a5f902e4c3f84f8ae892d87cc1435eb700aa2eeab`.
The build path is local and is not a durable artifact; the committed source,
lockfile, toolchain record, and commands are the reproducible inputs.

## Frozen external input

The standard input is SEC 2 version 1.0, §3.2.1:

- URL: <https://www.secg.org/SEC2-Ver-1.0.pdf>
- downloaded byte count checked on 2026-10-06: `149426`
- SHA-256: `d1b16728ad83888fd656d16b99dc71bcd5541d42d848ffd0de7c62c19010d8c3`

The temporary copy downloaded for this recorded check was removed after its
digest was checked; no standards PDF is committed. The certificate freezes the
curve parameters and this source digest.

## Producer and replay commands

Build and focused tests:

```sh
CARGO_TARGET_DIR=/workspace/crypto-target-sect-v2 \
  cargo test --offline --release --lib sect113r1 -- --nocapture
CARGO_TARGET_DIR=/workspace/crypto-target-sect-v2 \
  cargo build --offline --release --bin sect113r1_audit
```

Producer and same-implementation replay:

```sh
/workspace/crypto-target-sect-v2/release/sect113r1_audit run \
  --output research/sect113r1-weak-curve-20261006/diagnostics/pre-admission-certificate.json
/workspace/crypto-target-sect-v2/release/sect113r1_audit verify \
  --input research/sect113r1-weak-curve-20261006/diagnostics/pre-admission-certificate.json
```

Replay output:

```text
verified diagnostic=PRE_ADMISSION_DIAGNOSTIC_EXACT_CERTIFICATES
sha256=1cf4a3733e98c6acdc2219c24f29c15771fc67578bdb33dd426ab643e12572d8
```

The embedded value is a semantic digest of the typed report with the digest
field cleared. The raw JSON byte hash is separate. All deserialized certificate
structs reject unknown fields, and tests reject both unknown-field additions
and a semantic mutation even when the attacker recomputes the embedded digest.

## Diagnostic artifacts

| path | bytes | SHA-256 | role |
|---|---:|---|---|
| `diagnostics/pre-admission-certificate.json` | 14,110 | `8f77c0a21d5a51df372ce0a3a69784c58c3e404e5f1f52ad953c5e5714eeed24` | authoritative typed certificate; embedded semantic digest `1cf4a3733e98c6acdc2219c24f29c15771fc67578bdb33dd426ab643e12572d8` |
| `degree5_curve_records.json` | see Git object | `621846313f27b2d7b30e9c7430674937b93ce8cc232c85adadd37a0668c73d26` | derived endpoint identity view |
| `isogeny-routes.json` | see Git object | `fa83ed9cf56ae3c83febe7d8749e72aac989c73eb6d132a7efa2cdd6a50c6ee3` | derived ordered-route view |
| `README.md` | see Git object | `8c5e2f2aaa8ed8cce513cb63c20f945311b2f2e72bdb3d644e3dfef64418a826` | prospective protocol and diagnostic summary |
| `REPORT.md` | see Git object | `586b7ca6707f8f0688a283d780506e105798c2441a5eead2fe7b60d863c6938e` | source report |
| `REPORT.pdf` | 120,350 | `9b2a33d0e43162bf17d82c35c3f02bd2e3192bb22a86d0f9e8d1024527aadd3a` | rendered nine-page report |
| `evidence-flow.dot` | see Git object | `e15ddbe4e407619f65de0f9f6e3458b147e0c3529322298313deca2eb37f3f0e` | editable diagram source |
| `evidence-flow.svg` | see Git object | `e4dad52aaac50eca418343fb2cf3ca14bae6ec916811c179b360ebaa2d687be9` | rendered diagram |
| `evidence-flow.pdf` | 29,572 | `fc1abdf2ecf0dbc14d7526af924b8e429678ee6be37abbbbacf97c17de323465` | rendered one-page diagram |
| `pdf-images.lua` | see Git object | `4b7e01eff969511457213e213094dbc912ecb378c090bce4c3b82747e35d2f2b` | local SVG-to-PDF render filter |

This index does not hash itself; the final Git commit and pull-request head bind
its bytes. Any edit to the source Markdown or DOT requires regenerating the
rendered PDF/SVG and updating this table.

## Registry-derived views

The two valid degree-5 endpoint identities were added to the authoritative
registry, then every dependent view was rebuilt after merging the current
`origin/main`:

```sh
python3 scripts/build_curve_registry.py
python3 scripts/build_curve_registry.py --check
cargo run --offline --release --bin curve_cover_check --
cargo run --offline --release --bin curve_cover_check -- --check
python3 scripts/build_ic_leaderboard.py
python3 scripts/build_ic_leaderboard.py --check
python3 scripts/build_lab_browser.py
python3 scripts/build_lab_browser.py --check
python3 scripts/test_cover_catalog.py
python3 scripts/site/test_build.py
```

| path | SHA-256 |
|---|---|
| `docs/curves/registry.json` | `2342a3d22dd3a947dfc8ac2a791e2f42061dad4e7b723db6bd0abe7b851b64a6` |
| `src/cryptanalysis/curve_aliases.json` | `1de174e3da19076f02b1a0c57369f007b6a517a81051a3f2561c53288237941e` |
| `docs/curves/covers.json` | `4300d27564ed72b21623ea4e68bebd36fc76d33714b6f144a4b6ef8b6d4e4865` |
| `docs/ic/leaderboard.json` | `616d04d06339b69c8adc5203594f583e4c331a344e18f882429528c3ee78a097` |
| `docs/ic/LEADERBOARD.md` | `ad6313f1ab13a645bf342d7714fc1fea4bb6d27b99b128a19e51e0083ae6e8d4` |
| `docs/ic-leaderboard.html` | `c6f3d37dfc045e11dacce087f1730d824b8206cb77cff3208c3e9c4fc95a1075` |
| `docs/browser/data.json` | `91379fe1038d22919cb5130862e97fe34cb1d83c79d0dbb87991a2af03f8964d` |

Regeneration reported 124 curves; 124 verified cover records and zero invalid or
unsupported records; 46 methods; 29 factor bases; 258 candidates; 30 rounds; 27
sessions; and 1,956 yields. Each new endpoint occurs once in every authoritative
roster. The generated same-field degree-3 cover records have
`subgroup_transfer=not_tested` and `dlp_advantage=null`; they are not weakness
findings.

## Validation receipts

- Final focused native tests after the upstream merge: **11 passed, 0 failed,
  3,696 filtered out**.
- Direct `rustfmt` of `sect113r1_audit.rs` and repository `git diff --check`:
  passed.
- Registry builder and `--check`: passed.
- Native cover generator and `--check`: 124 verified, 0 unsupported, 0 invalid.
- IC leaderboard builder and `--check`: passed.
- Lab-browser builder and `--check`: passed.
- Certificate-to-record/route cross-view check: both digest links, both endpoint
  identity tuples, and both kernel/dual/verdict tuples matched (2/2).
- Cover-catalog tests: 2/2 passed.
- Site/HTML tests: 47/47 passed.
- Curve-name/registry policy check: 0 problems.
- Optional Playwright browser test: not run because its bundled Chromium shell
  is absent; this does not validate the arithmetic or registry joins.
- `REPORT.pdf`: nine US-letter pages; inspected at rendered resolution.
- `evidence-flow.pdf`: one portrait page; inspected at rendered resolution.

An earlier full library run on the initial diagnostic branch reported 3,580
passed, 16 failed, and 109 ignored. The failures were pre-existing or
environmental: missing sparse fixtures, sandbox-denied collaboration networking,
missing CryptoPro chain fixtures, and one existing BIKE mismatch. The focused
final module and generated-view checks above are the validation attributable to
this change; a new all-library clean pass is not claimed.

## Claim and reconstruction boundaries

- `implemented_diagnostic_gates_passed=true` means only that the explicitly
  implemented bounded gates replayed. It is not admission, complete attack
  coverage, production exploitability, or independent validation.
- The certificate recomputes source and endpoint ICV1/EC1/UID identities, the
  twist order, forward and dual kernel certificates, subgroup images, dual
  composition on `G` and the dependent planted target `[d]G`, and exact target
  pullback. The endpoint/route JSON files are derived views, not separate native
  certificates.
- The diagnostic does not serialize dual-kernel coefficients. An outside checker
  must rederive them from the codomain; no independent degree-5 oracle was run.
  Global morphism validation over the broader protocol set of kernel points,
  infinity, and sign partners remains an admitted-run obligation.
- The diagnostic proves only `ord_113(2)=28` for the GHS lane. Magic number,
  type, genus, and attack applicability are deliberately uncomputed, so no GHS
  competitiveness conclusion is asserted.
- The singular companion is a nodal cubic, not an elliptic curve and not an
  isogenous representative. Its bounded recovery is conditional on an
  attacker-selected `BinaryPoint` reaching unchecked `scalar_mul` with a
  distinguishable full-point output. No production wrapper exposure was shown.
- The full-width singular fixture recovered a residue only. Its remaining
  59-bit prime component was not executed; `29.710003` bits is a modeled rho
  cost, not an observed solve.
- No codomain discrete logarithm, randomized rho timing, representative-specific
  speedup, prevalence estimate, or exhaustive isogeny-class search was executed.

## Additive 2026-10-07 GHS structural screen

The frozen certificate, parent report, and parent evidence-flow diagram above
remain unchanged. A separate additive pre-admission diagnostic screens the
registered source model and the two registered degree-5 codomains named by
their ICV1 slugs. Its native same-codebase result is magic number 113, type I,
genus `2^112-1`, and zero candidates within the frozen genus bound 64 at each
of those three nodes.

- [`addenda/ghs-three-node-20261007/REPORT.md`](addenda/ghs-three-node-20261007/REPORT.md)
  gives the bounded interpretation and open transfer obligations.
- [`addenda/ghs-three-node-20261007/manifest.json`](addenda/ghs-three-node-20261007/manifest.json)
  records the exact producer revision, toolchain, binary and source hashes,
  commands, parameters, raw output hashes, byte-identical replays, and focused
  verifier receipt.
- [`addenda/ghs-three-node-20261007/SHA256SUMS`](addenda/ghs-three-node-20261007/SHA256SUMS)
  binds the complete additive bundle.
- [`tests/sect113r1_ghs_evidence.rs`](../../tests/sect113r1_ghs_evidence.rs)
  is a same-codebase replay and custody verifier, not independent arithmetic
  confirmation.

This is a three-node negative structural result only. It does not revise the
parent certificate in place, exhaust the isogeny class, construct a descended
Jacobian, transfer the subgroup to one, solve a DLP, establish a speedup, or
support a prevalence estimate.
