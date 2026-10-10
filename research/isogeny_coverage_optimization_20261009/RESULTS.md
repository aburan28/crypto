# Isogeny native optimization and curve coverage

Dated 2026-10-09. The implementation is Rust plus assembly. This round extends the dedicated binary beyond P192/P224 and measures native ARM Linux instruction counts under the local Docker VM.

| Command | Baseline median instructions | Final candidate median instructions | Baseline / candidate | Reduction |
|---|---:|---:|---:|---:|
| screen-p224 | 2471560507 | 78735743 | 31.3906 | 96.814% |
| verify-p224 | 50199483972 | 47033190946 | 1.0673 | 6.307% |
| search-p192 | 12331509351 | 12477584697 | 0.9883 | -1.185% |

Five rounds per panel. Candidate lists, replay records and map bytes match the frozen baseline. The screening improvement includes setup/output; verification excludes construction; the full-search panel includes the constructor child and replay. The small full-search case regresses slightly. No native Mac wall-time claim is made. Candidate 1's small verifier improvement and Candidate 2's larger instruction reduction remain separately recorded in PROFILE_SUMMARY.json.

## Coverage

The catalogue originally held 345 models: 200 prime, 72 binary, 60 Koblitz, 4 subfield and 9 extension models. It includes 170 short-Weierstrass prime models with construction support through 640 bits, of which 142 have registered complete public subgroup generators. The four newly verified P384/P521 codomains bring the canonical model count to 349.

Montgomery/Edwards entries and binary/extension models are visible and can be structurally screened from their known orders, but their construction/independent-map adapters remain unresolved. Catalogue completeness means the pinned repository inputs; the importer retains seven unresolved standards records and disclaims a worldwide census.

## Implemented and verified

Screening now resolves the curve and trace once outside the prime loop. Independent replay uses exact-size fields, binary extended-GCD inversion, a fast path for one, and allocation-free small-integer conversion. A recorded profile identified field multiplication as 78.65% of verifier instructions; forced inlining was tested as a second candidate. Producer and verifier arithmetic remain independent. No exact kernel, rational-map or subgroup check was removed.

The CLI suite passed 58 tests, field reference checks passed 3 tests across widths through P521, and the existing standalone suite passed. P384/19 and P521/7 each produced two independently certified maps, including 20 scalar transport checks per map. P224/1471 replay passed the optimized implementation.

## Scope and remaining gates

The baseline is c87d96eda, Candidate 1 is 6e3147ad9 and Candidate 2 is 3413c5329. Rust 1.93.1, identical release flags and the same 345-model source snapshot were used for paired measurements. Callgrind 3.19.0 records guest ARM Linux user-space instructions. This is not native Mac throughput and no end-to-end solver/security claim is made.

Root library compilation still has 643 pre-existing errors. Hosted CI and release publication remain gated. The old study sources/evidence remain frozen at their recorded revisions. No IC performance ledger row is promoted.

[Protocol](PROTOCOL.md), [profile spec](PROFILE_SPEC.json), [summary](PROFILE_SUMMARY.json), [models](VERIFIED_MODELS.json), [diagram](COVERAGE.svg), [PDF](REPORT.pdf), and the raw profile archives in evidence/ carry the receipts. The original candidate's regression and all failed test panels are retained.
