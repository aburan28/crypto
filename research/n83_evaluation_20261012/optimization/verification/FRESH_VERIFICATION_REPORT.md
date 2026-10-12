# Independent native replay of the optimization outputs

2026-10-12. Correct exact model:
`icv1-f2m83-t6151469093347-9eba6cd5`, field low terms `[0,2,4,7]`,
curve `y²+xy=x³+x²+1`, subgroup prime `8569786107849059`.
This report describes arithmetic and file-provenance verification of
completed outputs. It adds no timing measurement or algebraic discovery.

![Completed answers were replayed after all timed children finished. Repeated training and rho controls are counts of answers, not independent-target success-rate samples.](fresh_verification.pdf)

## Recorded arithmetic result

The parent ran the prepared standalone C++ checker after timings completed.
`final_receipt.json` reports **PASS**: **156 rho scalar answers, five IC
scalar answers and five IC four-point witnesses** checked independently,
with zero invalid rows.

| Output group | Distinct frozen points | Scalar-answer replays | Four-point witnesses |
|---|---:|---:|---:|
| Rho training/kernel/policy | 3 | 126/126 | Not applicable |
| Rho holdout/default/selected | 5 | 30/30 | Not applicable |
| IC holdout | The same 5 | 5/5 | 5/5 |

All 157 process dispositions are preserved as `COMPLETED_VERIFIED`. There
were 156 rho children and one IC child with five target rows. Each of the
three training points appears in 42 rho records; each holdout point appears
in six rho records plus its one IC record. The eight-point freeze and
generation law are recorded in the preserved `GENERATION_REPORT.md/PDF`.
These counts describe retained output checks; they do not certify an IID
target population, a timing ratio or a statistical algorithm winner.

The final checker reproduced the eight-fixture SplitMix/rejection stream,
subgroup membership, `[d]G=Q` fixture certificates, disjoint signed-orbit
exclusions, individual point-file hashes and five-row holdout file. The
fixture-manifest and validation-sidecar hashes match the earlier freeze.
No point, scalar sidecar, retained base or old receipt was rewritten.

## What each answer check proves

The native implementation validates the exact retained field, G, subgroup
and 99,600-point factor base with its 600 signed Frobenius orbit columns.
For every recovered rho or IC scalar it recomputes `[d]G` with standalone
polynomial-field affine arithmetic and compares the result to that row's
exact frozen public point. The scalar sidecar is compared only after the
answer's point replay. Measured-child inputs remain point-only.

For all five IC answers it checks the four actual point indices against the
retained base, compares their x codes and both intermediate sum abscissae,
and independently verifies `P1+P2+P3+P4=Q`. Orbit labels and signs were
checked against the full base before using those indices. Caps, failures or
missing answers would remain explicit records rather than certified answers;
none were present in these completed outputs.

The earlier `cdcc5432` historical model alias describes another modulus.
The exact-field checker and these outputs use `9eba6cd5` and `[0,2,4,7]`.
The parent study retains the historical identity erratum. No new identity
or library-arithmetic result is asserted here: these are applications of
the existing independent arithmetic to the recorded answers.

## File provenance and remaining limits

The post-run native SHA-256 audit also passed. It audited 157 execution
receipts and hashed 157 producer files, retaining 656 distinct file hashes.
It compared every recorded executable, existing argument input/config,
semantic output, stdout and stderr hash with actual bytes. Two executable
SHA-256 identities covered the rho and IC children. All 157 recorded
known-answer environment flags were false. The audit is preserved as
`final_file_hash_audit.json` with source `verify_file_provenance.cpp`.

This hash audit establishes consistency of the retained file provenance;
it does not independently observe historical process behavior or certify
every possible branch of library field arithmetic. The final answer replay
validates the recorded target outputs, not a universal arithmetic identity.

The full 600-row rank/factor-log transcript remains unverified. The retained
IC target and summary outputs do not contain every setup query, admitted
relation witness, echelon matrix and solved factor-log vector. Independent
reproduction of operation counts, timings, host isolation and second-host
promotion gates also remain separate. The parent evaluation owns the
statistical report and interpretation; this report makes no speedup claim.

## Evidence and reproduction

Final arithmetic receipt SHA-256:

```text
6640a9f69e8ef66f3e31a3eab70c18925bdc8140a623aac3559c01fdeb57e294
```

Native verifier sources are `check_optimization.cpp`, `native_helpers.hpp`
and `check.cpp`. Source/receipt and executable hashes are preserved beside
the report; binaries are excluded from the repository. Use the compile and
replay instructions in `README.md` and create a new output receipt path.
Run verifiers only after timed children finish. `final_receipt.json` and
all timed records were read without modification during this review.

The editable diagram is `fresh_verification.svg`, with its vector PDF
rendering. The original fixture-generation report and hash manifest are
preserved; the final manifest adds completed receipts and documentation.
