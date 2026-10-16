# Prime Degree Frontier implementation and admission controls

2026-10-15. Requested: build a challenge to find the largest prime-degree
isogenies, with all registered standard prime curves ranked separately.

| Deliverable | Acceptance |
|---|---|
| Complete registered standard prime roster | Native policy includes every pinned standard prime-field row with exact ICV1 and operation capability; unsupported rows remain visible |
| Comparable records | One leaderboard per exact registered source; degree is an exact decimal integer; no composite product or cross-curve aggregate |
| Native execution | Commands embedded in the dedicated ReleaseMe isogeny executable; no runtime Cargo or repository dependency |
| Judge | Fresh process-isolated independent replay from a copied explicit map; no stored PASS or score is trusted |
| Ladder | Descending low-order shortlist above degree 10453 by default; screening remains candidate-only |
| Failure evidence | Invalid proofs rejected; resource failures censored; fresh output directories and retained raw/resource receipts |
| Controls | Real P192/13 and P521/7 positives; composite, inflated degree, source/kernel/map mutation, compact-only, traversal, missing hash and forged-score negatives |
| Publication | Frozen native source and rules, complete native validation, reviewed report/diagram/PDF, ready PR |

Hypothesis: the new front door can admit mathematically valid public maps,
reject corrupted artifacts even after their receipt hash is recomputed, and
rank two exact source curves independently. A higher prime-degree discovery
is an optional subsequent rung, not a prerequisite for delivering the challenge.

Freeze the source and policy before measured construction/judgement examples.
Retain binary/source/registry/rules digests, commands, stdout/stderr, child exits,
raw maps, fresh certificates and process-resource receipts. Tests are controls,
not additional-map discoveries or speed comparisons. The main curve and degree
scope remain in docs/challenges/isogeny-prime-degree/README.md and the embedded
rules. No claim of an absolute mathematical maximum or current world record.

Resources: replay at most 3600 seconds/8 GiB sampled RSS per submission;
bundled construction at most 1800 seconds/8 GiB per invocation; three bounded
artifact files. No cloud allocation. External construction costs are not
independently attested by the judge. Stop and preserve output on a failed
check, process/resource failure or artifact-limit violation. A shortlist or
an old certificate cannot promote a challenge score.

The initial examples will record the complete standard-curve roster and candidate ladder,
P192/13 and P521/7 admission, negative controls and a fresh replay of the existing
P192/10453 artifact where feasible. Every submitted map is replayed again when
a leaderboard is generated. Earlier performance and 2000-additional-map work
remain separate. Canonical DLP/rho and curve-edge charts change only with new
verified evidence relevant to their own contracts; a challenge launch creates
no such claim.

After the first full-suite run observed a valid short replay censored with zero
RSS samples, require a supervisor acknowledgement before replay starts. The
installed control packet is rerun from the fixed clean source. Add twenty fresh
P192/13 judgements of the same known map as a scheduling-race stress control;
every accepted judgement must have a positive RSS sample. These are repeated
admission controls, not distinct discoveries. Also judge the historical
P192/10453 artifact with a one-second replay budget before its full-budget
leaderboard admission; retain the actual result, including censorship, without
changing its mathematical validity classification.

The final starting leaderboard also includes the committed P224/1471 artifact
alongside P192/10453 and the P521/7 control. Each is freshly replayed from its
raw map under the final rules. The historical map sources are pinned to commit
3515af45efc81beadc1ce213fa88a7eb7419487f. These re-admissions seed separate curve
rows and do not count as additional discoveries.

The final adversarial review adds a duplicate-kernel control: copying one valid
map into both positions must not substantiate COMPLETE coverage. The ff71f5772
baseline accepted this claim (degree 13 remained valid; completeness was false),
and its failing native regression is retained. The corrected judge compares
kernels modulo the source characteristic and checks distinct Frobenius roots
before accepting a COMPLETE two-line claim. Partial-map degrees remain eligible;
repeated-root totals are left unestablished. Rerun the entire native suite,
nine negative cases, twenty resource controls and three starting entries from
the next clean source freeze before publication.
