# Weak-curve conjectures and discovery methodology

This addition provides **seven open hypotheses** about class-support laws,
additional covering support, primitive intersections of field towers,
higher-degree quotient geometry, torus-lift degree rigidity and polarization
profiles. Each has a domain, prior-work boundary, first experiment, logical
falsifier and Lean obligation. One exploratory boundary hypothesis is
already refuted by a fixed native review of **87,486 complete census rows**.

- [Conjectures and distinctive properties](CONJECTURES.md)
- [Machine-readable hypothesis ledger](conjectures.json)
- [Discovery, verification and scaling methodology](METHODOLOGY.md)
- [Primary-source audit and published-class backlog](LITERATURE.md)
- [Lean proof map](PROOF_MAP.md)
- [Requested scope and preregistered review](PROTOCOL.md)

The candidate geometries and open conjectures have not been promoted to the
verified family catalogue. Published Type I/II coverings, Scholten forms,
genus-3 isogenies and odd-prime cover mechanisms keep their attribution.
The 2013 Iijima–Momose–Chao classification identifies additional published
cases to normalize and implement. Exact novelty of the proposed statements
remains unresolved against the bounded reviewed corpus.

## Exploration completed

The native review reads only the seven already certified class CSVs. It
checks the union labels, duplicates, complete signed ordinary trace domain,
the class count `p^2*(p-1)` and the cubic depth-at-least-two condition.
These results preserve the original family-specific scope.

| p | Classes | Cubic only | Quadratic only | Both | Union zeros | Zeros outside gap <=64p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 7 | 294 | 38 | 122 | 88 | 46 | 10 |
| 11 | 1,210 | 100 | 530 | 442 | 138 | 58 |
| 13 | 2,028 | 128 | 918 | 800 | 182 | 106 |
| 17 | 4,624 | 220 | 2,156 | 1,978 | 270 | 168 |
| 19 | 6,498 | 238 | 3,118 | 2,852 | 290 | 150 |
| 29 | 23,548 | 360 | 11,544 | 11,150 | 494 | 302 |
| 37 | 49,284 | 466 | 24,278 | 23,886 | 654 | 390 |

There are 2,074 union-zero rows, of which 1,184 violate the proposed strip.
At p7, t=±22 is a verified counterexample: both branch labels are zero and
`2*7^3-22=664>64*7=448`. Thus a global explanation must accommodate interior
gaps. The seven zero counts motivate the separate growth conjecture, while
the refuted geometric screen stays in the ledger. [REVIEW.tsv](REVIEW.tsv)
retains the full descriptive counts and [FALSIFIERS.tsv](FALSIFIERS.tsv)
retains every violation.

Three Lean theorems check the universal inclusion-exclusion and partition
identities used in this table. They do not establish the imported curve
labels or prove the seven open hypotheses. The earlier 814-model inventory,
property records, class oracles and proof additions remain unchanged.

## Reproduction

Run in this directory using the pinned Lean toolchain:

```sh
mkdir -p .lake/build/bin
rustc --edition=2021 review.rs -o .lake/build/bin/review
.lake/build/bin/review
rustc --edition=2021 --test review.rs -o .lake/build/bin/review-tests
.lake/build/bin/review-tests
jq -e -f validate.jq conjectures.json
lake build
lake env lean AxiomAudit.lean
lake env leanchecker --fresh --verbose CensusCounting
```

The review binds all fixed input CSV hashes in `INPUTS_SHA256.tsv`. Saved
output comparisons, exact verification results and retained failed attempts
are documented in `VALIDATION.md`. New geometric pilots and 192–252-bit
property runs are planned, not executed by this initial review.
