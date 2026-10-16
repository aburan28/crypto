# Index Calculus for the Elliptic-Curve Discrete Logarithm, Measured

`ic_program.tex` (built: [`ic_program.pdf`](ic_program.pdf), 21 pages) is the
programme record of the repository's ECDLP index-calculus work from
2026-07-01 to 2026-10-14: the routes pursued, the measurement infrastructure
(ecbench, the boundary ledger, claim hygiene), what improved and on which
boundary, the negative results, and where the frontier stands.

**Status**: research draft, anonymized, not peer-reviewed, not submitted,
not cleared. A synthesis of committed evidence: **no new measurement** was
run for it. Every number cites a committed artefact (Appendix B of the PDF;
plotted values with their sources in [`figures/data.json`](figures/data.json)).
Restore authors with `make deanonymize` after creating an `AUTHORS` file
(see `../AUTHORS.example`).

**Change class** (per `AGENTS.md`): none of its own; it reports results
already classed by their studies. With one modelled exception (the `k = 4`
crossover of §8.4, an extrapolation fed by a measured constant) every
improvement it reports is *engineering* or *accounting*.

## What it says

- **Measurement** (§3): the unit `S = GAE / √r`, cold one-target cost
  beside the online window, sealed records and exact replays, isolation
  levels L0–L3, paired cluster-bootstrap ratios of totals, the frontier
  comparability key, and seven corrections the accounting forced on earlier
  claims (exploratory online ratios, a Buchberger baseline labelled F4,
  unpriced F4 word XORs, cross-host counter drift, an `f4_fp` soundness
  bug, F6-IC's charged cost, the ECC2K-130 constant).
- **Chronology** (§4): commit activity by ISO week and theme, four phases
  from 6 September to 14 October, with the July–August prime-field prelude.
- **Root-finding decomposition** (§5): how one suggestion in an audited
  manuscript (the roots of `S₃(x_R, x_F, X)` are `x(R ± F)`) became a
  proved identity (Lean), prime and binary oracles, `m ≤ 4` pair tables,
  the PDP family, a CUDA kernel, Kummer-tower and dense-index variants, and
  four ecbench improvement loops; stage factors 383–1,535× (binary vs
  matrix F4), 3,021–5,802× (NIST vs Buchberger), end to end 50–105× (vs
  msolve); best complete Koblitz `n = 13` arm `S` 38.82 → 22.97, still
  1.55× strong rho.
- **Gröbner engines** (§6), **Koblitz curves** (§7), **prime and extension
  fields** (§8: PKM tower 29–141× rho; FGHR `D₄` cuts `C₄` 47.5–48.4×,
  modelled crossover `2^150 → 2^135`), **transfer routes** (§9: the
  transport-conservation lemma and the negative results).
- **Frontier** (§10): a table of every measured improvement by boundary and
  class; cold one-target frontier still rho; one admitted online point
  (compact-orbit IC at Koblitz `n = 61`, 13.9× [9.8, 19.1] over strong rho,
  64 targets, set-up excluded; 83× rho with set-up charged).
- **Appendix A**: how to organise the material for publication (below).

## What it does not say

- No ECDLP speed-up on any deployed curve, no key recovery, no asymptotic
  sub-rho claim, no change to the Gaudry–Diem picture for prime fields, no
  ECC2K-130 or P-256 security statement.
- Stage ratios are not end-to-end speed-ups; L0 wall ratios are
  descriptive; lower-bound costs are labelled; exploratory online ratios
  (n = 61/71/73/83) are not admitted results.
- The commit-activity figure counts keyword matches on commit subjects; it
  indicates where effort went, not how much of it was index calculus.

## How to organise the publication

The PDF is a programme record and deliberately long. Appendix A proposes
splitting it into five self-contained papers, each citing this record for
provenance:

| paper | content | owed before submission |
|:--|:--|:--|
| P1. Root-finding decomposition | §5; extends `../rootfinding-decomposition/` with the PDP family, GPU, Kummer tower, ecbench loops and the audit corrections | novelty audit; Magma/msolve cross-check on every row; a GPU run on hardware; `m ≥ 3` prime measurements |
| P2. ecbench: one unit for every ECDLP method | §3, Definition 1, corrections, isolation rubric, frontier synthesis | an L3 cross-host replay certificate; inherited-F4 counter drift resolved; a four-size bound per incumbent |
| P3. Index calculus on Koblitz curves, measured | §7, compact orbit, online ladder, admitted `n = 61` panel, ECC2K-130 extrapolation | isolated panels at `n = 71, 73, 83`; the constant `C` reconciled; a cold `n = 83` session |
| P4. Symmetry at `k = 4` | §8.4 | a second prime and more seeds; the degenerate target explained; a measured end-to-end point |
| P5. Negative results | §9, lit2026 negatives, SAT scaling, symmetrised bases, fibre-aware IC, Lemma 1 | literature check of Lemma 1; the conductor-gap yield experiment |

Each paper should follow: the question stated as a boundary (Definition 1)
→ the mathematics with proofs or certificates → implementation and
fail-closed verification → results on the stated boundary first, stage
results second → what is not established → reproduction (commit, binary
hash, spec, session id, replay certificate). P1, P3 and P4 suit a
cryptanalysis venue via IACR ePrint (`iacr-eprint` skill); P2 an artefact or
methodology track.

## Build

```sh
make          # pdflatex + bibtex
make check    # rebuild; fail on LaTeX errors, unresolved refs/cites, overfull boxes
```

Figures are rendered from `figures/data.json` (matplotlib, recorded values
only; PDF for LaTeX, SVG and 150 dpi PNG for browsing):

```sh
python3 figures/make_figures.py
```

## Reproduction of the commit-activity tally

From a full (unshallow) clone at the paper's base commit:

```sh
git log --since=2026-07-01 --format='%ad|%s' --date=format:'%G-W%V'
```

Themes are case-insensitive regular expressions on the subject:

| theme | regex |
|:--|:--|
| Koblitz / binary IC | `koblitz\|n83\|n53\|n61\|n71\|n73\|orbit\|frobenius\|ggmp\|binary ic\|boundary` |
| Gröbner F4/F5/F6 | `\bf4\b\|\bf5\b\|\bf6\|groebner\|gröbner\|macaulay\|m4ri\|inherited` |
| root finding / PDP | `root.?find\|rootfind\|semaev\|pdp\|pair.?table\|subtract\|summation` |
| ecbench / measurement | `ecbench\|isolab\|replay\|audit\|icms\|bound` |
| isogeny / transfer | `isogen\|conductor\|fiber\|ghs\|weil\|transfer\|glv\|gls\|volcano\|jacobian\|genus` |
| prime-field IC | `prime\|gaudry\|pkm\|kummer\|p-256\|p256\|nist` |

## Before submission

- [ ] Split per the table above; each split paper gets its own directory here.
- [ ] Verify every bibliography entry against the publisher's record.
- [ ] Re-read every cited artefact at the commit of submission (numbers can
      be superseded by later studies; frozen evidence is not rewritten).
- [ ] `make deanonymize`, `make check`, then the `iacr-eprint` workflow.
