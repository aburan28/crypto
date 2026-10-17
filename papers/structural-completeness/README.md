# LaTeX paper drafts

This directory holds the structural-completeness paper (below) and, in
[`ecc2k130/`](ecc2k130/), the ECC2K-130 campaign paper with its own
Makefile and README.

`structural_completeness.tex` is a paper-quality LaTeX version of
[`research/notes/cm-isogeny/PAPER_STRUCTURAL_COMPLETENESS.md`](../../research/notes/cm-isogeny/PAPER_STRUCTURAL_COMPLETENESS.md).

**Status**: anonymized for double-blind review.  Restore real
authors via `make deanonymize` after creating an `AUTHORS` file
(see `../AUTHORS.example`).

## Building

Requires a TeX distribution (TeX Live, MacTeX, MiKTeX) with
`pdflatex` and standard packages (`amsmath`, `amssymb`, `hyperref`,
`booktabs`).

Use the provided `Makefile` (the rules are the shared `../paper.mk`):

```
make           # build the PDF (two passes for cross-refs)
make arxiv     # produce arxiv-submission.tar.gz
make eprint    # produce eprint-submission.zip
make check     # sanity check unresolved references + page count
make clean     # remove build artefacts
```

Manually (without Makefile):

```
cd papers/structural-completeness
pdflatex structural_completeness.tex
pdflatex structural_completeness.tex    # second pass for cross-refs
```

Output: `structural_completeness.pdf` (~15 pages).

## Contents

- Abstract: structural-completeness theorem stated.
- §1: introduction and motivation.
- §2: setup — isogeny graphs and CM.
- §3: attack surface taxonomy (the seven families).
- §4: the seven structural blocks (B1–B7).
- §5: novel proposals — CGA-HNC and VFCG-ρ, both null.
- §6: the Howe-gluing structural fact (existence, not explicit
  construction; errata regarding the naive construction).
- §7: empirical verification across secp256k1, P-256,
  brainpoolP256r1, Curve25519.
- §8: the structural-completeness theorem with proof sketches.
- §9: open questions and out-of-scope directions.
- §10: reproducibility — full list of companion PARI scripts.
- Bibliography.

## Corrections

- **2026-10-17** (sheriff case SHF-20261017-e4d7f6). P-N6: the abstract
  (`.tex`, `abstract.txt`, Markdown twin) said "The theorem is proven by case
  analysis"; §8 holds proof sketches, so it now says "argued by case analysis
  (proof sketches in Section 8, not yet independently reviewed)". The Markdown
  twin's stale "proof is *unconditional*" paragraph is deleted. P-N4: the
  twist ℓ-rank remark now says both scripts skip `r = 0` (dropping
  `(7,7,5)`), says "hypothesis classes", and gives the distinct `(ℓ,p,#E)`
  check counts (26 278 instances = 336 distinct checks; 9 864 = 249),
  recounted independently. P-N5: `recount/README.md` provenance corrected and
  `recount/exp_zdoubleprime_lemma_extended.txt` regenerated.

## Companion artefacts

Everything cited in the paper is in this repository:

- PARI scripts: `../../secp256k1_cm_audit/*.gp`
- Rust integration: `../../tests/curve_audit.rs`
- Research notes: `../../RESEARCH_*.md`
- Markdown paper (kept in sync with this LaTeX): `../../research/notes/cm-isogeny/PAPER_STRUCTURAL_COMPLETENESS.md`

## Status

Research draft, *not* peer-reviewed.  Intended as the starting
point for a submission to e.g. IACR ePrint, AsiaCrypt, ANTS, or
the journal *Designs, Codes and Cryptography*.

Suggested next steps before submission:
1. Independent verification of the seven structural blocks via
   another reviewer / re-implementation.
2. Tighten the proof sketches in §8 with explicit citations
   and inline lemmas.
3. Add a section on the GLV-HNP research direction (currently
   in `../RESEARCH_GLV_HNP.md`) as motivated future work.
4. Consider adding an Appendix with the explicit Mestre-algorithm
   description (currently in `../RESEARCH_MESTRE_HOWE.md`).
