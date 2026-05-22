# LaTeX paper draft

`structural_completeness.tex` is a paper-quality LaTeX version of
[`PAPER_STRUCTURAL_COMPLETENESS.md`](../PAPER_STRUCTURAL_COMPLETENESS.md).

**Status**: anonymized for double-blind review.  Restore real
authors via `make deanonymize` after creating an `AUTHORS` file
(see `AUTHORS.example`).

## Building

Requires a TeX distribution (TeX Live, MacTeX, MiKTeX) with
`pdflatex` and standard packages (`amsmath`, `amssymb`, `hyperref`,
`booktabs`).

Use the provided `Makefile`:

```
make           # build the PDF (two passes for cross-refs)
make arxiv     # produce arxiv-submission.tar.gz
make eprint    # produce eprint-submission.zip
make check     # sanity check unresolved references + page count
make clean     # remove build artefacts
```

Manually (without Makefile):

```
cd paper
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
- §8: the structural-completeness theorem with proof.
- §9: open questions and out-of-scope directions.
- §10: reproducibility — full list of companion PARI scripts.
- Bibliography.

## Companion artefacts

Everything cited in the paper is in this repository:

- PARI scripts: `../secp256k1_cm_audit/*.gp`
- Rust integration: `../tests/curve_audit.rs`
- Research notes: `../RESEARCH_*.md`
- Markdown paper (kept in sync with this LaTeX): `../PAPER_STRUCTURAL_COMPLETENESS.md`

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
