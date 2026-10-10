# ECC2K-130 on Commodity Hardware

`ecc2k130.tex` is a paper-quality LaTeX account of the ECC2K-130
Pollard-rho campaign documented in [`../../ecc2k130/`](../../ecc2k130/),
[`../../hdl/ecc2k130/`](../../hdl/ecc2k130/) and the research notes this
repository holds, with the structural alternatives priced in
[`../../RESEARCH_ECC2K130_IC_FEASIBILITY.md`](../../RESEARCH_ECC2K130_IC_FEASIBILITY.md)
and the sibling library [aburan28/cryptanalysis](https://github.com/aburan28/cryptanalysis).

**Status**: research draft, anonymized for double-blind review. Restore
real authors via `make deanonymize` after creating an `AUTHORS` file
(see `AUTHORS.example`).

## Building

Requires a TeX distribution with `pdflatex`, `bibtex`, and the packages
`amsmath`, `amssymb`, `amsthm`, `mathtools`, `hyperref`, `booktabs`,
`tabularx`, `graphicx`, `xcolor`, `microtype`, `enumitem`, `caption`,
`tikz`, `pgfplots`, `natbib`. On Debian/Ubuntu:

```
sudo apt-get install texlive-latex-extra texlive-science texlive-fonts-recommended \
                     texlive-bibtex-extra poppler-utils
```

```
make           # PDF (pdflatex + bibtex + two more passes)
make check     # unresolved refs, overfull boxes, page count
make arxiv     # arxiv-submission.tar.gz
make eprint    # eprint-submission.zip
make deanonymize
make clean
```

Output: `ecc2k130.pdf`. A built copy is committed so a reader without
TeX can open the paper.

## Contents

- Abstract.
- §1 Introduction: the challenge, contributions, method, the two
  repositories, what is not claimed.
- §2 The curve, the type-II optimal normal basis, the golden model,
  ECC2K-95 as the end-to-end control.
- §3 The walk: the σ-walk, the r-adding table walk, the x-only bridge
  map, the measured walk constant, cost per solve.
- §4 Field arithmetic: bitsliced and packed backends, LOP3 fusion,
  bit-for-bit compatibility.
- §5 The GPU client: the carry-less floor and roofline, acceptance
  rules, the 0.85→14.64 B/s ledger, the table walk at 20.08, the RTX
  PRO 4500, the twelve-SKU survey, every declined candidate.
- §6 CPU clients.
- §7 The FPGA engine: pre-board estimate, design, measured F2 images.
- §8 The campaign: configuration, the 2^28.41 interval, fleet
  economics, storage and certification, state.
- §9 Structural alternatives: index calculus on the Koblitz curve,
  point decomposition against the budget.
- §10 Related work. §11 Conclusion.
- Appendix A: reproduction commands.

Every number is footnoted to a repository artefact (`\art{path}`);
paths prefixed `cryptanalysis:` are in the sibling library.

## Known caveats carried into the text

- The RTX PRO 6000 σ-walk headline has moved past the audited 14.64 B/s
  (15.44 fused, 15.88 with nibble tables); the campaign preset still
  ships 14.64 and the paper says which is which.
- The 20.08 B/s table-walk rate is a rate, not a cost per solve; the
  cycle-rule v2 cost projection (0.81–0.86× σ) is unmeasured.
- The RTX PRO 4500 bridge map's DP34 collection rate, the L40S
  collection rate and the G7e family are unmeasured.
- The IC feasibility note's fitted constant `C` has three readings in
  the source; the paper reports the range and uses the pessimistic fit.
