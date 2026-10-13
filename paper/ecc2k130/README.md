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
make eprint    # eprint-submission.zip; the form, fields and rules are in EPRINT.md
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
- §7 The other GPU kernels: an evidence-class inventory, 256-bit
  secp256k1 arithmetic and its instruction ledger, rho with order-six
  folding and baby-step giant-step, kangaroo, P-256/P-384, the ECC2K-95
  kernel and pair table, the polynomial-system and linear-algebra
  solvers, and the F4 elimination kernels that did run on a device.
- §8 The FPGA engine: pre-board estimate, design, measured F2 images.
- §9 The campaign: configuration, the 2^28.41 interval, fleet
  economics, storage and certification, state.
- §10 Structural alternatives: index calculus on the Koblitz curve,
  point decomposition against the budget.
- §11 Related work. §12 Conclusion.
- Appendix A: reproduction commands.

Every number is footnoted to a repository artefact (`\art{path}`);
paths prefixed `cryptanalysis:` are in the sibling library.

## Revision history

- Draft 1: first complete text.
- Draft 2: every section fact-checked line by line against its cited
  files (55 corrections, among them the two-chain B200 verdict, the
  batch-18 control, the cycle-rule v2/v3 validation attribution, the
  inversion chain the client actually uses, and the FPGA host's
  `--verify N` sampling); bibliography checked against Crossref; the
  reporting rule restated in full so the paper does not depend on the
  companion draft; an availability paragraph naming the three public
  repositories; the campaign state refreshed from the published feed
  (`data/status_20261010T154157Z.json`); related work extended to
  ECC2K-108 and the 2016 binary-curve record.

- Draft 3: a new section on the repositories' other GPU arithmetic
  (prime-field secp256k1 / NIST and binary-field ECC2K-95 / solver
  kernels), written from their READMEs, optimisation notes and
  experiment records, with an evidence-class table up front. Apart
  from the F4 elimination kernels and one Macaulay scheduling run, none
    of that code has run on a device; the section quotes static
  instruction counts and emulated walk constants and says so.

- Draft 4: the NIST subsection rewritten after the P-256/P-384 rho
  kernel's first device run (2026-10-12, RTX PRO 6000 via Modal):
  the two first-run bugs, the distinguished-point correctness check,
  the field and walk table (1.340 / 0.438 G steps/s, 1.25x / 2.09x
  over the first working cells, classed as engineering), and the
  fifteen-curve CPU solver's audit and measured interval, Frobenius
  folding and Smart-attack results. Curve25519 added to the
  secp256k1 row. Abstract, introduction, conclusion and the
  evidence-class table updated to say which kernels have now run.

## Known caveats carried into the text

- The RTX PRO 6000 σ-walk headline has moved past the audited 14.64 B/s
  (15.44 fused, 15.88 with nibble tables); the campaign preset still
  ships 14.64 and the paper leads with that number, naming the others.
- The 20.08 B/s table-walk rate is a rate, not a cost per solve; the
  cycle-rule v2 cost projection (0.81–0.86× σ) is unmeasured, and the
  v3 rule's throughput and walk constant are not yet measured at all.
- The RTX PRO 4500 bridge map's DP34 collection rate, the L40S
  collection rate and the G7e family are unmeasured.
- The packed-walk witness cost (34.7%) is stated in `CAIRN-WITNESS.md`
  but only the job script is frozen under `benchmarks/witness-cost/`.
- The IC feasibility note's fitted constant `C` has three readings in
  the source; the paper carries the values recomputed from the exact
  subgroup orders (0.33, 1.51, 1.27) and says the extrapolation uses
  `C = 1`.
- Whether the storage bucket's copy of `campaign.json` had its iteration
  guard raised to 2^32 before the fleet stopped is not recorded.
- The other-GPU section's numbers are static PTX counts (clang,
  `ptxas` 12.9) and host-emulation operation counts, not device rates,
  except the NIST table (`gpu/nist/RESULTS.md`), a first device run at
  one block size that has not yet been through the pre-registered
  three-repeat protocol;
  the `gpu/ecc` occupancy and stack tables predate the order-six fold
  and the cycle ring and are quoted as such.
- The author block is restored (Adam Buran, Kevin Nguyen, Isogeny Labs) for the ePrint release; `make deanonymize` with an `AUTHORS`
  file restores it.
