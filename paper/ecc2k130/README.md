# Pollard Rho on NVIDIA Blackwell: ECC2K-130, Carry-less Multiplication and the NIST Curves

`ecc2k130.tex` is a paper-quality LaTeX account of Pollard rho on
current NVIDIA hardware, with ECC2K-130 as the yardstick: the
optimisation ledger, the PTX `clmad` instruction across dies, twelve
GPUs from Turing to Blackwell, the same optimisations on NIST P-256 and
P-384, and the 2009–2012 rates for comparison. The campaign it reports
is documented in [`../../ecc2k130/`](../../ecc2k130/),
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

All rates are in billions of iterations per second on one device.

- Abstract.
- §1 Introduction: rho on Blackwell, contributions, the reporting rule,
  availability, what is not claimed.
- §2 Where 2009 left off: the 2009–2012 rates per platform (Core 2,
  Cell, GTX 295, Spartan-3, ASIC estimate) in one unit, and the same
  yardstick today (about 230× per GPU card, 74× per FPGA).
- §3 The curve, the type-II optimal normal basis, the golden model,
  ECC2K-95 as the end-to-end control.
- §4 The walk: the σ-walk, the r-adding table walk, the x-only bridge
  map, the measured walk constant, cost per solve.
- §5 Field arithmetic: bitsliced and packed backends, LOP3 fusion,
  bit-for-bit compatibility.
- §6 What current NVIDIA hardware changes: PTX 9.3 `clmad` (sm_80+,
  CUDA 13.3), its issue rate on the RTX PRO 6000, B200 and H100, why the
  same change wins on one die and loses on another, the integer pipe.
- §7 The GPU client: floor and roofline, acceptance rules, the
  0.85 → 14.64 ledger, the table walk's arithmetic, σ-walk against table
  walk (cycle rules v1/v2/v3), the RTX PRO 4500, twelve GPUs per card,
  per SM and per dollar, every declined candidate.
- §8 CPU clients.
- §9 The same ladder on the NIST curves: P-256 1.075 → 1.632 (1.52×),
  P-384 0.210 → 0.548 (2.61×), the fifteen-curve CPU solver.
- §10 The other GPU kernels: evidence-class inventory.
- §11 The FPGA engine. §12 The campaign. §13 Structural alternatives.
- §14 Open questions and next optimisations: exact table-walk speed,
  σ-walk integer slots, pipe and `clmad`-rate probes, other Blackwell
  parts, the NIST walk, method.
- §15 Related work. §16 Conclusion.
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
- Draft 6: refocused on rho on Blackwell and retitled. New §2 (2009–2012
  rates from Bailey et al., the INDOCRYPT 2010 GPU paper, the Cell
  papers and Niederhagen's thesis, against today's), new §6 on `clmad`
  and the integer pipe across dies, a per-GPU matrix splitting σ-walk
  and table walk, a σ-walk-versus-table-walk subsection (the 20.08 table
  rate ran superseded rule v2; exact v3 is 5.10), and a new §14 on open
  questions. Every rate is now in billions of iterations per second.
  Bibliography: INDOCRYPT 2010 GPU paper, Niederhagen thesis,
  Wenger–Wolfger SAC 2014 and JCE 2016 added; SHARCS 2009 author list
  and pages corrected; CRYPTO 2009 volume added.
- Draft 5: a prose pass over every section for readability (first
  person plural, shorter sentences, fewer stock phrases; every number,
  footnote, table and label kept), and a new §7 on the NIST curves (replacing the Draft 4 NIST subsection of §8):
  the `gpu/nist` rho walk's first device runs on the RTX PRO 6000
  (`gpu/nist/RESULTS.md`), set against the ECC2K-130 levers they
  borrow from. The P-256/P-384 row of the other-GPU inventory moves
  from "never" to "measured", and the abstract, introduction,
  related work, conclusion and reproduction appendix are updated to
  match. Graph and report copies are in
  `benchmarks/20261013-nist-gpu-rho-ladder/`.

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
- The NIST round-2 rows (W=32, register caps; 1.632 / 0.548 Gstep/s,
  1.52× / 2.61× cumulative) come from aburan28/crypto PR 1778, which
  was open when Draft 4 was written; round 1 (1.340 / 0.438) is on
  `main`. The NIST sweep cells are single runs, not interleaved
  repeated pairs, and the paper says so in the table caption. A
  Solinas-field walk has not been built; only the field benchmark
  compares the two representations.
- The author block is anonymized; `make deanonymize` with an `AUTHORS`
  file restores it.
