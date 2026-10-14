# Accelerating Pollard Rho with GPU Field Arithmetic: ECC2K-130 on NVIDIA Blackwell

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

- §1 Introduction. §2 History: 2009–2012 rates against today.
- §3 Curve and field. §4 Walks. §5 Field arithmetic (with a figure of
  the packed carry-less product).
- §6 Instruction sets: NVIDIA (PRMT, SHF, LOP3, carry chains, PTX 9.3
  `clmad`), x86 (PCLMULQDQ, MULX/ADX, VPTERNLOGD, IFMA, VPCLMULQDQ/GFNI,
  AVX10.2, APX), Arm (PMULL, EOR3, SVE2), RISC-V (Zbc, Zvbc), with a
  timeline figure.
- §7 Carry-less multiply on NVIDIA GPUs: `clmad` rates on three dies.
- §8 GPU: floor and roofline, the ledger (staircase figure), speedup
  breakdown (no algorithmic change), kernel anatomy, step data flow and
  memory-layout figures, table walk vs. σ-walk (cycle-rule figure),
  GPU comparison.
- §9 CPU (with a per-ISA chart).
- §10 NIST curves. §11 Other GPU kernels. §12 FPGA. §13 Campaign.
- §14 Alternatives to rho. §15 Upper bounds (figure).
- §16 Time to solve: expected time and cost for every measured device at
  one, 100 and 1,000 devices, variance, the cutoff-32 detection delay,
  the campaign corpus credit, a fleet-size figure, and the NIST curves.
- §17 Open problems. §18 Related work. §19 Conclusion.
- Appendix A: reproduction (companion repository and Docker image).
- Appendix B: evidence records.

Measured numbers cite evidence records [E*n*] (Appendix B) instead of
file paths. `app-evidence.tex` lists them, and `evidence-map.json` beside
the source maps each record to the files it came from. The companion
repository `rho-blackwell`, which will hold a copy of each record under its
number together with the code and a Docker image, is assembled but **not
yet published**, and its Docker image has not yet been built or run on a
GPU; the paper says so in its availability paragraph and Appendix A.

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
- Draft 7: no file paths in the text: every footnoted path became a
  numbered evidence record (88 records, Appendix B), with the companion
  repository `rho-blackwell` holding copies and the Docker image that
  starts the ECC2K-130 walk. New §6 on instruction sets (GPU and CPU),
  §9 rewritten around CPU instruction sets (including the packed PMULL
  client), §15 upper bounds, and a "what changed and what did not"
  subsection. Thirteen new figures: history chart, ledger staircase,
  speedup breakdown, kernel anatomy, step data flow, carry-less product,
  memory layout, roofline, cycle rules, per-GPU throughput and per-dollar,
  CPU per-ISA rates, instruction timeline, upper bounds. Corrections:
  the Frobenius networks replaced repeated squarings (not basis
  conversions); 385,024 workers run as four waves of the 96,256 that fit
  at once. Consecutive citations now render separately, and long URLs
  break.
- Draft 8: new §16, the time to a solve. A table gives expected
  device-years, calendar days at 100 and 1,000 devices, and list-price
  cost for every measured rate (eleven GPUs, the RTX PRO 4500 bridge
  map, the F2 FPGA, one CPU core, the RTX PRO 6000 ceiling and the 2010
  GTX 295). The model adds rho's spread (median 0.94x, 90% 1.71x of the
  expected work), the detection delay at cutoff 32 (about 1.6 days on an
  RTX PRO 6000, whatever the fleet size) and the campaign corpus credit
  (0.949 of the work remains). A fleet-size figure is new. The NIST
  subsection puts P-256 at about 6e21 card-years. That corrects the
  earlier "about 1e30 seconds" to 1.9e29 here and in the NIST benchmark
  report. Abstract, introduction and conclusion updated.
  The NIST round-2 results are no longer described as unmerged: PR 1778
  merged with `gpu/nist` identical to the cited commit, so E49 now points
  at `gpu/nist/RESULTS.md` on `main`.
  Retitled "Accelerating Pollard Rho with GPU Field Arithmetic: ECC2K-130 on NVIDIA Blackwell" to lead with the
  arithmetic speedup.
  The exact v3 table-walk rate is updated from 5.10 to the 2026-10-12
  screen's 8.24 (new evidence record E89), and §8 now says the 20.08
  receipt ran rule v2.
  New §8 "Profiling" subsection (rates, static and dynamic SASS counts,
  pipe-rate probes, clock64 phase builds, where Nsight Compute did and did
  not read counters) and a "Use of AI" paragraph in the introduction;
  evidence records E90 (THROUGHPUT-29B, Nsight counters) and E91 (probes
  and SASS cost tools).
  Section and subsection headings shortened to terse labels (for example
  "NIST curves", "FPGA", "Table walk vs. σ-walk", "Walk constant"), with the remaining rhetorical questions
  and informal phrases in the text replaced.
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
  cycle-rule v2 cost projection (0.81–0.86× σ) is unmeasured. The 20.08
  build ran rule v2 (116 registers, no frame; v3 needs 128 and a 400-byte
  frame). Exact v3 measured 5.10 in its first paired panel and 8.24 in an
  unpaired 2026-10-12 screen of the same knob set (E89), whose raw logs are
  not in the tree; v3's walk constant is not yet measured.
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
  1.52× / 2.61× cumulative) come from aburan28/crypto PR 1778, merged
  on 2026-10-14 with `gpu/nist` unchanged from the commit the paper
  cited (c0fc85eb7); round 1 (1.340 / 0.438) was already on `main`. The
  NIST sweep cells are single runs, not interleaved
  repeated pairs, and the paper says so in the table caption. A
  Solinas-field walk has not been built; only the field benchmark
  compares the two representations.
- The author block is anonymized; `make deanonymize` with an `AUTHORS`
  file restores it.
