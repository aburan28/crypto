# Pre-submission Checklist (IACR ePrint / conference)

A self-review checklist before submitting `structural_completeness.tex`
to IACR ePrint (https://eprint.iacr.org/) or a conference venue.

## Before submission

### Content review

- [ ] Read the paper end-to-end one more time
- [ ] Verify Theorem 8.1 (Structural Completeness) is stated precisely
      and the seven blocks (B1–B7) are individually correct
- [ ] Cross-check Section 6 (Howe-gluing): the errata about the naive
      construction (`y² = f_1·f_2` not being the Howe-glue) is
      properly stated
- [ ] Cross-check Section 7 (Empirical verification): the cross-curve
      table values match the actual `cargo test` output
- [ ] Verify all references compile (no `??` in PDF output)

### Build verification

- [ ] `pdflatex structural_completeness.tex` produces no errors
- [ ] Second pass resolves all forward references
- [ ] PDF page count is in the expected range (12–18 pages for a
      research note this length)
- [ ] Tables render correctly (no overfull/underfull boxes)

### Anonymization (if double-blind)

- [ ] Remove author block (`\author{...}`) — replace with "Anonymous"
- [ ] Remove the "audit programme conducted on the `crypto`
      repository" attribution
- [ ] Search for self-citations and rephrase to avoid de-anonymising

### Companion artefact disclosure

- [ ] Link to the public repository in the abstract / introduction
- [ ] Verify all PARI scripts are committed and reproducible
- [ ] Verify the Rust `cargo test --test curve_audit` runs cleanly
      and produces the table in §7
- [ ] Add a "How to verify" section pointing to specific scripts

### Citations

- [ ] Replace placeholder bibliography entries with proper venues
      and DOIs
- [ ] Verify Howe (1996), Tate (1966), Diem (2011) are the correct
      primary references for the blocks they support
- [ ] Add a citation for the Mestre algorithm (currently §6 references
      it without a primary citation in the running text)
- [ ] Add citation for Gallant-Lambert-Vanstone 2001 (GLV) since
      Section 5 discusses VFCG-ρ at j=0 curves

### Submission package

For IACR ePrint:

- [ ] Compile to PDF; archive PDF + .tex sources + bibliography
- [ ] Fill the ePrint submission form:
      - Title: from `abstract.txt`
      - Authors: real names + affiliations
      - Abstract: from `abstract.txt`
      - Keywords: from `abstract.txt`
      - Category: `Implementation` or `Foundations` (likely
        `Foundations` for the theorem, `Implementation` for the
        audit infrastructure)
- [ ] Choose appropriate ePrint metadata (license, version)

For a conference (e.g., AsiaCrypt, CRYPTO, ANTS):

- [ ] Format to the conference's LaTeX template (LNCS, IEEE, etc.)
- [ ] Trim or expand to the target page limit
- [ ] Add anonymization metadata if double-blind

## Suggested venues (ranked by fit)

1. **IACR ePrint** — fastest, no peer review.  Preprint server,
   universal acceptance, gets DOI.  Recommended FIRST step regardless
   of subsequent venue.

2. **ANTS (Algorithmic Number Theory Symposium)** — best fit for
   the audit-style results.  Welcomes empirical work with
   computational artefacts.  Biennial, deadline ~6 months before
   conference.

3. **Designs, Codes and Cryptography** — journal venue, long-form,
   accepts comprehensive audits.

4. **Journal of Cryptology** — journal venue, more theoretical
   slant; the Theorem 8.1 itself could anchor a J. Crypto. paper
   with the empirical work as supplementary material.

5. **CRYPTO/AsiaCrypt** — top conferences; competition would be
   stiff.  A "structural-completeness theorem" framing might
   resonate but the work is more incremental than these venues
   prefer.

6. **Cryptanalysis Workshops** (CHES, FDTC, COSADE) — better fit
   for the side-channel/HNP companion direction
   (`RESEARCH_GLV_HNP.md`) than for the main structural-completeness
   theorem.

## Risks and weaknesses to address

1. **Limited novelty**: each block (B1–B7) is individually well-
   known.  The novelty is the *combination* into a closed-form
   theorem.  Some reviewers may see this as "expository" rather
   than "new research".  Counter-argument: nobody has previously
   stated the closure as a theorem, and the empirical verification
   (§7) is genuinely new.

2. **Limited scope**: the theorem covers ordinary prime-field
   curves; doesn't address supersingular curves, twisted curves
   over extension fields, or pairing-based protocols.  Add a
   "What's not in this theorem" section explicit about scope.

3. **No explicit Howe cover**: §6 ends with the existence statement
   only; constructive recovery is left open.  Some reviewers will
   want the explicit C.  Counter-argument: the theorem doesn't
   depend on the explicit construction (B5 applies to any cover);
   the constructive piece is a separate paper.

4. **Conditional on standing assumptions**: B5 depends on Gaudry's
   index-calculus bounds being tight.  If a future paper improves
   Gaudry's bounds (e.g., a sub-exponential genus-2 prime-field
   DLP), then B5 weakens.  Note this honestly.

## Post-submission

- [ ] Add a "Published in" line to README.md once accepted
- [ ] Update PAPER_STRUCTURAL_COMPLETENESS.md with a link to the
      ePrint / DOI
- [ ] Continue developing the companion directions
      (Mestre algorithm, GLV-HNP) as follow-up papers

## Co-author note

The work as committed is the audit programme.  Real-name
attribution depends on:

- Who actually executed and reviewed the work
- Acknowledgement of code-assistance from the audit conversation
- Disclosure of which structural blocks were independently checked
  by which author

Standard cryptanalysis-paper practice: state "AI-assisted research"
where relevant, with the human authors verifying every theorem and
empirical claim.
