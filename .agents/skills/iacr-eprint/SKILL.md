---
name: iacr-eprint
description: Prepare a LaTeX paper for the IACR Cryptology ePrint Archive - restore the author block, build and check the PDF, generate the form-ready abstract and a field-by-field EPRINT.md for https://eprint.iacr.org/submit, and enforce the archive's submission rules (no anonymous papers, author emails in the PDF, abstract and keyword format, A4/Letter, license). Use when asked to publish, post, submit or revise a paper on ePrint / IACR eprint. Stops before the submit button; an author submits.
---

# Publishing a paper on IACR ePrint

ePrint has no account and no peer review: one web form, a PDF upload, and
an editor's accept/reject by email. Acceptance mails a **password**, the
only way to revise or withdraw. This skill takes a paper from the repository
to a filled-in form; it never submits.

## Rules the archive enforces

From https://eprint.iacr.org/operations.html and the form at
https://eprint.iacr.org/submit (re-read both if anything below seems off;
the script's header records the date they were last read).

1. **No anonymous papers.** Page 1 shows the title, the author names and a
   contact address or affiliation, and the PDF contains the author email(s).
   Papers in this repository are built anonymized for double-blind venues,
   so the author block must be restored first.
2. **Only the authors submit, and all of them must agree.** The license
   granted to IACR is non-exclusive and **irrevocable**; an accepted paper
   can never be fully removed (withdrawal keeps title, abstract and old
   versions, and a withdrawn paper cannot be resubmitted). Never press
   Submit for the user, never submit without every named author's consent,
   and never invent an author, affiliation or email: ask.
3. **Scope and claims.** Cryptology, self-contained, "proofs or convincing
   arguments for any claims". Repository claim hygiene applies doubly here:
   no key recovery, sub-rho asymptotic or deployed-curve weakness unless the
   evidence gates in AGENTS.md were met.
4. **Abstract**: at least 64 characters, printable UTF-8 only (text pasted
   from a PDF often is not). LaTeX *math* in `$…$` renders with MathJax;
   text commands (`\emph`, `\texttt`, `\cite`, `\ref`, footnotes) are
   rejected; HTML is not rendered. Paper macros (`\FF`) must be expanded
   (`\mathbb{F}`). The generator unwraps formatting commands (`\texttt{clmad}`
   becomes `clmad`) and expands macros; it still fails on `\cite`, `\ref`,
   footnotes and the like, which have no plain-text form and must be
   reworded in the paper's abstract.
5. **Title**: plain UTF-8, math in `$…$` allowed, no HTML, avoid `&`.
6. **Authors**: in paper order, UTF-8 names (not TeX), at least one email;
   ORCID / ROR optional.
7. **Keywords**: comma-separated, each at most 40 characters, no LaTeX.
8. **Category** (one): Applications, Cryptographic protocols, Foundations,
   Implementation, Secret-key cryptography, Public-key cryptography,
   Attacks and cryptanalysis.
9. **License** (one): CC BY, CC BY-SA, CC BY-NC, CC BY-NC-SA, CC BY-NC-ND,
   CC0. Default to CC BY; check it does not conflict with a venue's
   copyright agreement.
10. **PDF**: fits A4 or US Letter. Posting a paper under review at an
    anonymous-submission venue is allowed.
11. **Revisions** go on the same ePrint number via
    https://eprint.iacr.org/revise, never as a new entry.

## Steps

1. **Authors.** Get names (order), affiliation and contact email from the
   user. Write the paper's git-ignored `AUTHORS` file (see the paper's
   `AUTHORS.example`): the `\author{…}` / `\date{…}` block plus a
   `%pdfauthor=Name, Name` line for the PDF metadata. Then
   `make deanonymize` (papers under `paper/` share this Makefile pattern).
   Commit the restored `\author` block; the `AUTHORS` file stays ignored.
2. **Build and check.** `make && make check`: no unresolved references or
   citations, no overfull boxes. Look at page 1 yourself.
3. **Generate the form material.**

   ```sh
   python3 .agents/skills/iacr-eprint/scripts/eprint_form.py paper/<name> \
       --category "Attacks and cryptanalysis" \
       --keywords "kw one, kw two, ..." --license "CC BY"
   ```

   It reads the title, authors and abstract (following `\input` from the
   abstract environment), writes `eprint-abstract.txt` (paste-ready) and
   `EPRINT.md` (every form field with its value, plus check results) into
   the paper directory, and **exits 1** on any hard-rule failure: anonymous
   author or page 1, no email in the PDF, rejected text commands, unbalanced
   `$`, short or non-printable abstract, keywords over 40 chars or with
   LaTeX, a page size that is neither A4 nor Letter. Fix every FAIL; read
   every warn (leftover commands outside math, author names not found
   verbatim on page 1, empty PDF author metadata).
4. **Bundle.** `make eprint` (PDF plus sources) for the record; the form
   itself takes only the PDF.
5. **Land it.** Commit the restored author block, rebuilt PDF,
   `EPRINT.md` and `eprint-abstract.txt`, open the PR ready for review, and
   stop (AGENTS.md "Default workflow").
6. **Hand over.** Tell the user the PDF and `EPRINT.md` paths and that an
   author submits at https://eprint.iacr.org/submit. After acceptance they
   keep the password and record the ePrint number in the paper's README.

## Choosing the category

Attacks on a scheme, curve or problem instance, including measured
campaigns: *Attacks and cryptanalysis*. Kernels, hardware or software
speed without an attack claim: *Implementation*. Pick one; editors may
move it.
