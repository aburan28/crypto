# Submitting to the IACR ePrint archive

ePrint (https://eprint.iacr.org/submit) is author-operated: there is no
account; the submitter fills one form, uploads the PDF, and receives the
decision and a revision password by email. Two of its rules bind this paper:

1. **No anonymous papers.** The first page must carry the title, the
   author names and a contact address or affiliation, and the PDF must
   include the author email address(es). The committed PDF is anonymized,
   so the author block has to be restored first (step 1 below).
2. **The license is irrevocable** and every named author must agree to
   the submission. Only an author should press Submit.

## 1. Restore the author block and rebuild

```sh
cd paper/ecc2k130
cp AUTHORS.example AUTHORS      # fill in names, affiliation, email;
                                 # keep the %pdfauthor= line (PDF metadata)
make deanonymize                 # rewrites \author{...} and pdfauthor; keeps a .anon-backup
make && make check               # PDF with authors; no bad refs / overfull boxes
make eprint                      # eprint-submission.zip (PDF + sources)
```

`AUTHORS` is git-ignored. Check the first page of the PDF before uploading.

## 2. The form

| field | value |
|---|---|
| Title | `Pollard Rho on NVIDIA Blackwell: ECC2K-130, Carry-less Multiplication and the NIST Curves` |
| Authors | as in `AUTHORS`, in paper order; at least one email; ORCID optional |
| Abstract | paste `eprint-abstract.txt` (UTF-8 with `$…$` math only; MathJax renders it; `\emph`, `\cite`, HTML are rejected) |
| Category | **Attacks and cryptanalysis** (alternative: Implementation) |
| Keywords | `ECC2K-130, Pollard rho, Blackwell, GPU, carry-less multiplication, Koblitz curves, NIST P-256, elliptic curve discrete logarithm, Certicom challenge, optimal normal basis, FPGA, binary fields` |
| Publication info | Published nowhere else (the repositories are not a publication) |
| Content | `ecc2k130.pdf` from step 1 (A4 or US Letter; author emails on page 1) |
| Additional notes | optional, public; leave empty for a first version |
| License | **CC BY** unless a venue agreement says otherwise |
| Confirmation | tick after reading the acceptance, license and withdrawal conditions |

Nothing in the paper claims a key recovery, a sub-rho asymptotic or a
weakness of a deployed curve, and the abstract says so; that is what the
editors' "convincing arguments for any claims" criterion will look for.

## 3. After acceptance

Keep the password from the acceptance email; revisions (`/revise`) go on
the same ePrint number, never as a new entry. Record the ePrint number in
`README.md` and in `cryptanalysis/paper/` where the companion cites this
paper.
