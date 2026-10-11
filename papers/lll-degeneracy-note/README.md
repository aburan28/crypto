# An Empirical Note on LLL Degeneracy of the Boneh–Venkatesan Lattice for secp256k1

`lll_degeneracy_note.tex` is a short empirical note on how the Boneh–Venkatesan
hidden-number-problem lattice degenerates under LLL and BKZ on secp256k1 and
comparison curves, with the probe tables it reports.

**Status**: research note, not peer-reviewed, not submitted, not posted.
Anonymized; restore authors with `make deanonymize` after creating an `AUTHORS`
file (see `../AUTHORS.example`).

## What it says

Measured LLL/BKZ behaviour of the Boneh–Venkatesan lattice at the parameters in
its tables, on secp256k1 against comparison curves. Read the tables for the
exact probes; the note claims nothing beyond them.

## What it does not say

No key-recovery claim against any deployed system, no new lattice algorithm,
and no statement about parameters it did not probe.

## Building

```
make           # PDF (pdflatex, no bibliography)
make check     # unresolved refs, LaTeX errors, overfull boxes, page count
make clean
```

The rules are the shared `../paper.mk`. The built PDF is committed.

## Before submission

1. Re-run every probe in the tables from the companion scripts under
   `../../` named in the note and confirm the numbers.
2. Resolve the citation of the companion audit (the structural-completeness
   paper in `../structural-completeness/`) with its ePrint number once it has one.
3. Independent read of the interpretation of the degeneracy pattern.
