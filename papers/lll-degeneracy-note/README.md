# LLL Failures on the Boneh–Venkatesan Lattice Were Floating-Point Artifacts: A Corrected Note

`lll_degeneracy_note.tex` corrects this directory's earlier note, "An Empirical
Note on LLL Degeneracy of the Boneh–Venkatesan Lattice for secp256k1" (PDF of
2026-10-11).

**Status**: corrected 2026-10-16. Research note, not peer-reviewed, not
submitted, not posted. Anonymized; restore the authors with `make deanonymize`
after creating an `AUTHORS` file (see `../AUTHORS.example`).

## Correction

The earlier version reported two failure modes as curve properties:

- secp256k1 failing LLL at 256 bits where P-256 and brainpoolP256r1 succeed;
- a uniform failure at 384 bits and above, blamed on the iteration cap.

Both were bugs in the `f64` Gram–Schmidt step, which
`research/notes/lattice-hnp/RESEARCH_LLL_GS_ANALYSIS.md` §10 had already
diagnosed and fixed in May 2026:

- the squared norms overflow (`n⁴ > f64::MAX`) for secp256k1 and for every
  384-bit curve;
- P-521 additionally suffers catastrophic cancellation.

A sheriff audit found the stale claim: crypto-autoresearcher
`research/sheriff/SHF-20261016-3caa6f`, blocker P1. The curve-specificity
claim and everything built on it are withdrawn. The pre-fix tables are kept,
labelled as superseded.

## What it says

- The two failure modes and their numerical causes.
- A 2026-10-16 re-run of the committed probes on the fixed code
  (`research/lll_degeneracy_rerun_20261016/`). Every probe at 256 and 384 bits
  recovers the key, secp256k1 included (6/6). P-521 recovers 3/3 with 2048-bit
  fixed-point Gram–Schmidt.

## What it does not say

- No security claim about any curve.
- No LLL cost comparison between curves.
- No statement beyond `m = 8`, `k_bits = 3/4 · n_bits`.

## Building

```
make           # PDF (pdflatex, no bibliography)
make check     # unresolved refs, LaTeX errors, overfull boxes, page count
make clean
```

The rules are the shared `../paper.mk`. The built PDF is committed.

## Before submission

1. Decide whether a correction note is worth posting at all. The earlier
   version was never posted.
2. Novelty audit: the `f64` overflow of LLL Gram–Schmidt at large entries is a
   known numerical hazard. Cite the standard treatments (Nguyen–Stehlé's
   floating-point LLL) before claiming any lesson as ours.
3. Measurements owed: more seeds and larger `m` if a cost statement is ever
   wanted.
