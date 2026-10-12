# Port exact F5 rank certificate to current mainline

## Scope and hypothesis

The row-basis rank candidate at `236205132` was measured on a branch that
predates mainline's exact GF(2) prefix/resume kernel. Port the selected
column certificate and eight-table rank-only kernel onto mainline
`51926e7d5` without deleting or changing the prefix/resume API or its
tests. Keep the default `matrix_f5_f2` RREF output unchanged; the sparse
original-row result remains explicit opt-in.

The required invariant is exact rank and canonical row space on every
frozen F5 case, including dependent smaller cases, and exact rank on the
GF(2) kernel's dense, sparse, zero, partial-word, and dependent matrices.
Default route and prefix/resume tests must still pass. A benchmark of
the earlier branch does not establish a result on this source snapshot.

## Measurements and stop rule

Build one release binary on the current Apple host and run frozen seeds
`0` and `badc0de1`, comparing selective echelon with the ported
row-basis path. Record all seven F5 cases, phase costs, counted work,
fingerprints, source/binary hashes, A/A noise and failures. This local
screen is nonpromoting. The requested further-2× result needs the
original four-seed isolated Linux x86-64 one-thread gate: every paired
complete-call median and exact bootstrap lower bound must exceed 2.00×
against selective echelon, with no smaller-case or two-thread regression.
Stop or revise the port if exactness fails; do not use the old branch's
timings as a substitute for a current-mainline measurement.

The comparison is a matrix-F5 solver-stage diagnostic, not a one-target
IC online or matched-rho measurement.
