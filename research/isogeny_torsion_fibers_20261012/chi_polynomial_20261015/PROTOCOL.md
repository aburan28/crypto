# Frozen chi-polynomial decomposition gate (2026-10-15)

Status: preregistered before running this implementation. This gate extends the
[extension-field protocol](../extension_field_20261014/PROTOCOL.md) and its
[next-test list](../pairing_transfer_20261015/NEXT_TESTS.md).

## Inputs and boundary

Keep the four exact `FIXTURES.jsonl` rows, SHA-256
`1a63cc22a22e0fe27d44359325d502038d26d3ac616374fc5558c041939218b7`:
`(p,r)=(61,19),(271,61),(787,211),(3319,811)`. Keep their public `G,Q`,
eight/eight/fifteen/twenty-nine projected columns, original line, cofactor,
seed `0x20261014 XOR p` plus `round * 0x9e3779b9`, three rounds, and
`min(1_000_000,20r)` relation-query cap. The known scalar is available only
to a separate verifier. The generic two-class floor is
`S=sqrt(pi/4)=0.886226925`. Same-target negation and six-class rho are the
references; the earlier V4 pair-table and chi-index results are controls.

## Polynomial and solution contract

Use `z=chi(P)=(y+1)/(y-1)`, `d=z-1`, `n=z+1` on `E:y^2=x^3+1` over the
specified `F_(p^2)`. The explicit polynomial constraints for each summand are
`x^3*d^2-4z=0` and the frozen factor-base root polynomial in `z`, with
the factor-base coordinate/lift filter. For a finite query `R=(x_R,y_R)`,
the additional polynomial is the cleared collinearity determinant for
`P_1,P_2,-R`. Equal summands use the cleared tangent equation. Infinity,
`T=(0,1)` and `2T=(0,-1)` are separate branches. Solve the finite system
without a precomputed pair-sum table or chi-key lookup. Enumerative
root-solving is allowed but must charge every candidate polynomial
evaluation and lift; it is not an optimized resultant solver. Translate any
canonical chi-orbit representative back by its recorded torsion correction,
then verify each returned pair by full source-curve addition.

First compare the **exact** returned pair set, for every point of the order-19
subgroup, with an independent native enumeration of all signed projected
factor-base pairs. Abort DLP timing on a mismatch and retain the mismatch.
Then run the four frozen target solves, plus direct `S3`, quotient-invariant,
projected-pair, chi-index, matched rho and A/A controls as implemented. Keep
same target and source query streams, equal memory cap, and one-target cold
setup accounting. A run stops after a verified scalar, `10^6` relation
queries, `10^7` polynomial/root operations, or the frozen memory cap.

Record construction, unsuccessful evaluations, roots/lifts, curve operations,
relation collection, rank, recovery, independent scalar replay, peak RSS,
thread CPU and wall time, including failed attempts. `S` is complete only if
all native work has a calibrated conversion; otherwise report whole-process
CPU and wall separately as L0 shared-host descriptive data. A four-size
scaling claim needs L2 isolation, interleaving, A/A admission and independent
host replay. Preserve all earlier runs as historical records.
