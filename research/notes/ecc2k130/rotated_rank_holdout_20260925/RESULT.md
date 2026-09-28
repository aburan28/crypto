# Exact modular rank and point-only holdouts on frozen rotated toy rows

**Decision: PASS for the bounded n=13 rank/recovery semantics.** The
preregistered source/input head `c1ae12eec1469241da560cd329184e1843476915`
selected 16 projected-positive nonzero subgroup Q points per arm by a fixed
SHA-256 order, preserving the 4/0/2/0 preceding negative candidates. It
removed **every** archived Q+T row for each holdout Q before exact
modulo-2003 elimination. The remaining #766 certified rows have rank equal
to their full 2/2/2/3 canonical term-column sets, with first full rank at
training row 3/3/3/5. The vast majority of valid group identities are
linearly dependent. The independent bit-serial/Fermat replay reconstructed
the split and rank, checked all training equations and actual group points for
each solved base log, rebuilt the archive-oracle lookup, then verified each
recovered scalar by `[k]H=Q` before opening sealed k labels.

| Rotated arm | Frozen positive rows | Excluded for 16 Q | Retained rows | Zero rows | Dependent incl. zero | Independent / columns | First full rank | Point-only recovered |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| beta 3, m5 | 1,799 | 19 | 1,780 | 1 | 1,778 | 2 / 2 | 3 | 16 / 16 |
| beta 3, m6 | 6,307 | 46 | 6,261 | 1 | 6,259 | 2 / 2 | 3 | 16 / 16 |
| beta 7, m5 | 1,799 | 18 | 1,781 | 1 | 1,779 | 2 / 2 | 3 | 16 / 16 |
| beta 7, m6 | 8,012 | 64 | 7,948 | 1 | 7,945 | 3 / 3 | 5 | 16 / 16 |

The rank denominator is the set of columns **appearing in archived terms**
for that arm, fixed before the holdout exclusion. The per-arm ranks cannot
be added into one global matrix rank because their column sets differ. All
64 holdout Qs were selected for positive support by rule; no natural-target
hit rate follows. The point-only child was passed Q coordinates, solved
column logs and an oracle response with only signed coefficients and a
torsion index. It was not passed k or witness indices; this is audited
source/argument isolation, not OS filesystem confinement. The archive-oracle
builder scanned all 17,917 old row records (26,194,782 uncompressed bytes)
and selected the first saved T-index row for each Q. Its 64/64 hits are
**precomputed witness lookups**, not a fresh implicit PDP solver or new
relation yield. The original #762 support census and #766 row generation
remain prerequisite work outside this tiny rank-stage receipt.

The frozen local run occupied 12:38:13–12:38:15 UTC. Input construction took
0.051 s wall and peaked at 33,587,200 bytes RSS. The producer's complete
training, oracle scan and point-only recovery took 0.681 s wall: child stages
were 0.429, 0.211 and 0.041 s. Independent sealed-label replay took 0.441 s.
The rank trainer's own peak RSS was 176,963,584 bytes, the oracle's was
80,134,144 bytes, and the verifier's was 181,403,648 bytes. Exact modular,
field and curve operation counters and all stage receipts are archived.
These are local-host toy costs and do not include a new PDP search or a
same-workload rho reference, so there is no attack-speed ratio.

This answers the immediate ambiguity left by #766: 17,917 correct n=13 row
identities provide only 2/2/2/3 independent equations in their respective
tiny column spaces. It also checks that, *given* a saved decomposition row,
solved column logs recover held-out scalars by the signed Frobenius row rule.
The next n=131 admission needs exact or bounded canonical projected-base
cardinality after sign/torsion duplicates, a solver that finds new rows for
natural targets, measured independent-rank yield, charged memory and a
matched automorphism-aware rho comparison. The n=131 SHA-fixed planted rows
in #766 do not supply any of those quantities.
