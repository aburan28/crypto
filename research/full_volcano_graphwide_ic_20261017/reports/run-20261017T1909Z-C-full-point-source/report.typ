#set page(paper:"a4",margin:13mm)
#set text(size:9pt)
= Original C: complete selected point orders
Run source. Started 2026-10-17T19:09:32.476381+00:00.
Host: INTEL(R) XEON(R) PLATINUM 8573C, x86_64, kernel 6.18.55.
Population: original source model only. Primary observations; replays are kept separately.

#image("pipeline.svg",width:100%)

Native accounting: 4 method records; 2 matched verified targets.
Point orders M range from 545986327049404416 to 545986327049404416.
Both arms recover modulo the complete verified M. A selected point need not
generate the full curve group; maximality is not asserted.

The field is p=1637958981148213249, and N=p-1=2^26*3^20*7.
Every selected point order is smooth over 2,3,7. Its large numeric value is
not a large-prime security level. Both arms exploit the known factorization.
Solver arguments contain the public curve, G, Q and M, never scalar truth.
The runner verifies the recovered point and the separate public known answer.

The finite-field route varies the first argument of reduced Tate/Miller,
checks exact prime-power image order, then performs multiplicative PH.
Initially unknown projected small-prime base logs are recovered by checked
one-/two-factor relations and modular Gaussian elimination. No generic-log
fallback is used. This is a smooth-order diagnostic, not NFS/integer-smoothness
index calculus, a large-k backend or a new asymptotic result.

Recorded finite-field setup equations/queries:
14 / 14.
Target digit equations/queries:
92 / 92.
Identity and dependent setup equations remain counted; per-component ranks,
matrix work and pair tests are retained. Tiny prime digit groups explain the
dense relation coverage; these yields do not describe a large-order oracle.

Fixture construction, exact order finding and known-target generation are
recorded separately, outside both solver intervals. Cold solver work includes
input-order checks, projections, method setup, auxiliary searches, pairing,
relations/matrix work, digit recovery, CRT and final EC point verification.
Raw end-to-end observations remain in the native CSV, with min/max/mean and n.
No admitted runtime comparison is made: one round, no A/A, PSI unavailable.
CPU 4 was reserved; wrapper observed contended=false. Replays do not
contribute to primary means. The full-prefix replay is explicitly contended.

The immutable C prefix contains 1,869 vertices and 4,903 replayed maps;
1,024 vertices are expanded. Local labels certify 517 vertices / 24 observed
strata and keep 1,352 unknown. Original full C remains 1,934,901,248 vertices
and 165 strata. This does not complete original A--E, the full eight-target /
five-round matrix, transported-target comparisons, representation registration
or independent arithmetic replay.

Raw records, fixtures, scalar truth, native CSV/statistics, isolation journals,
exact binaries and both source/replay populations accompany the delivery.
Native implementation is source-pinned at d954e80d9; linked-library source
2bf7700b7. Source hash #text(size:7pt)[#raw("bc90c72be7042e6ff355ed529302f2fcbc48d5e0e4ce15c4bd6dd241cd179234")].

#pagebreak()
= Recorded primitive work across local strata
#image("operations.svg",width:100%)
#set text(size:7pt)
#table(columns:(1.2fr,1.8fr,.5fr,.9fr,1fr,.9fr),inset:3pt,
[Arm], [Metric], [n], [Min], [Mean], [Max],
[EC PH / BSGS], [Curve add + double], [2], [1,619], [1,622.00], [1,625],
[EC PH / BSGS], [Explicit Fp mul + square], [2], [0], [0.00], [0],
[EC PH / BSGS], [Curve-field inverse], [2], [1,561], [1,569.00], [1,577],
[Tate + FF PH relations], [Curve add + double], [2], [993], [997.00], [1,001],
[Tate + FF PH relations], [Explicit Fp mul + square], [2], [72,622], [72,881.00], [73,140],
[Tate + FF PH relations], [Curve-field inverse], [2], [992], [996.00], [1,000],)
#set text(size:8pt)
These are native descriptive summaries. Curve additions/doublings, curve-field
inversions and explicit Fp multiplication/square counts are different units;
they are not a priced total. EC inversion and square-root internals remain
unexpanded. Explicit Fp inversions already contribute to the multiplication /
square ledger and must not be priced twice. Allocation, comparisons, scalar
integer work and Gaussian inverse details remain unpriced with counters where
available. All per-phase counters and present-value n are in the native CSV.
Observed ranges are not confidence intervals. Unknown strata remain separate.
Source and replay plots may match because every non-clock field reproduces;
they are not extra independent samples or an independent arithmetic proof.
