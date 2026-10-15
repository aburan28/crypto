#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= Native pairing batch V3/V4 phase ledger
*Runs:* V3 on order 19 and 61, V4 on order 67; three rounds each. *Host:* shared Darwin arm64, L0. *Specs:* `RUN_SPEC_V3.json`, `RUN_SPEC_V4_401.json`.

All target scalars agreed with charged negation rho and passed independent source-curve replay: 54/54 on order 19, 180/180 on order 61 and 198/198 on order 67. The degree-two order-67 pairing cell supplied 66 distinct nonidentity targets, including the requested 64-target prefix.

Median coarse pairing/rho CPU ratios at the stated prefixes were 411.05 for `(p,r,N)=(61,19,16)`, 37.59 for `(271,61,60)` and 3.42 for `(401,67,64)`. Six-class rho ratios were 199.31 and 31.00 on the first two cells; the automorphism does not preserve the prime-field order-67 subgroup.

#image("GRAPH.svg", width: 100%)

At order 67, the coarse process residual includes per-target JSON reporting, while the phase sum omits input and between-phase work. Their gap motivates the frozen V5 compact run. Both are descriptive L0 diagnostics. Raw SHA-256 values: order-19 V3 `f77508edb018433f04cdbc0fd059c725785cf34077cb5db8c8a1dc11088a1b1b`; order-61 V3 `35c1d0d4e4c34c0ba0f9479ae2c95e7b7027a2840811e749e301652c9cb1af0f`; order-67 V4 `c4ac091859c4d09dd40cdb462c46d38255007f41beef930d362dbe060e406c5f`. See `REPORT.md` for raw links and limits.
