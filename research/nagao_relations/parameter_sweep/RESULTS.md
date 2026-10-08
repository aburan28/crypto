# Joint dimension / decomposition-length results

Frozen run: `results_v1/`; source and inputs pinned in its provenance. This is an accounting/coverage pilot, not a calibrated ECDLP speedup.

1584 generic stage cells, 504 algebraic controls, 132 cold IC runs, 8 matched rho runs, 33 exact coverage cells. 86 IC scalars verified; 46 incomplete runs retained. All 4195 full-run decomposition calls were independently audited, including empty results.

The full previous benchmark replay retains 240 stage cells and 48 verified DLP runs; `regression_comparison.json` checks exact completed-baseline operation counters and 384 oracle calls.

## Exact coverage, including failures

Each probability covers **all affine targets** in the actual distinct-abscissa domain, not just successful sampled targets. M is the actual rational abscissa count; B is the number of distinct nonzero columns after cofactor projection. The paper expression is an intensity, not a probability above one.

| Profile | ell | k | M | B | Covered / affine | Exact coverage | Paper intensity | Poisson(actual mean) |
|:--|--:|--:|--:|--:|:--|--:|--:|--:|
| binary5 | 2 | 2 | 2 | 2 | 4/43 | 0.0930 | 0.2500 | 0.0888 |
| binary5 | 2 | 3 | 2 | 2 | 0/43 | 0.0000 | 0.3333 | 0.0000 |
| binary5 | 2 | 4 | 2 | 2 | 0/43 | 0.0000 | 0.3333 | 0.0000 |
| binary5 | 3 | 2 | 5 | 4 | 28/43 | 0.6512 | 1.0000 | 0.5671 |
| binary5 | 3 | 3 | 5 | 4 | 33/43 | 0.7674 | 2.6667 | 0.8126 |
| binary5 | 3 | 4 | 5 | 4 | 29/43 | 0.6744 | 5.3333 | 0.7943 |
| binary5 | 4 | 2 | 10 | 5 | 43/43 | 1.0000 | 4.0000 | 0.9808 |
| binary5 | 4 | 3 | 10 | 5 | 43/43 | 1.0000 | 21.3333 | 1.0000 |
| binary5 | 4 | 4 | 10 | 5 | 43/43 | 1.0000 | 85.3333 | 1.0000 |
| binary9 | 2 | 2 | 3 | 3 | 12/507 | 0.0237 | 0.0156 | 0.0234 |
| binary9 | 2 | 3 | 3 | 3 | 8/507 | 0.0158 | 0.0208 | 0.0157 |
| binary9 | 2 | 4 | 3 | 3 | 0/507 | 0.0000 | 0.0208 | 0.0000 |
| binary9 | 3 | 2 | 5 | 5 | 40/507 | 0.0789 | 0.0625 | 0.0759 |
| binary9 | 3 | 3 | 5 | 5 | 80/507 | 0.1578 | 0.1667 | 0.1460 |
| binary9 | 3 | 4 | 5 | 5 | 80/507 | 0.1578 | 0.3333 | 0.1460 |
| binary9 | 4 | 2 | 11 | 10 | 170/507 | 0.3353 | 0.2500 | 0.3520 |
| binary9 | 4 | 3 | 11 | 10 | 440/507 | 0.8679 | 1.3333 | 0.9251 |
| binary9 | 4 | 4 | 11 | 10 | 507/507 | 1.0000 | 5.3333 | 1.0000 |
| binary9 | 5 | 2 | 17 | 16 | 360/507 | 0.7101 | 1.0000 | 0.6580 |
| binary9 | 5 | 3 | 17 | 16 | 507/507 | 1.0000 | 10.6667 | 1.0000 |
| binary9 | 5 | 4 | 17 | 16 | 507/507 | 1.0000 | 85.3333 | 1.0000 |
| subfield8-k0 | 1 | 2 | 1 | 0 | 0/507 | 0.0000 | 0.0625 | 0.0000 |
| subfield8-k0 | 1 | 3 | 1 | 0 | 0/507 | 0.0000 | 0.1667 | 0.0000 |
| subfield8-k0 | 1 | 4 | 1 | 0 | 0/507 | 0.0000 | 0.3333 | 0.0000 |
| subfield8-k0 | 2 | 2 | 29 | 26 | 498/507 | 0.9822 | 4.0000 | 0.9589 |
| subfield8-k0 | 2 | 3 | 29 | 26 | 507/507 | 1.0000 | 85.3333 | 1.0000 |
| subfield8-k0 | 2 | 4 | 29 | 26 | 507/507 | 1.0000 | 1365.3333 | 1.0000 |
| subfield8-b22 | 1 | 2 | 5 | 0 | 11/467 | 0.0236 | 0.0625 | 0.0702 |
| subfield8-b22 | 1 | 3 | 5 | 0 | 11/467 | 0.0236 | 0.1667 | 0.0860 |
| subfield8-b22 | 1 | 4 | 5 | 0 | 1/467 | 0.0021 | 0.3333 | 0.0128 |
| subfield8-b22 | 2 | 2 | 25 | 6 | 437/467 | 0.9358 | 4.0000 | 0.9218 |
| subfield8-b22 | 2 | 3 | 25 | 6 | 467/467 | 1.0000 | 85.3333 | 1.0000 |
| subfield8-b22 | 2 | 4 | 25 | 6 | 467/467 | 1.0000 | 1365.3333 | 1.0000 |

## Full cold pipeline: separate units, all attempts

Each row sums **both fixed seeds**, including incomplete work. Compare cost only where both workloads complete. A/M/Q are binary-field API additions/multiplications/squarings; scalar A/M/I are separate modulo-subgroup operations. Inversions in the binary field are already expanded into field operations. No conversion prices Python control, lookups, field operations and scalar arithmetic together. Consequently all S/rho/floor ratios are N/A; do not turn the vector into a single operation score.

| Profile | ell/k | Variant | Verified | Field A | Field M | Field Q | Scalar A/M/I | Cold seconds | S | /rho | /floor | Class |
|:--|:--|:--|:--|--:|--:|--:|:--|--:|:--|:--|:--|:--|
| binary5 | 2/2 | enumerate-last | 2/2 | 7,011 | 14,321 | 532 | 36/32/6 | 0.026844 | N/A | N/A | N/A | Reference |
| binary5 | 2/2 | pair-mitm | 2/2 | 7,011 | 14,321 | 532 | 36/32/6 | 0.026735 | N/A | N/A | N/A | Engineering candidate |
| binary5 | 2/3 | enumerate-last | 0/2 | 941 | 1,785 | 532 | 0/0/0 | 0.003677 | N/A | N/A | N/A | Reference; incomplete |
| binary5 | 2/3 | pair-mitm | 0/2 | 997 | 1,881 | 532 | 0/0/0 | 0.003961 | N/A | N/A | N/A | Engineering candidate; incomplete |
| binary5 | 2/4 | enumerate-last | 0/2 | 941 | 1,785 | 532 | 0/0/0 | 0.003771 | N/A | N/A | N/A | Reference; incomplete |
| binary5 | 2/4 | pair-mitm | 0/2 | 997 | 1,881 | 532 | 0/0/0 | 0.004188 | N/A | N/A | N/A | Engineering candidate; incomplete |
| binary5 | 3/2 | enumerate-last | 2/2 | 2,281 | 4,466 | 596 | 76/92/10 | 0.008869 | N/A | N/A | N/A | Reference |
| binary5 | 3/2 | pair-mitm | 2/2 | 2,281 | 4,466 | 596 | 76/92/10 | 0.009090 | N/A | N/A | N/A | Engineering candidate |
| binary5 | 3/3 | enumerate-last | 2/2 | 3,480 | 6,361 | 596 | 97/107/10 | 0.013632 | N/A | N/A | N/A | Reference |
| binary5 | 3/3 | pair-mitm | 2/2 | 2,422 | 4,689 | 596 | 97/107/10 | 0.009323 | N/A | N/A | N/A | Engineering candidate |
| binary5 | 3/4 | enumerate-last | 2/2 | 6,837 | 11,920 | 596 | 102/102/10 | 0.024561 | N/A | N/A | N/A | Reference |
| binary5 | 3/4 | pair-mitm | 2/2 | 3,312 | 6,134 | 596 | 102/102/10 | 0.012582 | N/A | N/A | N/A | Engineering candidate |
| binary5 | 4/2 | enumerate-last | 2/2 | 3,577 | 7,488 | 716 | 452/434/12 | 0.015589 | N/A | N/A | N/A | Reference |
| binary5 | 4/2 | pair-mitm | 2/2 | 3,577 | 7,488 | 716 | 452/434/12 | 0.016892 | N/A | N/A | N/A | Engineering candidate |
| binary5 | 4/3 | enumerate-last | 2/2 | 2,872 | 5,695 | 716 | 245/260/12 | 0.011031 | N/A | N/A | N/A | Reference |
| binary5 | 4/3 | pair-mitm | 2/2 | 4,780 | 9,018 | 716 | 245/260/12 | 0.017079 | N/A | N/A | N/A | Engineering candidate |
| binary5 | 4/4 | enumerate-last | 2/2 | 3,543 | 7,034 | 716 | 420/392/12 | 0.018598 | N/A | N/A | N/A | Reference |
| binary5 | 4/4 | pair-mitm | 2/2 | 5,349 | 10,178 | 716 | 420/392/12 | 0.018818 | N/A | N/A | N/A | Engineering candidate |
| binary5 | — | rho | 2/2 | 1,113 | 2,158 | 486 | 46/2/2 | 0.004685 | N/A | N/A | N/A | Reference |
| binary9 | 2/2 | enumerate-last | 0/2 | 40,842 | 114,157 | 13,364 | 6/12/3 | 0.336095 | N/A | N/A | N/A | Reference; incomplete |
| binary9 | 2/2 | pair-mitm | 0/2 | 40,842 | 114,157 | 13,364 | 6/12/3 | 0.361970 | N/A | N/A | N/A | Engineering candidate; incomplete |
| binary9 | 2/3 | enumerate-last | 0/2 | 73,915 | 204,400 | 13,364 | 6/8/2 | 0.518121 | N/A | N/A | N/A | Reference; incomplete |
| binary9 | 2/3 | pair-mitm | 0/2 | 41,009 | 114,637 | 13,364 | 6/8/2 | 0.324749 | N/A | N/A | N/A | Engineering candidate; incomplete |
| binary9 | 2/4 | enumerate-last | 0/2 | 14,456 | 23,571 | 13,364 | 0/0/0 | 0.083719 | N/A | N/A | N/A | Reference; incomplete |
| binary9 | 2/4 | pair-mitm | 0/2 | 14,624 | 24,051 | 13,364 | 0/0/0 | 0.086530 | N/A | N/A | N/A | Engineering candidate; incomplete |
| binary9 | 3/2 | enumerate-last | 1/2 | 42,178 | 112,770 | 13,468 | 106/136/10 | 0.319376 | N/A | N/A | N/A | Reference; incomplete |
| binary9 | 3/2 | pair-mitm | 1/2 | 42,178 | 112,770 | 13,468 | 106/136/10 | 0.325358 | N/A | N/A | N/A | Engineering candidate; incomplete |
| binary9 | 3/3 | enumerate-last | 2/2 | 75,365 | 194,820 | 13,468 | 215/230/12 | 0.492769 | N/A | N/A | N/A | Reference |
| binary9 | 3/3 | pair-mitm | 2/2 | 28,517 | 68,711 | 13,468 | 215/230/12 | 0.201161 | N/A | N/A | N/A | Engineering candidate |
| binary9 | 3/4 | enumerate-last | 2/2 | 96,060 | 249,338 | 13,468 | 182/194/12 | 0.675562 | N/A | N/A | N/A | Reference |
| binary9 | 3/4 | pair-mitm | 2/2 | 35,355 | 82,054 | 13,468 | 182/194/12 | 0.245358 | N/A | N/A | N/A | Engineering candidate |
| binary9 | 4/2 | enumerate-last | 2/2 | 35,794 | 86,804 | 13,708 | 432/596/22 | 0.253034 | N/A | N/A | N/A | Reference |
| binary9 | 4/2 | pair-mitm | 2/2 | 35,794 | 86,804 | 13,708 | 432/596/22 | 0.253903 | N/A | N/A | N/A | Engineering candidate |
| binary9 | 4/3 | enumerate-last | 2/2 | 39,995 | 94,192 | 13,708 | 745/893/22 | 0.267920 | N/A | N/A | N/A | Reference |
| binary9 | 4/3 | pair-mitm | 2/2 | 21,184 | 44,407 | 13,708 | 745/893/22 | 0.143760 | N/A | N/A | N/A | Engineering candidate |
| binary9 | 4/4 | enumerate-last | 2/2 | 53,471 | 130,058 | 13,708 | 886/1014/22 | 0.380653 | N/A | N/A | N/A | Reference |
| binary9 | 4/4 | pair-mitm | 2/2 | 22,926 | 49,497 | 13,708 | 1623/1707/22 | 0.150142 | N/A | N/A | N/A | Engineering candidate |
| binary9 | 5/2 | enumerate-last | 2/2 | 35,893 | 86,287 | 14,092 | 2315/2739/34 | 0.249588 | N/A | N/A | N/A | Reference |
| binary9 | 5/2 | pair-mitm | 2/2 | 35,893 | 86,287 | 14,092 | 2315/2739/34 | 0.252956 | N/A | N/A | N/A | Engineering candidate |
| binary9 | 5/3 | enumerate-last | 2/2 | 39,022 | 92,871 | 14,092 | 3449/3861/34 | 0.254968 | N/A | N/A | N/A | Reference |
| binary9 | 5/3 | pair-mitm | 2/2 | 27,088 | 62,288 | 14,092 | 3259/3674/34 | 0.177028 | N/A | N/A | N/A | Engineering candidate |
| binary9 | 5/4 | enumerate-last | 2/2 | 36,837 | 87,021 | 14,092 | 4033/4405/34 | 0.255691 | N/A | N/A | N/A | Reference |
| binary9 | 5/4 | pair-mitm | 2/2 | 27,814 | 64,490 | 14,092 | 4593/4949/34 | 0.182233 | N/A | N/A | N/A | Engineering candidate |
| binary9 | — | rho | 2/2 | 15,552 | 27,342 | 13,262 | 138/58/2 | 0.100603 | N/A | N/A | N/A | Reference |
| subfield8-b22 | 1/2 | enumerate-last | 0/2 | 20,801 | 53,037 | 16,220 | 0/0/0 | 0.167704 | N/A | N/A | N/A | Reference; incomplete |
| subfield8-b22 | 1/2 | pair-mitm | 0/2 | 20,801 | 53,037 | 16,220 | 0/0/0 | 0.169201 | N/A | N/A | N/A | Engineering candidate; incomplete |
| subfield8-b22 | 1/3 | enumerate-last | 0/2 | 20,801 | 53,037 | 16,220 | 0/0/0 | 0.161158 | N/A | N/A | N/A | Reference; incomplete |
| subfield8-b22 | 1/3 | pair-mitm | 0/2 | 21,361 | 54,637 | 16,220 | 0/0/0 | 0.168195 | N/A | N/A | N/A | Engineering candidate; incomplete |
| subfield8-b22 | 1/4 | enumerate-last | 0/2 | 20,801 | 53,037 | 16,220 | 0/0/0 | 0.163472 | N/A | N/A | N/A | Reference; incomplete |
| subfield8-b22 | 1/4 | pair-mitm | 0/2 | 21,361 | 54,637 | 16,220 | 0/0/0 | 0.172987 | N/A | N/A | N/A | Engineering candidate; incomplete |
| subfield8-b22 | 2/2 | enumerate-last | 2/2 | 29,970 | 80,939 | 17,548 | 523/541/14 | 0.262011 | N/A | N/A | N/A | Reference |
| subfield8-b22 | 2/2 | pair-mitm | 2/2 | 29,970 | 80,939 | 17,548 | 523/541/14 | 0.266015 | N/A | N/A | N/A | Engineering candidate |
| subfield8-b22 | 2/3 | enumerate-last | 2/2 | 29,067 | 77,626 | 17,548 | 440/464/14 | 0.234732 | N/A | N/A | N/A | Reference |
| subfield8-b22 | 2/3 | pair-mitm | 2/2 | 42,191 | 115,785 | 17,548 | 440/464/14 | 0.340318 | N/A | N/A | N/A | Engineering candidate |
| subfield8-b22 | 2/4 | enumerate-last | 2/2 | 29,612 | 79,171 | 17,548 | 490/492/14 | 0.246695 | N/A | N/A | N/A | Reference |
| subfield8-b22 | 2/4 | pair-mitm | 2/2 | 42,479 | 116,632 | 17,548 | 490/492/14 | 0.349002 | N/A | N/A | N/A | Engineering candidate |
| subfield8-b22 | — | rho | 2/2 | 20,765 | 52,999 | 12,942 | 58/2/2 | 0.162447 | N/A | N/A | N/A | Reference |
| subfield8-k0 | 1/2 | enumerate-last | 0/2 | 14,454 | 23,375 | 16,476 | 0/0/0 | 0.088782 | N/A | N/A | N/A | Reference; incomplete |
| subfield8-k0 | 1/2 | pair-mitm | 0/2 | 14,454 | 23,375 | 16,476 | 0/0/0 | 0.090666 | N/A | N/A | N/A | Engineering candidate; incomplete |
| subfield8-k0 | 1/3 | enumerate-last | 0/2 | 14,454 | 23,375 | 16,476 | 0/0/0 | 0.089498 | N/A | N/A | N/A | Reference; incomplete |
| subfield8-k0 | 1/3 | pair-mitm | 0/2 | 14,454 | 23,375 | 16,476 | 0/0/0 | 0.093007 | N/A | N/A | N/A | Engineering candidate; incomplete |
| subfield8-k0 | 1/4 | enumerate-last | 0/2 | 14,454 | 23,375 | 16,476 | 0/0/0 | 0.091463 | N/A | N/A | N/A | Reference; incomplete |
| subfield8-k0 | 1/4 | pair-mitm | 0/2 | 14,454 | 23,375 | 16,476 | 0/0/0 | 0.088461 | N/A | N/A | N/A | Engineering candidate; incomplete |
| subfield8-k0 | 2/2 | enumerate-last | 2/2 | 32,181 | 78,199 | 17,932 | 8047/9236/54 | 0.249952 | N/A | N/A | N/A | Reference |
| subfield8-k0 | 2/2 | pair-mitm | 2/2 | 32,181 | 78,199 | 17,932 | 8047/9236/54 | 0.241134 | N/A | N/A | N/A | Engineering candidate |
| subfield8-k0 | 2/3 | enumerate-last | 2/2 | 64,174 | 173,570 | 17,932 | 45601/46388/54 | 0.497496 | N/A | N/A | N/A | Reference |
| subfield8-k0 | 2/3 | pair-mitm | 2/2 | 59,817 | 168,679 | 17,932 | 53286/53975/54 | 0.522557 | N/A | N/A | N/A | Engineering candidate |
| subfield8-k0 | 2/4 | enumerate-last | 2/2 | 79,704 | 217,587 | 17,932 | 41277/41744/54 | 0.611704 | N/A | N/A | N/A | Reference |
| subfield8-k0 | 2/4 | pair-mitm | 2/2 | 61,502 | 172,704 | 17,932 | 41304/41771/54 | 0.509566 | N/A | N/A | N/A | Engineering candidate |
| subfield8-k0 | — | rho | 2/2 | 15,552 | 27,342 | 13,262 | 138/58/2 | 0.094895 | N/A | N/A | N/A | Reference |

## Predeclared neighboring-parameter gate

For each directed neighboring change, both individual seeds must verify and each field component must not rise, with at least one decreasing. Scalar nonincrease is reported separately; neither flag establishes a common-unit end-to-end gain. All passing cases are listed; no selected winner is promoted as universal.

| Profile | Variant | Baseline ell/k | Candidate ell/k | Field-vector gate | Scalar nonincrease |
|:--|:--|:--|:--|:--|:--|
| binary5 | enumerate-last | 3/3 | 3/2 | Pass on both seeds | True |
| binary5 | enumerate-last | 3/4 | 3/3 | Pass on both seeds | False |
| binary5 | enumerate-last | 4/2 | 3/2 | Pass on both seeds | True |
| binary5 | enumerate-last | 4/4 | 4/3 | Pass on both seeds | True |
| binary5 | pair-mitm | 3/4 | 3/3 | Pass on both seeds | False |
| binary5 | pair-mitm | 4/2 | 3/2 | Pass on both seeds | True |
| binary5 | pair-mitm | 4/3 | 3/3 | Pass on both seeds | True |
| binary5 | pair-mitm | 4/3 | 4/2 | Pass on both seeds | False |
| binary5 | pair-mitm | 4/4 | 3/4 | Pass on both seeds | True |
| binary5 | pair-mitm | 4/4 | 4/3 | Pass on both seeds | True |
| binary9 | enumerate-last | 3/4 | 3/3 | Pass on both seeds | False |
| binary9 | enumerate-last | 4/4 | 4/3 | Pass on both seeds | True |
| binary9 | pair-mitm | 3/4 | 3/3 | Pass on both seeds | False |
| binary9 | pair-mitm | 4/2 | 4/3 | Pass on both seeds | False |
| binary9 | pair-mitm | 4/4 | 4/3 | Pass on both seeds | True |
| binary9 | pair-mitm | 5/2 | 5/3 | Pass on both seeds | False |
| binary9 | pair-mitm | 5/3 | 4/3 | Pass on both seeds | True |
| binary9 | pair-mitm | 5/4 | 4/4 | Pass on both seeds | True |
| subfield8-k0 | enumerate-last | 2/3 | 2/2 | Pass on both seeds | True |
| subfield8-k0 | pair-mitm | 2/3 | 2/2 | Pass on both seeds | True |
| subfield8-b22 | pair-mitm | 2/3 | 2/2 | Pass on both seeds | False |

## Retained limits

{
  "fewer admissible abscissae than summands": 12,
  "rank deficient": 10,
  "zero projected factor columns": 24
}

Single stage repetition; two tiny DLP seeds. Wall times are descriptive, not a statistically established runtime gain. Generic controls cover k=2,3,4 and GF(8) supports; cached pullback/S4/S3 are measured only at q=2,k=3. No general-k Nagao, GF(8) polynomial adapter, F4 expansion, large-prime variation, prime-field result or asymptotic exponent fit is claimed.

Detailed stage timings, operation vectors, paired cells, phase totals, ranks, certificates and all failed attempts are in `comparison.json` and `results_v1/raw.jsonl`. See `README.md` for the interpretation and rerun commands.
