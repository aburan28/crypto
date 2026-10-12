# IC autolab development checkpoint — 2026-09-24

The portable loop completed six frozen, audited development screens. **888/960 jobs verified**; every rejection is retained. No performance winner is promoted. The first 192-job run completed its solves but failed during rho null-rank reporting; its original evaluator and receipts remain unchanged.

All jobs include curve/target generation, base construction, unsuccessful attempts, solving, relation processing, scalar-field linear algebra, descent, verification and output. The independent Python audit is outside the child. Operation counts, S and floor ratios are unmeasured, not zero. Seconds do not substitute for them.

## Sequence and findings

| Panel | Cells | Verified / jobs | A/A paired wall ratio |
|---|---|---:|---:|
| implementation_v2 | 13a0, 9a0 | 192 / 192 | 0.8546 |
| algebra | 9a0, 9a1 | 216 / 216 | 0.8363 |
| factor_base | 13a0, 9a0 | 156 / 204 | 0.9175 |
| union_followup | 13a0, 9a0 | 84 / 108 | 0.9468 |
| interactions | 13a0, 17a1 | 120 / 120 | 0.9142 |
| interaction_replication | 13a0, 17a1 | 120 / 120 | 0.8617 |

- The implementation grid tests batch/window/linear-algebra interactions as complete jobs.
- All seven PDP backends, each with dense and sparse relation linear algebra, passed the degree-9 panel.
- Initial two-generator Frobenius unions failed the existing nontrivial-matrix admission rule. That is a protocol rejection, not evidence that a one-column IC method is impossible.
- Enlarging the union to seed masks `[1,2,4,8]` passed both cells. Other unions remained partly or wholly rejected.
- That union, its batch-1 and dense-LA combinations, and their individual parents passed the degree-13/17 interaction panel and its separate-seed development replication.
- Rho controls request 1, 8 and 32 walks; the implementation can cap effective widths on tiny groups. The strongest general reference is not established by this small timing panel.
- The native A/A controls expose uncontrolled timing variation. The local rho test overlap in `implementation_v2` is recorded explicitly. No native speedup claim follows from these screens.

## Complete static measurements

One unit: complete cold child milliseconds. Wall-time columns are diagnostics; ratios pair identical cases and weight curve cells equally. Repetitions are not new targets. Each rho ratio is reported separately, avoiding a post-hoc claim about the best reference. Incomplete arms have no aggregate cost or ratio. Their per-case admitted ranks and failure reasons remain in `RESULTS.json`.

### implementation_v2

| Variant | Cold ms | Paired / incumbent | / rho 1 | / rho 8 | / rho 32 | Verified |
|---|---:|---:|---:|---:|---:|---:|
| incumbent | 21.8719 | 1.0000 | 3.0034 | 2.8149 | 2.9752 | 12/12 |
| implementation_01 | 14.0861 | 0.5871 | 1.7632 | 1.6526 | 1.7467 | 12/12 |
| implementation_02 | 11.3661 | 0.4663 | 1.4006 | 1.3127 | 1.3874 | 12/12 |
| implementation_03 | 17.1487 | 0.5512 | 1.6556 | 1.5517 | 1.6401 | 12/12 |
| implementation_04 | 14.0995 | 0.6559 | 1.9698 | 1.8462 | 1.9513 | 12/12 |
| implementation_05 | 19.3113 | 0.7642 | 2.2953 | 2.1513 | 2.2738 | 12/12 |
| implementation_06 | 14.7324 | 0.6011 | 1.8053 | 1.6920 | 1.7884 | 12/12 |
| implementation_07 | 25.8580 | 1.0254 | 3.0798 | 2.8865 | 3.0509 | 12/12 |
| implementation_08 | 12.4010 | 0.5899 | 1.7718 | 1.6606 | 1.7551 | 12/12 |
| implementation_09 | 15.7062 | 0.6501 | 1.9525 | 1.8299 | 1.9341 | 12/12 |
| implementation_10 | 35.1935 | 1.1015 | 3.3083 | 3.1006 | 3.2772 | 12/12 |
| implementation_11 | 18.5262 | 0.8513 | 2.5567 | 2.3962 | 2.5327 | 12/12 |
| aa_control | 19.7760 | 0.8546 | 2.5668 | 2.4057 | 2.5427 | 12/12 |
| rho_w1 | 8.0088 | 0.3330 | 1.0000 | 0.9372 | 0.9906 | 12/12 |
| rho_w8 | 8.4347 | 0.3553 | 1.0670 | 1.0000 | 1.0569 | 12/12 |
| rho_w32 | 7.9121 | 0.3361 | 1.0095 | 0.9461 | 1.0000 | 12/12 |

### algebra

| Variant | Cold ms | Paired / incumbent | / rho 1 | / rho 8 | / rho 32 | Verified |
|---|---:|---:|---:|---:|---:|---:|
| incumbent | 17.5033 | 1.0000 | 2.4006 | 2.8468 | 2.0742 | 12/12 |
| algebra_01 | 14.9381 | 0.8749 | 2.1003 | 2.4905 | 1.8147 | 12/12 |
| algebra_02 | 14.6655 | 0.8676 | 2.0828 | 2.4699 | 1.7996 | 12/12 |
| algebra_03 | 14.4003 | 0.9063 | 2.1758 | 2.5801 | 1.8799 | 12/12 |
| algebra_04 | 26.4680 | 1.5454 | 3.7100 | 4.3994 | 3.2055 | 12/12 |
| algebra_05 | 25.9082 | 1.8775 | 4.5072 | 5.3447 | 3.8943 | 12/12 |
| algebra_06 | 25.2754 | 1.1423 | 2.7424 | 3.2520 | 2.3694 | 12/12 |
| algebra_07 | 26.2983 | 1.5179 | 3.6440 | 4.3212 | 3.1485 | 12/12 |
| algebra_08 | 20.8496 | 1.1332 | 2.7205 | 3.2260 | 2.3505 | 12/12 |
| algebra_09 | 21.6427 | 1.2295 | 2.9515 | 3.5000 | 2.5502 | 12/12 |
| algebra_10 | 18.0390 | 1.1084 | 2.6608 | 3.1552 | 2.2989 | 12/12 |
| algebra_11 | 16.4471 | 0.9716 | 2.3326 | 2.7660 | 2.0154 | 12/12 |
| algebra_12 | 33.6458 | 2.0596 | 4.9444 | 5.8632 | 4.2720 | 12/12 |
| algebra_13 | 27.8681 | 1.5900 | 3.8170 | 4.5263 | 3.2979 | 12/12 |
| aa_control | 12.9528 | 0.8363 | 2.0078 | 2.3809 | 1.7347 | 12/12 |
| rho_w1 | 7.4675 | 0.4166 | 1.0000 | 1.1858 | 0.8640 | 12/12 |
| rho_w8 | 6.1080 | 0.3513 | 0.8433 | 1.0000 | 0.7286 | 12/12 |
| rho_w32 | 7.3074 | 0.4821 | 1.1574 | 1.3725 | 1.0000 | 12/12 |

### factor_base

| Variant | Cold ms | Paired / incumbent | / rho 1 | / rho 8 | / rho 32 | Verified |
|---|---:|---:|---:|---:|---:|---:|
| incumbent | 23.3162 | 1.0000 | 3.0396 | 2.7568 | 2.9965 | 12/12 |
| factor_base_01 | 29.1259 | 1.2987 | 3.9474 | 3.5801 | 3.8915 | 12/12 |
| factor_base_02 | 28.2757 | 1.2220 | 3.7142 | 3.3687 | 3.6616 | 12/12 |
| factor_base_03 | 26.6037 | 1.0804 | 3.2841 | 2.9786 | 3.2376 | 12/12 |
| factor_base_04 | 27.0244 | 1.1541 | 3.5080 | 3.1816 | 3.4583 | 12/12 |
| factor_base_05 | 29.1676 | 1.1071 | 3.3651 | 3.0520 | 3.3174 | 12/12 |
| factor_base_06 | 28.1824 | 1.2344 | 3.7521 | 3.4031 | 3.6990 | 12/12 |
| factor_base_07 | 642.5599 | 8.9549 | 27.2189 | 24.6867 | 26.8337 | 12/12 |
| factor_base_08 | 599.7235 | 6.7610 | 20.5505 | 18.6387 | 20.2597 | 12/12 |
| factor_base_09 | — | — | — | — | — | 0/12 |
| factor_base_10 | — | — | — | — | — | 0/12 |
| factor_base_11 | — | — | — | — | — | 0/12 |
| factor_base_12 | — | — | — | — | — | 0/12 |
| aa_control | 24.2085 | 0.9175 | 2.7886 | 2.5292 | 2.7492 | 12/12 |
| rho_w1 | 7.2721 | 0.3290 | 1.0000 | 0.9070 | 0.9858 | 12/12 |
| rho_w8 | 8.1195 | 0.3627 | 1.1026 | 1.0000 | 1.0870 | 12/12 |
| rho_w32 | 8.7662 | 0.3337 | 1.0144 | 0.9200 | 1.0000 | 12/12 |

### union_followup

| Variant | Cold ms | Paired / incumbent | / rho 1 | / rho 8 | / rho 32 | Verified |
|---|---:|---:|---:|---:|---:|---:|
| incumbent | 22.8911 | 1.0000 | 2.6607 | 2.8439 | 2.5425 | 12/12 |
| union_1 | — | — | — | — | — | 6/12 |
| union_2 | — | — | — | — | — | 0/12 |
| union_3 | — | — | — | — | — | 6/12 |
| union_4 | 16.6916 | 0.8995 | 2.3933 | 2.5581 | 2.2870 | 12/12 |
| aa_control | 19.3925 | 0.9468 | 2.5192 | 2.6926 | 2.4073 | 12/12 |
| rho_w1 | 8.0198 | 0.3758 | 1.0000 | 1.0689 | 0.9556 | 12/12 |
| rho_w8 | 6.9927 | 0.3516 | 0.9356 | 1.0000 | 0.8940 | 12/12 |
| rho_w32 | 9.1493 | 0.3933 | 1.0465 | 1.1185 | 1.0000 | 12/12 |

### interactions

| Variant | Cold ms | Paired / incumbent | / rho 1 | / rho 8 | / rho 32 | Verified |
|---|---:|---:|---:|---:|---:|---:|
| incumbent | 36.2327 | 1.0000 | 3.8099 | 4.0160 | 3.8757 | 12/12 |
| combination_01 | 23.7385 | 0.6503 | 2.4777 | 2.6117 | 2.5204 | 12/12 |
| combination_02 | 26.7312 | 0.7665 | 2.9205 | 3.0784 | 2.9708 | 12/12 |
| union_4 | 38.6873 | 0.8036 | 3.0618 | 3.2274 | 3.1146 | 12/12 |
| implementation_01 | 28.2180 | 0.7802 | 2.9727 | 3.1334 | 3.0239 | 12/12 |
| implementation_10 | 35.7428 | 0.9813 | 3.7386 | 3.9408 | 3.8031 | 12/12 |
| aa_control | 27.4043 | 0.9142 | 3.4829 | 3.6713 | 3.5430 | 12/12 |
| rho_w1 | 10.1807 | 0.2625 | 1.0000 | 1.0541 | 1.0172 | 12/12 |
| rho_w8 | 9.4395 | 0.2490 | 0.9487 | 1.0000 | 0.9651 | 12/12 |
| rho_w32 | 8.1268 | 0.2580 | 0.9830 | 1.0362 | 1.0000 | 12/12 |

### interaction_replication

| Variant | Cold ms | Paired / incumbent | / rho 1 | / rho 8 | / rho 32 | Verified |
|---|---:|---:|---:|---:|---:|---:|
| incumbent | 38.4235 | 1.0000 | 4.4269 | 3.8557 | 4.7951 | 12/12 |
| combination_01 | 27.1512 | 0.8077 | 3.5757 | 3.1143 | 3.8731 | 12/12 |
| combination_02 | 33.2033 | 0.8791 | 3.8916 | 3.3894 | 4.2152 | 12/12 |
| union_4 | 38.6317 | 0.9238 | 4.0894 | 3.5618 | 4.4296 | 12/12 |
| implementation_01 | 28.0912 | 0.7138 | 3.1600 | 2.7522 | 3.4228 | 12/12 |
| implementation_10 | 44.8494 | 1.1014 | 4.8759 | 4.2467 | 5.2814 | 12/12 |
| aa_control | 37.3434 | 0.8617 | 3.8149 | 3.3226 | 4.1322 | 12/12 |
| rho_w1 | 10.3586 | 0.2259 | 1.0000 | 0.8710 | 1.0832 | 12/12 |
| rho_w8 | 10.6090 | 0.2594 | 1.1482 | 1.0000 | 1.2436 | 12/12 |
| rho_w32 | 8.5442 | 0.2085 | 0.9232 | 0.8041 | 1.0000 | 12/12 |

## Validation and continuation

Release tests: 62 descent-related tests passed; 36 rho-related tests passed and 3 expensive tests remained ignored. The Python suites have 50 tournament/autolab tests (2 Linux-only capability skips) and 14 boundary tests. All four repo skills validate. See `validation/` for exact output.

Restore raw evidence with `python3 research/ic_candidate_tournament_20260915/evidence/restore.py --archive autolab-20260924`; audit each completed panel using its frozen `evaluator/autolab.py verify --round PATH`. The restored artifact checksum does not require executing its archived binary.

Next: restore the archived optimized IC winner; qualify a strong matched rho reference; run a new Linux-amd64/Valgrind instruction tournament with fresh, disjoint confirmation points, the retained portfolio and the combination candidates. For new F4/F5/SAT implementation changes, also run the applicable frozen algebra regression and fresh algebra holdouts. Keep the one-column policies in the archive for a separately validated admission protocol.
