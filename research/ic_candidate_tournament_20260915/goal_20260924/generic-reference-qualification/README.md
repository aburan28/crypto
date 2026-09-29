# Generic/reference comparison and observer evidence

[Frozen protocol](PROTOCOL.md) · [Full result export](RESULTS.json) · [Exposed fixtures](fixtures.json) · [Completed workflow](https://github.com/aburan28/crypto/actions/runs/36290704597)

## Measured development comparison; no promoted winner

The registered run 36290704597 completed all 1,350 native/profile pairs (30 A/A, 330 smoke and 990 development), with zero research failures. The whole-mode observer study completed all 360 pairs, with zero incomplete or failed pairs. This is accounting and reference qualification on public synthetic toy groups. It consumes no improvement round and does not change the accepted reference binding. One of three improvement rounds is complete; two remain. No global optimum or held-out improvement is established.

## Decision and comparison scope

The optimized pairinv incumbent remains the cold IC leader. Prepared both is the online point-estimate leader, at 0.9770794 times incumbent time, with a descriptive 95% interval [0.9176026, 1.043898]; this does not establish an IC improvement. Generic dense and sparse IC are roughly 3.7 times slower online, ten times slower in cold native time and 37 times more costly in cold instructions than the optimized incumbent. Future rounds must retain the optimized reference. Choosing the generic worker as a weaker baseline would not establish progress. The frozen selectors choose rho_incumbent_8 for cold instructions and rho_generic_dense_16 for online time. These are development selections; a reviewed, versioned binding is still required.

## Primary metric and uncertainty

Every job solves one supplied public point. Native online time begins after reusable factor-base, index and log preparation, includes all target-dependent attempts, and ends after scalar replay. Fixture construction, launch and input loading remain outside this interval. Three process repetitions are medianed per point; equal-cell geometric means combine 15 development points across five cells. The repetitions are not 45 independent targets and this is not shared-target amortization. Intervals are descriptive paired development bootstraps, with the frozen estimator and seed; they do not correct for selecting leaders on this panel. IC/rho intervals resample their paired observations directly, rather than dividing marginal confidence endpoints.

## Single-target online wall time: all 22 configurations

| Variant | Online ms | / incumbent | Descriptive 95% interval | Selected online rho / variant | Paired rho/IC 95% interval | Verified / scheduled |
| --- | --- | --- | --- | --- | --- | --- |
| incumbent | 0.03564707 | 1 | [1, 1] | 7.407728 | [6.127355, 8.983053] | 45/45 |
| prepared_both | 0.03483002 | 0.9770794 | [0.9176026, 1.043898] | 7.5815 | [6.447833, 9.08608] | 45/45 |
| generic_dense | 0.1317613 | 3.696273 | [3.147132, 4.194848] | 2.004107 | [1.751399, 2.36352] | 45/45 |
| generic_sparse | 0.1331706 | 3.735807 | [3.158101, 4.305252] | 1.982899 | [1.758373, 2.330959] | 45/45 |
| rho_incumbent_1 | 0.3244606 | 9.102027 | [7.171945, 12.1651] | 0.8138547 | unknown / not applicable | 45/45 |
| rho_incumbent_2 | 0.303651 | 8.518258 | [7.064248, 10.17534] | 0.8696294 | unknown / not applicable | 45/45 |
| rho_incumbent_4 | 0.2976156 | 8.348948 | [6.967493, 9.902724] | 0.8872648 | unknown / not applicable | 45/45 |
| rho_incumbent_8 | 0.2920471 | 8.192736 | [6.775364, 9.869981] | 0.9041824 | unknown / not applicable | 45/45 |
| rho_incumbent_16 | 0.2949671 | 8.274651 | [6.818359, 9.967512] | 0.8952314 | unknown / not applicable | 45/45 |
| rho_incumbent_32 | 0.2978373 | 8.355169 | [6.918814, 9.947265] | 0.8866041 | unknown / not applicable | 45/45 |
| rho_prepared_both_1 | 0.3236894 | 9.080393 | [7.170246, 12.09378] | 0.8157937 | unknown / not applicable | 45/45 |
| rho_prepared_both_2 | 0.3045554 | 8.543631 | [7.103577, 10.19957] | 0.8670468 | unknown / not applicable | 45/45 |
| rho_prepared_both_4 | 0.2971817 | 8.336778 | [6.961106, 9.859495] | 0.88856 | unknown / not applicable | 45/45 |
| rho_prepared_both_8 | 0.2986676 | 8.378461 | [7.03723, 9.9969] | 0.8841394 | unknown / not applicable | 45/45 |
| rho_prepared_both_16 | 0.2983169 | 8.368624 | [6.85242, 10.14203] | 0.8851787 | unknown / not applicable | 45/45 |
| rho_prepared_both_32 | 0.2896356 | 8.125087 | [6.726716, 9.812233] | 0.9117106 | unknown / not applicable | 45/45 |
| rho_generic_dense_1 | 0.2975753 | 8.34782 | [6.755088, 10.82661] | 0.8873847 | unknown / not applicable | 45/45 |
| rho_generic_dense_2 | 0.2775978 | 7.787394 | [6.550785, 9.379188] | 0.951246 | unknown / not applicable | 45/45 |
| rho_generic_dense_4 | 0.2815038 | 7.896969 | [6.557484, 9.483427] | 0.938047 | unknown / not applicable | 45/45 |
| rho_generic_dense_8 | 0.2711357 | 7.606114 | [6.403084, 8.861318] | 0.9739175 | unknown / not applicable | 45/45 |
| rho_generic_dense_16 | 0.2640638 | 7.407728 | [6.127355, 8.983053] | 1 | unknown / not applicable | 45/45 |
| rho_generic_dense_32 | 0.2646889 | 7.425264 | [6.469343, 8.496569] | 0.9976383 | unknown / not applicable | 45/45 |

## Paired supplied points

The next table shows the selected online IC and rho on each of the 15 development points. Every interval is one supplied point after reusable preparation through scalar replay. Full IC1 candidate IDs, RHO1 reference IDs, workload IDs, execution IDs and all four IC comparisons are retained in RESULTS.json; aliases below map to those exact manifests. Speedup is rho_online_ms / IC_online_ms, limited to this charged instrumented-worker comparison.

## Selected online IC and rho: point-level observations

| Case | Public Q (x,y) | Workload | IC alias | IC ms | Rho alias | Rho ms | Rho / IC | Verified |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| n17a1-000 | 13339,66814 | 9a33222ec2d9 | prepared_both | 0.030817 | rho_generic_dense_16 | 0.184645 | 5.99166 | True |
| n17a1-001 | 56873,110715 | 07c397c45d01 | prepared_both | 0.028533 | rho_generic_dense_16 | 0.202739 | 7.105422 | True |
| n17a1-002 | 51397,57521 | c32b90bc5602 | prepared_both | 0.030837 | rho_generic_dense_16 | 0.182851 | 5.929598 | True |
| n19a0-000 | 296334,340408 | f31b422d6fae | prepared_both | 0.026941 | rho_generic_dense_16 | 0.181931 | 6.752942 | True |
| n19a0-001 | 448223,478555 | b3b34806a9a1 | prepared_both | 0.027351 | rho_generic_dense_16 | 0.174265 | 6.371431 | True |
| n19a0-002 | 123350,229349 | 11877dbd59fd | prepared_both | 0.029425 | rho_generic_dense_16 | 0.235541 | 8.004792 | True |
| n23a0-000 | 3078340,7150380 | 382eb968aa38 | prepared_both | 0.044463 | rho_generic_dense_16 | 0.257622 | 5.794076 | True |
| n23a0-001 | 3481084,224419 | 3cfca02b926f | prepared_both | 0.038532 | rho_generic_dense_16 | 0.398495 | 10.34192 | True |
| n23a0-002 | 2882272,5658981 | e7824ff19518 | prepared_both | 0.040696 | rho_generic_dense_16 | 0.275496 | 6.769609 | True |
| n23a1-000 | 6292485,1304930 | da9c44af019d | prepared_both | 0.037029 | rho_generic_dense_16 | 0.305391 | 8.247347 | True |
| n23a1-001 | 4304213,1360375 | df590c7f4d1e | prepared_both | 0.034955 | rho_generic_dense_16 | 0.481079 | 13.76281 | True |
| n23a1-002 | 4316435,4581027 | c1af3f76b644 | prepared_both | 0.058129 | rho_generic_dense_16 | 0.265036 | 4.559445 | True |
| n31a0-000 | 2137032792,1243931738 | bedcd6bd70d3 | prepared_both | 0.035236 | rho_generic_dense_16 | 0.284933 | 8.086417 | True |
| n31a0-001 | 645710776,1570337940 | 7803cf52aac1 | prepared_both | 0.034905 | rho_generic_dense_16 | 0.394889 | 11.31325 | True |
| n31a0-002 | 1069130512,1602516399 | aabc9ccee14a | prepared_both | 0.035506 | rho_generic_dense_16 | 0.334165 | 9.411508 | True |

## Supplementary complete cold accounting

Cold native time is the full cold child process, including launch, input and reporting tail. It is not the primary single-target online metric. Cold Ir is Valgrind 3.22 amd64 user-space guest instructions with exclusive complete phase closure; it is not curve additions. S = Ir / sqrt(r) uses subgroup order r, not field size. The existing K-instruction floor applies only to these full-rank relation collectors and is not a universal IC bound; it is inapplicable to rho. No observer or witness overhead is subtracted.

## Supplementary cold native wall time

| Variant | Cold process ms | / incumbent | Descriptive 95% interval | Verified / scheduled |
| --- | --- | --- | --- | --- |
| incumbent | 0.9671105 | 1 | [1, 1] | 45/45 |
| prepared_both | 0.9671934 | 1.000086 | [0.9437302, 1.053397] | 45/45 |
| generic_dense | 9.718014 | 10.0485 | [8.437789, 12.73373] | 45/45 |
| generic_sparse | 9.698364 | 10.02819 | [8.467366, 12.7153] | 45/45 |
| rho_incumbent_1 | 0.9732892 | 1.006389 | [0.9330855, 1.106445] | 45/45 |
| rho_incumbent_2 | 0.9460796 | 0.9782538 | [0.9236538, 1.032655] | 45/45 |
| rho_incumbent_4 | 0.9318276 | 0.9635172 | [0.908026, 1.015199] | 45/45 |
| rho_incumbent_8 | 0.9318193 | 0.9635086 | [0.9021839, 1.026448] | 45/45 |
| rho_incumbent_16 | 0.9296367 | 0.9612517 | [0.9014169, 1.02424] | 45/45 |
| rho_incumbent_32 | 0.9297516 | 0.9613706 | [0.9002262, 1.024885] | 45/45 |
| rho_prepared_both_1 | 0.9873451 | 1.020923 | [0.9383698, 1.126567] | 45/45 |
| rho_prepared_both_2 | 0.9469785 | 0.9791833 | [0.9273201, 1.030762] | 45/45 |
| rho_prepared_both_4 | 0.9326889 | 0.9644077 | [0.9120017, 1.013289] | 45/45 |
| rho_prepared_both_8 | 0.9372449 | 0.9691187 | [0.8976228, 1.040794] | 45/45 |
| rho_prepared_both_16 | 0.9407592 | 0.9727525 | [0.9064158, 1.036223] | 45/45 |
| rho_prepared_both_32 | 0.9286504 | 0.9602319 | [0.8992078, 1.023721] | 45/45 |
| rho_generic_dense_1 | 2.03265 | 2.101776 | [1.949268, 2.268944] | 45/45 |
| rho_generic_dense_2 | 2.001914 | 2.069995 | [1.91672, 2.240392] | 45/45 |
| rho_generic_dense_4 | 2.02344 | 2.092253 | [1.932762, 2.285182] | 45/45 |
| rho_generic_dense_8 | 1.997169 | 2.065088 | [1.89729, 2.253677] | 45/45 |
| rho_generic_dense_16 | 2.000877 | 2.068923 | [1.913883, 2.249989] | 45/45 |
| rho_generic_dense_32 | 1.996225 | 2.064112 | [1.908056, 2.243763] | 45/45 |

## Supplementary complete cold instructions and boundaries

| Variant | Cold Ir | S = Ir/sqrt(r) | / incumbent | Descriptive 95% interval | / selected cold rho | / K-instruction floor | Verified / scheduled |
| --- | --- | --- | --- | --- | --- | --- | --- |
| incumbent | 1826648 | 2280.558 | 1 | [1, 1] | 0.7281677 | 228330.9 | 45/45 |
| prepared_both | 1859774 | 2321.915 | 1.018135 | [1.009134, 1.026871] | 0.7413729 | 232471.7 | 45/45 |
| generic_dense | 6.810338e+07 | 85026.63 | 37.28326 | [29.65719, 46.23664] | 27.14847 | 8512922 | 45/45 |
| generic_sparse | 6.807298e+07 | 84988.67 | 37.26662 | [29.65709, 46.23805] | 27.13635 | 8509122 | 45/45 |
| rho_incumbent_1 | 2822509 | 3523.885 | 1.545185 | [1.333969, 1.85461] | 1.125154 | unknown / not applicable | 45/45 |
| rho_incumbent_2 | 2615622 | 3265.587 | 1.431925 | [1.322881, 1.537357] | 1.042682 | unknown / not applicable | 45/45 |
| rho_incumbent_4 | 2571092 | 3209.991 | 1.407547 | [1.239192, 1.547416] | 1.02493 | unknown / not applicable | 45/45 |
| rho_incumbent_8 | 2508553 | 3131.913 | 1.37331 | [1.205196, 1.528233] | 1 | unknown / not applicable | 45/45 |
| rho_incumbent_16 | 2508557 | 3131.918 | 1.373312 | [1.205226, 1.528208] | 1.000002 | unknown / not applicable | 45/45 |
| rho_incumbent_32 | 2508572 | 3131.936 | 1.37332 | [1.205224, 1.528203] | 1.000007 | unknown / not applicable | 45/45 |
| rho_prepared_both_1 | 2822736 | 3524.168 | 1.54531 | [1.334079, 1.854718] | 1.125245 | unknown / not applicable | 45/45 |
| rho_prepared_both_2 | 2615828 | 3265.845 | 1.432038 | [1.322983, 1.537509] | 1.042764 | unknown / not applicable | 45/45 |
| rho_prepared_both_4 | 2571307 | 3210.26 | 1.407664 | [1.239309, 1.547556] | 1.025016 | unknown / not applicable | 45/45 |
| rho_prepared_both_8 | 2508725 | 3132.127 | 1.373404 | [1.205294, 1.52833] | 1.000068 | unknown / not applicable | 45/45 |
| rho_prepared_both_16 | 2508754 | 3132.164 | 1.37342 | [1.20534, 1.528258] | 1.00008 | unknown / not applicable | 45/45 |
| rho_prepared_both_32 | 2508806 | 3132.228 | 1.373448 | [1.205366, 1.528315] | 1.000101 | unknown / not applicable | 45/45 |
| rho_generic_dense_1 | 4565777 | 5700.344 | 2.499539 | [2.134634, 2.886798] | 1.820084 | unknown / not applicable | 45/45 |
| rho_generic_dense_2 | 4392212 | 5483.649 | 2.404521 | [2.027831, 2.848576] | 1.750894 | unknown / not applicable | 45/45 |
| rho_generic_dense_4 | 4366872 | 5452.011 | 2.390648 | [2.029696, 2.839163] | 1.740793 | unknown / not applicable | 45/45 |
| rho_generic_dense_8 | 4317113 | 5389.887 | 2.363408 | [1.969635, 2.833795] | 1.720957 | unknown / not applicable | 45/45 |
| rho_generic_dense_16 | 4317109 | 5389.882 | 2.363405 | [1.969619, 2.833833] | 1.720955 | unknown / not applicable | 45/45 |
| rho_generic_dense_32 | 4317119 | 5389.895 | 2.363411 | [1.969642, 2.833814] | 1.72096 | unknown / not applicable | 45/45 |

## Actual factor bases and clipped rho widths

All four IC arms have the following independently checked usable base counts B before sign/Frobenius folding and effective column counts K. The requested points = 6n is a construction parameter, never B. Factor-base policy was explicitly a comparison variable. The raw receipts retain rank, verified yield, unsuccessful PDP attempts, query histories, matrix work, descent and certificates. Larger unresolved PDP queries remain unresolved; they do not become unsatisfiability proofs.

## Actual base inventory, shared counts across four IC arms

| Cell | Usable points B | Folded columns K |
| --- | --- | --- |
| n17a1 | 272 | 8 |
| n19a0 | 304 | 8 |
| n23a0 | 368 | 8 |
| n23a1 | 368 | 8 |
| n31a0 | 496 | 8 |

## Effective rho widths, in n17a1 / n19a0 / n23a0 / n23a1 / n31a0 order

| Variant | Observed widths by cell |
| --- | --- |
| rho_incumbent_1 | n17a1: [1]; n19a0: [1]; n23a0: [1]; n23a1: [1]; n31a0: [1] |
| rho_incumbent_2 | n17a1: [1]; n19a0: [1]; n23a0: [2]; n23a1: [2]; n31a0: [2] |
| rho_incumbent_4 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [4]; n31a0: [2] |
| rho_incumbent_8 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [5]; n31a0: [2] |
| rho_incumbent_16 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [5]; n31a0: [2] |
| rho_incumbent_32 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [5]; n31a0: [2] |
| rho_prepared_both_1 | n17a1: [1]; n19a0: [1]; n23a0: [1]; n23a1: [1]; n31a0: [1] |
| rho_prepared_both_2 | n17a1: [1]; n19a0: [1]; n23a0: [2]; n23a1: [2]; n31a0: [2] |
| rho_prepared_both_4 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [4]; n31a0: [2] |
| rho_prepared_both_8 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [5]; n31a0: [2] |
| rho_prepared_both_16 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [5]; n31a0: [2] |
| rho_prepared_both_32 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [5]; n31a0: [2] |
| rho_generic_dense_1 | n17a1: [1]; n19a0: [1]; n23a0: [1]; n23a1: [1]; n31a0: [1] |
| rho_generic_dense_2 | n17a1: [1]; n19a0: [1]; n23a0: [2]; n23a1: [2]; n31a0: [2] |
| rho_generic_dense_4 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [4]; n31a0: [2] |
| rho_generic_dense_8 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [5]; n31a0: [2] |
| rho_generic_dense_16 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [5]; n31a0: [2] |
| rho_generic_dense_32 | n17a1: [1]; n19a0: [1]; n23a0: [4]; n23a1: [5]; n31a0: [2] |

## Interpretation of width comparisons

Requested widths 8, 16 and 32 clip to the same widths 1, 1, 4, 5 and 2 on this panel. They are retained as scheduled configurations, but do not constitute independent algorithms or additional target samples. A selected width cannot support a general batching advantage.

## Cold instruction phase shares: stage diagnostics only

| Phase | incumbent (%) | prepared_both (%) | generic_dense (%) | generic_sparse (%) |
| --- | --- | --- | --- | --- |
| factor_base | 9.19431 | 9.03158 | 2.83991 | 2.84105 |
| isogeny | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| matrix_build | 0.34291 | 0.33800 | 2.12701 | 2.16075 |
| pdp | 21.41665 | 22.55614 | 3.73324 | 3.74131 |
| precompute | 22.28751 | 21.91274 | 76.97399 | 76.98531 |
| queries | 1.13315 | 1.11416 | 0.26451 | 0.26501 |
| recovery_check | 9.81296 | 9.67559 | 0.14503 | 0.14529 |
| relation_check | 5.65615 | 5.56964 | 2.76250 | 2.36337 |
| relation_la | 1.56821 | 1.53935 | 0.17364 | 0.13795 |
| setup | 24.93305 | 24.56466 | 10.11048 | 10.48977 |
| target_descent | 3.65510 | 3.69815 | 0.86969 | 0.87018 |

## What the phase ledger supports

Shares are means of per-run fractions of complete cold Ir, with equal target and cell weight. Generic precomputation accounts for about 77%; final relation LA accounts for less than 0.2%. The optimized incumbent spends about 1.6% in final relation LA. Eliminating that stage alone therefore cannot supply a 20% complete-cold improvement on these measurements. This is prioritization evidence, not attribution of cost to timers or witness generation. F4/F5 internal matrix reduction and final relation LA remain distinct stages.

## Whole-mode observer effects

The registered enabled/legacy study checks semantic agreement and retains both outputs. It changes more than timers: collection mode, tracing and matrix/batch reporting also differ. Both modes retain witness reporting. Legacy has OBS1 observation identifiers and null scientific admission and phase costs. The common outer interval omits enabled phase-boundary bookkeeping, so it cannot replace the fully charged scientific online interval. Ratios near one do not establish zero or low timer overhead. No overhead is subtracted. The following 95% intervals are descriptive, from 10,000 paired within-cell resamples; they do not supply promotion evidence.

## Observer: Common supplied-point outer wall interval

| Variant | Enabled / legacy | Descriptive 95% interval | Complete pairs |
| --- | --- | --- | --- |
| generic_dense | 1.002486 | [0.9722461, 1.036954] | 45 |
| generic_sparse | 0.9796783 | [0.9516532, 1.008251] | 45 |
| rho_generic_dense_1 | 1.072294 | [1.02128, 1.13056] | 45 |
| rho_generic_dense_16 | 0.9947425 | [0.9581826, 1.034325] | 45 |
| rho_generic_dense_2 | 1.000835 | [0.9519597, 1.053193] | 45 |
| rho_generic_dense_32 | 1.016135 | [0.9725968, 1.063648] | 45 |
| rho_generic_dense_4 | 1.033987 | [0.9952481, 1.073919] | 45 |
| rho_generic_dense_8 | 1.02511 | [0.9769169, 1.079501] | 45 |

## Observer: Whole-process wall time

| Variant | Enabled / legacy | Descriptive 95% interval | Complete pairs |
| --- | --- | --- | --- |
| generic_dense | 0.9991448 | [0.9878167, 1.011534] | 45 |
| generic_sparse | 1.018569 | [1.00936, 1.027498] | 45 |
| rho_generic_dense_1 | 1.13504 | [1.095328, 1.176311] | 45 |
| rho_generic_dense_16 | 1.060653 | [1.003871, 1.124465] | 45 |
| rho_generic_dense_2 | 1.066009 | [1.001955, 1.125561] | 45 |
| rho_generic_dense_32 | 1.083342 | [1.02697, 1.14215] | 45 |
| rho_generic_dense_4 | 1.030169 | [0.9734824, 1.089787] | 45 |
| rho_generic_dense_8 | 1.062823 | [1.008795, 1.115971] | 45 |

## Noise control and environment

The 30-run A/A instruction gate passed. A/A online ratio is 0.9558727, interval [0.8785708, 1.005275]; cold process ratio is 0.9752561, interval [0.9240123, 1.025244]. This shows material wall-time noise. The experiment used the same Linux x86-64 Azure CI host, Rust 1.94.1, Valgrind 3.22.0, CPU affinity 3, one Rayon thread, 8 GiB and a 180-second child cap. The generic worker reported pclmulqdq dispatch. This is a virtualized CI environment; no bare-metal, ARM, GPU or FPGA result is implied. Host/runtime/build records are retained. The host manifest does not record the physical CPU model or concurrent host load; neither is inferred.

## Retained failures, portability and delivery gate

Research failures are zero, but the prerequisite controls retain intentional incomplete runs, including 12 incomplete mixed-adapter profiles and nine incomplete observer-control pairs. Nine prerequisite audits and the 360-pair observer replay passed locally without workers. Strict macOS replay of the main comparison failed on a one-ULP statistic: development/rho_comparisons/rho_prepared_both_8/per_cell/n23a0 was stored as 1.1843664724946852 and recomputed as 1.1843664724946854. Original Linux bytes and this failure are preserved. The dedicated Linux evidence CI must freshly restore the archive, reconstruct every export and pass all eleven frozen audits without worker execution. The verifier has not been weakened to accommodate the local difference.

## Next bounded research step

Do not redispatch this qualification, sealed round one, or sealed round two. The reviewed observer evidence and separate cold/online IC and rho leaders remain bound in the version-two reference contract. Round two completed under seed 2026092552 and retained the incumbent; selected stop5_word failed the frozen familywise and complete-cost gates. Exclude all exposed points from that round and earlier panels, including the readiness, qualification and adapter-control corpora. A third attempt needs a new registered panel and seed; never retune on confirmation or replay. Promotion still requires fresh confirmation and replay under the predeclared familywise rule, at least 20% lower complete cold Ir and cold native time, no online regression, no cell regression above 10%, and independently verified answers. One round remains.

## Durable archive and reproduction

The archive is committed in this repository, with its identity in [the evidence manifest](../../evidence/manifest.json). It contains the comparison, observer study, frozen evaluators, raw profiles, source and controlled build, prerequisite artifacts, GitHub artifact provenance, workflow record and local audit receipts. Compression preserves original profile bytes; no measurement is rewritten.

Archive: ic-generic-reference-qualification-20260926.tar.zst; 25,623,677 bytes, 95,496 retained files, 348,522,674 reconstructed file bytes.

SHA-256: `680369a9b822dcc3321380bff751bd2b94f0892c915f201a67916fc7bd49fc5a`.

From the repository root, using Python 3.12 and zstd:

```sh
python3.12 -m unittest discover \
  -s research/ic_candidate_tournament_20260915/goal_20260924/generic-reference-qualification \
  -p test_archive.py -v
```

Run on Linux for the complete strict replay; the Linux-only audit test is explicitly skipped on macOS. No test executes a research worker. The separate evidence workflow avoids repeating the large historical replay in every producer job.

The immutable fixture export has SHA-256 `b881798bfc56acdd8b8cc52a14b7501572c181a4e25e6bab1c8e7e61a742c7eb`. RESULTS.json retains all 990 canonical development run records, the original qualification and observer summaries, exact identities and 60 paired IC/rho point rows. The raw archive retains the complete 1,350-slot comparison and 360 observer pairs.

Regenerate this note and the canonical scoreboard with `render_report.py`; `render_report.py --check` verifies that both remain exact renderings of the exported evidence.
