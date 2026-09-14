# Incremental pullback benchmark results

Engineering candidate; calibrated normalized-cost classification remains pending. No end-to-end speedup is established. Counts below come from [comparison.json](comparison.json); the derivation, frozen gate and limitations are in [README.md](README.md).

The run contains 240 solver cells and 48 verified cold ECDLP runs. The audit checked 384 decomposition calls across 102 distinct targets, including empty answers. It reproduced 10 completed prior stage-counter records and all 6 prior cold baseline-counter records exactly. 1 formerly completed stage replay(s) were censored by this run's deadline; they are retained and excluded from equal-output ratios.

| Variant / reference | Full-DLP S | Cost/rho | Cost/floor | Class/status |
|---|---|---|---|---|
| coefficient-pullback | null | null | null | Reference/control |
| normalized-pullback | null | null | null | Engineering candidate; pending calibration |
| cached-pullback | null | null | null | Engineering candidate; pending calibration |
| s4-symmetric | null | null | null | Reference/control |
| chained-s3 | null | null | null | Reference/control |
| rho | null | null | null | Reference/control |

## Complete eleven-bit enumeration gate

Each corpus has four targets. These pairs complete with identical outputs. The coordinate cache preserves the normalized field vector exactly. The gate concerns multiplications; the API sum is an uncalibrated diagnostic.

| Corpus | Variant | Additions | Multiplications | Squarings | Field-API sum |
|---|---|---:|---:|---:|---:|
| frozen | coefficient-pullback | 97,804 | 90,822 | 11,545 | 200,171 |
| frozen | normalized-pullback | 70,028 | 47,770 | 11,165 | 128,963 |
| frozen | cached-pullback | 70,028 | 47,770 | 11,165 | 128,963 |
| fresh | coefficient-pullback | 97,024 | 90,446 | 11,509 | 198,979 |
| fresh | normalized-pullback | 69,486 | 47,394 | 11,129 | 128,009 |
| fresh | cached-pullback | 69,486 | 47,394 | 11,129 | 128,009 |

Frozen: 47.40% fewer multiplications; 4/4 complete equal-output targets. Frozen 10% diagnostic gate: PASS.

Fresh: 47.60% fewer multiplications; 4/4 complete equal-output targets. Frozen 10% diagnostic gate: PASS.

## Resolution at the fixed budget

Each row includes four targets in first mode plus the same four in enumeration mode. Resolved includes proved-empty cases. Relations found during incomplete enumeration do not make that cell resolved. Supported targets and uniform targets are distinct strata in the raw file; their mixture is not a natural relation-yield estimate.

| Corpus | n / d | Variant | Resolved / 8 | Resolved within 3 s / 8 |
|---|---|---|---:|---:|
| fresh | 11 / 5 | cached-pullback | 8 | 8 |
| fresh | 11 / 5 | chained-s3 | 8 | 8 |
| fresh | 11 / 5 | coefficient-pullback | 8 | 8 |
| fresh | 11 / 5 | normalized-pullback | 8 | 8 |
| fresh | 11 / 5 | s4-symmetric | 8 | 8 |
| fresh | 23 / 7 | cached-pullback | 8 | 8 |
| fresh | 23 / 7 | chained-s3 | 0 | 0 |
| fresh | 23 / 7 | coefficient-pullback | 1 | 1 |
| fresh | 23 / 7 | normalized-pullback | 8 | 8 |
| fresh | 23 / 7 | s4-symmetric | 0 | 0 |
| fresh | 29 / 8 | cached-pullback | 1 | 1 |
| fresh | 29 / 8 | chained-s3 | 0 | 0 |
| fresh | 29 / 8 | coefficient-pullback | 1 | 1 |
| fresh | 29 / 8 | normalized-pullback | 1 | 1 |
| fresh | 29 / 8 | s4-symmetric | 0 | 0 |
| frozen | 11 / 5 | cached-pullback | 8 | 8 |
| frozen | 11 / 5 | chained-s3 | 8 | 8 |
| frozen | 11 / 5 | coefficient-pullback | 8 | 8 |
| frozen | 11 / 5 | normalized-pullback | 8 | 8 |
| frozen | 11 / 5 | s4-symmetric | 8 | 8 |
| frozen | 23 / 7 | cached-pullback | 8 | 8 |
| frozen | 23 / 7 | chained-s3 | 0 | 0 |
| frozen | 23 / 7 | coefficient-pullback | 2 | 2 |
| frozen | 23 / 7 | normalized-pullback | 8 | 8 |
| frozen | 23 / 7 | s4-symmetric | 0 | 0 |
| frozen | 29 / 8 | cached-pullback | 2 | 2 |
| frozen | 29 / 8 | chained-s3 | 0 | 0 |
| frozen | 29 / 8 | coefficient-pullback | 0 | 0 |
| frozen | 29 / 8 | normalized-pullback | 1 | 1 |
| frozen | 29 / 8 | s4-symmetric | 0 | 0 |

## Cold ECDLP field-operation vectors

Each row sums three cold runs on the indicated corpus and group. Every run rebuilds setup and verifies its scalar. Scalar modular counts and phase times are separate in the comparison file. These partial vectors are not calibrated total operations. Times are descriptive; one repetition cannot establish a paired runtime-confidence claim.

| Corpus | Field bits / subgroup | Variant | Additions | Multiplications | Squarings | Cold seconds |
|---|---|---|---:|---:|---:|---:|
| frozen | 5 / 11 | coefficient-pullback | 16,066 | 22,476 | 4,595 | 0.050780 |
| frozen | 5 / 11 | normalized-pullback | 12,026 | 16,714 | 4,332 | 0.039678 |
| frozen | 5 / 11 | cached-pullback | 12,026 | 16,714 | 4,332 | 0.039599 |
| frozen | 5 / 11 | rho | 2,037 | 3,996 | 717 | 0.008465 |
| frozen | 9 / 127 | coefficient-pullback | 115,830 | 153,397 | 39,495 | 0.500559 |
| frozen | 9 / 127 | normalized-pullback | 88,449 | 112,424 | 38,690 | 0.394685 |
| frozen | 9 / 127 | cached-pullback | 88,449 | 112,424 | 38,690 | 0.391632 |
| frozen | 9 / 127 | rho | 25,866 | 45,389 | 19,869 | 0.156617 |
| fresh | 5 / 11 | coefficient-pullback | 16,033 | 22,630 | 4,653 | 0.050829 |
| fresh | 5 / 11 | normalized-pullback | 12,131 | 17,120 | 4,395 | 0.039511 |
| fresh | 5 / 11 | cached-pullback | 12,131 | 17,120 | 4,395 | 0.039501 |
| fresh | 5 / 11 | rho | 1,978 | 3,901 | 717 | 0.007597 |
| fresh | 9 / 127 | coefficient-pullback | 104,339 | 141,353 | 37,779 | 0.472088 |
| fresh | 9 / 127 | normalized-pullback | 81,302 | 106,582 | 37,091 | 0.381069 |
| fresh | 9 / 127 | cached-pullback | 81,302 | 106,582 | 37,091 | 0.375997 |
| fresh | 9 / 127 | rho | 25,691 | 45,023 | 19,869 | 0.156243 |

24/24 candidate/reference pairs followed identical full attempt paths. All costs remain in the totals if a valid first relation changes the downstream path. Fresh seeds can repeat Q in tiny groups: five-bit fresh seeds provide one new distinct target beyond the frozen panel; nine-bit fresh seeds provide two. The raw runs and distinctness audit preserve these collisions.

Rho is the same-group, same-target control. No field/scalar/control conversion or calibrated full-cost floor has been measured, so S and boundary ratios remain null. The branch exponent and support-counting ceiling are unchanged.

## Evidence and next gate

The initial harness import failure is preserved in `results/`; it occurred before any solver measurement. The accepted campaign is `results_v2/`. Its raw records include targets, partial outputs, failures, phase costs, scalar certificates and source hashes. Neither the contract nor candidate arithmetic changed after the first measurement.

Next iterations must rerun this entire matched suite, preserve failures and add fresh holdouts. Promote a candidate only after calibrated full-pipeline costs improve; larger-group scaling and runtime confidence intervals remain separate outstanding gates.
