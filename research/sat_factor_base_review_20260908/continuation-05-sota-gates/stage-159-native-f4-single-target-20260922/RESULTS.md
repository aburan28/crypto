# Stage 159: native F4 on one frozen n=59 Phase-B target

Strong internal engineering and one finite public toy-PDP native-F4 result. It is not licensed Magma F4, a full 160-input native-F4 panel, an end-to-end index-calculus speedup, or a Koblitz index-calculus SOTA.

The selected fixed-X1 F4 arm is class-blind during execution. It verifies the historical source export, defines the base as the public polynomial subspace, skips only x-values with no rational curve lift, and accepts SAT only after an exact three-point group check. Truth scoring occurs after the run seal.

| Backend | Host | Terminal | Wall s | Core s | Peak RSS B | Conflicts / native work |
|:--|:--|:--|--:|--:|--:|:--|
| native-f4-selected | local arm64 | sat | 100.179380 | 99.920352 | 921681920 | {'word_xors_elimination_only': 87513949370, 'f4_calls': 65} |
| native-xor | local arm64 | unknown_inconclusive | 8.906472 | 8.890455 | 168361984 | 100000 |
| direct-mitm | local arm64 | sat | 2.999661 | 2.991745 | 45842432 | {'group_additions': 117015, 'pair_entries': 116635} |
| wdsat | GitHub Linux x86_64 | timeout_inconclusive | 120.003056 | 88.060134 | 19632128 | null |
| cryptominisat | GitHub Linux x86_64 | timeout_inconclusive | 120.005939 | 87.983526 | 69353472 | null |
| magma-f4 | not run | not_run_licensed_tool_missing | null | null | null | null |
| ggmp | not run | not_same_instance | null | null | null | null |

Rational-X1 membership reduced F4 wall to 0.487 of the complete fixed-X1 baseline, core to 0.487, and elimination word XORs to 0.479; the exact witness is unchanged.

On the same local host and target, selected native F4 is 33.40 times direct MITM wall, 33.40 times its CPU, and 20.11 times its RSS. The trace constraint is rejected for speed: 1.182 wall and 1.588 word XORs, despite lower RSS.

The WDSat and CryptoMiniSat rows are exact-target frozen Linux receipts, but their host differs from the local F4 host; their seconds are descriptive and are not used as paired speed ratios. Licensed Magma remains unexecuted. GGMP has no same n=59 cell. Conflict count is null for F4 and MITM because neither exposes SAT conflicts.

The measured campaign total is a lower bound because some early developer checks lacked an outer meter. Full end-to-end IC cost and the ratio to automorphism-optimized rho remain unchanged and false. This does not establish a SOTA.
