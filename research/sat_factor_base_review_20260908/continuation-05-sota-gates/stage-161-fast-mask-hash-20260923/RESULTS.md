# Stage 161: fast hashing for trusted F4 monomial masks

Strong internal engineering and one finite public toy-PDP native-F4 hashing improvement. It is not licensed Magma F4, a full native-F4 panel, an end-to-end index-calculus speedup, a direct-MITM improvement, or a Koblitz index-calculus SOTA.

A charged 20-second sample of the Stage 160 binary put 1,143 of 3,334 top-of-stack samples in dense echelon elimination, 533 in critical-pair installation, 238 in column packing, and 225 directly in SipHash writes. Stage 161 retains exact-key hash maps but replaces SipHash only for solver-constructed `u64` monomial masks with deterministic SplitMix64 hashing.

| Backend | Host | Terminal | Wall s | Core s | Peak RSS B | Conflicts / native work |
|:--|:--|:--|--:|--:|--:|:--|
| native-f4-selected | local arm64 | sat | 60.645599 | 60.498754 | 964263936 | {'word_xors_elimination_only': 87513949370, 'f4_calls': 65} |
| native-xor | local arm64 | unknown_inconclusive | 8.906472 | 8.890455 | 168361984 | 100000 |
| direct-mitm | local arm64 | sat | 2.999661 | 2.991745 | 45842432 | {'group_additions': 117015, 'pair_entries': 116635} |
| wdsat | GitHub Linux x86_64 | timeout_inconclusive | 120.003056 | 88.060134 | 19632128 | null |
| cryptominisat | GitHub Linux x86_64 | timeout_inconclusive | 120.005939 | 87.983526 | 69353472 | null |
| magma-f4 | not run | not_run_licensed_tool_missing | null | null | null | null |
| ggmp | not run | not_same_instance | null | null | null | null |

The clean selected F4 process takes 60.645599 wall seconds and 60.498754 core-seconds. Symbolic/matrix build falls from 30.523303 to 17.204156 seconds, a 1.77x improvement. Whole-target wall improves 1.22x with identical equations, calls, matrices, pair counters, word XORs, roots, and verified witness.

Peak RSS increases 5.85% to 964263936 bytes. Direct MITM still wins: selected F4 uses 20.22x its wall, 20.22x its CPU, and 21.03x its RSS. Clean build plus run costs 228.348 wall seconds and 224.536 core-seconds.

Licensed Magma, a full native-F4 panel, end-to-end IC/rho cost, and independent external reproduction remain open. This one-target engineering result does not establish a SOTA.
