# Stage 160: constant-linear fixed-X1 S4 construction

Strong internal engineering and one finite public toy-PDP native-F4 construction improvement. It is not licensed Magma F4, a full native-F4 panel, an end-to-end index-calculus speedup, a direct-MITM improvement, or a Koblitz index-calculus SOTA.

The selected constructor forms the shared symbolic `X2*X3` product once and applies the fixed `X1` and target powers as exact linear maps in the polynomial basis. A direct unit gate compares its complete Boolean polynomial vectors with the generic S4 expansion over three fields, five fixed values, and four targets per field.

| Backend | Host | Terminal | Wall s | Core s | Peak RSS B | Conflicts / native work |
|:--|:--|:--|--:|--:|--:|:--|
| native-f4-selected | local arm64 | sat | 73.757990 | 73.585054 | 910934016 | {'word_xors_elimination_only': 87513949370, 'f4_calls': 65} |
| native-xor | local arm64 | unknown_inconclusive | 8.906472 | 8.890455 | 168361984 | 100000 |
| direct-mitm | local arm64 | sat | 2.999661 | 2.991745 | 45842432 | {'group_additions': 117015, 'pair_entries': 116635} |
| wdsat | GitHub Linux x86_64 | timeout_inconclusive | 120.003056 | 88.060134 | 19632128 | null |
| cryptominisat | GitHub Linux x86_64 | timeout_inconclusive | 120.005939 | 87.983526 | 69353472 | null |
| magma-f4 | not run | not_run_licensed_tool_missing | null | null | null | null |
| ggmp | not run | not_same_instance | null | null | null | null |

The clean selected F4 run takes 73.757990 wall seconds and 73.585054 core-seconds. Construction falls from 26.323585 s to 0.973558 s, a 27.04x improvement. Whole-target wall improves 1.36x with identical equations, F4 calls, word XORs, roots, and verified witness.

Direct MITM still wins decisively on the same host: selected F4 uses 24.59x its wall time, 24.60x its CPU, and 19.87x its RSS. The clean build-plus-run cost is 241.262 wall seconds and 236.554 core-seconds.

WDSat and CryptoMiniSat remain exact-target frozen Linux timeouts on a different host. Licensed Magma remains unexecuted, GGMP has no same n=59 cell, full IC/rho cost is unchanged, and independent reproduction is still outstanding. This does not establish a SOTA.
