# Stage 26: one-CPU Phase-B matrix

Stage 26 reruns the exact 160 truth-free Stage 20 source systems as four parallel cells. Each cell contains 40 inputs and executes native XOR SAT, WDSat, and CryptoMiniSat sequentially under a Linux singleton CPU affinity mask.

The packet is bound to inventory SHA-256 `c937afb9b172d114768b0b96a4b5ccf66e91fbf278b77c58ed49f0fd74af37f7`. Its verifier confirms that it contains no truth labels or known witnesses. The algebraic factor-base manifests state that the target subgroup was not enumerated and discrete-log labels were not used.

The run records elapsed time, total CPU time, per-process time and memory, solver conflicts, and a sampled live process-tree RSS peak. Pinned tool acquisition and metered build receipts are retained separately. Truth scoring remains withheld until every cell is sealed and downloaded.

This is a public synthetic PDP benchmark. It is not an end-to-end index-calculus result, does not pass the full-cost gate, and does not establish a Koblitz index-calculus SOTA.
