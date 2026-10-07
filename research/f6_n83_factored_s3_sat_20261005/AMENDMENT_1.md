# Amendment 1: native XOR rows for the factored S3 circuit

Registered after the original planted Kissat timeout was frozen in
commit `49bccdd80` and before changing the emitter or rerunning a solver.
The 30-second Kissat control fixed all 90 source bits but did not infer
the 249 intermediate bits. It yielded no model, so the original runner
correctly stopped before natural T001 offsets. This amendment changes
only the **representation and solver**: retain the exact same S3 field
circuit, but write each XOR gate as a native extended-DIMACS parity row
instead of four CNF clauses. AND gates keep their three exact Tseitin
clauses. Use CryptoMiniSat 5.14.7 at
`/opt/homebrew/opt/cryptominisat/bin/cryptominisat5`, binary SHA-256
`a3f85c3709b5e2a040bf82a4a604d1c7b9f10219bbf180a9e0f72319a2e892ac`.
Freeze its default Gaussian settings and one thread. Record all exact
extended-DIMACS inputs, gate/row/variable counts and hashes.

First run a correctness control with **all 339** source and
intermediate input bits fixed to the planted `[0,2,4,6,8]` assignment.
This tests the emitted parity-row convention, solver model parser and
full-group witness verifier; it does not estimate search cost. Require a
verified planted model within 30 seconds. Then run a second planted
control with only the 90 source bits fixed, 30 seconds, to measure
whether native XOR reasoning closes the intermediate bits. Preserve a
timeout and continue to the ordinary gate if the fully pinned control
passed. For ordinary public T001, run all four torsion offsets in order
0,1,2,3, 120 seconds each. Use the exact same factor base, target lifts,
expanded-polynomial cross-checks and full-group verification as the
original protocol. One ordinary verified witness is structural success;
SAT without a full-group witness is algebraic-only evidence. Timeouts
are inconclusive, UNSAT is limited to that offset's exact encoding.

Keep one solver process at a time, `--threads=1`, a live 7-GiB RSS
sampler/kill, exact stdout/stderr and exit statuses. Preserve the
initial Kissat failure untouched. CPU times on this contended host are
feasibility diagnostics, never a controlled F6/F4/F5 or IC/rho speed
claim. The complete candidate ID remains unset until a complete
pipeline exists.
