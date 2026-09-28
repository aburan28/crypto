# Instruction accounting, frozen before instrumented runs

Use the identical generated `rho_parity_cost` harness with each revision.
The generator reads the immutable original Weil benchmark, independent of
candidate edits. Valgrind 3.22.0 Callgrind on x86-64 counts executed user-mode
instructions (`Ir`), including library allocation, field arithmetic, failed
queries, validation of candidate witnesses and final scalar verification.
All threads are counted. Four exclusive scopes cover curve setup, target
generation, factor-base/plan/index setup, and driver/final verification.
Rho uses the same first two scopes and its complete solver/verification scope.
The sum is the cold single-target total. Diagnostic JSON, independent exhaustive
truth tables, and external oracle replay remain outside these attack scopes.
No timing under Valgrind is used as a native runtime measurement.

Report raw instruction counts and paired candidate/rho instruction ratios.
For a documented conversion, measure 4,096 additions of public, seeded random
nonzero subgroup point pairs on each curve and a matching loop control. The
mean incremental instructions per addition is a machine-specific conversion,
not an assertion that every group operation has the same cost. Report
`S_I = total_I / (measured_I_per_addition * sqrt(N))` separately from prior
algebraic-operation scoreboards. These counts include all scoped phases but
exclude kernel execution, physical memory latency and network/cloud costs
(caches are off). They establish implementation cost under this model, not
a hardware-independent cryptanalytic advance. The counting-yield bound does
not supply a calibrated instruction floor; its cost ratio remains null.

Run all ten cases and all eight frozen/fresh seed/log pairs for baseline and
each retained candidate, with pair-table and rho controls. Keep failed DLPs
and every profile. Instrumented verified workloads must agree with native
ones. Numerical parity requires candidate/rho <=1 on every useful case,
with no missing answers; also retain the stricter native runtime gate in
README.md. Neither gate implies degree-131 or asymptotic parity.

Profiler semantics: [Callgrind manual](https://valgrind.org/docs/manual/cl-manual.html).
Raw profiles use exclusive counts; summing inclusive function costs would
double count callees and is explicitly prohibited.
