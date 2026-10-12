# Pre-run backend clarification

The available `pycryptosat==5.16.0` binding exposes one-thread construction,
native XOR and per-call time limits, but no random-seed argument. Its backend
uses its version's default seed. The extension's intended seed 13 is applied
to Z3; CryptoMiniSat reports `seed_configurable=false`, exact version and
single-thread configuration. This is fixed before any CryptoMiniSat run and
precludes a seed-controlled claim across the SAT backends. All instances,
resource limits, correctness checks and comparison groups remain frozen.

An executable WDSat/Sage/msolve backend is unavailable in this environment.
Do not call the Boolean Macaulay implementation `msolve` or F5B Magma F5;
they are explicitly the pinned Python prototypes. Mark missing engines as
`unavailable`, with no fabricated runtime or speedup.
