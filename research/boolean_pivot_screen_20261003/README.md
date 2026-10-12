# Sibling Boolean pivot-trace feasibility screen

`PROTOCOL.md` and `protocol.json` freeze a structural question before any
candidate code or results. The inputs are generated public quadratic Boolean
systems. The study will compare exact degree-3 Macaulay matrices after the two
values of one branch variable, measuring the rank of their difference and the
survival of a deterministic pivot trace. The Boolean derivative is affine, but
that does not imply the matrix difference is low rank.

The native Rust worker and replay verifier are implemented in the follow-on
branch. The source-pinned [discovery bundle](discovery_01/manifest.json)
completed all 80 public branch pairs and rejected full pivot-trace replay:
0/64 eligible n>=12 pairs passed the frozen screen. The result and its exact
mathematical limit are in [CONCLUSION.md](CONCLUSION.md), with hashes in
[RUNS.json](RUNS.json). The unused holdout seeds remain sealed in the protocol.
No speedup was measured. This work does not accept curve or key inputs.

`run.sh discovery NEW_DIRECTORY` builds and tests the worker, copies exact
sources and binary into a new bundle, writes the full packed matrix corpus and
structural records, then seals and replays every member. It refuses to overwrite
an attempt. A later `holdout` run requires the complete qualified discovery
bundle and is rejected by the worker if its discovery screen failed. The
structural stop rule is a feasibility filter, not a timing or cryptanalytic
claim.
