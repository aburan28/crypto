# Sibling Boolean pivot-trace feasibility screen

`PROTOCOL.md` and `protocol.json` freeze a structural question before any
candidate code or results. The inputs are generated public quadratic Boolean
systems. The study will compare exact degree-3 Macaulay matrices after the two
values of one branch variable, measuring the rank of their difference and the
survival of a deterministic pivot trace. The Boolean derivative is affine, but
that does not imply the matrix difference is low rank.

The native Rust worker and replay verifier are implemented in the follow-on
branch. Until a source-pinned discovery bundle is admitted, the status remains
**implementation pending measurement**: no discovery or holdout result has
passed and no speedup has been measured. The unused holdout seeds remain sealed
in the protocol. This work does not accept curve or key inputs.

`run.sh discovery NEW_DIRECTORY` builds and tests the worker, copies exact
sources and binary into a new bundle, writes the full packed matrix corpus and
structural records, then seals and replays every member. It refuses to overwrite
an attempt. A later `holdout` run requires the complete qualified discovery
bundle and is rejected by the worker if its discovery screen failed. The
structural stop rule is a feasibility filter, not a timing or cryptanalytic
claim.
