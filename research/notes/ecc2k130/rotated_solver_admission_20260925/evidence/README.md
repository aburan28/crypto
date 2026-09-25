# Static solver-admission receipts

`final/` is the accepted frozen interface audit. Its `receipt.json` records
source commit, `FROZEN.json` hash, child command, UTC interval, exit code,
child peak-RSS upper bound and SHA-256 of every raw file. `result.json` records
all pinned corpus/source facts, missing prerequisites and the local binary
inventory. The stdout/stderr files are retained even though empty. The result
SHA-256 is `7d1acc98e428af3674e3f3d875f223a825dbdd66e0b8a8702e0bb93c6857dc8b`;
receipt SHA-256 is `7e2eeb25931d8bb4aca2f2289e3fd4de1445c91b8d945f4ee69dbff3a646ce17`.

`path_only_initial/` preserves the first successful static audit from draft
head `c901d44`. Its inventory listed executable names on PATH but omitted
the #764 WDSat fixture binary stored outside PATH. No solver result changed;
the original source/receipt remain reachable at that head. The corrected
source was committed and its hash-only CI passed before the new final child.
The original directory is historical raw evidence, not a final inventory.

Reproduce the accepted archive replay from a checkout containing the pinned
sources and evidence:

```sh
python3 research/notes/ecc2k130/rotated_solver_admission_20260925/ci_replay.py --evidence research/notes/ecc2k130/rotated_solver_admission_20260925/evidence/final
```

This reruns only the static source/corpus audit without querying installed
binaries. The final receipt provides the binary path/version/hash observations
from the original host. Neither the initial nor final child invoked a solver
on a PDP instance, so these files cannot be used as solver timing data.
