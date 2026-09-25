# Retained pre-outcome checker-control failure

Before freezing the draft protocol or running any selected #785 negative query, the first local `ci_replay.py` preflight exited 1 at `common.py`'s mutated-proof assertion. The fixed two-variable XOR CNF was `p cnf 2 4` with all four sign combinations. The valid DRAT proof was `1 0`, `2 0`, `0`; the attempted mutant deleted only the final `0`. DRAT-trim still verified that mutant because the earlier derived unit clauses already exposed a unit-propagation conflict. The assertion therefore tested an invalid expectation, not a checker defect or a selected PDP outcome.

The corrected pre-outcome control keeps the claimed empty-clause line but deletes **both derived unit steps**, yielding the mutant `0` alone. On the same non-unit XOR CNF, the pinned checker returns exit 1 and `s NOT VERIFIED` because the claimed conflict is unjustified. The valid proof still returns exit 0 and `s VERIFIED`. This correction was made before `FROZEN.json`, draft PR creation, CI, or any selected proof-generation run. The failed preflight command was:

```sh
/Library/Frameworks/Python.framework/Versions/3.13/bin/python3 research/notes/ecc2k130/n13_oaware_unsat_proof_20260925/ci_replay.py
```

Its recorded Python exception was `AssertionError` at the mutated-control assertion in `checker_selftest`. No selected query or proof was produced by that invocation.
