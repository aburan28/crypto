# Post-registration source feasibility audit

This is a source-level finding made after the one permitted v2 dispatch, not a
measurement or a change to its frozen protocol. [Actions run 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479)
ended [operationally censored](RESULT.md) without a retained bundle; the static
audit does not substitute for that missing empirical record. Its seed must not
be retried. The
machine-readable [audit](STATIC-FEASIBILITY.json) has SHA-256
`a492ef84c38cf77a1d12d2e489d29fab18a27d8c7a94cf9aa645d2ea4342cbb2`.

The pinned worker source at `765c3c5f19032bd852163805f257c56babef2040`
materializes `subgroup_orbits` with an **ambient basis of length `ell=n`**.
All three algebraic engines (`f4`, `f5`, `inherited_f4`) enter
`groebner_decompose`, whose reused Semaev template needs
`m*ell+(m-2)*n` Boolean variables. For the registered `m=3`, this is `4n`.
`koblitz_groebner::MAX_VARS` is 64, and the template returns `None` above
that limit. `groebner_decompose` turns `None` into an `unsupported` attempt
before F4/F5 solving. The five cells therefore all fail this deterministic
encoder condition, regardless of their target, query, solver budget, or LA:

| Cell | Basis length | Variables | Cap | F4/F5 layout |
| --- | ---: | ---: | ---: | --- |
| n17a1 | 17 | 68 | 64 | unsupported |
| n19a0 | 19 | 76 | 64 | unsupported |
| n23a0 | 23 | 92 | 64 | unsupported |
| n23a1 | 23 | 92 | 64 | unsupported |
| n31a0 | 31 | 124 | 64 | unsupported |

The SAT arms use a separate wide S4 encoding, so this check does not decide
their outcome. The running campaign's actual receipts still determine which
processes completed, failed, timed out, or were censored. No family speedup or
qualification follows from this static audit, and its results must not be
filled in as measured zero yield.

[`generic_solver_feasibility.py`](../../generic_solver_feasibility.py) checks
the exact reviewed source objects and must run against the worker checkout
*before* a future panel generates any fresh fixture. With `--require-pass`,
it exits nonzero if an algebraic arm exceeds the cap. A pass means only that
the layout fits; it does not certify dispatch, relation yield, correctness,
runtime, or scientific admission. For a new source revision, review and update
the source objects and bound before using the gate.

## Replacement base control, not a selected candidate

The pinned source also supports `standard_subspace`, whose actual algebraic
basis has a chosen length. At dimension 6, the five layouts need 35, 37, 41,
41, and 49 variables, respectively. The independent Python constructor now
rebuilds the exact polynomial-basis abscissae and point order. A local worker
inventory on one **previously disclosed** target per cell agrees with that
constructor and its independently computed subgroup census. The retained
[raw reports and replay receipts](standard-subspace-d6-inventory-control.json)
record the worker checkout and binary digest. This is a base-construction
control only; the worker build is not a qualified measured build.

| Cell | Geometric points | Usable subgroup points `B` | Folded columns |
| --- | ---: | ---: | ---: |
| n17a1 | 63 | 62 | 29 |
| n19a0 | 65 | 62 | 27 |
| n23a0 | 75 | 72 | 33 |
| n23a1 | 53 | 52 | 23 |
| n31a0 | 69 | 66 | 27 |

`standard_subspace` is generally not Frobenius-closed, whereas the incumbent
sampled-orbit base is. Any complete-pipeline comparison that changes to this
base is a **factor-base-policy comparison**; it cannot isolate the algebraic
solver as the cause of a cost difference. Before registering another fresh
panel, use a source-bound worker on disclosed points to check actual F4/F5
dispatch, ordinary-query outcome mix including failures, novel rank and
bounded complete recovery. Retain zero-yield cells and limits. If that pilot
justifies another registration, use a new seed, exclude every v2 target that
could have been generated, and keep the three sealed confirmation sets closed.
