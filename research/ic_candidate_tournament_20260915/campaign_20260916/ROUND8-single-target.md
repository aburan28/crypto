# Round 0008 pre-registration: the single-target win under the merged protocol

Written before any measured stage of this round was run.

## Why a further round

Rounds 0006 and 0007 on the single-target panel were frozen and audited with
the evaluator of their day. Meanwhile round 0006-batch16 landed on `main` with
a stricter admission rule: every target's logarithm must carry the descent
relation `[a]G + [b]Q = Σ P_i` it was derived from, and the checker verifies
that relation in the group and the scalar as its consequence under the
verified column logs (`certified_descents`). The round-0007 winner is an
index-calculus solver in exactly that sense, but its frozen output did not
carry the relation, so its result cannot be re-verified by the merged checker.
This round re-measures it with the certificate, under the merged evaluator and
checker, so that the single-target win stands under the current protocol.

## Objective and gates

As round 0007: `--objective rho`. Promotion requires no regression against
the incumbent and, on confirmation and replay, candidate/rho upper paired 95%
limits and every cell below one in both instructions and native process wall
(`beats_rho_strict`). Since the only challenger is a configuration control, the
expected outcome is *retained*; the decision then records the incumbent's own
`winner_over_rho`, `beats_rho_strict` and `rho_parity`, which is the
measurement this round exists for. A promotion of the control would be reported
as such.

## Parent, seed, budget

Incumbent: the round-0007 winner `tiny2` (frozen
`runs/round-0007/source_candidates/tiny2/source`) plus the certificate:
`solve_target` returns the witness it used and the worker writes it as
`"relation": {"a", "b", "points"}` per solution, the same shape as the
round-0006-batch16 baseline. Nothing else changes ([patch](round8-tiny2_cert.patch)).
A module test checks each witness in the general arithmetic and checks the
scalar as its consequence. Configuration `batch_trials: 1`. Challenger
`tiny2_cert_batch4`: the same executable with `batch_trials: 4`. Fresh seed
2026091608; target count 1; pilot profile; 1,356 paired jobs of an 1,800 budget;
one pinned CPU; 8 GiB cap; 60-second watchdog; Valgrind 3.22.0 `Ir`; native
progress recorded. Rho is the shipped per-target signed-Frobenius solver on the
incumbent's executable.

Development evidence before freezing (not a claim): on the round-0007
confirmation fixtures all 60 certified outputs pass the merged checker with
`certified_descents = 1` and the incumbent's factor-base fingerprints.

## Boundary, floor, class, honesty

Unchanged from round 0007: same unit and boundary, base support, `m = 3`, no
direct relations, every verification obligation, worker phase dumps, fixture
construction, rho branch and general-arithmetic final check. Class:
engineering. Rho is the shipped implementation; the native ratio is small and
specific to this virtual machine's process and CPUID costs; no
arithmetic-complexity, family-wide or cryptographic-size claim follows. Fresh
fixtures from the new seed; every failure retained; panels stay separate.
