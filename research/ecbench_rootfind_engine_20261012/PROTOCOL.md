# Protocol: the root-finding decomposition engine in ecbench

Written before any session ran. The question, the arms, the sizes, the
unit and the reading rules below are fixed; the sessions, audits and
tables that follow them are in this directory's `README.md`.

## Question

`RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md` (PR #1729) replaces the
Gröbner-basis solve of the ECDLP index-calculus decomposition step by root
finding: on binary curves a precomputed pair-sum table
(`DecompositionStrategy::PairTable`), on prime curves a direct-subtraction
sweep with one Montgomery batched inversion a target. That note measured the
*oracle* in isolation, on its own clock. The question here is what the
same engines read as **whole cold methods in ecbench's one unit**, beside
the Gröbner engines already in the harness and beside the matched rho, on
the same targets and seeds.

New engine ids (`docs/ecbench/README.md` §7):

- `pdp3-koblitz:m=2|3,engine=rootfind` on `koblitz-standard-subspace`
  (the base of `pdp3-koblitz:engine=inherited-f4`, the Gröbner control).
- `rootfind-subtract:m=2|3` on a prime curve (the base of `subtract` and
  `mitm`; ecbench has no prime Gröbner oracle).

## Boundary, reference, unit

- Floor `√(π/2A)` in `S`: `A = 2n` on the Koblitz curves, `A = 2` on the
  prime curves, recorded per run.
- Reference: `rho.signed_frobenius_strong` (Koblitz), `rho.negation`
  (prime). `rho.signed_frobenius` is the operations-only Koblitz baseline.
- Unit: `S = charged group-addition equivalents / √r`, cold, every phase
  charged. A root-finding arm's `S` is a complete count. A Gröbner arm's
  `S` is a **lower bound**: its Boolean word XORs are counted, not charged.
  **The ratio of a complete count to a lower bound is not a speedup and
  will not be reported as one.** The engines are compared on
  (i) verification and query-stream identity, (ii) the resource vector
  (`ecbench resources`: additions, lookups, word XORs, wall, RSS) and
  (iii) wall time, which on this macOS host is L0 and therefore
  descriptive only.

## Hypotheses

- H1 (correctness): every root-finding run verifies (`[k]G = Q`, and
  `k = planted` where planted) and every replay is `identical`.
- H2 (same stream): on a base and seed where the inherited-F4 arm does not
  exhaust its node budget, the root-finding arm at the same `m` has the same
  `targets_tried`, `relations_found` and `matrix_rank`; on the prime curves
  `rootfind-subtract:m=2` has the same three counters as `subtract` and the
  same `lookups`, and its relations phase charges exactly `relations_found`
  more additions.
- H3 (constant): the root-finding arms' solve wall time is below the
  Gröbner arms' on every completed workload at `m = 3`; descriptive only.

Stop conditions: a `wrong_answer` anywhere stops the round and is reported
as a defect. `error`, `exhausted` and `timeout` rows are kept and listed.

## Sizes and sessions

| spec | curves | targets × rounds | role |
|---|---|---|---|
| `docs/ecbench/specs/smoke-koblitz-rootfind.json` | `koblitz a=0 n=13`, d=5 | 2 × 2 | plumbing |
| `docs/ecbench/specs/smoke-prime-rootfind.json` | `prime_search 20 bits` | 2 × 2 | plumbing |
| `docs/ecbench/specs/claim-koblitz-rootfind-n13.json` | `koblitz a=0 n=13`, d=5, public targets | 8 × 3 | head-to-head |
| `docs/ecbench/specs/claim-koblitz-rootfind-n19.json` | `koblitz a=0 n=19`, d=7, public targets | 8 × 3 | head-to-head |
| `docs/ecbench/specs/panel-prime-rootfind.json` | `prime_search 18, 20, 22, 24 bits` | 8 × 3 | head-to-head (prime) |

Two Koblitz sizes are not a scaling claim (AGENTS.md §5 needs four); the
n = 23 and n = 31 sizes of the F6 ladder are out of this round's budget
(inherited F4 took 128 s a run at n = 23 on an idle host) and are listed as
future work. The prime panel reuses the curves, targets and seeds of
`research/ecbench_all_candidates_20261003/spec-prime.json` so its
`subtract` and `mitm` rows reproduce that session's.

## Inadmissible

Changing the workloads, seeds, unit or accounting after seeing results;
dropping a failed run; reading a cold wall ratio as an online speedup;
inferring anything about a deployed curve or about ECC2K-130.

## Classification

By AGENTS.md §3 from the ratio to the floor and to the reference, not from
the headline. Expected: *engineering* against `mitm` (the same work in a
faster constant) and *accounting* against the Gröbner arms (a complete count
beside a lower bound); no *advance* is possible from this round because the
asymptotic group-operation count of the decomposition is unchanged.
