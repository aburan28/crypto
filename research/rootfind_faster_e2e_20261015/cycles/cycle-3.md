# Cycle 3 — ecbench smoke + prime e2e after probe wins

**Change.** No new engine id (still `rootfind-negfold`). Added smoke spec
`docs/ecbench/specs/smoke-koblitz-rootfind-faster-e2e-20261015.json` and
measured the post-cycle-1/2 binary under fail-closed ecbench + prime e2e.

**Hypothesis.** Complete negfold still ≪ XOR-priced F4 in `S`; cold vs
rho-strong remains ~1.5–1.7× (hard stop). Prime rootfind still ≫ msolve
in CPU time on 16-bit.

## Measure

### ecbench smoke (n=13, 2 public × 2 rounds, L0)

Session `ECBS1hd5339a8cb067`, audit OK (18/18 verified, 8/8 replays identical).
Receipt SHA-256 `5723fd6322be993fa37183c7712a2638a4530e1be72ced6349ee7d92f93ac088`.

| arm | mean S | S / rho | lower bound | med wall ms |
|---|---:|---:|---|---:|
| rho-strong | 14.542 | 1.000 | yes | 0.128 |
| **ic-negfold-m2-d5** | **24.609** | **1.692** | no | 1.615 |
| ic-f4-m3-d5 | 418.599 | 28.786 | no | 126.078 |

Paired: `negfold / rho` = **1.692** [1.404, 2.230] (unpriced; diagnostic).
`negfold / F4` = **0.0588** [0.046, 0.070] (complete beside complete).

### Prime e2e (16-bit, msolve paired CPU)

Root-finding **31.71 ms** measured vs msolve composed **1.72 s** → **54×**.
At each method's best `|F|`: root `|F|=120` 8.60 ms vs msolve `|F|=12` 1.38 s → **161×**.

## Outcome

**Confirmed.** Probe-path engineering did not break ecbench identity or
completeness. Hard stop vs rho stands (~1.69× on this smoke; inherited
~1.54× on larger claim panels). Negfold remains ≪ F4. Classification:
**engineering**.
