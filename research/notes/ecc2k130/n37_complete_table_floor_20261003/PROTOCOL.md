# n37 complete pair-table cold setup floor: preregistered gate

## Question and scope

Can the **specific complete at-most-three-summand oracle** used by the four
equal-support n37 policies remain plausible for a cold, one-public-target
ECDLP solve? This is a necessary-cost gate, not a complete IC run or an attack
speed claim. It does not test sparse, streaming, algebraic, or larger-field
PDP oracles. The four policies are original source, transported leaf,
descendant-native leaf, and pullback source; exactly one policy is used per
cold solve. Their inputs are fixed by the merged support and PDP gates.

The source curve is `icv1-f2m37-tm534059-32aad96b`, `K_0/F_(2^37)`, with
prime subgroup order `r = 230603167` and signed-Frobenius class size 74. The
generic rho floor is `sqrt(pi/(2*74))` in `S = GAE/sqrt(r)`; the matched
reference is `rho.signed_frobenius_strong` on the same source subgroup. This
gate will compare **charged group operations only**. The host cannot admit
wall-clock figures at L0; no measured speedup is proposed.

## Frozen inputs and calculations

- Four-policy support:
  `research/notes/ecc2k130/n37_four_policy_support_20261003/RESULT.json.gz`,
  SHA-256 `8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5`.
  Its source 42-column base is pinned by SHA-256
  `0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75`;
  native seed manifest SHA-256 is
  `bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c`.
- Every policy has 42 log columns, 1,554 signed classes and 3,108 physical
  subgroup points. The complete pair table enumerates all nondecreasing
  pairs, including diagonal doublings: `M = P(P+1)/2 = 4,831,386` operations
  for `P = 3,108`. One addition or doubling is one GAE. The frozen PDP
  producer used precisely this enumeration for each policy; the new Rust
  counter must independently replay the additions with `CountedGroup`,
  validate the support hash and point count, and give the same `M` for all
  four policies. All omitted factor-base generation, isogeny construction,
  sorting, relation collection, rank, descent, and recovery costs are
  nonnegative, so `M` is a lower bound for a **cold implementation that
  constructs this complete pair table**.
- The rho spec is [`rho-public-spec.json`](rho-public-spec.json): 8 distinct
  public targets derived by `ecbench` from target seed `202610030037`, one
  warmup and 2 measured rounds, algorithm seed `373737`, strong rho defaults
  (`lanes=32`, `dp_bits=4`, `step_cap_factor=2000`). Each workload is a
  separate one-target DLP. The spec requests L2 to prohibit wall-time
  admission on this host; operation counts remain valid at L0. No planted
  scalar is passed to rho. `ecbench` independently checks `[d]G = Q`.

## Decision rule fixed before measuring

First require that `ecbench plan` resolves the registered source curve and
the same subgroup order/generator as the frozen n37 source; that every
measured rho run verifies; that `ecbench verify --replay-all` reproduces
every deterministic measured run; and that the counter's four policy
receipts are identical under rerun. Keep failures and exhausted runs.

For each verified measured rho run, divide the fixed `M` by that run's
**total charged GAE**, including its setup and internal verification. If
the minimum ratio over all 16 runs is at least 100, reject this exact
complete-table oracle as a rational cold single-target implementation at
n37: its mandatory table construction alone uses at least 100 times the
matched rho's charged operations on every sampled target. If any run fails
or the minimum is below 100, report the gate as inconclusive. Report the
per-run ratios, not a pooled IC/rho speedup. The result is an operation-count
diagnostic; rho's unpriced native work and both methods' wall time still
prevent a runtime comparison.

The 8 workloads, 2 measured rounds, seeds, policy inputs, counted unit,
ratio threshold, and treatment of failed runs cannot change after the
session starts. Do not claim a scaling law from this single curve size.

## Execution and evidence

Open the protocol PR before running the new counter or rho session. Then
build the release Rust producer and `ecbench`, run the counter twice to
separate output files, run the spec into a new sealed `sessions/` directory,
and audit it with `ecbench verify --replay-all`. Commit the producer, exact
commands, source hashes, both counter outputs, session, audit receipt,
per-run comparison, host facts, result and decision. Update the canonical
scoreboard in that same PR. Retain failures and timeouts. Any larger cold
four-policy DLP implementation remains a separate experiment and must use
the native `ecbench` harness for a full method comparison.
