# Preregistered full-Frobenius-orbit relation-stage control

The solver probes above fix `(0,1,2)` or restrict each of three phase
selectors to three choices. Neither searches every Frobenius phase. To
measure this omitted factor on identical m83 quotient-base instances,
preregister an exhaustive meet-in-the-middle (MITM) relation-stage control:

- Reuse the original `s=2` quotient factor-base construction, E0 curve,
  prime subgroup, modulus, seeds `260938/260939`, deterministic generator,
  first two proper planted triples and four natural subgroup targets in their
  original order. Compare target coordinates and input SHA-256 to the first
  run before interpreting outcomes. Do not use target scalars in the search.
- For every subgroup lift of every nonzero payload, include both signs and
  all 83 Frobenius phases; deduplicate exact point collisions, preserving
  the corresponding representative, coefficient and phase. Verify every
  orbit point and Frobenius eigenvalue in the subgroup.
- Build all unordered pair sums, including repeated equal points; index by
  group element. For each target test every signed-orbit third point and
  look up the matching pair; reject opposite-point cancellations, verify
  group addition, collect relation rows, identify duplicate rows and compute
  natural-target rank over the prime subgroup. Record planted relations as
  positive controls and do not count them as natural rank. Count every
  group addition during orbit setup, pair building, target search and final
  verification separately, including failed trials. Preserve setup time,
  memory, timeouts, zero yields and verification outcomes.
- One child process per seed, 1 GiB address-space limit and 70 s wall limit.
  If it cannot finish, retain its partial progress receipt and stop. This
  is a tiny-factor-base, three-point, full-orbit relation-stage control;
  do not infer actual DLP success, end-to-end speedup, matched-rho `S`, or
  ECC2K-130 scaling from its completion.

The full rho reference would require an unfeasible number of group additions
on this ~81-bit prime subgroup; keep its measured cost and every full-DLP
speedup null. This control can quantify phase coverage and relation yield on
its sampled natural targets, and can falsify an overbroad fixed-phase claim.
