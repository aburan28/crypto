# Frozen n37 four-policy natural-target PDP and rank gate

Status: **preregistered, unmeasured**. Commit this protocol and open its PR
before running the new PDP producer. This is the next gate after the merged
[equal-useful-support census](../n37_four_policy_support_20261003/RESULT.md),
not an IC/rho performance comparison or an ECC2K-130 attack-speed claim.
The source is `icv1-f2m37-tm534059-32aad96b`, with prime subgroup order
`r = 230603167`; its archived degree-73 image is the fixed descendant.

## Question, boundary, and frozen inputs

At the **same 42 log columns and 1,554 signed classes (3,108 physical
points)**, does selecting seeds on the descendant change exact natural-target
three-summand decomposition yield or rank relative to source-selected seeds?
Original↔transported and descendant-native↔pullback are algebraically paired
controls. Their hit/miss flags, canonical witness indices, coefficient rows,
and rank trajectories must agree after mapping their target and generator.
Only original-versus-pullback (equivalently transported-versus-native) compares
different seed selections.

- Read the immutable four-policy support manifest only from
  `research/notes/ecc2k130/n37_four_policy_support_20261003/RESULT.json.gz`,
  gzip SHA-256
  `8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5`,
  expanded SHA-256
  `05664103a6dab090a8b4298f34c0964e6f253f625cd73d661562d85a2fce87af`.
  Its independent support replay receipt has SHA-256
  `758159afc85d989d2cb95008c8b0c0e2ac4c411ca9ae14cfed0105b40dc80dd2`.
  Recheck the manifest's four 3,108-point arrays, coefficient labels, point
  membership, and pairing before a query; never substitute a seed or support
  point after seeing target outcomes.
- Read the **point-only**, ordered 1,024-line public targets in block b03 then
  b04 from `research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/`:
  `n37_L1024_b03.points.jsonl` has SHA-256
  `84400a2914f06e4a001d0f113f0195952f692634d4ff8285e23ade2599e2bde2`;
  `n37_L1024_b04.points.jsonl` has SHA-256
  `72c5b5361a0b846c9626269ce3928820cb7d420aae479d31dc26f5cbf7e78aa6`.
  Check all 2,048 points for curve and subgroup membership, then map each
  source target once to the leaf. Neither producer nor replay may open a
  fixture-scalar file. These blocks were held out from support construction,
  although older source-policy work used them; any later speed claim needs a
  fresh disjoint public-Q stream.
- Rebuild the field bridge and oriented isogeny from the source, native, and
  degree-73 archive inputs pinned by the support manifest. Use the frozen
  `research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock`, SHA-256
  `b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365`.
  The source generator is the registered Koblitz generator; the leaf
  generator is its complete degree-73 image.

The exact **at-most-two-summand** population counting ceiling for uniformly
random subgroup targets is at most
`(1 + 3108 + 3108·3109/2)/r = 4834495/230603167 ≈ 2.09646%`.
Collisions and cancellations can lower it. It is a population ceiling, not a
deterministic upper bound on these fixed 2,048 targets. The formal
three-multiset count is `C(3110,3) = 5008536820 ≈ 21.72r`, so its analogous
counting bound clips to one and cannot predict yield. The attack reference is
cold same-Q signed-Frobenius rho in a later common-unit comparison; this gate
sets complete cost `S`, rho ratio, and n131 transfer to **null**.

## Exact oracle and immutable probe schedule

Use all 3,108 **full** point entries, in the manifest's original order, for
each policy. The primary oracle admits zero to three factors; an empty factor
list is the identity. Build a complete table by enumerating the identity,
every singleton `i`, then every nondecreasing pair `(i,j)` with `i≤j` in
lexicographic order. There are exactly 4,834,495 raw candidates per policy.
Key by the full point, retain the first candidate for each sum, and record
both raw and distinct counts. For a target `Q`, scan the empty singleton
then each factor `k` in manifest order and look up `Q−P_k`; accept the first
full-point match after recomputing the at-most-three-factor sum. This fixes a
canonical witness. A **proved miss** requires a complete table and all 3,109
residuals checked. Any allocation failure, timeout, invalid point, mismatch,
or truncated scan is an **unknown/error**, never a miss.

As a secondary control, use the same fixed point set to decide at-most-two
factors exactly. Query identity and each singleton against the complete
identity/singleton set. Record separate witness, hit and probe counts. The
two-summand control is not a replacement for the primary three-summand
comparison, and neither metric licenses seed selection on these targets.
Signed cancellation and repeated factors are allowed; the emitted relation
row reduces all factors by their frozen `(column, coefficient)` labels modulo
`r`, including repeated or negated entries.

For rank, use the same target-blind SplitMix64 stream on all four policies:
start state `0x6e33375f34627064`, advance by the standard SplitMix64 step,
take `1 + output mod (r−1)`, and discard repeated scalars. The first 256
distinct values are the cap. Query `[a]G_source` and its leaf image with the
primary oracle. Verify each witness as a full point before adding its row
`Σ c_{i,k} log(B_i) = a (mod r)` to a 42-column incremental rank tracker.
Stop an arm at rank 42 or after 256 probes; retain misses, dependent rows,
every row, right-hand side, witness and rank transition. If rank 42 is
reached, solve the independent system and verify each of the 42 base logs by
scalar multiplication. Then derive and verify every public target log with
a primary witness on **both** the source and leaf generators. A missing
witness or incomplete rank leaves that target's recovered log null.

## Audit, accounting, and decisions

The producer may use a sorted complete pair table; the independent native
replay must use a separately implemented full-point lookup for the two
source-curve supports and infer each paired leaf result only after checking
the isogeny, every mapped witness, and the prime-subgroup bijection. Replay
all proved misses, witnesses, coefficient rows, rank transitions, base logs,
and recovered scalar checks. Mutating a target point, witness index, row, or
hit flag must make replay fail. Preserve raw outputs, hashes, code hashes,
host facts, failures and an independent replay receipt in this PR. Commit the
raw manifest compressed, with its compressed and expanded digests; CI must
inflate and replay it on Linux. Never overwrite a previous run.

Record each policy's pair additions, table cardinality/collisions, build
wall time and retained bytes; target and rank residual subtractions/lookups,
verification additions, scalar multiplications, map calls, rank operations,
and phase times separately. Charge the native `[lambda]` support closure,
source-base scan, leaf-seed selection, transport, and exact pullbacks only in
the later **fresh cold** process; archived support is an input to this gate.
All wall times here are descriptive host-specific stage diagnostics. Keep
`S`, rho ratio, and end-to-end speedup null.

The correctness gate passes only if all four complete tables, all 2,048
target decisions per arm, all paired equalities, every emitted witness, and
every reported rank/log transition replay. Report source-vs-native paired
discordances and an exact conditional binomial (McNemar) p-value; call a
three-summand **selection lead on these blocks** only for an absolute yield
gap of at least `21/2048` and `p<0.01`, with full verification. If both
selection policies hit all 2,048 targets, classify three-summand n37 yield
as **saturated/non-discriminating**, regardless of different first-hit probe
counts; preserve the two-summand control and move the next yield test to a
larger n. If rank fails to reach 42 by the fixed cap, label full-rank and
log-recovery incomplete rather than tuning the probe stream. If resources
prevent complete tables or scans, preserve the partial record and report
unknown, with no exact-yield or speed claim. A measured cost advantage still
requires the later cold ecbench and matched rho gate, even if this fixed
sample passes the selection-lead rule.
