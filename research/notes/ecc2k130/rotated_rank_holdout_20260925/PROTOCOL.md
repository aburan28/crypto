# Preregistered modular-rank and point-only holdout gate for rotated rows

Status: protocol only; no rank, row-selection or recovered-scalar outcome has
been run. This follows merged [PR #766](https://github.com/aburan28/crypto/pull/766)
(commit `e762852591ec6cb147707daac9127b22e5150143`). Its immutable
compressed-row archive SHA-256 is
`863f50ea872ce01b1999d0ef67bc0341fa0c764f2ad4187d16736768b0ab9b3d`
and its frozen source manifest SHA-256 is
`e08fd485952cc8a3674964ed9b9d95f0f8fa1ba6978b9c0c8cc2713974eb5ae4`.
The source curve is E:y²+xy=x³+1 over F_(2^13), polynomial `0x201b`,
subgroup order q=2003, cofactor four, H=(4793,2429), Frobenius scalar
λ=89. Four rotated arms are beta 3/7 × m 5/6. PR #766 established row
semantics on 17,917 positive Q+T witnesses but did **not** count their
modular rank or recover held-out logs. This gate asks those two separate
toy questions, without a PDP solver or n=131 speed claim.

## Frozen split and complete cost

Before any rank or solve outcome, generate one literal input manifest from
the #766 archive. For each arm, score each nonzero k∈{1,…,2002} by
`SHA256("ECC2K130-ROTATED-RANK-20260925-v1/beta/m/k")`, sort by digest
bytes then k, inspect the exact #766 four-coset roster in that order, and
retain the **first 16 projected-positive k** (at least one saved Q+T row).
Preserve every scanned negative and positive candidate in the selection
transcript; no favorable target may be added after outcomes. Freeze the
selected k labels in a sealed sidecar and a distinct `point_only.json`
containing only opaque case IDs and full Q coordinates. Freeze source, all
literal inputs and SHA-256s, caps, commands and a draft PR before solving.
The selection rule targets positive support conditionally; its success rate
is not a natural-target yield estimate. It uses the complete old archive as
an oracle input and charges archive loading and scan.

Training excludes **all** archived rows whose Q equals any selected holdout
point, across all four torsion cosets. Each remaining row gives
`Σ_j c_j log_H(C_j)=4k (mod 2003)`, where C_j is its canonical [4]F0
column and the signed c_j were independently certified in #766. Define
the full candidate column set from the frozen arm summary's
`term_column_coordinates` (not just columns seen after selection), sorted
lexicographically. Stream archived rows in their recorded `(k,T)` order;
include zero rows and count dependent rows. Perform exact modular Gaussian
elimination with fixed leftmost pivot order, verify every row is consistent,
record first rank attainment, total rank, nullity and pivot transcript.
If rank is deficient, retain the failure and do **not** invent unique base
logs or call a target recovered. Full rank permits a unique base-log vector;
check every training equation independently before any holdout query.

The point-only recovery child receives `point_only.json`, solved base logs,
and a separately built **archive-oracle** response keyed by Q coordinates.
It receives no k, training witness index or sealed label. The oracle builder
scans all #766 row records, and for each Q picks the first available
T-index row, strips k, source indices and all labels, and records its signed
column coefficients, Q coordinates and torsion index. Charge the full scan,
index construction, lookup and misses. This is a precomputed archive witness
lookup, **not** an implicit PDP solver or a new relation yield. For an
available oracle row with known columns, recover
`k̂ = 4^(−1) Σ_j c_j log_H(C_j) (mod 2003)`. Independently verify
`[k̂]H=Q` by full group law for every returned scalar, then open the sealed
labels and compare `k̂=k`. Preserve all 16 cases per arm, including lookup
misses, rank-deficient skips and verification failures. Never use a saved
row's k as a recovery input.

## Independent replay, caps and decision

Use PR #766's Euclid-inverse arithmetic for the producer and the separate
bit-serial/Fermat field/group arithmetic for independent replay. The replay
must rederive the SHA holdout ordering and complete exclusion set from the
immutable archive, independently recompute rank/nullity and a solution,
check every retained row equation, rebuild each oracle response from the
original source row, verify each returned scalar by group law, and compare
sealed k only after group validation. Archive raw inputs/results, stages,
zero/dependent/independent row counts, exact source/input hashes and a failure
receipt. Enforce 120 s/512 MiB for the producer and 300 s/512 MiB for the
independent replay; wall and CPU/RSS are host diagnostics. Record field and
point operations, cold archive/oracle setup, misses, rank and recovery
separately. Preserve partial files and a quantitative failure receipt on a
cap or assertion failure.

Pass only if source/input hashes, full holdout exclusion, modular row
consistency, full-rank solution (if one exists), point-only recovery process
isolation, every `[k̂]H=Q` and sealed-label match all replay independently.
A rank-deficient arm is reported as a censored or failed recovery, not
selectively discarded. This is a tiny-field semantic/end-to-end toy gate;
it cannot establish n=131 base cardinality, PDP search cost, relation yield,
independent n=131 rank, a Certicom challenge log or an IC/rho crossover.
The next n=131 admission needs a separate exact compressed-base cardinality
bound, including projected duplicates and x=1 torsion, and a fully charged
solver/memory experiment. Update the canonical scoreboard and decision
ledger only after an accepted archive-only replay.
