# Preregistered n37 Frobenius-folded pair-table gate

Status: **preregistered, unmeasured**. Commit this protocol and open its PR before
running the new producer. The curve is the registered
`icv1-f2m37-tm534059-32aad96b`, with prime subgroup order `r = 230603167`.
This gate tests a specific way to reduce the table cost exposed by the
complete-table cold floor in [PR #1280](https://github.com/aburan28/crypto/pull/1280).
PR #1280 remains pending CI and is not treated as accepted evidence here.

## Question and frozen inputs

Can the existing `FrobeniusPairTable` represent the **same exact m≤2 and
m≤3 point-sum decisions** as the complete pair table on each frozen source
support, while making only one pair addition per signed-Frobenius pair orbit?
The two supports are `original` and `pullback` from the merged
[equal-useful-support gate](../n37_four_policy_support_20261003/RESULT.md).
Each has 42 log columns, 1,554 signed classes and 3,108 physical points.
The paired leaf supports are excluded: their fast action is a scalar
`[lambda]`, not the source curve's geometric Frobenius, so applying this
table to them would conflate distinct operation costs.

- Read only the immutable support manifest
  `research/notes/ecc2k130/n37_four_policy_support_20261003/RESULT.json.gz`,
  gzip SHA-256
  `8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5`,
  expanded SHA-256
  `05664103a6dab090a8b4298f34c0964e6f253f625cd73d661562d85a2fce87af`.
  Reject changed schema, curve, order, array lengths, point or column labels.
  Verify uniqueness, negation closure, signed-Frobenius closure, curve
  membership and prime-subgroup membership before any target query.
- Use the ordered **point-only** b03 and b04 n37/L1,024 target files, with
  SHA-256 `84400a2914f06e4a001d0f113f0195952f692634d4ff8285e23ade2599e2bde2`
  and `72c5b5361a0b846c9626269ce3928820cb7d420aae479d31dc26f5cbf7e78aa6`.
  Neither producer nor replay may read fixture-scalar files. These 2,048
  targets were used in the earlier PDP decision and are **regression
  controls**, not a fresh policy-selection or speed sample.
- Compare hit/miss for both arities against the independently replayed
  [complete-table manifest](../n37_four_policy_pdp_20261003/RESULT.md),
  `RESULT_CI_FIXED.json.gz` SHA-256
  `c84066a521137cfa61a2f490ed18a3a957b79da6e27f587ad1579b350cb9764a`.
  The folded table may return a different valid witness; verify every
  returned full-point sum and its `(column, coefficient)` relation row.
  The m≤2 decisions include many proved misses and are the main exactness
  control; the old m≤3 decisions all hit on this sample.
- Build against frozen `Cargo.lock` from
  `research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock`, SHA-256
  `b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365`.

## Frozen method and accounting

For each support, construct `FactorBase::from_column_map` in manifest order,
then call the existing native `FrobeniusPairTable::build`. It must report
exactly 42 representatives and **66,822 counted group-addition requests**:
`(2·37)·42·43/2`. A complete nondecreasing pair table requests 4,831,386
additions, so the predicted setup-addition reduction is
`4831386/66822 = 72.30232558`. This is a testable count, not a predicted
wall speedup. Preserve the actual unique table entries, setup wall time,
minimum retained bytes, canonicalisations, Frobenius maps and lookups.
The latter operations and field work must be reported as **unpriced native
work** until a matched calibration charges them.

For an at-most-`m` decision, test the identity and direct factor-base lookup
first; for `m=2` then call `decompose_mitm_frobenius` with `m=2`; for `m=3`
try the same at-most-two path before calling it with `m=3`. This handles
zero, one and two summands explicitly, including identity pair sums omitted
from the folded table. Do not consult the historical manifest during search.
Count query
group additions, canonicalisations, lookups and mismatch counters separately
by support and arity. A returned witness must contain at most `m` indices,
sum to the target under the general source group law, and yield a relation
row from the frozen labels. A missing witness is called a **proved miss on
this query only if the folded-table construction completed and the oracle
performed its entire prescribed scan**. Allocation failure, closure
failure, timeout and counter mismatch are errors, not misses.

The producer writes a new JSON receipt without overwriting prior output.
An independent native replay must rebuild each support from the pinned
manifest, independently recompute all reported witness sums and rows, check
all 8,192 hit/miss decisions against the complete-table reference, and
recompute table and query counters. It must fail after mutation of a
witness index, hit flag, support coordinate or setup-addition count. Commit
the producer, replay, compact raw receipt, input and output hashes, host
facts, mutation receipts, analysis and decision to this PR. CI replays the
receipt on Linux. Do not erase failed or superseded runs.

This uses a **frozen selected factor base**. The original source-base scan,
descendant leaf-seed selection, isogeny transport, exact pullbacks and
their validation are not charged by this stage. Table memory, native
canonicalisation, field operations, relation collection, linear algebra,
target-log recovery, and matched rho are also outside a complete cost
comparison. Record `S`, IC/rho ratio and end-to-end speedup as **null**;
there is no ECC2K-130 attack-speed claim from this gate.

## Decision rule and next experiment

Pass only if both supports build with the predicted count and zero closure
or pair-image mismatch, all m≤2 and m≤3 decisions agree with the complete
reference, every hit has an independently verified full-point witness and
row, both runs agree bit-for-bit on non-timing evidence, and independent
Linux replay succeeds. Failure or unknown remains a retained negative result.
If it passes, register the fixed support as an explicitly uncharged-input
`ic.pipeline` factor-base plugin and run a **fresh**, point-only n37 Q panel
through full-rank recovery against same-Q signed-Frobenius rho. Then
replace the frozen input with a cold native support constructor and charge
all selection and native operations before any crossover claim. Test n41
and n53 only after the n37 accounting is sound. The leaf needs a separately
priced scalar-action orbit table or transport to the source before a fair
comparison. A failed n37 single-target cost gate does not by itself rule
out batch amortisation or a larger-n crossover.
