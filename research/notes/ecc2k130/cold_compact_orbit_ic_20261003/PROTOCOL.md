# Preregistered cold compact-orbit source pipeline gate

Status: **preregistered, unmeasured**. Commit this protocol and open its PR
before a new `ecbench` session. It follows the fixed-support Frobenius-table
exactness result in [PR #1281](https://github.com/aburan28/crypto/pull/1281),
which is pending CI and is not yet accepted as a merged result. The purpose
here is to replace that experiment's free archived source support with a
fresh, counted constructor and to recover new public logarithms end to end.

## Boundary and question

For `E_0: y²+xy=x³+1` over `GF(2^37)`, use registered ICV1 slug
`icv1-f2m37-tm534059-32aad96b`, prime subgroup order `r=230603167`, and
signed-Frobenius order `A=74`. The generic-group floor in the repo's unit is
`S_rho,floor = sqrt(pi/(2A)) = 0.1456948091`. The **reference** is
`rho.signed_frobenius_strong` on each exact public Q in the same native
`ecbench` spec/session. A 42-column base has 3,108 signed points and its
orbit-folded pair table still has the unavoidable setup-addition count
`37·42·43 = 66,822`, hence `S_setup,adds >= 4.4003461`, or 30.2025 times
the generic rho floor before base selection, queries or rank. This floor
does not predict any individual rho walk or a batch speedup.

Question: can a **fresh cold source base**, counted folded table, relation
collection, full-rank solve and target recovery be reproduced and charged
beside strong same-Q rho, and does the table setup already rule out the
single-target route at n37? No wall or full-speed claim is allowed if any
native field operation remains unpriced.

## Frozen construction and correctness controls

Implement an `ic_framework` `compact-orbit-scan:columns=42,raw_x_cap=1000000`
factor-base plug-in and expose it to `ecbench fb` and `ic.pipeline`.
Its source algorithm is the merged compact-orbit producer's deterministic
scan: start raw x at one; lift both points on the Koblitz curve in the
constructor's order; clear the subgroup cofactor; discard infinity and
projected x≤1; key a candidate by its least Frobenius-conjugate x; accept
the first 42 distinct keys; for each accepted seed append its 37
coordinate-Frobenius images, each followed by its negative, with labels
`(column, ±lambda^j mod r)`. Stop with an error at the cap, with no
substitution of a later seed or target-dependent choice. Count every
cofactor multiplication and order check through `CountedGroup`; retain
the raw-x scan, lifts, inversions, Artin–Schreier solves, field squares,
Frobenius maps, hash operations and memory as explicit native counters.
Price only units with a pinned curve calibration. Mark the rest
`*_uncharged` and preserve the run as a lower bound. Refuse duplicated,
off-curve, wrong-order or nonclosed points.

Before any fresh-Q measurement, the n37 constructor must reproduce the
**entire** point/label array and all 42 representatives in
`research/notes/ecc2k130/n37_four_policy_support_20261003/SOURCE42.jsonl`,
SHA-256 `0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75`,
internal base hash
`8423b135df3515b0e284126d9eeb95f71e3901ae12908b3cad8eb03d1f31d4bc`.
That source base was selected before this experiment. The immutable
four-policy support manifest has gzip SHA-256
`8eac2ae4b8d8f0fc452b7cd7c8edbd3183558e55454abbd67349665bd64610d5`.
Use the archived `Cargo.lock` SHA-256
`b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365`.
If exact parity fails, retain the failure and repair the constructor using
only the archived source algorithm; do not inspect target outcomes first.

Give the new counted oracle a distinct plug-in name
`mitm-frobenius-counted:m=3`. It uses the existing folded table and
decomposition logic, but records the table's native setup counters in the
oracle-setup phase. Its control is the existing
`mitm-frobenius:m=3` on the **same freshly built base**; the two arms must
recover identical scalars and have identical group-operation, relation,
rank and query traces. A change beyond native setup accounting is a
correctness failure, not an optimisation. The strong-rho arm is a separate
same-Q reference.

## Frozen first measurement and decisions

Use `ecbench.spec/v1` on n37 with `target_kind=public`,
`targets_per_curve=8`, `target_seed=202610030641`, alternate arm order,
`warmup=1`, `rounds=2`, deterministic session seed `202610030642`, an
L0 operation-count session and a 600-second per-run timeout. Each arm gets
the same eight Q; the IC relation source is `walk` with
`max_trials=1000000`. Store the exact spec and `ecbench plan` output before
running. Do not give either IC arm planted scalars or fixture-scalar files.
Run `verify --replay-all`, require all 16 measured answers **per arm**
(48 in total) to satisfy `[d]G=Q`, every replay to be identical and no IC rank
or target exhaustion. Save the full session, audit, per-Q paired comparison,
phase costs, source/input hashes, host record, failures and decision in this
PR. Use the pinned calibration and disclose every unpriced counter; any
incomplete or nondeterministic pair prevents an IC/rho verdict.

The n37 **single-target cost rejection** is allowed if the minimum counted
IC group-addition-equivalent lower bound across all 16 measured IC runs
exceeds the maximum fully recorded strong-rho GAE on those same Q and
rounds. Report the actual per-Q ratios and range, and name the rejected
route precisely: an eager folded m=3 table with this 42-column source
selection. If that condition fails, classify cost as inconclusive and
inspect which phase changed; never retune the support, target seed, unit
or trial cap after looking. Even if it passes, it says nothing by itself
about a shared-setup batch run, another arity, n41/n53, degree-263 leaves
or ECC2K-130 at n131. `S` and speedup remain bounded or null wherever
field work is unpriced, and L0 wall times are descriptive only.

After this admission gate, retain the same source algorithm at n41/K255
(`r=549756390943`, table-addition floor `2,676,480`) and n53/K440
(`r=21044858204113`, floor `10,284,120`). Their formal three-multiset
count/r ratios are 2.7721 and 0.8035, respectively; those are counting
capacity checks, not observed PDP yields. A separate preregistered
shared-setup batch process must amortize the cold table over new disjoint
public Q and compare with matched automorphism-aware **batched** rho,
including full rank, every recovered scalar, memory and native work.
Descendant-native and transported/pullback arms need their own charged
isogeny and scalar-action construction. No n131 transfer conclusion may be
drawn without the repository's n83 confidence gate and explicit scaling
assumptions.
