# n37 four-policy equal-useful-support gate

Decision: **`EQUAL_USEFUL_SUPPORT_ADMITTED_FOR_N37_K42`**. The [preregistered
protocol](PROTOCOL.md) was committed as `bc19274f0e5a20c3727006de0beddfec4ac7297c`
before the new census. Four target-blind policies now have the same 42 log
columns, 1,554 distinct signed subgroup classes, and 3,108 physical points.
An independent native replay recomputed every point, label, action and
exceptional control. This admits the fixed inputs to a natural-target PDP
comparison. **PDP yield, relation rank, complete cost, `S` and IC/rho speedup
remain unset.** No target file or fixture scalar was read by this gate.

The source base is the exact n37/L1,024 compact-orbit base from the archived
v2 b00 process, not a new base selected to fit the leaf. Its 127,620-byte
[header](SOURCE42.jsonl) has SHA-256
`0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75`
and internal base hash
`8423b135df3515b0e284126d9eeb95f71e3901ae12908b3cad8eb03d1f31d4bc`.
The fixed 42 leaf seeds come from the prior
[`NATIVE42.json`](../n37_native_basis_bridge_20261002/NATIVE42.json), and
the degree-73 map is rebuilt from the prior archived certificate. All input
and code digests, including the earlier preserved runs, are in
[`EVIDENCE.json`](EVIDENCE.json).

| Policy | Curve | Log columns | Signed classes | Physical points | Candidate construction action recorded by this gate |
| --- | --- | ---: | ---: | ---: | --- |
| Original | source | 42 | 1,554 | 3,108 | 1,554 Frobenius steps, including closure checks |
| Transported | leaf | 42 | 1,554 | 3,108 | 1,554 degree-73 maps of positive source members; negatives are derived by sign |
| Descendant-native | leaf | 42 | 1,554 | 3,108 | 1,512 scalar actions `[lambda^k]` beyond the 42 seeds |
| Pullback | source | 42 | 1,554 | 3,108 | 42 exact rational inverses, then 1,554 source Frobenius steps |

The leaf's scalar action is **not** being credited as a cheap geometric
Frobenius. The producer explicitly evaluates it 1,512 times. The audit also
performs 1,554 extra maps to check transported negatives and 3,150 maps to
check every native/pullback positive and negative plus each seed; those
cross-checks are not candidate construction costs. The producer checked all
12,432 physical points for curve and order membership. Its independent
[`REPLAY.json`](REPLAY.json) verified all 6,216 full-point map/sign pairs,
6,216 signed classes, 42 exact seed pullbacks and 1,554 native scalar-action
pairs. The replay derives source and pullback orbits with the general
Koblitz group law, then derives native leaf members by mapping the pullback
and compares them with scalar multiplication. It does not trust the
producer's point or label arrays. Infinity and the rational 2-torsion point
passed forward and inverse exception controls and are excluded from all
bases.

The final [`RESULT.json`](RESULT.json) is 2,611,829 bytes, SHA-256
`05664103a6dab090a8b4298f34c0964e6f253f625cd73d661562d85a2fce87af`.
The replay receipt is SHA-256
`758159afc85d989d2cb95008c8b0c0e2ac4c411ca9ae14cfed0105b40dc80dd2`.
Changing the first original point's x-coordinate in a scratch copy made
the replay exit nonzero with `original: bad full point or label at 0/0/1`.
Two prior PASS receipts are preserved compressed: the first preceded the
explicit negative-map checks, and the second preceded retained-memory
accounting. Their canonical point-and-label arrays have the **same SHA-256**
as the final array, `f1e8cb50eb8f191291d979777afde410bd055390f7ed42909158e0a26068b438`.
Neither earlier receipt is used for the decision.

The local [host](HOST.json) was an Apple M4 Pro on macOS, at L0. Descriptive
phase readings in the final manifest include 1,512 leaf scalar actions in
4.72 ms, 1,554 transported positive-point maps plus 1,554 sign audit maps
in 51.43 ms, and 42 exact inverses plus full pullback audits in 994.20 ms.
They are not isolated wall-time measurements, not pure operation prices, and
do not include a cold source-base scan or leaf-seed selection. The point
arrays retain at least 149,184 bytes and producer RSS before JSON output
peaked at 18,546,688 bytes; neither is the final attack's memory budget.
The generic floor for the matched source curve is
`sqrt(pi/(2*74)) = 0.145695` in group-addition equivalents per `sqrt(r)`, but
this stage has no calibrated `S` with which to divide it. The reference in
the next full attack panel remains same-Q signed-Frobenius rho.

The decisive follow-on is a frozen natural-target three-summand PDP test on
the point-only b03 and b04 blocks named in the protocol, followed by rank
and cold logarithm recovery. The isogeny is bijective on this rational
subgroup, so original and transported must have identical exact hit/miss
and rank trajectories after target mapping; native and pullback must too.
That pairwise equality is a strong implementation control. Only the
original-versus-native seed selection can change exact support. A failure
to find a witness when one exists is a solver failure, not an algebraic
advantage of either curve. The existing blocks are held out from this base
gate, although prior source compact work has used them; any speed claim
needs a newly frozen disjoint Q stream and complete cold accounting through
the native ecbench harness. The changed endomorphism order on a degree-263
descendant motivates this action-cost test, but this n37 result does not
predict degree-263 PDP yield or an ECC2K-130 crossover.

Reproduce the arithmetic audit on a checkout containing the archived
inputs and frozen lockfile:

```sh
cp research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock Cargo.lock
cargo build --release --locked --example n37_four_policy_support --example n37_four_policy_support_replay
target/release/examples/n37_four_policy_support_replay \
  research/notes/ecc2k130/n37_four_policy_support_20261003/RESULT.json \
  /tmp/n37-four-policy-replay-fresh.json
cmp /tmp/n37-four-policy-replay-fresh.json \
  research/notes/ecc2k130/n37_four_policy_support_20261003/REPLAY.json
```

The [CI workflow](../../../../.github/workflows/n37-four-policy-support.yml)
performs that replay on Linux. To rerun the producer, give it a new output
path; it refuses overwrites. The policy array digest can be checked with
`jq -c '.policies' RESULT.json | shasum -a 256` from this directory. Its
phase times and host field will differ across runs, while point arrays and
labels must agree exactly.
