# Frozen counting admission for the native three-summand chain

Status: preregistered before the calculator and its output. This is an exact
integer **necessary-bound** gate for the full-rank route in [the parent paired
protocol](../../../ecc2k130_leaf_native_m3_gate_20260926/PROTOCOL.md). It
does not enumerate an n=131 factor base, submit an n=131 CNF to a solver,
measure natural PDP yield, or estimate an end-to-end ECDLP cost.

## Hypothesis and fixed inputs

The [native m3 chain admitted in PR #976](../native_m3_chain_20260929/DECISION.md)
uses three ordered disjoint normal-basis x spaces of dimensions `[44,44,43]`
over GF(2^131). A three-slot base with the same **total physical slot-choice
budget** as the admitted balanced m10 structural screen will have a negligible
necessary hit ceiling on uniform challenge-subgroup targets. Conversely, even
a 1% necessary hit ceiling at m3 will demand far more physical choices and
potential log columns. This is a comparison of counting envelopes, not of
solvers or measured support. The m10 screen's tuple count is likewise only a
ceiling and does not certify positive yield or rank.

Freeze parent `origin/main` commit
`ce79b74645c8e564ec361918dfb72bd990cc56d6` and these SHA-256 inputs:

| Input | SHA-256 | Required fields |
| --- | --- | --- |
| `research/notes/ecc2k130/m10_export_capacity_20260925/inputs/balanced_m10_result.json` | `1bf1dc09dec347dcd421266bc4ece4b529931f21874fc3e6a26b5620f7a2cc25` | `q=680564733841876926932320129493409985129`, `physical_f0_points=7977`, `m=10` |
| `research/notes/ecc2k130/native_m3_chain_20260929/chain.py` | `03b3e97b74e35a100ab60a7c5be0cea09108405e185cccda74c6942a94b21944` | `LEAF_DIMS=(44,44,43)` |
| `research/ecc2k130_leaf_native_m3_gate_20260926/admission.py` | `2a03a30c337d5dfb4d0d30ddff102293db6c8461ea1520138ef268ea5c8f0ece` | existing ordered-slot bound for cross-check |

There is no random stream. Freeze `N=1,000,000` uniformly marginal target
probes as an illustrative bounded panel, one relation row per probe maximum,
and necessary support-ceiling thresholds `1/100` and `1/2`. These are fixed
before calculating the output. If the eventual experiment probes torsion
translates or uses a different target distribution, it must separately count
all probes and use its correct uniform coset size. Here `q` is the most
optimistic one-coset denominator; a uniform full-group target uses `4q`.

## Exact arithmetic and checks

For ordered physical slot sizes `s0,s1,s2`, at most `C=s0*s1*s2` distinct
target sums exist. With total physical slot-choice budget `S`, the maximum
product is the balanced integer partition
`M(S)=a^(3-r)*(a+1)^r`, where `a=S//3` and `r=S%3`. Therefore the probability
that a uniformly marginal subgroup target lies in support is at most
`min(1,M(S)/q)`. The expected number of supported probes among `N` targets
is at most `N*min(1,M(S)/q)`, whether or not the targets are correlated, and
the probability of any hit is at most the same value capped at one. These
are upper bounds, never predicted success rates.

Set the fixed budget `S=10*7977`, calculate the balanced slots and exact
fractions for support, expected hits and any-hit probability. For each frozen
threshold `u/v`, compute the **smallest integer** `S` with
`M(S) >= ceil(u*q/v)`. Record `M(S-1)` and `M(S)` as a checkable minimality
certificate. Do not optimize a physical factor set using outcomes.

For the fixed `[44,44,43]` implicit domains, each slot has at most
`2^(d+1)-1` physical points: x=0 has one point and each other x has at most
two. Report their combined maximum separately; it is not an observed count.
The three x spaces have pairwise intersection `{0}` by the admitted
131-rank normal-basis partition, so their physical-point union has `S-2`
points. On the 4q-order curve, `[4]` has at most four rational preimages of
any projected point, including O. Even giving the leaf a free transported
signed Frobenius action with orbit size 262, retaining logs for **every**
nonzero projected factor point requires at least
`ceil(max(0,S-6)/(4*262))` distinct orbit log classes. Apply this only at
the two minimum-budget thresholds. This conditional log-column floor does
not prove that all columns must be solved by a particular future algorithm;
it is an optimistic screen for the presently specified full-rank factor-base
route. The actual transported action is not free to compute.

The calculator must reject input/hash drift and emit exact integers and
fractions, with decimals only as clearly marked display values. A separate
verifier must brute-check the balanced-product identity for small budgets,
check both threshold minimality brackets, and run the parent `admission.py`
on the fixed-budget balanced slots to compare its hit and rank-one bounds.
STOP on any mismatch. Archive commands and both receipts. An exact bound
under the frozen inputs is `PASS_BOUND`; a mismatch is `FAIL_VERIFICATION`.
Do not infer solver difficulty, a matrix runtime, rho crossover or a
challenge-log result from either outcome. The decision must state which
*bounded route* the envelope rules out and which higher-arity or compressed
route remains unmeasured.
