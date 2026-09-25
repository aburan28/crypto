# Preregistered branch-complete rational recursive-S3 exporter gate

Status: protocol and implementation frozen before the first new exporter or
oracle outcome. [PR #777](https://github.com/aburan28/crypto/pull/777)
proved that finite S3 roots on complete rational x fibres correspond exactly
to rational signed sums, while [PR #774](https://github.com/aburan28/crypto/pull/774)
recorded one n13 witness whose first two factors sum to O and are invisible
to the affine-only chain. This gate exports a concrete one-hot DIMACS CNF
with O states and checks every clause and model path against a separate
full-point group-law oracle. It does not run a SAT solver or price n131.

## Exact transition and target semantics

A factor variable chooses an x from a **complete rational point fibre**
`L(x)`. A prefix state is either `O` or a finite rational x. The first
transition is `L(a)+L(b) -> state`; subsequent transitions are
`state+L(a) -> next_state`. The intended local table is:

- `O+a -> a`, with no O outcome;
- finite `u+a -> c` for every finite root `S3(u,a,c)=0`;
- finite `u+a -> O` iff `u=a`, including `u=a=0`.

For `u=a!=0`, the finite root is the doubling outcome and O is the
inverse-sign outcome. For `u=a=0`, only O is allowed. For `u!=a` with one
zero, the unique finite root is retained. Every finite terminal x accepts
both exact rational target signs because the entire signed tuple may be
negated without changing its factor-x or prefix-state path; terminal O
accepts only the identity. The induction also works across O prefixes:
negating all earlier factors leaves O at O, and `O+L(a)` may take either
allowed sign. If a factor x does not lift, or if a fibre omits a sign,
export is rejected before CNF creation.

Each factor group and prefix-state group is exactly one hot. For every
locally **forbidden** triple, the exporter emits one ternary negative clause;
there are no auxiliary or denominator-clearing variables. For each target
point, `schema.json` gives the final-state assumption literal. The emitted
`base.cnf` plus that unit assumption is the exact point-target query under
the sign-completeness premise. The independent verifier must parse DIMACS,
reconstruct every permitted/forbidden triple using *point sums*, compare the
complete canonical clause multiset, and independently enumerate all signed
factor tuples. It may read only the frozen variable mapping, not the
producer's transition predicate or root cache. It then compares every
`(factor-x tuple, full prefix-state path, terminal state)` from CNF with
the oracle, and checks that each finite path realizes both full terminal
point signs. It separately checks all frozen target assumptions, including
O, and records exceptional paths and masks; absence of a model is called
UNSAT only for these exhaustively checked toy domains.

## Frozen domains and controls

The toy panels are exactly `(n,m,poly)=(2,4,0x7),(3,5,0xb),(4,4,0x13)`
for `E:y²+xy=x³+1`. Every slot admits **all** rational x fibres in that
field. Every full curve point, including O, is an exact target label. All
x-domain tuples and all signed point tuples are enumerated, without
sampling, stopping early, or using a solver. These panels exercise
x=0, doubling, inverse-to-O, O-prefix regeneration, target-O and finite
target signs. The field moduli are checked irreducible, and the point
census is checked against the Koblitz order recurrence.

The fixed regression panel is n13,m5,poly `0x201b`, with the five signed
factor slots from `raw/n13-m5/factors.json` in the #767 raw corpus archive.
It uses all 32 exact Q+T full-point labels in #774's frozen
`evidence/n13-m5/producer/result.json`, plus an explicit O target. This
panel must reproduce the #774 exceptional-only factor-x mask
`[0,0,0,2,1]` for target index 12, `Q+O=(7256,3272)`, whose saved signed
witness starts `(0,1),(0,1)` and passes through O. Its factor-x coordinates
are converted to mask labels using #770's frozen n13 rotated basis only
for that regression assertion. The independent oracle instead reads the
archived signed factor points and computes every full tuple anew. No #774
candidate or oracle output is used as the truth set.

Four deterministic negative controls must be rejected: a nonliftable n3
factor x=2, a sign-incomplete n2 x=1 fibre, a mutated CNF clause that
forbids the n2 `(0,0)->O` transition, and a wrong n2 target-O assumption
literal replaced by a finite x=0 literal. Record the first rejecting check
and its class. These controls distinguish invalid inputs and corrupt
exports from a valid toy UNSAT result.

## Freeze, resource cap, evidence and decision

The producer uses bit-serial GF arithmetic and algebraic quadratic roots
to build the transition tables and DIMACS. The independent verifier uses
the #762 polynomial-product/reduction plus Euclid-inverse group law to
construct local transition tables and a complete signed-tuple oracle.
Freeze exact protocol, producer/verifier/runner/CI bytes, parent source,
#767 raw corpus, #774 result and #770 n13 input/basis bytes by SHA-256 in
`FROZEN.json`. Open a draft PR and pass hash-only CI before outcomes.
Run one cold producer and one independent verifier sequentially, each under
180 s wall and 512 MiB peak RSS, with an external 195 s child stop.
Archive exact DIMACS, schema, all raw path/target rows, complete verifier
result, stdout/stderr, UTC/exits/wall/CPU/RSS, native operation counters,
source hashes, failures/censor receipts and a byte/hash manifest. A timeout,
source drift, mismatch, failed negative control or missing O witness fails
this gate; retain any failed attempt rather than overwrite it.

Pass only when the independent verifier accepts every CNF clause, every
complete model path and target label, the exceptional n13 witness and all
negative controls within caps. Passing admits a **branch-complete toy
semantic exporter** for these rational sign-complete domains. It is not a
direct S6/S7 resultant equivalence, a generic Boolean ANF exporter, a SAT
solver result, n131 relation rank/yield, ECDLP cost S or matched-rho ratio.
The next gate would encode the same #767 fixed full-point systems in the
intended direct/chain solver format, verify returned models and complete
negative proofs independently, then admit bounded solvers.
