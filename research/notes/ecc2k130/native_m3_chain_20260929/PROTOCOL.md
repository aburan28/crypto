# Native degree-263 three-summand implicit-domain chain

Status: preregistered in the first PR commit before candidate source or held
evaluation. The initial frozen panel passed its declared controls, but its
selected SAT witness was the x=0 triple, exercising only inverse and copy.
The expanded semantic control below was specified **after inspecting that
result and before its own source revision or run**; the initial raw archive
is retained in `evidence/initial`. This is the next semantic and
representation gate after the generic full-point edge in
[#973](https://github.com/aburan28/crypto/pull/973). It does not estimate
natural PDP yield, relation rank, a full logarithm cost, or crossover with rho.
The conditional paired native/pullback experiment remains governed by
[`ecc2k130_leaf_native_m3_gate_20260926/PROTOCOL.md`](../../../ecc2k130_leaf_native_m3_gate_20260926/PROTOCOL.md).

## Hypothesis and exact construction

The complete `E_(a,b): y²+xy=x³+a x²+b` point edge from #973 can be bound into
a three-factor prefix chain without enumerating an infeasible native base.
Every factor is a *finite* physical curve point `(x_i,y_i)`, with a free
`y_i` word checked by the edge and an `x_i` word wired from independent
Boolean selectors and a fixed GF(2)-linear basis. Bind `F0+F1=S2` and
`S2+F2=SUM` through two copies of the same generic edge, using shared full
point roles, independent slope words and a single hash-consed DAG. Require
both edges and a canonical pinned `SUM` target. Include all curve equation,
infinity, inverse, doubling and ordinary-addition constraints. The source
edge and prior source-curve chain builder are immutable dependencies, not
modified inputs.

A positive solver model is valid only after (i) every exported CNF clause is
checked against the complete assignment, (ii) the input bits and factor
selectors are decoded with no missing or contradictory Boolean values,
(iii) the two edge outputs hold in the native DAG, and (iv) an independent
bit-polynomial field and group law confirm each factor is in its declared
physical domain, `S2=F0+F1`, and `SUM=S2+F2=target`. A solver status or
internal DAG result alone is insufficient.

## Fixed small-field solver panel

Use GF(2^5) modulus `0x25`, curve `a=0,b=1`, and the three ordered one-bit
x-domain bases `[[1],[2],[4]]`. Independently enumerate each physical factor
set, all ordered factor triples and their full point sums. Select the first
non-infinity supported target in canonical `(infinity,x,y)` order as the
positive; select the first rational point absent from support as the
negative. If either does not exist, STOP; do not alter bases or target rule.
Export one complete Tseitin CNF per target using the frozen parent exporter.
Run a single-thread CryptoMiniSat 5 child with 120-second external cap for
each target, recording the exact binary SHA-256, version, command, stdout,
stderr, exit, wall, CPU, RSS and CNF SHA-256. The positive must return SAT
with a fully lifted point witness. The negative must return UNSAT, and the
independent exhaustive group oracle must confirm it is absent. Because the
solver emits no independently checked UNSAT proof here, label the negative
`SOLVER_UNSAT_ORACLE_CONFIRMED_TOY`, not general checked UNSAT.
A timeout, unsupported solver output, incomplete assignment or replay
mismatch is STOP, not a negative PDP result.
The containing toy producer and verifier have 300-second reported wall
acceptance caps and 330-second external watchdogs, each with a 1-GiB peak-RSS
acceptance cap.

The expanded panel must evaluate the native three-factor DAG for **every
ordered physical factor triple**: the independent exact `S2` and `SUM` with
exact slopes must be accepted. Replacing `SUM` with the first different
rational point in canonical order while retaining all other inputs must be
rejected. Archive every row and a stream SHA-256, not only aggregate counts.
Count the branch type at each edge and require the initial panel's 27 triples,
including 24 generic/generic triples, to be independently replayed. Select
one additional solver target: the canonical first supported point for which
**every** witness tuple uses generic addition at both edges. If none exists,
STOP. Its complete SAT model must lift to a generic/generic witness. Keep the
original first-positive and first-negative targets and the same binary/caps.
These three toy solver queries are disclosed controls, not a natural-target
PDP yield sample. Only the expanded panel can support
`PASS_SEMANTIC_AND_CAPACITY`.

## Frozen n131 capacity panel

Read the literal `a=0,b` values for both nonconjugate degree-263 lines
`[1,0]` and `[1,4]` from the archived exact leaf smoke, and require its
independent replay PASS. In GF(2^131) with polynomial modulus
`0x800000000000000000000000000002007`, take normal element beta `3` and
three basis dimensions `[44,44,43]`. Slot `i` has vectors
`beta^(2^(3j+i))` for `j<d_i`; require all 131 conjugates have GF(2)
rank 131, each slot has its declared rank, and the union has rank 131.
The x domain is implicit: do not enumerate its approximately 2^43 members.
Build each leaf chain in a separate process. Count the exact generic and
canonical-infinity-target DIMACS size under a 2,000,000-node DAG cap and
256-MiB target CNF cap. Each child has a 600-second reported wall acceptance cap and a
2-GiB peak-RSS acceptance cap; an external watchdog kills it after 630 seconds.
RSS is checked after exit, not represented
as an OS-enforced limit. Preserve full or partial counts and errors. A cap
is `CAPACITY_CENSORED`, not evidence of PDP hardness. Do not submit n131
CNF to a solver in this gate or infer attack speed from representation size.

Commit the candidate, producer and independent verifier in a second revision;
commit their SHA-256 source/input lock in a separate third revision before
any held panel. Archive the exact raw
toy CNFs, model and solver streams; two leaf count receipts; frozen source and
input hashes; child process resource records; all failures; and an
independent replay. CI must validate the archive without depending on the
local solver executable. The decision is `PASS_SEMANTIC_AND_CAPACITY` only
when the toy SAT/negative controls and both n131 capacity arms pass all
checks within caps. Any mismatch is FAIL; timeout, missing evidence, loader
failure, or resource cap is STOP or censored as specified, never silently
omitted. Even PASS authorizes only a separately frozen natural-target
solver/yield experiment followed by the exact native/pullback parity and
full-rank cost comparisons in the parent m≥3 protocol.
