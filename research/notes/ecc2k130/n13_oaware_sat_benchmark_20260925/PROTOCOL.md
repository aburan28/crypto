# Preregistered n13-m5 O-aware dense-CNF SAT solver-stage control

Status: freeze-only draft. No timed solver outcomes existed when this protocol,
code, input pins and binary pins were committed. This is a **toy solver-stage**
measurement on [PR #781](https://github.com/aburan28/crypto/pull/781)'s
independently validated rational, complete-fibre, O-aware recursive-S3
one-hot DIMACS. It is not an S6 resultant benchmark, relation generator,
n131 feasibility result, full ECDLP cost, or rho comparison. The dense base
is 1,263 variables, 1,195,344 clauses and 17,630,779 bytes. It is archived
once by #781; this PR pins its SHA-256 and archives only query metadata and
raw child output.

## Hypothesis, truth and admission

A current native CNF engine can solve at least some of the fixed 32 exact
n13-m5 Q+T instances within 15 seconds and 2 GiB sampled child-tree RSS,
returning a model whose chosen factor x values admit a rational signed
tuple summing to the exact full point. Four planted and four complete
point-oracle-negative subgroup Q from [PR #767](https://github.com/aburan28/crypto/pull/767)
are frozen. Every Q has torsion branches T=O,(0,1),(1,0),(1,1), in that
order; all 32 branch targets and final-state assumption literals are exactly
#781's merged schema order, no adaptive labels or target choice. The #767
complete full-point tuple oracle supplies an independent *toy expected
answer*; the runner derives exact Q+T points using an independently coded
field/group law. A SAT model is accepted only after deriving complete
rational point fibres from its five factor x choices and enumerating at
most 2^5 signs to find a tuple whose group sum is exactly Q+T. No saved
#767 witness is used as a SAT certificate. An engine UNSAT response that
agrees with the complete point oracle is **not a separately proof-checked
UNSAT certificate**. No DRAT/LRAT checker is run.

## Frozen solver arms and execution

Three locally installed CNF interfaces are pinned by absolute binary path,
version and SHA-256 in `FROZEN.json`: CryptoMiniSat 5.14.7, Kissat 4.0.4 and
CaDiCaL 3.0.1. MiniSat, WDSat, msolve, F4/F5, crossbred and FES are not
interchangeable CNF arms for this fixed file. The local runner uses Python 3.12.8 with psutil 7.2.2; the preserved pre-outcome Python 3.9 import failure is documented in `evidence/preoutcome_setup_failure.txt`. No `m!` permutation ordering
is added: the five rotated factor slots have distinct x domains. The engine
list is fixed; each arm must first pass two tiny SAT/UNSAT child-interface
smokes showing exit 10/20 and parseable `s`/`v` output. If one fails, retain
its raw failure, exclude it from the matched comparison and report the
missing arm. Never replace the binary after observing a corpus outcome.

For every target, construct one byte-exact query CNF by replacing the base
`p cnf` clause count with count+1 and appending its positive one-unit
assumption. Record its SHA-256, byte count, literal and construction wall;
delete the temporary query file after all three children. Each of the
96 solver invocations is a fresh process that reads the full derived
17.6 MB file. Solver order rotates deterministically by target ordinal
through the frozen engine list; target order is Q0T0...Q7T3. Each child has
15.0 s monotonic wall and 2,147,483,648-byte sampled process-tree RSS caps,
polled every 20 ms with psutil; on cap, kill the process group. Preserve
all stdout/stderr bytes, exact command, UTC interval, exit, cap reason,
child CPU, sampled RSS and Darwin `ru_maxrss` before/after. Measure common
source/schema/oracle preflight, base load, every query-file construction,
child wall, and model parse/independent certificate separately. A killed,
missing-status, malformed, contradictory or incomplete factor model is
censored/invalid, never converted to UNSAT. Preserve any failed attempt
under a new directory; do not overwrite receipts.

The parseable solver convention is an exact `s SATISFIABLE` with exit 10
and a `v` signed-literal model, or exact `s UNSATISFIABLE` with exit 20.
Other statuses are unknown. A SAT model must select exactly one positive
literal in every factor group; the independent fibre/group calculation
must certify the exact point. The child status itself is not a proof of
UNSAT, even when it agrees with #767's exhaustive point oracle.

## Cost and decision gate

Report per branch cold child wall, CPU/RSS, setup and verification. For
Q0..Q3, report the fixed T-order cumulative cost through the first
**verified** SAT witness and the separate cost of running all four
branches. For Q4..Q7, a projected-negative decision requires all four
branch UNSAT statuses and oracle agreement; report all-four cost, or
censor if any branch is unknown/invalid. Include common preflight and
query assembly in cost views, and report the complete eight-Q, 32-branch
per-engine portfolio and the actual 96-child process wall separately.
A per-query maximum wall is not a portfolio cost. Never fit n131 cost from
these toy values. A completed stage comparison requires all 96 queries
and all 32 fixed target statuses for each included engine; otherwise
show the raw censored panel with no favorable ranking. Classify the gate
as `ADMITTED_TOY_SOLVER_STAGE` only if all three smokes pass, no source
drift, no malformed model/contradiction, and every target returns SAT or
UNSAT inside caps with SAT witnesses independently certified and UNSAT
statuses oracle-consistent. Otherwise classify `CENSORED_OR_BLOCKED` and
name every failure. Even an admitted panel says nothing about a full
index-calculus/rho crossover.

Freeze source, merged corpus, exact target order, solver binary hashes,
versions, commands, cap rules, runner, independent verifier, analysis and
archive replay in this draft PR; pass its hash-only CI before the first
child. Host scheduling: run after #784's unequal census and the short
sparse-CNF semantic gate release. Preserve the tiny smoke, every raw
child attempt and failed run, result/analysis receipts and archive-only
CI. The prior #781 base CNF remains a single source archive, linked by
hash in each query receipt.
