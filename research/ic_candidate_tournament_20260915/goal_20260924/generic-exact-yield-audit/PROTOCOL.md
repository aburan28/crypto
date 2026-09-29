# Exact natural-yield audit of the disclosed recovery pilot

This protocol is registered **before** analyzing the raw three-summand
ordinary-query stream from [PR #965](https://github.com/aburan28/crypto/pull/965).
The immutable [panel](panel.json) pins its 4,914,932-byte source-bound
`evidence.tar.gz` by SHA-256
`e0ce19fc28c58e2dc1cae9649a16af74099016fff8183e5ce69f65b39c804f02`.
The original source commit, five disclosed public points, base recipe, query
seed and worker limits remain unchanged. No target, worker or solver is rerun.
This audit does not use any sealed confirmation set.

**Question.** For every natural ordinary query actually recorded by F4, F5,
native-XOR SAT and CNF SAT, did a three-point decomposition into the exact
factor base exist? The worker's `incomplete` verdict means its node or
conflict budget was exhausted, so the previous relation-yield table cannot
distinguish a hard solvable query from a query with no relation. The exact
group oracle provides that distinction for the recorded finite workloads.

Before classification, independently replay the original source/build,
factor-base, query-law, group relation, solver-dispatch, matrix and phase
audits for all sixteen reports. The four jobs that timed out without a report
retain unknown query counts and unknown yield. For each cell, decode the
verified **geometric** base in order and insert every pair `i <= j` into an
index keyed by the independently computed group sum. For each recorded
ordinary query `[a]G`, check all base points `P_k`: if
`[a]G - P_k` is in the pair index, independently re-add the three points and
record a valid witness. This exhausts all three-summand index triples,
including repeated points, without trusting either the worker's row solver
or its PDP verdict. Keep the first third-index/pair-index witness only;
the classification is existence, not a rank-maximizing oracle.

The analyzer must classify **all** reported attempts in the frozen
cell-by-solver order and compare each observed witness or proved-UNSAT verdict
against exact group existence. It records a solvable-but-bounded-incomplete
count for each arm, plus the exact yes/no labels for every query. It fails
closed on a mismatch rather than discarding a row. The same ordered query
stream across reported arms is checked, but a worker stopped early may have
only a prefix. No wall time, amortized throughput or speedup is inferred from
the Python oracle. A planted witness is not introduced. The output is a
diagnostic on previously disclosed points, not a fresh-target yield estimate
or SAT qualification.

The runner, panel hash, archive hash and this protocol are committed in a
follow-on PR before the analysis is executed. Run the registered analyzer
once to a new output directory, retain its complete query labels and digest,
then report every arm including censored timeouts. If the exact audit shows
SAT missed feasible queries, the next solver experiment needs a separately
frozen SAT backend or encoding and complete one-target admission before fresh
competitive points are spent.
